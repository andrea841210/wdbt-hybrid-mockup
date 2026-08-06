
import csv
import json
import re
from io import StringIO
from pathlib import Path

import pandas as pd
import streamlit as st


# -----------------------------
# Basic settings
# -----------------------------
st.set_page_config(
    page_title="WDBT Hybrid Workflow Mockup",
    layout="wide",
)

DATA_DIR = Path("data")
ASSET_DIR = Path("assets")


# -----------------------------
# Data loading helpers
# -----------------------------
def read_csv_with_encoding_fallback(path: Path, **read_csv_kwargs) -> pd.DataFrame:
    """Read a local CSV with a deterministic encoding fallback."""
    raw = path.read_bytes()
    for encoding in ("utf-8-sig", "utf-8", "cp950", "big5"):
        try:
            text = raw.decode(encoding)
            return pd.read_csv(StringIO(text), **read_csv_kwargs)
        except UnicodeDecodeError:
            continue
    raise ValueError(f"Unable to decode CSV: {path.name}")


@st.cache_data
def load_pan_syndrome_master() -> pd.DataFrame:
    path = DATA_DIR / "pan_syndrome_master.csv"
    df = read_csv_with_encoding_fallback(
        path,
        sep=None,
        engine="python",
        dtype=str,
        keep_default_na=False,
    )
    expected_cols = {
        "pan_syndrome_id",
        "display_order",
        "section_code",
        "section_name",
        "pan_syndrome_zh",
        "primary_phenotype_category",
        "secondary_phenotype_category",
    }
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(f"pan_syndrome_master.csv missing columns: {missing}")
    df["secondary_phenotype_category"] = df["secondary_phenotype_category"].fillna("")
    df["display_order"] = pd.to_numeric(df["display_order"], errors="raise")
    return df.sort_values("display_order")


@st.cache_data
def load_category_to_hpo() -> pd.DataFrame:
    path = DATA_DIR / "phenotype_category_to_hpo.csv"
    df = pd.read_csv(path, sep=None, engine="python")
    expected_cols = {"phenotype_category", "hpo_id", "hpo_term"}
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(f"phenotype_category_to_hpo.csv missing columns: {missing}")
    return df


@st.cache_data
def load_gene_disease_phenotype_map() -> pd.DataFrame:
    path = DATA_DIR / "gene_disease_phenotype_map.csv"
    df = read_csv_with_encoding_fallback(
        path,
        sep=",",
        engine="python",
        dtype=str,
        keep_default_na=False,
    )
    canonical_cols = [
        "disease_id",
        "gene",
        "disease_name",
        "omim_id",
        "phenotype_concern_list",
        "inheritance",
    ]
    expected_cols = set(canonical_cols)
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(f"gene_disease_phenotype_map.csv missing columns: {missing}")
    df = df[canonical_cols].copy()
    for col in canonical_cols:
        df[col] = df[col].astype(str).str.strip()
    df["gene"] = df["gene"].str.upper()

    required_nonblank = [
        "disease_id",
        "gene",
        "disease_name",
        "omim_id",
        "phenotype_concern_list",
    ]
    blank_counts = {col: int((df[col] == "").sum()) for col in required_nonblank}
    blank_counts = {col: count for col, count in blank_counts.items() if count}
    if blank_counts:
        raise ValueError(f"gene_disease_phenotype_map.csv blank required values: {blank_counts}")
    return df


@st.cache_data
def load_function2_static() -> pd.DataFrame:
    path = DATA_DIR / "function2_risk_strata_static.csv"
    if not path.exists():
        return pd.DataFrame(
            [
                {
                    "Level": 1,
                    "Color": "🔴 Red",
                    "Variant": "GJB2 c.109G>A",
                    "Concern": "Hearing loss",
                    "Disease ID": "D001",
                    "Disease term": "Xxxx",
                },
                {
                    "Level": 3,
                    "Color": "🟡 Yellow",
                    "Variant": "SLC26A4 c.919-2A>G",
                    "Concern": "Hearing loss",
                    "Disease ID": "D050",
                    "Disease term": "Xxxx",
                },
            ]
        )
    return pd.read_csv(path, sep=None, engine="python")


@st.cache_data
def load_reporting_package_mock() -> dict:
    path = DATA_DIR / "reporting_package_mock.json"
    if not path.exists():
        return {
            "sample_id": "TS260508003",
            "status": "static preview",
            "outputs": [
                "Pan-syndrome recall",
                "Risk strata ranking",
                "Reporting package preview",
            ],
        }
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


# -----------------------------
# Function 1 helpers
# -----------------------------
def split_categories(primary: str, secondary: str) -> list[str]:
    """
    Convert primary and secondary phenotype category fields into a clean list.
    Secondary categories may be separated by semicolons or commas.
    """
    categories = []

    if isinstance(primary, str) and primary.strip():
        raw_primary = primary.replace("；", ";")
        for item in raw_primary.replace(",", ";").split(";"):
            item = item.strip()
            if item:
                categories.append(item)

    if isinstance(secondary, str) and secondary.strip():
        raw_secondary = secondary.replace("；", ";")
        for item in raw_secondary.replace(",", ";").split(";"):
            item = item.strip()
            if item:
                categories.append(item)

    # Deduplicate while preserving order
    seen = set()
    clean = []
    for cat in categories:
        if cat not in seen:
            clean.append(cat)
            seen.add(cat)
    return clean


def split_semicolon_list(value: str) -> list[str]:
    """
    Split PM mapping phenotype_concern_list.
    """
    if not isinstance(value, str) or not value.strip():
        return []
    raw = value.replace("；", ";").replace(",", ";")
    items = [item.strip() for item in raw.split(";") if item.strip()]
    seen = set()
    clean = []
    for item in items:
        if item not in seen:
            clean.append(item)
            seen.add(item)
    return clean


def build_function1_outputs(
    selected_ids: list[str],
    sample_id: str,
    pan_df: pd.DataFrame,
    hpo_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selected_pan = pan_df[pan_df["pan_syndrome_id"].isin(selected_ids)].copy()

    category_rows = []
    hpo_rows = []

    for _, row in selected_pan.iterrows():
        ps_id = row["pan_syndrome_id"]
        ps_text = row["pan_syndrome_zh"]
        categories = split_categories(
            row["primary_phenotype_category"],
            row["secondary_phenotype_category"],
        )

        for cat in categories:
            category_rows.append(
                {
                    "pan_syndrome_id": ps_id,
                    "pan_syndrome_zh": ps_text,
                    "phenotype_category": cat,
                }
            )

            matched = hpo_df[hpo_df["phenotype_category"] == cat]
            for _, hpo_row in matched.iterrows():
                hpo_rows.append(
                    {
                        "pan_syndrome_id": ps_id,
                        "pan_syndrome_zh": ps_text,
                        "phenotype_category": cat,
                        "hpo_id": hpo_row["hpo_id"],
                        "hpo_term": hpo_row["hpo_term"],
                    }
                )

    category_expansion = pd.DataFrame(category_rows)
    hpo_mapping = pd.DataFrame(hpo_rows)

    if hpo_mapping.empty:
        geneyx_ready = pd.DataFrame(
            columns=[
                "sample_id",
                "hpo_id",
                "hpo_term",
                "source_pan_syndrome_ids",
                "source_phenotype_categories",
            ]
        )
    else:
        geneyx_ready = (
            hpo_mapping.groupby(["hpo_id", "hpo_term"], as_index=False)
            .agg(
                source_pan_syndrome_ids=(
                    "pan_syndrome_id",
                    lambda x: "; ".join(sorted(set(x))),
                ),
                source_phenotype_categories=(
                    "phenotype_category",
                    lambda x: "; ".join(sorted(set(x))),
                ),
            )
        )
        geneyx_ready.insert(0, "sample_id", sample_id)

    return selected_pan, category_expansion, hpo_mapping, geneyx_ready


def make_csv_download(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False, encoding="utf-8-sig").encode("utf-8-sig")


# -----------------------------
# Function 2 helpers
# -----------------------------
RELEVANCE_HIGH_BUCKETS = [
    "Dominant HET",
    "Recessive Compound HET",
    "Recessive HOM",
    "Mitochondria",
]
P_LP_CLASSIFICATIONS = {"pathogenic", "likelypathogenic"}


def decode_uploaded_tsv(uploaded_file) -> str:
    """
    Decode uploaded TSV bytes with common encodings.
    """
    raw = uploaded_file.getvalue()
    for encoding in ["utf-8-sig", "utf-8", "cp950", "big5", "latin1"]:
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", errors="replace")


def parse_geneyx_tsv(uploaded_file, inheritance_bucket: str) -> tuple[dict, pd.DataFrame]:
    """
    Parse one Geneyx TSV file.

    Expected structure:
    #Main Sample:...
    #Analysis:...
    #Patient:...
    Relevance    Pathogenic    Note    Location    Gene ...
    variant rows...
    """
    text = decode_uploaded_tsv(uploaded_file)
    lines = text.splitlines()

    metadata = {
        "sample_id": "",
        "analysis_id": "",
        "patient_id": "",
        "inheritance_bucket": inheritance_bucket,
    }

    header_idx = None

    for idx, line in enumerate(lines):
        clean_line = line.strip()

        if clean_line.startswith("#Main Sample:"):
            metadata["sample_id"] = clean_line.replace("#Main Sample:", "", 1).strip()
        elif clean_line.startswith("#Analysis:"):
            metadata["analysis_id"] = clean_line.replace("#Analysis:", "", 1).strip()
        elif clean_line.startswith("#Patient:"):
            metadata["patient_id"] = clean_line.replace("#Patient:", "", 1).strip()

        if clean_line.startswith("Relevance\tPathogenic\t"):
            header_idx = idx
            break

    if header_idx is None:
        return metadata, pd.DataFrame()

    body_text = "\n".join(lines[header_idx:]).strip()
    if not body_text:
        return metadata, pd.DataFrame()

    df = pd.read_csv(
        StringIO(body_text),
        sep="\t",
        dtype=str,
        keep_default_na=False,
        quoting=csv.QUOTE_NONE,
        engine="python",
    )

    df = df.dropna(how="all")
    if df.empty:
        return metadata, df

    df.insert(0, "inheritance_bucket", inheritance_bucket)
    df.insert(0, "patient_id", metadata["patient_id"])
    df.insert(0, "analysis_id", metadata["analysis_id"])
    df.insert(0, "sample_id", metadata["sample_id"])

    return metadata, df


def normalize_pathogenicity(value: str) -> str:
    """Normalize Geneyx pathogenicity labels for exact P/LP matching."""
    return re.sub(r"[^a-z]", "", str(value).lower())


def is_geneyx_interpreted_candidate(df: pd.DataFrame) -> pd.Series:
    """
    Geneyx already performs an AI/HPO-based interpretation layer.

    Function 2 decision tree should only use rows where both:
    - Relevance is not blank
    - Pathogenic is exactly Pathogenic or Likely Pathogenic after normalization

    Other values such as VUS/UncertainSignificance, Benign, or a merely
    non-blank value must not enter risk stratification.
    """
    if df.empty or "Relevance" not in df.columns or "Pathogenic" not in df.columns:
        return pd.Series(False, index=df.index)

    relevance_filled = df["Relevance"].astype(str).str.strip() != ""
    pathogenic_is_plp = df["Pathogenic"].map(normalize_pathogenicity).isin(
        P_LP_CLASSIFICATIONS
    )
    return relevance_filled & pathogenic_is_plp


def parse_all_geneyx_files(
    uploaded_files: dict,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Parse all required Geneyx TSV files.

    Returns:
    bucket_summary
    all_candidates
    interpreted_plp_candidates
    """
    all_tables = []
    summary_rows = []

    for bucket, uploaded_file in uploaded_files.items():
        metadata, df = parse_geneyx_tsv(uploaded_file, bucket)

        raw_variant_count = len(df) if not df.empty else 0

        if raw_variant_count > 0:
            interpreted_mask = is_geneyx_interpreted_candidate(df)
            interpreted_count = int(interpreted_mask.sum())
        else:
            interpreted_count = 0

        summary_rows.append(
            {
                "inheritance_bucket": bucket,
                "sample_id": metadata.get("sample_id", ""),
                "analysis_id": metadata.get("analysis_id", ""),
                "patient_id": metadata.get("patient_id", ""),
                "raw_variant_count": raw_variant_count,
                "interpreted_plp_count": interpreted_count,
                "status": "has interpreted P/LP" if interpreted_count > 0 else "no interpreted P/LP",
            }
        )

        if raw_variant_count > 0:
            all_tables.append(df)

    bucket_summary = pd.DataFrame(summary_rows)

    if not all_tables:
        all_candidates = pd.DataFrame()
        interpreted_candidates = pd.DataFrame()
        return bucket_summary, all_candidates, interpreted_candidates

    all_candidates = pd.concat(all_tables, ignore_index=True)

    interpreted_mask = is_geneyx_interpreted_candidate(all_candidates)
    interpreted_candidates = all_candidates[interpreted_mask].copy()

    return bucket_summary, all_candidates, interpreted_candidates


def split_mapping_coverage(
    interpreted_candidates: pd.DataFrame,
    mapping_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split P/LP candidates by whether the active disease map contains the gene."""
    if interpreted_candidates.empty or "Gene" not in interpreted_candidates.columns:
        return pd.DataFrame(), interpreted_candidates.copy()

    active_genes = set(mapping_df["gene"].astype(str).str.upper().str.strip())
    covered_mask = (
        interpreted_candidates["Gene"].astype(str).str.upper().str.strip().isin(active_genes)
    )
    return (
        interpreted_candidates[covered_mask].copy(),
        interpreted_candidates[~covered_mask].copy(),
    )


def simplify_geneyx_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Reduce Geneyx interpreted P/LP candidates to traceable input columns.

    The inheritance bucket comes from the uploaded Geneyx TSV filename and is
    the inheritance-mode input used by the decision tree.
    """
    display_cols = [
        "sample_id",
        "analysis_id",
        "patient_id",
        "inheritance_bucket",
        "Relevance",
        "Pathogenic",
        "Gene",
        "HGVSC",
        "HGVSP",
        "Zygosity",
        "ACMG",
    ]

    if df.empty:
        return pd.DataFrame(columns=display_cols)

    out = df.copy()
    for col in display_cols:
        if col not in out.columns:
            out[col] = ""

    return out[display_cols]


def build_selected_category_context(
    selected_pan: pd.DataFrame,
) -> tuple[list[str], dict[str, list[str]]]:
    """
    Merge primary and secondary categories with equal decision-tree weight.

    Source labels are retained for QA trace only; they never change the tier.
    """
    selected_categories = []
    category_sources: dict[str, list[str]] = {}

    ordered_pan = selected_pan.sort_values("display_order")
    for _, row in ordered_pan.iterrows():
        ps_id = str(row.get("pan_syndrome_id", "")).strip()
        for field, source_type in (
            ("primary_phenotype_category", "primary"),
            ("secondary_phenotype_category", "secondary"),
        ):
            for category in split_semicolon_list(row.get(field, "")):
                if category not in category_sources:
                    selected_categories.append(category)
                    category_sources[category] = []
                source = f"{ps_id}.{source_type}"
                if source not in category_sources[category]:
                    category_sources[category].append(source)

    return selected_categories, category_sources


def aggregate_gene_disease_maps(gene_maps: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse multiple OMIM provenance rows into one disease-level map.

    Risk identity is variant + inheritance bucket + disease_id. OMIM IDs and
    mapping inheritance remain visible as reference-only provenance.
    """
    if gene_maps.empty:
        return pd.DataFrame()

    rows = []
    for disease_id, group in gene_maps.groupby("disease_id", sort=False):
        disease_names = list(dict.fromkeys(group["disease_name"].astype(str)))
        omim_ids = list(dict.fromkeys(group["omim_id"].astype(str)))
        inheritances = list(dict.fromkeys(group["inheritance"].astype(str)))

        concerns = []
        for concern_list in group["phenotype_concern_list"]:
            for concern in split_semicolon_list(concern_list):
                if concern not in concerns:
                    concerns.append(concern)

        rows.append(
            {
                "disease_id": disease_id,
                "gene": str(group["gene"].iloc[0]),
                "disease_name": " / ".join(disease_names),
                "omim_ids": "; ".join(omim_ids),
                "mapping_inheritance_reference": "; ".join(inheritances),
                "phenotype_concern_list": "; ".join(concerns),
            }
        )

    return pd.DataFrame(rows)


def build_phenotype_relevance_and_risk_preview(
    mapped_candidates: pd.DataFrame,
    mapping_df: pd.DataFrame,
    selected_categories: list[str],
    category_sources: dict[str, list[str]] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Connect Geneyx P/LP candidates to the active disease mapping,
    then apply the fixed four-tier rules.
    """
    if mapped_candidates.empty:
        columns = [
            "risk_level",
            "color",
            "risk_tier",
            "concern",
            "phenotype_relevance",
            "matched_phenotype_categories",
            "matched_phenotype_sources",
            "gene",
            "variant",
            "disease_id",
            "disease_name",
            "omim_ids",
            "mapping_inheritance_reference",
            "inheritance_bucket",
            "relevance",
            "pathogenic",
            "rule_trace",
        ]
        return pd.DataFrame(columns=columns), pd.DataFrame(columns=columns)

    category_sources = category_sources or {}
    selected_set = set(selected_categories)
    rows = []

    for _, variant in mapped_candidates.iterrows():
        gene = str(variant.get("Gene", "")).upper().strip()
        gene_maps = mapping_df[mapping_df["gene"] == gene].copy()
        gene_maps = aggregate_gene_disease_maps(gene_maps)
        if gene_maps.empty:
            continue

        for _, mapping in gene_maps.iterrows():
            phenotype_concerns = split_semicolon_list(mapping.get("phenotype_concern_list", ""))
            matched_categories = sorted(selected_set.intersection(set(phenotype_concerns)))
            phenotype_relevance = "Yes" if matched_categories else "No"
            matched_sources = []
            for category in matched_categories:
                for source in category_sources.get(category, []):
                    source_label = f"{source}: {category}"
                    if source_label not in matched_sources:
                        matched_sources.append(source_label)

            inheritance_bucket = str(variant.get("inheritance_bucket", ""))
            relevance = str(variant.get("Relevance", "")).strip()

            # Only Geneyx-interpreted P/LP rows reach this function. Primary
            # and secondary phenotype categories have identical match weight.
            if inheritance_bucket in RELEVANCE_HIGH_BUCKETS:
                if phenotype_relevance == "Yes":
                    risk_level = 1
                    color = "🔴 Red"
                    risk_tier = "High Attention"
                    rule_trace = "P/LP compatible bucket + phenotype relevance Yes"
                else:
                    risk_level = 2
                    color = "🟠 Orange"
                    risk_tier = "Targeted Follow-up"
                    rule_trace = "P/LP compatible bucket + phenotype relevance No"
            elif inheritance_bucket == "Recessive HET":
                if phenotype_relevance == "Yes":
                    risk_level = 3
                    color = "🟡 Yellow"
                    risk_tier = "Routine Monitoring"
                    rule_trace = "P/LP Recessive HET + phenotype relevance Yes"
                else:
                    risk_level = 4
                    color = "🔵 Blue"
                    risk_tier = "General Awareness"
                    rule_trace = "P/LP Recessive HET + phenotype relevance No"
            else:
                risk_level = 4
                color = "🔵 Blue"
                risk_tier = "General Awareness"
                rule_trace = "No MVP bucket rule matched"

            hgvsc = str(variant.get("HGVSC", "")).strip()
            hgvsp = str(variant.get("HGVSP", "")).strip()
            aa = str(variant.get("AA", "")).strip()
            variant_label = hgvsc or hgvsp or aa or str(variant.get("Location", "")).strip()

            rows.append(
                {
                    "risk_level": risk_level,
                    "color": color,
                    "risk_tier": risk_tier,
                    "concern": "; ".join(phenotype_concerns),
                    "phenotype_relevance": phenotype_relevance,
                    "matched_phenotype_categories": "; ".join(matched_categories),
                    "matched_phenotype_sources": "; ".join(matched_sources),
                    "gene": gene,
                    "variant": variant_label,
                    "disease_id": mapping.get("disease_id", ""),
                    "disease_name": mapping.get("disease_name", ""),
                    "omim_ids": mapping.get("omim_ids", ""),
                    "mapping_inheritance_reference": mapping.get(
                        "mapping_inheritance_reference", ""
                    ),
                    "inheritance_bucket": inheritance_bucket,
                    "relevance": relevance,
                    "pathogenic": variant.get("Pathogenic", ""),
                    "rule_trace": rule_trace,
                }
            )

    relevance_table = pd.DataFrame(rows)

    if relevance_table.empty:
        return relevance_table, relevance_table

    risk_preview = relevance_table.sort_values(
        by=["risk_level", "gene", "disease_id", "variant"],
        ascending=[True, True, True, True],
    ).reset_index(drop=True)

    return relevance_table, risk_preview


# -----------------------------
# Sidebar
# -----------------------------
st.sidebar.title("WDBT Mockup v0.5")
mode = st.sidebar.radio(
    "Mode",
    ["Function 1｜HPO Translation", "Function 2｜Risk Strata Preview"],
)

st.sidebar.caption(
    "Hybrid mockup only. Not a production system, not a clinical decision support tool."
)
st.sidebar.link_button(
    "View app.py on GitHub",
    "https://github.com/andrea841210/wdbt-hybrid-mockup/blob/main/app.py",
    use_container_width=True,
)

pdf_path = ASSET_DIR / "WDBT_UI_all-2.pdf"
if pdf_path.exists():
    with open(pdf_path, "rb") as f:
        st.sidebar.download_button(
            "Download UI reference PDF",
            data=f,
            file_name="WDBT_UI_all-2.pdf",
            mime="application/pdf",
        )


# -----------------------------
# App title
# -----------------------------
st.markdown("## WDBT Hybrid Workflow Mockup")
st.caption(
    "A hybrid workflow mockup for WDBT input translation and post-Geneyx risk-strata preview."
)


# -----------------------------
# Function 1
# -----------------------------
if mode == "Function 1｜HPO Translation":
    st.markdown("### Function 1｜兒基安泛性症狀轉譯層")
    st.caption("PGSafe DB_v1 · HPO database_v1")
    st.caption(
        "Function 1 performs input preparation only: "
        "pan-syndrome → phenotype category → HPO / Geneyx-ready input. "
        "No interpretation, disease inference, or risk stratification is performed."
    )

    pan_df = load_pan_syndrome_master()
    hpo_df = load_category_to_hpo()

    left, right = st.columns([1, 2])

    with left:
        st.subheader("Input")

        display_options = {
            f"[{row.section_code}] {row.pan_syndrome_id}｜{row.pan_syndrome_zh}": row.pan_syndrome_id
            for _, row in pan_df.iterrows()
        }

        default_labels = [
            label
            for label, ps_id in display_options.items()
            if ps_id in ["PS001", "PS003", "PS005"]
        ]

        selected_labels = st.multiselect(
            "Pan-syndrome ID",
            options=list(display_options.keys()),
            default=default_labels,
        )
        selected_ids = [display_options[label] for label in selected_labels]

        physician_note = st.text_input("醫師補充描述", value="Headache")
        top_n = st.number_input("Top N to show HPO", min_value=1, max_value=500, value=200)
        sample_id = st.text_input("Sample ID", value="TS260508003")

        run_translation = st.button("Run Translation", type="primary")

    with right:
        if not selected_ids:
            st.warning("Please select at least one pan-syndrome item.")
        elif run_translation:
            selected_pan, category_expansion, hpo_mapping, geneyx_ready = build_function1_outputs(
                selected_ids=selected_ids,
                sample_id=sample_id,
                pan_df=pan_df,
                hpo_df=hpo_df,
            )

            st.subheader("1. Pan-syndrome list")
            st.dataframe(
                selected_pan[
                    [
                        "pan_syndrome_id",
                        "section_code",
                        "pan_syndrome_zh",
                        "primary_phenotype_category",
                        "secondary_phenotype_category",
                    ]
                ],
                use_container_width=True,
                hide_index=True,
            )

            st.subheader("2. Phenotype category expansion")
            st.dataframe(
                category_expansion,
                use_container_width=True,
                hide_index=True,
            )

            st.subheader(f"3. HPO mapping preview (n={min(top_n, len(hpo_mapping))})")
            st.dataframe(
                hpo_mapping.head(int(top_n)),
                use_container_width=True,
                hide_index=True,
            )

            st.subheader("4. Geneyx-ready HPO input package")
            st.dataframe(
                geneyx_ready,
                use_container_width=True,
                hide_index=True,
            )

            st.download_button(
                "Download Geneyx-ready HPO input package",
                data=make_csv_download(geneyx_ready),
                file_name=f"{sample_id}_geneyx_ready_hpo_input.csv",
                mime="text/csv",
                type="primary",
            )


# -----------------------------
# Function 2
# -----------------------------
else:
    st.markdown("### Function 2｜分子證據與風險分層預覽")
    st.caption(
        "Active scope: all genes present in gene_disease_phenotype_map.csv. "
        "Only Geneyx-interpreted P/LP rows are eligible: Relevance must be filled and "
        "Pathogenic must equal Pathogenic or Likely Pathogenic. Primary and secondary "
        "pan-syndrome categories have equal matching weight."
    )


    pan_df = load_pan_syndrome_master()
    mapping_df = load_gene_disease_phenotype_map()

    left, right = st.columns([1, 2])

    with left:
        st.subheader("Input")

        st.markdown("**Upload Geneyx TSV files**")

        dominant_het_file = st.file_uploader(
            "Dominant HET.tsv",
            type=["tsv", "txt"],
            key="dominant_het_file",
        )
        mitochondria_file = st.file_uploader(
            "Mitochondria.tsv",
            type=["tsv", "txt"],
            key="mitochondria_file",
        )
        recessive_compound_het_file = st.file_uploader(
            "Recessive Compound HET.tsv",
            type=["tsv", "txt"],
            key="recessive_compound_het_file",
        )
        recessive_het_file = st.file_uploader(
            "Recessive HET.tsv",
            type=["tsv", "txt"],
            key="recessive_het_file",
        )
        recessive_hom_file = st.file_uploader(
            "Recessive HOM.tsv",
            type=["tsv", "txt"],
            key="recessive_hom_file",
        )

        st.markdown("---")

        display_options = {
            f"[{row.section_code}] {row.pan_syndrome_id}｜{row.pan_syndrome_zh}": row.pan_syndrome_id
            for _, row in pan_df.iterrows()
        }

        default_labels = [
            label
            for label, ps_id in display_options.items()
            if ps_id in ["PS001", "PS003", "PS005"]
        ]

        selected_labels = st.multiselect(
            "Pan-syndrome ID",
            options=list(display_options.keys()),
            default=default_labels,
            key="function2_ps",
        )
        selected_ids = [display_options[label] for label in selected_labels]

        physician_note = st.text_input("醫師補充描述", value="Headache", key="function2_note")
        sample_id = st.text_input("Sample ID", value="TS260508003", key="function2_sample")

        run_preview = st.button("Run Risk Preview", type="primary")

    with right:
        uploaded_files = {
            "Dominant HET": dominant_het_file,
            "Mitochondria": mitochondria_file,
            "Recessive Compound HET": recessive_compound_het_file,
            "Recessive HET": recessive_het_file,
            "Recessive HOM": recessive_hom_file,
        }

        missing_files = [bucket for bucket, file in uploaded_files.items() if file is None]

        if run_preview:
            if missing_files:
                st.error(
                    "Please upload all five Geneyx TSV files before running MVP risk preview. "
                    f"Missing: {', '.join(missing_files)}"
                )
            elif not selected_ids:
                st.warning("Please select at least one pan-syndrome item.")
            else:
                selected_pan = pan_df[pan_df["pan_syndrome_id"].isin(selected_ids)].copy()
                selected_categories, category_sources = build_selected_category_context(
                    selected_pan
                )

                bucket_summary, all_candidates, interpreted_candidates = parse_all_geneyx_files(
                    uploaded_files
                )
                mapped_candidates, unmapped_candidates = split_mapping_coverage(
                    interpreted_candidates,
                    mapping_df,
                )
                interpreted_display = simplify_geneyx_candidates(interpreted_candidates)
                mapped_display = simplify_geneyx_candidates(mapped_candidates)
                relevance_table, risk_preview = build_phenotype_relevance_and_risk_preview(
                    mapped_candidates=mapped_candidates,
                    mapping_df=mapping_df,
                    selected_categories=selected_categories,
                    category_sources=category_sources,
                )

                st.subheader("1. Selected pan-syndrome context")
                st.dataframe(
                    selected_pan[
                        [
                            "pan_syndrome_id",
                            "section_code",
                            "pan_syndrome_zh",
                            "primary_phenotype_category",
                            "secondary_phenotype_category",
                        ]
                    ],
                    use_container_width=True,
                    hide_index=True,
                )



                st.subheader("2. Geneyx bucket summary")
                st.dataframe(
                    bucket_summary,
                    use_container_width=True,
                    hide_index=True,
                )

                st.subheader("3. Interpreted P/LP Geneyx candidates")
                st.caption(
                    "Relevance is filled and Pathogenic is Pathogenic/Likely Pathogenic. "
                    "VUS and benign rows are excluded."
                )
                st.dataframe(
                    interpreted_display,
                    use_container_width=True,
                    hide_index=True,
                )

                if not unmapped_candidates.empty:
                    unmapped_genes = sorted(
                        set(
                            unmapped_candidates["Gene"]
                            .astype(str)
                            .str.upper()
                            .str.strip()
                        )
                    )
                    st.info(
                        "P/LP candidates outside the active disease map are retained for trace "
                        f"but not assigned a risk tier: {', '.join(unmapped_genes)}"
                    )

                st.subheader("4. P/LP candidates covered by active disease mapping")
                st.dataframe(
                    mapped_display,
                    use_container_width=True,
                    hide_index=True,
                )

                st.subheader("5. Phenotype relevance mapping")
                st.dataframe(
                    relevance_table,
                    use_container_width=True,
                    hide_index=True,
                )

                st.subheader("6. Risk strata preview")
                st.dataframe(
                    risk_preview,
                    use_container_width=True,
                    hide_index=True,
                )

                st.download_button(
                    "Download risk strata preview",
                    data=make_csv_download(risk_preview),
                    file_name=f"{sample_id}_risk_strata_preview.csv",
                    mime="text/csv",
                    type="primary",
                )
