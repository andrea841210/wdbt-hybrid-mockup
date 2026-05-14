
import csv
import json
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
@st.cache_data
def load_pan_syndrome_master() -> pd.DataFrame:
    path = DATA_DIR / "pan_syndrome_master.csv"
    df = pd.read_csv(path, sep=None, engine="python")
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
    df = pd.read_csv(path, sep=None, engine="python")
    expected_cols = {
        "disease_id",
        "gene",
        "disease_name",
        "phenotype_concern_list",
        "inheritance",
        "severity",
        "referral",
        "possible_tiers",
    }
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(f"gene_disease_phenotype_map.csv missing columns: {missing}")
    df["gene"] = df["gene"].astype(str).str.upper().str.strip()
    df["phenotype_concern_list"] = df["phenotype_concern_list"].fillna("")
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
MVP_SCOPE_GENES = ["DMD", "FOLR1", "GJB2", "GRIN2B", "SCN1A"]
RELEVANCE_HIGH_BUCKETS = [
    "Dominant HET",
    "Recessive Compound HET",
    "Recessive HOM",
    "Mitochondria",
]


def decode_uploaded_tsv(uploaded_file) -> str:
    """
    Decode uploaded TSV bytes with common encodings.
    """
    raw = uploaded_file.getvalue()
    for encoding in ["utf-8-sig", "utf-8", "big5", "latin1"]:
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


def parse_all_geneyx_files(uploaded_files: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Parse all required Geneyx TSV files.

    Returns:
    bucket_summary
    all_candidates
    mvp_in_scope_candidates
    """
    all_tables = []
    summary_rows = []

    for bucket, uploaded_file in uploaded_files.items():
        metadata, df = parse_geneyx_tsv(uploaded_file, bucket)

        variant_count = len(df) if not df.empty else 0
        summary_rows.append(
            {
                "inheritance_bucket": bucket,
                "sample_id": metadata.get("sample_id", ""),
                "analysis_id": metadata.get("analysis_id", ""),
                "patient_id": metadata.get("patient_id", ""),
                "variant_count": variant_count,
                "status": "has variants" if variant_count > 0 else "no call-out",
            }
        )

        if variant_count > 0:
            all_tables.append(df)

    bucket_summary = pd.DataFrame(summary_rows)

    if not all_tables:
        all_candidates = pd.DataFrame()
        mvp_candidates = pd.DataFrame()
        return bucket_summary, all_candidates, mvp_candidates

    all_candidates = pd.concat(all_tables, ignore_index=True)

    if "Gene" not in all_candidates.columns:
        mvp_candidates = pd.DataFrame()
    else:
        mvp_candidates = all_candidates[
            all_candidates["Gene"].astype(str).str.upper().isin(MVP_SCOPE_GENES)
        ].copy()

    return bucket_summary, all_candidates, mvp_candidates


def simplify_geneyx_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Reduce Geneyx candidate table to Function 2 display columns.
    Missing columns are kept as empty strings.
    """
    display_cols = [
        "sample_id",
        "analysis_id",
        "patient_id",
        "inheritance_bucket",
        "Relevance",
        "Pathogenic",
        "Gene",
        "Location",
        "Omim",
        "Omim Inheritance",
        "AA",
        "HGVSC",
        "HGVSP",
        "Zygosity",
        "RefSeq",
        "dbSNP",
        "ACMG",
        "ClinVar",
        "Effect",
        "Severity",
        "Max AF",
        "CADD Phred",
        "SpliceAi Score",
        "Phenotype",
        "Matched Phenotypes",
        "Filter",
    ]

    if df.empty:
        return pd.DataFrame(columns=display_cols)

    out = df.copy()
    for col in display_cols:
        if col not in out.columns:
            out[col] = ""

    return out[display_cols]


def build_phenotype_relevance_and_risk_preview(
    mvp_candidates: pd.DataFrame,
    mapping_df: pd.DataFrame,
    selected_categories: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Connect Geneyx MVP candidates to PM gene-disease-phenotype mapping,
    then apply MVP risk tier rules.
    """
    if mvp_candidates.empty:
        columns = [
            "risk_level",
            "color",
            "risk_tier",
            "concern",
            "phenotype_relevance",
            "matched_phenotype_categories",
            "gene",
            "variant",
            "disease_id",
            "disease_name",
            "inheritance_bucket",
            "relevance",
            "pathogenic",
            "acmg",
            "clinvar",
            "severity",
            "referral",
            "rule_trace",
        ]
        return pd.DataFrame(columns=columns), pd.DataFrame(columns=columns)

    selected_set = set(selected_categories)
    rows = []

    for _, variant in mvp_candidates.iterrows():
        gene = str(variant.get("Gene", "")).upper().strip()
        gene_maps = mapping_df[mapping_df["gene"] == gene].copy()

        # If no PM mapping exists, keep the candidate but mark as out-of-map.
        if gene_maps.empty:
            gene_maps = pd.DataFrame(
                [
                    {
                        "disease_id": "",
                        "gene": gene,
                        "disease_name": "",
                        "phenotype_concern_list": "",
                        "inheritance": "",
                        "severity": "",
                        "referral": "",
                        "possible_tiers": "",
                    }
                ]
            )

        for _, mapping in gene_maps.iterrows():
            phenotype_concerns = split_semicolon_list(mapping.get("phenotype_concern_list", ""))
            matched_categories = sorted(selected_set.intersection(set(phenotype_concerns)))
            phenotype_relevance = "Yes" if matched_categories else "No"

            inheritance_bucket = str(variant.get("inheritance_bucket", ""))
            relevance = str(variant.get("Relevance", "")).strip()
            relevance_high = relevance.lower() == "high"

            # MVP v0.4 rule calibration:
            # Geneyx "Relevance" is not consistently populated across buckets in the test files.
            # Therefore, phenotype relevance is used as the primary MVP trigger,
            # while Relevance is retained as a displayed evidence field rather than a hard gate.
            if inheritance_bucket in RELEVANCE_HIGH_BUCKETS:
                if phenotype_relevance == "Yes":
                    risk_level = 1
                    color = "🔴 Red"
                    risk_tier = "High Attention"
                    rule_trace = "Compatible bucket + phenotype relevance Yes"
                else:
                    risk_level = 2
                    color = "🟠 Orange"
                    risk_tier = "Targeted Follow-up"
                    rule_trace = "Compatible bucket + phenotype relevance No"
            elif inheritance_bucket == "Recessive HET":
                if phenotype_relevance == "Yes":
                    risk_level = 3
                    color = "🟡 Yellow"
                    risk_tier = "Routine Monitoring"
                    rule_trace = "Recessive HET + phenotype relevance Yes"
                else:
                    risk_level = 4
                    color = "🔵 Blue"
                    risk_tier = "General Awareness"
                    rule_trace = "Recessive HET + phenotype relevance No"
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
                    "gene": gene,
                    "variant": variant_label,
                    "disease_id": mapping.get("disease_id", ""),
                    "disease_name": mapping.get("disease_name", ""),
                    "inheritance_bucket": inheritance_bucket,
                    "relevance": relevance,
                    "pathogenic": variant.get("Pathogenic", ""),
                    "acmg": variant.get("ACMG", ""),
                    "clinvar": variant.get("ClinVar", ""),
                    "severity": mapping.get("severity", ""),
                    "referral": mapping.get("referral", ""),
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
st.sidebar.title("WDBT Mockup v0.4")
mode = st.sidebar.radio(
    "Mode",
    ["Function 1｜HPO Translation", "Function 2｜Risk Strata Preview"],
)

st.sidebar.caption(
    "Hybrid mockup only. Not a production system, not a clinical decision support tool."
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
        "MVP scope: DMD, FOLR1, GJB2, GRIN2B, SCN1A only. "
        "This page uses uploaded Geneyx TSV files and PM-defined phenotype mapping "
        "to generate a rule-based risk-strata preview."
    )
    st.caption(
        "本頁僅示意 MVP 範圍內之基因：DMD、FOLR1、GJB2、GRIN2B、SCN1A。"
        "系統依上傳之 Geneyx TSV 與 PM 設定之 phenotype mapping 產生風險分層預覽；"
        "本頁不代表正式解讀或臨床決策報告。"
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

        run_preview = st.button("Run MVP Risk Preview", type="primary")

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
                selected_categories = []

                for _, row in selected_pan.iterrows():
                    selected_categories.extend(
                        split_categories(
                            row["primary_phenotype_category"],
                            row["secondary_phenotype_category"],
                        )
                    )

                selected_categories = sorted(set(selected_categories))

                bucket_summary, all_candidates, mvp_candidates = parse_all_geneyx_files(uploaded_files)
                mvp_display = simplify_geneyx_candidates(mvp_candidates)
                relevance_table, risk_preview = build_phenotype_relevance_and_risk_preview(
                    mvp_candidates=mvp_candidates,
                    mapping_df=mapping_df,
                    selected_categories=selected_categories,
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

                st.caption("Selected phenotype categories")
                st.write(", ".join(selected_categories) if selected_categories else "None")

                st.subheader("2. Geneyx bucket summary")
                st.dataframe(
                    bucket_summary,
                    use_container_width=True,
                    hide_index=True,
                )

                st.subheader("3. MVP in-scope candidate variants")
                st.caption("MVP scope genes: DMD, FOLR1, GJB2, GRIN2B, SCN1A.")
                st.dataframe(
                    mvp_display,
                    use_container_width=True,
                    hide_index=True,
                )

                st.subheader("4. Phenotype relevance mapping")
                st.dataframe(
                    relevance_table,
                    use_container_width=True,
                    hide_index=True,
                )

                st.subheader("5. Risk strata preview")
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
        else:
            st.caption("Upload all five Geneyx TSV files, select pan-syndrome context, then run MVP risk preview.")
