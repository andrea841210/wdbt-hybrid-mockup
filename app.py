
import json
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
# Sidebar
# -----------------------------
st.sidebar.title("WDBT Mockup v0.3")
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
    "Function 1: pan-syndrome → phenotype category → HPO / Geneyx-ready input. "
    "Function 2: static preview only."
)


# -----------------------------
# Function 1
# -----------------------------
if mode == "Function 1｜HPO Translation":
    st.markdown("## Function 1｜兒基安泛性症狀轉譯層")
    st.info(
        "This page performs one-layer input preparation only: "
        "pan-syndrome → phenotype category → HPO list. "
        "It does not perform interpretation, disease inference, or risk stratification."
    )

    pan_df = load_pan_syndrome_master()
    hpo_df = load_category_to_hpo()

    left, right = st.columns([1, 2])

    with left:
        st.subheader("Input")

        module = st.selectbox("Module", ["Neuro"], index=0)

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
            )

        else:
            st.warning("Select pan-syndrome IDs and click Run Translation.")


# -----------------------------
# Function 2
# -----------------------------
else:
    st.header("Function 2｜分子證據與風險分層預覽")
    st.warning(
        "Preview only. Formal risk-stratification logic will be added after "
        "PM / bioinformatics confirmation.\n\n"
        "本頁僅為示意預覽；正式風險分層邏輯待 PM / 生資確認後導入。"
    )

    pan_df = load_pan_syndrome_master()
    risk_df = load_function2_static()
    reporting_mock = load_reporting_package_mock()

    left, right = st.columns([1, 2])

    with left:
        st.subheader("Input")

        uploaded_file = st.file_uploader(
            "Geneyx output file",
            type=["csv", "xlsx", "xls"],
            help="Static preview only. Uploaded file is not processed in v0.3.",
        )

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
        show_preview = st.button("Show Static Preview", type="primary")

    with right:
        if show_preview:
            selected_pan = pan_df[pan_df["pan_syndrome_id"].isin(selected_ids)].copy()

            st.subheader("1. Pan-syndrome list recall")
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

            st.subheader("2. Risk strata ranking")
            st.dataframe(
                risk_df,
                use_container_width=True,
                hide_index=True,
            )

            st.subheader("3. Reporting package output")
            st.json(reporting_mock)

            st.download_button(
                "Download reporting package mock JSON",
                data=json.dumps(reporting_mock, ensure_ascii=False, indent=2).encode("utf-8"),
                file_name=f"{sample_id}_reporting_package_mock.json",
                mime="application/json",
            )
        else:
            st.warning("Click Show Static Preview to display the static completed-flow preview.")
