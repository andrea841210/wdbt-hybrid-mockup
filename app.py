from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import pandas as pd
import streamlit as st

APP_DIR = Path(__file__).resolve().parent
DATA_DIR = APP_DIR / "data"
ASSET_DIR = APP_DIR / "assets"

PAN_MASTER_PATH = DATA_DIR / "pan_syndrome_master.csv"
PAN_TO_HPO_PATH = DATA_DIR / "pan_syndrome_to_hpo.csv"
FUNCTION2_STATIC_PATH = DATA_DIR / "function2_risk_strata_static.csv"
REPORTING_PACKAGE_PATH = DATA_DIR / "reporting_package_mock.json"
REFERENCE_PDF_PATH = ASSET_DIR / "WDBT_UI_all-2.pdf"


st.set_page_config(
    page_title="WDBT Hybrid Workflow Mockup v0.2",
    page_icon="🧬",
    layout="wide",
)


@st.cache_data
def load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        st.error(f"Missing data file: {path.name}")
        return pd.DataFrame()
    return pd.read_csv(path)


@st.cache_data
def load_json(path: Path) -> dict:
    if not path.exists():
        st.error(f"Missing data file: {path.name}")
        return {}
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def ordered_unique(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys([str(v) for v in values]))


def build_geneyx_input_package(
    sample_id: str,
    module: str,
    selected_pan_ids: list[str],
    physician_note: str,
    hpo_df: pd.DataFrame,
) -> pd.DataFrame:
    package = hpo_df[["HPO_ID", "HPO_term"]].drop_duplicates().copy()
    package.insert(0, "Sample_ID", sample_id)
    package.insert(1, "Module", module)
    package.insert(2, "Selected_pan_syndrome_IDs", ", ".join(selected_pan_ids))
    package["Physician_note"] = physician_note
    package["Input_status"] = "Geneyx-ready mock input"
    return package


def filter_function1_outputs(
    module: str,
    selected_pan_ids: list[str],
    top_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pan_master = load_csv(PAN_MASTER_PATH)
    pan_to_hpo = load_csv(PAN_TO_HPO_PATH)

    pan_list = (
        pan_master[
            (pan_master["module"] == module)
            & (pan_master["pan_syndrome_id"].isin(selected_pan_ids))
        ]
        .sort_values("display_order")
        .reset_index(drop=True)
    )

    hpo_preview = (
        pan_to_hpo[pan_to_hpo["pan_syndrome_id"].isin(selected_pan_ids)]
        .sort_values(["pan_syndrome_id", "HPO_ID"])
        .head(top_n)
        .reset_index(drop=True)
    )

    return pan_list, hpo_preview


def render_pan_syndrome_list(df: pd.DataFrame, title: str) -> None:
    st.subheader(title)
    if df.empty:
        st.info("No pan-syndrome records found for current selection.")
        return
    for _, row in df.iterrows():
        st.markdown(
            f"**{row['pan_syndrome_id']}**　{row['pan_syndrome_label_zh']}"
        )


def render_scope_statement() -> None:
    st.info(
        "This mockup is not a production system, not a clinical decision support tool, "
        "not connected to Geneyx or TGIA base, and does not generate formal reports.\n\n"
        "本工具不是正式系統，不是臨床決策工具，未連接 Geneyx 或 TGIA base，也不產生正式報告。"
    )


def function1_page() -> None:
    st.title("WDBT Function 1｜兒基安泛性症狀轉譯層")
    st.caption("Pan-syndrome → HPO / Geneyx-ready input")

    pan_master = load_csv(PAN_MASTER_PATH)
    available_modules = ordered_unique(pan_master["module"].tolist()) or ["Neuro"]

    left, right = st.columns([0.32, 0.68], gap="large")

    with left:
        st.markdown("### Input")
        module = st.selectbox("Module", options=available_modules, index=0)
        module_pan = pan_master[pan_master["module"] == module]
        pan_options = ordered_unique(module_pan["pan_syndrome_id"].tolist())
        default_pan = [p for p in ["PS001", "PS003", "PS005"] if p in pan_options]
        selected_pan_ids = st.multiselect(
            "Pan_syndrome ID",
            options=pan_options,
            default=default_pan or pan_options[:3],
        )
        physician_note = st.text_input("醫師補充描述", value="Headache")
        st.text_input("HPO database_version", value="PGSafe DB_v1", disabled=True)
        top_n = st.number_input("Top N to show_HPO", min_value=1, max_value=500, value=200, step=10)
        sample_id = st.text_input("Sample ID", value="TS260508003")
        run_translation = st.button("Run Translation", type="primary", use_container_width=True)

    with right:
        if not run_translation:
            st.warning("Select pan-syndrome IDs and click **Run Translation** to generate the preview.")
            render_scope_statement()
            return

        if not selected_pan_ids:
            st.error("Please select at least one pan-syndrome ID.")
            return

        pan_list, hpo_preview = filter_function1_outputs(
            module=module,
            selected_pan_ids=selected_pan_ids,
            top_n=int(top_n),
        )
        geneyx_package = build_geneyx_input_package(
            sample_id=sample_id,
            module=module,
            selected_pan_ids=selected_pan_ids,
            physician_note=physician_note,
            hpo_df=hpo_preview,
        )

        render_pan_syndrome_list(pan_list, "Pan_syndrome list")

        st.subheader(f"HPO mapping preview (n={len(hpo_preview)})")
        st.dataframe(
            hpo_preview[["pan_syndrome_id", "HPO_ID", "HPO_term", "source_note"]],
            use_container_width=True,
            hide_index=True,
        )


        st.subheader("Geneyx-ready HPO input package")
        st.dataframe(geneyx_package, use_container_width=True, hide_index=True)
        st.download_button(
            "Download Geneyx-ready HPO input package",
            data=geneyx_package.to_csv(index=False).encode("utf-8-sig"),
            file_name=f"{sample_id}_geneyx_ready_hpo_input.csv",
            mime="text/csv",
            use_container_width=True,
        )

        render_scope_statement()


def function2_page() -> None:
    st.title("WDBT Function 2｜分子證據與風險分層預覽")
    st.caption("Geneyx evidence + original pan-syndrome context → combined risk strata")

    st.warning(
        "Demo preview only. Risk strata logic pending PM / bioinformatics confirmation.\n\n"
        "本頁為示意預覽；正式風險分層邏輯待 PM / 生資確認後導入。"
    )

    pan_master = load_csv(PAN_MASTER_PATH)
    risk_static = load_csv(FUNCTION2_STATIC_PATH)
    reporting_package = load_json(REPORTING_PACKAGE_PATH)
    pan_options = ordered_unique(pan_master["pan_syndrome_id"].tolist())
    default_pan = [p for p in ["PS001", "PS003", "PS005"] if p in pan_options]

    left, right = st.columns([0.32, 0.68], gap="large")

    with left:
        st.markdown("### Input")
        uploaded = st.file_uploader(
            "Geneyx output file",
            type=["csv", "tsv", "xlsx"],
            help="Optional only. This v0.2 page does not calculate from uploaded files.",
        )
        if uploaded is None:
            st.caption("Demo file: hpo_variation_merged_demo.csv")
        selected_pan_ids = st.multiselect(
            "Pan_syndrome ID",
            options=pan_options,
            default=default_pan or pan_options[:3],
        )
        physician_note = st.text_input("醫師補充描述", value="Headache")
        st.text_input("Risk logic version", value="Geneyx_v1", disabled=True)
        sample_id = st.text_input("Sample ID", value="TS260508003")
        run_risk = st.button("Run Risk Strata", type="primary", use_container_width=True)

    with right:
        if not run_risk:
            st.warning("Click **Run Risk Strata** to show the static completed-flow preview.")
            return

        pan_list = (
            pan_master[pan_master["pan_syndrome_id"].isin(selected_pan_ids)]
            .sort_values("display_order")
            .reset_index(drop=True)
        )
        render_pan_syndrome_list(pan_list, "Pan_syndrome list (recall)")

        st.subheader("Risk strata ranking")
        st.dataframe(risk_static, use_container_width=True, hide_index=True)

        st.subheader("Reporting package output")
        reporting_package = {
            **reporting_package,
            "sample_id": sample_id,
            "pan_syndrome_recall": selected_pan_ids,
            "physician_note": physician_note,
        }
        st.json(reporting_package, expanded=True)
        st.download_button(
            "Download reporting package output",
            data=json.dumps(reporting_package, ensure_ascii=False, indent=2).encode("utf-8"),
            file_name=f"{sample_id}_reporting_package_static_preview.json",
            mime="application/json",
            use_container_width=True,
        )

        render_scope_statement()


def main() -> None:
    st.sidebar.title("WDBT Hybrid Workflow Mockup v0.2")
    mode = st.sidebar.radio(
        "Mode",
        options=[
            "Function 1｜Translation demo",
            "Function 2｜Static completed-flow preview",
        ],
    )
    st.sidebar.markdown("---")
    st.sidebar.markdown("**Reference layout**")
    if REFERENCE_PDF_PATH.exists():
        st.sidebar.download_button(
            "Download WDBT_UI_all-2.pdf",
            data=REFERENCE_PDF_PATH.read_bytes(),
            file_name="WDBT_UI_all-2.pdf",
            mime="application/pdf",
            use_container_width=True,
        )
    st.sidebar.markdown(
        "Function 1 performs only one-layer translation. "
        "Function 2 is static until PM / bioinformatics confirms risk strata logic."
    )

    if mode.startswith("Function 1"):
        function1_page()
    else:
        function2_page()


if __name__ == "__main__":
    main()
