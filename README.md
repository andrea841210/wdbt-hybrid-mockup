# WDBT Hybrid Workflow Mockup v0.2

This Streamlit app is a hybrid WDBT mockup.
Function 1 demonstrates a one-layer pan-syndrome to HPO translation workflow.
Function 2 is a static preview of the post-Geneyx risk strata screen.

This is not a production system, not a clinical decision support tool,
not connected to Geneyx or TGIA base, and does not generate formal reports.

---

## 中文說明

本 Streamlit app 為 WDBT hybrid mockup。  
Function 1 示意泛性症狀至 HPO / Geneyx-ready input 的單層轉譯流程。  
Function 2 僅為 Geneyx 後風險分層畫面的靜態預覽。

本工具不是正式系統，不是臨床決策工具，未連接 Geneyx 或 TGIA base，也不產生正式報告。

---

## v0.2 updates

- Function 1 removes the original disease scope preview table.
- Function 2 adds two visible fields to the static risk strata ranking table:
  - `Color`: visual risk color label. Level 1 = red, Level 2 = orange, Level 3 = yellow, Level 4 = blue.
  - `Concern`: simplified external-facing risk strata concern.
- Reference layout updated to `assets/WDBT_UI_all-2.pdf`.

---

## What this mockup does

### Function 1｜Translation demo

User selects:

- Module
- Pan_syndrome ID
- Physician supplemental note
- Top N HPO terms to show
- Sample ID

Then the app reads:

- `data/pan_syndrome_master.csv`
- `data/pan_syndrome_to_hpo.csv`

and generates:

1. Pan_syndrome list
2. HPO mapping preview
3. Geneyx-ready HPO input package download

Logic boundary:

> selected pan-syndrome ID → linked HPO terms / IDs → Geneyx-ready input package

It does not perform interpretation, risk stratification, disease inference, or report generation.

### Function 2｜Static completed-flow preview

Function 2 displays a fixed preview only:

1. Pan_syndrome list recall
2. Risk strata ranking, including `Color` and `Concern`
3. Reporting package output

Risk strata logic is not implemented in v0.2.

Displayed disclaimer:

> Demo preview only. Risk strata logic pending PM / bioinformatics confirmation.  
> 本頁為示意預覽；正式風險分層邏輯待 PM / 生資確認後導入。

---

## Repository structure

```text
wdbt-hybrid-mockup/
├── app.py
├── data/
│   ├── pan_syndrome_master.csv
│   ├── pan_syndrome_to_hpo.csv
│   ├── function2_risk_strata_static.csv
│   └── reporting_package_mock.json
├── assets/
│   └── WDBT_UI_all-2.pdf
├── requirements.txt
└── README.md
```

---

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

---

## Deployment note

For Streamlit Community Cloud:

1. Push this folder to GitHub.
2. Create a new Streamlit app from the GitHub repository.
3. Set the main file path to `app.py`.
4. Deploy.
