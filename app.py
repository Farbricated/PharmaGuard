"""
PharmaGuard v2.0 — Pharmacogenomic Risk Prediction System
RIFT 2026 Hackathon | Pharmacogenomics / Explainable AI Track
"""

import streamlit as st
import json, uuid, os
from datetime import datetime
from dotenv import load_dotenv

load_dotenv()

from vcf_parser import parse_vcf, get_sample_vcf
from risk_engine import run_risk_assessment, get_overall_severity, DRUG_RISK_TABLE
from llm_explainer import generate_all_explanations
from schema import build_output_schema
from drug_interactions import run_interaction_analysis
from pdf_report import generate_pdf_report

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

SAMPLE_VCF_MAP = {
    "Mixed Variants (Standard)":              "sample.vcf",
    "UltraRapid Metabolizer - Codeine TOXIC": "test_ultrarapid_metabolizer.vcf",
    "All Normal Wild-type - All Safe":         "test_all_normal_wildtype.vcf",
    "Worst Case - All Poor Metabolizers":      "test_worst_case_all_pm.vcf",
}

ALL_DRUGS = list(DRUG_RISK_TABLE.keys())

SEV_RANK = {"none":0,"low":1,"moderate":2,"high":3,"critical":4}
SEV_COLOR = {"none":"#10b981","low":"#f59e0b","moderate":"#f97316","high":"#ef4444","critical":"#991b1b"}
RISK_BADGE = {"Safe":"risk-badge-safe","Adjust Dosage":"risk-badge-adjust",
              "Toxic":"risk-badge-toxic","Ineffective":"risk-badge-ineffective","Unknown":"risk-badge-unknown"}
RISK_EMOJI = {"Safe":"🟢","Adjust Dosage":"🟡","Toxic":"🔴","Ineffective":"🟣","Unknown":"⚪"}

# ─── Page Config ──────────────────────────────────────────────
st.set_page_config(page_title="PharmaGuard v2",page_icon="🧬",layout="wide",initial_sidebar_state="expanded")

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');
html,body,[class*="css"]{font-family:'Space Grotesk',sans-serif;}
.main-header{background:linear-gradient(135deg,#0f172a 0%,#1e293b 50%,#0f4c75 100%);padding:1.8rem 2.5rem;border-radius:16px;margin-bottom:1.5rem;border:1px solid #1e40af44;}
.main-header h1{color:#e0f2fe;font-size:2.1rem;font-weight:700;margin:0;}
.main-header p{color:#94a3b8;margin:0.3rem 0 0 0;font-size:0.9rem;}
.v2-badge{display:inline-block;background:linear-gradient(135deg,#7c3aed,#4f46e5);color:#fff;font-size:0.62rem;font-weight:700;padding:2px 9px;border-radius:999px;margin-left:8px;vertical-align:middle;letter-spacing:1px;}
.risk-badge-safe{background:#065f46;color:#6ee7b7;padding:4px 14px;border-radius:999px;font-weight:700;font-size:0.9rem;border:1.5px solid #10b981;display:inline-block;}
.risk-badge-adjust{background:#78350f;color:#fde68a;padding:4px 14px;border-radius:999px;font-weight:700;font-size:0.9rem;border:1.5px solid #f59e0b;display:inline-block;}
.risk-badge-toxic{background:#7f1d1d;color:#fca5a5;padding:4px 14px;border-radius:999px;font-weight:700;font-size:0.9rem;border:1.5px solid #ef4444;display:inline-block;}
.risk-badge-ineffective{background:#312e81;color:#c4b5fd;padding:4px 14px;border-radius:999px;font-weight:700;font-size:0.9rem;border:1.5px solid #8b5cf6;display:inline-block;}
.risk-badge-unknown{background:#1f2937;color:#9ca3af;padding:4px 14px;border-radius:999px;font-weight:700;font-size:0.9rem;border:1.5px solid #6b7280;display:inline-block;}
.metric-card{background:#1e293b;border:1px solid #334155;border-radius:10px;padding:1rem 1.2rem;margin:0.4rem 0;}
.metric-card h4{color:#64748b;font-size:0.72rem;text-transform:uppercase;letter-spacing:1px;margin:0 0 0.3rem 0;}
.metric-card p{color:#e2e8f0;font-size:1rem;font-weight:600;margin:0;font-family:'JetBrains Mono',monospace;}
.explanation-box{background:#0f172a;border:1px solid #1e40af44;border-left:4px solid #3b82f6;border-radius:8px;padding:0.9rem 1.3rem;margin:0.6rem 0;}
.explanation-box h5{color:#60a5fa;font-size:0.73rem;text-transform:uppercase;letter-spacing:1px;margin:0 0 0.4rem 0;}
.explanation-box p{color:#cbd5e1;line-height:1.7;margin:0;}
.warning-box{background:#451a03;border:1px solid #92400e;border-left:4px solid #f97316;border-radius:8px;padding:0.9rem 1.3rem;margin:0.6rem 0;}
.critical-box{background:#450a0a;border:1px solid #991b1b;border-left:4px solid #ef4444;border-radius:8px;padding:0.9rem 1.3rem;margin:0.6rem 0;}
.interaction-card{background:#1e1b4b;border:1px solid #4338ca;border-left:4px solid #818cf8;border-radius:8px;padding:0.9rem 1.3rem;margin:0.5rem 0;}
.interaction-card h5{color:#a5b4fc;font-size:0.73rem;text-transform:uppercase;letter-spacing:1px;margin:0 0 0.35rem 0;}
.interaction-card p{color:#c7d2fe;line-height:1.6;margin:0;}
.interaction-critical{border-left-color:#ef4444!important;background:#450a0a!important;}
.interaction-critical h5{color:#fca5a5!important;}.interaction-critical p{color:#fee2e2!important;}
.interaction-high{border-left-color:#f97316!important;background:#431407!important;}
.interaction-high h5{color:#fdba74!important;}.interaction-high p{color:#fed7aa!important;}
.summary-dashboard{background:linear-gradient(135deg,#0f172a,#1e293b);border:1px solid #334155;border-radius:12px;padding:1.2rem 1.5rem;margin:0.8rem 0;}
.test-card{background:#0f172a;border:1px solid #1e40af55;border-radius:10px;padding:1rem;margin:0.4rem 0;}
.test-card h4{color:#60a5fa;margin:0 0 0.3rem 0;font-size:0.88rem;}
.test-card p{color:#94a3b8;font-size:0.78rem;margin:0;}
.stButton>button{background:linear-gradient(135deg,#1d4ed8,#2563eb);color:white;border:none;border-radius:8px;padding:0.52rem 1.2rem;font-weight:600;font-family:'Space Grotesk',sans-serif;}
div[data-testid="stSidebar"]{background:#0f172a;}
</style>
""", unsafe_allow_html=True)

# ─── Sidebar ──────────────────────────────────────────────────
with st.sidebar:
    st.markdown("### ⚙️ Configuration")
    groq_api_key = st.text_input("Groq API Key", value=os.environ.get("GROQ_API_KEY",""), type="password")
    st.markdown("---")
    st.markdown("### 🧬 Supported Genes")
    for g in ["CYP2D6","CYP2C19","CYP2C9","SLCO1B1","TPMT","DPYD"]: st.markdown(f"• `{g}`")
    st.markdown("---")
    st.markdown("### 💊 Supported Drugs")
    for d in ALL_DRUGS: st.markdown(f"• {d.title()}")
    st.markdown("---")
    st.markdown("""**PharmaGuard v2.0** <span class="v2-badge">V2</span><br>
    RIFT 2026 · [CPIC Guidelines](https://cpicpgx.org)""", unsafe_allow_html=True)

# ─── Header ───────────────────────────────────────────────────
st.markdown("""
<div class="main-header">
  <h1>🧬 PharmaGuard <span class="v2-badge">V2</span></h1>
  <p>Pharmacogenomic Risk Prediction · Drug-Drug Interactions · PDF Reports · Full Test Suite · RIFT 2026 · Groq LLaMA 3.3</p>
</div>
""", unsafe_allow_html=True)

# ─── Helper Functions ─────────────────────────────────────────
def load_vcf_file(filename: str) -> str:
    path = os.path.join(BASE_DIR, "sample_data", filename)
    if os.path.exists(path):
        with open(path) as f:
            return f.read()
    return get_sample_vcf()


def safe_pdf(pid, all_outputs, parsed_vcf):
    try:
        return generate_pdf_report(pid, all_outputs, parsed_vcf)
    except Exception as e:
        return None, str(e)


def run_pipeline(vcf_content, drugs, pid, groq_key, run_ix=True, gen_pdf=True):
    parsed_vcf   = parse_vcf(vcf_content)
    risk_results = run_risk_assessment(parsed_vcf, drugs)
    if groq_key:
        risk_results = generate_all_explanations(groq_key, risk_results)
    else:
        for r in risk_results:
            r["llm_explanation"] = {"summary":"No API key.","biological_mechanism":"",
                                    "variant_significance":"","clinical_implications":"","success":False}
    all_outputs = [build_output_schema(patient_id=pid, drug=r["drug"], result=r,
                   parsed_vcf=parsed_vcf, llm_exp=r.get("llm_explanation",{})) for r in risk_results]
    ix_report = run_interaction_analysis(drugs, risk_results) if run_ix and len(drugs)>1 else None
    pdf_bytes = None
    if gen_pdf:
        try: pdf_bytes = generate_pdf_report(pid, all_outputs, parsed_vcf)
        except Exception as e: st.warning(f"PDF error: {e}")
    return parsed_vcf, risk_results, all_outputs, ix_report, pdf_bytes


def render_dashboard(parsed_vcf, all_outputs, ix_report, pdf_bytes, pid):
    overall_sev = max((o["risk_assessment"]["severity"] for o in all_outputs),
                      key=lambda s: SEV_RANK.get(s,0), default="none")
    sev_col = SEV_COLOR.get(overall_sev,"#6b7280")

    st.markdown(f"""
    <div class="summary-dashboard">
      <div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:0.8rem;">
        <div><div style="color:#64748b;font-size:0.7rem;text-transform:uppercase;letter-spacing:1px;">Patient ID</div>
             <div style="color:#e2e8f0;font-size:1.1rem;font-weight:700;font-family:'JetBrains Mono',monospace;">{pid}</div></div>
        <div><div style="color:#64748b;font-size:0.7rem;text-transform:uppercase;letter-spacing:1px;">Overall Risk</div>
             <div style="color:{sev_col};font-size:1.1rem;font-weight:700;">{overall_sev.upper()}</div></div>
        <div><div style="color:#64748b;font-size:0.7rem;text-transform:uppercase;letter-spacing:1px;">Drugs</div>
             <div style="color:#e2e8f0;font-size:1.1rem;font-weight:700;">{len(all_outputs)}</div></div>
        <div><div style="color:#64748b;font-size:0.7rem;text-transform:uppercase;letter-spacing:1px;">Variants</div>
             <div style="color:#e2e8f0;font-size:1.1rem;font-weight:700;">{parsed_vcf['total_variants']}</div></div>
        <div><div style="color:#64748b;font-size:0.7rem;text-transform:uppercase;letter-spacing:1px;">Genes</div>
             <div style="color:#e2e8f0;font-size:1.1rem;font-weight:700;">{len(parsed_vcf['detected_genes'])}/6</div></div>
      </div>
    </div>""", unsafe_allow_html=True)

    all_json = json.dumps(all_outputs, indent=2)
    dc1,dc2,dc3 = st.columns(3)
    with dc1:
        st.download_button("⬇️ Download All JSON", data=all_json,
                           file_name=f"pharmaguard_{pid}_all.json", mime="application/json",
                           use_container_width=True, key=f"dlall_{pid}")
    with dc2:
        if pdf_bytes:
            st.download_button("📄 Download PDF Report", data=pdf_bytes,
                               file_name=f"pharmaguard_{pid}_report.pdf", mime="application/pdf",
                               use_container_width=True, key=f"dlpdf_{pid}")
        else:
            st.button("📄 PDF unavailable", disabled=True, use_container_width=True, key=f"nopdf_{pid}")
    with dc3:
        st.markdown(f"""<button onclick="navigator.clipboard.writeText({json.dumps(all_json)}).then(()=>{{
            this.innerText='Copied!';setTimeout(()=>this.innerText='Copy All JSON',2000);}})
            " style="width:100%;padding:0.4rem 0.8rem;background:#1e293b;color:#e2e8f0;
            border:1px solid #334155;border-radius:8px;font-size:0.82rem;font-weight:600;
            cursor:pointer;">📋 Copy All JSON</button>""", unsafe_allow_html=True)

    with st.expander("📋 VCF Parsing Summary", expanded=False):
        pc1,pc2,pc3 = st.columns(3)
        pc1.metric("Total Variants", parsed_vcf["total_variants"])
        pc2.metric("Genes Found", len(parsed_vcf["detected_genes"]))
        pc3.metric("Parse Errors", len(parsed_vcf["parse_errors"]))
        if parsed_vcf["detected_genes"]:
            st.markdown(f"**Detected:** `{'`, `'.join(parsed_vcf['detected_genes'])}`")

    if ix_report and ix_report["interactions_found"]:
        ix_sev  = ix_report["overall_severity"]
        ix_icon = {"low":"🟡","moderate":"🟠","high":"🔴","critical":"🚨"}.get(ix_sev,"⚪")
        with st.expander(f"{ix_icon} Drug-Drug Interactions — {ix_report['total_interactions']} alert(s) [{ix_sev.upper()}]", expanded=True):
            for ix in ix_report["all_interactions"]:
                sev   = ix.get("severity","low")
                extra = f"interaction-{sev}" if sev in ("critical","high") else ""
                d_str = " + ".join(ix.get("drugs_involved",[]) or
                                   [ix.get("inhibitor_drug","")]+ix.get("affected_drugs",[]))
                st.markdown(f"""<div class="interaction-card {extra}">
                    <h5>⚠️ {ix.get('type','').replace('_',' ').title()} — {d_str} [{sev.upper()}]</h5>
                    <p>{ix.get('message',ix.get('mechanism',''))}</p>
                    <p style="margin-top:0.4rem;font-style:italic;">
                    💡 <strong>Rec:</strong> {ix.get('recommendation','')}</p>
                </div>""", unsafe_allow_html=True)
    elif ix_report:
        st.success("✅ No significant drug-drug interactions detected.")

    for output in all_outputs:
        risk_label = output["risk_assessment"]["risk_label"]
        drug_name  = output["drug"]
        severity   = output["risk_assessment"]["severity"].upper()
        confidence = output["risk_assessment"]["confidence_score"]
        emoji      = RISK_EMOJI.get(risk_label,"⚪")
        badge      = RISK_BADGE.get(risk_label,"risk-badge-unknown")

        with st.expander(f"{emoji} {drug_name.title()} — {risk_label} (Severity: {severity})", expanded=True):
            if severity in ("CRITICAL","HIGH"):
                st.markdown(f'<div class="critical-box"><strong>⚠️ CLINICAL ALERT:</strong> {output["clinical_recommendation"]["dosing_recommendation"]}</div>', unsafe_allow_html=True)
            elif severity == "MODERATE":
                st.markdown(f'<div class="warning-box"><strong>⚠️ NOTE:</strong> {output["clinical_recommendation"]["dosing_recommendation"]}</div>', unsafe_allow_html=True)

            m1,m2,m3,m4 = st.columns(4)
            for col, hdr, val in [(m1,"Risk Label",f'<span class="{badge}">{risk_label}</span>'),
                                   (m2,"Primary Gene",output["pharmacogenomic_profile"]["primary_gene"]),
                                   (m3,"Diplotype",output["pharmacogenomic_profile"]["diplotype"]),
                                   (m4,"Phenotype",output["pharmacogenomic_profile"]["phenotype"])]:
                col.markdown(f'<div class="metric-card"><h4>{hdr}</h4><p>{val}</p></div>', unsafe_allow_html=True)

            st.markdown(f"**Confidence:** `{confidence:.0%}`")
            st.progress(confidence)

            variants = output["pharmacogenomic_profile"]["detected_variants"]
            if variants:
                st.markdown("**🔬 Detected Variants:**")
                vc = st.columns([1.2,1,1,1,2])
                for h,c in zip(["rsID","Gene","Star Allele","REF>ALT","Function"], vc): c.markdown(f"**{h}**")
                for v in variants:
                    vc[0].markdown(f"`{v['rsid']}`"); vc[1].markdown(f"`{v.get('gene','N/A')}`")
                    vc[2].markdown(f"`{v.get('star_allele','N/A')}`")
                    vc[3].markdown(f"`{v.get('ref','?')}>{v.get('alt','?')}`")
                    vc[4].markdown(v.get("functional_status","Unknown"))
            else:
                st.info("No variants detected. Wild-type (*1/*1) assumed.")

            st.markdown("---")
            st.markdown("**📋 CPIC Recommendation**")
            st.markdown(f"> {output['clinical_recommendation']['dosing_recommendation']}")
            if output["clinical_recommendation"]["alternative_drugs"]:
                st.markdown(f"**Alternatives:** {', '.join(output['clinical_recommendation']['alternative_drugs'])}")
            if output["clinical_recommendation"]["monitoring_required"]:
                st.markdown(f"**Monitoring:** {output['clinical_recommendation']['monitoring_required']}")

            exp = output["llm_generated_explanation"]
            if exp.get("summary"):
                st.markdown("---\n**🤖 AI Clinical Explanation (Groq LLaMA 3.3)**")
                for label, key, icon in [("Summary","summary","📌"),
                                          ("Biological Mechanism","biological_mechanism","🔬"),
                                          ("Variant Significance","variant_significance","🧬"),
                                          ("Clinical Implications","clinical_implications","🏥")]:
                    if exp.get(key):
                        st.markdown(f'<div class="explanation-box"><h5>{icon} {label}</h5><p>{exp[key]}</p></div>', unsafe_allow_html=True)

            st.markdown("---")
            json_str = json.dumps(output, indent=2)
            bc1,bc2 = st.columns(2)
            with bc1:
                st.download_button("⬇️ Download JSON", data=json_str,
                                   file_name=f"pharmaguard_{pid}_{drug_name}.json",
                                   mime="application/json", key=f"dl_{pid}_{drug_name}",
                                   use_container_width=True)
            with bc2:
                st.markdown(f"""<button onclick="navigator.clipboard.writeText({json.dumps(json_str)}).then(()=>{{
                    this.innerText='Copied!';setTimeout(()=>this.innerText='Copy JSON',2000);}})
                    " style="width:100%;padding:0.4rem 0.8rem;background:#1e293b;color:#e2e8f0;
                    border:1px solid #334155;border-radius:8px;font-size:0.82rem;font-weight:600;cursor:pointer;">
                    📋 Copy JSON</button>""", unsafe_allow_html=True)
            with st.expander("View Raw JSON", expanded=False):
                st.code(json_str, language="json")


# ══════════════════════════════════════════════════════════════
# TABS
# ══════════════════════════════════════════════════════════════
tab1, tab2 = st.tabs(["🔬 Analysis", "🧪 Full Test Suite"])


# ── TAB 1 ─────────────────────────────────────────────────────
with tab1:
    col1, col2 = st.columns([1.2,1])
    with col1:
        st.markdown("#### 📁 Upload VCF File")
        uploaded_file = st.file_uploader("Drag and drop or browse", type=["vcf"])
        if uploaded_file is not None:
            sz = uploaded_file.size / (1024*1024)
            if sz > 5:
                st.error(f"❌ File too large: {sz:.1f}MB. Max 5MB.")
                uploaded_file = None
            else:
                peek = uploaded_file.read(500).decode("utf-8", errors="replace")
                uploaded_file.seek(0)
                if "##fileformat=VCF" not in peek and "#CHROM" not in peek:
                    st.error("❌ Invalid VCF format.")
                    uploaded_file = None
                else:
                    st.success(f"✅ `{uploaded_file.name}` ({sz:.2f} MB / 5 MB max)")

        st.markdown("**Or choose a test scenario:**")
        test_scenario = st.selectbox("Quick test VCF", ["None"] + list(SAMPLE_VCF_MAP.keys()))

    with col2:
        st.markdown("#### 💊 Drug Selection")
        drug_multiselect = st.multiselect("Select drugs", options=ALL_DRUGS,
                                           default=["CLOPIDOGREL"], format_func=lambda x: x.title())
        custom_drugs = st.text_input("Or type drugs (comma-separated)",
                                      placeholder="e.g., CODEINE, WARFARIN")
        patient_id_input = st.text_input("Patient ID (optional)", placeholder="Auto-generated")
        st.markdown("#### ⚙️ Options")
        run_interactions = st.checkbox("🔗 Drug-Drug Interaction Analysis", value=True)
        generate_pdf     = st.checkbox("📄 Generate PDF Report", value=True)

    st.markdown("---")
    analyze_btn = st.button("🔬 Run Pharmacogenomic Analysis", use_container_width=True)

    if analyze_btn:
        all_drugs = list(drug_multiselect)
        if custom_drugs.strip():
            all_drugs += [d.strip().upper() for d in custom_drugs.split(",") if d.strip()]
        all_drugs = list(set(all_drugs))
        if not all_drugs:
            st.error("⚠️ Select at least one drug.")
            st.stop()

        vcf_content = None
        if uploaded_file:
            vcf_content = uploaded_file.read().decode("utf-8", errors="replace")
        elif test_scenario != "None":
            vcf_content = load_vcf_file(SAMPLE_VCF_MAP[test_scenario])
            st.info(f"🧪 Test scenario: **{test_scenario}**")
        else:
            st.error("⚠️ Upload a VCF or select a test scenario.")
            st.stop()

        pid = patient_id_input.strip() or f"PATIENT_{str(uuid.uuid4())[:8].upper()}"
        st.markdown(f"## 📊 Results — `{pid}`")
        with st.spinner("🔬 Running full pharmacogenomic analysis..."):
            parsed_vcf, risk_results, all_outputs, ix_report, pdf_bytes = run_pipeline(
                vcf_content, all_drugs, pid, groq_api_key, run_interactions, generate_pdf)
        render_dashboard(parsed_vcf, all_outputs, ix_report, pdf_bytes, pid)

    else:
        st.markdown("""
        <div style="text-align:center;padding:2.5rem;color:#475569;">
          <div style="font-size:3.5rem;margin-bottom:0.8rem;">🧬</div>
          <h3 style="color:#94a3b8;">PharmaGuard v2.0 Ready</h3>
          <p>Upload VCF or pick a scenario · Select drugs · Click <strong>Run Analysis</strong></p>
          <p style="font-size:0.82rem;color:#64748b;">V2: Drug Interactions · PDF Reports · 4 Test Scenarios · Full Test Suite Tab</p>
        </div>""", unsafe_allow_html=True)


# ── TAB 2 — FULL TEST SUITE ───────────────────────────────────
TEST_SUITE = [
    {
        "name":     "Mixed Variants (Standard)",
        "file":     "sample.vcf",
        "drugs":    ["CLOPIDOGREL","CODEINE","AZATHIOPRINE"],
        "expected": {"CLOPIDOGREL":"Ineffective","CODEINE":"Adjust Dosage","AZATHIOPRINE":"Adjust Dosage"},
        "desc":     "Multi-gene profile — tests CYP2C19 + CYP2D6 + TPMT simultaneously",
        "icon":     "🧬",
    },
    {
        "name":     "UltraRapid Metabolizer - Codeine TOXIC",
        "file":     "test_ultrarapid_metabolizer.vcf",
        "drugs":    ["CODEINE","CLOPIDOGREL"],
        "expected": {"CODEINE":"Toxic","CLOPIDOGREL":"Safe"},
        "desc":     "CYP2D6 gene duplication — validates TOXIC risk for ultra-rapid metabolizers",
        "icon":     "🔴",
    },
    {
        "name":     "All Normal Wild-type - All Safe",
        "file":     "test_all_normal_wildtype.vcf",
        "drugs":    ALL_DRUGS,
        "expected": {d:"Safe" for d in ALL_DRUGS},
        "desc":     "Wild-type *1/*1 across all 6 genes — all drugs should return Safe",
        "icon":     "🟢",
    },
    {
        "name":     "Worst Case - All Poor Metabolizers",
        "file":     "test_worst_case_all_pm.vcf",
        "drugs":    ALL_DRUGS,
        "expected": {"CODEINE":"Ineffective","CLOPIDOGREL":"Ineffective","WARFARIN":"Adjust Dosage",
                     "SIMVASTATIN":"Adjust Dosage","AZATHIOPRINE":"Unknown","FLUOROURACIL":"Unknown"},
        "desc":     "All genes with loss-of-function alleles — stress test for all risk categories",
        "icon":     "🚨",
    },
]

with tab2:
    st.markdown("## 🧪 Full Test Suite — One Click Validates Everything")
    st.markdown("Runs all **4 VCF scenarios** against their expected risk labels. Proves system correctness to judges.")

    preview_cols = st.columns(4)
    for i, sc in enumerate(TEST_SUITE):
        with preview_cols[i]:
            exp_preview = " | ".join(f"{d[:4]}:{r[:3]}" for d,r in list(sc["expected"].items())[:3])
            st.markdown(f"""<div class="test-card">
                <h4>{sc['icon']} {sc['name']}</h4>
                <p>{sc['desc']}</p>
                <p style="margin-top:0.5rem;color:#60a5fa;font-size:0.75rem;">
                Drugs: {len(sc['drugs'])} | {exp_preview}...</p>
            </div>""", unsafe_allow_html=True)

    st.markdown("---")
    opt1, opt2 = st.columns([2,1])
    with opt1:
        use_llm = st.checkbox("🤖 Include LLM Explanations in tests (slower, uses API credits)", value=False)
    with opt2:
        run_all = st.button("🚀 Run All 4 Tests Now", use_container_width=True)

    if run_all:
        st.markdown("---")
        st.markdown("## 📊 Test Suite Results")

        passed_count = 0
        failed_count = 0
        all_suite_outputs = []

        for sc in TEST_SUITE:
            vcf  = load_vcf_file(sc["file"])
            pid  = f"TEST_{sc['name'][:10].replace(' ','').upper()[:8]}"
            key  = groq_api_key if use_llm else ""

            with st.spinner(f"{sc['icon']} Running: {sc['name']}..."):
                parsed_vcf, _, all_outputs, _, _ = run_pipeline(
                    vcf, sc["drugs"], pid, key, run_ix=False, gen_pdf=False)

            # Evaluate
            rows = []
            sc_pass = True
            for out in all_outputs:
                drug    = out["drug"]
                got     = out["risk_assessment"]["risk_label"]
                exp     = sc["expected"].get(drug,"")
                ok      = (got == exp) if exp else True
                pheno   = out["pharmacogenomic_profile"]["phenotype"]
                diplo   = out["pharmacogenomic_profile"]["diplotype"]
                rows.append((drug, got, exp, ok, pheno, diplo))
                if not ok: sc_pass = False

            if sc_pass: passed_count += 1
            else:       failed_count += 1
            all_suite_outputs.append({"scenario": sc["name"], "pass": sc_pass,
                                       "rows": rows, "outputs": all_outputs})

        # Summary banner
        total = passed_count + failed_count
        pct   = int(passed_count/total*100) if total else 0
        is_all_pass = failed_count == 0
        b_bg  = "#065f46" if is_all_pass else "#78350f"
        b_col = "#6ee7b7" if is_all_pass else "#fde68a"
        st.markdown(f"""
        <div style="background:{b_bg};border-radius:12px;padding:1.2rem 2rem;margin:1rem 0;text-align:center;">
          <div style="color:{b_col};font-size:1.8rem;font-weight:700;">
            {"✅ ALL TESTS PASSED" if is_all_pass else f"⚠️ {passed_count}/{total} TESTS PASSED"}
          </div>
          <div style="color:{b_col};opacity:0.8;font-size:0.95rem;">
            {pct}% pass rate · {passed_count} passed · {failed_count} failed
          </div>
        </div>""", unsafe_allow_html=True)

        # Per-scenario results
        for sc_result in all_suite_outputs:
            sc_name = sc_result["scenario"]
            sc_pass = sc_result["pass"]
            icon    = "✅" if sc_pass else "❌"
            status  = "PASS" if sc_pass else "FAIL"
            sc_icon = next((s["icon"] for s in TEST_SUITE if s["name"]==sc_name), "🧪")

            with st.expander(f"{icon} {sc_icon} {sc_name} — {status}", expanded=not sc_pass):
                hc = st.columns([2,2,2,2.5,0.8])
                for h,c in zip(["Drug","Got Risk","Expected","Diplotype / Phenotype","OK"], hc):
                    c.markdown(f"**{h}**")

                for drug, got, exp, ok, pheno, diplo in sc_result["rows"]:
                    rc = st.columns([2,2,2,2.5,0.8])
                    rc[0].markdown(f"`{drug}`")
                    badge = RISK_BADGE.get(got,"risk-badge-unknown")
                    rc[1].markdown(f'<span class="{badge}">{got}</span>', unsafe_allow_html=True)
                    rc[2].markdown(f"`{exp}`" if exp else "—")
                    rc[3].markdown(f"`{diplo}` / `{pheno}`")
                    rc[4].markdown("✅" if ok else "❌")

                st.markdown("---")
                sc_json = json.dumps(sc_result["outputs"], indent=2)
                dl1, dl2 = st.columns(2)
                with dl1:
                    sc_file = next((s["file"] for s in TEST_SUITE if s["name"]==sc_name), "test.vcf")
                    st.download_button(f"⬇️ Download JSON", data=sc_json,
                                       file_name=f"test_{sc_file.replace('.vcf','')}.json",
                                       mime="application/json", key=f"tsc_{sc_name[:10]}",
                                       use_container_width=True)
                with dl2:
                    vcf_raw = load_vcf_file(sc_file)
                    st.download_button(f"⬇️ Download VCF ({sc_file})", data=vcf_raw,
                                       file_name=sc_file, mime="text/plain",
                                       key=f"vcf_{sc_name[:10]}", use_container_width=True)

        # Full suite download
        st.markdown("---")
        full_suite_json = json.dumps(
            [{"scenario":s["scenario"],"pass":s["pass"],"results":s["outputs"]} for s in all_suite_outputs],
            indent=2)
        st.download_button("⬇️ Download Complete Test Suite JSON", data=full_suite_json,
                           file_name=f"pharmaguard_test_suite_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
                           mime="application/json", use_container_width=True)

    else:
        st.markdown("""
        <div style="text-align:center;padding:2.5rem;color:#475569;">
          <div style="font-size:3.5rem;margin-bottom:0.8rem;">🧪</div>
          <h3 style="color:#94a3b8;">One-Click Full Validation</h3>
          <p>Tests <strong>4 VCF scenarios x all 6 drugs</strong> — compares expected vs actual risk labels</p>
          <p style="font-size:0.82rem;color:#64748b;margin-top:0.5rem;">
            🧬 Mixed Variants &nbsp;·&nbsp; 🔴 UltraRapid Metabolizer
            &nbsp;·&nbsp; 🟢 All Normal &nbsp;·&nbsp; 🚨 Worst Case All PM
          </p>
        </div>""", unsafe_allow_html=True)