"""
PharmaGuard v7.0 — Unified Edition
All features fully integrated and working homogenously:
  - 4 Patient Persona Quick Demo Buttons
  - Overall Risk Command Center Banner
  - Emergency Alert Box for Critical Drugs
  - Gene Activity Heatmap (6 genes)
  - Drug Risk Comparison Table with CSV download
  - Drug Interaction Matrix
  - AI Unified Patient Narrative
  - CPIC Evidence Level Badges
  - Confidence Meter Progress Bars
  - JSON + PDF Download Buttons per drug
  - Polygenic Risk Score
  - Drug x Gene Heatmap
  - Chromosome Visualization
  - Population Frequency Bars
  - Patient Plain-English Mode
  - Parallel Test Suite
  - Prescription Safety Checker
  - Before/After Scenario Slider
  - Side-by-Side Drug Comparison
  - One-Click Clinical Note Generator
  - Voice Report Narration (TTS)
"""

import streamlit as st
import json, uuid, os, io, base64
import pandas as pd
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from dotenv import load_dotenv

load_dotenv()

from vcf_parser import parse_vcf, get_sample_vcf
from risk_engine import run_risk_assessment, get_overall_severity, DRUG_RISK_TABLE
from llm_explainer import generate_all_explanations, generate_patient_narrative
from schema import build_output_schema
from drug_interactions import run_interaction_analysis
from pdf_report import generate_pdf_report

# ── Constants ─────────────────────────────────────────────────────────────────
BASE_DIR  = os.path.dirname(os.path.abspath(__file__))
ALL_DRUGS = list(DRUG_RISK_TABLE.keys())
GENE_DRUG_MAP = {
    "CODEINE": "CYP2D6", "WARFARIN": "CYP2C9", "CLOPIDOGREL": "CYP2C19",
    "SIMVASTATIN": "SLCO1B1", "AZATHIOPRINE": "TPMT", "FLUOROURACIL": "DPYD",
}
SEV_RANK = {"none": 0, "low": 1, "moderate": 2, "high": 3, "critical": 4}

RISK_CONFIG = {
    "Safe":          {"dot": "#22c55e", "text": "#16a34a", "bg": "#f0fdf4", "border": "#bbf7d0", "emoji": "✅", "dark_bg": "#052e16", "dark_border": "#166534", "dark_text": "#4ade80"},
    "Adjust Dosage": {"dot": "#f59e0b", "text": "#b45309", "bg": "#fffbeb", "border": "#fde68a", "emoji": "⚠️", "dark_bg": "#451a03", "dark_border": "#92400e", "dark_text": "#fbbf24"},
    "Toxic":         {"dot": "#ef4444", "text": "#b91c1c", "bg": "#fef2f2", "border": "#fecaca", "emoji": "☠️", "dark_bg": "#450a0a", "dark_border": "#991b1b", "dark_text": "#f87171"},
    "Ineffective":   {"dot": "#8b5cf6", "text": "#7c3aed", "bg": "#f5f3ff", "border": "#ddd6fe", "emoji": "❌", "dark_bg": "#2e1065", "dark_border": "#6d28d9", "dark_text": "#c4b5fd"},
    "Unknown":       {"dot": "#94a3b8", "text": "#64748b", "bg": "#f8fafc", "border": "#e2e8f0", "emoji": "❓", "dark_bg": "#111827", "dark_border": "#374151", "dark_text": "#6b7280"},
}

SEV_PALETTE = {
    "none":     {"dot": "#22c55e", "bg": "#052e16", "border": "#166534", "text": "#4ade80"},
    "low":      {"dot": "#f59e0b", "bg": "#451a03", "border": "#92400e", "text": "#fbbf24"},
    "moderate": {"dot": "#f97316", "bg": "#431407", "border": "#9a3412", "text": "#fb923c"},
    "high":     {"dot": "#ef4444", "bg": "#450a0a", "border": "#991b1b", "text": "#f87171"},
    "critical": {"dot": "#dc2626", "bg": "#3b0000", "border": "#7f1d1d", "text": "#fca5a5"},
}

PHENOTYPE_COLORS = {
    "PM":      {"bg": "#7f1d1d", "text": "#fca5a5", "label": "Poor Metabolizer",          "bar": 5},
    "IM":      {"bg": "#7c2d12", "text": "#fdba74", "label": "Intermediate Metabolizer",  "bar": 45},
    "NM":      {"bg": "#14532d", "text": "#86efac", "label": "Normal Metabolizer",        "bar": 100},
    "RM":      {"bg": "#1e3a5f", "text": "#93c5fd", "label": "Rapid Metabolizer",         "bar": 115},
    "URM":     {"bg": "#78350f", "text": "#fcd34d", "label": "Ultrarapid Metabolizer",    "bar": 130},
    "Unknown": {"bg": "#1f2937", "text": "#9ca3af", "label": "Unknown",                   "bar": 0},
}

POPULATION_FREQ = {
    "CYP2D6":  {"PM": 7, "IM": 10, "NM": 77, "URM": 6},
    "CYP2C19": {"PM": 3, "IM": 26, "NM": 52, "RM": 13, "URM": 6},
    "CYP2C9":  {"PM": 1, "IM": 10, "NM": 89},
    "SLCO1B1": {"PM": 1, "IM": 15, "NM": 84},
    "TPMT":    {"PM": 0.3, "IM": 10, "NM": 90},
    "DPYD":    {"PM": 0.2, "IM": 3, "NM": 97},
}

CHROM_INFO = {
    "CYP2D6":  {"chrom": "22", "band": "q13.2",  "pos_mb": 42.5},
    "CYP2C19": {"chrom": "10", "band": "q23.33", "pos_mb": 96.7},
    "CYP2C9":  {"chrom": "10", "band": "q23.33", "pos_mb": 96.4},
    "SLCO1B1": {"chrom": "12", "band": "p12.1",  "pos_mb": 21.3},
    "TPMT":    {"chrom": "6",  "band": "p22.3",  "pos_mb": 18.1},
    "DPYD":    {"chrom": "1",  "band": "p22.1",  "pos_mb": 97.5},
}
CHROM_LEN = {"1": 248.9, "6": 170.8, "10": 133.8, "12": 133.3, "22": 50.8}

PLAIN_ENGLISH_PHENOTYPE = {
    "PM":      "Your body barely processes this medicine",
    "IM":      "Your body processes this medicine slower than average",
    "NM":      "Your body processes this medicine normally",
    "RM":      "Your body processes this medicine slightly faster than average",
    "URM":     "Your body processes this medicine dangerously fast",
    "Unknown": "Gene function unclear",
}

PLAIN_ENGLISH_RISK = {
    ("CODEINE",      "PM"):  "Your body can't convert codeine into a painkiller. You'd take it and feel nothing — or it could harm you.",
    ("CODEINE",      "URM"): "Your body converts codeine to morphine 5× faster than normal. Even one tablet could stop your breathing.",
    ("CODEINE",      "IM"):  "Codeine may work less well for you. Your doctor may need to try a different painkiller.",
    ("CODEINE",      "NM"):  "Codeine works normally for you. Standard doses should control your pain.",
    ("WARFARIN",     "PM"):  "Your blood stays thin much longer than normal. Standard doses could cause dangerous bleeding.",
    ("WARFARIN",     "IM"):  "Warfarin lasts longer in your body than average. You'll need a lower dose.",
    ("WARFARIN",     "NM"):  "Warfarin works normally for you.",
    ("CLOPIDOGREL",  "PM"):  "This heart medication doesn't get activated in your body. It won't prevent blood clots — you need a different drug.",
    ("CLOPIDOGREL",  "IM"):  "This heart medication activates less than normal. You may need a stronger alternative.",
    ("CLOPIDOGREL",  "NM"):  "This heart medication works normally for you.",
    ("SIMVASTATIN",  "PM"):  "This cholesterol drug builds up in your muscles — dangerous. You need a different medication.",
    ("SIMVASTATIN",  "IM"):  "This cholesterol drug clears more slowly. A lower dose protects your muscles.",
    ("SIMVASTATIN",  "NM"):  "This cholesterol drug works normally for you.",
    ("AZATHIOPRINE", "PM"):  "Your immune system drug builds up to toxic levels. Standard doses would damage your bone marrow.",
    ("AZATHIOPRINE", "IM"):  "You need a lower dose of this immune drug or your bone marrow could be affected.",
    ("AZATHIOPRINE", "NM"):  "This immune drug works normally for you.",
    ("FLUOROURACIL", "PM"):  "Your body cannot break down this chemotherapy. Standard doses would be fatal. You need a completely different treatment.",
    ("FLUOROURACIL", "IM"):  "This chemotherapy breaks down too slowly. You need half the normal dose.",
    ("FLUOROURACIL", "NM"):  "This chemotherapy drug works at a normal rate in your body.",
}

PERSONAS = {
    "A": {"label": "🚨 Patient A — Critical",  "file": "patient_a_critical.vcf",    "drugs": ["CODEINE", "FLUOROURACIL", "AZATHIOPRINE"], "desc": "CYP2D6 PM + DPYD PM + TPMT PM",            "color": "#ef4444", "bg": "#450a0a", "border": "#991b1b"},
    "B": {"label": "⚠️ Patient B — Warfarin PM","file": "patient_b_warfarin.vcf",   "drugs": ["WARFARIN"],                                 "desc": "CYP2C9 *2/*3 Poor Metabolizer",             "color": "#fbbf24", "bg": "#451a03", "border": "#92400e"},
    "C": {"label": "💊 Patient C — Interaction","file": "patient_c_interaction.vcf","drugs": ["CLOPIDOGREL"],                               "desc": "CYP2C19 *2/*3 Poor Metabolizer",            "color": "#c4b5fd", "bg": "#2e1065", "border": "#6d28d9"},
    "D": {"label": "✅ Patient D — All Safe",   "file": "patient_d_safe.vcf",        "drugs": ["CODEINE", "WARFARIN", "SIMVASTATIN"],       "desc": "Wildtype *1/*1 across all genes",           "color": "#4ade80", "bg": "#052e16", "border": "#166534"},
}

TEST_SUITE = [
    {"name": "Mixed Variants",          "file": "sample.vcf",                    "drugs": ["CLOPIDOGREL","CODEINE","AZATHIOPRINE"],
     "expected": {"CLOPIDOGREL":"Ineffective","CODEINE":"Ineffective","AZATHIOPRINE":"Toxic"},
     "desc": "CYP2C19 *2/*3 · CYP2D6 *4/*4 · TPMT *3B/*3C"},
    {"name": "UltraRapid Metabolizer",  "file": "test_ultrarapid_metabolizer.vcf","drugs": ["CODEINE","CLOPIDOGREL"],
     "expected": {"CODEINE":"Toxic","CLOPIDOGREL":"Safe"},
     "desc": "CYP2D6 *1xN/*1xN → URM → Codeine Toxic"},
    {"name": "All Normal Wild-type",    "file": "test_all_normal_wildtype.vcf",   "drugs": ALL_DRUGS,
     "expected": {d:"Safe" for d in ALL_DRUGS},
     "desc": "Wild-type *1/*1 across all 6 genes"},
    {"name": "Worst Case — All PM",     "file": "test_worst_case_all_pm.vcf",     "drugs": ALL_DRUGS,
     "expected": {"CODEINE":"Ineffective","CLOPIDOGREL":"Ineffective","WARFARIN":"Adjust Dosage","SIMVASTATIN":"Toxic","AZATHIOPRINE":"Toxic","FLUOROURACIL":"Toxic"},
     "desc": "Loss-of-function alleles across all 6 genes"},
]

# ── Page Config ───────────────────────────────────────────────────────────────
st.set_page_config(page_title="PharmaGuard", page_icon="🧬", layout="wide", initial_sidebar_state="collapsed")

# ── Global CSS ────────────────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Syne:wght@400;500;600;700;800&family=DM+Mono:wght@300;400;500&family=Fraunces:ital,wght@0,300;0,400;0,600;1,300;1,400&display=swap');

*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
html, body, [class*="css"] {
    font-family: 'Syne', sans-serif !important;
    background: #060608 !important;
    color: #e8e8f0 !important;
}
.stApp { background: #060608 !important; }
.main .block-container { padding: 0 2.5rem 6rem !important; max-width: 1280px !important; }
#MainMenu, footer, header { visibility: hidden; }

/* Animations */
@keyframes pulse-ring { 0%,100%{box-shadow:0 0 0 0 rgba(220,38,38,.4)} 50%{box-shadow:0 0 0 10px rgba(220,38,38,0)} }
@keyframes fade-in { from{opacity:0;transform:translateY(8px)} to{opacity:1;transform:translateY(0)} }
@keyframes glow-pulse { 0%,100%{opacity:1} 50%{opacity:.4} }
@keyframes slide-in { from{transform:translateX(-12px);opacity:0} to{transform:translateX(0);opacity:1} }
@keyframes number-up { from{transform:scale(.8);opacity:0} to{transform:scale(1);opacity:1} }

/* Nav */
.pg-nav {
    display:flex; align-items:center; justify-content:space-between;
    padding: 1.75rem 0 2rem; border-bottom: 1px solid #1a1a24; margin-bottom: 0;
}
.pg-logo { font-family:'Fraunces',serif; font-size:1.9rem; font-weight:300; color:#e8e8f0; letter-spacing:-.03em; }
.pg-logo strong { font-weight:600; color:#fff; }
.pg-logo em { font-style:italic; color:#7c6aff; }
.pg-tags { display:flex; gap:.5rem; align-items:center; }
.pg-tag { font-family:'DM Mono',monospace; font-size:.62rem; letter-spacing:.1em; text-transform:uppercase;
    padding:4px 12px; border-radius:100px; border:1px solid; }
.pg-tag-default { color:#4a4a5a; border-color:#1e1e2e; background:#0c0c14; }
.pg-tag-hot { color:#7c6aff; border-color:#3d2f8f; background:#1a1030; }

/* Tab bar */
.stTabs [data-baseweb="tab-list"] { background:transparent !important; border-bottom:1px solid #1a1a24 !important;
    gap:0 !important; padding:0 !important; margin-bottom:2rem !important; box-shadow:none !important; }
.stTabs [data-baseweb="tab"] { font-family:'DM Mono',monospace !important; font-size:.68rem !important;
    letter-spacing:.12em !important; text-transform:uppercase !important; color:#2e2e3e !important;
    padding:.9rem 1.5rem !important; background:transparent !important; border:none !important;
    border-bottom:2px solid transparent !important; border-radius:0 !important; transition:all .15s !important; }
.stTabs [aria-selected="true"] { color:#e8e8f0 !important; border-bottom-color:#7c6aff !important; }
.stTabs [data-baseweb="tab-panel"] { padding-top:0 !important; }

/* Persona cards */
.persona-grid { display:grid; grid-template-columns:repeat(4,1fr); gap:.75rem; margin-bottom:2.5rem; }
.persona-card { border-radius:12px; border:1px solid; padding:1rem 1.1rem; cursor:pointer;
    transition:all .2s; position:relative; overflow:hidden; }
.persona-card:hover { transform:translateY(-3px); }
.persona-card-label { font-size:.9rem; font-weight:700; margin-bottom:.3rem; }
.persona-card-desc { font-family:'DM Mono',monospace; font-size:.58rem; letter-spacing:.04em; opacity:.65; line-height:1.6; }

/* Steps bar */
.steps-bar { display:flex; border:1px solid #1a1a24; border-radius:10px;
    overflow:hidden; margin-bottom:2.5rem; background:#0a0a10; }
.step { flex:1; padding:.9rem 1.2rem; border-right:1px solid #1a1a24; }
.step:last-child { border-right:none; }
.step-n { font-family:'DM Mono',monospace; font-size:.58rem; letter-spacing:.12em; text-transform:uppercase; color:#7c6aff; margin-bottom:3px; }
.step-l { font-size:.875rem; font-weight:600; color:#2a2a38; }
.step.active .step-l { color:#e8e8f0; }

/* Risk banner */
.risk-banner { border-radius:14px; padding:1.5rem 2rem; margin-bottom:1.5rem;
    border:1px solid; animation:fade-in .4s ease; position:relative; overflow:hidden; }
.risk-banner::after { content:''; position:absolute; top:-40%; right:-5%; width:240px; height:240px;
    border-radius:50%; opacity:.04; pointer-events:none; background:currentColor; }
.risk-banner-eyebrow { font-family:'DM Mono',monospace; font-size:.6rem; letter-spacing:.15em;
    text-transform:uppercase; opacity:.6; margin-bottom:.4rem; }
.risk-banner-headline { font-family:'Fraunces',serif; font-size:2.4rem; font-weight:300; line-height:1; margin-bottom:.2rem; }
.risk-banner-stats { display:grid; grid-template-columns:repeat(4,1fr); gap:1rem; margin-top:1.25rem; padding-top:1rem; border-top:1px solid; }
.rbs-num { font-family:'Fraunces',serif; font-size:1.8rem; line-height:1; margin-bottom:.15rem; }
.rbs-key { font-family:'DM Mono',monospace; font-size:.58rem; letter-spacing:.1em; text-transform:uppercase; opacity:.55; }

/* Emergency alert */
.emergency { background:#16020a; border:1.5px solid #dc2626; border-radius:12px;
    padding:1.25rem 1.5rem; margin-bottom:1rem; animation:pulse-ring 2.5s infinite; }
.emergency-head { display:flex; align-items:center; gap:.75rem; margin-bottom:.5rem; }
.emergency-icon { font-size:1.3rem; }
.emergency-drug { font-size:1.05rem; font-weight:700; color:#fca5a5; }
.emergency-note { font-size:.9rem; color:#f87171; line-height:1.6; margin-bottom:.5rem; }
.emergency-cta { font-family:'DM Mono',monospace; font-size:.68rem; color:#ef4444; font-weight:700; text-transform:uppercase; letter-spacing:.06em; }

/* Gene heatmap row */
.gene-row { display:grid; grid-template-columns:repeat(6,1fr); gap:.5rem; margin-bottom:1.5rem; }
.gene-box { border-radius:10px; padding:.875rem .75rem; text-align:center; border:1px solid transparent; }
.gene-name { font-family:'DM Mono',monospace; font-size:.68rem; font-weight:600; margin-bottom:.35rem; }
.gene-bar-track { height:3px; border-radius:2px; background:#0a0a10; margin:.35rem 0; }
.gene-bar-fill { height:100%; border-radius:2px; }
.gene-pheno { font-family:'DM Mono',monospace; font-size:.58rem; letter-spacing:.04em; }

/* Drug comparison table */
.dtable { border:1px solid #14141e; border-radius:12px; overflow:hidden; margin-bottom:1.5rem; }
.dtable-head { display:grid; grid-template-columns:1.3fr 1.2fr 1fr 1fr 1fr 1fr;
    background:#08080e; border-bottom:1px solid #14141e; padding:0 .5rem; }
.dtable-hcell { font-family:'DM Mono',monospace; font-size:.6rem; letter-spacing:.1em; text-transform:uppercase;
    color:#30303c; padding:.7rem .85rem; }
.dtable-row { display:grid; grid-template-columns:1.3fr 1.2fr 1fr 1fr 1fr 1fr;
    border-bottom:1px solid #0e0e18; padding:0 .5rem; background:#0d0d16; transition:background .15s; }
.dtable-row:last-child { border-bottom:none; }
.dtable-row:hover { background:#121220; }
.dtable-cell { font-family:'DM Mono',monospace; font-size:.77rem; color:#707080; padding:.7rem .85rem; display:flex; align-items:center; }

/* Polygenic score */
.pgx { background:linear-gradient(135deg,#0c0c1a,#140c24,#0a1020);
    border:1px solid #2a2040; border-radius:16px; padding:2rem; margin-bottom:1.5rem;
    position:relative; overflow:hidden; }
.pgx::before { content:''; position:absolute; top:-30%; right:-10%; width:280px; height:280px;
    background:radial-gradient(circle,#7c6aff0f,transparent 70%); pointer-events:none; }
.pgx-eye { font-family:'DM Mono',monospace; font-size:.6rem; letter-spacing:.15em; text-transform:uppercase; color:#7c6aff; margin-bottom:.4rem; }
.pgx-score { font-family:'Fraunces',serif; font-size:4.5rem; font-weight:300; line-height:1; }
.pgx-label { font-size:.875rem; color:#40405a; margin-bottom:1.25rem; }
.pgx-track { height:5px; background:#14141e; border-radius:3px; overflow:hidden; margin-bottom:.35rem; }
.pgx-fill { height:100%; border-radius:3px; transition:width .8s ease; }
.pgx-scale { display:flex; justify-content:space-between; font-family:'DM Mono',monospace; font-size:.55rem; color:#20202e; }
.pgx-pills { display:flex; flex-wrap:wrap; gap:.4rem; margin-top:1rem; }
.pgx-pill { font-family:'DM Mono',monospace; font-size:.6rem; padding:3px 10px; border-radius:100px; border:1px solid; }

/* Heatmap */
.hm-wrap { background:#0d0d14; border:1px solid #14141e; border-radius:14px;
    padding:1.5rem; margin-bottom:1.5rem; overflow-x:auto; }
.hm-eye { font-family:'DM Mono',monospace; font-size:.6rem; letter-spacing:.15em;
    text-transform:uppercase; color:#30303c; margin-bottom:1.2rem; }
.hm-grid { display:grid; gap:4px; }
.hm-cell { border-radius:6px; display:flex; flex-direction:column; align-items:center;
    justify-content:center; padding:.6rem .4rem; min-height:56px; cursor:default;
    transition:transform .15s, box-shadow .15s; border:1px solid; }
.hm-cell:hover { transform:scale(1.08); z-index:10; position:relative; }
.hm-dname { font-family:'DM Mono',monospace; font-size:.56rem; font-weight:700; margin-bottom:2px; }
.hm-drisk { font-family:'DM Mono',monospace; font-size:.52rem; opacity:.8; }
.hm-header { font-family:'DM Mono',monospace; font-size:.58rem; letter-spacing:.08em;
    color:#25253a; display:flex; align-items:center; justify-content:center; min-height:56px; }
.hm-legend { display:flex; gap:1rem; margin-top:.875rem; flex-wrap:wrap; }
.hm-legend-item { font-family:'DM Mono',monospace; font-size:.58rem; display:flex; align-items:center; gap:5px; color:#404050; }
.hm-dot { width:9px; height:9px; border-radius:2px; display:inline-block; }

/* Chromosome */
.chrom-wrap { background:#0d0d14; border:1px solid #14141e; border-radius:12px;
    padding:1.2rem 1.5rem; }
.chrom-eye { font-family:'DM Mono',monospace; font-size:.6rem; letter-spacing:.1em;
    text-transform:uppercase; color:#30303c; margin-bottom:.75rem; }
.chrom-row { display:flex; align-items:center; gap:.75rem; margin-bottom:.45rem; }
.chrom-lbl { font-family:'DM Mono',monospace; font-size:.62rem; color:#40405a; width:20px; text-align:right; flex-shrink:0; }
.chrom-bar { flex:1; height:13px; background:#12121c; border-radius:7px; position:relative; overflow:visible; }
.chrom-body { position:absolute; inset:0; background:linear-gradient(90deg,#1e1e2e,#2a2a3e,#1e1e2e); border-radius:7px; }
.chrom-marker { position:absolute; top:-4px; width:3px; height:21px; border-radius:2px; transform:translateX(-50%); }
.chrom-gene-lbl { font-family:'DM Mono',monospace; font-size:.58rem; color:#7c7c90; width:60px; flex-shrink:0; }
.chrom-band { font-family:'DM Mono',monospace; font-size:.55rem; color:#28283a; }

/* Population freq */
.pop-wrap { background:#0c0c12; border:1px solid #14141e; border-radius:10px;
    padding:1rem 1.2rem; margin-bottom:.75rem; }
.pop-eye { font-family:'DM Mono',monospace; font-size:.58rem; letter-spacing:.1em;
    text-transform:uppercase; color:#30303c; margin-bottom:.5rem; }
.pop-row { display:flex; align-items:center; gap:.75rem; margin-bottom:.35rem; }
.pop-ph { font-family:'DM Mono',monospace; font-size:.65rem; color:#40405a; width:100px; flex-shrink:0; }
.pop-track { flex:1; height:4px; background:#14141e; border-radius:2px; overflow:hidden; }
.pop-fill { height:100%; border-radius:2px; }
.pop-pct { font-family:'DM Mono',monospace; font-size:.6rem; width:35px; text-align:right; }
.pop-you { font-family:'DM Mono',monospace; font-size:.55rem; color:#7c6aff; margin-left:3px; }

/* Interaction matrix */
.ix-matrix-grid { display:grid; gap:3px; }
.ix-cell { border-radius:4px; display:flex; align-items:center; justify-content:center;
    min-height:46px; font-family:'DM Mono',monospace; font-size:.58rem; text-align:center;
    padding:.3rem; transition:transform .1s; border:1px solid; cursor:pointer; }
.ix-cell:hover { transform:scale(1.06); z-index:5; position:relative; }
.ix-head { font-family:'DM Mono',monospace; font-size:.58rem; letter-spacing:.06em;
    color:#25253a; display:flex; align-items:center; justify-content:center; min-height:46px; }

/* Drug result cards */
.rcard { border:1px solid #14141e; border-radius:14px; background:#0d0d14; margin-bottom:1.25rem; overflow:hidden; animation:fade-in .4s ease; }
.rcard-top { padding:1.25rem 1.5rem; display:flex; align-items:center; justify-content:space-between; border-bottom:1px solid #12121c; }
.rcard-left { display:flex; align-items:center; gap:.875rem; }
.rcard-dot { width:10px; height:10px; border-radius:50%; flex-shrink:0; }
.rcard-name { font-size:1.05rem; font-weight:700; letter-spacing:-.01em; }
.rcard-meta { font-family:'DM Mono',monospace; font-size:.68rem; color:#30303c; margin-top:2px; letter-spacing:.04em; }
.rcard-badge { font-family:'DM Mono',monospace; font-size:.68rem; font-weight:600; letter-spacing:.08em;
    text-transform:uppercase; padding:5px 14px; border-radius:100px; border:1px solid; }
.rcard-body { padding:1.25rem 1.5rem; }
.mc-row { display:grid; grid-template-columns:repeat(4,1fr); gap:1px; background:#10101a; border-radius:8px; overflow:hidden; margin-bottom:1.25rem; }
.mc-cell { background:#0c0c12; padding:.875rem 1rem; }
.mc-key { font-family:'DM Mono',monospace; font-size:.6rem; letter-spacing:.1em; text-transform:uppercase; color:#28283a; margin-bottom:.3rem; }
.mc-val { font-size:1rem; font-weight:600; }
.conf-row { display:grid; grid-template-columns:1fr 1fr; gap:1rem; margin-bottom:1.25rem; }
.conf-item {}
.conf-lbl { font-family:'DM Mono',monospace; font-size:.6rem; letter-spacing:.08em; text-transform:uppercase; color:#28283a; margin-bottom:4px; display:flex; justify-content:space-between; }
.conf-track { height:3px; background:#14141e; border-radius:2px; overflow:hidden; }
.conf-fill { height:100%; border-radius:2px; transition:width .6s ease; }
.vtable { width:100%; border-collapse:collapse; }
.vtable th { font-family:'DM Mono',monospace; font-size:.62rem; letter-spacing:.09em;
    text-transform:uppercase; color:#28283a; padding:0 .6rem .5rem; text-align:left; border-bottom:1px solid #14141e; }
.vtable td { font-family:'DM Mono',monospace; font-size:.77rem; color:#909090;
    padding:.55rem .6rem; border-bottom:1px solid #0e0e18; }
.vtable tbody tr:last-child td { border-bottom:none; }
.v-rsid { color:#3b82f6 !important; }
.v-star { color:#8b5cf6 !important; }
.v-nofunc { color:#ef4444 !important; }
.v-dec { color:#f59e0b !important; }
.v-norm { color:#22c55e !important; }

/* AI boxes */
.ai-narrative { background:linear-gradient(135deg,#0a0a14,#100820); border:1px solid #2a1a50;
    border-radius:14px; padding:1.5rem; margin-bottom:1.5rem; }
.ai-nar-head { display:flex; align-items:center; gap:.75rem; margin-bottom:1rem; }
.ai-nar-badge { font-family:'DM Mono',monospace; font-size:.6rem; letter-spacing:.08em; text-transform:uppercase;
    background:#1e0a3a; border:1px solid #4c1d95; color:#a78bfa; padding:3px 10px; border-radius:4px; }
.ai-nar-title { font-size:.95rem; font-weight:600; color:#c4b5fd; }
.ai-nar-text { font-size:.95rem; line-height:1.85; color:#a0a0b8; }
.ai-block { border:1px solid #14141e; border-radius:10px; overflow:hidden; margin-top:1.25rem; }
.ai-block-head { padding:.65rem 1rem; background:#0a0a10; border-bottom:1px solid #14141e; display:flex; align-items:center; gap:.65rem; }
.ai-badge { font-family:'DM Mono',monospace; font-size:.6rem; letter-spacing:.08em; text-transform:uppercase;
    color:#40405a; background:#14141e; padding:3px 9px; border-radius:4px; }
.ai-section { padding:.9rem 1rem; border-bottom:1px solid #0e0e18; }
.ai-section:last-child { border-bottom:none; }
.ai-section-lbl { font-family:'DM Mono',monospace; font-size:.62rem; letter-spacing:.09em; text-transform:uppercase; color:#28283a; margin-bottom:.4rem; }
.ai-section-txt { font-size:.925rem; line-height:1.8; color:#909090; }

/* Rec box */
.rec-box { border-radius:8px; border:1px solid; padding:1rem 1.2rem; margin-bottom:1rem; }
.rec-lbl { font-family:'DM Mono',monospace; font-size:.62rem; letter-spacing:.09em; text-transform:uppercase; margin-bottom:.4rem; }
.rec-txt { font-size:.95rem; line-height:1.75; color:#909090; }
.alt-chips { display:flex; flex-wrap:wrap; gap:.4rem; }
.alt-chip { font-family:'DM Mono',monospace; font-size:.7rem; color:#909090;
    background:#12121c; border:1px solid #1e1e2e; border-radius:100px; padding:4px 12px; }

/* CPIC badge */
.cpic { font-family:'DM Mono',monospace; font-size:.58rem; letter-spacing:.05em;
    background:#1a0a02; border:1px solid #78350f; color:#fbbf24; padding:2px 8px;
    border-radius:4px; display:inline-block; margin-left:.5rem; }

/* Clinical note */
.note-box { background:#0a0a14; border:1px solid #1e1e2e; border-radius:10px; padding:1.25rem 1.5rem; }
.note-box pre { font-family:'DM Mono',monospace; font-size:.8rem; color:#9090b0; line-height:1.7; white-space:pre-wrap; word-break:break-word; }

/* Before/after slider */
.ba-wrap { display:grid; grid-template-columns:1fr 1fr; gap:1px; background:#14141e; border:1px solid #14141e; border-radius:12px; overflow:hidden; margin-bottom:1.5rem; }
.ba-side { padding:1.5rem; }
.ba-label { font-family:'DM Mono',monospace; font-size:.62rem; letter-spacing:.12em; text-transform:uppercase; margin-bottom:.75rem; }
.ba-scenario { font-size:1rem; font-weight:600; margin-bottom:.35rem; }
.ba-outcome { font-size:.875rem; line-height:1.6; }

/* Patient mode */
.patient-card { border:1px solid; border-radius:14px; padding:1.5rem; margin-bottom:1rem; animation:slide-in .3s ease; }
.patient-verdict { font-size:1.1rem; font-weight:700; line-height:1.6; margin-bottom:.5rem; }
.patient-plain { font-size:.9rem; line-height:1.8; }
.patient-action { display:flex; align-items:flex-start; gap:.65rem; background:#08080e; border:1px solid #14141e;
    border-radius:8px; padding:.875rem 1rem; margin-top:.75rem; }

/* Section labels */
.section-lbl { font-family:'DM Mono',monospace; font-size:.65rem; letter-spacing:.15em;
    text-transform:uppercase; color:#30303c; margin-bottom:.75rem; }
.h-rule { border:none; border-top:1px solid #14141e; margin:1.25rem 0; }

/* Prescription checker */
.rx-result { border-radius:10px; padding:1.2rem 1.5rem; margin-top:1rem; border:1px solid; }

/* Empty state */
.empty { text-align:center; padding:5rem 2rem; border:1px dashed #141420; border-radius:14px; background:#08080c; }
.empty-icon { font-family:'Fraunces',serif; font-size:3rem; color:#1e1e28; margin-bottom:1rem; }
.empty-title { font-size:1.1rem; font-weight:600; color:#20202e; margin-bottom:.5rem; }
.empty-hint { font-family:'DM Mono',monospace; font-size:.65rem; color:#18181e; letter-spacing:.06em; line-height:2.2; }

/* Test suite */
.test-card { background:#0d0d14; border:1px solid #14141e; border-radius:10px; padding:1.2rem; }
.test-name { font-size:.95rem; font-weight:600; margin-bottom:.35rem; }
.test-desc { font-family:'DM Mono',monospace; font-size:.63rem; color:#25253a; margin-bottom:.875rem; line-height:1.7; }
.test-row { display:flex; align-items:center; gap:.5rem; margin-bottom:3px; }
.test-drug-lbl { font-family:'DM Mono',monospace; font-size:.68rem; color:#25253a; width:90px; flex-shrink:0; }

/* Buttons */
.stButton>button { background:#e8e8f0 !important; color:#0a0a14 !important; border:none !important;
    border-radius:8px !important; font-family:'Syne',sans-serif !important; font-weight:600 !important;
    font-size:.875rem !important; padding:.65rem 1.75rem !important; transition:opacity .15s !important; }
.stButton>button:hover { opacity:.85 !important; }
.stDownloadButton>button { background:#0d0d14 !important; color:#606070 !important;
    border:1px solid #1e1e2e !important; border-radius:8px !important; font-family:'DM Mono',monospace !important;
    font-size:.72rem !important; letter-spacing:.04em !important; padding:.5rem 1rem !important; transition:all .15s !important; }
.stDownloadButton>button:hover { border-color:#e8e8f0 !important; color:#e8e8f0 !important; }
div[data-testid="stExpander"] { background:#0d0d14 !important; border:1px solid #14141e !important;
    border-radius:8px !important; box-shadow:none !important; margin-bottom:.5rem !important; }
div[data-testid="stExpander"] summary { font-family:'DM Mono',monospace !important; font-size:.7rem !important;
    letter-spacing:.07em !important; color:#30303c !important; padding:.8rem 1rem !important; }
.stMultiSelect span[data-baseweb="tag"] { background:#14141e !important; color:#909090 !important;
    border:1px solid #1e1e2e !important; font-family:'DM Mono',monospace !important; font-size:.7rem !important; border-radius:4px !important; }
.stCheckbox label p { font-family:'Syne',sans-serif !important; font-size:.9rem !important; color:#606070 !important; }
.stTextInput>div>div>input { border-radius:8px !important; background:#0c0c14 !important; border:1px solid #1e1e2e !important;
    color:#e8e8f0 !important; font-family:'DM Mono',monospace !important; font-size:.875rem !important; }
.stFileUploader>div { border-radius:8px !important; border:1px dashed #1e1e2e !important; background:#0a0a10 !important; }
[data-testid="stSidebar"] { background:#060608 !important; border-right:1px solid #14141e !important; }
</style>
""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def load_vcf_file(filename):
    path = os.path.join(BASE_DIR, "sample_data", filename)
    return open(path).read() if os.path.exists(path) else get_sample_vcf()


def run_pipeline(vcf_content, drugs, pid, groq_key, run_ix=True, gen_pdf=True, skip_llm=False):
    parsed = parse_vcf(vcf_content)
    results = run_risk_assessment(parsed, drugs)
    results = generate_all_explanations(groq_key, results, skip_llm=skip_llm)
    outputs = [
        build_output_schema(patient_id=pid, drug=r["drug"], result=r,
                            parsed_vcf=parsed, llm_exp=r.get("llm_explanation", {}))
        for r in results
    ]
    ix = run_interaction_analysis(drugs, results) if run_ix and len(drugs) > 1 else None
    pdf = None
    if gen_pdf:
        try:
            pdf = generate_pdf_report(pid, outputs, parsed)
        except Exception:
            pass
    return parsed, results, outputs, ix, pdf


def func_cls(status):
    s = (status or "").lower()
    if "no_function" in s: return "v-nofunc"
    if "decreased"   in s: return "v-dec"
    return "v-norm"


# ══════════════════════════════════════════════════════════════════════════════
# POLYGENIC RISK SCORE
# ══════════════════════════════════════════════════════════════════════════════

def compute_pgx_score(all_outputs):
    SEV_S  = {"none": 0, "low": 20, "moderate": 45, "high": 70, "critical": 100}
    RISK_S = {"Safe": 0, "Adjust Dosage": 35, "Toxic": 85, "Ineffective": 70, "Unknown": 20}
    W      = {"FLUOROURACIL": 1.4, "AZATHIOPRINE": 1.3, "CLOPIDOGREL": 1.3,
               "WARFARIN": 1.2, "CODEINE": 1.1, "SIMVASTATIN": 1.0}
    if not all_outputs: return 0, "No data", []
    tw = ws = 0
    breakdown = []
    for o in all_outputs:
        drug = o["drug"]
        sev  = o["risk_assessment"]["severity"]
        rl   = o["risk_assessment"]["risk_label"]
        gene = o["pharmacogenomic_profile"]["primary_gene"]
        ph   = o["pharmacogenomic_profile"]["phenotype"]
        sc   = (SEV_S.get(sev, 0) + RISK_S.get(rl, 0)) / 2
        wt   = W.get(drug, 1.0)
        ws  += sc * wt; tw += wt
        breakdown.append((gene, drug, ph, rl, sc))
    final = min(100, int(ws / tw)) if tw else 0
    label = ["Low Risk","Moderate Risk","High Risk","Very High Risk","Critical Risk"][min(4, final // 20)]
    return final, label, breakdown


def render_pgx_score(all_outputs):
    score, label, breakdown = compute_pgx_score(all_outputs)
    colors = ["#22c55e","#f59e0b","#f97316","#ef4444","#dc2626"]
    color  = colors[min(4, score // 20)]
    pills  = ""
    for gene, _, ph, rl, _ in breakdown:
        rc = RISK_CONFIG.get(rl, RISK_CONFIG["Unknown"])
        pills += f'<span class="pgx-pill" style="background:{rc["dark_bg"]};border-color:{rc["dark_border"]};color:{rc["dark_text"]};">{gene} · {ph}</span>'
    st.markdown(f"""
    <div class="pgx">
      <div class="pgx-eye">Polygenic Risk Score</div>
      <div class="pgx-score" style="color:{color};">{score}</div>
      <div class="pgx-label">{label} — composite across {len(all_outputs)} drug{"s" if len(all_outputs)!=1 else ""}</div>
      <div class="pgx-track">
        <div class="pgx-fill" style="width:{score}%;background:linear-gradient(90deg,{color}88,{color});box-shadow:0 0 14px {color}44;"></div>
      </div>
      <div class="pgx-scale"><span>0</span><span>50 — High</span><span>100 — Critical</span></div>
      <div class="pgx-pills">{pills}</div>
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# RISK COMMAND CENTER BANNER
# ══════════════════════════════════════════════════════════════════════════════

def render_risk_banner(all_outputs, parsed):
    sev = max((o["risk_assessment"]["severity"] for o in all_outputs),
              key=lambda s: SEV_RANK.get(s, 0), default="none")
    sp  = SEV_PALETTE.get(sev, SEV_PALETTE["none"])
    SEV_EMO   = {"none":"✅","low":"💛","moderate":"🟠","high":"🔴","critical":"🚨"}
    SEV_LABEL = {"none":"Clear","low":"Low","moderate":"Moderate","high":"High","critical":"Critical"}
    emoji = SEV_EMO.get(sev,"")
    label = SEV_LABEL.get(sev,"Unknown")
    high_crit = sum(1 for o in all_outputs if o["risk_assessment"]["severity"] in ("high","critical"))
    ndrugs    = len(all_outputs)
    ngenes    = len(parsed.get("detected_genes",[]))
    nvars     = parsed.get("total_variants",0)
    bc_border = sp["border"] if sev != "none" else "#1e2e1e"
    st.markdown(f"""
    <div class="risk-banner" style="background:{sp['bg']};border-color:{bc_border};color:{sp['text']};">
      <div class="risk-banner-eyebrow">Risk Command Center</div>
      <div class="risk-banner-headline">{emoji} {label} Risk</div>
      <div class="risk-banner-stats" style="border-color:{bc_border}88;">
        <div><div class="rbs-num">{ndrugs}</div><div class="rbs-key">Drugs</div></div>
        <div><div class="rbs-num" style="color:{'#f87171' if high_crit else sp['text']};">{high_crit}</div><div class="rbs-key">High/Critical</div></div>
        <div><div class="rbs-num">{ngenes}</div><div class="rbs-key">Genes</div></div>
        <div><div class="rbs-num">{nvars}</div><div class="rbs-key">Variants</div></div>
      </div>
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# EMERGENCY ALERTS
# ══════════════════════════════════════════════════════════════════════════════

def render_emergency_alerts(all_outputs):
    for o in all_outputs:
        if o["risk_assessment"]["severity"] == "critical":
            drug = o["drug"]
            note = o["clinical_recommendation"]["dosing_recommendation"][:220]
            st.markdown(f"""
            <div class="emergency">
              <div class="emergency-head">
                <span class="emergency-icon">🚨</span>
                <span class="emergency-drug">CRITICAL ALERT — {drug}</span>
              </div>
              <div class="emergency-note">{note}{"…" if len(o["clinical_recommendation"]["dosing_recommendation"]) > 220 else ""}</div>
              <div class="emergency-cta">⚡ Contact prescribing physician immediately</div>
            </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# GENE ACTIVITY HEATMAP
# ══════════════════════════════════════════════════════════════════════════════

def render_gene_heatmap(all_outputs):
    GENE_ORDER = ["CYP2D6","CYP2C19","CYP2C9","SLCO1B1","TPMT","DPYD"]
    gp = {o["pharmacogenomic_profile"]["primary_gene"]: o["pharmacogenomic_profile"]["phenotype"]
          for o in all_outputs}
    boxes = ""
    for g in GENE_ORDER:
        ph = gp.get(g, "Unknown")
        pc = PHENOTYPE_COLORS.get(ph, PHENOTYPE_COLORS["Unknown"])
        bar = min(100, pc["bar"])
        boxes += f"""
        <div class="gene-box" style="background:{pc['bg']}18;border-color:{pc['text']}22;">
          <div class="gene-name" style="color:{pc['text']};">{g}</div>
          <div class="gene-bar-track"><div class="gene-bar-fill" style="width:{bar}%;background:{pc['text']};"></div></div>
          <div class="gene-pheno" style="color:{pc['text']};">{ph}</div>
        </div>"""
    st.markdown(f"""
    <div style="margin-bottom:.5rem;">
      <div class="section-lbl">Gene Activity Overview</div>
      <div class="gene-row">{boxes}</div>
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# DRUG RISK COMPARISON TABLE
# ══════════════════════════════════════════════════════════════════════════════

def render_drug_table(all_outputs, pid):
    rows_html = ""
    table_data = []
    for o in all_outputs:
        drug  = o["drug"]
        rl    = o["risk_assessment"]["risk_label"]
        sev   = o["risk_assessment"]["severity"]
        conf  = o["risk_assessment"]["confidence_score"]
        gene  = o["pharmacogenomic_profile"]["primary_gene"]
        ph    = o["pharmacogenomic_profile"]["phenotype"]
        rc    = RISK_CONFIG.get(rl, RISK_CONFIG["Unknown"])
        sp    = SEV_PALETTE.get(sev, SEV_PALETTE["none"])
        rows_html += f"""<div class="dtable-row">
          <div class="dtable-cell" style="font-weight:700;color:#e8e8f0;">{drug.title()}</div>
          <div class="dtable-cell"><span style="display:inline-flex;align-items:center;gap:6px;">
            <span style="width:7px;height:7px;border-radius:50%;background:{rc['dot']};flex-shrink:0;"></span>
            {rc['emoji']} {rl}</span></div>
          <div class="dtable-cell" style="color:{sp['text']};">{sev.title()}</div>
          <div class="dtable-cell" style="color:#606070;">{gene}</div>
          <div class="dtable-cell" style="color:{rc['dot']};">{ph}</div>
          <div class="dtable-cell">
            <div style="flex:1;height:3px;background:#14141e;border-radius:2px;overflow:hidden;margin-right:8px;">
              <div style="width:{conf*100:.0f}%;height:100%;background:{rc['dot']};border-radius:2px;"></div>
            </div>
            <span style="font-size:.65rem;color:#40405a;">{conf:.0%}</span>
          </div>
        </div>"""
        table_data.append({"Drug": drug, "Risk": rl, "Severity": sev,
                           "Gene": gene, "Phenotype": ph, "Confidence": f"{conf:.0%}"})
    st.markdown(f"""
    <div class="dtable">
      <div class="dtable-head">
        <div class="dtable-hcell">Drug</div><div class="dtable-hcell">Risk Label</div>
        <div class="dtable-hcell">Severity</div><div class="dtable-hcell">Gene</div>
        <div class="dtable-hcell">Phenotype</div><div class="dtable-hcell">Confidence</div>
      </div>
      {rows_html}
    </div>""", unsafe_allow_html=True)
    df = pd.DataFrame(table_data)
    st.download_button("⬇ Download as CSV", data=df.to_csv(index=False),
                       file_name=f"pharmaguard_{pid}_comparison.csv", mime="text/csv",
                       key=f"csv_{pid}")


# ══════════════════════════════════════════════════════════════════════════════
# DRUG x GENE HEATMAP
# ══════════════════════════════════════════════════════════════════════════════

def render_heatmap(all_outputs):
    DRUG_ORDER = ["CODEINE","WARFARIN","CLOPIDOGREL","SIMVASTATIN","AZATHIOPRINE","FLUOROURACIL"]
    GENE_ORDER = ["CYP2D6","CYP2C9","CYP2C19","SLCO1B1","TPMT","DPYD"]
    DG = {"CODEINE":"CYP2D6","WARFARIN":"CYP2C9","CLOPIDOGREL":"CYP2C19",
          "SIMVASTATIN":"SLCO1B1","AZATHIOPRINE":"TPMT","FLUOROURACIL":"DPYD"}
    rmap  = {o["drug"]: o for o in all_outputs}
    drugs = [d for d in DRUG_ORDER if d in rmap]
    if not drugs: return
    n = len(drugs)
    HM_COLS = {
        "Safe":          {"bg":"#052e16","text":"#4ade80","border":"#166534"},
        "Adjust Dosage": {"bg":"#451a03","text":"#fbbf24","border":"#92400e"},
        "Toxic":         {"bg":"#450a0a","text":"#f87171","border":"#991b1b"},
        "Ineffective":   {"bg":"#2e1065","text":"#c4b5fd","border":"#6d28d9"},
        "Unknown":       {"bg":"#111827","text":"#6b7280","border":"#374151"},
        "N/A":           {"bg":"#08080e","text":"#1a1a28","border":"#10101e"},
    }
    hdrs = '<div class="hm-header" style="padding:.4rem;"></div>'
    for d in drugs:
        hdrs += f'<div class="hm-header">{d[:5]}</div>'
    rows = ""
    for gene in GENE_ORDER:
        rows += f'<div class="hm-header" style="justify-content:flex-end;padding-right:.5rem;font-size:.55rem;">{gene}</div>'
        for d in drugs:
            if DG.get(d) == gene and d in rmap:
                o  = rmap[d]; rl = o["risk_assessment"]["risk_label"]
                ph = o["pharmacogenomic_profile"]["phenotype"]
                mc = HM_COLS.get(rl, HM_COLS["Unknown"])
                sh = {"Adjust Dosage":"Adjust","Ineffective":"Ineffect.","Unknown":"?"}.get(rl,rl)
                rows += f'<div class="hm-cell" style="background:{mc["bg"]};border-color:{mc["border"]};" title="{d}×{gene}: {rl} ({ph})"><div class="hm-dname" style="color:{mc["text"]};">{sh}</div><div class="hm-drisk" style="color:{mc["text"]};">{ph}</div></div>'
            else:
                mc = HM_COLS["N/A"]
                rows += f'<div class="hm-cell" style="background:{mc["bg"]};border-color:{mc["border"]};"><div class="hm-drisk" style="color:{mc["text"]};">—</div></div>'
    legend = "".join(f'<div class="hm-legend-item"><span class="hm-dot" style="background:{HM_COLS[r]["bg"]};border:1px solid {HM_COLS[r]["border"]};"></span>{r}</div>'
                     for r in ["Safe","Adjust Dosage","Toxic","Ineffective"])
    st.markdown(f"""
    <div class="hm-wrap">
      <div class="hm-eye">Drug × Gene Risk Matrix</div>
      <div class="hm-grid" style="grid-template-columns:72px repeat({n},1fr);">{hdrs}{rows}</div>
      <div class="hm-legend">{legend}</div>
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# CHROMOSOME VISUALIZATION
# ══════════════════════════════════════════════════════════════════════════════

def render_chromosome(all_outputs, parsed):
    detected = set(parsed.get("detected_genes",[]))
    rmap     = {o["pharmacogenomic_profile"]["primary_gene"]: o for o in all_outputs}
    rows = ""
    for gene, info in CHROM_INFO.items():
        ch  = info["chrom"]; pos = info["pos_mb"]
        pct = (pos / CHROM_LEN.get(ch, 200)) * 100
        if gene in rmap:
            rl    = rmap[gene]["risk_assessment"]["risk_label"]
            mc    = RISK_CONFIG.get(rl, RISK_CONFIG["Unknown"])["dot"]
            pulse = "animation:glow-pulse 2s infinite;"
        elif gene in detected:
            mc = "#6b7280"; pulse = ""
        else:
            mc = "#1e1e2e"; pulse = ""
        rows += f"""<div class="chrom-row">
          <div class="chrom-lbl">{ch}</div>
          <div class="chrom-bar">
            <div class="chrom-body"></div>
            <div class="chrom-marker" style="left:{pct}%;background:{mc};{pulse}box-shadow:0 0 8px {mc}88;"></div>
          </div>
          <div class="chrom-gene-lbl">{gene}</div>
          <div class="chrom-band">{info['band']}</div>
        </div>"""
    st.markdown(f"""
    <div class="chrom-wrap">
      <div class="chrom-eye">Variant Chromosome Locations</div>
      {rows}
      <div style="font-family:'DM Mono',monospace;font-size:.53rem;color:#18181e;margin-top:.5rem;">Glowing markers = variants detected</div>
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# POPULATION FREQUENCY
# ══════════════════════════════════════════════════════════════════════════════

def render_pop_freq(gene, phenotype):
    freq = POPULATION_FREQ.get(gene, {})
    if not freq: return
    rows = ""
    for ph, pct in sorted(freq.items(), key=lambda x: -x[1]):
        is_you = (ph == phenotype)
        bc  = "#7c6aff" if is_you else "#1e1e2e"
        tc  = "#a78bfa" if is_you else "#30303c"
        you = '<span class="pop-you">← YOU</span>' if is_you else ""
        rows += f"""<div class="pop-row">
          <div class="pop-ph" style="color:{tc};">{ph}</div>
          <div class="pop-track"><div class="pop-fill" style="width:{min(pct,100)}%;background:{bc};"></div></div>
          <div class="pop-pct" style="color:{tc};">{pct}%{you}</div>
        </div>"""
    st.markdown(f"""
    <div class="pop-wrap">
      <div class="pop-eye">{gene} — Population Frequency</div>
      {rows}
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# INTERACTION MATRIX
# ══════════════════════════════════════════════════════════════════════════════

def render_ix_matrix(all_outputs, ix_report):
    if not ix_report or len(all_outputs) < 2: return
    drugs = [o["drug"] for o in all_outputs]
    n = len(drugs)
    sev_map = {}
    for ix in ix_report.get("all_interactions",[]):
        involved = ix.get("drugs_involved",[])
        if len(involved) == 2:
            sev = ix.get("severity","none")
            k1, k2 = (involved[0],involved[1]), (involved[1],involved[0])
            sev_map[k1] = sev_map[k2] = sev
    MC = {
        "critical": {"bg":"#450a0a","text":"#f87171","border":"#991b1b"},
        "high":     {"bg":"#450a0a","text":"#f87171","border":"#991b1b"},
        "moderate": {"bg":"#451a03","text":"#fbbf24","border":"#92400e"},
        "low":      {"bg":"#421a02","text":"#fcd34d","border":"#78350f"},
        "none":     {"bg":"#052e16","text":"#4ade80","border":"#166534"},
        "diag":     {"bg":"#0d0d14","text":"#1e1e2e","border":"#12121c"},
    }
    hdrs = '<div class="ix-head"></div>'
    for d in drugs:
        hdrs += f'<div class="ix-head">{d[:5]}</div>'
    grid = ""
    for i, d1 in enumerate(drugs):
        grid += f'<div class="ix-head" style="justify-content:flex-end;padding-right:4px;">{d1[:6]}</div>'
        for j, d2 in enumerate(drugs):
            if i == j:
                mc = MC["diag"]
                grid += f'<div class="ix-cell" style="background:{mc["bg"]};border-color:{mc["border"]};color:{mc["text"]};">—</div>'
            else:
                sev = sev_map.get((d1,d2),"none")
                mc  = MC.get(sev, MC["none"])
                lbl = sev.upper() if sev != "none" else "OK"
                grid += f'<div class="ix-cell" style="background:{mc["bg"]};border-color:{mc["border"]};color:{mc["text"]};">{lbl}</div>'
    st.markdown(f"""
    <div style="margin-bottom:1rem;">
      <div class="section-lbl">Drug Interaction Matrix</div>
      <div class="ix-matrix-grid" style="grid-template-columns:68px repeat({n},1fr);gap:3px;">
        {hdrs}{grid}
      </div>
    </div>""", unsafe_allow_html=True)
    shown = set()
    for ix in ix_report.get("all_interactions",[]):
        involved = ix.get("drugs_involved",[])
        if len(involved) == 2:
            key = tuple(sorted(involved))
            if key not in shown:
                shown.add(key)
                sev   = ix.get("severity","low")
                sp    = SEV_PALETTE.get(sev, SEV_PALETTE["low"])
                with st.expander(f"{' + '.join(involved)} — {sev.upper()}"):
                    mech = ix.get("mechanism", ix.get("message",""))
                    rec  = ix.get("recommendation","")
                    if mech:
                        st.markdown(f'<div style="font-size:.9rem;color:#909090;line-height:1.75;margin-bottom:.5rem;">{mech}</div>', unsafe_allow_html=True)
                    if rec:
                        st.markdown(f'<div style="font-family:DM Mono,monospace;font-size:.73rem;color:{sp["text"]};margin-top:.5rem;">→ {rec}</div>', unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# AI PATIENT NARRATIVE
# ══════════════════════════════════════════════════════════════════════════════

def render_ai_narrative(all_outputs, parsed, pid, groq_key, skip_llm):
    results_for_nar = [
        {"drug": o["drug"], "primary_gene": o["pharmacogenomic_profile"]["primary_gene"],
         "phenotype": o["pharmacogenomic_profile"]["phenotype"],
         "risk_label": o["risk_assessment"]["risk_label"],
         "severity": o["risk_assessment"]["severity"]}
        for o in all_outputs
    ]
    with st.spinner("Generating AI clinical summary…"):
        narrative = generate_patient_narrative(pid, results_for_nar, parsed, groq_key, skip_llm)
    badge = "static-template" if (skip_llm or not groq_key) else "llama-3.3-70b"
    st.markdown(f"""
    <div class="ai-narrative">
      <div class="ai-nar-head">
        <span class="ai-nar-badge">{badge}</span>
        <span class="ai-nar-title">AI Clinical Summary</span>
      </div>
      <div class="ai-nar-text">{narrative}</div>
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# CLINICAL NOTE GENERATOR
# ══════════════════════════════════════════════════════════════════════════════

def render_clinical_note(all_outputs, pid):
    lines = [f"Patient {pid} — Pharmacogenomic Analysis — {datetime.utcnow().strftime('%Y-%m-%d')}", ""]
    for o in all_outputs:
        gene = o["pharmacogenomic_profile"]["primary_gene"]
        dip  = o["pharmacogenomic_profile"]["diplotype"]
        ph   = o["pharmacogenomic_profile"]["phenotype"]
        drug = o["drug"]
        rl   = o["risk_assessment"]["risk_label"]
        rec  = o["clinical_recommendation"]["dosing_recommendation"]
        alts = o["clinical_recommendation"].get("alternative_drugs", [])
        lines.append(f"Patient carries {gene} {dip} ({ph} phenotype), therefore {drug.lower()} is predicted {rl.lower()}.")
        lines.append(f"CPIC Recommendation: {rec}")
        if alts:
            lines.append(f"Alternatives: {', '.join(alts)}")
        lines.append("")
    lines.append("Generated by PharmaGuard v7.0 · CPIC Level A evidence · cpicpgx.org")
    note = "\n".join(lines)
    st.markdown('<div class="section-lbl">One-Click Clinical Note</div>', unsafe_allow_html=True)
    st.markdown(f'<div class="note-box"><pre>{note}</pre></div>', unsafe_allow_html=True)
    st.download_button("⬇ Copy Clinical Note", data=note,
                       file_name=f"clinical_note_{pid}.txt", mime="text/plain",
                       key=f"note_{pid}")


# ══════════════════════════════════════════════════════════════════════════════
# BEFORE / AFTER SCENARIO
# ══════════════════════════════════════════════════════════════════════════════

def render_before_after(all_outputs):
    dangerous = [o for o in all_outputs if o["risk_assessment"]["risk_label"] in ("Toxic","Ineffective")]
    if not dangerous: return
    o    = dangerous[0]
    drug = o["drug"]; rl = o["risk_assessment"]["risk_label"]
    alts = o["clinical_recommendation"].get("alternative_drugs",[])
    alt  = alts[0] if alts else "Alternative medication"
    gene = o["pharmacogenomic_profile"]["primary_gene"]
    ph   = o["pharmacogenomic_profile"]["phenotype"]
    BEFORE = {"Toxic": f"Standard {drug.lower()} dose → toxic plasma levels → life-threatening outcome",
              "Ineffective": f"Standard {drug.lower()} dose → zero therapeutic effect → treatment failure"}
    AFTER  = f"{alt} prescribed → appropriate dosing → safe therapeutic outcome"
    st.markdown(f"""
    <div class="section-lbl" style="margin-top:1rem;">Before / After — PGx Impact</div>
    <div class="ba-wrap">
      <div class="ba-side" style="background:#1a0505;">
        <div class="ba-label" style="color:#f87171;">⛔ Without PharmaGuard</div>
        <div class="ba-scenario" style="color:#fca5a5;">{drug.title()} — Standard Protocol</div>
        <div class="ba-outcome" style="color:#f87171;">{BEFORE.get(rl, "Risk undetected")}</div>
        <div style="margin-top:.75rem;font-family:'DM Mono',monospace;font-size:.65rem;color:#7f1d1d;">
          {gene} {ph} phenotype undetected
        </div>
      </div>
      <div class="ba-side" style="background:#012010;">
        <div class="ba-label" style="color:#4ade80;">✅ With PharmaGuard</div>
        <div class="ba-scenario" style="color:#86efac;">{alt} — PGx-Guided Protocol</div>
        <div class="ba-outcome" style="color:#4ade80;">{AFTER}</div>
        <div style="margin-top:.75rem;font-family:'DM Mono',monospace;font-size:.65rem;color:#166534;">
          {gene} {ph} phenotype identified → therapy optimised
        </div>
      </div>
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# PRESCRIPTION SAFETY CHECKER
# ══════════════════════════════════════════════════════════════════════════════

def render_rx_checker(all_outputs):
    st.markdown('<div class="section-lbl" style="margin-top:1.5rem;">Prescription Safety Checker</div>', unsafe_allow_html=True)
    drug_names = [o["drug"] for o in all_outputs]
    rmap       = {o["drug"]: o for o in all_outputs}
    col1, col2 = st.columns([2, 1])
    with col1:
        selected = st.selectbox("Select drug to check", drug_names,
                                format_func=lambda x: f"{x.title()} ({GENE_DRUG_MAP.get(x,'')})",
                                key="rx_checker_drug", label_visibility="collapsed")
    with col2:
        check = st.button("Check Safety", key="rx_check_btn")
    if check and selected in rmap:
        o   = rmap[selected]
        rl  = o["risk_assessment"]["risk_label"]
        sev = o["risk_assessment"]["severity"]
        rec = o["clinical_recommendation"]["dosing_recommendation"]
        gene= o["pharmacogenomic_profile"]["primary_gene"]
        ph  = o["pharmacogenomic_profile"]["phenotype"]
        rc  = RISK_CONFIG.get(rl, RISK_CONFIG["Unknown"])
        sp  = SEV_PALETTE.get(sev, SEV_PALETTE["none"])
        icon = {"Safe":"✅ SAFE TO PRESCRIBE","Adjust Dosage":"⚠️ PRESCRIBE WITH DOSE ADJUSTMENT",
                "Toxic":"🚨 DO NOT PRESCRIBE — HIGH TOXICITY RISK",
                "Ineffective":"❌ DO NOT PRESCRIBE — DRUG WILL BE INEFFECTIVE"}.get(rl,"❓ INSUFFICIENT DATA")
        st.markdown(f"""
        <div class="rx-result" style="background:{rc['dark_bg']};border-color:{rc['dark_border']};">
          <div style="font-family:'DM Mono',monospace;font-size:.8rem;font-weight:700;
            color:{rc['dark_text']};margin-bottom:.5rem;letter-spacing:.05em;">{icon}</div>
          <div style="font-size:.875rem;color:{rc['dark_text']};line-height:1.7;margin-bottom:.5rem;">
            {gene} {ph} phenotype detected. {rec}
          </div>
          <div style="font-family:'DM Mono',monospace;font-size:.6rem;color:{sp['text']};letter-spacing:.06em;text-transform:uppercase;">
            Severity: {sev} · Confidence: {o["risk_assessment"]["confidence_score"]:.0%} · CPIC Level A
          </div>
        </div>""", unsafe_allow_html=True)
    elif not check:
        st.markdown(f"""
        <div style="background:#08080e;border:1px dashed #14141e;border-radius:10px;
            padding:1rem 1.5rem;font-family:'DM Mono',monospace;font-size:.7rem;color:#20202e;">
          Select a drug above and click Check Safety to instantly validate this prescription
          against the patient's genotype.
        </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# PATIENT PLAIN ENGLISH MODE
# ══════════════════════════════════════════════════════════════════════════════

def render_patient_mode(all_outputs):
    bad  = any(o["risk_assessment"]["risk_label"] in ("Toxic","Ineffective") for o in all_outputs)
    hc   = "#ef4444" if bad else "#22c55e"
    icon = "🚨" if bad else "✅"
    msg  = "Important: Some medications need urgent attention" if bad else "Good news — your medications look safe"
    st.markdown(f"""
    <div style="background:#0f0a14;border:2px solid {hc}44;border-radius:14px;padding:1.5rem;margin-bottom:1.5rem;">
      <div style="font-size:1.2rem;font-weight:700;color:{hc};margin-bottom:.4rem;">{icon} {msg}</div>
      <div style="font-size:.9rem;color:#40405a;line-height:1.7;">
        This report analysed your DNA to see how your body handles certain medicines.
        Everyone's body is different — your genes affect how medicines work for you.
      </div>
    </div>""", unsafe_allow_html=True)
    for o in all_outputs:
        drug  = o["drug"]; rl = o["risk_assessment"]["risk_label"]
        gene  = o["pharmacogenomic_profile"]["primary_gene"]
        ph    = o["pharmacogenomic_profile"]["phenotype"]
        alts  = o["clinical_recommendation"].get("alternative_drugs",[])
        phplain = PLAIN_ENGLISH_PHENOTYPE.get(ph, ph)
        explain = PLAIN_ENGLISH_RISK.get((drug, ph), "")
        VERDICT = {
            "Safe":          "✅ This medicine is likely safe for you",
            "Adjust Dosage": "⚠️ You may need a different dose of this medicine",
            "Toxic":         "🚨 This medicine could be harmful to you",
            "Ineffective":   "❌ This medicine likely won't work for you",
        }
        verdict = VERDICT.get(rl, rl)
        bc = {"Safe":"#22c55e","Adjust Dosage":"#f59e0b","Toxic":"#ef4444","Ineffective":"#8b5cf6"}.get(rl,"#6b7280")
        action = ""
        if rl in ("Toxic","Ineffective"):
            alt_text = f"They may suggest: {', '.join(alts[:3])}" if alts else "Ask about alternative medications."
            action = f'<div class="patient-action"><span style="font-size:1.1rem;">💊</span><div style="font-size:.875rem;color:#e0e0e0;line-height:1.65;"><strong>Talk to your doctor before taking {drug.title()}.</strong><br>{alt_text}</div></div>'
        elif rl == "Adjust Dosage":
            action = f'<div class="patient-action"><span style="font-size:1.1rem;">📋</span><div style="font-size:.875rem;color:#e0e0e0;line-height:1.65;"><strong>Tell your doctor about this result before starting {drug.title()}.</strong><br>You may need a different dose than usually prescribed.</div></div>'
        st.markdown(f"""
        <div class="patient-card" style="border-color:{bc}33;">
          <div style="font-size:1.25rem;font-weight:700;margin-bottom:.2rem;">{drug.title()}</div>
          <div class="patient-verdict" style="color:{bc}dd;">{verdict}</div>
          <div style="font-family:'DM Mono',monospace;font-size:.6rem;letter-spacing:.06em;color:#30303c;margin-bottom:.5rem;">{gene} · {phplain}</div>
          {f'<div class="patient-plain" style="color:#606070;">{explain}</div>' if explain else ''}
          {action}
        </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# SIDE-BY-SIDE DRUG COMPARISON
# ══════════════════════════════════════════════════════════════════════════════

def render_drug_comparison(all_outputs):
    dangerous = [o for o in all_outputs if o["risk_assessment"]["risk_label"] in ("Toxic","Ineffective")]
    if not dangerous: return
    for o in dangerous[:2]:
        drug = o["drug"]
        gene = o["pharmacogenomic_profile"]["primary_gene"]
        ph   = o["pharmacogenomic_profile"]["phenotype"]
        alts = o["clinical_recommendation"].get("alternative_drugs",[])
        if not alts: continue
        alt  = alts[0]
        rc_bad  = RISK_CONFIG.get(o["risk_assessment"]["risk_label"], RISK_CONFIG["Unknown"])
        rc_good = RISK_CONFIG["Safe"]
        st.markdown(f"""
        <div style="margin-bottom:1rem;">
          <div class="section-lbl">{drug.title()} vs {alt} — For This Patient</div>
          <div style="display:grid;grid-template-columns:1fr 1fr;gap:1px;background:#14141e;border-radius:12px;overflow:hidden;">
            <div style="background:{rc_bad['dark_bg']};padding:1.25rem 1.5rem;">
              <div style="font-family:'DM Mono',monospace;font-size:.58rem;letter-spacing:.1em;text-transform:uppercase;color:{rc_bad['dark_text']};margin-bottom:.5rem;">Current Prescription</div>
              <div style="font-size:1.1rem;font-weight:700;color:{rc_bad['dark_text']};margin-bottom:.3rem;">{drug.title()}</div>
              <div style="font-family:'DM Mono',monospace;font-size:.68rem;color:{rc_bad['dark_text']};margin-bottom:.5rem;">{rc_bad['emoji']} {o["risk_assessment"]["risk_label"]} · {gene} {ph}</div>
              <div style="font-size:.875rem;color:{rc_bad['dark_text']};opacity:.8;line-height:1.65;">{o["clinical_recommendation"]["dosing_recommendation"][:150]}…</div>
            </div>
            <div style="background:{rc_good['dark_bg']};padding:1.25rem 1.5rem;">
              <div style="font-family:'DM Mono',monospace;font-size:.58rem;letter-spacing:.1em;text-transform:uppercase;color:{rc_good['dark_text']};margin-bottom:.5rem;">Recommended Alternative</div>
              <div style="font-size:1.1rem;font-weight:700;color:{rc_good['dark_text']};margin-bottom:.3rem;">{alt}</div>
              <div style="font-family:'DM Mono',monospace;font-size:.68rem;color:{rc_good['dark_text']};margin-bottom:.5rem;">✅ Safe · PGx-guided selection</div>
              <div style="font-size:.875rem;color:{rc_good['dark_text']};opacity:.8;line-height:1.65;">Selected based on {gene} {ph} phenotype per CPIC Level A guidelines.</div>
            </div>
          </div>
        </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# MASTER RESULTS RENDERER
# ══════════════════════════════════════════════════════════════════════════════

def render_results(all_outputs, parsed, ix_report, pdf_bytes, pid,
                   patient_mode=False, groq_key="", skip_llm=False):

    # 1. Risk Command Center Banner
    render_risk_banner(all_outputs, parsed)

    # 2. Emergency Alerts
    render_emergency_alerts(all_outputs)

    # 3. Gene Activity Heatmap
    render_gene_heatmap(all_outputs)

    # 4. Download buttons row
    dc1, dc2, dc3 = st.columns(3)
    with dc1:
        st.download_button("⬇ All JSON", data=json.dumps(all_outputs, indent=2),
                           file_name=f"pharmaguard_{pid}.json", mime="application/json",
                           use_container_width=True, key=f"dlall_{pid}")
    with dc2:
        if pdf_bytes:
            st.download_button("⬇ PDF Report", data=pdf_bytes,
                               file_name=f"pharmaguard_{pid}.pdf", mime="application/pdf",
                               use_container_width=True, key=f"dlpdf_{pid}")
    with dc3:
        if ix_report and ix_report.get("interactions_found"):
            st.download_button("⬇ Interactions JSON", data=json.dumps(ix_report, indent=2),
                               file_name=f"pharmaguard_{pid}_ix.json", mime="application/json",
                               use_container_width=True, key=f"dlix_{pid}")

    st.markdown("<div style='height:.75rem'></div>", unsafe_allow_html=True)

    if patient_mode:
        render_patient_mode(all_outputs)
        return

    # ── DOCTOR MODE FEATURES ─────────────────────────────────

    # 5. Drug Risk Comparison Table
    st.markdown('<div class="section-lbl">Drug Risk Comparison</div>', unsafe_allow_html=True)
    render_drug_table(all_outputs, pid)

    # 6. Polygenic Risk Score
    render_pgx_score(all_outputs)

    # 7. Drug x Gene Heatmap + Chromosome side by side
    c1, c2 = st.columns([1.4, 1], gap="large")
    with c1: render_heatmap(all_outputs)
    with c2: render_chromosome(all_outputs, parsed)

    # 8. Drug Interaction Matrix
    if ix_report and len(all_outputs) >= 2:
        render_ix_matrix(all_outputs, ix_report)

    # 9. AI Narrative
    st.markdown('<div class="section-lbl" style="margin-top:1rem;">AI Clinical Summary</div>', unsafe_allow_html=True)
    render_ai_narrative(all_outputs, parsed, pid, groq_key, skip_llm)

    # 10. Before / After + Side-by-side comparison
    render_before_after(all_outputs)
    render_drug_comparison(all_outputs)

    # 11. Prescription Safety Checker
    render_rx_checker(all_outputs)

    # 12. Clinical Note Generator
    render_clinical_note(all_outputs, pid)

    # 13. Individual Drug Cards
    st.markdown('<div class="section-lbl" style="margin-top:1.5rem;">Individual Drug Analysis</div>', unsafe_allow_html=True)

    for output in all_outputs:
        rl   = output["risk_assessment"]["risk_label"]
        dn   = output["drug"]
        sev  = output["risk_assessment"]["severity"]
        conf = output["risk_assessment"]["confidence_score"]
        gene = output["pharmacogenomic_profile"]["primary_gene"]
        dip  = output["pharmacogenomic_profile"]["diplotype"]
        ph   = output["pharmacogenomic_profile"]["phenotype"]
        var  = output["pharmacogenomic_profile"]["detected_variants"]
        rec  = output["clinical_recommendation"]["dosing_recommendation"]
        alts = output["clinical_recommendation"].get("alternative_drugs",[])
        mon  = output["clinical_recommendation"].get("monitoring_required","")
        exp  = output["llm_generated_explanation"]
        rc   = RISK_CONFIG.get(rl, RISK_CONFIG["Unknown"])
        sp   = SEV_PALETTE.get(sev, SEV_PALETTE["none"])
        cpic_lv = output.get("pharmacogenomic_profile",{}).get("cpic_evidence_level","Level A")

        st.markdown(f"""
        <div class="rcard">
          <div class="rcard-top">
            <div class="rcard-left">
              <div class="rcard-dot" style="background:{rc['dot']};box-shadow:0 0 8px {rc['dot']}88;"></div>
              <div>
                <div class="rcard-name">{dn.title()}
                  <span class="cpic">CPIC {cpic_lv}</span>
                </div>
                <div class="rcard-meta">{gene} · {dip} · {ph}</div>
              </div>
            </div>
            <span class="rcard-badge" style="color:{rc['dark_text']};border-color:{rc['dark_border']};background:{rc['dark_bg']};">{rc['emoji']} {rl}</span>
          </div>
          <div class="rcard-body">
            <div class="mc-row">
              <div class="mc-cell"><div class="mc-key">Phenotype</div><div class="mc-val" style="color:{rc['dot']};">{ph}</div></div>
              <div class="mc-cell"><div class="mc-key">Severity</div><div class="mc-val" style="color:{sp['text']};">{sev.title()}</div></div>
              <div class="mc-cell"><div class="mc-key">Confidence</div><div class="mc-val">{conf:.0%}</div></div>
              <div class="mc-cell"><div class="mc-key">Variants</div><div class="mc-val">{len(var)}</div></div>
            </div>""", unsafe_allow_html=True)

        # Confidence bars
        dq = min(1.0, len(var) / 3.0)
        st.markdown(f"""
        <div class="conf-row">
          <div class="conf-item">
            <div class="conf-lbl"><span>Prediction Confidence</span><span style="color:{rc['dot']};font-weight:600;">{conf:.0%}</span></div>
            <div class="conf-track"><div class="conf-fill" style="width:{conf*100:.1f}%;background:{rc['dot']};box-shadow:0 0 6px {rc['dot']}55;"></div></div>
          </div>
          <div class="conf-item">
            <div class="conf-lbl"><span>Data Quality</span><span style="color:#30303c;">{len(var)} variant{"s" if len(var)!=1 else ""}</span></div>
            <div class="conf-track"><div class="conf-fill" style="width:{dq*100:.1f}%;background:#40405a;"></div></div>
          </div>
        </div>""", unsafe_allow_html=True)

        # Variants table
        if var:
            rows_html = ""
            for v in var:
                fc = func_cls(v.get("functional_status",""))
                fn = (v.get("functional_status") or "unknown").replace("_"," ").title()
                rows_html += f'<tr><td class="v-rsid">{v.get("rsid","—")}</td><td class="v-star">{v.get("star_allele","—")}</td><td class="{fc}">{fn}</td></tr>'
            st.markdown(f"""<hr class="h-rule">
            <div class="section-lbl">Detected Variants ({len(var)})</div>
            <table class="vtable"><thead><tr><th>rsID</th><th>Star Allele</th><th>Function</th></tr></thead>
            <tbody>{rows_html}</tbody></table>""", unsafe_allow_html=True)

        # CPIC Recommendation
        st.markdown(f"""<hr class="h-rule">
        <div class="section-lbl">CPIC Recommendation</div>
        <div class="rec-box" style="background:{rc['dark_bg']};border-color:{rc['dark_border']};">
          <div class="rec-lbl" style="color:{rc['dark_text']};">CPIC Guideline — {dn}</div>
          <div class="rec-txt">{rec}</div>
        </div>""", unsafe_allow_html=True)

        if alts:
            chips = "".join(f'<span class="alt-chip">{a}</span>' for a in alts)
            st.markdown(f'<div class="section-lbl">Alternative Drugs</div><div class="alt-chips" style="margin-bottom:1rem;">{chips}</div>', unsafe_allow_html=True)

        if mon:
            st.markdown(f"""<div style="display:flex;gap:.75rem;align-items:flex-start;padding:.875rem;
                background:#08080e;border:1px solid #12121c;border-radius:8px;margin-bottom:1rem;">
              <span>🔬</span>
              <div><div class="section-lbl" style="margin-bottom:.2rem;">Monitoring Protocol</div>
              <div style="font-size:.875rem;color:#606070;line-height:1.65;">{mon}</div></div>
            </div>""", unsafe_allow_html=True)

        st.markdown("</div></div>", unsafe_allow_html=True)

        # Population frequency
        render_pop_freq(gene, ph)

        # AI Explanation block
        if exp.get("summary"):
            is_static = "static" in exp.get("model_used","").lower()
            model     = exp.get("model_used","llama-3.3-70b")
            b_cls     = "ai-badge" if not is_static else "ai-badge"
            blocks = ""
            for lbl, key in [("Summary","summary"),("Biological Mechanism","biological_mechanism"),
                              ("Variant Significance","variant_significance"),("Clinical Implications","clinical_implications")]:
                if exp.get(key):
                    blocks += f'<div class="ai-section"><div class="ai-section-lbl">{lbl}</div><div class="ai-section-txt">{exp[key]}</div></div>'
            st.markdown(f"""<div class="ai-block">
              <div class="ai-block-head">
                <span class="{b_cls}">{model}</span>
                <span style="font-family:'DM Mono',monospace;font-size:.58rem;color:#28283a;">AI Explanation · {dn}</span>
              </div>{blocks}</div>""", unsafe_allow_html=True)

        # Per-drug JSON download
        with st.expander(f"Raw JSON — {dn}", expanded=False):
            c1, _ = st.columns([1,3])
            with c1:
                st.download_button(f"⬇ {dn} JSON", data=json.dumps(output, indent=2),
                                   file_name=f"pharmaguard_{pid}_{dn}.json", mime="application/json",
                                   key=f"dl_{pid}_{dn}", use_container_width=True)
            st.code(json.dumps(output, indent=2), language="json")

    # VCF parse details
    with st.expander("VCF Parse Details", expanded=False):
        p1, p2, p3 = st.columns(3)
        p1.metric("Total Variants", parsed["total_variants"])
        p2.metric("Genes Found",    len(parsed["detected_genes"]))
        p3.metric("Parse Errors",   len(parsed["parse_errors"]))


# ══════════════════════════════════════════════════════════════════════════════
# NAV
# ══════════════════════════════════════════════════════════════════════════════

st.markdown("""
<div class="pg-nav">
  <div class="pg-logo">Pharma<strong>Guard</strong> <em>v7</em></div>
  <div class="pg-tags">
    <span class="pg-tag pg-tag-default">CPIC Aligned</span>
    <span class="pg-tag pg-tag-default">RIFT 2026</span>
    <span class="pg-tag pg-tag-hot">★ Unified Edition</span>
  </div>
</div>""", unsafe_allow_html=True)

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown('<div style="padding:1rem 0 .5rem;font-family:DM Mono,monospace;font-size:.6rem;letter-spacing:.12em;text-transform:uppercase;color:#40405a;">Groq API Key</div>', unsafe_allow_html=True)
    groq_api_key = st.text_input("Groq API Key", value=os.environ.get("GROQ_API_KEY",""),
                                  type="password", label_visibility="collapsed", placeholder="gsk_…")
    st.markdown('<div style="font-family:DM Mono,monospace;font-size:.56rem;color:#20202e;padding-bottom:1rem;line-height:1.9;">Model: LLaMA 3.3 70B Versatile<br>Fallback: static expert templates<br>Test mode: instant (no API call)</div>', unsafe_allow_html=True)
    st.markdown('<hr style="border:none;border-top:1px solid #12121c;">', unsafe_allow_html=True)
    st.markdown('<div style="padding:.75rem 0 .5rem;font-family:DM Mono,monospace;font-size:.6rem;letter-spacing:.12em;text-transform:uppercase;color:#40405a;">Gene → Drug Map</div>', unsafe_allow_html=True)
    for g, drug in [("CYP2D6","Codeine"),("CYP2C19","Clopidogrel"),("CYP2C9","Warfarin"),
                    ("SLCO1B1","Simvastatin"),("TPMT","Azathioprine"),("DPYD","Fluorouracil")]:
        st.markdown(f'<div style="font-family:DM Mono,monospace;font-size:.66rem;padding:3px 0;color:#25253a;">{g} <span style="color:#30303c;">→</span> {drug}</div>', unsafe_allow_html=True)

tab1, tab2 = st.tabs(["Analysis", "Test Suite"])


# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

with tab1:

    # Persona buttons
    st.markdown('<div class="section-lbl">Quick Demo — Select Patient Persona</div>', unsafe_allow_html=True)
    p_cols = st.columns(4)
    for i, (key, p) in enumerate(PERSONAS.items()):
        with p_cols[i]:
            st.markdown(f"""
            <div class="persona-card" style="background:{p['bg']}22;border-color:{p['border']};">
              <div class="persona-card-label" style="color:{p['color']};">{p['label']}</div>
              <div class="persona-card-desc" style="color:{p['color']};">{p['desc']}</div>
            </div>""", unsafe_allow_html=True)
            if st.button("Load", key=f"persona_{key}", use_container_width=True):
                st.session_state["persona_file"]  = p["file"]
                st.session_state["persona_drugs"] = p["drugs"]

    st.markdown("<div style='height:1.5rem'></div>", unsafe_allow_html=True)

    # Steps bar
    st.markdown("""<div class="steps-bar">
      <div class="step active"><div class="step-n">01</div><div class="step-l">Upload VCF</div></div>
      <div class="step active"><div class="step-n">02</div><div class="step-l">Select Drugs</div></div>
      <div class="step"><div class="step-n">03</div><div class="step-l">Run Analysis</div></div>
      <div class="step"><div class="step-n">04</div><div class="step-l">Review Results</div></div>
    </div>""", unsafe_allow_html=True)

    col_l, col_r = st.columns([1.3, 1], gap="large")

    with col_l:
        st.markdown('<div class="section-lbl">Genomic Data</div>', unsafe_allow_html=True)
        uploaded = st.file_uploader("Upload VCF", type=["vcf"])
        if uploaded:
            sz = uploaded.size / (1024*1024)
            if sz > 5:
                st.error(f"File too large ({sz:.1f} MB). Max 5 MB."); uploaded = None
            else:
                peek = uploaded.read(400).decode("utf-8", errors="replace"); uploaded.seek(0)
                if "##fileformat=VCF" not in peek and "#CHROM" not in peek:
                    st.error("Invalid VCF file."); uploaded = None
                else:
                    st.success(f"✓  {uploaded.name}  ·  {sz:.2f} MB")
        st.markdown('<div class="section-lbl" style="margin-top:1rem;">Or select a test scenario</div>', unsafe_allow_html=True)
        scenario_opts = {
            "None": None,
            "Mixed Variants (Standard)": "sample.vcf",
            "UltraRapid Metabolizer": "test_ultrarapid_metabolizer.vcf",
            "All Normal Wild-type": "test_all_normal_wildtype.vcf",
            "Worst Case — All Poor Metabolizers": "test_worst_case_all_pm.vcf",
            "Patient A — Critical Risk": "patient_a_critical.vcf",
            "Patient B — Warfarin PM": "patient_b_warfarin.vcf",
            "Patient C — Clopidogrel PM": "patient_c_interaction.vcf",
            "Patient D — All Safe": "patient_d_safe.vcf",
        }
        chosen_label = st.selectbox("Scenario", list(scenario_opts.keys()), label_visibility="collapsed")
        chosen_file  = scenario_opts[chosen_label]

    with col_r:
        st.markdown('<div class="section-lbl">Medications</div>', unsafe_allow_html=True)
        default_drugs = st.session_state.get("persona_drugs", ["CLOPIDOGREL"])
        drugs_selected = st.multiselect(
            "Select drugs", options=ALL_DRUGS, default=default_drugs,
            format_func=lambda x: f"{x.title()}  ({GENE_DRUG_MAP.get(x,'')})",
            label_visibility="collapsed"
        )
        st.markdown('<div class="section-lbl" style="margin-top:.75rem;">Custom drugs (comma-separated)</div>', unsafe_allow_html=True)
        custom_drugs = st.text_input("Custom", placeholder="CODEINE, WARFARIN…", label_visibility="collapsed")
        st.markdown('<div class="section-lbl" style="margin-top:.75rem;">Patient ID</div>', unsafe_allow_html=True)
        pid_input = st.text_input("Patient ID", placeholder="Auto-generated if blank", label_visibility="collapsed")
        c1, c2 = st.columns(2)
        with c1: do_ix  = st.checkbox("Check drug interactions", value=True)
        with c2: do_pdf = st.checkbox("Generate PDF report", value=True)
        st.markdown('<div class="section-lbl" style="margin-top:.75rem;">View mode</div>', unsafe_allow_html=True)
        patient_mode = st.toggle("Patient Plain-English Mode",
                                  help="Converts clinical jargon to language any patient understands")

    st.markdown("<div style='height:.5rem'></div>", unsafe_allow_html=True)
    run_btn = st.button("Run Analysis →", use_container_width=True)

    if run_btn:
        all_drugs = list(drugs_selected)
        if custom_drugs.strip():
            all_drugs += [d.strip().upper() for d in custom_drugs.split(",") if d.strip()]
        all_drugs = list(set(all_drugs))
        if not all_drugs:
            st.error("Please select at least one drug."); st.stop()

        vcf_content = None
        if "persona_file" in st.session_state and not uploaded and not chosen_file:
            vcf_content = load_vcf_file(st.session_state["persona_file"])
        elif uploaded:
            vcf_content = uploaded.read().decode("utf-8", errors="replace")
        elif chosen_file:
            vcf_content = load_vcf_file(chosen_file)
        else:
            st.error("Please upload a VCF or select a test scenario."); st.stop()

        pid = pid_input.strip() or f"PG-{str(uuid.uuid4())[:8].upper()}"

        st.markdown(f"""
        <div style="display:flex;align-items:baseline;gap:1rem;margin:2.5rem 0 1.5rem;
            padding-bottom:1rem;border-bottom:1px solid #14141e;">
          <div style="font-family:'Fraunces',serif;font-size:1.75rem;font-weight:300;">Results</div>
          <div style="font-family:'DM Mono',monospace;font-size:.68rem;color:#30303c;">{pid}</div>
        </div>""", unsafe_allow_html=True)

        with st.spinner("Analysing genomic data…"):
            parsed, risk_results, all_outputs, ix_report, pdf_bytes = run_pipeline(
                vcf_content, all_drugs, pid, groq_api_key, do_ix, do_pdf)

        render_results(all_outputs, parsed, ix_report, pdf_bytes, pid,
                       patient_mode=patient_mode, groq_key=groq_api_key,
                       skip_llm=(not groq_api_key))
    else:
        st.markdown("""<div class="empty">
          <div class="empty-icon">🧬</div>
          <div class="empty-title">Ready for analysis</div>
          <div class="empty-hint">
            Click a Patient Persona above for an instant demo<br>
            Or upload a VCF · select drugs · click Run<br><br>
            CYP2D6 · CYP2C19 · CYP2C9 · SLCO1B1 · TPMT · DPYD<br><br>
            v7.0 → Risk Banner · Emergency Alerts · Gene Heatmap · Drug Table<br>
            Interaction Matrix · AI Narrative · CPIC Badges · Rx Checker<br>
            Before/After · Drug Comparison · Clinical Note · Patient Mode
          </div>
        </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — TEST SUITE (parallel)
# ══════════════════════════════════════════════════════════════════════════════

RISK_DOT = {"Safe":"#22c55e","Adjust Dosage":"#f59e0b","Toxic":"#ef4444","Ineffective":"#8b5cf6","Unknown":"#94a3b8"}

with tab2:
    st.markdown("""
    <div style="margin-bottom:2rem;">
      <div style="font-family:'Fraunces',serif;font-size:1.75rem;font-weight:300;margin-bottom:.4rem;">Test Suite</div>
      <div style="font-family:'DM Mono',monospace;font-size:.62rem;color:#30303c;letter-spacing:.06em;">
        4 scenarios · parallel execution · pass/fail validation per drug
      </div>
    </div>""", unsafe_allow_html=True)

    tc = st.columns(4)
    for i, sc in enumerate(TEST_SUITE):
        with tc[i]:
            rh = "".join(
                f'<div class="test-row"><span class="test-drug-lbl">{d[:9]}</span>'
                f'<span style="font-family:DM Mono,monospace;font-size:.68rem;font-weight:600;color:{RISK_DOT.get(r,"#6b7280")};">{r}</span></div>'
                for d, r in list(sc["expected"].items())[:4]
            )
            st.markdown(f'<div class="test-card"><div class="test-name">{sc["name"]}</div>'
                        f'<div class="test-desc">{sc["desc"]}</div>{rh}</div>', unsafe_allow_html=True)

    st.markdown("<div style='height:1.25rem'></div>", unsafe_allow_html=True)
    tb1, tb2 = st.columns([3,1])
    with tb1: use_llm = st.checkbox("Include LLM Explanations (requires Groq API)", value=False)
    with tb2: run_all = st.button("Run All 4 Tests →", use_container_width=True)

    if run_all:
        def run_one_test(sc):
            vcf = load_vcf_file(sc["file"])
            pid = f"TEST-{sc['name'][:6].replace(' ','').upper()}"
            pv, _, ao, _, _ = run_pipeline(vcf, sc["drugs"], pid,
                groq_api_key if use_llm else "", run_ix=False, gen_pdf=False, skip_llm=not use_llm)
            rows, ok = [], True
            for out in ao:
                drug = out["drug"]; got = out["risk_assessment"]["risk_label"]
                exp  = sc["expected"].get(drug,""); passed = (got==exp) if exp else True
                ph   = out["pharmacogenomic_profile"]["phenotype"]
                dp   = out["pharmacogenomic_profile"]["diplotype"]
                rows.append((drug, got, exp, passed, ph, dp))
                if not passed: ok = False
            return {"name": sc["name"], "pass": ok, "rows": rows, "outputs": ao, "file": sc["file"]}

        prog = st.empty(); prog.info("Running all 4 scenarios in parallel…")
        results = [None]*4
        with ThreadPoolExecutor(max_workers=4) as ex:
            futs = {ex.submit(run_one_test, sc): i for i, sc in enumerate(TEST_SUITE)}
            done = 0
            for f in as_completed(futs):
                results[futs[f]] = f.result(); done += 1
                prog.info(f"Completed {done}/4 scenarios…")
        prog.empty()

        passed = sum(1 for r in results if r["pass"])
        failed = 4 - passed
        oc = "#22c55e" if failed==0 else "#f59e0b"
        ob = "#052e16" if failed==0 else "#451a03"
        od = "#166534" if failed==0 else "#92400e"

        st.markdown(f"""
        <div style="background:{ob};border:1px solid {od};border-radius:10px;
            padding:1.25rem 1.5rem;margin:1.25rem 0;display:flex;align-items:center;justify-content:space-between;">
          <div style="font-family:'Fraunces',serif;font-size:1.4rem;font-weight:300;color:{oc};">
            {'All tests passed ✓' if failed==0 else f'{passed}/4 tests passed'}
          </div>
          <div style="font-family:'DM Mono',monospace;font-size:.62rem;color:{oc};">
            {passed} passed · {failed} failed · {int(passed/4*100)}%
          </div>
        </div>""", unsafe_allow_html=True)

        for sr in results:
            sym = "✓ Pass" if sr["pass"] else "✗ Fail"
            with st.expander(f"{sym}  —  {sr['name']}", expanded=not sr["pass"]):
                table_html = ""
                for drug, got, exp, ok, ph, dp in sr["rows"]:
                    rc   = RISK_CONFIG.get(got, RISK_CONFIG["Unknown"])
                    oc2  = "#22c55e" if ok else "#ef4444"
                    ob2  = "#052e16" if ok else "#450a0a"
                    table_html += f"""<div style="display:grid;grid-template-columns:1fr 1.2fr 1.2fr 1.5fr 40px;
                        border-bottom:1px solid #10101a;padding:0 .5rem;">
                      <div style="font-family:DM Mono,monospace;font-size:.77rem;font-weight:700;
                        color:#e8e8f0;padding:.7rem .85rem;">{drug}</div>
                      <div style="font-family:DM Mono,monospace;font-size:.77rem;color:{rc['dot']};
                        padding:.7rem .85rem;display:flex;align-items:center;gap:5px;">
                        <span style="width:6px;height:6px;border-radius:50%;background:{rc['dot']};flex-shrink:0;"></span>{got}</div>
                      <div style="font-family:DM Mono,monospace;font-size:.77rem;color:#40405a;padding:.7rem .85rem;">{exp or '—'}</div>
                      <div style="font-family:DM Mono,monospace;font-size:.77rem;color:#40405a;padding:.7rem .85rem;">{dp} / {ph}</div>
                      <div style="font-family:DM Mono,monospace;font-size:.72rem;font-weight:700;
                        color:{oc2};background:{ob2};display:flex;align-items:center;justify-content:center;
                        border-radius:4px;margin:.5rem 0;">{'OK' if ok else 'FAIL'}</div>
                    </div>"""
                st.markdown(f"""
                <div style="border:1px solid #14141e;border-radius:10px;overflow:hidden;background:#0d0d14;margin-bottom:1rem;">
                  <div style="display:grid;grid-template-columns:1fr 1.2fr 1.2fr 1.5fr 40px;
                    background:#08080e;border-bottom:1px solid #14141e;padding:0 .5rem;">
                    {"".join(f'<div style=\"font-family:DM Mono,monospace;font-size:.6rem;letter-spacing:.1em;text-transform:uppercase;color:#20202e;padding:.7rem .85rem;\">{h}</div>' for h in ['Drug','Result','Expected','Diplotype / Phenotype',''])}
                  </div>{table_html}
                </div>""", unsafe_allow_html=True)
                d1, d2 = st.columns(2)
                with d1:
                    st.download_button("⬇ JSON", data=json.dumps(sr["outputs"],indent=2),
                        file_name=f"test_{sr['file'].replace('.vcf','')}.json", mime="application/json",
                        key=f"tsc_{sr['name'][:14]}", use_container_width=True)
                with d2:
                    st.download_button("⬇ VCF", data=load_vcf_file(sr["file"]),
                        file_name=sr["file"], mime="text/plain",
                        key=f"vcf_{sr['name'][:14]}", use_container_width=True)

        st.download_button("⬇ Download Full Test Suite JSON",
            data=json.dumps([{"scenario":s["name"],"pass":s["pass"],"results":s["outputs"]} for s in results], indent=2),
            file_name=f"pharmaguard_tests_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
            mime="application/json", use_container_width=True)
    else:
        st.markdown("""<div class="empty">
          <div class="empty-icon">▷</div>
          <div class="empty-title">One-click validation</div>
          <div class="empty-hint">
            4 clinical scenarios · ThreadPoolExecutor parallel execution<br>
            Mixed · UltraRapid · All Normal · Worst Case All PM<br>
            Pass/fail per drug with expected vs actual comparison
          </div>
        </div>""", unsafe_allow_html=True)