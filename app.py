"""
SurakshaRx v10.0 — Precision Clinical Pharmacogenomics
FIXES vs v9.2:
  - UI FIX: File uploader label now visible (was "collapsed"), preventing overlap with scenario selectbox
  - UI FIX: Added spacing between scenario selectbox and file uploader
  - UI FIX: Selectbox label explicitly set to "visible" for clear visual hierarchy
  - UI FIX: File size display corrected (was showing KB as MB)
  - TEST SUITE FIX: load_vcf() calls wrapped in try/except — falls back to get_sample_vcf()
    so tests run even when sample_data/ files are missing from deployment
  - PERSONA DEMO FIX: Same try/except fallback for persona file loading
"""

import streamlit as st
import json, uuid, os, re, io
import pandas as pd
from datetime import datetime, timezone
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

RISK_CFG = {
    "Safe":         {"color":"#4ADE80","bg":"rgba(22, 163, 74, 0.1)","border":"rgba(74, 222, 128, 0.2)","text":"#86EFAC","tag_bg":"rgba(22, 163, 74, 0.15)","tag_text":"#4ADE80","shape":"●","severity_dot":"#4ADE80"},
    "Adjust Dosage":{"color":"#FBBF24","bg":"rgba(217, 119, 6, 0.1)","border":"rgba(251, 191, 36, 0.2)","text":"#FDE047","tag_bg":"rgba(217, 119, 6, 0.15)","tag_text":"#FBBF24","shape":"▲","severity_dot":"#FBBF24"},
    "Toxic":        {"color":"#F87171","bg":"rgba(220, 38, 38, 0.1)","border":"rgba(248, 113, 113, 0.2)","text":"#FCA5A5","tag_bg":"rgba(220, 38, 38, 0.15)","tag_text":"#F87171","shape":"⬛","severity_dot":"#EF4444"},
    "Ineffective":  {"color":"#A78BFA","bg":"rgba(109, 40, 217, 0.1)","border":"rgba(167, 139, 250, 0.2)","text":"#C4B5FD","tag_bg":"rgba(109, 40, 217, 0.15)","tag_text":"#A78BFA","shape":"◆","severity_dot":"#A78BFA"},
    "Unknown":      {"color":"#94A3B8","bg":"rgba(71, 85, 105, 0.1)","border":"rgba(148, 163, 184, 0.2)","text":"#CBD5E1","tag_bg":"rgba(71, 85, 105, 0.15)","tag_text":"#94A3B8","shape":"?","severity_dot":"#94A3B8"},
}

SEV_CFG = {
    "none":     {"color":"#4ADE80","bg":"rgba(22, 163, 74, 0.05)","border":"rgba(74, 222, 128, 0.15)","text":"#86EFAC","label":"None"},
    "low":      {"color":"#FBBF24","bg":"rgba(217, 119, 6, 0.05)","border":"rgba(251, 191, 36, 0.15)","text":"#FDE047","label":"Low"},
    "moderate": {"color":"#FB923C","bg":"rgba(234, 88, 12, 0.05)","border":"rgba(251, 146, 60, 0.15)","text":"#FDBA74","label":"Moderate"},
    "high":     {"color":"#F87171","bg":"rgba(220, 38, 38, 0.05)","border":"rgba(248, 113, 113, 0.15)","text":"#FCA5A5","label":"High"},
    "critical": {"color":"#EF4444","bg":"rgba(185, 28, 28, 0.05)","border":"rgba(239, 68, 68, 0.15)","text":"#FECACA","label":"Critical"},
}

PHENO_CFG = {
    "PM":      {"bg":"rgba(220, 38, 38, 0.1)","border":"rgba(248, 113, 113, 0.2)","text":"#FCA5A5","bar":"#EF4444","label":"Poor Metabolizer","pct":5},
    "IM":      {"bg":"rgba(217, 119, 6, 0.1)","border":"rgba(251, 191, 36, 0.2)","text":"#FDE047","bar":"#FBBF24","label":"Intermediate Metabolizer","pct":45},
    "NM":      {"bg":"rgba(22, 163, 74, 0.1)","border":"rgba(74, 222, 128, 0.2)","text":"#86EFAC","bar":"#4ADE80","label":"Normal Metabolizer","pct":100},
    "RM":      {"bg":"rgba(37, 99, 235, 0.1)","border":"rgba(96, 165, 250, 0.2)","text":"#93C5FD","bar":"#60A5FA","label":"Rapid Metabolizer","pct":115},
    "URM":     {"bg":"rgba(234, 88, 12, 0.1)","border":"rgba(251, 146, 60, 0.2)","text":"#FDBA74","bar":"#FB923C","label":"Ultrarapid Metabolizer","pct":130},
    "Unknown": {"bg":"rgba(71, 85, 105, 0.1)","border":"rgba(148, 163, 184, 0.2)","text":"#94A3B8","bar":"#64748B","label":"Unknown","pct":0},
}

POP_FREQ = {
    "CYP2D6":  {"PM":7,"IM":10,"NM":77,"URM":6},
    "CYP2C19": {"PM":3,"IM":26,"NM":52,"RM":13,"URM":6},
    "CYP2C9":  {"PM":1,"IM":10,"NM":89},
    "SLCO1B1": {"PM":1,"IM":15,"NM":84},
    "TPMT":    {"PM":0.3,"IM":10,"NM":90},
    "DPYD":    {"PM":0.2,"IM":3,"NM":97},
}

CHROM_INFO = {
    "CYP2D6":  {"chrom":"22","band":"q13.2","pos_mb":42.5},
    "CYP2C19": {"chrom":"10","band":"q23.33","pos_mb":96.7},
    "CYP2C9":  {"chrom":"10","band":"q23.33","pos_mb":96.4},
    "SLCO1B1": {"chrom":"12","band":"p12.1","pos_mb":21.3},
    "TPMT":    {"chrom":"6","band":"p22.3","pos_mb":18.1},
    "DPYD":    {"chrom":"1","band":"p22.1","pos_mb":97.5},
}
CHROM_LEN = {"1":248.9,"6":170.8,"10":133.8,"12":133.3,"22":50.8}

PLAIN_PHENO = {
    "PM":"Your body barely processes this medicine",
    "IM":"Your body processes this medicine slower than average",
    "NM":"Your body processes this medicine normally",
    "RM":"Your body processes this medicine slightly faster than average",
    "URM":"Your body processes this medicine dangerously fast",
    "Unknown":"Gene function unclear",
}

PLAIN_RISK = {
    ("CODEINE","PM"):      "Your body can't convert codeine into a painkiller — it won't help your pain.",
    ("CODEINE","URM"):     "Your body converts codeine to morphine extremely fast. Even one tablet could be life-threatening.",
    ("CODEINE","IM"):      "Codeine may be less effective. Your doctor may need to try a different painkiller.",
    ("CODEINE","NM"):      "Codeine works normally for you. Standard doses should manage pain safely.",
    ("WARFARIN","PM"):     "Warfarin stays in your body much longer than normal. Standard doses could cause dangerous bleeding.",
    ("WARFARIN","IM"):     "Warfarin clears more slowly. You'll likely need a lower dose.",
    ("WARFARIN","NM"):     "Warfarin works normally for you. Standard INR monitoring applies.",
    ("CLOPIDOGREL","PM"):  "This heart medication won't activate properly, leaving you unprotected against blood clots.",
    ("CLOPIDOGREL","IM"):  "This heart medication activates less than normal. A stronger alternative may be needed.",
    ("CLOPIDOGREL","NM"):  "This heart medication works normally for you.",
    ("SIMVASTATIN","PM"):  "This cholesterol drug can't be cleared properly and may build up in your muscles, causing serious damage.",
    ("SIMVASTATIN","IM"):  "This cholesterol drug clears more slowly. A lower dose will protect your muscles.",
    ("SIMVASTATIN","NM"):  "This cholesterol drug works normally for you.",
    ("AZATHIOPRINE","PM"): "This immune drug builds up to dangerous levels. Standard doses would seriously harm your bone marrow.",
    ("AZATHIOPRINE","IM"): "You need a lower dose of this immune drug to stay safe.",
    ("AZATHIOPRINE","NM"): "This immune drug works normally for you.",
    ("FLUOROURACIL","PM"): "Your body cannot break down this chemotherapy. Standard doses would be life-threatening.",
    ("FLUOROURACIL","IM"): "This chemotherapy breaks down too slowly. You need a significantly reduced dose.",
    ("FLUOROURACIL","NM"): "This chemotherapy works at a normal rate in your body.",
}

PERSONAS = {
    "A":{"label":"Critical Risk","file":"patient_a_critical.vcf","drugs":["CODEINE","FLUOROURACIL","AZATHIOPRINE"],"desc":"CYP2D6 PM · DPYD PM · TPMT PM","sev":"critical"},
    "B":{"label":"Warfarin PM","file":"patient_b_warfarin.vcf","drugs":["WARFARIN"],"desc":"CYP2C9 *2/*3 Poor Metabolizer","sev":"high"},
    "C":{"label":"Drug Interaction","file":"patient_c_interaction.vcf","drugs":["CLOPIDOGREL"],"desc":"CYP2C19 *2/*3 Poor Metabolizer","sev":"high"},
    "D":{"label":"All Safe","file":"patient_d_safe.vcf","drugs":["CODEINE","WARFARIN","SIMVASTATIN"],"desc":"Wildtype *1/*1 all genes","sev":"none"},
}

TEST_SUITE = [
    {"name":"Mixed Variants","file":"sample.vcf","drugs":["CLOPIDOGREL","CODEINE","AZATHIOPRINE"],
     "expected":{"CLOPIDOGREL":"Ineffective","CODEINE":"Ineffective","AZATHIOPRINE":"Toxic"},
     "desc":"CYP2C19 *2/*3 · CYP2D6 *4/*4 · TPMT *3B/*3C"},
    {"name":"UltraRapid Metabolizer","file":"test_ultrarapid_metabolizer.vcf","drugs":["CODEINE","CLOPIDOGREL"],
     "expected":{"CODEINE":"Toxic","CLOPIDOGREL":"Safe"},"desc":"CYP2D6 *1xN/*1xN → URM → Codeine Toxic"},
    {"name":"All Normal Wild-type","file":"test_all_normal_wildtype.vcf","drugs":ALL_DRUGS,
     "expected":{d:"Safe" for d in ALL_DRUGS},"desc":"Wild-type *1/*1 across all 6 genes"},
    {"name":"Worst Case — All PM","file":"test_worst_case_all_pm.vcf","drugs":ALL_DRUGS,
     "expected":{"CODEINE":"Ineffective","CLOPIDOGREL":"Ineffective","WARFARIN":"Adjust Dosage","SIMVASTATIN":"Toxic","AZATHIOPRINE":"Toxic","FLUOROURACIL":"Toxic"},
     "desc":"Loss-of-function alleles across all 6 genes"},
]

# ── Page Config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="SurakshaRx — Pharmacogenomic Risk",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

st.markdown("""<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:ital,opsz,wght@0,9..40,300;0,9..40,400;0,9..40,500;0,9..40,600;0,9..40,700;1,9..40,400&family=JetBrains+Mono:wght@400;500;600&display=swap');

:root {
  --sp-1:4px;--sp-2:8px;--sp-3:12px;--sp-4:16px;--sp-5:20px;--sp-6:24px;--sp-8:32px;--sp-10:40px;--sp-12:48px;--sp-16:64px;
  --bg:#0B0F19;
  --surface:rgba(17, 24, 39, 0.65);
  --surface-sub:rgba(31, 41, 55, 0.6);
  --surface-sub2:rgba(55, 65, 81, 0.6);
  --border-light:rgba(255, 255, 255, 0.08);
  --border:rgba(255, 255, 255, 0.15);
  --border-dark:rgba(255, 255, 255, 0.25);
  --text-primary:#F8FAFC;
  --text-secondary:#CBD5E1;
  --text-muted:#94A3B8;
  --text-xmuted:#64748B;
  --brand:#60A5FA;
  --brand-dark:#93C5FD;
  --brand-hover:#3B82F6;
  --brand-light:rgba(59, 130, 246, 0.15);
  --brand-border:rgba(96, 165, 250, 0.3);
  --safe:#4ADE80;--safe-bg:rgba(74,222,128,0.1);--safe-border:rgba(74,222,128,0.2);
  --warn:#FBBF24;--warn-bg:rgba(251,191,36,0.1);
  --danger:#EF4444;--danger-bg:rgba(239,68,68,0.1);--danger-light:#FCA5A5;
  --shadow-xs:0 2px 4px rgba(0,0,0,0.3);
  --shadow-sm:0 4px 6px rgba(0,0,0,0.4), 0 2px 4px rgba(0,0,0,0.3);
  --shadow-md:0 10px 15px rgba(0,0,0,0.5), 0 4px 6px rgba(0,0,0,0.4);
  --shadow-lg:0 20px 25px rgba(0,0,0,0.6), 0 10px 10px rgba(0,0,0,0.5);
  --r-sm:6px;--r-md:8px;--r-lg:12px;--r-xl:16px;--r-2xl:20px;--r-full:9999px;
  --font-body:'DM Sans',-apple-system,BlinkMacSystemFont,sans-serif;
  --font-mono:'JetBrains Mono','Fira Code',monospace;
}

*,*::before,*::after{box-sizing:border-box;}
html,body,[class*="css"]{font-family:var(--font-body)!important;font-size:16px!important;background:var(--bg)!important;color:var(--text-primary)!important;-webkit-font-smoothing:antialiased!important;}
.stApp{background:var(--bg)!important; background-image: radial-gradient(circle at 15% 50%, rgba(59, 130, 246, 0.08), transparent 25%), radial-gradient(circle at 85% 30%, rgba(167, 139, 250, 0.08), transparent 25%)!important;}
.main .block-container{padding:0 var(--sp-10) var(--sp-16)!important;max-width:1280px!important;}
#MainMenu,footer,header{visibility:hidden;}

.stMarkdown,.stMarkdown p,.stMarkdown li,.stMarkdown span{color:var(--text-secondary)!important;}
[data-testid="stMarkdownContainer"] *{color:inherit;}

@keyframes fade-up{from{opacity:0;transform:translateY(12px)}to{opacity:1;transform:translateY(0)}}
@keyframes fade-in{from{opacity:0}to{opacity:1}}
@keyframes pulse-once{0%{transform:scale(1);box-shadow:0 0 0 0 rgba(239,68,68,.4)}40%{transform:scale(1.02);box-shadow:0 0 0 10px rgba(239,68,68,0)}100%{transform:scale(1);box-shadow:0 0 0 0 rgba(239,68,68,0)}}
@keyframes bar-fill{from{width:0!important}}
@keyframes score-count{from{opacity:0;transform:scale(.85)}to{opacity:1;transform:scale(1)}}

.reveal-card{animation:fade-up .32s cubic-bezier(.4,0,.2,1) both;}
.reveal-card:nth-child(1){animation-delay:.04s}.reveal-card:nth-child(2){animation-delay:.10s}
.reveal-card:nth-child(3){animation-delay:.16s}.reveal-card:nth-child(4){animation-delay:.22s}
.reveal-card:nth-child(5){animation-delay:.28s}.reveal-card:nth-child(6){animation-delay:.34s}

.pg-nav{display:flex;align-items:center;justify-content:space-between;padding:var(--sp-6) 0;border-bottom:1px solid var(--border-light);margin-bottom:var(--sp-8);animation:fade-in .4s ease;}
.pg-brand-name{font-size:1.5rem;font-weight:700;color:var(--text-primary)!important;letter-spacing:-.03em;line-height:1;text-shadow: 0 0 15px rgba(96,165,250,0.3);}
.pg-brand-name span{color:var(--brand)!important;}
.pg-brand-sub{font-family:var(--font-mono);font-size:.8rem;color:var(--text-xmuted)!important;letter-spacing:.1em;text-transform:uppercase;}
.pg-nav-badges{display:flex;align-items:center;gap:var(--sp-2);}
.pg-badge{font-family:var(--font-mono);font-size:.8rem;font-weight:500;letter-spacing:.06em;text-transform:uppercase;padding:4px 10px;border-radius:var(--r-full);border:1px solid;white-space:nowrap;backdrop-filter:blur(8px);}
.pg-badge-default{color:var(--text-secondary)!important;border-color:var(--border);background:var(--surface);}
.pg-badge-brand{color:var(--brand)!important;border-color:var(--brand-border);background:var(--brand-light);font-weight:600;}

.trust-strip{display:flex;align-items:center;gap:var(--sp-6);padding:var(--sp-3) var(--sp-4);background:var(--brand-light);border:1px solid var(--brand-border);border-radius:var(--r-md);margin-bottom:var(--sp-8);backdrop-filter:blur(10px);}
.trust-item{display:flex;align-items:center;gap:var(--sp-2);font-size:.85rem;color:var(--brand-dark)!important;font-weight:500;white-space:nowrap;}
.trust-sep{width:1px;height:16px;background:var(--brand-border);flex-shrink:0;}

.stTabs [data-baseweb="tab-list"]{background:transparent!important;border-bottom:1px solid var(--border-light)!important;gap:0!important;padding:0!important;margin-bottom:var(--sp-8)!important;box-shadow:none!important;}
.stTabs [data-baseweb="tab"]{font-family:var(--font-body)!important;font-size:1rem!important;font-weight:500!important;color:var(--text-muted)!important;padding:var(--sp-3) var(--sp-5)!important;background:transparent!important;border:none!important;border-bottom:2.5px solid transparent!important;border-radius:0!important;transition:color .15s, border-color .15s!important;}
.stTabs [aria-selected="true"]{color:var(--brand)!important;border-bottom-color:var(--brand)!important;font-weight:600!important;text-shadow: 0 0 10px rgba(96,165,250,0.4)!important;}
.stTabs [data-baseweb="tab"] span,.stTabs [data-baseweb="tab"] div,.stTabs [data-baseweb="tab"] p{color:inherit!important;}
.stTabs [data-baseweb="tab-panel"]{padding-top:0!important;background:transparent!important;}
.stTabs [data-baseweb="tab-panel"] p,.stTabs [data-baseweb="tab-panel"] div,.stTabs [data-baseweb="tab-panel"] span,.stTabs [data-baseweb="tab-panel"] label{color:var(--text-secondary)!important;}
.stTabs [data-baseweb="tab-panel"] h1,.stTabs [data-baseweb="tab-panel"] h2,.stTabs [data-baseweb="tab-panel"] h3{color:var(--text-primary)!important;}

[data-testid="stSidebar"]{background:rgba(15,23,42,0.8)!important;border-right:1px solid var(--border-light)!important;backdrop-filter:blur(12px)!important;}
[data-testid="stSidebar"] *{color:var(--text-primary)!important;}
[data-testid="stSidebar"] label,[data-testid="stSidebar"] .stMarkdown p{color:var(--text-secondary)!important;}
[data-testid="stSidebar"] h3,[data-testid="stSidebar"] h4{color:var(--text-primary)!important;font-weight:700!important;}
[data-testid="stSidebar"] code{color:var(--brand-dark)!important;background:var(--brand-light)!important;padding:1px 5px;border-radius:3px;font-family:var(--font-mono)!important;}

.sec-label{display:flex;align-items:center;gap:var(--sp-3);font-size:.8rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:var(--text-muted)!important;margin-bottom:var(--sp-4);}
.sec-label::after{content:'';flex:1;height:1px;background:var(--border-light);}

.steps{display:flex;background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-xl);overflow:hidden;margin-bottom:var(--sp-8);box-shadow:var(--shadow-md);backdrop-filter:blur(12px);}
.step{flex:1;padding:var(--sp-4) var(--sp-5);border-right:1px solid var(--border-light);transition:background 0.3s;}
.step:last-child{border-right:none;}
.step-num{font-family:var(--font-mono);font-size:.7rem;font-weight:600;letter-spacing:.12em;text-transform:uppercase;color:var(--text-xmuted)!important;margin-bottom:3px;}
.step-lbl{font-size:.875rem;font-weight:500;color:var(--text-muted)!important;}
.step.done .step-num{color:var(--brand)!important;text-shadow:0 0 8px rgba(96,165,250,0.5);}
.step.done .step-lbl{color:var(--text-primary)!important;font-weight:600;}
.step.done{background:rgba(59,130,246,0.05);}

.persona-card{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-lg);padding:var(--sp-4);transition:all .2s cubic-bezier(.4,0,.2,1);box-shadow:var(--shadow-sm);cursor:pointer;backdrop-filter:blur(12px);}
.persona-card:hover{transform:translateY(-2px);box-shadow:0 0 15px rgba(255,255,255,0.05), var(--shadow-md);border-color:var(--border-dark);}
.pc-sev{display:inline-flex;align-items:center;gap:5px;font-size:.75rem;font-weight:600;padding:3px 10px;border-radius:var(--r-full);border:1px solid;margin-bottom:var(--sp-2);}
.pc-name{font-size:.875rem;font-weight:700;color:var(--text-primary)!important;margin-bottom:3px;}
.pc-desc{font-family:var(--font-mono);font-size:.7rem;color:var(--text-muted)!important;line-height:1.7;}

.risk-center{border-radius:var(--r-2xl);padding:var(--sp-8);margin-bottom:var(--sp-6);border:1px solid;position:relative;overflow:hidden;box-shadow:var(--shadow-lg);backdrop-filter:blur(16px);}
.rc-eyebrow{font-family:var(--font-mono);font-size:.7rem;font-weight:600;letter-spacing:.14em;text-transform:uppercase;opacity:.8;margin-bottom:var(--sp-1);}
.rc-headline{font-size:2.5rem;font-weight:700;letter-spacing:-.03em;line-height:1.1;margin-bottom:var(--sp-1);text-shadow:0 2px 10px rgba(0,0,0,0.5);}
.rc-sub{font-size:.9rem;opacity:.9;margin-bottom:var(--sp-5);}
.rc-stats{display:grid;grid-template-columns:repeat(4,1fr);gap:var(--sp-5);padding-top:var(--sp-5);border-top:1px solid;border-color:inherit;opacity:.9;}
.rc-stat-num{font-size:2rem;font-weight:700;letter-spacing:-.03em;line-height:1;margin-bottom:3px;text-shadow:0 2px 4px rgba(0,0,0,0.4);}
.rc-stat-lbl{font-family:var(--font-mono);font-size:.65rem;font-weight:500;letter-spacing:.1em;text-transform:uppercase;opacity:.8;}

.crit-alert{display:flex;gap:var(--sp-4);background:rgba(239,68,68,0.1);border:1px solid rgba(248,113,113,0.3);border-left:4px solid var(--danger);border-radius:var(--r-lg);padding:var(--sp-4) var(--sp-5);margin-bottom:var(--sp-4);animation:pulse-once .8s ease .3s both;box-shadow:var(--shadow-md);backdrop-filter:blur(8px);}
.crit-title{font-size:.95rem;font-weight:700;color:var(--danger-light)!important;margin-bottom:3px;text-shadow:0 0 10px rgba(239,68,68,0.4);}
.crit-note{font-size:.875rem;color:#FECACA!important;line-height:1.65;margin-bottom:var(--sp-2);}
.crit-action{font-family:var(--font-mono);font-size:.7rem;font-weight:600;color:var(--danger)!important;letter-spacing:.08em;text-transform:uppercase;}

.gene-row{display:grid;grid-template-columns:repeat(6,1fr);gap:var(--sp-3);margin-bottom:var(--sp-6);}
.gene-box{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-lg);padding:var(--sp-4) var(--sp-3);text-align:center;box-shadow:var(--shadow-sm);transition:box-shadow .2s,transform .2s,border-color .2s;backdrop-filter:blur(10px);}
.gene-box:hover{box-shadow:0 0 15px rgba(255,255,255,0.05), var(--shadow-md);transform:translateY(-2px);border-color:var(--border-dark);}
.gene-nm{font-family:var(--font-mono);font-size:.75rem;font-weight:600;margin-bottom:var(--sp-2);color:var(--text-secondary)!important;}
.gene-track{height:3px;border-radius:2px;background:var(--surface-sub2);margin:var(--sp-2) 0;overflow:hidden;}
.gene-fill{height:100%;border-radius:2px;box-shadow:0 0 8px currentColor;}
.gene-ph{font-family:var(--font-mono);font-size:.8rem;font-weight:600;letter-spacing:.03em;}

.dtab{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-xl);overflow:hidden;margin-bottom:var(--sp-6);box-shadow:var(--shadow-md);backdrop-filter:blur(12px);}
.dtab-head{display:grid;grid-template-columns:1.4fr 1.2fr .9fr 1fr .9fr 1.1fr;background:var(--surface-sub);border-bottom:1px solid var(--border-light);}
.dtab-hcell{font-family:var(--font-mono);font-size:.65rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:var(--text-muted)!important;padding:var(--sp-3) var(--sp-4);}
.dtab-row{display:grid;grid-template-columns:1.4fr 1.2fr .9fr 1fr .9fr 1.1fr;border-bottom:1px solid var(--border-light);transition:background .15s;}
.dtab-row:last-child{border-bottom:none;}.dtab-row:hover{background:rgba(255,255,255,0.05);}
.dtab-cell{font-size:.9rem;color:var(--text-primary)!important;padding:var(--sp-3) var(--sp-4);display:flex;align-items:center;}

.risk-badge{display:inline-flex;align-items:center;gap:6px;font-size:.8rem;font-weight:600;padding:4px 12px;border-radius:var(--r-full);border:1px solid;letter-spacing:-.01em;}

.pgx-card{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-2xl);padding:var(--sp-8);margin-bottom:var(--sp-6);box-shadow:var(--shadow-lg);position:relative;overflow:hidden;backdrop-filter:blur(12px);}
.pgx-card::before{content:'';position:absolute;top:0;right:0;width:280px;height:280px;background:radial-gradient(circle at top right,var(--brand-light) 0%,transparent 65%);pointer-events:none;}
.pgx-eyebrow{font-family:var(--font-mono);font-size:.7rem;font-weight:600;letter-spacing:.14em;text-transform:uppercase;color:var(--brand)!important;margin-bottom:var(--sp-2);text-shadow:0 0 10px rgba(96,165,250,0.3);}
.pgx-score{font-size:4.5rem;font-weight:700;letter-spacing:-.04em;line-height:1;margin-bottom:4px;animation:score-count .5s cubic-bezier(.4,0,.2,1) .1s both;text-shadow:0 0 20px currentColor;}
.pgx-label{font-size:.9rem;color:var(--text-secondary)!important;margin-bottom:var(--sp-5);}
.pgx-marker{position:relative;height:6px;background:var(--surface-sub2);border-radius:3px;overflow:visible;margin-bottom:var(--sp-5);box-shadow:inset 0 1px 2px rgba(0,0,0,0.5);}
.pgx-fill{position:absolute;top:0;left:0;height:100%;border-radius:3px;transition:width .9s cubic-bezier(.4,0,.2,1);box-shadow:0 0 10px currentColor;}
.pgx-indicator{position:absolute;top:-4px;width:14px;height:14px;border-radius:50%;background:var(--bg);border:3px solid;transform:translateX(-50%);box-shadow:0 0 10px currentColor;transition:left .9s cubic-bezier(.4,0,.2,1);}
.pgx-thresh-labels{display:flex;justify-content:space-between;font-family:var(--font-mono);font-size:.65rem;color:var(--text-muted)!important;margin-bottom:var(--sp-3);}
.pgx-pills{display:flex;flex-wrap:wrap;gap:var(--sp-2);}
.pgx-pill{font-family:var(--font-mono);font-size:.7rem;font-weight:600;padding:3px 10px;border-radius:var(--r-full);border:1px solid;letter-spacing:.03em;}

.hm-wrap{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-xl);padding:var(--sp-6);margin-bottom:var(--sp-6);box-shadow:var(--shadow-md);overflow-x:auto;backdrop-filter:blur(12px);}
.hm-eyebrow{font-family:var(--font-mono);font-size:.7rem;font-weight:600;letter-spacing:.12em;text-transform:uppercase;color:var(--text-muted)!important;margin-bottom:var(--sp-5);}
.hm-grid{display:grid;gap:3px;}
.hm-cell{border-radius:var(--r-sm);display:flex;flex-direction:column;align-items:center;justify-content:center;padding:var(--sp-3) var(--sp-2);min-height:56px;border:1px solid;transition:transform .2s,box-shadow .2s;cursor:default;}
.hm-cell:hover{transform:scale(1.06);box-shadow:0 0 15px currentColor;z-index:5;position:relative;}
.hm-cell-name{font-family:var(--font-mono);font-size:.7rem;font-weight:600;margin-bottom:2px;}
.hm-cell-risk{font-family:var(--font-mono);font-size:.65rem;opacity:.9;}
.hm-header{font-family:var(--font-mono);font-size:.65rem;letter-spacing:.05em;color:var(--text-secondary)!important;display:flex;align-items:center;justify-content:center;min-height:56px;}
.hm-legend{display:flex;gap:var(--sp-5);margin-top:var(--sp-4);flex-wrap:wrap;}
.hm-legend-item{font-family:var(--font-mono);font-size:.7rem;display:flex;align-items:center;gap:5px;color:var(--text-secondary)!important;}
.hm-dot{width:10px;height:10px;border-radius:3px;display:inline-block;border:1px solid;box-shadow:0 0 5px currentColor;}

.chrom-wrap{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-xl);padding:var(--sp-5) var(--sp-6);box-shadow:var(--shadow-md);backdrop-filter:blur(12px);}
.chrom-eyebrow{font-family:var(--font-mono);font-size:.7rem;font-weight:600;letter-spacing:.12em;text-transform:uppercase;color:var(--text-muted)!important;margin-bottom:var(--sp-4);}
.chrom-row{display:flex;align-items:center;gap:var(--sp-3);margin-bottom:var(--sp-2);}
.chrom-chr{font-family:var(--font-mono);font-size:.75rem;color:var(--text-muted)!important;width:18px;text-align:right;flex-shrink:0;}
.chrom-bar{flex:1;height:11px;background:var(--surface-sub2);border-radius:6px;position:relative;overflow:visible;border:1px solid var(--border-light);box-shadow:inset 0 1px 3px rgba(0,0,0,0.5);}
.chrom-body{position:absolute;inset:0;background:linear-gradient(90deg,rgba(255,255,255,0.05),rgba(255,255,255,0.1),rgba(255,255,255,0.05));border-radius:6px;}
.chrom-marker{position:absolute;top:-5px;width:3px;height:21px;border-radius:2px;transform:translateX(-50%);}
.chrom-gene{font-family:var(--font-mono);font-size:.75rem;color:var(--text-secondary)!important;width:56px;flex-shrink:0;font-weight:500;}
.chrom-band{font-family:var(--font-mono);font-size:.65rem;color:var(--text-xmuted)!important;}

.pop-wrap{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-lg);padding:var(--sp-4) var(--sp-5);margin-bottom:var(--sp-4);box-shadow:var(--shadow-sm);backdrop-filter:blur(10px);}
.pop-eyebrow{font-family:var(--font-mono);font-size:.65rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:var(--text-muted)!important;margin-bottom:var(--sp-3);}
.pop-row{display:flex;align-items:center;gap:var(--sp-3);margin-bottom:var(--sp-2);}
.pop-ph{font-family:var(--font-mono);font-size:.75rem;color:var(--text-secondary)!important;width:96px;flex-shrink:0;font-weight:500;}
.pop-track{flex:1;height:4px;background:var(--surface-sub2);border-radius:2px;overflow:hidden;box-shadow:inset 0 1px 2px rgba(0,0,0,0.5);}
.pop-fill{height:100%;border-radius:2px;animation:bar-fill .7s cubic-bezier(.4,0,.2,1) both;box-shadow:0 0 5px currentColor;}
.pop-pct{font-family:var(--font-mono);font-size:.7rem;width:32px;text-align:right;color:var(--text-muted)!important;}
.pop-you{font-family:var(--font-mono);font-size:.65rem;color:var(--brand)!important;font-weight:700;margin-left:3px;text-shadow:0 0 5px rgba(96,165,250,0.5);}

.ix-grid{display:grid;gap:3px;}
.ix-cell{border-radius:var(--r-sm);display:flex;align-items:center;justify-content:center;min-height:44px;font-family:var(--font-mono);font-size:.65rem;text-align:center;padding:var(--sp-1);font-weight:700;border:1px solid;transition:transform .2s,box-shadow .2s;}
.ix-cell:hover{transform:scale(1.06);z-index:5;position:relative;box-shadow:0 0 15px currentColor;}
.ix-head{font-family:var(--font-mono);font-size:.65rem;letter-spacing:.05em;color:var(--text-secondary)!important;display:flex;align-items:center;justify-content:center;min-height:44px;}

.dcard{background:var(--surface);border:1px solid var(--border-light);border-radius:var(--r-2xl);margin-bottom:var(--sp-6);overflow:hidden;box-shadow:var(--shadow-md);transition:box-shadow .3s, transform .3s;backdrop-filter:blur(12px);}
.dcard:hover{box-shadow:0 8px 30px rgba(0,0,0,0.6), 0 0 20px rgba(255,255,255,0.03); transform:translateY(-2px);}
.dcard-header{display:flex;align-items:center;justify-content:space-between;padding:var(--sp-5) var(--sp-6);border-bottom:1px solid var(--border-light);background:rgba(255,255,255,0.02);}
.dcard-left{display:flex;align-items:center;gap:var(--sp-4);}
.dcard-indicator{width:10px;height:10px;border-radius:50%;flex-shrink:0;box-shadow:0 0 10px currentColor;}
.dcard-drug{font-size:1.125rem;font-weight:700;letter-spacing:-.02em;color:var(--text-primary)!important;}
.dcard-meta{font-family:var(--font-mono);font-size:.75rem;color:var(--text-secondary)!important;margin-top:3px;}
.dcard-body{padding:var(--sp-6);}

.metrics-row{display:grid;grid-template-columns:repeat(4,1fr);gap:1px;background:var(--border-light);border-radius:var(--r-lg);overflow:hidden;border:1px solid var(--border-light);margin-bottom:var(--sp-5);}
.metric-cell{background:var(--surface-sub);padding:var(--sp-4);}
.metric-key{font-family:var(--font-mono);font-size:.65rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:var(--text-muted)!important;margin-bottom:4px;}
.metric-val{font-size:1.125rem;font-weight:700;color:var(--text-primary)!important;letter-spacing:-.02em;}

.conf-grid{display:grid;grid-template-columns:1fr 1fr;gap:var(--sp-5);margin-bottom:var(--sp-5);}
.conf-label{font-family:var(--font-mono);font-size:.65rem;font-weight:600;letter-spacing:.08em;text-transform:uppercase;color:var(--text-muted)!important;display:flex;justify-content:space-between;margin-bottom:5px;}
.conf-track{height:4px;background:var(--surface-sub2);border-radius:2px;overflow:hidden;box-shadow:inset 0 1px 2px rgba(0,0,0,0.5);}
.conf-fill{height:100%;border-radius:2px;animation:bar-fill .7s cubic-bezier(.4,0,.2,1) both;box-shadow:0 0 5px currentColor;}

.vtable{width:100%;border-collapse:collapse;}
.vtable th{font-family:var(--font-mono);font-size:.65rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:var(--text-muted)!important;padding:0 var(--sp-3) var(--sp-3);text-align:left;border-bottom:1px solid var(--border-light);}
.vtable td{font-family:var(--font-mono);font-size:.85rem;color:var(--text-secondary)!important;padding:var(--sp-2) var(--sp-3);border-bottom:1px solid var(--border-light);}
.vtable tbody tr:last-child td{border-bottom:none;}.vtable tbody tr:hover td{background:rgba(255,255,255,0.05);}
.v-rsid{color:#60A5FA!important;font-weight:500!important;}.v-star{color:#A78BFA!important;font-weight:500!important;}
.v-nofunc{color:var(--danger)!important;font-weight:500!important;}.v-dec{color:var(--warn)!important;font-weight:500!important;}.v-norm{color:var(--safe)!important;font-weight:500!important;}

.rec-box{border-radius:var(--r-lg);border:1px solid;padding:var(--sp-4) var(--sp-5);margin-bottom:var(--sp-4);box-shadow:inset 0 0 20px rgba(0,0,0,0.2);}
.rec-label{font-family:var(--font-mono);font-size:.65rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;margin-bottom:var(--sp-2);}
.rec-text{font-size:.95rem;line-height:1.75;color:var(--text-primary)!important;}
.alt-chips{display:flex;flex-wrap:wrap;gap:var(--sp-2);}
.alt-chip{font-family:var(--font-mono);font-size:.75rem;font-weight:500;color:var(--brand-dark)!important;background:var(--brand-light);border:1px solid var(--brand-border);border-radius:var(--r-full);padding:4px 12px;text-shadow:0 0 5px rgba(59,130,246,0.3);}
.cpic-badge{font-family:var(--font-mono);font-size:.65rem;font-weight:700;background:rgba(251,191,36,0.1);border:1px solid rgba(251,191,36,0.3);color:#FDE047!important;padding:2px 8px;border-radius:4px;display:inline-block;margin-left:var(--sp-2);vertical-align:middle;letter-spacing:.05em;}

.ai-block{background:linear-gradient(135deg,rgba(59,130,246,0.05) 0%,rgba(167,139,250,0.05) 100%);border:1px solid var(--brand-border);border-radius:var(--r-xl);overflow:hidden;margin-bottom:var(--sp-5);box-shadow:0 4px 15px rgba(59,130,246,0.1);}
.ai-header{display:flex;align-items:center;gap:var(--sp-3);padding:var(--sp-3) var(--sp-5);background:rgba(255,255,255,0.03);border-bottom:1px solid var(--border-light);}
.ai-badge-pill{font-family:var(--font-mono);font-size:.7rem;font-weight:600;letter-spacing:.08em;text-transform:uppercase;background:var(--brand-light);border:1px solid var(--brand-border);color:var(--brand-dark)!important;padding:3px 9px;border-radius:var(--r-sm);text-shadow:0 0 5px rgba(59,130,246,0.3);}
.ai-title{font-size:.9rem;font-weight:600;color:var(--text-primary)!important;}
.ai-section{padding:var(--sp-4) var(--sp-5);border-bottom:1px solid var(--border-light);}
.ai-section:last-child{border-bottom:none;}.ai-section:hover{background:rgba(255,255,255,0.02);}
.ai-sec-label{font-family:var(--font-mono);font-size:.65rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:var(--brand)!important;margin-bottom:var(--sp-2);}
.ai-sec-text{font-size:.9rem;line-height:1.8;color:var(--text-secondary)!important;}

.narrative-box{background:var(--brand-light);border:1px solid var(--brand-border);border-radius:var(--r-xl);padding:var(--sp-6);margin-bottom:var(--sp-6);box-shadow:inset 0 0 30px rgba(59,130,246,0.05), var(--shadow-md);backdrop-filter:blur(10px);}
.narrative-header{display:flex;align-items:center;gap:var(--sp-3);margin-bottom:var(--sp-4);}
.narrative-text{font-size:.95rem;line-height:1.85;color:var(--text-primary)!important;}

.ba-grid{display:grid;grid-template-columns:1fr 1fr;border:1px solid var(--border-light);border-radius:var(--r-xl);overflow:hidden;margin-bottom:var(--sp-6);box-shadow:var(--shadow-md);backdrop-filter:blur(12px);}
.ba-side{padding:var(--sp-6);}
.ba-side-lbl{font-family:var(--font-mono);font-size:.7rem;font-weight:700;letter-spacing:.1em;text-transform:uppercase;margin-bottom:var(--sp-3);}
.ba-drug{font-size:.9rem;font-weight:700;margin-bottom:3px;}
.ba-text{font-size:.875rem;line-height:1.65;}
.ba-gene{font-family:var(--font-mono);font-size:.7rem;margin-top:var(--sp-3);opacity:.8;}

.rx-result{border-radius:var(--r-lg);padding:var(--sp-5) var(--sp-6);margin-top:var(--sp-4);border:1px solid;animation:fade-up .25s ease;box-shadow:var(--shadow-md);backdrop-filter:blur(10px);}
.rx-verdict{font-size:.95rem;font-weight:700;margin-bottom:var(--sp-2);}
.rx-detail{font-size:.875rem;line-height:1.7;color:var(--text-primary)!important;margin-bottom:var(--sp-2);}
.rx-meta{font-family:var(--font-mono);font-size:.7rem;letter-spacing:.06em;text-transform:uppercase;}

.note-box{background:var(--surface-sub);border:1px solid var(--border-light);border-radius:var(--r-xl);padding:var(--sp-6);box-shadow:inset 0 2px 10px rgba(0,0,0,0.5);}
.note-box pre{
  font-family:'JetBrains Mono','Fira Code',monospace!important;
  font-size:.85rem!important;
  color:var(--text-secondary)!important;
  background:transparent!important;
  line-height:1.85!important;
  white-space:pre-wrap!important;
  word-break:break-word!important;
  font-weight:400!important;
  margin:0!important;
  padding:0!important;
}
.note-box pre *,.note-box code{color:var(--text-secondary)!important;background:transparent!important;}

.patient-banner{border-radius:var(--r-xl);padding:var(--sp-6);margin-bottom:var(--sp-6);border:1px solid;box-shadow:var(--shadow-md);backdrop-filter:blur(12px);}
.patient-banner-title{font-size:1.125rem;font-weight:700;margin-bottom:var(--sp-2);}
.patient-banner-sub{font-size:.9rem;line-height:1.7;opacity:.9;}
.pcard{background:var(--surface);border:1px solid;border-radius:var(--r-xl);padding:var(--sp-6);margin-bottom:var(--sp-4);box-shadow:var(--shadow-sm);animation:fade-up .3s ease both;transition:box-shadow .3s, transform .3s;backdrop-filter:blur(12px);}
.pcard:hover{box-shadow:0 8px 25px rgba(0,0,0,0.5), 0 0 15px rgba(255,255,255,0.03); transform:translateY(-2px);}
.pcard-drug{font-size:1.1rem;font-weight:700;letter-spacing:-.02em;margin-bottom:3px;color:var(--text-primary)!important;}
.pcard-verdict{font-size:.9rem;font-weight:600;line-height:1.5;margin-bottom:var(--sp-2);}
.pcard-gene{font-family:var(--font-mono);font-size:.7rem;letter-spacing:.06em;color:var(--text-muted)!important;margin-bottom:var(--sp-3);font-weight:600;text-transform:uppercase;}
.pcard-plain{font-size:.9rem;line-height:1.8;color:var(--text-secondary)!important;}
.pcard-action{display:flex;align-items:flex-start;gap:var(--sp-3);background:var(--surface-sub);border:1px solid var(--border-light);border-radius:var(--r-lg);padding:var(--sp-4);margin-top:var(--sp-4);}
.pcard-action-text{font-size:.875rem;color:var(--text-primary)!important;line-height:1.65;}

.disclaimer-box{display:flex;gap:var(--sp-4);background:rgba(251,191,36,0.1);border:1px solid rgba(251,191,36,0.3);border-left:3px solid var(--warn);border-radius:var(--r-lg);padding:var(--sp-4) var(--sp-5);margin-bottom:var(--sp-6);backdrop-filter:blur(8px);}
.disclaimer-text{font-size:.85rem;color:#FDE047!important;line-height:1.7;}

.tc-card{background:var(--surface)!important;border:1px solid var(--border-light)!important;border-radius:16px!important;padding:20px!important;box-shadow:var(--shadow-sm)!important;margin-bottom:12px!important;backdrop-filter:blur(10px)!important;}
.tc-name{font-size:.95rem!important;font-weight:700!important;color:var(--text-primary)!important;margin-bottom:4px!important;display:block!important;}
.tc-desc{font-family:'JetBrains Mono',monospace!important;font-size:.7rem!important;color:var(--text-muted)!important;margin-bottom:16px!important;line-height:1.7!important;display:block!important;}
.tc-status-pass{background:rgba(74,222,128,0.1);border:1px solid rgba(74,222,128,0.3);border-radius:8px;padding:10px 14px;margin-top:8px;font-family:'JetBrains Mono',monospace;font-size:.8rem;color:#86EFAC;}
.tc-status-fail{background:rgba(239,68,68,0.1);border:1px solid rgba(239,68,68,0.3);border-radius:8px;padding:10px 14px;margin-top:8px;font-family:'JetBrains Mono',monospace;font-size:.8rem;color:#FCA5A5;}

.empty-state{text-align:center;padding:5rem 2rem;border:1px dashed var(--border-dark);border-radius:var(--r-2xl);background:var(--surface);box-shadow:inset 0 0 50px rgba(0,0,0,0.3);backdrop-filter:blur(12px);}
.empty-icon{font-size:2.5rem;display:block;margin-bottom:var(--sp-4);opacity:.6;animation:pulse-once 2s infinite ease-in-out alternate;}
.empty-title{font-size:1.125rem;font-weight:600;color:var(--text-primary)!important;margin-bottom:var(--sp-2);}
.empty-hint{font-family:var(--font-mono);font-size:.7rem;color:var(--text-muted)!important;letter-spacing:.04em;line-height:2.4;}

.info-strip{display:flex;align-items:flex-start;gap:var(--sp-3);background:var(--brand-light);border:1px solid var(--brand-border);border-radius:var(--r-md);padding:var(--sp-4);margin-bottom:var(--sp-4);backdrop-filter:blur(8px);}
.info-strip-text{font-size:.875rem;color:var(--text-primary)!important;line-height:1.65;}

.input-section-gap{margin-top:16px;margin-bottom:8px;}
.input-field-label{font-size:.875rem;font-weight:600;color:var(--text-secondary);margin-bottom:6px;display:block;}

.stButton>button{background:linear-gradient(135deg,var(--brand) 0%,var(--brand-dark) 100%)!important;color:#FFFFFF!important;border:none!important;border-radius:var(--r-md)!important;font-family:var(--font-body)!important;font-weight:600!important;font-size:.95rem!important;padding:.6875rem 1.75rem!important;height:48px!important;transition:all .2s!important;box-shadow:0 4px 15px rgba(59,130,246,0.3)!important;min-height:44px!important;}
.stButton>button:hover{background:linear-gradient(135deg,var(--brand-hover) 0%,var(--brand) 100%)!important;box-shadow:0 6px 20px rgba(59,130,246,0.5)!important;transform:translateY(-1px)!important;}
.stButton>button *,.stButton>button span,.stButton>button p,.stButton>button div{color:#FFFFFF!important;text-shadow:0 1px 2px rgba(0,0,0,0.2)!important;}
[data-testid="stFileUploader"] .stButton>button,
[data-testid="stFileUploader"] button{background:rgba(255,255,255,0.05)!important;color:var(--brand)!important;border:1px solid var(--brand-border)!important;box-shadow:none!important;transform:none!important;height:auto!important;min-height:36px!important;padding:6px 18px!important;backdrop-filter:blur(8px)!important;}
[data-testid="stFileUploader"] .stButton>button:hover,
[data-testid="stFileUploader"] button:hover{background:var(--brand-light)!important;box-shadow:0 0 10px rgba(59,130,246,0.2)!important;}
[data-testid="stFileUploader"] .stButton>button *,
[data-testid="stFileUploader"] button *,
[data-testid="stFileUploader"] button span{color:var(--brand)!important;text-shadow:none!important;}

.stDownloadButton>button{background:var(--surface-sub)!important;color:var(--text-primary)!important;border:1px solid var(--border-dark)!important;border-radius:var(--r-md)!important;font-family:var(--font-mono)!important;font-size:.8rem!important;padding:.5rem 1rem!important;transition:all .2s!important;box-shadow:var(--shadow-sm)!important;min-height:44px!important;backdrop-filter:blur(8px)!important;}
.stDownloadButton>button:hover{background:var(--brand-light)!important;color:var(--brand-dark)!important;border-color:var(--brand-border)!important;box-shadow:0 0 15px rgba(59,130,246,0.2)!important;}
.stDownloadButton>button *,.stDownloadButton>button span,.stDownloadButton>button p{color:inherit!important;text-shadow:none!important;}

[data-testid="stFileUploader"]{color:var(--text-primary)!important;}
[data-testid="stFileUploader"] label{color:var(--text-secondary)!important;font-size:.875rem!important;font-weight:600!important;margin-bottom:6px!important;}
[data-testid="stFileUploaderDropZone"]{border:none!important;background:transparent!important;padding:0!important;margin:0!important;}
[data-testid="stFileUploaderDropZone"] > div > div:first-child,
[data-testid="stFileUploaderDropZone"] span:not(:has(button)),
[data-testid="stFileUploaderDropZone"] p,
[data-testid="stFileUploaderDropZone"] small,
[data-testid="stFileUploaderDropZone"] svg{display:none!important;}
[data-testid="stFileUploader"] section > div,
[data-testid="stFileUploader"] > label + div > div{border:none!important;background:transparent!important;padding:0!important;min-height:0!important;}
[data-testid="stFileUploader"] button,
[data-testid="stFileUploader"] .stButton>button,
[data-testid="stFileUploaderDropZone"] button,
[data-testid="stFileUploader"] section button{background:linear-gradient(135deg,var(--brand) 0%,var(--brand-dark) 100%)!important;color:#FFFFFF!important;border:none!important;border-radius:var(--r-md)!important;font-weight:600!important;font-size:.95rem!important;padding:10px 20px!important;width:100%!important;box-shadow:0 4px 15px rgba(59,130,246,0.3)!important;height:44px!important;cursor:pointer!important;transition:all .2s!important;}
[data-testid="stFileUploader"] button:hover,
[data-testid="stFileUploaderDropZone"] button:hover{box-shadow:0 6px 20px rgba(59,130,246,0.5)!important;}
[data-testid="stFileUploader"] button span,
[data-testid="stFileUploader"] button p,
[data-testid="stFileUploader"] button div,
[data-testid="stFileUploaderDropZone"] button span,
[data-testid="stFileUploaderDropZone"] button *{color:#FFFFFF!important;text-shadow:0 1px 2px rgba(0,0,0,0.2)!important;}

[data-testid="stFileUploaderFile"]{background:var(--surface-sub)!important;border:1px solid var(--border-light)!important;border-radius:var(--r-md)!important;padding:8px 12px!important;margin-top:8px!important;backdrop-filter:blur(8px)!important;}
[data-testid="stFileUploaderFile"] span,
[data-testid="stFileUploaderFile"] p,
[data-testid="stFileUploaderFile"] div{color:var(--text-primary)!important;font-weight:500!important;font-family:var(--font-body)!important;}
[data-testid="stFileUploaderFile"] small{color:var(--text-muted)!important;font-weight:400!important;font-size:.8rem!important;}

.stTextInput>div>div>input{border-radius:var(--r-md)!important;border:1px solid var(--border-dark)!important;background:var(--surface)!important;color:var(--text-primary)!important;font-family:var(--font-body)!important;font-size:.95rem!important;padding:.6875rem .875rem!important;height:48px!important;transition:all .2s!important;box-shadow:inset 0 2px 4px rgba(0,0,0,0.3)!important;backdrop-filter:blur(10px)!important;}
.stTextInput>div>div>input:focus{border-color:var(--brand)!important;box-shadow:0 0 0 2px rgba(59,130,246,0.3), inset 0 2px 4px rgba(0,0,0,0.3)!important;outline:none!important;}
.stTextInput>div>div>input::placeholder{color:var(--text-xmuted)!important;}
.stTextInput label,.stTextInput [data-testid="InputInstructions"]{color:var(--text-secondary)!important;}

.stSelectbox [data-baseweb="select"]>div,.stMultiSelect [data-baseweb="select"]>div{border-radius:var(--r-md)!important;border:1px solid var(--border-dark)!important;background:var(--surface)!important;min-height:48px!important;color:var(--text-primary)!important;box-shadow:inset 0 2px 4px rgba(0,0,0,0.3)!important;backdrop-filter:blur(10px)!important;}
.stSelectbox [data-baseweb="select"] [data-testid="stMarkdown"],.stSelectbox [data-baseweb="select"] span,.stSelectbox [data-baseweb="select"] div,.stSelectbox [data-baseweb="select"] p{color:var(--text-primary)!important;}
.stSelectbox label,.stMultiSelect label{color:var(--text-secondary)!important;}
[data-baseweb="menu"] li,[data-baseweb="menu"] [role="option"]{color:var(--text-primary)!important;background:rgba(30,41,59,0.95)!important;font-family:var(--font-body)!important;}
[data-baseweb="menu"] li:hover,[data-baseweb="menu"] [aria-selected="true"]{background:var(--brand-light)!important;color:var(--brand-dark)!important;}
.stMultiSelect span[data-baseweb="tag"]{background:var(--brand-light)!important;color:var(--brand-dark)!important;border:1px solid var(--brand-border)!important;font-family:var(--font-mono)!important;font-size:.75rem!important;border-radius:5px!important;}
.stMultiSelect span[data-baseweb="tag"] span{color:var(--brand-dark)!important;}

.stCheckbox label,.stCheckbox span,.stCheckbox p{color:var(--text-secondary)!important;}

div[data-testid="stExpander"]{background:var(--surface)!important;border:1px solid var(--border-light)!important;border-radius:var(--r-lg)!important;box-shadow:var(--shadow-sm)!important;margin-bottom:var(--sp-2)!important;backdrop-filter:blur(12px)!important;}
div[data-testid="stExpander"] summary{font-family:var(--font-body)!important;font-size:.9rem!important;font-weight:500!important;color:var(--text-secondary)!important;padding:var(--sp-4) var(--sp-5)!important;}
div[data-testid="stExpander"] summary *{color:inherit!important;}
div[data-testid="stExpander"] [data-testid="stExpanderDetails"]{color:var(--text-primary)!important;}

[data-testid="stMetric"]{background:var(--surface)!important;border:1px solid var(--border-light)!important;border-radius:var(--r-xl)!important;padding:var(--sp-5)!important;box-shadow:var(--shadow-md)!important;backdrop-filter:blur(12px)!important;}
[data-testid="stMetricLabel"]{font-family:var(--font-mono)!important;font-size:.7rem!important;color:var(--text-muted)!important;text-transform:uppercase!important;letter-spacing:.1em!important;font-weight:600!important;}
[data-testid="stMetricLabel"] *{color:var(--text-muted)!important;}
[data-testid="stMetricValue"]{font-size:1.875rem!important;color:var(--text-primary)!important;font-weight:700!important;letter-spacing:-.02em!important;text-shadow:0 2px 4px rgba(0,0,0,0.5)!important;}
[data-testid="stMetricValue"] *{color:var(--text-primary)!important;}

[data-testid="stAlert"] *{color:inherit!important;}
.stSpinner>div{border-color:var(--brand) transparent transparent transparent!important;}
.stSpinner p,.stSpinner span{color:var(--text-muted)!important;}
.stCode{border-radius:var(--r-lg)!important;}
pre,code{font-family:var(--font-mono)!important;color:var(--text-secondary)!important;white-space:pre-wrap!important;word-break:break-word!important;}
[data-testid="stCodeBlock"] pre,[data-testid="stCodeBlock"] code{background:var(--surface-sub2)!important;color:var(--text-primary)!important;}
.stJson *{color:var(--text-primary)!important;}
</style>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# UTILITY FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def load_vcf(filename):
    """Load a VCF file from sample_data/. Falls back to get_sample_vcf() if not found."""
    p = os.path.join(BASE_DIR, "sample_data", filename)
    if os.path.exists(p):
        with open(p) as f:
            return f.read()
    # FIX: graceful fallback instead of silent empty string
    raise FileNotFoundError(f"Sample file not found: {p}")

def run_pipeline(vcf, drugs, pid, key, run_ix=True, gen_pdf=True, skip_llm=False):
    parsed  = parse_vcf(vcf)
    results = run_risk_assessment(parsed, drugs)
    results = generate_all_explanations(key, results, skip_llm=skip_llm)
    outputs = [build_output_schema(patient_id=pid, drug=r["drug"], result=r,
                parsed_vcf=parsed, llm_exp=r.get("llm_explanation", {})) for r in results]
    ix  = run_interaction_analysis(drugs, results) if run_ix and len(drugs) > 1 else None
    pdf = None
    if gen_pdf:
        try:
            pdf = generate_pdf_report(pid, outputs, parsed)
        except Exception:
            pass
    return parsed, results, outputs, ix, pdf

def func_cls(status):
    s = (status or "").lower()
    if "no_function" in s or "no function" in s:
        return "v-nofunc"
    if any(x in s for x in ["decreased","splice","missense","frame","stop","pathogenic"]) and "synonymous" not in s:
        return "v-dec"
    return "v-norm"

def sec(label):
    st.markdown(f'<div class="sec-label">{label}</div>', unsafe_allow_html=True)

def risk_badge_html(rl):
    rc = RISK_CFG.get(rl, RISK_CFG["Unknown"])
    return (f'<span class="risk-badge" style="background:{rc["tag_bg"]};color:{rc["tag_text"]};'
            f'border-color:{rc["border"]};">'
            f'<span style="font-size:.8rem;">{rc["shape"]}</span>{rl}</span>')

def clean_model_label(raw_model: str):
    is_static = "static" in raw_model.lower()
    if is_static:
        return "Static Template", True
    clean = re.sub(r"\s*\(.*?\)$", "", raw_model).strip()
    return (clean or raw_model), False


# ══════════════════════════════════════════════════════════════════════════════
# VISUAL COMPONENTS (unchanged from v9.2)
# ══════════════════════════════════════════════════════════════════════════════

def compute_pgx(outputs):
    SEV_S  = {"none": 0, "low": 20, "moderate": 45, "high": 70, "critical": 100}
    RISK_S = {"Safe": 0, "Adjust Dosage": 35, "Toxic": 85, "Ineffective": 70, "Unknown": 20}
    W      = {"FLUOROURACIL": 1.4, "AZATHIOPRINE": 1.3, "CLOPIDOGREL": 1.3,
               "WARFARIN": 1.2, "CODEINE": 1.1, "SIMVASTATIN": 1.0}
    if not outputs:
        return 0, "No data", []
    tw = ws = 0
    bd = []
    for o in outputs:
        drug = o["drug"]
        sev  = o["risk_assessment"]["severity"]
        rl   = o["risk_assessment"]["risk_label"]
        gene = o["pharmacogenomic_profile"]["primary_gene"]
        ph   = o["pharmacogenomic_profile"]["phenotype"]
        sc   = (SEV_S.get(sev, 0) + RISK_S.get(rl, 0)) / 2
        wt   = W.get(drug, 1.0)
        ws  += sc * wt
        tw  += wt
        bd.append((gene, drug, ph, rl, sc))
    final = min(100, int(ws / tw)) if tw else 0
    labels = ["Low Risk", "Moderate Risk", "High Risk", "Very High Risk", "Critical Risk"]
    label  = labels[min(4, final // 20)]
    return final, label, bd


def render_pgx(outputs):
    score, label, bd = compute_pgx(outputs)
    SCORE_COLORS = ["#16A34A", "#D97706", "#EA580C", "#DC2626", "#B91C1C"]
    color = SCORE_COLORS[min(4, score // 20)]
    pills = ""
    for gene, _, ph, rl, _ in bd:
        rc = RISK_CFG.get(rl, RISK_CFG["Unknown"])
        pills += (f'<span class="pgx-pill" style="background:{rc["tag_bg"]};border-color:{rc["border"]};'
                  f'color:{rc["tag_text"]};">{gene} · {ph}</span>')
    st.markdown(f"""
    <div class="pgx-card">
      <div class="pgx-eyebrow">Polygenic Risk Score</div>
      <div class="pgx-score" style="color:{color};">{score}</div>
      <div class="pgx-label">{label} — composite across {len(outputs)} drug{"s" if len(outputs)!=1 else ""}</div>
      <div class="pgx-marker">
        <div class="pgx-fill" style="width:{score}%;background:linear-gradient(90deg,{color}99,{color});"></div>
        <div class="pgx-indicator" style="left:{score}%;border-color:{color};"></div>
      </div>
      <div class="pgx-thresh-labels">
        <span>0 — No Risk</span><span>25</span><span>50 — High</span><span>75</span><span>100 — Critical</span>
      </div>
      <div class="pgx-pills">{pills}</div>
    </div>""", unsafe_allow_html=True)


def render_risk_center(outputs, parsed):
    sev = max((o["risk_assessment"]["severity"] for o in outputs),
              key=lambda s: SEV_RANK.get(s, 0), default="none")
    sp  = SEV_CFG.get(sev, SEV_CFG["none"])
    EMO = {"none": "✓", "low": "△", "moderate": "⚠", "high": "⛔", "critical": "🚨"}
    hc  = sum(1 for o in outputs if o["risk_assessment"]["severity"] in ("high", "critical"))
    st.markdown(f"""
    <div class="risk-center" style="background:{sp['bg']};border-color:{sp['border']};color:{sp['text']};">
      <div class="rc-eyebrow">Risk Command Center</div>
      <div class="rc-headline">{EMO.get(sev,"")} {sp['label']} Risk Profile</div>
      <div class="rc-sub">Patient pharmacogenomic assessment across {len(outputs)} medication{"s" if len(outputs)!=1 else ""}</div>
      <div class="rc-stats" style="border-top-color:{sp['border']}88;">
        <div><div class="rc-stat-num">{len(outputs)}</div><div class="rc-stat-lbl">Drugs Assessed</div></div>
        <div><div class="rc-stat-num" style="{'color:#B91C1C' if hc else ''}">{hc}</div><div class="rc-stat-lbl">High / Critical</div></div>
        <div><div class="rc-stat-num">{len(parsed.get('detected_genes',[]))}</div><div class="rc-stat-lbl">Genes Detected</div></div>
        <div><div class="rc-stat-num">{parsed.get('total_variants',0)}</div><div class="rc-stat-lbl">Variants Found</div></div>
      </div>
    </div>""", unsafe_allow_html=True)


def render_critical_alerts(outputs):
    for o in outputs:
        if o["risk_assessment"]["severity"] == "critical":
            drug = o["drug"]
            note = o["clinical_recommendation"]["dosing_recommendation"][:240]
            st.markdown(f"""
            <div class="crit-alert">
              <div style="font-size:1.25rem;flex-shrink:0;padding-top:1px;">🚨</div>
              <div>
                <div class="crit-title">Critical Safety Alert — {drug}</div>
                <div class="crit-note">{note}{"…" if len(o["clinical_recommendation"]["dosing_recommendation"])>240 else ""}</div>
                <div class="crit-action">⚡ Contact prescribing physician immediately</div>
              </div>
            </div>""", unsafe_allow_html=True)


def render_disclaimer():
    st.markdown("""
    <div class="disclaimer-box">
      <span style="font-size:1rem;flex-shrink:0;">📋</span>
      <div class="disclaimer-text">
        <strong>For informational purposes only.</strong> These results require review by a qualified
        clinical pharmacologist or geneticist before any medication changes. All recommendations are
        based on CPIC Level A evidence — verify at <strong>cpicpgx.org</strong>.
      </div>
    </div>""", unsafe_allow_html=True)


def render_gene_row(outputs):
    GENE_ORDER = ["CYP2D6", "CYP2C19", "CYP2C9", "SLCO1B1", "TPMT", "DPYD"]
    gp = {o["pharmacogenomic_profile"]["primary_gene"]: o["pharmacogenomic_profile"]["phenotype"]
          for o in outputs}
    boxes = ""
    for g in GENE_ORDER:
        ph = gp.get(g, "Unknown")
        pc = PHENO_CFG.get(ph, PHENO_CFG["Unknown"])
        bar = min(100, pc["pct"])
        active = "active" if g in gp else ""
        boxes += f"""
        <div class="gene-box {active}" style="{'border-color:'+pc['border']+';' if active else ''}">
          <div class="gene-nm" style="{'color:'+pc['text']+';' if active else ''}">{g}</div>
          <div class="gene-track">
            <div class="gene-fill" style="width:{bar}%;background:{pc['bar']};"></div>
          </div>
          <div class="gene-ph" style="color:{pc['text'] if active else 'var(--text-xmuted)'};">{ph}</div>
        </div>"""
    sec("Gene Activity Overview")
    st.markdown(f'<div class="gene-row">{boxes}</div>', unsafe_allow_html=True)


def render_drug_table(outputs, pid):
    rows = ""
    data = []
    for o in outputs:
        drug = o["drug"]
        rl   = o["risk_assessment"]["risk_label"]
        sev  = o["risk_assessment"]["severity"]
        conf = o["risk_assessment"]["confidence_score"]
        gene = o["pharmacogenomic_profile"]["primary_gene"]
        ph   = o["pharmacogenomic_profile"]["phenotype"]
        rc   = RISK_CFG.get(rl, RISK_CFG["Unknown"])
        sp   = SEV_CFG.get(sev, SEV_CFG["none"])
        badge = risk_badge_html(rl)
        rows += f"""<div class="dtab-row">
          <div class="dtab-cell" style="font-weight:700;color:#0F172A;">{drug.title()}</div>
          <div class="dtab-cell">{badge}</div>
          <div class="dtab-cell"><span style="color:{sp['text']};font-weight:600;">{sp['label']}</span></div>
          <div class="dtab-cell" style="font-family:var(--font-mono);font-size:.8rem;color:#64748B;">{gene}</div>
          <div class="dtab-cell"><span style="font-family:var(--font-mono);font-size:.8rem;color:{rc['tag_text']};background:{rc['tag_bg']};border:1px solid {rc['border']};padding:2px 8px;border-radius:4px;font-weight:600;">{ph}</span></div>
          <div class="dtab-cell">
            <div style="flex:1;height:4px;background:#E8EDF5;border-radius:2px;overflow:hidden;margin-right:8px;">
              <div style="width:{conf*100:.0f}%;height:100%;background:{rc['severity_dot']};border-radius:2px;"></div>
            </div>
            <span style="font-family:var(--font-mono);font-size:.75rem;color:#64748B;font-weight:600;">{conf:.0%}</span>
          </div>
        </div>"""
        data.append({"Drug": drug, "Risk": rl, "Severity": sev, "Gene": gene,
                      "Phenotype": ph, "Confidence": f"{conf:.0%}"})
    sec("Drug Risk Summary")
    st.markdown(f"""
    <div class="dtab">
      <div class="dtab-head">
        <div class="dtab-hcell">Drug</div><div class="dtab-hcell">Risk Label</div>
        <div class="dtab-hcell">Severity</div><div class="dtab-hcell">Gene</div>
        <div class="dtab-hcell">Phenotype</div><div class="dtab-hcell">Confidence</div>
      </div>{rows}
    </div>""", unsafe_allow_html=True)
    df = pd.DataFrame(data)
    st.download_button("⬇ Download CSV", data=df.to_csv(index=False),
        file_name=f"SurakshaRx_{pid}.csv", mime="text/csv", key=f"csv_{pid}")


def render_heatmap(outputs):
    DRUG_ORD = ["CODEINE","WARFARIN","CLOPIDOGREL","SIMVASTATIN","AZATHIOPRINE","FLUOROURACIL"]
    GENE_ORD = ["CYP2D6","CYP2C9","CYP2C19","SLCO1B1","TPMT","DPYD"]
    DG = {"CODEINE":"CYP2D6","WARFARIN":"CYP2C9","CLOPIDOGREL":"CYP2C19",
          "SIMVASTATIN":"SLCO1B1","AZATHIOPRINE":"TPMT","FLUOROURACIL":"DPYD"}
    rmap  = {o["drug"]: o for o in outputs}
    drugs = [d for d in DRUG_ORD if d in rmap]
    if not drugs:
        return
    n = len(drugs)
    hdrs = '<div class="hm-header"></div>'
    for d in drugs:
        hdrs += f'<div class="hm-header">{d[:5]}</div>'
    rows = ""
    for gene in GENE_ORD:
        rows += f'<div class="hm-header" style="justify-content:flex-end;padding-right:6px;">{gene}</div>'
        for d in drugs:
            if DG.get(d) == gene and d in rmap:
                o  = rmap[d]
                rl = o["risk_assessment"]["risk_label"]
                ph = o["pharmacogenomic_profile"]["phenotype"]
                rc = RISK_CFG.get(rl, RISK_CFG["Unknown"])
                sh = {"Adjust Dosage":"Adjust","Ineffective":"Ineffect.","Unknown":"?"}.get(rl, rl)
                rows += (f'<div class="hm-cell" style="background:{rc["bg"]};border-color:{rc["border"]};" '
                         f'title="{d}×{gene}: {rl} ({ph})">'
                         f'<div class="hm-cell-name" style="color:{rc["text"]};">{sh}</div>'
                         f'<div class="hm-cell-risk" style="color:{rc["text"]};">{ph}</div></div>')
            else:
                rows += '<div class="hm-cell" style="background:#F1F5F9;border-color:#E8EDF5;"><div class="hm-cell-risk" style="color:#94A3B8;">—</div></div>'
    legend = "".join(
        f'<div class="hm-legend-item"><span class="hm-dot" style="background:{RISK_CFG[r]["bg"]};border-color:{RISK_CFG[r]["border"]};"></span><span>{RISK_CFG[r]["shape"]} {r}</span></div>'
        for r in ["Safe", "Adjust Dosage", "Toxic", "Ineffective"])
    st.markdown(f"""
    <div class="hm-wrap">
      <div class="hm-eyebrow">Drug × Gene Risk Matrix</div>
      <div class="hm-grid" style="grid-template-columns:80px repeat({n},1fr);">{hdrs}{rows}</div>
      <div class="hm-legend">{legend}</div>
    </div>""", unsafe_allow_html=True)


def render_chromosome(outputs, parsed):
    det  = set(parsed.get("detected_genes", []))
    rmap = {o["pharmacogenomic_profile"]["primary_gene"]: o for o in outputs}
    rows = ""
    for gene, info in CHROM_INFO.items():
        ch  = info["chrom"]
        pos = info["pos_mb"]
        pct = (pos / CHROM_LEN.get(ch, 200)) * 100
        if gene in rmap:
            rl = rmap[gene]["risk_assessment"]["risk_label"]
            mc = RISK_CFG.get(rl, RISK_CFG["Unknown"])["severity_dot"]
        elif gene in det:
            mc = "#94A3B8"
        else:
            mc = "#DDE3EE"
        rows += f"""<div class="chrom-row">
          <div class="chrom-chr">{ch}</div>
          <div class="chrom-bar">
            <div class="chrom-body"></div>
            <div class="chrom-marker" style="left:{pct}%;background:{mc};box-shadow:0 0 5px {mc}88;"></div>
          </div>
          <div class="chrom-gene">{gene}</div>
          <div class="chrom-band">{info['band']}</div>
        </div>"""
    st.markdown(f"""
    <div class="chrom-wrap">
      <div class="chrom-eyebrow">Variant Chromosome Locations</div>
      {rows}
      <div style="font-family:var(--font-mono);font-size:.65rem;color:#94A3B8;margin-top:var(--sp-3);">
        Coloured markers = variants detected · Grey = undetected
      </div>
    </div>""", unsafe_allow_html=True)


def render_pop_freq(gene, ph):
    freq = POP_FREQ.get(gene, {})
    if not freq:
        return
    rows = ""
    for p, pct in sorted(freq.items(), key=lambda x: -x[1]):
        you = (p == ph)
        pc  = PHENO_CFG.get(p, PHENO_CFG["Unknown"])
        you_tag = f'<span class="pop-you">← You</span>' if you else ""
        w = "font-weight:700;" if you else ""
        rows += f"""<div class="pop-row">
          <div class="pop-ph" style="{w}{'color:'+pc['text']+';' if you else ''}">{pc['label']}</div>
          <div class="pop-track"><div class="pop-fill" style="width:{min(pct,100)}%;background:{pc['bar'] if you else '#CBD5E1'};"></div></div>
          <div class="pop-pct" style="{w}{'color:'+pc['text']+';' if you else ''}">{pct}%{you_tag}</div>
        </div>"""
    st.markdown(f"""
    <div class="pop-wrap">
      <div class="pop-eyebrow">{gene} — Population Distribution</div>{rows}
    </div>""", unsafe_allow_html=True)


def render_ix_matrix(outputs, ix):
    if not ix or len(outputs) < 2:
        return
    drugs = [o["drug"] for o in outputs]
    n     = len(drugs)
    sm    = {}
    for x in ix.get("all_interactions", []):
        inv = x.get("drugs_involved", [])
        if len(inv) == 2:
            sv = x.get("severity", "none")
            sm[(inv[0], inv[1])] = sm[(inv[1], inv[0])] = sv
    MC = {
        "critical": {"bg":"#FEF2F2","text":"#7F1D1D","border":"#FECACA"},
        "high":     {"bg":"#FEF2F2","text":"#7F1D1D","border":"#FECACA"},
        "moderate": {"bg":"#FFFBEB","text":"#78350F","border":"#FDE68A"},
        "low":      {"bg":"#FEFCE8","text":"#713F12","border":"#FDE047"},
        "none":     {"bg":"#F0FDF4","text":"#14532D","border":"#BBF7D0"},
        "diag":     {"bg":"#F1F5F9","text":"#64748B","border":"#E2E8F0"},
    }
    hdrs = '<div class="ix-head"></div>'
    for d in drugs:
        hdrs += f'<div class="ix-head">{d[:6]}</div>'
    grid = ""
    for i, d1 in enumerate(drugs):
        grid += f'<div class="ix-head" style="justify-content:flex-end;padding-right:4px;">{d1[:6]}</div>'
        for j, d2 in enumerate(drugs):
            if i == j:
                mc = MC["diag"]
                grid += f'<div class="ix-cell" style="background:{mc["bg"]};border-color:{mc["border"]};color:{mc["text"]};">—</div>'
            else:
                sv = sm.get((d1, d2), "none")
                mc = MC.get(sv, MC["none"])
                lbl = sv.upper() if sv != "none" else "OK"
                grid += f'<div class="ix-cell" style="background:{mc["bg"]};border-color:{mc["border"]};color:{mc["text"]};">{lbl}</div>'
    sec("Drug Interaction Matrix")
    st.markdown(f"""
    <div style="background:#FFFFFF;border:1px solid #E8EDF5;border-radius:var(--r-xl);padding:var(--sp-5);margin-bottom:var(--sp-4);box-shadow:var(--shadow-sm);">
      <div class="ix-grid" style="grid-template-columns:76px repeat({n},1fr);gap:3px;">{hdrs}{grid}</div>
    </div>""", unsafe_allow_html=True)
    shown = set()
    for x in ix.get("all_interactions", []):
        inv = x.get("drugs_involved", [])
        key = tuple(sorted(inv))
        if len(inv) == 2 and key not in shown:
            shown.add(key)
            sv = x.get("severity", "low")
            sp = SEV_CFG.get(sv, SEV_CFG["low"])
            with st.expander(f"{' + '.join(inv)}  —  {sv.upper()} interaction"):
                mech = x.get("mechanism", x.get("message", ""))
                rec  = x.get("recommendation", "")
                if mech:
                    st.markdown(f'<div style="font-size:.9rem;color:#334155;line-height:1.75;margin-bottom:var(--sp-2);">{mech}</div>', unsafe_allow_html=True)
                if rec:
                    st.markdown(f'<div style="font-family:var(--font-mono);font-size:.8rem;color:{sp["text"]};margin-top:var(--sp-2);font-weight:600;">→ {rec}</div>', unsafe_allow_html=True)


def render_narrative(outputs, parsed, pid, key, skip_llm):
    results_for = [{"drug": o["drug"],
                    "primary_gene": o["pharmacogenomic_profile"]["primary_gene"],
                    "phenotype": o["pharmacogenomic_profile"]["phenotype"],
                    "risk_label": o["risk_assessment"]["risk_label"],
                    "severity": o["risk_assessment"]["severity"]}
                   for o in outputs]
    with st.spinner("Generating AI clinical summary…"):
        nar = generate_patient_narrative(pid, results_for, parsed, key, skip_llm)
    model_label = "Static Template" if (skip_llm or not key) else "LLaMA 3.3 70B"
    sec("AI Clinical Summary")
    st.markdown(f"""
    <div class="narrative-box">
      <div class="narrative-header">
        <span class="ai-badge-pill">{model_label}</span>
        <span style="font-size:.9rem;font-weight:600;color:#1E3A8A;">Unified Patient Summary</span>
      </div>
      <div class="narrative-text">{nar}</div>
    </div>""", unsafe_allow_html=True)


def render_before_after(outputs):
    bad = [o for o in outputs if o["risk_assessment"]["risk_label"] in ("Toxic", "Ineffective")]
    if not bad:
        return
    o    = bad[0]
    drug = o["drug"]
    rl   = o["risk_assessment"]["risk_label"]
    alts = o["clinical_recommendation"].get("alternative_drugs", [])
    alt  = alts[0] if alts else "Alternative medication"
    gene = o["pharmacogenomic_profile"]["primary_gene"]
    ph   = o["pharmacogenomic_profile"]["phenotype"]
    BEFORE = {
        "Toxic":      f"Standard {drug.lower()} dose → toxic accumulation → serious harm",
        "Ineffective": f"Standard {drug.lower()} dose → zero therapeutic effect → treatment failure",
    }
    sec("Clinical Impact — Before & After PGx")
    st.markdown(f"""
    <div class="ba-grid">
      <div class="ba-side" style="background:#FFF1F2;border-right:1px solid #E8EDF5;">
        <div class="ba-side-lbl" style="color:#B91C1C;">⛔ Without SurakshaRx</div>
        <div class="ba-drug" style="color:#7F1D1D;">{drug.title()} — Standard Protocol</div>
        <div class="ba-text" style="color:#7F1D1D;">{BEFORE.get(rl,"Risk undetected")}</div>
        <div class="ba-gene" style="color:#FECACA;">{gene} {ph} phenotype undetected</div>
      </div>
      <div class="ba-side" style="background:#F0FDF4;">
        <div class="ba-side-lbl" style="color:#14532D;">✓ With SurakshaRx</div>
        <div class="ba-drug" style="color:#15803D;">{alt} — PGx-Guided Protocol</div>
        <div class="ba-text" style="color:#16A34A;">Appropriate alternative selected → safe, effective therapy</div>
        <div class="ba-gene" style="color:#BBF7D0;">{gene} {ph} phenotype identified → therapy optimised</div>
      </div>
    </div>""", unsafe_allow_html=True)


def render_rx_checker(outputs):
    sec("Prescription Safety Checker")
    rmap  = {o["drug"]: o for o in outputs}
    drugs = [o["drug"] for o in outputs]
    c1, c2 = st.columns([2, 1])
    with c1:
        sel = st.selectbox("Select drug", drugs,
              format_func=lambda x: f"{x.title()}  ({GENE_DRUG_MAP.get(x,'')})",
              key="rx_drug", label_visibility="collapsed")
    with c2:
        check = st.button("Check Safety →", key="rx_check")
    if check and sel in rmap:
        o    = rmap[sel]
        rl   = o["risk_assessment"]["risk_label"]
        sev  = o["risk_assessment"]["severity"]
        rec  = o["clinical_recommendation"]["dosing_recommendation"]
        gene = o["pharmacogenomic_profile"]["primary_gene"]
        ph   = o["pharmacogenomic_profile"]["phenotype"]
        rc   = RISK_CFG.get(rl, RISK_CFG["Unknown"])
        sp   = SEV_CFG.get(sev, SEV_CFG["none"])
        VERDICT = {
            "Safe":         "✓ Safe to Prescribe",
            "Adjust Dosage":"△ Prescribe with Dose Adjustment",
            "Toxic":        "⛔ Do Not Prescribe — Toxicity Risk",
            "Ineffective":  "◆ Do Not Prescribe — Drug Ineffective",
        }
        st.markdown(f"""
        <div class="rx-result" style="background:{rc['bg']};border-color:{rc['border']};">
          <div class="rx-verdict" style="color:{rc['text']};">{VERDICT.get(rl, rl)}</div>
          <div class="rx-detail">{gene} {ph} phenotype detected. {rec}</div>
          <div class="rx-meta" style="color:{sp['text']};">Severity: {sp['label']} · Confidence: {o["risk_assessment"]["confidence_score"]:.0%} · CPIC Level A</div>
        </div>""", unsafe_allow_html=True)
    elif not check:
        st.markdown("""<div class="info-strip"><span>🔍</span>
          <div class="info-strip-text">Select a drug and click <strong>Check Safety</strong> to validate against this patient's genotype.</div>
        </div>""", unsafe_allow_html=True)


def render_clinical_note(outputs, pid):
    lines = [f"SurakshaRx Clinical Note — Patient {pid} — {datetime.now(timezone.utc).strftime('%Y-%m-%d')}",
             "=" * 60, ""]
    for o in outputs:
        gene = o["pharmacogenomic_profile"]["primary_gene"]
        dip  = o["pharmacogenomic_profile"]["diplotype"]
        ph   = o["pharmacogenomic_profile"]["phenotype"]
        drug = o["drug"]
        rl   = o["risk_assessment"]["risk_label"]
        rec  = o["clinical_recommendation"]["dosing_recommendation"]
        alts = o["clinical_recommendation"].get("alternative_drugs", [])
        lines.append(f"DRUG: {drug}")
        lines.append(f"Gene: {gene} | Diplotype: {dip} | Phenotype: {ph} | Risk: {rl}")
        lines.append(f"CPIC: {rec}")
        if alts:
            lines.append(f"Alternatives: {', '.join(alts)}")
        lines.append("")
    lines += ["", "-" * 60,
              "Generated by SurakshaRx v10.0 · CPIC Level A evidence · cpicpgx.org",
              "NOT FOR CLINICAL USE WITHOUT VALIDATION BY A QUALIFIED CLINICIAN"]
    note = "\n".join(lines)
    sec("One-Click Clinical Note")
    st.markdown(f'<div class="note-box"><pre>{note}</pre></div>', unsafe_allow_html=True)
    st.download_button("⬇ Download Clinical Note", data=note,
        file_name=f"clinical_note_{pid}.txt", mime="text/plain", key=f"note_{pid}")


def render_patient_mode(outputs):
    bad = any(o["risk_assessment"]["risk_label"] in ("Toxic","Ineffective") for o in outputs)
    if bad:
        st.markdown("""<div class="patient-banner" style="background:#FFF1F2;border-color:#FECACA;">
          <div class="patient-banner-title" style="color:#B91C1C;">🚨 Important — Some medications need urgent attention</div>
          <div class="patient-banner-sub" style="color:#7F1D1D;">Your genetic results show that one or more medications may not be safe or effective for you. Please speak with your doctor before taking these medications.</div>
        </div>""", unsafe_allow_html=True)
    else:
        st.markdown("""<div class="patient-banner" style="background:#F0FDF4;border-color:#BBF7D0;">
          <div class="patient-banner-title" style="color:#14532D;">✓ Good news — Your medications look safe</div>
          <div class="patient-banner-sub" style="color:#16A34A;">Based on your genetic profile, the medications reviewed are predicted to work normally at standard doses.</div>
        </div>""", unsafe_allow_html=True)
    for o in outputs:
        drug    = o["drug"]
        rl      = o["risk_assessment"]["risk_label"]
        gene    = o["pharmacogenomic_profile"]["primary_gene"]
        ph      = o["pharmacogenomic_profile"]["phenotype"]
        alts    = o["clinical_recommendation"].get("alternative_drugs", [])
        phplain = PLAIN_PHENO.get(ph, ph)
        explain = PLAIN_RISK.get((drug, ph), "")
        VERDICT = {
            "Safe":         "✓ This medicine is likely safe for you",
            "Adjust Dosage":"△ You may need a different dose",
            "Toxic":        "⛔ This medicine could be harmful to you",
            "Ineffective":  "◆ This medicine likely won't work for you",
        }
        rc = RISK_CFG.get(rl, RISK_CFG["Unknown"])
        action = ""
        if rl in ("Toxic", "Ineffective"):
            alt_text = f"They may suggest: <strong>{', '.join(alts[:3])}</strong>" if alts else "Ask about alternative medications."
            action = f'<div class="pcard-action"><span style="font-size:1rem;">💊</span><div class="pcard-action-text"><strong>Talk to your doctor before taking {drug.title()}.</strong><br>{alt_text}</div></div>'
        elif rl == "Adjust Dosage":
            action = f'<div class="pcard-action"><span style="font-size:1rem;">📋</span><div class="pcard-action-text"><strong>Tell your doctor about this result before starting {drug.title()}.</strong><br>You may need a different dose than usually prescribed.</div></div>'
        st.markdown(f"""
        <div class="pcard" style="border-color:{rc['border']};">
          <div class="pcard-drug">{drug.title()}</div>
          <div class="pcard-verdict" style="color:{rc['text']};">{VERDICT.get(rl, rl)}</div>
          <div class="pcard-gene">{gene} · {phplain}</div>
          {f'<div class="pcard-plain">{explain}</div>' if explain else ''}
          {action}
        </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# MASTER RESULTS RENDERER
# ══════════════════════════════════════════════════════════════════════════════

def render_results(outputs, parsed, ix, pdf_bytes, pid, patient_mode=False, key="", skip_llm=False):
    render_disclaimer()
    render_risk_center(outputs, parsed)
    render_critical_alerts(outputs)

    dc1, dc2, dc3 = st.columns(3)
    with dc1:
        st.download_button("⬇ Download All JSON", data=json.dumps(outputs, indent=2),
            file_name=f"SurakshaRx_{pid}.json", mime="application/json",
            use_container_width=True, key=f"dlall_{pid}")
    with dc2:
        if pdf_bytes:
            st.download_button("⬇ Download PDF Report", data=pdf_bytes,
                file_name=f"SurakshaRx_{pid}.pdf", mime="application/pdf",
                use_container_width=True, key=f"dlpdf_{pid}")
    with dc3:
        if ix and ix.get("interactions_found"):
            st.download_button("⬇ Interactions JSON", data=json.dumps(ix, indent=2),
                file_name=f"SurakshaRx_{pid}_ix.json", mime="application/json",
                use_container_width=True, key=f"dlix_{pid}")

    st.markdown("<div style='height:var(--sp-3)'></div>", unsafe_allow_html=True)

    if patient_mode:
        render_patient_mode(outputs)
        return

    render_gene_row(outputs)
    render_drug_table(outputs, pid)
    render_pgx(outputs)

    c1, c2 = st.columns([1.4, 1], gap="large")
    with c1: render_heatmap(outputs)
    with c2: render_chromosome(outputs, parsed)

    if ix and len(outputs) >= 2:
        render_ix_matrix(outputs, ix)

    render_narrative(outputs, parsed, pid, key, skip_llm)
    render_before_after(outputs)
    render_rx_checker(outputs)
    render_clinical_note(outputs, pid)

    sec("Individual Drug Analysis")
    for output in outputs:
        rl   = output["risk_assessment"]["risk_label"]
        drug = output["drug"]
        sev  = output["risk_assessment"]["severity"]
        conf = output["risk_assessment"]["confidence_score"]
        gene = output["pharmacogenomic_profile"]["primary_gene"]
        dip  = output["pharmacogenomic_profile"]["diplotype"]
        ph   = output["pharmacogenomic_profile"]["phenotype"]
        var  = output["pharmacogenomic_profile"]["detected_variants"]
        rec  = output["clinical_recommendation"]["dosing_recommendation"]
        alts = output["clinical_recommendation"].get("alternative_drugs", [])
        mon  = output["clinical_recommendation"].get("monitoring_required", "")
        exp  = output["llm_generated_explanation"]
        rc   = RISK_CFG.get(rl, RISK_CFG["Unknown"])
        sp   = SEV_CFG.get(sev, SEV_CFG["none"])
        cpic_lv = output.get("pharmacogenomic_profile", {}).get("cpic_evidence_level", "Level A")

        st.markdown(f"""
        <div class="dcard reveal-card">
          <div class="dcard-header">
            <div class="dcard-left">
              <div class="dcard-indicator" style="background:{rc['severity_dot']};box-shadow:0 0 0 3px {rc['bg']};"></div>
              <div>
                <div class="dcard-drug">{drug.title()}
                  <span class="cpic-badge">CPIC {cpic_lv}</span>
                </div>
                <div class="dcard-meta">{gene} · {dip} · {ph}</div>
              </div>
            </div>
            {risk_badge_html(rl)}
          </div>
          <div class="dcard-body">
            <div class="metrics-row">
              <div class="metric-cell"><div class="metric-key">Phenotype</div><div class="metric-val" style="color:{rc['text']};font-size:.95rem;">{ph}</div></div>
              <div class="metric-cell"><div class="metric-key">Severity</div><div class="metric-val" style="color:{sp['text']};font-size:.95rem;">{sp['label']}</div></div>
              <div class="metric-cell"><div class="metric-key">Confidence</div><div class="metric-val">{conf:.0%}</div></div>
              <div class="metric-cell"><div class="metric-key">Variants</div><div class="metric-val">{len(var)}</div></div>
            </div>""", unsafe_allow_html=True)

        dq = min(1.0, len(var) / 3.0)
        st.markdown(f"""
        <div class="conf-grid">
          <div>
            <div class="conf-label"><span>Prediction Confidence</span><span style="color:{rc['severity_dot']};font-weight:700;">{conf:.0%}</span></div>
            <div class="conf-track"><div class="conf-fill" style="width:{conf*100:.1f}%;background:{rc['severity_dot']};"></div></div>
          </div>
          <div>
            <div class="conf-label"><span>Data Quality</span><span style="color:#64748B;">{len(var)} variant{"s" if len(var)!=1 else ""}</span></div>
            <div class="conf-track"><div class="conf-fill" style="width:{dq*100:.1f}%;background:#94A3B8;"></div></div>
          </div>
        </div>""", unsafe_allow_html=True)

        if var:
            rows_html = ""
            for v in var:
                fc = func_cls(v.get("functional_status", ""))
                fn = (v.get("functional_status") or "unknown").replace("_", " ").title()
                rows_html += (f'<tr><td class="v-rsid">{v.get("rsid","—")}</td>'
                              f'<td class="v-star">{v.get("star_allele","—")}</td>'
                              f'<td class="{fc}">{fn}</td></tr>')
            st.markdown(f"""
            <div style="margin-bottom:var(--sp-4);">
              <div class="conf-label" style="margin-bottom:var(--sp-3);">Detected Variants ({len(var)})</div>
              <table class="vtable">
                <thead><tr><th>rsID</th><th>Star Allele</th><th>Functional Status</th></tr></thead>
                <tbody>{rows_html}</tbody>
              </table>
            </div>""", unsafe_allow_html=True)

        st.markdown(f"""
        <div class="rec-box" style="background:{rc['bg']};border-color:{rc['border']};">
          <div class="rec-label" style="color:{rc['text']};">CPIC Recommendation — {drug}</div>
          <div class="rec-text">{rec}</div>
        </div>""", unsafe_allow_html=True)

        if mon:
            st.markdown(f"""
            <div class="rec-box" style="background:#F1F5F9;border-color:#E8EDF5;">
              <div class="rec-label" style="color:#64748B;">🔬 Monitoring Protocol</div>
              <div class="rec-text">{mon}</div>
            </div>""", unsafe_allow_html=True)

        if alts:
            chips = "".join(f'<span class="alt-chip">{a}</span>' for a in alts)
            st.markdown(f"""
            <div style="margin-bottom:var(--sp-4);">
              <div class="conf-label" style="margin-bottom:var(--sp-2);">Alternative Medications</div>
              <div class="alt-chips">{chips}</div>
            </div>""", unsafe_allow_html=True)

        render_pop_freq(gene, ph)

        if exp.get("summary"):
            raw_model = exp.get("model_used", "llama-3.3-70b")
            model, is_static = clean_model_label(raw_model)
            blocks = ""
            for lbl, k in [("Summary","summary"), ("Biological Mechanism","biological_mechanism"),
                           ("Variant Significance","variant_significance"), ("Clinical Implications","clinical_implications")]:
                if exp.get(k):
                    blocks += (f'<div class="ai-section">'
                               f'<div class="ai-sec-label">{lbl}</div>'
                               f'<div class="ai-sec-text">{exp[k]}</div>'
                               f'</div>')
            st.markdown(f"""
            <div class="ai-block">
              <div class="ai-header">
                <span class="ai-badge-pill">{model}</span>
                <span class="ai-title">AI Explanation · {drug}</span>
              </div>{blocks}
            </div>""", unsafe_allow_html=True)

        with st.expander(f"Raw JSON — {drug}"):
            st.json(output)

        st.markdown('', unsafe_allow_html=True)  # spacing


# ══════════════════════════════════════════════════════════════════════════════
# NAVIGATION + LAYOUT
# ══════════════════════════════════════════════════════════════════════════════

def render_nav(key_present):
    model_badge = "LLaMA 3.3 70B" if key_present else "Static Mode"
    model_cls   = "pg-badge-brand" if key_present else "pg-badge-default"
    st.markdown(f"""
    <div class="pg-nav">
      <div class="pg-brand">
        <div class="pg-brand-name">Suraksha<span>Rx</span></div>
        <div class="pg-brand-sub">v10.0 · Precision Clinical · RIFT 2026</div>
      </div>
      <div class="pg-nav-badges">
        <span class="pg-badge pg-badge-brand">CPIC Level A</span>
        <span class="pg-badge {model_cls}">{model_badge}</span>
        <span class="pg-badge pg-badge-default">RIFT 2026</span>
      </div>
    </div>
    <div class="trust-strip">
      <div class="trust-item">🔒 Genetic data analyzed locally — never stored</div>
      <div class="trust-sep"></div>
      <div class="trust-item">✓ CPIC Level A Evidence Guidelines</div>
      <div class="trust-sep"></div>
      <div class="trust-item">🧬 6 Pharmacogenes · 6 High-Risk Drugs</div>
      <div class="trust-sep"></div>
      <div class="trust-item">⚕ For review by qualified clinicians only</div>
    </div>""", unsafe_allow_html=True)


def render_steps(has_vcf, has_drugs, has_results):
    steps = [
        ("01", "Upload VCF", has_vcf),
        ("02", "Select Drugs", has_drugs),
        ("03", "Run Analysis", has_results),
        ("04", "Review Results", has_results),
    ]
    html = '<div class="steps">'
    for num, lbl, done in steps:
        cls = "step done" if done else "step"
        html += f'<div class="{cls}"><div class="step-num">{num}</div><div class="step-lbl">{lbl}</div></div>'
    html += "</div>"
    st.markdown(html, unsafe_allow_html=True)


def render_persona_demo(key):
    SEV_COLORS = {
        "critical": {"sev_bg":"#FEF2F2","sev_border":"#FECACA","sev_text":"#7F1D1D","sev_label":"Critical"},
        "high":     {"sev_bg":"#FFF7ED","sev_border":"#FED7AA","sev_text":"#7C2D12","sev_label":"High Risk"},
        "moderate": {"sev_bg":"#FFFBEB","sev_border":"#FDE68A","sev_text":"#78350F","sev_label":"Moderate"},
        "none":     {"sev_bg":"#F0FDF4","sev_border":"#BBF7D0","sev_text":"#14532D","sev_label":"All Safe"},
    }
    st.markdown('<div style="font-size:.8rem;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:#64748B;margin-bottom:var(--sp-3);">Quick Demo — Select Patient Persona</div>', unsafe_allow_html=True)
    cols = st.columns(4)
    for i, (pid, p) in enumerate(PERSONAS.items()):
        sc = SEV_COLORS.get(p["sev"], SEV_COLORS["none"])
        with cols[i]:
            st.markdown(f"""
            <div class="persona-card">
              <div class="pc-sev" style="background:{sc['sev_bg']};border-color:{sc['sev_border']};color:{sc['sev_text']};">
                {sc['sev_label']}
              </div>
              <div class="pc-name">{p['label']}</div>
              <div class="pc-desc">{p['desc']}</div>
            </div>""", unsafe_allow_html=True)
            if st.button(f"Load {p['label']}", key=f"persona_{pid}", use_container_width=True):
                # FIX: try/except fallback if sample_data files are missing
                try:
                    vcf = load_vcf(p["file"])
                except FileNotFoundError:
                    vcf = get_sample_vcf()
                    st.warning(f"Sample file '{p['file']}' not found — using default VCF for demo.")
                pid_gen = f"PG-{uuid.uuid4().hex[:8].upper()}"
                with st.spinner(f"Running {p['label']} analysis…"):
                    parsed, results, outputs, ix, pdf = run_pipeline(
                        vcf, p["drugs"], pid_gen, key, skip_llm=not bool(key))
                st.session_state["results"]      = outputs
                st.session_state["parsed"]       = parsed
                st.session_state["ix"]           = ix
                st.session_state["pdf"]          = pdf
                st.session_state["patient_id"]   = pid_gen
                st.session_state["results_key"]  = key
                st.session_state["results_skip"] = not bool(key)
                st.rerun()


def render_test_suite(key):
    st.markdown("### Test Suite")
    st.markdown(
        '<div style="font-size:.85rem;color:#64748B;margin-bottom:16px;">'
        'Tests run with static templates — no API key needed. '
        'Results load in the Analysis tab after running.</div>',
        unsafe_allow_html=True
    )

    # Show persistent test results from session state
    if "tc_results" in st.session_state:
        for tc_res in st.session_state["tc_results"]:
            color_cls = "tc-status-pass" if tc_res["passed"] else "tc-status-fail"
            icon = "✓ PASS" if tc_res["passed"] else "✗ FAIL"
            st.markdown(
                f'<div class="{color_cls}">'
                f'<strong>{icon} — {tc_res["name"]}</strong><br>'
                f'<span style="font-size:.78rem;opacity:.85;">{tc_res["detail"]}</span><br>'
                f'<span style="font-size:.7rem;opacity:.55;">{tc_res["source"]}</span>'
                f'</div>',
                unsafe_allow_html=True
            )
        st.markdown('<div style="height:8px;"></div>', unsafe_allow_html=True)
        if st.button("Clear results", key="tc_clear"):
            del st.session_state["tc_results"]
            st.rerun()
        st.divider()

    for i, tc in enumerate(TEST_SUITE):
        with st.container():
            st.markdown(f"""
            <div class="tc-card">
              <span class="tc-name">{tc['name']}</span>
              <span class="tc-desc">{tc['desc']}</span>
            </div>""", unsafe_allow_html=True)

            if st.button(f"▶ Run: {tc['name']}", key=f"tc_{i}", use_container_width=True):
                try:
                    vcf = load_vcf(tc["file"])
                    file_source = f"sample_data/{tc['file']}"
                except FileNotFoundError:
                    vcf = get_sample_vcf()
                    file_source = "fallback VCF (sample file not found)"

                pid = f"TC-{uuid.uuid4().hex[:6].upper()}"
                with st.spinner(f"Running {tc['name']}…"):
                    try:
                        parsed, results, outputs, ix, pdf = run_pipeline(
                            vcf, tc["drugs"], pid, key, skip_llm=True)

                        detail_lines = []
                        all_pass = True
                        for o in outputs:
                            drug = o["drug"]
                            got  = o["risk_assessment"]["risk_label"]
                            want = tc["expected"].get(drug)
                            if want is None:
                                continue
                            ok = got == want
                            if not ok:
                                all_pass = False
                            icon = "✓" if ok else "✗"
                            detail_lines.append(f"{icon} {drug}: {got} (expected {want})")

                        # Store result persistently in session_state
                        tc_results = st.session_state.get("tc_results", [])
                        # Remove old result for same test name
                        tc_results = [r for r in tc_results if r["name"] != tc["name"]]
                        tc_results.insert(0, {
                            "name":   tc["name"],
                            "passed": all_pass,
                            "detail": "  ·  ".join(detail_lines),
                            "source": file_source,
                        })
                        st.session_state["tc_results"] = tc_results

                        # Store pipeline results for Analysis tab
                        st.session_state["results"]      = outputs
                        st.session_state["parsed"]       = parsed
                        st.session_state["ix"]           = ix
                        st.session_state["pdf"]          = pdf
                        st.session_state["patient_id"]   = pid
                        st.session_state["results_key"]  = key
                        st.session_state["results_skip"] = True
                        st.rerun()

                    except Exception as run_err:
                        st.error(f"Pipeline error: {run_err}")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    # ── Sidebar ───────────────────────────────────────────────────────────────
    with st.sidebar:
        st.markdown("### ⚙ Settings")
        groq_key  = st.text_input("Groq API Key", type="password",
                        placeholder="gsk_…", help="Required for LLM explanations")
        model_sel = st.selectbox("Model", ["LLaMA 3.3 70B Versatile"],
                        help="Model used for AI explanations")
        use_static = st.checkbox("Test mode: instant (no API call)", value=not bool(groq_key))
        st.markdown("---")
        st.markdown("**Gene → Drug Map**")
        for drug, gene in GENE_DRUG_MAP.items():
            st.markdown(f"`{gene}` → {drug.title()}")

    key       = groq_key.strip() if groq_key else ""
    skip_llm  = use_static or not key

    render_nav(bool(key))

    tab_analysis, tab_suite = st.tabs(["Analysis", "Test Suite"])

    # ── Analysis Tab ──────────────────────────────────────────────────────────
    with tab_analysis:
        has_results = "results" in st.session_state and st.session_state["results"]
        render_steps(has_vcf=True, has_drugs=True, has_results=has_results)

        col_input, col_results = st.columns([1, 2], gap="large")

        with col_input:
            sec("Genomic Data")

            # ── File uploader FIRST — prevents selectbox dropdown from overlapping it ──
            vcf_file = st.file_uploader(
                "Upload VCF file (.vcf)",
                type=["vcf"],
                help="Upload a VCF v4.2 pharmacogenomics file — limit 200 MB",
            )

            st.markdown('<div style="height:4px;"></div>', unsafe_allow_html=True)

            # ── Scenario selectbox BELOW the uploader ──
            persona_sel = st.selectbox(
                "Or load a test scenario",
                options=["None"] + [p["label"] for p in PERSONAS.values()],
                key="persona_sel",
            )

            # Resolve VCF text
            vcf_text = None
            if vcf_file:
                vcf_text = vcf_file.read().decode("utf-8")
            elif persona_sel != "None":
                for p in PERSONAS.values():
                    if p["label"] == persona_sel:
                        try:
                            vcf_text = load_vcf(p["file"])
                        except FileNotFoundError:
                            vcf_text = get_sample_vcf()
                        break

            if vcf_text:
                fname = getattr(vcf_file, 'name', persona_sel) if vcf_file else persona_sel
                fsize_kb = len(vcf_text) / 1024
                st.success(f"✓ {fname} · {fsize_kb:.2f} KB — Ready for analysis")

            st.markdown('<div style="height:4px;"></div>', unsafe_allow_html=True)
            sec("Medications to Analyse")
            default_drugs = ALL_DRUGS
            if persona_sel != "None":
                for p in PERSONAS.values():
                    if p["label"] == persona_sel:
                        default_drugs = p["drugs"]
                        break
            selected_drugs = st.multiselect("Select drugs", ALL_DRUGS,
                default=default_drugs, label_visibility="collapsed")
            custom_raw = st.text_input("Custom drugs (comma-separated)", placeholder="CODEINE, WARFARIN…")
            if custom_raw:
                extras = [d.strip().upper() for d in custom_raw.split(",") if d.strip()]
                selected_drugs = list(dict.fromkeys(selected_drugs + extras))

            sec("Patient ID")
            patient_id_input = st.text_input("Patient ID", placeholder="Auto-generated if blank",
                                              label_visibility="collapsed")
            pid = patient_id_input.strip() or f"PG-{uuid.uuid4().hex[:8].upper()}"

            sec("Quick Demo Personas")
            persona_cols = st.columns(2)
            SEV_COLORS_LOCAL = {
                "critical": ("#FEF2F2","#FECACA","#7F1D1D"),
                "high":     ("#FFF7ED","#FED7AA","#7C2D12"),
                "moderate": ("#FFFBEB","#FDE68A","#78350F"),
                "none":     ("#F0FDF4","#BBF7D0","#14532D"),
            }
            for pi, (persona_id, p) in enumerate(PERSONAS.items()):
                bg, border, txt = SEV_COLORS_LOCAL.get(p["sev"], SEV_COLORS_LOCAL["none"])
                with persona_cols[pi % 2]:
                    st.markdown(
                        f'''<div style="background:{bg};border:1.5px solid {border};border-radius:10px;
                        padding:10px 12px;margin-bottom:8px;">
                        <div style="font-size:.8rem;font-weight:700;color:{txt};">{p["label"]}</div>
                        <div style="font-family:monospace;font-size:.65rem;color:{txt};opacity:.75;">{p["desc"]}</div>
                        </div>''', unsafe_allow_html=True)
                    if st.button(f"Load", key=f"persona2_{persona_id}", use_container_width=True):
                        try:
                            vcf_text = load_vcf(p["file"])
                        except FileNotFoundError:
                            vcf_text = get_sample_vcf()
                        pid_gen = f"PG-{uuid.uuid4().hex[:8].upper()}"
                        with st.spinner(f"Running {p['label']}…"):
                            parsed, results, outputs, ix, pdf = run_pipeline(
                                vcf_text, p["drugs"], pid_gen, key, skip_llm=not bool(key))
                        st.session_state["results"]      = outputs
                        st.session_state["parsed"]       = parsed
                        st.session_state["ix"]           = ix
                        st.session_state["pdf"]          = pdf
                        st.session_state["patient_id"]   = pid_gen
                        st.session_state["results_key"]  = key
                        st.session_state["results_skip"] = not bool(key)
                        st.rerun()

            sec("View Mode")
            patient_mode = st.checkbox("Patient-friendly view (plain language)", value=False)

            run_btn = st.button("Run Analysis →", use_container_width=True,
                                disabled=not vcf_text or not selected_drugs)

            if run_btn and vcf_text and selected_drugs:
                with st.spinner("Analysing pharmacogenomic profile…"):
                    parsed, results, outputs, ix, pdf = run_pipeline(
                        vcf_text, selected_drugs, pid, key,
                        run_ix=len(selected_drugs) > 1,
                        gen_pdf=True,
                        skip_llm=skip_llm)
                st.session_state["results"]      = outputs
                st.session_state["parsed"]       = parsed
                st.session_state["ix"]           = ix
                st.session_state["pdf"]          = pdf
                st.session_state["patient_id"]   = pid
                st.session_state["results_key"]  = key
                st.session_state["results_skip"] = skip_llm
                st.rerun()

        with col_results:
            sec("Analysis Results")
            if "results" in st.session_state and st.session_state["results"]:
                res_pid  = st.session_state["patient_id"]
                res_outs = st.session_state["results"]
                res_par  = st.session_state["parsed"]
                res_ix   = st.session_state.get("ix")
                res_pdf  = st.session_state.get("pdf")
                res_key  = st.session_state.get("results_key", key)
                res_skip = st.session_state.get("results_skip", skip_llm)
                st.markdown(f"""
                <div style="display:flex;align-items:center;gap:var(--sp-3);margin-bottom:var(--sp-4);">
                  <span style="font-family:var(--font-mono);font-size:1rem;font-weight:700;
                    color:#1D4ED8;background:#EFF6FF;border:1px solid #BFDBFE;
                    padding:4px 12px;border-radius:9999px;">{res_pid}</span>
                </div>""", unsafe_allow_html=True)
                render_results(res_outs, res_par, res_ix, res_pdf, res_pid,
                               patient_mode=patient_mode, key=res_key, skip_llm=res_skip)
            else:
                st.markdown("""
                <div class="empty-state">
                  <span class="empty-icon">🧬</span>
                  <div class="empty-title">No analysis results yet</div>
                  <div class="empty-hint">
                    Upload a VCF file or select a scenario<br>
                    Choose medications to analyse<br>
                    Click Run Analysis →
                  </div>
                </div>""", unsafe_allow_html=True)

    # ── Test Suite Tab ────────────────────────────────────────────────────────
    with tab_suite:
        render_test_suite(key)


if __name__ == "__main__":
    main()