# 🧬 SurakshaRx — Pharmacogenomic Risk Prediction System

> RIFT 2026 Hackathon | Pharmacogenomics / Explainable AI Track | **v10.0**

[![Live Demo](https://img.shields.io/badge/Live-Demo-brightgreen)](https://pharmaguard-mkmkwjxgblr9fx4mwx9jdc.streamlit.app/)
[![Demo Video](https://img.shields.io/badge/Demo-Video-blue)](https://drive.google.com/file/d/1JGpvu2KgMW0KmkMa3ot3NsQAi3XPAP3L/view?usp=sharing)
[![GitHub](https://img.shields.io/badge/GitHub-SurakshaRx-black)](https://github.com/Farbricated/SurakshaRx)
[![CPIC](https://img.shields.io/badge/CPIC-Aligned-orange)](https://cpicpgx.org)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://python.org)
[![Groq](https://img.shields.io/badge/LLM-LLaMA_3.3_70B-purple)](https://groq.com)

## 🔗 Links
- **Live Demo:** https://pharmaguard-mkmkwjxgblr9fx4mwx9jdc.streamlit.app/
- **Demo Video:** https://drive.google.com/file/d/1JGpvu2KgMW0KmkMa3ot3NsQAi3XPAP3L/view?usp=sharing
- **GitHub:** https://github.com/Farbricated/SurakshaRx

---

## 📋 What is SurakshaRx?

Adverse drug reactions kill **over 100,000 Americans annually** — and many are preventable. Most people don't know their genes affect how medicines work in their body. SurakshaRx changes that.

Upload a genetic data file (VCF), select medications, and SurakshaRx instantly tells you:
- Which drugs are **safe** for you
- Which need a **different dose**
- Which are **toxic** or **won't work** based on your DNA
- **Why** — with full AI-generated clinical explanations

All powered by CPIC guidelines (the gold standard in clinical pharmacogenomics) and LLaMA 3.3 70B via Groq.

---

## ✨ Current Features (v10.0)

### Core Analysis Engine
| Feature | Description |
|---------|-------------|
| **VCF v4.2 Parsing** | Genotype-aware parser — reads GT field, skips 0/0 (homozygous reference) variants; handles heterozygous and homozygous alt correctly |
| **Diplotype → Phenotype Mapping** | Full CPIC star-allele tables for all 6 genes; reverse diplotype lookup; xN duplication handling for URM |
| **Risk Assessment** | Per-drug risk label (Safe / Adjust Dosage / Toxic / Ineffective), severity, and confidence score |
| **Drug Interaction Checker** | Shared-gene risk detection, known dangerous combos, CYP inhibitor phenocopying |
| **PDF Report Generation** | Clinical-grade PDF via fpdf2; Unicode-safe (ASCII fallback for Helvetica compatibility) |

### AI & LLM
| Feature | Description |
|---------|-------------|
| **LLaMA 3.3 70B via Groq** | Per-drug clinical explanations: Summary, Biological Mechanism, Variant Significance, Clinical Implications |
| **Thread-safe RLock cache** | Same drug/phenotype pair served instantly; zero duplicate API calls in parallel mode |
| **Rate-limit handling** | 3 retries with 1s/3s/8s backoff; sets `_RATE_LIMITED` flag to skip downstream calls |
| **Static expert templates** | 18 hand-written fallback templates covering all 6 drugs × key phenotypes; all cite rsIDs |
| **Unified Patient Narrative** | Holistic paragraph summary across all drugs; has its own cache key; skips API if rate-limited |

### Visualisations
| Feature | Description |
|---------|-------------|
| **Polygenic Risk Score** | Composite 0–100 score weighted by drug criticality (FLUOROURACIL 1.4×, AZATHIOPRINE 1.3×, etc.) |
| **Drug × Gene Heatmap** | Color-coded matrix of all selected drugs × all genes with hover details |
| **Chromosome Visualisation** | Animated variant markers positioned on chromosome diagrams by genomic band |
| **Population Frequency** | Bar chart comparing patient phenotype to global population frequencies per gene |
| **Patient Plain-English Mode** | Toggle converts all clinical output to plain-language patient cards |
| **Before & After PGx Panel** | Side-by-side clinical impact: without vs. with pharmacogenomic guidance |
| **Prescription Safety Checker** | Per-drug safety verdict widget with gene/phenotype context |
| **One-Click Clinical Note** | Downloadable plain-text clinical note for EHR pasting |

### UI & UX
| Feature | Description |
|---------|-------------|
| **4 Quick Demo Personas** | Pre-loaded patient scenarios (Critical Risk, Warfarin PM, Drug Interaction, All Safe) with graceful fallback if sample files are missing |
| **Test Suite tab** | 4 validated scenarios run with static templates (no API key required); results persist in session state |
| **Parallel test execution** | `concurrent.futures.ThreadPoolExecutor` — all 4 scenarios in ~3s |
| **File uploader UX fix** | Drop zone hidden; clean "Browse" button only; file size displays correctly in KB |
| **Input layout fix** | File uploader rendered before selectbox to prevent dropdown overlap |
| **Custom CSS design system** | DM Sans + JetBrains Mono; CSS variables; animated cards; responsive grid |

---

## 🏗️ Architecture

```
User uploads VCF file + selects drug(s)
        ↓
VCF Parser (vcf_parser.py)
  • Parses VCF v4.2 format
  • Reads GT field — skips 0/0 (homozygous ref) variants
  • Handles heterozygous (0/1) and homozygous alt (1/1) correctly
  • Extracts diplotype data per gene via star allele accumulation
  • Falls back to rsID-based gene inference if GENE= INFO tag is absent
        ↓
Risk Engine (risk_engine.py)
  • Diplotype → Phenotype lookup (CPIC tables, with reverse lookup)
  • Phenotype × Drug → Risk Label + Severity + Confidence
  • Returns: risk_label, severity, confidence_score, cpic_recommendation
        ↓
┌─────────────────────────────────────────────────────────────┐
│  LLM Explainer (llm_explainer.py) — Thread-safe v5.1        │
│  • Thread-safe RLock cache: same drug/phenotype served      │
│    instantly; dedicated cache key for patient narrative     │
│  • LLaMA 3.3 70B via Groq, 3 retries + 1s/3s/8s backoff   │
│  • _RATE_LIMITED flag: if hit during drug calls, narrative  │
│    skips API immediately (no spinner hang)                  │
│  • 18 static expert templates as fallback (rsID-cited)      │
│  • skip_llm=True path for instant test mode                 │
└─────────────────────────────────────────────────────────────┘
        ↓
Schema Builder (schema.py)
  • Assembles fully validated JSON output
  • All fields populated — no nulls
  • SLCO1B1 phenotype keys use standard enum (NM/IM/PM)
  • Variant output: rsid, gene, star_allele, functional_status only
        ↓
Drug Interactions (drug_interactions.py)
  • Shared-gene risk detection by phenotype
  • 6 known dangerous combination alerts
  • CYP inhibitor phenocopying detection
        ↓
PDF Report (pdf_report.py)
  • fpdf2 clinical report with Unicode → ASCII sanitisation
  • Per-drug: risk, variants table, CPIC recommendation, AI explanation
        ↓
Streamlit UI (app.py v10.0) — All visual features + UX fixes
  • Polygenic Risk Score gauge
  • Drug × Gene Heatmap
  • Chromosome Visualisation
  • Population Frequency bars
  • Patient Plain-English Mode
  • Before & After PGx panel
  • Prescription Safety Checker
  • Clinical Note generator
  • PDF + CSV + JSON download
  • Parallel Test Suite
```

---

## 🧬 Supported Genes & Drugs

| Gene | Drug | Clinical Significance |
|------|------|-----------------------|
| **CYP2D6** | Codeine | URM = respiratory death risk; PM = zero pain relief |
| **CYP2C19** | Clopidogrel | PM = heart attack risk — drug won't work at all |
| **CYP2C9** | Warfarin | PM = fatal bleeding risk at standard doses |
| **SLCO1B1** | Simvastatin | PM = rhabdomyolysis (muscle destruction) |
| **TPMT** | Azathioprine | PM = bone marrow destruction at standard doses |
| **DPYD** | Fluorouracil | PM = fatal chemotherapy toxicity |

---

## 📊 Test Suite — All 4 Scenarios

| Scenario | Drugs | Key Variants | Expected Results |
|----------|-------|-------------|-----------------|
| **Mixed Variants** | Clopidogrel, Codeine, Azathioprine | CYP2C19 *2/*3, CYP2D6 *4/*4, TPMT *3B/*3C | Ineffective, Ineffective, Toxic ✅ |
| **UltraRapid Metabolizer** | Codeine, Clopidogrel | CYP2D6 *1xN/*1xN | Toxic, Safe ✅ |
| **All Normal Wild-type** | All 6 drugs | *1/*1 across all genes | All Safe ✅ |
| **Worst Case — All PM** | All 6 drugs | LOF alleles all 6 genes | Ineffective×2, Adjust×1, Toxic×3 ✅ |

All tests run with static templates (no API key required). Results persist in session state and load into the Analysis tab automatically.

---

## 🛠️ Tech Stack

| Layer | Technology |
|-------|-----------|
| Frontend + Backend | Streamlit (Python) |
| LLM | Groq — LLaMA 3.3 70B Versatile |
| LLM Fallback | Static expert templates (rsID-cited, 18 templates) |
| Parallel Execution | `concurrent.futures.ThreadPoolExecutor` |
| Thread Safety | `threading.RLock` for LLM cache + `threading.Event` for rate-limit flag |
| PDF Reports | fpdf2 (with Unicode → ASCII sanitiser) |
| Data Validation | Pydantic v2 |
| VCF Parsing | Pure Python — no external bioinformatics dependencies |
| Guidelines | CPIC (cpicpgx.org) — Level A evidence |
| Deployment | Streamlit Community Cloud |

---

## 🚀 Installation & Local Setup

```bash
# 1. Clone the repository
git clone https://github.com/Farbricated/SurakshaRx.git
cd SurakshaRx

# 2. Create virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate        # macOS/Linux
# venv\Scripts\activate         # Windows

# 3. Install dependencies
pip install -r requirements.txt

# 4. Set up environment variables
cp .env.example .env
# Edit .env and add: GROQ_API_KEY=your_key_here

# 5. Run the app
streamlit run app.py
```

The app will open at `http://localhost:8501`

---

## 🔑 Environment Variables

```env
GROQ_API_KEY=your_groq_api_key_here
```

Get your **free** Groq API key at: https://console.groq.com

> **No API key?** The app works fully without one — all features run using expert-written static templates that cite all relevant rsIDs and star alleles. Enable "Test mode: instant (no API call)" in the sidebar.

---

## 📤 Output JSON Schema

Every analysis produces a fully-validated JSON file:

```json
{
  "patient_id": "PG-A1B2C3D4",
  "drug": "CLOPIDOGREL",
  "timestamp": "2026-02-19T14:30:00Z",
  "risk_assessment": {
    "risk_label": "Ineffective",
    "confidence_score": 0.94,
    "severity": "high"
  },
  "pharmacogenomic_profile": {
    "primary_gene": "CYP2C19",
    "diplotype": "*2/*3",
    "phenotype": "PM",
    "detected_variants": [
      {
        "rsid": "rs4244285",
        "gene": "CYP2C19",
        "star_allele": "*2",
        "functional_status": "no_function"
      }
    ],
    "cpic_evidence_level": "Level A"
  },
  "clinical_recommendation": {
    "cpic_guideline": "CPIC Guideline for CLOPIDOGREL",
    "dosing_recommendation": "Use alternative antiplatelet: prasugrel or ticagrelor.",
    "alternative_drugs": ["Prasugrel", "Ticagrelor"],
    "monitoring_required": "Platelet function testing. Monitor for ischaemic events."
  },
  "llm_generated_explanation": {
    "summary": "...",
    "biological_mechanism": "...",
    "variant_significance": "...",
    "clinical_implications": "...",
    "model_used": "llama-3.3-70b-versatile"
  },
  "quality_metrics": {
    "vcf_parsing_success": true,
    "variants_detected": 2,
    "genes_analyzed": ["CYP2C19", "CYP2D6", "CYP2C9", "SLCO1B1", "TPMT", "DPYD"],
    "parse_errors": [],
    "explanation_generated": true
  }
}
```

---

## 📁 Project Structure

```
SurakshaRx/
├── app.py                          # Main Streamlit app (v9.3)
├── vcf_parser.py                   # VCF v4.2 parsing — genotype-aware (GT field)
├── risk_engine.py                  # CPIC-based risk prediction engine
├── llm_explainer.py                # Groq LLM explainer (thread-safe, rate-limit aware)
├── schema.py                       # Output JSON schema builder
├── drug_interactions.py            # Drug-drug interaction checker
├── pdf_report.py                   # Clinical PDF report generator (Unicode-safe)
├── requirements.txt                # Python dependencies
├── .env.example                    # Environment variable template
├── .gitignore                      # Git ignore (excludes .env, patient data)
├── patient_high_risk_multigene.vcf # Example multi-gene high-risk VCF
├── sample_data/
│   ├── sample.vcf                          # Mixed variants test case
│   ├── patient_a_critical.vcf              # Critical risk: CYP2D6 PM + DPYD PM + TPMT PM
│   ├── patient_b_warfarin.vcf              # Warfarin PM: CYP2C9 *2/*3
│   ├── patient_c_interaction.vcf           # Drug interaction: CYP2C19 *2/*3
│   ├── patient_d_safe.vcf                  # All safe: wildtype *1/*1
│   ├── test_ultrarapid_metabolizer.vcf     # CYP2D6 *1xN/*1xN URM
│   ├── test_all_normal_wildtype.vcf        # All NM baseline
│   └── test_worst_case_all_pm.vcf          # All PM: LOF alleles across 6 genes
└── README.md
```

---

## 📖 Usage Examples

### Example 1 — Run the sample analysis
1. Open the app at https://pharmaguard-mkmkwjxgblr9fx4mwx9jdc.streamlit.app/
2. Select **"Mixed Variants"** from the scenario dropdown or use a Quick Demo Persona
3. Select drugs: **Clopidogrel + Codeine + Azathioprine**
4. Click **Run Analysis →**
5. See: Polygenic Risk Score, Drug×Gene Heatmap, Chromosome Visualisation, AI explanations

### Example 2 — Patient Plain-English Mode
1. Run any analysis
2. Toggle **"Patient-friendly view (plain language)"** in the left panel
3. All clinical jargon converts to plain cards a patient can read directly

### Example 3 — Run all tests (no API key needed)
1. Go to the **Test Suite** tab
2. Click any scenario button — runs with static templates instantly
3. Results appear in the Test Suite tab and load into the Analysis tab automatically
4. Download full JSON from the Analysis tab

### Example 4 — API integration
```python
import json

result = json.load(open("SurakshaRx_PATIENT_001_WARFARIN.json"))
risk      = result["risk_assessment"]["risk_label"]            # "Adjust Dosage"
gene      = result["pharmacogenomic_profile"]["primary_gene"]  # "CYP2C9"
diplotype = result["pharmacogenomic_profile"]["diplotype"]     # "*2/*3"
rec       = result["clinical_recommendation"]["dosing_recommendation"]
alts      = result["clinical_recommendation"]["alternative_drugs"]
```

---

## 🆕 Changelog

### v10.0 (Current)
- 🎨 **UI/UX REDESIGN:** Complete overhaul to a modern dark glassmorphism theme (blurs, gradients, glowing borders).
- 🐛 **BUG FIX:** Fixed file handle leak in `load_vcf()` by implementing context managers.
- 🐛 **BUG FIX:** Fixed unclosed HTML divs in drug card rendering that caused DOM structure issues.
- 🐛 **BUG FIX:** Fixed `datetime.utcnow()` deprecation warnings in `app.py`, `schema.py`, and `pdf_report.py`.
- 📝 **BRAND FIX:** Corrected navigation header to display `SurakshaRx` instead of `PharmaGuard`.

### v9.3
- 🐛 **UI FIX:** File uploader label visible (was "collapsed"), preventing overlap with scenario selectbox
- 🐛 **UI FIX:** Selectbox label explicitly set to "visible" for clear visual hierarchy
- 🐛 **UI FIX:** File size display corrected (was showing KB as MB)
- 🐛 **UI FIX:** File uploader rendered before scenario selectbox in DOM order
- 🐛 **TEST SUITE FIX:** `load_vcf()` calls wrapped in try/except — falls back to `get_sample_vcf()` when `sample_data/` files are missing from deployment
- 🐛 **PERSONA DEMO FIX:** Same try/except fallback for persona file loading

### v9.2
- Drug × Gene Heatmap column/row rendering fixes
- Persona demo cards display in 4-column layout
- Test suite results persist across reruns via session state

### v5.1 (llm_explainer.py)
- `generate_patient_narrative()` checks dedicated narrative cache key before any API call
- Narrative API call uses `max_tokens=300` to reduce rate-limit pressure
- `_RATE_LIMITED` threading.Event flag: if hit during per-drug calls, narrative skips API entirely
- Eliminates infinite spinner after 6 prior Groq calls

### v5.0
- Polygenic Risk Score (0–100 composite gauge)
- Drug × Gene Heatmap
- Chromosome Visualisation with animated markers
- Population Frequency bars
- Patient Plain-English Mode toggle
- Parallel test execution via ThreadPoolExecutor
- Thread-safe RLock LLM cache

### v4.2
- Reduced backoff: 5s/15s/30s → 1s/3s/8s
- All static templates cite rsIDs

### v4.1
- Fixed TEST_SUITE expected values (CODEINE *4/*4 → Ineffective, AZATHIOPRINE *3B/*3C → Toxic)

### v4.0
- LLM retry + fallback to static templates
- Drug-drug interaction checker
- PDF report generation
- Drug-specific monitoring recommendations

### VCF Parser (v4.1)
- **GT field parsing:** Skips 0/0 (homozygous reference) variants — patients who do NOT carry a variant are no longer incorrectly flagged
- Handles 0/1 (heterozygous) and 1/1 (homozygous alt) zygosity correctly for diplotype construction
- rsID-based gene inference fallback for VCFs without GENE= INFO tags

### Risk Engine / Schema fixes
- SLCO1B1 diplotype table updated to standard phenotype enum: NM / IM / PM
- SIMVASTATIN risk table, CPIC recommendations, and schema keys aligned to NM / IM / PM
- `contraindicated` field removed from schema output (was causing Pydantic validation errors)
- Variant output in JSON limited to: rsid, gene, star_allele, functional_status (no chrom/pos/ref/alt)

---

## 👥 Team

| Name | Role |
|------|------|
| Member 1 | Lead Full-Stack Developer — Streamlit UI, pipeline architecture, VCF parser, deployment |
| Member 2 | Clinical & Pharmacogenomics Domain Expert — CPIC guidelines, risk engine, static templates, drug interaction logic |
| Member 3 | AI/ML Engineer — LLM integration, Groq API, prompt engineering, thread-safe caching, explainability layer |


---

## 🏆 Hackathon Submission

**Event:** RIFT 2026  
**Track:** Pharmacogenomics / Explainable AI  
**Hashtags:** `#RIFT2026` `#SurakshaRx` `#Pharmacogenomics` `#AIinHealthcare` `#ExplainableAI`

---

## ⚠️ Disclaimer

SurakshaRx is a hackathon prototype for educational purposes only. It is **NOT** a medical device and should **NOT** be used for clinical decision-making without validation by a qualified clinical pharmacologist or geneticist. All recommendations should be verified against current CPIC guidelines at [cpicpgx.org](https://cpicpgx.org).