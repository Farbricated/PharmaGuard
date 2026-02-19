"""
PDF Report Generator for PharmaGuard v2.0
Generates professional clinical pharmacogenomic reports
Uses fpdf2 library
"""

from fpdf import FPDF
from datetime import datetime
from typing import List, Dict
import io


RISK_COLORS = {
    "Safe":          (6,   95,  70),   # dark green
    "Adjust Dosage": (120, 53,  15),   # dark amber
    "Toxic":         (127, 29,  29),   # dark red
    "Ineffective":   (49,  46,  129),  # dark purple
    "Unknown":       (31,  41,  55),   # dark gray
}

SEVERITY_COLORS = {
    "none":     (16,  185, 129),
    "low":      (251, 191, 36),
    "moderate": (249, 115, 22),
    "high":     (239, 68,  68),
    "critical": (185, 28,  28),
}


class PharmaGuardPDF(FPDF):
    def __init__(self):
        super().__init__()
        self.set_margins(15, 15, 15)
        self.set_auto_page_break(auto=True, margin=20)

    def header(self):
        # Dark header bar
        self.set_fill_color(15, 23, 42)
        self.rect(0, 0, 210, 22, 'F')
        self.set_text_color(224, 242, 254)
        self.set_font("Helvetica", "B", 13)
        self.set_xy(10, 5)
        self.cell(0, 12, "PharmaGuard | Pharmacogenomic Risk Report", ln=False)
        self.set_font("Helvetica", "", 8)
        self.set_text_color(148, 163, 184)
        self.set_xy(10, 14)
        self.cell(0, 6, f"Generated: {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')} | RIFT 2026 | Powered by Groq LLaMA 3.3", ln=True)
        self.ln(8)

    def footer(self):
        self.set_y(-15)
        self.set_fill_color(15, 23, 42)
        self.rect(0, 282, 210, 15, 'F')
        self.set_text_color(100, 116, 139)
        self.set_font("Helvetica", "I", 7)
        self.set_xy(10, 284)
        self.cell(0, 6, "DISCLAIMER: PharmaGuard is a research tool. Not for clinical use without validation. Consult CPIC guidelines at cpicpgx.org", ln=False)
        self.set_xy(-30, 284)
        self.cell(0, 6, f"Page {self.page_no()}", align="R")

    def section_title(self, title: str, color=(30, 64, 175)):
        self.set_fill_color(*color)
        self.set_text_color(255, 255, 255)
        self.set_font("Helvetica", "B", 10)
        self.cell(0, 8, f"  {title}", ln=True, fill=True)
        self.ln(2)
        self.set_text_color(0, 0, 0)

    def key_value(self, key: str, value: str, bold_val=False):
        self.set_font("Helvetica", "B", 9)
        self.set_text_color(71, 85, 105)
        self.cell(55, 6, key + ":", ln=False)
        self.set_font("Helvetica", "B" if bold_val else "", 9)
        self.set_text_color(15, 23, 42)
        self.cell(0, 6, str(value), ln=True)

    def colored_badge(self, label: str, risk_label: str):
        color = RISK_COLORS.get(risk_label, (31, 41, 55))
        self.set_fill_color(*color)
        self.set_text_color(255, 255, 255)
        self.set_font("Helvetica", "B", 10)
        self.cell(50, 9, f"  {label}", ln=False, fill=True)
        self.ln(11)
        self.set_text_color(0, 0, 0)

    def info_box(self, text: str, bg=(241, 245, 249), border_color=(59, 130, 246)):
        x = self.get_x()
        y = self.get_y()
        # Left border accent
        self.set_fill_color(*border_color)
        self.rect(x, y, 3, 20, 'F')
        # Background
        self.set_fill_color(*bg)
        self.rect(x + 3, y, 177, 20, 'F')
        self.set_xy(x + 6, y + 2)
        self.set_font("Helvetica", "", 8)
        self.set_text_color(30, 41, 59)
        self.multi_cell(170, 4.5, text)
        self.ln(2)
        self.set_text_color(0, 0, 0)

    def alert_box(self, text: str, severity: str = "high"):
        colors = {
            "critical": ((69, 10, 10),   (239, 68, 68)),
            "high":     ((69, 10, 10),   (239, 68, 68)),
            "moderate": ((69, 26, 3),    (249, 115, 22)),
            "low":      ((120, 53, 15),  (251, 191, 36)),
        }
        bg, border = colors.get(severity, ((31, 41, 55), (107, 114, 128)))
        x, y = self.get_x(), self.get_y()
        self.set_fill_color(*bg)
        self.rect(x, y, 180, 14, 'F')
        self.set_fill_color(*border)
        self.rect(x, y, 4, 14, 'F')
        self.set_xy(x + 7, y + 3)
        self.set_font("Helvetica", "B", 8)
        self.set_text_color(254, 226, 226)
        self.multi_cell(170, 4.5, f"CLINICAL ALERT: {text}")
        self.ln(4)
        self.set_text_color(0, 0, 0)


def generate_pdf_report(patient_id: str, all_outputs: List[Dict], parsed_vcf: Dict) -> bytes:
    """Generate a complete PDF clinical report and return as bytes."""
    
    pdf = PharmaGuardPDF()
    pdf.add_page()

    # ── Patient Summary Header ─────────────────────────────────
    pdf.set_fill_color(248, 250, 252)
    pdf.rect(15, 30, 180, 28, 'F')
    pdf.set_xy(18, 32)
    pdf.set_font("Helvetica", "B", 16)
    pdf.set_text_color(15, 23, 42)
    pdf.cell(0, 10, "Pharmacogenomic Risk Assessment", ln=True)
    pdf.set_x(18)
    pdf.set_font("Helvetica", "", 10)
    pdf.set_text_color(71, 85, 105)
    pdf.cell(60, 7, f"Patient ID: {patient_id}", ln=False)
    pdf.cell(60, 7, f"Genes Analyzed: {len(parsed_vcf.get('detected_genes', []))}/6", ln=False)
    pdf.cell(0, 7, f"Drugs Evaluated: {len(all_outputs)}", ln=True)
    pdf.set_x(18)
    pdf.cell(0, 7, f"Variants Detected: {parsed_vcf.get('total_variants', 0)}  |  Report Date: {datetime.utcnow().strftime('%B %d, %Y')}", ln=True)
    pdf.ln(6)

    # ── Genomic Profile Summary ────────────────────────────────
    pdf.section_title("GENOMIC PROFILE SUMMARY", color=(15, 23, 42))
    
    genes_found = parsed_vcf.get("detected_genes", [])
    all_genes = ["CYP2D6", "CYP2C19", "CYP2C9", "SLCO1B1", "TPMT", "DPYD"]
    
    pdf.set_font("Helvetica", "B", 8)
    pdf.set_text_color(100, 116, 139)
    pdf.cell(0, 5, "GENE STATUS OVERVIEW", ln=True)
    pdf.ln(1)
    
    for gene in all_genes:
        found = gene in genes_found
        pdf.set_fill_color(6, 95, 70) if found else pdf.set_fill_color(51, 65, 85)
        pdf.set_text_color(255, 255, 255)
        pdf.set_font("Helvetica", "B", 8)
        pdf.cell(32, 7, f"  {gene}", fill=True, ln=False)
        pdf.set_fill_color(248, 250, 252)
        pdf.set_text_color(71, 85, 105)
        pdf.set_font("Helvetica", "", 8)
        status = "Variants Detected" if found else "Wild-type (*1/*1)"
        pdf.cell(58, 7, f"  {status}", fill=True, ln=False)
        pdf.ln(0)
        if all_genes.index(gene) % 2 == 1:
            pdf.ln(8)
        else:
            pdf.set_x(15 + 90)

    pdf.ln(10)

    # ── Per-Drug Risk Cards ─────────────────────────────────────
    for i, output in enumerate(all_outputs):
        if pdf.get_y() > 220:
            pdf.add_page()

        drug = output["drug"]
        risk_label = output["risk_assessment"]["risk_label"]
        severity = output["risk_assessment"]["severity"]
        confidence = output["risk_assessment"]["confidence_score"]
        gene = output["pharmacogenomic_profile"]["primary_gene"]
        diplotype = output["pharmacogenomic_profile"]["diplotype"]
        phenotype = output["pharmacogenomic_profile"]["phenotype"]
        exp = output["llm_generated_explanation"]

        # Drug section header
        color = RISK_COLORS.get(risk_label, (31, 41, 55))
        pdf.set_fill_color(*color)
        pdf.set_text_color(255, 255, 255)
        pdf.set_font("Helvetica", "B", 11)
        pdf.cell(0, 10, f"  {drug.title()}  -  {risk_label}  (Severity: {severity.upper()})", ln=True, fill=True)
        pdf.ln(2)

        # Alert box for high/critical
        if severity in ("high", "critical", "moderate"):
            rec = output["clinical_recommendation"]["dosing_recommendation"]
            pdf.alert_box(rec, severity)

        # Key metrics row
        pdf.set_text_color(15, 23, 42)
        pdf.key_value("Primary Gene", gene, bold_val=True)
        pdf.key_value("Diplotype", diplotype, bold_val=True)
        pdf.key_value("Phenotype", phenotype, bold_val=True)
        pdf.key_value("Confidence Score", f"{confidence:.0%}")
        pdf.ln(2)

        # Detected variants table
        variants = output["pharmacogenomic_profile"]["detected_variants"]
        if variants:
            pdf.set_font("Helvetica", "B", 8)
            pdf.set_fill_color(30, 41, 59)
            pdf.set_text_color(255, 255, 255)
            pdf.cell(35, 6, "  rsID", fill=True)
            pdf.cell(25, 6, "Star Allele", fill=True)
            pdf.cell(20, 6, "REF>ALT", fill=True)
            pdf.cell(0, 6, "Functional Status", fill=True, ln=True)
            
            for j, v in enumerate(variants[:6]):
                pdf.set_fill_color(248, 250, 252) if j % 2 == 0 else pdf.set_fill_color(241, 245, 249)
                pdf.set_text_color(15, 23, 42)
                pdf.set_font("Helvetica", "", 8)
                pdf.cell(35, 5.5, f"  {v.get('rsid','N/A')}", fill=True)
                pdf.cell(25, 5.5, v.get("star_allele","N/A"), fill=True)
                pdf.cell(20, 5.5, f"{v.get('ref','?')}>{v.get('alt','?')}", fill=True)
                pdf.cell(0, 5.5, v.get("functional_status","Unknown"), fill=True, ln=True)
            pdf.ln(3)

        # LLM Explanation sections
        if exp.get("summary"):
            pdf.set_font("Helvetica", "B", 8)
            pdf.set_text_color(59, 130, 246)
            pdf.cell(0, 5, "AI CLINICAL EXPLANATION (Groq LLaMA 3.3)", ln=True)
            
            sections = [
                ("Summary", exp.get("summary", "")),
                ("Biological Mechanism", exp.get("biological_mechanism", "")),
                ("Variant Significance", exp.get("variant_significance", "")),
                ("Clinical Implications", exp.get("clinical_implications", "")),
            ]
            for label, text in sections:
                if text:
                    pdf.set_font("Helvetica", "B", 8)
                    pdf.set_text_color(30, 64, 175)
                    pdf.cell(0, 5, label, ln=True)
                    pdf.set_font("Helvetica", "", 8)
                    pdf.set_text_color(51, 65, 85)
                    pdf.multi_cell(0, 4.5, text)
                    pdf.ln(1)

        # CPIC Recommendation
        pdf.set_font("Helvetica", "B", 8)
        pdf.set_text_color(6, 95, 70)
        pdf.cell(0, 5, "CPIC DOSING RECOMMENDATION", ln=True)
        pdf.set_font("Helvetica", "", 8)
        pdf.set_text_color(15, 23, 42)
        pdf.multi_cell(0, 4.5, output["clinical_recommendation"]["dosing_recommendation"])
        
        alts = output["clinical_recommendation"].get("alternative_drugs", [])
        if alts:
            pdf.set_font("Helvetica", "B", 8)
            pdf.set_text_color(71, 85, 105)
            pdf.cell(30, 5, "Alternatives:", ln=False)
            pdf.set_font("Helvetica", "", 8)
            pdf.set_text_color(15, 23, 42)
            pdf.cell(0, 5, ", ".join(alts), ln=True)

        pdf.ln(6)
        # Divider
        pdf.set_draw_color(203, 213, 225)
        pdf.line(15, pdf.get_y(), 195, pdf.get_y())
        pdf.ln(6)

    # ── Disclaimer page ──────────────────────────────────────────
    pdf.add_page()
    pdf.section_title("IMPORTANT CLINICAL DISCLAIMER", color=(127, 29, 29))
    pdf.set_font("Helvetica", "", 9)
    pdf.set_text_color(30, 41, 59)
    disclaimer = (
        "This report has been generated by PharmaGuard, an AI-powered pharmacogenomics research tool "
        "developed for the RIFT 2026 Hackathon. The risk assessments and recommendations provided herein "
        "are based on established CPIC (Clinical Pharmacogenomics Implementation Consortium) guidelines "
        "and LLM-generated explanations.\n\n"
        "THIS REPORT IS NOT A MEDICAL DEVICE AND SHOULD NOT BE USED FOR CLINICAL DECISION-MAKING "
        "WITHOUT VALIDATION BY A QUALIFIED CLINICAL PHARMACOLOGIST OR GENETICIST.\n\n"
        "All recommendations should be reviewed in the context of the patient's complete clinical picture, "
        "co-medications, comorbidities, and other relevant factors. Consult CPIC guidelines at cpicpgx.org "
        "for the most current evidence-based dosing recommendations."
    )
    pdf.multi_cell(0, 5.5, disclaimer)
    
    pdf.ln(5)
    pdf.section_title("REFERENCES", color=(30, 58, 138))
    refs = [
        "1. CPIC Guidelines - cpicpgx.org",
        "2. PharmGKB - pharmgkb.org",
        "3. Relling MV et al. CPIC: Clinical Pharmacogenomics Implementation Consortium guidelines. Clin Pharmacol Ther. 2011.",
        "4. Scott SA et al. CYP2C19 allele nomenclature. Pharmacogenet Genomics. 2012.",
        "5. Caudle KE et al. Standardizing terms for clinical pharmacogenomic test results. Genet Med. 2017.",
    ]
    pdf.set_font("Helvetica", "", 8)
    pdf.set_text_color(51, 65, 85)
    for ref in refs:
        pdf.cell(0, 6, ref, ln=True)

    return bytes(pdf.output())