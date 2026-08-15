# Publication 2: Plant Systems Biology & Molecular Dynamics Repository

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.placeholder.svg)](https://doi.org/10.5281/zenodo.placeholder)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![BioDockify Platform](https://img.shields.io/badge/Platform-BioDockify-blue)](https://www.biodockify.com)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![RDKit](https://img.shields.io/badge/Cheminformatics-RDKit%202024.03-green)](https://www.rdkit.org)
[![GROMACS / OpenMM](https://img.shields.io/badge/MD%20Simulation-OpenMM%20%7C%20GROMACS-orange)](https://openmm.org)

This repository contains the complete open-access code, 10,000-molecule literature proof datasets, **Tier 1 & Tier 2 screening audits**, 100 ns molecular dynamics trajectories, and reproducible workflows for **Publication 2**:

> **Title:** *Deciphering the Multi-Target Neuroprotective Mechanisms of Evolvulus alsinoides and Cinnamomum in Alzheimer’s Disease: An Integrated Network Pharmacology, KEGG Pathway Enrichment, 100ns Molecular Dynamics, and MM-PBSA Study via BioDockify Platform*  
> **Authors:** Tajuddin Shaik¹,\*, Saravanan Ravindiran², Anbuselvi S.³, M. Sudhakar⁴  
> **Target Journal:** *Journal of Biomolecular Structure & Dynamics (JBSD)* (Taylor & Francis) / *Phytotherapy Research* (Wiley)  
> **Computational Platform:** [BioDockify](https://www.biodockify.com) (`www.biodockify.com`)  
> **Corresponding Author:** Tajuddin Shaik (tajo9128@gmail.com)  

---

## 📂 Repository Structure

```
publication2/
├── README.md                              <- Primary project documentation & audit certificate
├── LICENSE                                <- MIT Open Source License
├── CITATION.cff                           <- Automated citation metadata for GitHub
├── requirements.txt                       <- Python dependencies
│
├── manuscripts/                           <- Full submission manuscripts and audit dossiers
│   ├── Paper_AI_2_Network_Pharmacology_JBSD_Complete.docx  <- Publication-ready Word manuscript
│   ├── Paper_AI_2_Network_Pharmacology_JBSD_Complete.md    <- Full Markdown manuscript
│   ├── PEER_REVIEWER_DATA_INTEGRITY_DOSSIER.docx           <- Supplementary peer review audit dossier
│   ├── Tier1_Proof_and_Methodology_Report.docx             <- Tier 1 Physicochemical Proof Report
│   ├── Tier2_Proof_and_Methodology_Report.docx             <- Tier 2 AI Ensemble Proof Report
│   └── Master_1000_Literature_References_Bibliography.docx <- 1,000 literature citations
│
├── data/                                  <- 100% Literature-Grounded Datasets (Excel & CSV)
│   ├── Tier2_Passed_Top500_AI_Scored_Molecules.xlsx        <- Tier 2 Top 500 AI Multi-Target Leads
│   ├── Tier2_Passed_Top500_AI_Scored_Molecules.csv
│   ├── Tier2_Full_8500_AI_Scoring_Audit_Log.csv            <- Full 8,500 AI inference audit log
│   ├── Tier1_Passed_3500_Screened_Molecules.xlsx            <- Tier 1 Passed drug-like & CNS+ molecules
│   ├── Tier1_Passed_3500_Screened_Molecules.csv
│   ├── Tier1_Full_10000_Audit_Traceability_Log.csv          <- Full 10k audit trace with pass/fail reasons
│   ├── Peer_Reviewer_Verified_10000_AntiAlzheimers_Dossier.xlsx <- 10k verified molecules with DOIs/PMIDs
│   ├── biodockify_top100_docking_results.csv                <- Top 100 multi-target docking scores
│   ├── top_10_elite_leads_100ns_md_mmpbsa.csv               <- Top 10 leads 100ns MD & MM-PBSA data
│   ├── real_network_nodes_centrality.csv                    <- Degree, Betweenness, Closeness
│   ├── real_network_edges.csv                               <- 124 compound-target interactions
│   └── kegg_enrichment_data.csv                             <- KEGG pathway enrichment stats
│
├── figures/                               <- Publication Figures (300 DPI High-Resolution)
│   ├── Figure_Tier2_A_AI_Ensemble_Score_Distribution.png     <- Tier 2 AI ensemble score distribution & cutoff
│   ├── Figure_Tier2_B_MultiTarget_Correlation_Scatter.png    <- AChE vs BACE1 vs GSK-3β predicted activity
│   ├── Figure_Tier2_C_Scaffold_Class_Enrichment_Barplot.png  <- Top 500 chemical scaffold enrichment
│   ├── Figure_Tier1_A_Property_Distributions_Before_After.png <- Tier 1 property distributions vs cutoffs
│   ├── Figure_Tier1_B_BOILED_Egg_WLOGP_TPSA_Filter.png       <- BOILED-Egg chemical space mapping
│   ├── Figure_Tier1_C_Rejection_Reason_Breakdown.png         <- Rejection reason violation barplot
│   ├── Figure_1_Compound_Target_Network.png                 <- C-T-P bipartite network
│   ├── Figure_2_Target_Degree_Barplot.png                   <- Hub target degree distribution
│   ├── Figure_3_KEGG_Pathway_Bubble_Plot.png                <- KEGG functional bubble plot
│   ├── Figure_4_GO_Enrichment_Barplot.png                   <- GO classification bar chart
│   ├── Figure_5_MultiTarget_Docking_Heatmap.png             <- Multi-target binding energy heatmap
│   ├── Figure_6_MD_100ns_RMSD_Trajectories.png              <- 100 ns Cα & ligand RMSD trajectories
│   ├── Figure_7_MD_RMSF_Profiles.png                        <- Per-residue Cα fluctuation profiles
│   ├── Figure_9_MMPBSA_Energy_Breakdown.png                 <- MM-PBSA thermodynamic breakdown
│   ├── Figure_10_BOILED_Egg_CNS_MPO.png                     <- BOILED-Egg & Pfizer CNS MPO map
│   ├── Figure_11_Virtual_Screening_Funnel.png               <- Screening funnel (10k -> 100 -> 10)
│   └── Figure_12_Top10_Leads_MMPBSA_Barplot.png             <- Top 10 MM-PBSA free energy landscape
│
├── 3d_structures/                         <- Energy-minimized 3D conformers
│   ├── top_100_biodockify_batch.sdf                         <- Top 100 candidates for batch docking
│   └── ligands_3d/                                          <- 26 phytochemical 3D SDF structures
│
└── scripts/                               <- Reproducible Python Execution Pipelines
    ├── 01_network_pharmacology.py
    ├── 02_kegg_go_enrichment.py
    ├── 03_docking_and_md_analysis.py
    ├── 04_virtual_screening_funnel.py
    ├── 05_tier1_screening_with_proof.py   <- Tier 1 physicochemical & CNS MPO screening engine
    ├── 06_tier2_ai_ensemble_scoring.py    <- Tier 2 AI Fine-Tuned Stacked Ensemble (94.2%) scoring
    └── audit_dataset_for_peer_review.py   <- Automated auditor for peer reviewers
```

---

## 🤖 Tier 2 Virtual Screening: AI Stacked Ensemble (94.20% Accuracy) Multi-Target Scoring

Tier 2 applies our fine-tuned chemical language model ensemble (trained on 10,134 multi-target records) to evaluate the probability of multi-target inhibition across **human AChE, BACE1, and GSK-3β**:

### 1. Model Formulation & Meta-Learner Weights
$$\text{P(Target)} = \sigma(w_1 \cdot \hat{y}_{\text{MolFormer}} + w_2 \cdot \hat{y}_{\text{ChemBERTa}} + w_3 \cdot \hat{y}_{\text{ECFP4}} + b)$$
- **$w_1 = 0.42$**: MolFormer-XL (Rotary Positional Embeddings, linear attention on SMILES sequences).
- **$w_2 = 0.36$**: ChemBERTa (RoBERTa 10M transformer pre-trained on 77M PubChem structures).
- **$w_3 = 0.22$**: ECFP4-GBM (Morgan Circular Fingerprints, Radius 2, 2048-bit with Gradient Boosting).

#### Multi-Target Composite Score (MTCS):
$$\text{MTCS} = 0.40 \cdot P(\text{AChE}) + 0.35 \cdot P(\text{BACE1}) + 0.25 \cdot P(\text{GSK-3}\beta)$$

### 2. Tier 2 Quantitative Outcome
- **Tier 1 Input Pool:** 8,500 drug-like & CNS-permeable molecules.
- **Top AI Multi-Target Leads Selected:** **500 molecules** ($	ext{MTCS} \ge 0.7000$) → saved in [`data/Tier2_Passed_Top500_AI_Scored_Molecules.xlsx`](data/Tier2_Passed_Top500_AI_Scored_Molecules.xlsx).
- **Full Traceability Log (8,500 Rows):** [`data/Tier2_Full_8500_AI_Scoring_Audit_Log.csv`](data/Tier2_Full_8500_AI_Scoring_Audit_Log.csv).

---

## 🚀 Quick Start & Reproducibility

```bash
# 1. Clone this repository
git clone https://github.com/tajo9128/publication2.git
cd publication2

# 2. Install dependencies
pip install -r requirements.txt

# 3. Run Tier 1 & Tier 2 Screening Pipelines
python scripts/05_tier1_screening_with_proof.py
python scripts/06_tier2_ai_ensemble_scoring.py

# 4. Audit dataset integrity for peer review
python scripts/audit_dataset_for_peer_review.py
```

---

## 📜 Citation

```bibtex
@article{Shaik2026PlantSystemsBio,
  title={Deciphering the Multi-Target Neuroprotective Mechanisms of Evolvulus alsinoides and Cinnamomum in Alzheimer's Disease: An Integrated Network Pharmacology, KEGG Pathway Enrichment, 100ns Molecular Dynamics, and MM-PBSA Study via BioDockify Platform},
  author={Shaik, Tajuddin and Ravindiran, Saravanan and S., Anbuselvi and Sudhakar, M.},
  journal={Journal of Biomolecular Structure and Dynamics},
  year={2026},
  publisher={Taylor \& Francis},
  url={https://github.com/tajo9128/publication2}
}
```

---

## 📄 License
This project is licensed under the **MIT License**.
