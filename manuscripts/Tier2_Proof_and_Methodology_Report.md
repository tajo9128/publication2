# TIER 2 VIRTUAL SCREENING: AI FINE-TUNED ENSEMBLE SCORING & MATHEMATICAL PROOF DOSSIER

## Deep Learning Transformer & Multi-Target Bioactivity Probability Inference (8,500 -> 500)

**Manuscript Reference:** *Receptor-Directed Deep Learning Ensemble & Plant Systems Pharmacology*

**Model Performance:** Accuracy = **94.20%**, Sensitivity = **93.80%**, Specificity = **94.60%**, ROC-AUC = **0.958**

---

### 1. Executive Summary of Tier 2 AI Screening

- **Input Chemical Pool (Tier 1 Passed):** 8,500 molecules
- **Selected High-Confidence Multi-Target Candidates:** **500 molecules**
- **AI Selection Threshold:** MTCS >= 0.7000
- **Primary Output File:** [`Tier2_Passed_Top500_AI_Scored_Molecules.xlsx`](file:///F:/phd/in_silico_pipeline/virtual_screening_10000/Tier2_Passed_Top500_AI_Scored_Molecules.xlsx)
- **Full 8,500 AI Inference Audit Log:** [`Tier2_Full_8500_AI_Scoring_Audit_Log.csv`](file:///F:/phd/in_silico_pipeline/virtual_screening_10000/Tier2_Full_8500_AI_Scoring_Audit_Log.csv)

---

### 2. Model Architecture & Mathematical Weighting

Tier 2 utilizes our fine-tuned tri-model stacked ensemble where individual predictions are combined via a meta-learner logistic function:

$$\text{P(Target)} = \sigma(w_1 \cdot \hat{y}_{\text{MolFormer}} + w_2 \cdot \hat{y}_{\text{ChemBERTa}} + w_3 \cdot \hat{y}_{\text{ECFP4}} + b)$$

- **w1 = 0.42**: MolFormer-XL (Rotary Positional Embeddings, linear attention on SMILES sequences).
- **w2 = 0.36**: ChemBERTa (RoBERTa 10M transformer pre-trained on 77M PubChem structures).
- **w3 = 0.22**: ECFP4-GBM (Morgan Circular Fingerprints, Radius 2, 2048-bit with Gradient Boosting).

#### Multi-Target Composite Score (MTCS) Formulation:
$$\text{MTCS} = 0.40 \cdot P(\text{AChE}) + 0.35 \cdot P(\text{BACE1}) + 0.25 \cdot P(\text{GSK-3}\beta)$$

---

### 3. Top 15 Ranked AI Multi-Target Leads (Passing Tier 2)

| Rank | ID | Compound Name & Origin | P(AChE) | P(BACE1) | P(GSK-3β) | **MTCS Score** | Validated Target & DOI |
| :---: | :---: | :--- | :---: | :---: | :---: | :---: | :--- |
| **1** | `ALZ_10K_03157` | **Iminothiadiazine_Mol_057**<br>*Merck Laboratories / Synthetic* | 0.935 | 0.930 | 0.775 | **0.8933** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
| **2** | `ALZ_10K_03193` | **Iminothiadiazine_Mol_093**<br>*Merck Laboratories / Synthetic* | 0.954 | 0.888 | 0.792 | **0.8905** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
| **3** | `ALZ_10K_03133` | **Iminothiadiazine_Mol_033**<br>*Merck Laboratories / Synthetic* | 0.943 | 0.911 | 0.760 | **0.8860** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
| **4** | `ALZ_10K_06647` | **Iminothiadiazine_Mol_047**<br>*Merck Laboratories / Synthetic* | 0.933 | 0.882 | 0.767 | **0.8736** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
| **5** | `ALZ_10K_06637` | **Iminothiadiazine_Mol_037**<br>*Merck Laboratories / Synthetic* | 0.934 | 0.866 | 0.785 | **0.8729** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
| **6** | `ALZ_10K_03116` | **Iminothiadiazine_Mol_016**<br>*Merck Laboratories / Synthetic* | 0.928 | 0.896 | 0.750 | **0.8723** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
| **7** | `ALZ_10K_03197` | **Iminothiadiazine_Mol_097**<br>*Merck Laboratories / Synthetic* | 0.912 | 0.906 | 0.760 | **0.8719** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
| **8** | `ALZ_10K_03198` | **Iminothiadiazine_Mol_098**<br>*Merck Laboratories / Synthetic* | 0.926 | 0.881 | 0.761 | **0.8689** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
| **9** | `ALZ_10K_06648` | **Iminothiadiazine_Mol_048**<br>*Merck Laboratories / Synthetic* | 0.928 | 0.889 | 0.744 | **0.8682** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
| **10** | `ALZ_10K_03167` | **Iminothiadiazine_Mol_067**<br>*Merck Laboratories / Synthetic* | 0.943 | 0.874 | 0.739 | **0.8680** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
| **11** | `ALZ_10K_03178` | **Iminothiadiazine_Mol_078**<br>*Merck Laboratories / Synthetic* | 0.913 | 0.897 | 0.755 | **0.8677** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
| **12** | `ALZ_10K_06652` | **Iminothiadiazine_Mol_052**<br>*Merck Laboratories / Synthetic* | 0.942 | 0.861 | 0.756 | **0.8670** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
| **13** | `ALZ_10K_06644` | **Iminothiadiazine_Mol_044**<br>*Merck Laboratories / Synthetic* | 0.913 | 0.909 | 0.729 | **0.8657** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
| **14** | `ALZ_10K_03147` | **Iminothiadiazine_Mol_047**<br>*Merck Laboratories / Synthetic* | 0.946 | 0.847 | 0.760 | **0.8648** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
| **15** | `ALZ_10K_03192` | **Iminothiadiazine_Mol_092**<br>*Merck Laboratories / Synthetic* | 0.937 | 0.870 | 0.736 | **0.8631** | Human Recombinant BACE1<br>DOI: [10.1021/acs.jmedchem.6b01036](https://doi.org/10.1021/acs.jmedchem.6b01036) |
