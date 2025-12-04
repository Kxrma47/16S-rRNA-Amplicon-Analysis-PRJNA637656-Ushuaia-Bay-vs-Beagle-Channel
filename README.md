# 16S rRNA Metagenomics Analysis — Ushuaia Bay vs Beagle Channel  
BioProject **PRJNA637656**

This repository contains a complete, reproducible workflow for processing and analyzing  
**16S rRNA V4 amplicon data** from seawater samples collected in **Ushuaia Bay** and the **Beagle Channel** (Argentina).  
The work is divided into two stages:

- **Part 1** — ASV inference using DADA2, taxonomy assignment (SILVA 138.1), generation of a phyloseq object.  
- **Part 2** — Genus-level filtering, visualization, alpha- and beta-diversity analyses.

The dataset is described in the article:

**“The 16S rRNA gene amplicon data set from seawater in Ushuaia Bay and the Beagle Channel (Argentina)”**  
ScienceDirect: https://www.sciencedirect.com/science/article/pii/S2352340920310659

---

## Dataset Links
- NCBI BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA637656  
- SRA runs (SRR11941686–SRR11941695): https://www.ncbi.nlm.nih.gov/sra?term=PRJNA637656  
- ENA FASTQ browser: https://www.ebi.ac.uk/ena/browser/view/PRJNA637656  

All FASTQs and metadata used in this analysis were automatically downloaded from ENA.

---

# Part 1 — ASV Pipeline (Short Report)

### 1. Data Acquisition
- Downloaded all **10 paired-end FASTQ samples**.
- Constructed metadata table containing sample_id, SRA/ENA accessions, platform, library layout, etc.
- Stored raw data under `data/` and cached metadata as `work/run_tbl.rds`.

### 2. Primer / Adapter Detection
- Extracted 2000 reads/sample using *ShortRead*.
- Detected Illumina + 16S V4 primers → trimming required.

### 3. Adapter / Primer Removal
- Used **cutadapt** with `-m 50` minimum length.
- Trimmed reads saved to `work/cutadapt/`.

### 4. Quality Profiles
- Generated `plotQualityProfile()` for forward/reverse reads.
- Standard Illumina curves observed: gradual decline toward read ends.

### 5. DADA2 Workflow
Parameters: `truncLen = c(150,150)`, `maxEE = c(3,3)`

Steps executed:
- `filterAndTrim()`
- `learnErrors()` for F/R
- `dada()` denoising
- `mergePairs()` (low merge rate expected due to short overlap)
- `makeSequenceTable()`
- `removeBimeraDenovo()`

**Final result:**  
➡️ **21 non-chimeric ASVs** across 10 samples.

Outputs saved:
- `work/read_tracking.tsv`  
- `work/reads_per_sample_final.tsv`  
- `work/asv_per_sample.tsv`

### 6. Taxonomy Assignment (SILVA 138.1)
- Used SILVA train + species files.
- Several ASVs remained **unclassified at genus level** — typical for marine datasets.

### 7. Phyloseq Object
The following were constructed and bundled:
- OTU table (**21 taxa × 10 samples**)  
- Taxonomy table (7 ranks)  
- Sample metadata (10 samples × 6 variables)  
- DNAStringSet reference sequences  

ASVs were renamed to **sp1–sp21**, exported as:
- `ASVs_sp.fasta`
- `work/phyloseq_object.rds`

### **Outcome of Part 1**
A complete ASV-processing pipeline was reproduced for PRJNA637656, yielding a high-confidence, SILVA-annotated dataset packaged as a **phyloseq** object ready for diversity analyses.

---

# Part 2 — Genus-Level Analyses (Representation, Diversity, Statistics)

## Task 1 — Selecting the Most Represented Genera
- Agglomerated ASVs → **Genus** level using `tax_glom`.
- Initial phyloseq (genus level): **11 genera × 10 samples**.
- Filtering rules (per assignment):
  - Present in **≥ 50%** of samples  
  - Has **≥ 5%** relative abundance in ≥ 1 sample  

After filtering:  
➡️ **5 genera** remained.

Both "before" and "after" phyloseq summaries were printed as required.

---

## Task 2 — Visualization  
### **Figure 1 – Relative Abundance**
Shows genera retained after filtering. The community is visually dominated by **Unclassified** genus, meaning genus-level ecological interpretation is limited.

### **Figure 2 – CLR Heatmap**
CLR calculation:  
`CLR = log(count + 1) – mean(log(count + 1))` per sample.

This accounts for compositional constraints and highlights genus-level variation.

---

## Task 3 — Alpha-Diversity + Statistical Testing
Computed per sample:
- Observed richness  
- Shannon diversity  
- Simpson diversity  

**Figure 3** shows boxplots by site (Beagle Channel vs Ushuaia Bay).

**Wilcoxon tests:**
| Metric   | p-value |
|----------|---------|
| Observed | 0.828   |
| Shannon  | 0.906   |
| Simpson  | 0.105   |

➡️ **No statistically significant differences** in alpha diversity.

---

## Task 4 — Beta-Diversity (Bray–Curtis)  
### Distances + PCoA  
**Figure 4**: PCoA of Bray–Curtis distances.  
→ Samples from the two sites overlap; no clear community separation.

### Heatmap  
**Figure 5**: Bray–Curtis matrix heatmap with hierarchical clustering.

### Dispersion Test (betadisper)
- Tests differences in within-site variability.
- Result: **p = 0.5347** → not significant.

### PERMANOVA (adonis2)
Model: `bray ~ site`  
Results:
- **p = 0.7333**  
- **R² = 0.19**  

➡️ No significant difference in community composition between sites.

### Interpretation of Methods
- **adonis2** tests *centroid separation* between groups.  
- **betadisper** tests *variance differences* within groups.  
Both are required to ensure correct interpretation.

---

# Summary & Biological Interpretation

- The dataset contains **21 ASVs**, collapsing to **11 genera**, and finally **5 genera** after filtering.  
- Genus-level taxonomic resolution is limited (many ASVs remain unclassified).  
- **No significant ecological differences** were observed between Ushuaia Bay and the Beagle Channel in:
  - α-diversity  
  - β-diversity dispersion  
  - Community composition (PERMANOVA)  

Given small sample size (*n* = 10) and the dominance of unclassified genera, subtle ecological signals might remain undetected.  
Nevertheless, all required analytical steps were performed and documented as specified.

---

# Reproducibility
- Part 1 code outputs `phyloseq_object.rds` for reuse.
- Part 2 uses the saved object to generate all visualizations and statistics.

---

# Citation
If using this analysis, cite:
- The Data in Brief article (Ushuaia Bay / Beagle Channel dataset)
- SILVA 138.1  
- DADA2 (Callahan et al.)  
- phyloseq (McMurdie & Holmes)

---

# License
Academic / educational use.
