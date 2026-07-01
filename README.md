# TOX TI-Treg — single-cell & ATAC-seq analysis code

[![DOI](https://zenodo.org/badge/1277810034.svg)](https://doi.org/10.5281/zenodo.20810031)

Analysis code for **"TOX enforces the immunosuppressive program of tumor-infiltrating regulatory T cells"**
(Park et al., *Nature Immunology*). This repository contains the single-cell RNA-seq, single-cell V(D)J,
and bulk ATAC-seq analysis used to generate the computational figures.

This is **analysis code built on standard, published packages** (Seurat, GSVA, clusterProfiler, ChIPseeker,
trackViewer, TOBIAS, MACS2, etc.) — **no previously unreported custom algorithm** is introduced. Each script
is a workflow that reproduces specific published figures from the deposited data.

> Mass-cytometry (CyTOF; Fig 1a,b) and flow-cytometry (FlowJo) analyses are not included here — they were
> performed by other authors. This repository covers the scRNA-seq, scVDJ-seq, and ATAC-seq analyses.

## Repository structure
```
scripts/
  01_scRNA_scVDJ_analysis.R        # Fig 6, Extended Data 7-8  (Seurat: HTO demux, FindMarkers,
                                   #   GSVA, GSEA, PD-1 GSE164033 integration, Venn, heatmap)
  02_public_CD4_GSE152022.R        # Fig 1c-f  (public CD4 T-cell reanalysis)
  03_ATAC_differential_csaw.R      # exploratory csaw differential-accessibility workflow
                                   #   (NOT used for any published figure; provided for completeness)
  04_ATAC_peak_annotation_tracks.R # Fig 7b,7c  (ChIPseeker annotation + trackViewer locus tracks)
  05_ATAC_TOBIAS_footprint.R       # Fig 7d, Extended Data 9  (TOBIAS BINDetect motif footprinting)
envs/
  scRNA_env_R-4.0.0.yml            # conda environment (scRNA/scVDJ)
  ATAC_env_ATACseq.yml             # conda environment (ATAC)
  KEY_VERSIONS.txt                 # key package versions (as cited in Methods)
```

## 1. System requirements
- **Operating system:** Linux (developed and tested on Ubuntu 20.04 / x86-64). Any OS that runs conda +
  R should work; not tested on Windows/macOS.
- **R + key packages** (exact versions in `envs/KEY_VERSIONS.txt`):
  - scRNA/scVDJ — R 4.0.5; Seurat 4.0.5, GSVA 1.38.2, clusterProfiler 3.18.1, ggvenn 0.1.10, pheatmap 1.0.12.
  - ATAC — R 4.4.1; ChIPseeker 1.40.0, trackViewer 1.40.0, csaw 1.38.0,
    TxDb.Mmusculus.UCSC.mm10.knownGene 3.10.0, clusterProfiler 4.12.6.
  - Command-line — Cell Ranger 7.0.1, Bowtie2 2.5.4, SAMtools 1.21, MACS2 2.2, TOBIAS 0.17.0 (JASPAR2024), mm10.
- **Hardware:** a standard multi-core workstation. No non-standard hardware (no GPU) is required. Some steps
  are memory-intensive (Seurat integration, per-cell GSVA with `parallel.sz = 28`); **≥ 32 GB RAM and
  ~8+ cores are recommended.**

## 2. Installation guide
Recreate the two conda environments:
```bash
conda env create -f envs/scRNA_env_R-4.0.0.yml   # R 4.0.5 — run scripts 01-02
conda env create -f envs/ATAC_env_ATACseq.yml    # R 4.4.1 — run scripts 03-05
```
**Typical install time:** ~10–20 min per environment on a normal desktop (dominated by conda/Bioconductor
package downloads and solving).

## 3. Demo
No separate simulated dataset is bundled; the scripts run on the study's deposited data and on public GEO data:
- **In-house data:** GEO SuperSeries **GSE336236** (processed matrices / peak & coverage files).
- **Public data:** **GSE152022** (Fig 1) and **GSE164033** (PD-1 WT/KO integration, Extended Data 7-8).

Minimal demo (Fig 1c-f): download the two `GSE152022` samples (GSM4598898, GSM4598899), set the input path at
the top of `02_public_CD4_GSE152022.R`, then:
```bash
conda activate <scRNA env>
Rscript scripts/02_public_CD4_GSE152022.R
```
**Expected output:** integrated UMAP, DotPlot, VlnPlot and a Seurat object (`.rds`) corresponding to Fig 1c-f.
**Expected run time:** ~10–20 min on a normal desktop. (Script 01 with per-cell GSVA is heavier: ~1–2 h;
script 04 ChIPseeker binning + trackViewer: ~20–40 min.)

## 4. Instructions for use
Scripts 01–02 use the scRNA environment; 03–05 use the ATAC environment. Input/output paths are set at the
top of each script — edit them to your local directory layout, then run with `Rscript scripts/<name>.R`.
Expected inputs:
- **scRNA/scVDJ:** Cell Ranger outputs (filtered_feature_bc_matrix; TCR `filtered_contig_annotations.csv`;
  Feature-Barcode / hashtag matrix; hashtags C0301 = YFP⁺ TOXKO, C0302 = YFP⁻ TOXWT).
- **ATAC:** aligned BAM → MACS2 broad peaks (blacklist-filtered) → ChIPseeker / trackViewer / TOBIAS.
- **Public:** GSE152022 (Fig 1) and GSE164033 (PD-1 WT/KO integration) downloaded from GEO.

### Reproduction instructions (script → figure)
| Script | Figures |
|--------|---------|
| `01_scRNA_scVDJ_analysis.R` | Fig 6; Extended Data 7–8 |
| `02_public_CD4_GSE152022.R` | Fig 1c–f |
| `04_ATAC_peak_annotation_tracks.R` | Fig 7b, 7c |
| `05_ATAC_TOBIAS_footprint.R` | Fig 7d; Extended Data 9 |

Running each script on the corresponding deposited/public data regenerates the listed figure panels
(byte-identical PDFs were confirmed for the ATAC locus tracks and motif volcano).

## Data availability
- In-house **scRNA-seq + scVDJ-seq + ATAC-seq**: GEO SuperSeries **GSE336236**.
- Reused public: **GSE152022**, **GSE164033** (PD-1 WT/KO; Kim et al., *Nat Immunol* 2023) — cited, not re-deposited.

## Notes / caveats
- ATAC-seq has **one library per condition** (Tox_FF = WT, Tox_cKO = KO; no biological replicate) →
  differential-accessibility / footprinting results are interpreted descriptively.
- scRNA WT-vs-KO differential expression is within a single hashtag-multiplexed library (cell-level test).
- Analysis logic in each script is as used for the published figures; only file paths are user-configurable.

## Citation
If you use this code, please cite Park et al., *Nature Immunology* (year), and this repository
(archived on Zenodo, DOI **10.5281/zenodo.20810031**).

## License
MIT (see `LICENSE`) — an [Open Source Initiative](https://opensource.org/licenses/MIT)-approved license.
