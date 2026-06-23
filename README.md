# TOX TI-Treg — single-cell & ATAC-seq analysis code

[![DOI](https://zenodo.org/badge/1277810034.svg)](https://doi.org/10.5281/zenodo.20810031)

Analysis code for **"TOX enforces the immunosuppressive program of tumor-infiltrating regulatory T cells"**
(Park et al., *Nature Immunology*). This repository contains the single-cell RNA-seq, single-cell V(D)J,
and bulk ATAC-seq analysis used to generate the computational figures.

> Mass-cytometry (CyTOF; Fig 1a,b) and flow-cytometry (FlowJo) analyses are not included here — they were
> performed by other authors. This repository covers the scRNA-seq, scVDJ-seq, and ATAC-seq analyses.

## Repository structure
```
scripts/
  01_scRNA_scVDJ_analysis.R        # Fig 6, Extended Data 7-8  (Seurat: HTO demux, FindMarkers,
                                   #   GSVA, GSEA, PD-1 GSE164033 integration, Venn, heatmap)
  02_public_CD4_GSE152022.R        # Fig 1c-f  (public CD4 T-cell reanalysis)
  03_ATAC_differential_csaw.R      # Fig 7 framework  (csaw differential accessibility)
  04_ATAC_peak_annotation_tracks.R # Fig 7b,7c  (ChIPseeker annotation + trackViewer locus tracks)
  05_ATAC_TOBIAS_footprint.R       # Fig 7d, Extended Data 9  (TOBIAS BINDetect motif footprinting)
envs/
  scRNA_env_R-4.0.0.yml            # conda environment (scRNA/scVDJ)
  ATAC_env_ATACseq.yml             # conda environment (ATAC)
  KEY_VERSIONS.txt                 # key package versions (as cited in Methods)
```

## Software environments
Two conda environments were used. Recreate with:
```bash
conda env create -f envs/scRNA_env_R-4.0.0.yml   # R 4.0.5 — run scripts 01-02
conda env create -f envs/ATAC_env_ATACseq.yml    # R 4.4.1 — run scripts 03-05
```
Key versions: Seurat 4.0.5, GSVA 1.38.2, clusterProfiler 3.18.1, ggvenn 0.1.10, pheatmap 1.0.12 (scRNA);
ChIPseeker 1.40.0, trackViewer 1.40.0, csaw 1.38.0, TxDb.Mmusculus.UCSC.mm10.knownGene 3.10.0 (ATAC);
MACS2 2.2.x, TOBIAS 0.17.0, JASPAR2024, genome mm10. 10x Cell Ranger (Single Cell 5′ v2) for demux.

## Running
Scripts 01–02 use the scRNA environment; 03–05 use the ATAC environment. Input/output paths are set at the
top of each script — edit them to your local directory layout. Expected inputs:
- **scRNA/scVDJ**: Cell Ranger outputs (filtered_feature_bc_matrix; TCR `filtered_contig_annotations.csv`;
  Feature-Barcode / hashtag matrix; hashtags C0301 = YFP⁺ TOXKO, C0302 = YFP⁻ TOXWT).
- **ATAC**: aligned BAM → MACS2 broad peaks (blacklist-filtered) → ChIPseeker / trackViewer / TOBIAS.
- **Public**: GSE152022 (Fig 1) and GSE164033 (PD-1 WT/KO integration, Extended Data 7-8) downloaded from GEO.

## Data availability
- In-house **scRNA-seq + scVDJ-seq + ATAC-seq**: GEO SuperSeries **GSE_XXXXXX** (assigned on deposition).
- Reused public: **GSE152022**, **GSE164033** (PD-1 WT/KO; Kim et al., *Nat Immunol* 2023) — cited, not re-deposited.

## Notes / caveats
- ATAC-seq has **one library per condition** (Tox_FF = WT, Tox_cKO = KO; no biological replicate) →
  differential-accessibility / footprinting results are interpreted descriptively.
- scRNA WT-vs-KO differential expression is within a single hashtag-multiplexed library (cell-level test).
- Analysis logic in each script is as used for the published figures; only file paths are user-configurable.

## Citation
If you use this code, please cite Park et al., *Nature Immunology* (year), and this repository
(archived on Zenodo, DOI listed in the release).

## License
MIT (see `LICENSE`).
