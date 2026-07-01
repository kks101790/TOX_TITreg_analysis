# Data & Code Availability — final statements (TOX TI-Treg, Nature Immunology)

Accessions final: GEO **GSE336236** (private; reviewer token issued), Zenodo concept DOI
**10.5281/zenodo.20810031**, GitHub **kks101790/TOX_TITreg_analysis**.
csaw confirmed NOT used in any figure (template only) → removed. Pipeline tools verified from
FASTQ headers / BAM @PG / Cell Ranger web_summary.

## Data Availability
The single-cell RNA-sequencing, single-cell V(D)J-sequencing, and bulk ATAC-sequencing data
generated in this study have been deposited in the Gene Expression Omnibus (GEO) under accession
GSE336236. Previously published single-cell RNA-sequencing data of PD-1-deficient
tumor-infiltrating regulatory T cells used for the integrative analysis are available under GEO
accession GSE164033. Publicly available single-cell RNA-sequencing data of CD4+ T cells from
normal and tumour-bearing lung are available under GEO accession GSE152022. Mass cytometry
(CyTOF) data (Fig. 1a,b) are available from the corresponding authors upon reasonable request.
Source data are provided with this paper.

## Code Availability
Custom analysis code for the single-cell RNA-seq, scVDJ-seq, and bulk ATAC-seq analyses is
available at GitHub (https://github.com/kks101790/TOX_TITreg_analysis) and archived at Zenodo
(https://doi.org/10.5281/zenodo.20810031). Sequencing reads were processed with Cell Ranger
v7.0.1 (10x Genomics; mm10/refdata-gex-mm10-2020-A) for the scRNA-seq/scVDJ-seq/Feature Barcode
libraries, and with Bowtie2 v2.5.4 and SAMtools v1.21 (alignment to mm10) followed by MACS2 v2.2
(broad peak calling, ENCODE blacklist filtering) for the ATAC-seq libraries. Downstream analyses
used Seurat v4.0.5, GSVA v1.38.2, clusterProfiler v3.18.1, ggvenn v0.1.10 and pheatmap v1.0.12
(scRNA-seq/scVDJ-seq; R 4.0.5); and ChIPseeker v1.40.0, trackViewer v1.40.0, clusterProfiler
v4.12.6 and TxDb.Mmusculus.UCSC.mm10.knownGene v3.10.0 (ATAC-seq; R 4.4.1), together with TOBIAS
v0.17.0 (ATACorrect/FootprintScores/BINDetect) using JASPAR2024 CORE non-redundant motifs, and
ggplot2/ggrepel for visualization. No previously unreported custom algorithms were used.

## Methods — data-processing sentences to ADD (these steps are currently missing from Methods)
scRNA/scVDJ (add near the start of the single-cell processing paragraph):
  "Raw sequencing reads were demultiplexed, aligned to the mouse genome (mm10,
   refdata-gex-mm10-2020-A) and quantified using Cell Ranger v7.0.1 (10x Genomics); the resulting
   filtered feature-barcode matrices were analysed in R with Seurat."
ATAC-seq (add before the MACS2 sentence):
  "ATAC-seq reads were aligned to mm10 with Bowtie2 v2.5.4 (--very-sensitive -X 1000); alignments
   were sorted and filtered with SAMtools v1.21 (retaining properly paired, non-mitochondrial
   reads), Tn5-shifted, and used for MACS2 broad peak calling with ENCODE-blacklist filtering."

## Notes (submission, NOT published text)
- GSE336236 private; give the GEO reviewer access token to the editor in the cover letter.
- Release GSE336236 publicly upon acceptance/publication.
- csaw removed from repo + Code Availability (template `ATAC-seq-master/csaw_workflow.R`; produced no
  figure, outputs never on disk, absent from Methods). REMOVE scripts/03 from the GitHub repo:
  `git rm scripts/03_ATAC_differential_csaw.R && git commit -m "remove unused csaw script" && git push`
- bigWig coverage-track tool not recoverable from scripts; described generically as coverage tracks.
- CyTOF = collaborators (S.P./Mt Sinai), "upon reasonable request".
