# Nature Research — Code and Software Submission Checklist (completed)

**Manuscript:** NI-A43611 — "TOX enforces the immunosuppressive program of tumor-infiltrating regulatory T cells"
**Corresponding author(s):** Kyungsoo Kim; Sang-Jun Ha

## Nature of the code
This study uses **analysis code built entirely on standard, published, open-source packages**
(Seurat, GSVA, clusterProfiler, ggvenn, pheatmap, ChIPseeker, trackViewer, TOBIAS, MACS2, Bowtie2, SAMtools,
Cell Ranger). **No previously unreported custom algorithm or standalone software tool** is introduced; the
scripts are workflows that reproduce specific figures from the deposited data.

## Access for editors and reviewers (single link)
All required content is openly accessible (no compiled binary; source code + documentation):
- **GitHub (source code + README):** https://github.com/kks101790/TOX_TITreg_analysis
- **Zenodo (archived, versioned; downloadable as a single archive):** https://doi.org/10.5281/zenodo.20810031
- **Associated data — GEO GSE336236** (private during review; reviewer token: `kbirycuuzvmfnsn`)

*Please confirm the GitHub repository visibility is set to **Public** (the Zenodo archive is public regardless).*

## Required content — checklist
- [x] **Source code** — provided (R scripts, `scripts/01–05`). No compiled/standalone binary (analysis code).
- [x] **A small dataset to demo the code** — no bundled simulated dataset; the demo runs on public GEO data
      **GSE152022** (GSM4598898, GSM4598899; Fig 1c–f). In-house data are at **GSE336236**. Demo command,
      expected output and expected run time are given in **README §3**.
- [x] **README file** that includes:
  - [x] **1. System requirements** — OS (Linux, x86-64), R 4.0.5 / 4.4.1, all package versions
        (`envs/KEY_VERSIONS.txt`), hardware (multi-core workstation, ≥32 GB RAM, no non-standard hardware/GPU). — README §1
  - [x] **2. Installation guide** — two `conda env create` commands; typical install time ~10–20 min per
        environment on a normal desktop. — README §2
  - [x] **3. Demo** — instructions to run, expected output, and expected run time (~10–20 min for Fig 1;
        ~1–2 h for script 01 per-cell GSVA; ~20–40 min for script 04). — README §3
  - [x] **4. Instructions for use** — how to run each script on your own data (edit paths at top, `Rscript`). — README §4
  - [x] **(Optional) Reproduction instructions** — script→figure map; running each script on the corresponding
        deposited/public data regenerates the listed panels (byte-identical PDFs confirmed for ATAC locus
        tracks and motif volcano). — README §4

## Additional information
- **License:** MIT — an Open Source Initiative-approved license (`LICENSE` file in the repository).
- **Link to code in an open-source repository:** GitHub (above); archived at Zenodo DOI 10.5281/zenodo.20810031.
- **Detailed description of the code's functionality (pseudocode):** provided in the manuscript **Methods**
  (single-cell processing paragraphs and "Bulk ATAC sequencing analysis"), and summarized as a script→figure
  table in the README. No novel algorithm/pseudocode is required, as only standard published packages were used
  (cited with versions in Methods and the Code Availability statement).
- **Code Ocean capsule:** not created (optional; not required for this submission).
