# --- explicit package loading (added for standalone reproducibility; analysis logic UNCHANGED;
#     original scripts assumed an interactive session with these already attached) ---
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(stringr))
suppressMessages(library(seqinr))

# =============================================================================
# ATAC-seq TF footprinting / motif (TOBIAS BINDetect)
# Manuscript: "TOX enforces the immunosuppressive program of tumor-infiltrating
#              regulatory T cells" (Park et al., Nature Immunology)
# Figures: Fig 7d; Extended Data 9
# Environment: ATAC env (R 4.4.1) + TOBIAS 0.17.0   (exact versions in ../envs/ and README.md)
# Provenance: original project script "Analyzing_ATACseq_3.R" — analysis logic UNCHANGED; only this
#             header added. Configure input/output paths per README directory layout.
# =============================================================================

# TOBIAS 

# TOBIAS BINDetect --motifs JASPAR2024_CORE_non-redundant_pfms_jaspar.txt --signals Tox_FF_footprints.bw Tox_cKO_footprints.bw --genome mm10.fa.gz --peaks Tox_merged_peaks.bed --peak_header merged_peaks_annotated_header.txt --outdir BINDetect_output --cond_names WT KO --cores 28

# Plot only significant genes
bindetect <- read.delim("BINDetect_output/bindetect_highlight.tsv")
bindetect$WT_KO_adjp <- p.adjust(bindetect$WT_KO_pvalue)
bindetect$WT_KO_highlighted[bindetect$WT_KO_change > 0 & bindetect$WT_KO_highlighted == "TRUE"] <- "WT_High"
bindetect$WT_KO_highlighted[bindetect$WT_KO_change < 0 & bindetect$WT_KO_highlighted == "TRUE"] <- "KO_High"

tmp <- bindetect[bindetect$WT_KO_highlighted != "FALSE",]
tmp <- tmp[order(abs(tmp$WT_KO_change), decreasing=T),]
tmp <- tmp[!duplicated(tmp$name),]

bindetect$WT_KO_highlighted[!(bindetect$output_prefix %in% tmp$output_prefix)] <- "FALSE"
bindetect <- bindetect[order(bindetect$WT_KO_highlighted),]

bindetect$WT_KO_highlighted[-log10(bindetect$WT_KO_adjp) < 26] <- "FALSE"
bindetect$WT_KO_highlighted[bindetect$name == "FOXP3"] <- "FALSE"

pdf("./../Results/BINDetect/VolcanoPlot.pdf", height=10, width=10)
ggplot(bindetect, aes(x=WT_KO_change, y=-log10(WT_KO_pvalue), color=WT_KO_highlighted)) +
  geom_point() +
  scale_color_manual(values=c("Red", "Blue"), limits = c('WT_High', 'KO_High')) +
	ggrepel::geom_text_repel(aes(x = WT_KO_change, y = -log10(WT_KO_pvalue), label = ifelse(WT_KO_highlighted != "FALSE", as.character(name),"")), max.overlaps=Inf, force=20, segment.color="grey") +
  guides(color=guide_legend(title="")) +
  theme_bw() +
  theme(legend.position = c(0.1, 0.1))
dev.off()

pdf("./../Results/BINDetect/VolcanoPlot_nofont.pdf", height=10, width=10)
ggplot(bindetect, aes(x=WT_KO_change, y=-log10(WT_KO_pvalue), color=WT_KO_highlighted)) +
  geom_point() +
  scale_color_manual(values=c("Red", "Blue"), limits = c('WT_High', 'KO_High')) +
	ggrepel::geom_text_repel(aes(x = WT_KO_change, y = -log10(WT_KO_pvalue), label = ifelse(WT_KO_highlighted != "FALSE", " ",""), size=0), max.overlaps=Inf, force=20, segment.color="grey") +
  guides(color=guide_legend(title="")) +
  theme_bw() +
  theme(legend.position = c(0.1, 0.1)) +
  theme(text = element_text(size=0))
dev.off()

pdf("./../Results/BINDetect/VolcanoPlot_nofont2.pdf", height=10, width=10)
ggplot(bindetect, aes(x=WT_KO_change, y=-log10(WT_KO_pvalue), color=WT_KO_highlighted)) +
  geom_point() +
  scale_color_manual(values=c("Red", "Blue"), limits = c('WT_High', 'KO_High')) +
	ggrepel::geom_text_repel(aes(x = WT_KO_change, y = -log10(WT_KO_pvalue), label = ifelse(WT_KO_highlighted != "FALSE", as.character(name),"")), max.overlaps=Inf, force=20, segment.color="grey") +
  guides(color=guide_legend(title="")) +
  theme_bw() +
  theme(legend.position = c(0.1, 0.1)) +
  theme(text = element_text(size=0))
dev.off()

pdf("./../Results/BINDetect/VolcanoPlot_nofont3.pdf", height=10, width=10)
for(sz in seq(1, 10, 1)){
  print(ggplot(bindetect, aes(x=WT_KO_change, y=-log10(WT_KO_pvalue), color=WT_KO_highlighted)) +
    geom_point() +
    scale_color_manual(values=c("Red", "Blue"), limits = c('WT_High', 'KO_High')) +
    ggrepel::geom_text_repel(aes(x = WT_KO_change, y = -log10(WT_KO_pvalue), label = ifelse(WT_KO_highlighted != "FALSE", as.character(name),"")), size=sz, max.overlaps=Inf, force=20, segment.color="grey") +
    guides(color=guide_legend(title="")) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(text = element_text(size=0)))
}
dev.off()


bindetect.dist <- read.delim("BINDetect_output/bindetect_distances.txt")

fasta <- read.fasta("mm10.fa.gz")

fasta$chrX


promoter <- "GCTTCAGATCCCTTCTTCTGTTCAACCCAGCGATCCTCCAACGTCTCACAAACAC"

i <- 1
cand <- 1:length(fasta$chrX)
for(chr in str_split(cns3, "")[[1]]){
  chr <- tolower(chr)

  cand <- cand[fasta$chrX[cand + i - 1] == chr]

  i <- i + 1
}

promoter <- "GCTTCAGATCCCTTCTTCTGTTCAACCCAGCGATCCTCCAACGTCTCACAAACACAATGCTGTCTCTACCTGCCTCGGGATGCCTTTGTGATTTGACTTATTTTCCCTCAGTTTTTTTTTTCTGACTCTACACACTTTTGTTTAAGAAATTGTGGTTTCTCATGAGCCCTGTTATCTCATTGATACCTTTTACCTCTGTGGTGAGGGGAAGAAATCATATTTTCAGATGACTTGTAAAGGGCAAAGAAAAAACCCAAAATTTCAAAATTTCCGTTTAAGTCTCATAAGAAAAGAATAAACAAAGTAAGAGAGCAAAGAAAAAAAAACTACAAGAACCCCCCCCCCACCCTGCAATTATCAGCACACACACTCATCAAAAAAAAATTGGATTATTAGAAGAGCGAGGTCTGCGGCTTCCAC"
# 7579204 ~ 7579623


cns1 <- "TAGATTACTCTTTTCTTGTGGGGCTTCTGTGTATGGTTTTGTGTTTTAAGTCTTTTGCACTTGAAAATGAGATAACTGTTCACCCCATGTTGGCTTCCAGTCTCCTTTATGGCTTCATTTTTTCCATTTACTGCAGAGGTCAAAAGTGTGGGTATGGGAGCCAGACTGTCTGGAACAACCTAGCCTCAACTCAAGTCATCTGTGTGAATTTTACCCAGGCTCTTAACCTCTCTGTACCTCCATTTCCTCGTATGTACTGTGATGATTATAACAGTACCTACCTCAGAGGATCTTTCTGAGGATTATTTTTATTAATGATGGTAGGTGCTCAGCACAAGGCC"
# 7581695 ~ 7582035

cns2 <- "TGGGTTTTGCATGGTAGCCAGATGGACGTCACCTACCACATCCGCTAGCACCCACATCACCCTACCTGGGCCTATCCGGCTACAGGATAGACTAGCCACTTCTCGGAACGAAACCTGTGGGGTAGATTATCTGCCCCCTTCTCTTCCTCCTTGTTGCCGATGAAGCCCAATGCATCCGGCCGCCATGACGTCAATGGCAGAAAAATCTGGCCAAGTTCAGGTTGTGACAACAGGGCCCAGATGTAGACCCCGATAGGAAAACATATTCTATGTCCCAGAAACAACCTCCATACAGCTTCTAAGAAACAGTCAAACAGGAACGCCCCAACAGACAGTGCAGGAAGCTGGCTGGCCAGCCCAGCCCTCCAGGTCCCTAGTACCACTAGACAGACCATATCCAATTCAGGTCCTCTTTCTGAGAATGTA"
# 7583960 ~ 7584385

cns3 <- "GTGAGGCCCGGGGCCCAGAATGGGGTAAGCAGGGTGGGGTACTTGGGCCTATAGGTGTCGACCTTTACTGTGGCATGTGGCGGGGGGGGGGGGGGGGGCTGGGGCACAGGAAGTGGTTTATGGGTCCCAGGCAAGTCTGACTTATGCAGATATTGCAGGGCCAAGAAAATCCCCACTCTCCAGGCTTCAGAGATTCAAGGCTTTCCCCACCCCTCCCAATCCTCATCCCGATAG"
# 7586562 ~ 7586795
