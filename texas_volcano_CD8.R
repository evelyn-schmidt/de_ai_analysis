# ============================================================
#  Texas-Themed Volcano Plot: CD8 T Cells — ICBdT vs ICB
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)

# ── 1. Add treatment column from orig.ident ──────────────────
#    orig.ident contains "Rep1_ICBdT", "Rep2_ICB", etc.
#    We extract the treatment label (ICBdT or ICB) after the last underscore.
merged$treatment <- sub(".*_(ICBdT|ICB)$", "\\1", merged$orig.ident)

# Quick check
table(merged$treatment)

# ── 2. Subset to CD8 T cells ─────────────────────────────────
#    Uses immgen_singler_main; adjust if you have a refined CD8 label
cd8 <- subset(merged, subset = immgen_singler_main == "T cells" &
                grepl("8", immgen_singler_fine))   # catches T.8* labels

# If the line above grabs too few/many cells, inspect and adjust:
# table(merged$immgen_singler_fine[merged$immgen_singler_main == "T cells"])
# Alternative: use a CD8 cluster number, e.g.:
# cd8 <- subset(merged, subset = seurat_clusters == YOUR_CD8_CLUSTER)

cat("CD8 T cells per treatment:\n")
print(table(cd8$treatment))

# ── 3. Pseudo-bulk DE: ICBdT vs ICB ──────────────────────────
Idents(cd8) <- "treatment"

de_results <- FindMarkers(
  cd8,
  ident.1    = "ICBdT",
  ident.2    = "ICB",
  test.use   = "wilcox",
  min.pct    = 0.1,
  logfc.threshold = 0,   # keep all genes; filter visually
  verbose    = FALSE
)

de_results$gene      <- rownames(de_results)
de_results$neg_log10p <- -log10(de_results$p_val_adj + 1e-300)

# Significance / direction labels for colouring
de_results <- de_results %>%
  mutate(significance = case_when(
    p_val_adj < 0.05 & avg_log2FC >  0.5 ~ "Up in ICBdT",
    p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Up in ICB",
    TRUE                                  ~ "NS"
  ))

# Top genes to label (by adj p-value, within significant hits)
top_up   <- de_results %>% filter(significance == "Up in ICBdT") %>%
              slice_min(p_val_adj, n = 10)
top_down <- de_results %>% filter(significance == "Up in ICB")   %>%
              slice_min(p_val_adj, n = 10)
label_genes <- bind_rows(top_up, top_down)

# ── 4. Texas-themed colour palette ───────────────────────────
#    Burnt orange (UT Austin) + midnight blue + warm sand background
tx_up    <- "#BF5700"   # burnt orange  — up in ICBdT
tx_down  <- "#003087"   # navy blue     — up in ICB
tx_ns    <- "#C4A882"   # sand / taupe  — not significant
tx_bg    <- "#FDF6EC"   # warm cream background
tx_grid  <- "#E8D5B7"   # light tan grid lines
tx_title <- "#2B1A0E"   # dark espresso text

# ── 5. Build the volcano plot ─────────────────────────────────
p <- ggplot(de_results, aes(x = avg_log2FC, y = neg_log10p,
                             colour = significance)) +

  # Reference lines
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = "dashed", linewidth = 0.4, colour = "#9E8060") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", linewidth = 0.4, colour = "#9E8060") +

  # Points
  geom_point(alpha = 0.65, size = 1.6, stroke = 0) +

  # Colour scale — Texas flag-inspired
  scale_colour_manual(
    values = c("Up in ICBdT" = tx_up,
               "Up in ICB"   = tx_down,
               "NS"           = tx_ns),
    name = NULL
  ) +

  # Gene labels
  ggrepel::geom_text_repel(
    data        = label_genes,
    aes(label   = gene),
    size        = 2.8,
    fontface    = "bold.italic",
    colour      = tx_title,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.colour = "#9E8060",
    segment.size   = 0.3,
    max.overlaps   = 20,
    seed           = 42
  ) +

  # Axes & labels
  labs(
    title    = "CD8 T Cells — ICBdT vs ICB",
    subtitle = "Mouse bladder tumour PBMCs  •  Wilcoxon rank-sum, min.pct = 0.1",
    x        = expression(Log[2]~Fold~Change~(ICBdT/ICB)),
    y        = expression(-Log[10]~Adjusted~italic(p)*"-value"),
    caption  = "Dashed lines: |log₂FC| = 0.5,  adj.p = 0.05"
  ) +

  # Texas-themed theme
  theme_minimal(base_size = 13) +
  theme(
    plot.background    = element_rect(fill = tx_bg,  colour = NA),
    panel.background   = element_rect(fill = tx_bg,  colour = NA),
    panel.grid.major   = element_line(colour = tx_grid, linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(colour = "#9E8060", fill = NA,
                                       linewidth = 0.8),

    plot.title         = element_text(colour = tx_title,  face = "bold",
                                       size = 16, family = "serif",
                                       margin = margin(b = 4)),
    plot.subtitle      = element_text(colour = "#6B4C2A", size = 10,
                                       family = "serif"),
    plot.caption       = element_text(colour = "#9E8060", size = 8,
                                       hjust = 0),
    plot.margin        = margin(20, 20, 15, 20),

    axis.title         = element_text(colour = tx_title, face = "bold",
                                       size = 11),
    axis.text          = element_text(colour = "#4A3520"),
    axis.ticks         = element_line(colour = "#9E8060"),

    legend.position    = "top",
    legend.text        = element_text(colour = tx_title, size = 10,
                                       face = "bold"),
    legend.key         = element_rect(fill = tx_bg, colour = NA),
    legend.background  = element_rect(fill = tx_bg, colour = NA)
  ) +

  # Lone-star decorative annotation
  annotate("text", x = Inf, y = Inf, label = "\u2605",
           hjust = 1.3, vjust = 1.5, size = 10,
           colour = "#BF5700", alpha = 0.25) +

  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))

# ── 6. Save ───────────────────────────────────────────────────
ggsave("texas_volcano_CD8_ICBdT_vs_ICB.pdf",
       plot   = p,
       width  = 9,
       height = 7,
       device = cairo_pdf)

ggsave("texas_volcano_CD8_ICBdT_vs_ICB.png",
       plot   = p,
       width  = 9,
       height = 7,
       dpi    = 300)

message("✔  Saved: texas_volcano_CD8_ICBdT_vs_ICB.pdf / .png")
print(p)

# ── 7. Summary table of top hits ──────────────────────────────
cat("\n── Top 10 Up in ICBdT ──\n")
print(top_up  %>% select(gene, avg_log2FC, p_val_adj) %>% arrange(p_val_adj))

cat("\n── Top 10 Up in ICB ──\n")
print(top_down %>% select(gene, avg_log2FC, p_val_adj) %>% arrange(p_val_adj))
