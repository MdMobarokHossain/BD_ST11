############################################################
## Refined Pangenome Figure for Publication
## Reviewer-friendly version
## Replaced misleading Venn diagram with
## biologically meaningful gene-frequency distribution
############################################################

library(tidyverse)
library(cowplot)
library(ggrepel)

setwd("C:/New_folder/Others/KLEB-2024")

############################################################
## 1. Read pangenome matrix
############################################################

rtab <- read_tsv("gene_presence_absence.Rtab")

# First column = gene ID
gene_mat <- as.matrix(rtab[, -1])

# Convert to numeric if needed
gene_mat <- apply(gene_mat, 2, as.numeric)

rownames(gene_mat) <- rtab[[1]]

n_genomes <- ncol(gene_mat)
n_genes   <- nrow(gene_mat)

############################################################
## 2. Define gene categories
############################################################

gene_freq <- rowSums(gene_mat > 0)
gene_prop <- gene_freq / n_genomes

gene_class <- case_when(
  gene_prop == 1                        ~ "Core",
  gene_prop >= 0.95 & gene_prop < 1    ~ "Soft core",
  gene_prop >= 0.15 & gene_prop < 0.95 ~ "Shell",
  TRUE                                  ~ "Cloud"
)

table(gene_class)

############################################################
## 3. Pan-genome rarefaction
############################################################

set.seed(123)

n_perm <- 100L

pan_mat <- matrix(
  NA_integer_,
  nrow = n_perm,
  ncol = n_genomes
)

for (p in seq_len(n_perm)) {
  
  order_idx <- sample(seq_len(n_genomes))
  
  seen <- rep(FALSE, n_genes)
  
  for (i in seq_len(n_genomes)) {
    
    g <- order_idx[i]
    
    seen <- seen | (gene_mat[, g] > 0)
    
    pan_mat[p, i] <- sum(seen)
  }
}

pan_mean <- colMeans(pan_mat)

pan_sd <- apply(pan_mat, 2, sd)

pan_df <- tibble(
  N        = 1:n_genomes,
  pan_mean = pan_mean,
  pan_sd   = pan_sd,
  pan_lo   = pan_mean - 1.96 * pan_sd / sqrt(n_perm),
  pan_hi   = pan_mean + 1.96 * pan_sd / sqrt(n_perm)
)

############################################################
## 4. Core-genome rarefaction
############################################################

core_mat <- matrix(
  NA_integer_,
  nrow = n_perm,
  ncol = n_genomes
)

for (p in seq_len(n_perm)) {
  
  order_idx <- sample(seq_len(n_genomes))
  
  present <- rep(TRUE, n_genes)
  
  for (i in seq_len(n_genomes)) {
    
    g <- order_idx[i]
    
    present <- present & (gene_mat[, g] > 0)
    
    core_mat[p, i] <- sum(present)
  }
}

core_mean <- colMeans(core_mat)

core_sd <- apply(core_mat, 2, sd)

core_df <- tibble(
  N         = 1:n_genomes,
  core_mean = core_mean,
  core_sd   = core_sd,
  core_lo   = core_mean - 1.96 * core_sd / sqrt(n_perm),
  core_hi   = core_mean + 1.96 * core_sd / sqrt(n_perm)
)

############################################################
## 5. PCA of gene presence/absence
############################################################

pca <- prcomp(
  t(gene_mat),
  center = TRUE,
  scale. = FALSE
)

pca_df <- as.data.frame(pca$x)

pca_df$genome <- colnames(gene_mat)

# Highlight study isolates
pca_df$highlight <- ifelse(
  grepl("^KPN", pca_df$genome),
  "KPN",
  "Other"
)

############################################################
## Remove PCA outliers
############################################################

pca_core <- pca_df[, c("PC1", "PC2")]

d <- mahalanobis(
  pca_core,
  colMeans(pca_core),
  cov(pca_core)
)

cutoff <- qchisq(0.975, df = 2)

pca_df$outlier <- d > cutoff

pca_clean <- subset(
  pca_df,
  outlier == FALSE
)

############################################################
## 6. Lancet-style theme
############################################################

theme_lancet <- function() {
  
  theme_classic(base_size = 12) +
    
    theme(
      axis.title = element_text(
        face = "bold",
        size = 12
      ),
      
      axis.text = element_text(
        color = "black",
        size = 10
      ),
      
      plot.title = element_text(
        face = "bold",
        hjust = 0,
        size = 13
      ),
      
      legend.title = element_text(
        face = "bold"
      ),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.6
      ),
      
      axis.line = element_line(
        colour = "black"
      ),
      
      plot.margin = margin(
        5, 5, 5, 5
      )
    )
}

############################################################
## 7A. Pan-genome plot
############################################################

p_pan <- ggplot(
  pan_df,
  aes(x = N, y = pan_mean)
) +
  
  geom_ribbon(
    aes(
      ymin = pan_lo,
      ymax = pan_hi
    ),
    fill = "#4575b4",
    alpha = 0.2
  ) +
  
  geom_line(
    linewidth = 1,
    color = "#4575b4"
  ) +
  
  labs(
    x = "Number of genomes",
    y = "Pan-genome size",
    title = "A. Pan-genome rarefaction"
  ) +
  
  theme_lancet()

############################################################
## 7B. Core-genome plot
############################################################

p_core <- ggplot(
  core_df,
  aes(x = N, y = core_mean)
) +
  
  geom_ribbon(
    aes(
      ymin = core_lo,
      ymax = core_hi
    ),
    fill = "#d73027",
    alpha = 0.2
  ) +
  
  geom_line(
    linewidth = 1,
    color = "#d73027"
  ) +
  
  labs(
    x = "Number of genomes",
    y = "Core-genome size",
    title = "B. Core-genome rarefaction"
  ) +
  
  theme_lancet()

############################################################
## 7C. PCA plot
############################################################

p_pca <- ggplot(
  pca_clean,
  aes(PC1, PC2)
) +
  
  geom_point(
    aes(color = highlight),
    size = 2.2,
    alpha = 0.85
  ) +
  
  scale_color_manual(
    values = c(
      "KPN" = "#1f78b4",
      "Other" = "grey70"
    )
  ) +
  
  stat_ellipse(
    level = 0.95,
    linewidth = 0.8,
    color = "black"
  ) +
  
  geom_text_repel(
    data = subset(
      pca_clean,
      highlight == "KPN"
    ),
    
    aes(label = genome),
    
    size = 3,
    fontface = "bold",
    color = "#1f78b4",
    
    box.padding = 0.5,
    point.padding = 0.3,
    
    segment.color = "grey40",
    segment.size = 0.4,
    
    max.overlaps = Inf
  ) +
  
  labs(
    x = paste0(
      "PC1 (",
      round(
        100 * summary(pca)$importance[2, 1],
        1
      ),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(
        100 * summary(pca)$importance[2, 2],
        1
      ),
      "%)"
    ),
    
    title = "C. PCA of gene presence–absence"
  ) +
  
  theme_lancet() +
  
  theme(
    legend.position = "right"
  )

############################################################
## 7D. Gene frequency distribution plot
## Replaces misleading Venn diagram
############################################################

freq_df <- data.frame(
  gene       = rownames(gene_mat),
  frequency  = gene_freq,
  proportion = gene_prop,
  class      = gene_class
)

freq_plot_df <- freq_df %>%
  count(frequency, class)

p_freq <- ggplot(
  freq_plot_df,
  aes(
    x = frequency,
    y = n,
    fill = class
  )
) +
  
  geom_col(
    width = 0.9,
    color = "black",
    linewidth = 0.2
  ) +
  
  scale_fill_manual(
    values = c(
      "Core"      = "#1b9e77",
      "Soft core" = "#66a61e",
      "Shell"     = "#e6ab02",
      "Cloud"     = "#d95f02"
    )
  ) +
  
  labs(
    x = "Number of genomes containing gene",
    y = "Number of genes",
    fill = "Gene category",
    title = "D. Gene frequency distribution"
  ) +
  
  theme_lancet() +
  
  theme(
    legend.position = "right"
  )

############################################################
## 8. Combine panels
############################################################

top_row <- plot_grid(
  p_pan,
  p_core,
  ncol = 2,
  rel_widths = c(1, 1)
)

bottom_row <- plot_grid(
  p_pca,
  p_freq,
  ncol = 2,
  rel_widths = c(1, 1)
)

combined_fig <- plot_grid(
  top_row,
  bottom_row,
  nrow = 2,
  rel_heights = c(1, 1)
)

############################################################
## 9. Save high-resolution figures
############################################################

ggsave(
  filename = "pangenome_refined_figure.tiff",
  plot = combined_fig,
  width = 8,
  height = 6,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "pangenome_refined_figure.pdf",
  plot = combined_fig,
  width = 8,
  height = 6
)

############################################################
## 10. Accessory gene analysis
############################################################

accessory_idx <- gene_class %in% c(
  "Shell",
  "Cloud"
)

accessory_genes <- rownames(gene_mat)[
  accessory_idx
]

############################################################
## Study isolates
############################################################

kpn_cols <- grepl(
  "^KPN",
  colnames(gene_mat)
)

study_mat <- gene_mat[
  ,
  kpn_cols,
  drop = FALSE
]

nonstudy_mat <- gene_mat[
  ,
  !kpn_cols,
  drop = FALSE
]

############################################################
## Accessory genes present in study isolates
############################################################

study_accessory_present <-
  rowSums(
    study_mat[
      accessory_idx,
      ,
      drop = FALSE
    ] > 0
  ) > 0

accessory_in_study <-
  accessory_genes[
    study_accessory_present
  ]

write.table(
  accessory_in_study,
  file = "accessory_genes_in_study_isolates.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

############################################################
## Accessory genes unique to study isolates
############################################################

if (ncol(nonstudy_mat) > 0) {
  
  study_present <-
    rowSums(
      study_mat[
        accessory_idx,
        ,
        drop = FALSE
      ] > 0
    ) > 0
  
  nonstudy_present <-
    rowSums(
      nonstudy_mat[
        accessory_idx,
        ,
        drop = FALSE
      ] > 0
    ) > 0
  
  accessory_unique_study <-
    accessory_genes[
      study_present &
        !nonstudy_present
    ]
  
  write.table(
    accessory_unique_study,
    file =
      "accessory_genes_unique_to_study_isolates.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

############################################################
## Accessory gene counts per isolate
############################################################

study_accessory_counts <-
  colSums(
    study_mat[
      accessory_idx,
      ,
      drop = FALSE
    ] > 0
  )

accessory_summary <- data.frame(
  genome_id = colnames(study_mat),
  n_accessory_genes =
    study_accessory_counts
) %>%
  arrange(
    desc(n_accessory_genes)
  )

write.csv(
  accessory_summary,
  file =
    "study_isolates_accessory_gene_counts.csv",
  row.names = FALSE
)

############################################################
## Summary
############################################################

table(gene_class)
