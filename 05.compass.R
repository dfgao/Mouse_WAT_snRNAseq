library(cowplot)
library(scales)
library(ggplot2)
library(dplyr)
set.seed(1234)

# adipo ----
adipo.compass <- read.csv('compass.adipo.csv',header = T,row.names = 1)

adipo.compass$change <- factor(ifelse(adipo.compass$wilcox_pval < 0.05, 
                                      ifelse(adipo.compass$cohens_d > 0, 'HP up','LP up'), 'No sig'))
table(adipo.compass$change)
unique(adipo.compass$subsystem)

ggplot(adipo.compass, aes(x = cohens_d, y = 1)) +
  geom_point(aes(color = change), position = position_jitter(height = 0.15, width = 0), size = 1, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black",size = 1) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(x = "Cohen's d", title = "Cohen's d distribution") +
  xlim(-2,2) + 
  scale_color_manual(values = c('#F28E2B','#4F81BD','#BFBFBF')) +
  theme_minimal(base_size = 12) +
  theme(axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5) )

adipo.compass[adipo.compass$subsystem == 'Valine, leucine, and isoleucine metabolism',]$subsystem <- 'BCAA cat.'

select.sub <- c('BCAA cat.','Triacylglycerol synthesis','Fatty acid synthesis','Fatty acid oxidation')
plot_df <- adipo.compass %>% filter(subsystem %in% select.sub) %>%mutate(subsystem = factor(subsystem, levels = rev(select.sub)))

ggplot(plot_df, aes(x = cohens_d, y = subsystem, color = change)) +
  # geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(position = position_jitter(height = 0.2, width = 0), size = 1, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black",size = 1) +
  labs(x = "Cohen's d", y = NULL, title = "Cohen's d by subsystem") +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c('#F28E2B','#4F81BD','#BFBFBF')) +
  xlim(-2,2) + 
  theme(
    axis.text.y = element_text(face = "bold",size = 12),
    panel.grid.major.y = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5) )

# socre 
library(tidyverse)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(KEGGREST)

nature.pathway <- rio::import('../nature.adipo.pathway.xlsx')
uniques <- lapply(nature.pathway, unique)
maxlen <- max(sapply(uniques, length))
uniq_df <- as.data.frame(lapply(uniques, 
                                function(x) c(x, rep(NA, maxlen - length(x)))), 
                         stringsAsFactors = F)

DNL.pathway <- keggGet('mmu00061')
DNL.pathway <- sapply(seq(2,38,2), function(x){
  gns <- unlist(strsplit(DNL.pathway[[1]][['GENE']][x],";"))[1]
})
uniques$DNL <- DNL.pathway

DefaultAssay(adipo.seu) <- "SCT"  
adipo.seu <- AddModuleScore(
  object  = adipo.seu,
  features = uniques,
  name     = "nat.score",
  assay    = "SCT",
  search   = F,        
  ctrl     = 100       
)
colnames(adipo.seu@meta.data) <- c(colnames(adipo.seu@meta.data)[-(44:55)], paste0(names(uniques),'_score'))

df <- adipo.seu@meta.data[,c(22,44:55)]
df_long <- df %>% pivot_longer(cols = -condition, names_to = "variable", values_to = "score")
df_long$condition <- factor(df_long$condition, levels = c("LP","HP"))

pairwise_res <- df_long %>%
  group_by(variable) %>%
  summarise(
    p.value = tryCatch(wilcox.test(score ~ condition)$p.value, error = function(e) NA_real_),
    ymax = max(score, na.rm = TRUE),
    ymin = min(score, na.rm = TRUE),
    .groups = "drop"
  )

pairwise_res <- pairwise_res %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    p.label = ifelse(is.na(p.adj), "NA",
                     ifelse(p.adj < 0.001, "<0.001", signif(p.adj, 3)))
  )

pairwise_res <- pairwise_res %>%
  mutate(
    y.position = ymax + 0.08 * (ymax - ymin + 1e-6),  
    x.position = 1.5  
  )

mycols <- c("LP" = "#0B996F", "HP" = "#D6570D")

p <- ggplot(df_long, aes(x = condition, y = score, fill = condition)) +
  geom_violin(trim = FALSE, width = 0.7, alpha = 0.85) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.7) +
  # geom_jitter(width = 0.12, size = 0.1, aes(color = condition), alpha = 0.6) +
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        axis.title.x = element_blank()) +
  labs(y = "Metabolic pathway scores")

p <- p + geom_text(
  data = pairwise_res,
  aes(x = x.position, y = y.position, label = p.label),
  inherit.aes = FALSE,
  size = 3
)
p

### ridge
library(ggridges)

df_long$condition <- factor(df_long$condition, levels = c('LP','HP'))
ggplot(df_long, aes(x = score, y = variable, fill = condition)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.5) +
  scale_fill_manual(values = mycols) +
  # ggtitle("SASP addscore") + xlab('SASP score') + ylab("") +
  # xlim(0.0,0.15) +
  # facet_wrap(~ variable, nrow = 3) + 
  
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6),
        legend.position ="none")

p1 <- ggplot(df, aes(x = Adap.therm_score, y = condition, fill = condition)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.6) +
  scale_fill_manual(values = mycols) +
  ggtitle("Adap.them. addscore") + xlab('Adap.them. score') + ylab("") +
  # xlim(0.0,0.15) +
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6),
        legend.position ="none")
p2 <- ggplot(df, aes(x = BCAA_score, y = condition, fill = condition)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.6) +
  scale_fill_manual(values = mycols) +
  ggtitle("BCAA addscore") + xlab('BCAA score') + ylab("") +
  # xlim(0.0,0.15) +
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6),
        legend.position ="none")
p3 <- ggplot(df, aes(x = BETA.OX_score, y = condition, fill = condition)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.6) +
  scale_fill_manual(values = mycols) +
  ggtitle("BETA.OX addscore") + xlab('BETA.OX score') + ylab("") +
  # xlim(0.0,0.15) +
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6),
        legend.position ="none")
p4 <- ggplot(df, aes(x = Lipid.dp_score, y = condition, fill = condition)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.6) +
  scale_fill_manual(values = mycols) +
  ggtitle("Lipid.dp. addscore") + xlab('Lipid.dp. score') + ylab("") +
  # xlim(0.0,0.15) +
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6),
        legend.position ="none")
p5 <- ggplot(df, aes(x = Lipolysis_score, y = condition, fill = condition)) +
  geom_density_ridges(quantile_lines = TRUE, quantiles = 2,alpha = 0.6) +
  scale_fill_manual(values = mycols) +
  ggtitle("Lipolysis addscore") + xlab('Lipolysis score') + ylab("") +
  # xlim(0.0,0.15) +
  theme_test(base_size = 15) +
  theme(panel.border=element_rect(linewidth = 1, color = "black"),
        strip.background = element_rect(linewidth = 1, fill = "white"),
        strip.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16, colour = "black"),
        axis.line = element_line(size = 0.6),
        legend.position ="none")
patchwork::wrap_plots(c(p1,p2,p3,p4,p5),nrow = 3)
(p1 + p2 + p3 + p4 + p5)

## GSVA
library(GSVA)
expr_mat <- GetAssayData(adipo.seu, slot="data",assay = 'SCT')
ssgsea_res <- gsva(as.matrix(expr_mat),
                   uniques,
                   method="ssgsea",
                   ssgsea.norm=TRUE,
                   parallel.sz=40)

meta_add <- t(ssgsea_res) %>% as.data.frame()
colnames(meta_add) <- c(paste0(colnames(meta_add), '_GSVA'))
adipo.seu <- AddMetaData(adipo.seu, metadata = meta_add)

df <- adipo.seu@meta.data[,c(22,56:67)]
df_long <- df %>% pivot_longer(cols = -condition, names_to = "variable", values_to = "score")
df_long$condition <- factor(df_long$condition, levels = c("LP","HP"))

pairwise_res <- df_long %>%
  group_by(variable) %>%
  summarise(
    p.value = tryCatch(wilcox.test(score ~ condition)$p.value, error = function(e) NA_real_),
    ymax = max(score, na.rm = TRUE),
    ymin = min(score, na.rm = TRUE),
    .groups = "drop"
  )

pairwise_res <- pairwise_res %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    p.label = ifelse(is.na(p.adj), "NA",
                     ifelse(p.adj < 0.001, "<0.001", signif(p.adj, 3)))
  )

pairwise_res <- pairwise_res %>%
  mutate(
    y.position = ymax + 0.08 * (ymax - ymin + 1e-6),  
    x.position = 1.5  
  )

mycols <- c("LP" = "#0B996F", "HP" = "#D6570D")

p <- ggplot(df_long, aes(x = condition, y = score, fill = condition)) +
  geom_violin(trim = FALSE, width = 0.7, alpha = 0.85) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.7) +
  # geom_jitter(width = 0.12, size = 0.1, aes(color = condition), alpha = 0.6) +
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        axis.title.x = element_blank()) +
  labs(y = "Metabolic pathway scores")

p <- p + geom_text(
  data = pairwise_res,
  aes(x = x.position, y = y.position, label = p.label),
  inherit.aes = FALSE,
  size = 3
)
p

# heatmap
## addscore
library(RColorBrewer)

df <- adipo.seu@meta.data[,c(27,56:67)]
score_cols <- setdiff(colnames(df), "SCT_snn_res.0.2")

df_long <- df %>%
  pivot_longer(cols = all_of(score_cols), names_to = "pathway", values_to = "score")
df_long$SCT_snn_res.0.2 <- factor(paste0('AD_',df_long$SCT_snn_res.0.2), levels = c(paste0('AD_',seq(0,4))))
colnames(df_long) <- c('cluster','pathway','score')

mean_df <- df_long %>%
  group_by(cluster, pathway) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

scaled_df <- mean_df %>%
  group_by(pathway) %>%
  mutate(scaled = (mean_score - min(mean_score, na.rm = TRUE)) /
           (max(mean_score, na.rm = TRUE) - min(mean_score, na.rm = TRUE) + 1e-9)) %>%
  ungroup()

cluster_levels <- sort(unique(scaled_df$cluster))
scaled_df$cluster <- factor(scaled_df$cluster, levels = cluster_levels)

pathway_order <- scaled_df %>%
  group_by(pathway) %>%
  summarise(overall = mean(mean_score, na.rm = TRUE)) %>%
  arrange(desc(overall)) %>%
  pull(pathway)

scaled_df$pathway <- factor(scaled_df$pathway, levels = rev(pathway_order))

colors <- colorRampPalette(c("#e6f5d0","#a1dab4","#41b6c4","#225ea8","#0c2c84"))(100)

ggplot(scaled_df, aes(x = cluster, y = pathway, fill = scaled)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = colors, limits = c(0,1), na.value = "grey90",
                       name = "Scaled score") +
  theme_minimal(base_size = 12) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 10),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm")
  ) +
  coord_fixed(ratio = 0.7) +
  scale_fill_gradientn(colors = colors, limits = c(0,1), na.value = "grey90", name = "Scaled score",
                       breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) + 
  guides(fill = guide_colorbar(
    title.position = "top",   
    title.hjust = 0.5,        
    barwidth = unit(4, "cm"),  
    barheight = unit(0.35, "cm")
  )) 


# macrop ----
macrop.compass <- read.csv('compass.macrop.csv',header = T,row.names = 1)

macrop.compass$change <- factor(ifelse(macrop.compass$wilcox_pval < 0.05, 
                                      ifelse(macrop.compass$cohens_d > 0, 'HP up','LP up'), 'No sig'))
table(macrop.compass$change)
unique(macrop.compass$subsystem)

ggplot(macrop.compass, aes(x = cohens_d, y = 1)) +
  geom_point(aes(color = change), position = position_jitter(height = 0.15, width = 0), size = 1, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black",size = 1) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(x = "Cohen's d", title = "Cohen's d distribution") +
  xlim(-2,2) + 
  scale_color_manual(values = c('#F28E2B','#4F81BD','#BFBFBF')) +
  theme_minimal(base_size = 12) +
  theme(axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5) )

macrop.compass[macrop.compass$subsystem == 'Valine, leucine, and isoleucine metabolism',]$subsystem <- 'BCAA cat.'

select.sub <- c('Pentose phosphate pathway','Oxidative phosphorylation','Glycolysis/gluconeogenesis','Fatty acid synthesis','Fatty acid oxidation')
plot_df <- macrop.compass %>% filter(subsystem %in% select.sub) %>%mutate(subsystem = factor(subsystem, levels = rev(select.sub)))

ggplot(plot_df, aes(x = cohens_d, y = subsystem, color = change)) +
  # geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(position = position_jitter(height = 0.2, width = 0), size = 1, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black",size = 1) +
  labs(x = "Cohen's d", y = NULL, title = "Cohen's d by subsystem") +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c('#F28E2B','#4F81BD','#BFBFBF')) +
  xlim(-2,2) + 
  theme(
    axis.text.y = element_text(face = "bold",size = 12),
    panel.grid.major.y = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5) )


# key pathway ----
## adipo
load('../../06.uterus/03.analysis/mus.uterus/GO_DATA.RData')
findGO(pattern = 'insulin')
ins.rep.act <- GO_DATA_BP$PATHID2EXTID$`GO:0030070`




