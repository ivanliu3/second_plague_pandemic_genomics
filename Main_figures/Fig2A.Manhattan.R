library(dplyr)
library(data.table)
library(qqman)
library(stringr)
library(qqman)
library(ggplot2)
library(ggbreak)
library(patchwork)
source("/maps/projects/ilab/people/pls394/plague/supplementaryM_scripts/calculate_lambda.R")

### A
### read pre vs post
### meta all sites
meta = fread("/maps/projects/ilab/people/pls394/plague/Association/GEMMA/output_2026April/meta_analysis_all/GEMMA_preVSpost_Lund-Trondheim.PC71.InverseVariancetbl",header=T)
dim(meta)


### read filtered single GWAS
### read filtered single GWAS                                                                                                                                                                                                                                                                                                                  
lund = fread("/maps/projects/ilab/people/pls394/plague/Association/GEMMA/output_2026April/filter/Lund.pre_vs_post.mm.7_pcs.maf005.rmMULTI-INDEL.rmLOWINFO095.rmLOWGP.rmQ.rmHWE1e-6.rmwithin5bp.tsv",header=T)
trond = fread("/maps/projects/ilab/people/pls394/plague/Association/GEMMA/output_2026April/filter/Trondheim.pre_vs_post.mm.7_pcs.maf005.rmMULTI-INDEL.rmLOWINFO095.rmLOWGP.rmQ.rmHWE1e-6.rmwithin5bp.tsv",header=T)
dim(lund)
dim(trond)

union_rs <- unique(c(lund$rs, trond$rs))
length(union_rs)


### get rid of close  within 5bp                                                                                                                                                                                                                                                                                                               
meta2.filter4= meta %>% filter(MarkerName %in% union_rs)
dim(meta2.filter4)

tmp <- strsplit(meta2.filter4$MarkerName, ":", fixed = TRUE)
meta2.filter4$chr <- as.integer(sub("^chr", "", sapply(tmp, `[`, 1)))
meta2.filter4$ps  <- as.integer(sapply(tmp, `[`, 2))
meta2.filter4 <- meta2.filter4[order(meta2.filter4$chr, meta2.filter4$ps), ]
meta2.filter4 = meta2.filter4 %>% arrange(chr,ps)
dim(meta2.filter4)

t_manhattan <- meta2.filter4[, .(CHR = as.numeric(chr),
                     BP = ps,
                     P = as.numeric(`P-value`),
                     SNP = MarkerName)]


t_plot <- t_manhattan %>%
  mutate(
    CHR = as.numeric(CHR),
    BP = as.numeric(BP),
    P = as.numeric(P),
    logp = -log10(P)
  ) %>%
    arrange(CHR, BP)

## Calculate cumulative chromosome positions
chr_df <- t_plot %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(
    chr_start = lag(cumsum(chr_len), default = 0)
  )

t_plot <- t_plot %>%
  left_join(chr_df, by = "CHR") %>%
  mutate(
    BPcum = BP + chr_start,
    CHR_factor = as.factor(CHR)
  )

## X-axis chromosome label positions
axis_df <- t_plot %>%
  group_by(CHR) %>%
  summarise(
    center = mean(range(BPcum)),
    .groups = "drop"
  )

##
##2. Plot
##-----------------------------
p1 <- ggplot(t_plot, aes(x = BPcum, y = logp, color = CHR_factor)) +
  geom_point(size = 1.3, alpha = 0.8) +
  geom_hline(
    yintercept = -log10(1.1e-8),
    linetype = "solid",
    linewidth = 0.5,
    color='red'
  ) +
    geom_hline(
    yintercept = -log10(1.2e-7),
    linetype = "dashed",
    linewidth = 0.5,
    color='darkgreen'
  ) + geom_hline(
    yintercept = -log10(5.8e-7),
    linetype = "dashed",
    linewidth = 0.5,
    color='purple3'
  ) +
 
  scale_color_manual(
    values = rep(c("blue4", "orange3"), length.out = length(unique(t_plot$CHR)))
  ) +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR,
    expand = c(0.01, 0.01)
  ) +
  scale_y_break(
    c(10, 40),
    scales = 0.2, ## relative height of upper to lower
    space = 0.15,
    symbol = "slash"
  ) +
  scale_y_continuous(
    limits = c(0, 43),
    breaks = c(0, 2, 4, 6, 8, 10, 40, 41, 42, 43)
  ) +
  labs(
    title = "Pre-BD (N=313) versus Post-BD (N=105)",
    ## subtitle = "Meta: Nsite=5,851,478\nLund: Nsite=5,637,038\nTrondheim: Nsite=5,633,253",
    subtitle = "Meta: Nsite=5,851,478",
    x = "Chromosome",
    y = expression(-log[10](P))
  ) +
    theme_classic(base_size = 14) +
    theme(
        legend.position = "none",
        
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),

        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y.right = element_blank(),

        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 19, hjust = 0)

    )

## p1


### B
###read pre vs present-day
### meta all sites
meta = fread("/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/meta_analysis_all/GEMMA_preVSmodern_Lund-Trondheim.PC7.all1.InverseVariancetbl",header=T)
dim(meta)


### read filtered single GWAS
lund = fread("/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/meta_analysis/Lund-danish_skane.pre_vs_present.mm.7_pcs.maf005.rmMULTI-INDEL.rmLOWINFO095.rmLOWGP.rmQ.rmHWE1e-6.rmwithin5bp.forMETAL.tsv",header=T)
trond = fread("/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/meta_analysis/Trondheim-norwegian.pre_vs_present.mm.7_pcs.maf005.rmMULTI-INDEL.rmLOWINFO095.rmLOWGP.rmQ.rmHWE1e-6.rmwithin5bp.forMETAL.tsv",header=T)
dim(lund)
dim(trond)

union_rs <- unique(c(lund$rs, trond$rs))
length(union_rs)




meta2.filter4= meta %>% filter(MarkerName %in% union_rs)
dim(meta2.filter4)

tmp <- strsplit(meta2.filter4$MarkerName, ":", fixed = TRUE)
meta2.filter4$chr <- as.integer(sub("^chr", "", sapply(tmp, `[`, 1)))
meta2.filter4$ps  <- as.integer(sapply(tmp, `[`, 2))
meta2.filter4 <- meta2.filter4[order(meta2.filter4$chr, meta2.filter4$ps), ]
meta2.filter4 = meta2.filter4 %>% arrange(chr,ps)
dim(meta2.filter4)

t_manhattan <- meta2.filter4[, .(CHR = as.numeric(chr),
                     BP = ps,
                     P = as.numeric(`P-value`),
                     SNP = MarkerName)]

t_plot <- t_manhattan %>%
  mutate(
    CHR = as.numeric(CHR),
    BP = as.numeric(BP),
    P = as.numeric(P),
    logp = -log10(P)
  ) %>%
    arrange(CHR, BP)

## Calculate cumulative chromosome positions
chr_df <- t_plot %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(
    chr_start = lag(cumsum(chr_len), default = 0)
  )

t_plot <- t_plot %>%
  left_join(chr_df, by = "CHR") %>%
  mutate(
    BPcum = BP + chr_start,
    CHR_factor = as.factor(CHR)
  )

## X-axis chromosome label positions
axis_df <- t_plot %>%
  group_by(CHR) %>%
  summarise(
    center = mean(range(BPcum)),
    .groups = "drop"
  )

##
##2. Plot
##-----------------------------
p2 <- ggplot(t_plot, aes(x = BPcum, y = logp, color = CHR_factor)) +
  geom_point(size = 1.3, alpha = 0.8) +
    geom_hline(
        yintercept = -log10(1.1e-8),
        linetype = "solid",
        linewidth = 0.5,
        color='red'
    ) +
    geom_hline(
        yintercept = -log10(1.2e-7),
        linetype = "dashed",
        linewidth = 0.5,
        color='darkgreen'
    ) + geom_hline(
            yintercept = -log10(5.8e-7),
            linetype = "dashed",
            linewidth = 0.5,
            color='purple3'
  ) +
  scale_color_manual(
    values = rep(c("blue4", "orange3"), length.out = length(unique(t_plot$CHR)))
  ) +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR,
    expand = c(0.01, 0.01)
  ) +
  scale_y_break(
    c(10, 40),
    scales = 0.2, ## relative height of upper to lower
    space = 0.15,
    symbol = "slash"
  ) +
  scale_y_continuous(
    limits = c(0, 43),
    breaks = c(0, 2, 4, 6, 8, 10, 40, 41, 42, 43)
  ) +
    labs(
        title = "Pre-BD (N=313) versus Present-day (N=1562)",
        ## subtitle = "Meta: Nsite=5,775,573\nLund: Nsite=5,652,669\nTrondheim: Nsite=5,659,336",
        subtitle = "Meta: Nsite=5,775,573",
        x = "Chromosome",
        y = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
        legend.position = "none",
        
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 14),
        axis.ticks.y = element_line(),
        axis.line.y = element_line(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_blank(),

        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.y.right = element_blank(),
        
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 19, hjust = 0)
    )
## p2









### Combined side-by-side plot
### Use the ggplot objects above so both panels use exactly the same y-axis scale.
### Only p1 keeps the visible y-axis; p2 has its y-axis removed.
combined_plot <- p1 + p2 +  plot_layout(ncol = 2, widths = c(1, 1))

## combined_plot

ggsave(
  filename = "/projects/ilab/people/pls394/plague/important_scripts/main_fig/Manhattan.GEMMA_side_by_side.single_shared_yaxis.2026Jun08.png",
  plot = combined_plot,
  height = 7,
  width = 25,
  dpi = 300,
  bg = "white"
)

