library(ggplot2)
library(dplyr)
library(data.table)
library(qqman)
library(stringr)
library(qqman)
source("/maps/projects/ilab/people/pls394/plague/supplementaryM_scripts/calculate_lambda.R")

filter_colors <- c(
  "Pass" = "grey40",
  "MAF < 0.05" = "#FF7F00",
  "info Score < 0.95" = "#377EB8",
  "deCODE quality concern" = "#4DAF4A",
  "Multiallelic INDEL" = "#984EA3",
  "HWE pvalue < 1e-6" =  "#E41A1C",
  "Variants within 5bp" = "#A65628"
)


### site information ###
info =fread("/maps/projects/ilab/people/pls394/plague/Site_filter/2026Feb/WWW_221011.infoQ.final.v2.tsv",header=T)
dim(info)

hwe = fread("/maps/projects/ilab/people/pls394/plague/HWE_new/output/pcangsd.pre-rm1strelative.default.hwe.2026Feb.tsv",header=T)

lund.beforeNearby = read.table("/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/filter/Lund-danish_skane.pre_vs_present.mm.7_pcs.maf005.rmMULTI-INDEL.rmLOWINFO095.rmLOWGP.rmQ.rmHWE1e-6.modernDenmarkHWE1e-6.tsv",sep="\t",header=T)
lund.afterNearby = fread("/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/filter/Lund-danish_skane.pre_vs_present.mm.7_pcs.maf005.rmMULTI-INDEL.rmLOWINFO095.rmLOWGP.rmQ.rmHWE1e-6.modernDenmarkHWE1e-6.rmwithin5bp.tsv",header=T)

nearby.bad.marker = lund.beforeNearby$rs[! lund.beforeNearby$rs %in% lund.afterNearby$rs]
nearby.df = data.frame(rs = nearby.bad.marker,removal_nearby= 'Y')

#### 2026 Jun 21, as Ida asked
## no filter lund
lund = fread("/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/Lund-danish_skane.99610455.2026-02-19.pre_vs_present.mm.7_pcs.assoc.min.txt.gz",header=T)
dim(lund)

lund$chr <- gsub(":.*", "", lund$rs)
lund$chr = gsub("chr","",lund$chr)
lund = lund %>% filter(chr != 'X')
lund$ps <- as.numeric(str_split_fixed(lund$rs, ":", 3)[, 2])
lund = lund %>% arrange(chr,ps)
dim(lund)

t_manhattan <- lund[, .(CHR = as.numeric(chr), 
                     BP = ps, 
                     P = p_score,
                     SNP = rs,
                     MAF = af)]
dim(t_manhattan)

merged = t_manhattan %>% left_join(info[,c(1,2,9:14)],by=c('SNP'= 'ID'))
dim(merged)
head(merged)

merged = merged %>% left_join(hwe %>% select(ID2,pvalue),by=c('ID2' = 'ID2'))
names(merged)[ncol(merged)] = 'HWEpvalue'


merged = merged %>% left_join(nearby.df,by=c('SNP'='rs'))





## chr1:225574884:SG
lead.chr =1
lead.pos = 225574884
step = 500000
lower.pos = lead.pos -step
upper.pos =lead.pos +step

locus_df <- merged %>%
  filter(
    CHR == lead.chr,
    BP > lower.pos,
    BP < upper.pos
  ) %>%
  mutate(
    info_filter  = infoScore_final < 0.95 |  Prop_lowMAXGP >0.05 ,
    seqq_filter  = str_detect(coalesce(SeqQFlags, ""), "Q"),
    indel_filter = MULTI_INDEL == "Y",
    pval_filter  = HWEpvalue < 1e-6,
    maf_filter   = MAF < 0.05 | MAF >0.95,
    nearby_filter = removal_nearby =='Y',
    
    filter_group = case_when(
      maf_filter   ~ "MAF < 0.05",
      info_filter  ~ "info Score < 0.95",
      seqq_filter  ~ "deCODE quality concern",
      indel_filter ~ "Multiallelic INDEL",
      pval_filter  ~ "HWE pvalue < 1e-6",
      nearby_filter ~ "Variants within 5bp",
      TRUE ~ "Pass"
    ),
    filter_group = factor(filter_group, levels = names(filter_colors)),

    logp = -log10(P),
    dot_size = if_else(maf_filter, 0.6, 1.0)
  )

p= ggplot(locus_df, aes(x = BP, y = logp, color = filter_group)) +
  geom_point(aes(size = dot_size), alpha = 0.8) +
  scale_size_identity(guide = "none") +
  scale_color_manual(values = filter_colors, breaks = names(filter_colors), drop = FALSE) +
  geom_vline(xintercept = lead.pos, linetype = "dashed") +
  labs(
    title = "LocusZoom-style plot",
    subtitle = paste0("chr", lead.chr, ":", lower.pos, "-", upper.pos),
    x = paste0("Chromosome ", lead.chr, " position"),
    y = expression(-log[10](P)),
    color = "Filter"
  ) +
    theme_classic()



pdf('/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/LocusZoom_colorBYfilter/Lund.preVpresent.nofilter.chr1.manhattan.2026Jun27.pdf',h=6,w=10)
plot(p)
dev.off()



## chr2:160798465:SG:1
lead.chr =2
lead.pos = 160798465
step = 500000
lower.pos = lead.pos -step
upper.pos =lead.pos +step

locus_df <- merged %>%
  filter(
    CHR == lead.chr,
    BP > lower.pos,
    BP < upper.pos
  ) %>%
  mutate(
    info_filter  = infoScore_final < 0.95 |  Prop_lowMAXGP >0.05 ,
    seqq_filter  = str_detect(coalesce(SeqQFlags, ""), "Q"),
    indel_filter = MULTI_INDEL == "Y",
    pval_filter  = HWEpvalue < 1e-6,
    maf_filter   = MAF < 0.05 | MAF >0.95,
    nearby_filter = removal_nearby =='Y',
    
    filter_group = case_when(
      maf_filter   ~ "MAF < 0.05",
      info_filter  ~ "info Score < 0.95",
      seqq_filter  ~ "deCODE quality concern",
      indel_filter ~ "Multiallelic INDEL",
      pval_filter  ~ "HWE pvalue < 1e-6",
      nearby_filter ~ "Variants within 5bp",
      TRUE ~ "Pass"
    ),
    filter_group = factor(filter_group, levels = names(filter_colors)),

    logp = -log10(P),
    dot_size = if_else(maf_filter, 0.6, 1.0)
  )


p= ggplot(locus_df, aes(x = BP, y = logp, color = filter_group)) +
  geom_point(aes(size = dot_size), alpha = 0.8) +
  scale_size_identity(guide = "none") +
  scale_color_manual(values = filter_colors, breaks = names(filter_colors), drop = FALSE) +
  geom_vline(xintercept = lead.pos, linetype = "dashed") +
  labs(
    title = "LocusZoom-style plot",
    subtitle = paste0("chr", lead.chr, ":", lower.pos, "-", upper.pos),
    x = paste0("Chromosome ", lead.chr, " position"),
    y = expression(-log[10](P)),
    color = "Filter"
  ) +
    theme_classic()



pdf('/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/LocusZoom_colorBYfilter/Lund.preVpresent.nofilter.chr2.manhattan.2026Jun27.pdf',h=6,w=10)
plot(p)
dev.off()


## chr3:171104061:IG
lead.chr =3 
lead.pos = 171104061
step = 500000
lower.pos = lead.pos -step
upper.pos =lead.pos +step


locus_df <- merged %>%
  filter(
    CHR == lead.chr,
    BP > lower.pos,
    BP < upper.pos
  ) %>%
  mutate(
    info_filter  = infoScore_final < 0.95 |  Prop_lowMAXGP >0.05 ,
    seqq_filter  = str_detect(coalesce(SeqQFlags, ""), "Q"),
    indel_filter = MULTI_INDEL == "Y",
    pval_filter  = HWEpvalue < 1e-6,
    maf_filter   = MAF < 0.05 | MAF >0.95,
    nearby_filter = removal_nearby =='Y',
    
    filter_group = case_when(
      maf_filter   ~ "MAF < 0.05",
      info_filter  ~ "info Score < 0.95",
      seqq_filter  ~ "deCODE quality concern",
      indel_filter ~ "Multiallelic INDEL",
      pval_filter  ~ "HWE pvalue < 1e-6",
      nearby_filter ~ "Variants within 5bp",
      TRUE ~ "Pass"
    ),
    filter_group = factor(filter_group, levels = names(filter_colors)),

    logp = -log10(P),
    dot_size = if_else(maf_filter, 0.6, 1.0)
  )

p= ggplot(locus_df, aes(x = BP, y = logp, color = filter_group)) +
  geom_point(aes(size = dot_size), alpha = 0.8) +
  scale_size_identity(guide = "none") +
  scale_color_manual(values = filter_colors, breaks = names(filter_colors), drop = FALSE) +
  geom_vline(xintercept = lead.pos, linetype = "dashed") +
  labs(
    title = "LocusZoom-style plot",
    subtitle = paste0("chr", lead.chr, ":", lower.pos, "-", upper.pos),
    x = paste0("Chromosome ", lead.chr, " position"),
    y = expression(-log[10](P)),
    color = "Filter"
  ) +
    theme_classic()



pdf('/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/LocusZoom_colorBYfilter/Lund.preVpresent.nofilter.chr3.manhattan.2026Jun27.pdf',h=6,w=10)
plot(p)
dev.off()

## chr4:139027159:SG
lead.chr =4
lead.pos = 139027159
step = 500000
lower.pos = lead.pos -step
upper.pos =lead.pos +step


locus_df <- merged %>%
  filter(
    CHR == lead.chr,
    BP > lower.pos,
    BP < upper.pos
  ) %>%
  mutate(
    info_filter  = infoScore_final < 0.95 |  Prop_lowMAXGP >0.05 ,
    seqq_filter  = str_detect(coalesce(SeqQFlags, ""), "Q"),
    indel_filter = MULTI_INDEL == "Y",
    pval_filter  = HWEpvalue < 1e-6,
    maf_filter   = MAF < 0.05 | MAF >0.95,
    nearby_filter = removal_nearby =='Y',
    
    filter_group = case_when(
      maf_filter   ~ "MAF < 0.05",
      info_filter  ~ "info Score < 0.95",
      seqq_filter  ~ "deCODE quality concern",
      indel_filter ~ "Multiallelic INDEL",
      pval_filter  ~ "HWE pvalue < 1e-6",
      nearby_filter ~ "Variants within 5bp",
      TRUE ~ "Pass"
    ),
    filter_group = factor(filter_group, levels = names(filter_colors)),

    logp = -log10(P),
    dot_size = if_else(maf_filter, 0.6, 1.0)
  )

p= ggplot(locus_df, aes(x = BP, y = logp, color = filter_group)) +
  geom_point(aes(size = dot_size), alpha = 0.8) +
  scale_size_identity(guide = "none") +
  scale_color_manual(values = filter_colors, breaks = names(filter_colors), drop = FALSE) +
  geom_vline(xintercept = lead.pos, linetype = "dashed") +
  labs(
    title = "LocusZoom-style plot",
    subtitle = paste0("chr", lead.chr, ":", lower.pos, "-", upper.pos),
    x = paste0("Chromosome ", lead.chr, " position"),
    y = expression(-log[10](P)),
    color = "Filter"
  ) +
    theme_classic()



pdf('/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/LocusZoom_colorBYfilter/Lund.preVpresent.nofilter.chr4.manhattan.2026Jun27.pdf',h=6,w=10)
plot(p)
dev.off()


## chr6:32655190:SG
lead.chr =6
lead.pos = 32655190
step = 500000
lower.pos = lead.pos -step
upper.pos =lead.pos +step

locus_df <- merged %>%
  filter(
    CHR == lead.chr,
    BP > lower.pos,
    BP < upper.pos
  ) %>%
  mutate(
    info_filter  = infoScore_final < 0.95 |  Prop_lowMAXGP >0.05 ,
    seqq_filter  = str_detect(coalesce(SeqQFlags, ""), "Q"),
    indel_filter = MULTI_INDEL == "Y",
    pval_filter  = HWEpvalue < 1e-6,
    maf_filter   = MAF < 0.05 | MAF >0.95,
    nearby_filter = removal_nearby =='Y',
    
    filter_group = case_when(
      maf_filter   ~ "MAF < 0.05",
      info_filter  ~ "info Score < 0.95",
      seqq_filter  ~ "deCODE quality concern",
      indel_filter ~ "Multiallelic INDEL",
      pval_filter  ~ "HWE pvalue < 1e-6",
      nearby_filter ~ "Variants within 5bp",
      TRUE ~ "Pass"
    ),
    filter_group = factor(filter_group, levels = names(filter_colors)),

    logp = -log10(P),
    dot_size = if_else(maf_filter, 0.6, 1.0)
  )


p= ggplot(locus_df, aes(x = BP, y = logp, color = filter_group)) +
  geom_point(aes(size = dot_size), alpha = 0.8) +
  scale_size_identity(guide = "none") +
  scale_color_manual(values = filter_colors, breaks = names(filter_colors), drop = FALSE) +
  geom_vline(xintercept = lead.pos, linetype = "dashed") +
  labs(
    title = "LocusZoom-style plot",
    subtitle = paste0("chr", lead.chr, ":", lower.pos, "-", upper.pos),
    x = paste0("Chromosome ", lead.chr, " position"),
    y = expression(-log[10](P)),
    color = "Filter"
  ) +
    theme_classic()



pdf('/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/LocusZoom_colorBYfilter/Lund.preVpresent.nofilter.chr6.manhattan.2026Jun27.pdf',h=6,w=10)
plot(p)
dev.off()


## chr14:22193127:SG
lead.chr =14
lead.pos = 22193127
step = 500000
lower.pos = lead.pos -step
upper.pos =lead.pos +step

locus_df <- merged %>%
  filter(
    CHR == lead.chr,
    BP > lower.pos,
    BP < upper.pos
  ) %>%
  mutate(
    info_filter  = infoScore_final < 0.95 |  Prop_lowMAXGP >0.05 ,
    seqq_filter  = str_detect(coalesce(SeqQFlags, ""), "Q"),
    indel_filter = MULTI_INDEL == "Y",
    pval_filter  = HWEpvalue < 1e-6,
    maf_filter   = MAF < 0.05 | MAF >0.95,
    nearby_filter = removal_nearby =='Y',
    
    filter_group = case_when(
      maf_filter   ~ "MAF < 0.05",
      info_filter  ~ "info Score < 0.95",
      seqq_filter  ~ "deCODE quality concern",
      indel_filter ~ "Multiallelic INDEL",
      pval_filter  ~ "HWE pvalue < 1e-6",
      nearby_filter ~ "Variants within 5bp",
      TRUE ~ "Pass"
    ),
    filter_group = factor(filter_group, levels = names(filter_colors)),

    logp = -log10(P),
    dot_size = if_else(maf_filter, 0.6, 1.0)
  )


p= ggplot(locus_df, aes(x = BP, y = logp, color = filter_group)) +
  geom_point(aes(size = dot_size), alpha = 0.8) +
  scale_size_identity(guide = "none") +
  scale_color_manual(values = filter_colors, breaks = names(filter_colors), drop = FALSE) +
  geom_vline(xintercept = lead.pos, linetype = "dashed") +
  labs(
    title = "LocusZoom-style plot",
    subtitle = paste0("chr", lead.chr, ":", lower.pos, "-", upper.pos),
    x = paste0("Chromosome ", lead.chr, " position"),
    y = expression(-log[10](P)),
    color = "Filter"
  ) +
    theme_classic()



pdf('/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/LocusZoom_colorBYfilter/Lund.preVpresent.nofilter.chr14.manhattan.2026Jun27.pdf',h=6,w=10)
plot(p)
dev.off()


### chr12:108886229:SG
## trond.beforeNearby = read.table("/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/filter/Trondheim-norwegian.pre_vs_present.mm.7_pcs.maf005.rmMULTI-INDEL.rmLOWINFO095.rmLOWGP.rmQ.rmHWE1e-6.modernNorwaHWE1e-6.tsv",header=T)
## trond.afterNearby = fread("/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/filter/Trondheim-norwegian.pre_vs_present.mm.7_pcs.maf005.rmMULTI-INDEL.rmLOWINFO095.rmLOWGP.rmQ.rmHWE1e-6.modernNorwaHWE1e-6.rmwithin5bp.tsv",header=T)

## nearby.bad.marker = trond.beforeNearby$rs[! trond.beforeNearby$rs %in% trond.afterNearby$rs]
## nearby.df = data.frame(rs = nearby.bad.marker,removal_nearby= 'Y')


## trond = fread("/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/Trondheim-norwegian.99610455.2026-02-19.pre_vs_present.mm.7_pcs.assoc.min.txt.gz",header=T)
## dim(trond)

## trond$chr <- gsub(":.*", "", trond$rs)
## trond$chr = gsub("chr","",trond$chr)
## trond = trond %>% filter(chr != 'X')
## trond$ps <- as.numeric(str_split_fixed(trond$rs, ":", 3)[, 2])
## trond = trond %>% arrange(chr,ps)
## dim(trond)

## t_manhattan2 <- trond[, .(CHR = as.numeric(chr), 
##                      BP = ps, 
##                      P = p_score,
##                      SNP = rs,
##                      MAF = af)]
## dim(t_manhattan2)

## merged2 = t_manhattan2 %>% left_join(info[,c(1,2,9:14)],by=c('SNP'= 'ID'))
## dim(merged2)
## head(merged)

## merged2 = merged2 %>% left_join(hwe %>% select(ID2,pvalue),by=c('ID2' = 'ID2'))
## names(merged2)[ncol(merged2)] = 'HWEpvalue'

## merged2 = merged2 %>% left_join(nearby.df,by=c('SNP'='rs'))


## lead.chr =12
## lead.pos = 108886229
## step = 500000
## lower.pos = lead.pos -step
## upper.pos =lead.pos +step

## locus_df <- merged %>%
##   filter(
##     CHR == lead.chr,
##     BP > lower.pos,
##     BP < upper.pos
##   ) %>%
##   mutate(
##     info_filter  = infoScore_final < 0.95 |  Prop_lowMAXGP >0.05 ,
##     seqq_filter  = str_detect(coalesce(SeqQFlags, ""), "Q"),
##     indel_filter = MULTI_INDEL == "Y",
##     pval_filter  = HWEpvalue < 1e-6,
##     maf_filter   = MAF < 0.05 | MAF >0.95,
##     nearby_filter = removal_nearby =='Y',
    
##     filter_group = case_when(
##       maf_filter   ~ "MAF < 0.05",
##       info_filter  ~ "info Score < 0.95",
##       seqq_filter  ~ "deCODE quality concern",
##       indel_filter ~ "Multiallelic INDEL",
##       pval_filter  ~ "HWE pvalue < 1e-6",
##       nearby_filter ~ "Variants within 5bp",
##       TRUE ~ "Pass"
##     ),
##     filter_group = factor(filter_group, levels = names(filter_colors)),

##     logp = -log10(P),
##     dot_size = if_else(maf_filter, 0.6, 1.0)
##   )


## p= ggplot(locus_df, aes(x = BP, y = logp, color = filter_group)) +
##   geom_point(aes(size = dot_size), alpha = 0.8) +
##   scale_size_identity(guide = "none") +
##   scale_color_manual(values = filter_colors, breaks = names(filter_colors), drop = FALSE) +
##   geom_vline(xintercept = lead.pos, linetype = "dashed") +
##   labs(
##     title = "LocusZoom-style plot",
##     subtitle = paste0("chr", lead.chr, ":", lower.pos, "-", upper.pos),
##     x = paste0("Chromosome ", lead.chr, " position"),
##     y = expression(-log[10](P)),
##     color = "Filter"
##   ) +
##     theme_classic()



## pdf('/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/LocusZoom_colorBYfilter/Trondheim.preVpresent.nofilter.chr12.manhattan.2026Jun27.pdf',h=5,w=15)
## plot(p)

## dev.off()
