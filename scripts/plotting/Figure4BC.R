library(dplyr)
library(data.table)
library(tidyr)
library(ComplexUpset)
library(ggplot2)
library(colorspace)
library(cowplot)#
library(scales)
fade_colours <- function(col, fade = 0.5) {
  scales::alpha(col, 1 - fade)
}

# masterFile <- "https://docs.google.com/spreadsheets/d/1V42OBU_UWXDCuF4P_75L8-GaSSlBh0pGvts0hVexTBY/edit"
## master <- as.data.frame(gsheet::gsheet2tbl(masterFile))
## df  =  master %>% filter(sampling_site == 'Lund' & included_analysis =='yes' )  %>% dplyr::select(id,Age_assignment)

## info = read.table("/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Summary.ancestry.tsv",header=T,sep='\t')
info = read.table("/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Summary.ancestry.2026April28.tsv",header=T,sep='\t')

### Lund
lund.info = info %>% filter(Site=='Lund')

lund.pca = read.table("/maps/projects/ilab/people/pls394/plague/ancestry_version1/PCA/output/Lund.pvalue.europe_more3_collapse.DF7.202603.csv",sep=",",header=T,row.names = 1)
colnames(lund.pca) <- gsub("_E"," E", colnames(lund.pca))
colnames(lund.pca) <- gsub("_I"," & I", colnames(lund.pca))
colnames(lund.pca) <- gsub("_B"," & B", colnames(lund.pca))
lund.outlier = lund.info %>% filter(Group == 'non-Scandinavian') %>% pull(ID)
lund.pca.outlier = lund.pca %>% filter(row.names(lund.pca) %in% lund.outlier)
lund.pca.outlier$Age <- lund.info$Age_assignment[match(rownames(lund.pca.outlier), lund.info$ID)]

df <- lund.pca.outlier
region_cols <- setdiff(colnames(df), "Age")
#region_cols = region_cols[-which(region_cols =='Scandinavia')]

# convert p-values to membership
df_sets <- df %>%
  mutate(across(all_of(region_cols), ~ . >= 0.01))

# keep only regions with at least one TRUE
valid_regions <- region_cols[
  colSums(df_sets[region_cols]) >= 1
]

df_sets_filtered <- df_sets %>%
    select(all_of(valid_regions), Age)



# Recode Age values
df_sets_filtered$Age <- dplyr::recode(
  df_sets_filtered$Age,
  "Pre" = "Pre-BD",
  "During" = "During BD",
  "Post" = "Post-BD"
)

df_sets_filtered$Age <- factor(df_sets_filtered$Age,
                               levels = c("Pre-BD", "During BD", "Post-BD"))

## check what inds are only consistent with Scandinavia
region_cols <- setdiff(names(df_sets_filtered), "Age")
df_sets_filtered[
  df_sets_filtered$Scandinavia == TRUE &
    rowSums(df_sets_filtered[region_cols]) == 1,  ]


## upset plot
base_col <- "#1122aa"

age_cols <- c(
    "Pre-BD" = fade_colours('#0b1874', 1-0.4),
    "During BD" = base_col,
    "Post-BD" = '#0b1874'
  ## "During BD" = base_col,
  ## "Post-BD" = "#0b1874"

)



upset.lund <- upset(
  df_sets_filtered,
  intersect = valid_regions,
  width_ratio = 0.35,
  height_ratio = 1.2,
  name = "Present-day references",
  guides = "over",

  base_annotations = list(
    "Intersection size" =
      intersection_size(
        counts = FALSE,
        mapping = aes(fill = Age),
        width = 0.7,
        position = position_dodge2(width = 0.7, preserve = "single")
      ) +
      scale_fill_manual(
        values = age_cols,
        breaks = c("Pre-BD", "During BD", "Post-BD"),
        drop = FALSE
      ) +
      theme(
        panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),

        ## axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        ## axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14),
      )
  ),

  set_sizes = upset_set_size(
      geom = geom_bar(
          aes(x = group, fill = Age),
          width = 0.8,
          position = position_dodge2(width = 0.8, preserve = "single")
      ),
      position = "left"
  ) +
    scale_fill_manual(
        values = age_cols,
        breaks = c("Pre-BD", "During BD", "Post-BD"),
        drop = FALSE
    ) +
      theme(
      axis.text.x = element_text(size = 14),
      ## axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      ## axis.title.y = element_text(size = 16)
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14),
    )
)  + # end upset
    theme(
        legend.position = "none",

        text = element_text(size = 14),
        ## axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        ## axis.title.y = element_text(size = 16),

        plot.title = element_text(size = 18),
        strip.text = element_text(size = 16),
    )


upset.lund


### Trondheim
trondheim.info = info %>% filter(Site=='Trondheim') 


trondheim.pca = read.table("/maps/projects/ilab/people/pls394/plague/ancestry_version1/PCA/output/Trondheim.pvalue.europe_more3_collapse.DF7.202603.csv",sep=",",header=T,row.names = 1)
colnames(trondheim.pca) <- gsub("_E"," E", colnames(lund.pca))
colnames(trondheim.pca) <- gsub("_I"," & I", colnames(lund.pca))
colnames(trondheim.pca) <- gsub("_B"," & B", colnames(lund.pca))

trondheim.outlier = trondheim.info %>% filter(Group == 'non-Scandinavian') %>% pull(ID)
trondheim.pca.outlier = trondheim.pca %>% filter(row.names(trondheim.pca) %in% trondheim.outlier)
trondheim.pca.outlier$Age <- trondheim.info$Age_assignment[match(rownames(trondheim.pca.outlier), trondheim.info$ID)]

df <- trondheim.pca.outlier
region_cols <- setdiff(colnames(df), "Age")
##region_cols = region_cols[-which(region_cols =='Scandinavia')]

# convert p-values to membership
df_sets <- df %>%
  mutate(across(all_of(region_cols), ~ . >= 0.01))

# keep only regions with at least one TRUE
valid_regions <- region_cols[
  colSums(df_sets[region_cols]) >= 1
]

df_sets_filtered <- df_sets %>%
    select(all_of(valid_regions), Age)

# Recode Age values
df_sets_filtered$Age <- dplyr::recode(
  df_sets_filtered$Age,
  "Pre" = "Pre-BD",
  "Post" = "Post-BD",
  "uncertain" = 'Uncertain'
)

df_sets_filtered$Age <- factor(df_sets_filtered$Age,
                               levels = c("Pre-BD","Post-BD",'Uncertain'))

## upset plot
base_col <- "#0b6b16"

age_cols <- c(
  "Pre-BD" = fade_colours("#0b6b16", 1-0.3),
  "Post-BD" = "#0b6b16",
  "Uncertain" = fade_colours("#0b6b16",1- 0.65)
)

upset.trondheim <- upset(
  df_sets_filtered,
  intersect = valid_regions,
  width_ratio = 0.35,
  height_ratio = 1.2,
  name = "Present-day references",
  guides = "over",

  base_annotations = list(
    "Intersection size" =
      intersection_size(
        counts = FALSE,
        mapping = aes(fill = Age),
        width = 0.7,
        position = position_dodge2(width = 0.7, preserve = "single")
      ) +
      scale_fill_manual(
        values = age_cols,
        breaks = c("Pre-BD","Post-BD","Uncertain"),
        drop = FALSE
      ) +
      theme(
        panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),

        ## axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        ## axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14),
      )
  ),

  set_sizes = upset_set_size(
      geom = geom_bar(
          aes(x = group, fill = Age),
          width = 0.8,
          position = position_dodge2(width = 0.8, preserve = "single")
      ),
      position = "left"
  ) +
    scale_fill_manual(
        values = age_cols,
        breaks = c("Pre-BD","Post-BD","Uncertain"),
        drop = FALSE
    ) +
      theme(
      axis.text.x = element_text(size = 14),
      ## axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      ## axis.title.y = element_text(size = 16)
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14),
    )
)  + # end upset
    theme(
        legend.position = "none",

        text = element_text(size = 14),
        ## axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        ## axis.title.y = element_text(size = 16),

        plot.title = element_text(size = 18),
        strip.text = element_text(size = 16),
    )


upset.trondheim


### Vilnius
vilnius.info = info %>% filter(Site=='Vilnius')
vilnius.info[vilnius.info$Age_assignment =='During','Age_assignment'] = 'Post'
vilnius.pca = read.table("/maps/projects/ilab/people/pls394/plague/ancestry_version1/PCA/output/Vilnius.pvalue.europe_more3_collapse.DF7.202603.csv",sep=",",header=T,row.names = 1)
colnames(vilnius.pca) <- gsub("_E"," E", colnames(vilnius.pca))
colnames(vilnius.pca) <- gsub("_I"," & I", colnames(vilnius.pca))
colnames(vilnius.pca) <- gsub("_B"," & B", colnames(vilnius.pca))

vilnius.outlier = vilnius.info %>% filter(Group == 'non-Baltic') %>% pull(ID)
vilnius.pca.outlier = vilnius.pca %>% filter(row.names(vilnius.pca) %in% vilnius.outlier)
vilnius.pca.outlier$Age <- vilnius.info$Age_assignment[match(rownames(vilnius.pca.outlier), vilnius.info$ID)]

df <- vilnius.pca.outlier
region_cols <- setdiff(colnames(df), "Age")
##region_cols = region_cols[-which(region_cols =='Baltic')] # it does not matter


# convert p-values to membership
df_sets <- df %>%
  mutate(across(all_of(region_cols), ~ . >= 0.01))

# keep only regions with at least one TRUE
valid_regions <- region_cols[
  colSums(df_sets[region_cols]) >= 1
]

df_sets_filtered <- df_sets %>%
    select(all_of(valid_regions), Age)

df_sets_filtered$Age <- factor(df_sets_filtered$Age,
                               levels = c("Pre", "Post"))
# Recode Age values
df_sets_filtered$Age <- dplyr::recode(
  df_sets_filtered$Age,
  "Pre" = "Pre-BD",
  "Post" = "Post-BD"
)

df_sets_filtered$Age <- factor(df_sets_filtered$Age,
                               levels = c("Pre-BD", "Post-BD"))

## upset plot
base_col <- "#7e0c3f"


age_cols <- c(
  "Pre-BD" = fade_colours("#7e0c3f", 0.8),
  "Post-BD" = "#7e0c3f"

)

upset.vilnius <- upset(
  df_sets_filtered,
  intersect = valid_regions,
  width_ratio = 0.35,
  height_ratio = 1.2,
  name = "Present-day references",
  guides = "over",

  base_annotations = list(
    "Intersection size" =
      intersection_size(
        counts = FALSE,
        mapping = aes(fill = Age),
        width = 0.7,
        position = position_dodge2(width = 0.7, preserve = "single")
      ) +
      scale_fill_manual(
        values = age_cols,
        breaks = c("Pre-BD","Post-BD"),
        drop = FALSE
      ) +
      theme(
        panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),

        ## axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        ## axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14),
      )
  ),

  set_sizes = upset_set_size(
      geom = geom_bar(
          aes(x = group, fill = Age),
          width = 0.8,
          position = position_dodge2(width = 0.8, preserve = "single")
      ),
      position = "left"
  ) +
    scale_fill_manual(
        values = age_cols,
        breaks = c("Pre-BD","Post-BD"),
        drop = FALSE
    ) +
      theme(
      axis.text.x = element_text(size = 14),
      ## axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      ## axis.title.y = element_text(size = 16)
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14),
    )
)  + # end upset
    theme(
        legend.position = "none",

        text = element_text(size = 14),
        ## axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        ## axis.title.y = element_text(size = 16),

        plot.title = element_text(size = 18),
        strip.text = element_text(size = 16),
    )


upset.vilnius




# Combine with extra space between them using rel_widths
# Add 0.1 width columns of empty space between each plot
plots <- list(
  upset.lund,
  upset.trondheim,
  upset.vilnius
)

# Combine with extra space between them using rel_widths
# Add 0.1 width columns of empty space between each plot
row1 <- plot_grid(
  plots[[1]], NULL, plots[[2]], NULL, plots[[3]],
  nrow = 1,
  rel_widths = c(1.3, 0.1/2, 1, 0.1/2, 1),  # 0.1 = 10% width for spacing
 # labels = c("A", "", "B", "", "C"),
  label_size = 14
)
row1

pdf('NonOutlier.PCA.upset.pdf',h=6,w=20)
row1
dev.off()


##### add row 2
library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggh4x)
library(pheatmap)
library(VennDiagram)
library(RColorBrewer)
library(openxlsx)
library(ComplexHeatmap)
library(patchwork)
library(scales)


masterFile <- "https://docs.google.com/spreadsheets/d/1HrTrz-77GykCpsILlXWiNyujMQRFFlSQnu0aFJUGCOs/edit"
master <- as.data.frame(gsheet::gsheet2tbl(masterFile))
samples.df <-master %>% filter((`Sampling Site` == 'Vilnius' | `Sampling Site` =='Trondheim' | `Sampling Site` == 'Lund') & `Included in Analyses?` =='yes' )  %>% select(Identity,`Sampling Site`,`Cohort Assignment`)
samples.df[which(samples.df$`Sampling Site` == 'Vilnius' & samples.df$`Cohort Assignment`=='during-BD'),'Cohort Assignment'] = 'post-BD'
samples <- samples.df %>% pull(Identity)

summary = read.table("/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Summary.ancestry.tsv",header=T,sep="\t")
samples <- summary %>% filter(Group == 'non-Scandinavian' | Group =='non-Baltic'| Group == 'nonlocal Scandinavian') %>% pull(ID)
## read UKbiobank ibd
t<-read.table('/maps/projects/ilab/people/pls394/plague/data_decode_2024Nov/IBD/UK_biobank_2026Fig/IBD_ALL_v_xbi_and_nonxbi_regrouped.masked.7cM.100maskedLociPerOrigCm.country_decaf_plus_changed_xbi.summary_to_send.txt',header=T)
nonxbi <- t %>%
    filter(totally_excluded == FALSE | is.na(totally_excluded)) %>%
    filter(dropped_for_figures == FALSE | is.na(dropped_for_figures))  

nonxbi$ID2G %>% unique()

regions = c("Lithuania","Estonia_and_Latvia","Belarus_Russia","Poland","Czechia_Hungary_Slovakia","France","Ashkenazi_ancestry","England","Wales","Sctoland","Ireland","Finland")
regions = c( "Ashkenazi_ancestry" ,"Baltic" ,"Belarus_Russia_Ukraine","England" ,"Finland" , "France"   ,  "Ireland"        ,        "Poland" ,                "Scotland" ,  "Southeast_Europe"   ,    "Wales" )
nonxbi <- nonxbi %>% filter(ID2Grp %in% regions)
nonxbi$ID2Grp <- factor(nonxbi$ID2Grp, levels = regions)

## plot for
samples.grp <- samples.df %>%
  mutate(
    grp = case_when(
      `Sampling Site` == "Lund" & `Cohort Assignment` == "pre-BD" ~ "Lund\nPre-BD",
      `Sampling Site` == "Lund" & `Cohort Assignment` == "during-BD" ~ "Lund\nDuring BD",
      `Sampling Site` == "Lund" & `Cohort Assignment` == "post-BD" ~ "Lund\nPost-BD",
      `Sampling Site` == "Lund" & `Cohort Assignment` == "uncertain" ~ "Lund\nUncertain",
      `Sampling Site` == "Trondheim" & `Cohort Assignment` =='pre-BD'  ~ "Trondheim\nPre-BD",
      `Sampling Site` == "Trondheim" & `Cohort Assignment` =='post-BD'  ~ "Trondheim\nPost-BD",
      `Sampling Site` == "Trondheim" & `Cohort Assignment` =='uncertain'  ~ "Trondheim\nUncertain",
      `Sampling Site` == "Vilnius" & `Cohort Assignment` =='pre-BD'  ~ "Vilnius\nPre-BD",
      `Sampling Site` == "Vilnius" & `Cohort Assignment` =='post-BD'  ~ "Vilnius\nPost-BD",
      `Sampling Site` == "Vilnius"& `Cohort Assignment` =='uncertain'  ~ "Vilnius\nUncertain",

      TRUE ~ NA_character_
    )
  ) %>%
  select(Identity, grp)

raw_plot_df <- nonxbi %>%
  left_join(samples.grp, by = c("ID" = "Identity")) %>%
  filter(!is.na(grp)) %>%
  mutate(
    grp = factor(
      grp,
      
      levels = c("Lund\nPre-BD","Lund\nDuring BD","Lund\nPost-BD","Lund\nUncertain","Trondheim\nPre-BD","Trondheim\nPost-BD","Trondheim\nUncertain","Vilnius\nPre-BD","Vilnius\nPost-BD","Vilnius\nUncertain")

    )
  )


plot_df <- raw_plot_df %>%
    ## filter(ID2Grp %in% c("Wales", "Ireland", "Ashkenazi_ancestry", "Finland")) %>%
    filter(ID2Grp %in% c("Ashkenazi_ancestry")) %>%
  mutate(
    label = case_when(
      ID2Grp == "Ashkenazi_ancestry" & ID == "LUN114" ~ "LUN114",
      TRUE ~ NA_character_
    )
  )


bar_fill_v2 <- c(
  `Lund\nPre-BD`        = fade_colours("#0b1874", 0.4),
  `Lund\nDuring BD`     = "#1122aa",
  `Lund\nPost-BD`       = "#0b1874",
  `Lund\nUncertain`     = fade_colours("#0b1874", 0.65),

  `Trondheim\nPre-BD`   = fade_colours("#0b6b16", 0.3),
  `Trondheim\nPost-BD`  = "#0b6b16",
  `Trondheim\nUncertain` = fade_colours("#0b6b16", 0.65),

  `Vilnius\nPre-BD`     = fade_colours("#7e0c3f", 0.2),
  `Vilnius\nPost-BD`    = "#7e0c3f"
)

plot_df$grp <- factor(
  plot_df$grp,
  levels = names(bar_fill_v2)
)

plot_df$ID2Grp <- gsub("_", " ", plot_df$ID2Grp)

ash.p <- ggplot(plot_df, aes(x = grp, y = shareID2Mean)) +
  geom_point(
    aes(color = grp),
    alpha = 0.7,
    position = position_jitter(width = 0.15, height = 0),
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = subset(plot_df, !is.na(label)),
    aes(label = label),
    size = 6,
    color = "black",
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  scale_color_manual(values = bar_fill_v2, drop = FALSE) +
  ggh4x::facet_wrap2(
    ~ ID2Grp,
    scales = "free_y",
    ncol = 3,
    axes = "all"
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      ID2Grp %in% c("Wales", "Ireland") ~
        scale_y_continuous(limits = c(0, 1.5))
    )
  ) +
  labs(
    x = NULL,
    y = "Mean IBD (cM)",
    title = ""
  ) +
  theme_bw() +
  theme(
      axis.text.x = element_text(
          size = 14,
          angle = 45,
          hjust = 1,
          vjust = 1
    ),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    panel.grid = element_blank()
  )


## read OminExpress
full = read.table("/maps/projects/ilab/people/pls394/plague/data_decode_2024Nov/IBD/omni_2026Fig/IBD_ALL_v_omni.masked.6cM.100maskedLociPerOrigCm.summary_to_send.txt",header=T)


raw_plot_df <- full %>%
  left_join(samples.grp, by = c("ID" = "Identity")) %>%
  filter(!is.na(grp)) %>%
  mutate(
    grp = factor(
      grp,
      
      levels = c("Lund\nPre-BD","Lund\nDuring BD","Lund\nPost-BD","Lund\nUncertain","Trondheim\nPre-BD","Trondheim\nPost-BD","Trondheim\nUncertain","Vilnius\nPre-BD","Vilnius\nPost-BD","Vilnius\nUncertain")

    )
  )

### a refined version of Iceland plot ###
plot_df <- raw_plot_df %>%
  filter(ID2Grp %in% c("Iceland")) %>%
  mutate(
    label = case_when(
##        county == "Iceland" & id %in% c("AV29","SK018","SK083","SK115","SK226","SK317","SK339","SK340","SK381") ~ id,
        ID2Grp == "Iceland" & ID %in% c("AV29","LUN138","SK083","SK115","SK226","SK317","SK339","SK340") ~ ID,

      TRUE ~ NA_character_
    )
  )

plot_df$grp <- factor(
  plot_df$grp,
  levels = names(bar_fill_v2)
)

plot_df$ID2Grp <- gsub("_", " ", plot_df$ID2Grp)

icelandic.p <- ggplot(plot_df, aes(x = grp, y = shareID2Mean)) +
  geom_point(
    aes(color = grp),
    alpha = 0.7,
    position = position_jitter(width = 0.15, height = 0),
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = subset(plot_df, !is.na(label)),
    aes(label = label),
    size = 6,
    color = "black",
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  scale_color_manual(values = bar_fill_v2, drop = FALSE) +
  ggh4x::facet_wrap2(
    ~ ID2Grp,
    scales = "free_y",
    ncol = 3,
    axes = "all"
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      ID2Grp %in% c("Wales", "Ireland") ~
        scale_y_continuous(limits = c(0, 1.5))
    )
  ) +
  labs(
    x = NULL,
    y = "Mean IBD (cM)",
    title = ""
  ) +
  theme_bw() +
  theme(
      axis.text.x = element_text(
          size = 14,
          angle = 45,
          hjust = 1,
          vjust = 1
    ),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    panel.grid = element_blank()
  )


### a refined version of NO_South plot ###
plot_df <- raw_plot_df %>%
  mutate(
    ID2Grp = if_else(ID2Grp == "NO_South", "NO:South", ID2Grp)
  ) %>%
  filter(ID2Grp %in% c("NO:South")) %>%
  mutate(
    label = case_when(
      ## from different ID2Grp
      ID2Grp == "NO:South" & ID %in% c("LUN391", "LUN269", "LUN210","LUN200","LUN359") ~ ID,
      ## from different regions within the same country
      ID2Grp == "NO:South" & ID%in% c("SK288", "SK106", "SK096", "SK047") ~ ID,
      TRUE ~ NA_character_
    )
  )

plot_df$grp <- factor(
  plot_df$grp,
  levels = names(bar_fill_v2)
)

plot_df$ID2Grp <- gsub("_", " ", plot_df$ID2Grp)

NOsouth.p <- ggplot(plot_df, aes(x = grp, y = shareID2Mean)) +
  geom_point(
    aes(color = grp),
    alpha = 0.7,
    position = position_jitter(width = 0.15, height = 0),
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = subset(plot_df, !is.na(label)),
    aes(label = label),
    size = 6,
    color = "black",
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  scale_color_manual(values = bar_fill_v2, drop = FALSE) +
  ggh4x::facet_wrap2(
    ~ ID2Grp,
    scales = "free_y",
    ncol = 3,
    axes = "all"
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      ID2Grp %in% c("Wales", "Ireland") ~
        scale_y_continuous(limits = c(0, 1.5))
    )
  ) +
  labs(
    x = NULL,
    y = "Mean IBD (cM)",
    title = ""
  ) +
  theme_bw() +
  theme(
      axis.text.x = element_text(
          size = 14,
          angle = 45,
          hjust = 1,
          vjust = 1
    ),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    panel.grid = element_blank()
  )

row2 <- plot_grid(
  ash.p, NULL, icelandic.p, NULL, NOsouth.p,
  nrow = 1,
  rel_widths = c(1, 0.05, 1, 0.05, 1)
)

combined_all <- plot_grid(
  row1,
  row2,
  ncol = 1,
  align = "none",
  rel_heights = c(1, .8)
)

ggsave(
  filename = "Fig4_B-C.2026Jun04.svg",
  plot = combined_all,
  width = 20,
  height = 11,
  units = "in",
  dpi = 300
)
