library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(grid)
library(GetoptLong)  
library(circlize)
library(scico)

############# Lund #################
plot.mat.2 <- read.table("/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Lund.IBD_standard.2026April28.tsv",header = T,row.names=1,sep="\t") %>% as.matrix
info <- read.table("/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Lund.IBD_info.2026April28.tsv",header=T,sep="\t")
table(info$Group,info$Age_assignment)
  
#### pre
###pre-nonScand
##pca
pre_nonScand_pca <- info %>% filter(Age_assignment =='Pre' & Group == 'non-Scandinavian' & Method == 'PCA') %>% pull(ID)
pre_nonScand_pca.data = plot.mat.2[,pre_nonScand_pca]
hm.pre_nonScand_pca = heatmap(pre_nonScand_pca.data)
##ibd
pre_nonScand_ibd <- info %>% filter(Age_assignment =='Pre' & Group == 'non-Scandinavian' & Method == 'IBD') %>% pull(ID)
pre_nonScand_ibd.data = plot.mat.2[,pre_nonScand_ibd]
hm.pre_nonScand_ibd = heatmap(pre_nonScand_ibd.data)

###pre-nonLocal
##pca
pre_nonLocal_pca <- info %>% filter(Age_assignment =='Pre' & Group == 'nonlocal Scandinavian'& Method =='PCA') %>% pull(ID)
pre_nonLocal_pca.data = plot.mat.2[,pre_nonLocal_pca]
hm.pre_nonLocal_pca = heatmap(pre_nonLocal_pca.data)
##ibd
pre_nonLocal_ibd <- info %>% filter(Age_assignment =='Pre' & Group == 'nonlocal Scandinavian'& Method=='IBD') %>% pull(ID)
pre_nonLocal_ibd.data = plot.mat.2[,pre_nonLocal_ibd]
hm.pre_nonLocal_ibd = heatmap(pre_nonLocal_ibd.data)

###pre-Local
##pca
pre_Local <- info %>% filter(Age_assignment =='Pre' & Group == 'local Scandinavian') %>% pull(ID)

pre_Local.data = plot.mat.2[,pre_Local]
hm.pre_Local = heatmap(pre_Local.data)

col_order.pre <- c(colnames(pre_nonScand_pca.data)[hm.pre_nonScand_pca$colInd],
                   colnames(pre_nonScand_ibd.data)[hm.pre_nonScand_ibd$colInd],
                   colnames(pre_nonLocal_pca.data)[hm.pre_nonLocal_pca$colInd],
                   colnames(pre_nonLocal_ibd.data)[hm.pre_nonLocal_ibd$colInd],
                   colnames(pre_Local.data)[hm.pre_Local$colInd])


#### during
###during-nonScand
##pca
during_nonScand_pca <- info %>% filter(Age_assignment =='During' & Group == 'non-Scandinavian' & Method == 'PCA') %>% pull(ID)
during_nonScand_pca.data = plot.mat.2[,during_nonScand_pca]
hm.during_nonScand_pca = heatmap(during_nonScand_pca.data)
##ibd
during_nonScand_ibd <- info %>% filter(Age_assignment =='During' & Group == 'non-Scandinavian' & Method == 'IBD') %>% pull(ID)
# one sample in during, skip heatmap
during_nonScand_ibd.data = plot.mat.2[,during_nonScand_ibd]
##hm.during_nonScand_ibd = heatmap(during_nonScand_ibd.data)

###during-nonLocal
##pca
during_nonLocal_pca <- info %>% filter(Age_assignment =='During' & Group == 'nonlocal Scandinavian'& Method=='PCA') %>% pull(ID)
# no sample, skip heatmap
during_nonLocal_pca.data = plot.mat.2[,during_nonLocal_pca]
##hm.during_nonLocal_pca = heatmap(during_nonLocal_pca.data)
##ibd
during_nonLocal_ibd <- info %>% filter(Age_assignment =='During' & Group == 'nonlocal Scandinavian'& Method=='IBD') %>% pull(ID)
during_nonLocal_ibd.data = plot.mat.2[,during_nonLocal_ibd]
hm.during_nonLocal_ibd = heatmap(during_nonLocal_ibd.data)

##during-Local
##pca
during_Local <- info %>% filter(Age_assignment =='During' & Group == 'local Scandinavian' ) %>% pull(ID)
during_Local.data = plot.mat.2[,during_Local]
hm.during_Local = heatmap(during_Local.data)

col_order.during <- c(colnames(during_nonScand_pca.data)[hm.during_nonScand_pca$colInd],
                      during_nonScand_ibd,

                      colnames(during_nonLocal_ibd.data)[hm.during_nonLocal_ibd$colInd],
                      colnames(during_Local.data)[hm.during_Local$colInd])


#### post
###post-nonScand
##pca
post_nonScand_pca <- info %>% filter(Age_assignment =='Post' & Group == 'non-Scandinavian' & Method == 'PCA') %>% pull(ID)
post_nonScand_pca.data = plot.mat.2[,post_nonScand_pca]
hm.post_nonScand_pca = heatmap(post_nonScand_pca.data)
##ibd
post_nonScand_ibd <- info %>% filter(Age_assignment =='Post' & Group == 'non-Scandinavian' & Method == 'IBD') %>% pull(ID)
#only one sample, so skip heatmap
post_nonScand_ibd.data = plot.mat.2[,post_nonScand_ibd]
##hm.post_nonScand_ibd = heatmap(post_nonScand_ibd.data)

###post-nonLocal
##pca
post_nonLocal_pca <- info %>% filter( Age_assignment=='Post' & Group == 'nonlocal Scandinavian'& Method=='PCA') %>% pull(ID)
post_nonLocal_pca.data = plot.mat.2[,post_nonLocal_pca]
hm.post_nonLocal_pca = heatmap(post_nonLocal_pca.data)
##ibd
post_nonLocal_ibd <- info %>% filter(Age_assignment =='Post' & Group == 'nonlocal Scandinavian'& Method=='IBD') %>% pull(ID)
post_nonLocal_ibd.data = plot.mat.2[,post_nonLocal_ibd]
hm.post_nonLocal_ibd = heatmap(post_nonLocal_ibd.data)

###post-Local
##pca
post_Local <- info %>% filter(Age_assignment =='Post' & Group == 'local Scandinavian' ) %>% pull(ID)
post_Local.data = plot.mat.2[,post_Local]
hm.post_Local = heatmap(post_Local.data)

col_order.post <- c(colnames(post_nonScand_pca.data)[hm.post_nonScand_pca$colInd],
                    post_nonScand_ibd,
                   colnames(post_nonLocal_pca.data)[hm.post_nonLocal_pca$colInd],
                   colnames(post_nonLocal_ibd.data)[hm.post_nonLocal_ibd$colInd],
                   colnames(post_Local.data)[hm.post_Local$colInd])



col_order <- c(col_order.pre,col_order.during,col_order.post)
plot.mat.lund <- plot.mat.2[,col_order]



## make heatmap
## annotation
annotate.df.row.lund <- info[match(col_order,info$ID),]
all.equal(annotate.df.row.lund$ID,col_order)
annotate.df.row.lund$Age_assignment <- factor(annotate.df.row.lund$Age_assignment,levels=unique(annotate.df.row.lund$Age_assignment))
annotate.df.row.lund$Group <- factor(annotate.df.row.lund$Group, levels = unique(annotate.df.row.lund$Group))
annotate.df.row.lund$Pop = 'Lund'

annotate.df.col.lund <- data.frame(IBDsource = c(rep("Local Scandinavian",2),rep("non-Local Scandinavian",12),rep("non-Scandinavian",5)))


rownames(plot.mat.lund) <- gsub("_",":", rownames(plot.mat.lund))
rownames(plot.mat.lund) <- c(
  "DK", "SE:Sk", "SE:So", "SE:Gt", "SE:St", "SE:Ce", "SE:No",
  "NO:Os", "NO:So", "NO:Ce", "NO:Ea", "NO:We", "NO:No", "NO:Tr",
  "NL", "DE", "UK:Sh", "UK:Or", "IS"
)

pheatmap_colors <- colorRampPalette(rev(c("#D73027", "#FC8D59", 
                                          "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(101)
ibd_color = colorRamp2(breaks = seq(-3,4,by=0.07), colors = pheatmap_colors)

## row annotation
row_ha.left = rowAnnotation(IBDsource=annotate.df.col.lund$IBDsource,col=list(IBDsource =c("non-Local Scandinavian"= "#E9C46A", "non-Scandinavian"="#005F73", "Local Scandinavian"="#E76F51")),
                            show_legend = FALSE,show_annotation_name = FALSE)
spacer_ha <- rowAnnotation(
    spacer = anno_empty(border = FALSE, width = unit(2, "mm"))
)
combined_left_annotation <- c(row_ha.left , spacer_ha)

## column annotation
col_ha = HeatmapAnnotation(Period=annotate.df.row.lund$Age_assignment, col=list(Period=c("Pre"="blue","During"="red","Post"="green")),na_col='white',show_legend = FALSE,show_annotation_name=FALSE)
col_ha2 =  HeatmapAnnotation("Result"=annotate.df.row.lund$Group, col=list("Result"=c("non-Scandinavian"="#005F73","nonlocal Scandinavian"="#E9C46A",'local Scandinavian'="#E76F51")),na_col='white',show_legend = FALSE,show_annotation_name = TRUE,simple_anno_size  = unit(12, "mm"),annotation_name_side = "left",
                             annotation_name_gp = gpar(fontsize = 24),
                             annotation_name_offset = unit(10, "mm") )
## col_ha3 =  HeatmapAnnotation(Method=annotate.df.row.lund$Method, col=list(Method=c("PCA"="beige","IBD"="bisque")),na_col='white',show_legend = FALSE,show_annotation_name = TRUE,simple_anno_size  = unit(12, "mm"),annotation_name_side = "left",
##                              annotation_name_gp = gpar(fontsize = 24),
##                              annotation_name_offset = unit(10, "mm"))

col_ha3 =  HeatmapAnnotation(Method=annotate.df.row.lund$Method, col=list(Method=c("PCA"="#55514D","IBD"="#A39F9A")),na_col='white',show_legend = FALSE,show_annotation_name = TRUE,simple_anno_size  = unit(12, "mm"),annotation_name_side = "left",
                             annotation_name_gp = gpar(fontsize = 24),
                             annotation_name_offset = unit(10, "mm"))



ha1 = HeatmapAnnotation(
    empty1 = anno_empty(border = FALSE, height  = unit(6, "mm")),
    gp = gpar(fontsize = 16)
)

ha2 = HeatmapAnnotation(
    empty2 = anno_empty(border = FALSE, height = unit(6, "mm")),
    gp = gpar(fontsize = 16)
)
spacer1 <- HeatmapAnnotation(
  empty3 = anno_empty(width = unit(22, "mm"),border=FALSE)
)
spacer2 <- HeatmapAnnotation(
  empty4 = anno_empty(width = unit(7, "mm"),border=FALSE)
)

combined_top_annotation <- c(ha1, ha2)
combined_bottom_annotation <- c(spacer1,col_ha2, spacer2,col_ha3)

ibd_lund <- Heatmap(plot.mat.lund,
        name='IBD_Lund',
        cluster_columns=FALSE,
        cluster_rows=FALSE,
        show_row_names = TRUE,
        show_column_names=FALSE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 24),
        row_title = "IBD sharing with",  # Add this line
        row_title_side = "left",            # Position the title
        row_title_gp = gpar(fontsize = 24,fontface = "bold"),
        ## row_split = factor(c('DK',rep('SE',6),rep('NO',7),rep('nonScand-Europe',5)),levels=c("DK","SE","NO","nonScand-Europe")),
        ## row_gap = unit(c(1.5,1.5,1.5,1.5),"mm"),
        ## row_title = NULL,
        ## col = pheatmap_colors,
        col = ibd_color,
        na_col = 'white',
        column_split = list(annotate.df.row.lund$Age_assignment,annotate.df.row.lund$Group),
        column_gap = unit(c(1.5,1.5,4,1.5,1.5,4,1.5,1.5), "mm"),
        column_title = NULL,
        left_annotation = combined_left_annotation,

        ## top_annotation = combined_top_annotation,
        bottom_annotation = combined_bottom_annotation,
        ##heatmap_legend_param = list(title = NULL)  ,
        show_heatmap_legend = FALSE
        ## width = unit(ncol(plot.mat.lund)*0.8, "mm")
        )


group_block_anno = function(group, empty_anno, gp = gpar(), 
    label = NULL, label_gp = gpar()) {

    seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
    loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
    seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
    loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))

    seekViewport("global")
    grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
        just = c("left", "bottom"), gp = gp)
    if(!is.null(label)) {
        grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5, gp = label_gp,rot  = 0)
    }
}

draw(ibd_lund)

group_block_anno(7:9, "empty2", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),  # This controls the BOX
                 label = "Post",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))  # This controls the TEXT

group_block_anno(4:6, "empty2", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),
                 label = "During", 
                 label_gp = gpar(fontsize = 22, fontface = "bold"))

group_block_anno(1:3, "empty2", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),
                 label = "Pre",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))

group_block_anno(1:9, "empty1", 
                 gp = gpar(fill = NA, col = NA, lwd = 2),
                 label = "Lund",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))




### lund PCA ###
lund = read.table("/maps/projects/ilab/people/pls394/plague/ancestry_version1/PCA/output/Lund.mahadist.europe_more3.DF7.202507.csv", sep=",",header=T,row.names=1) %>% t() %>% as.data.frame()
pop_target <- c("Danish","Swedish","Norwegian",
                "German", "Dutch", "Belgian", # Germany and Benelux
                "English","Welsh","Scottish","Orcadian","Northern_Irish","Irish", # Britain and Ireland
                "French", # France
                "Polish","Russian", # East Europe
                "Lithuanian","Estonian", # Baltic
                "Italian", # Italy
                "Finnish", # Uralic
                "Croatian","Greek","Turkish", # Southeast Europe
                "Icelandic" # Iceland)
)


lund.pca <- lund %>% filter( row.names(.) %in% pop_target)
pca_ordered <- lund.pca[match(pop_target, rownames(lund.pca)),match(colnames(plot.mat.lund),colnames(lund.pca)) ] %>% as.matrix()
## rownames(pca_ordered) <- recode(rownames(pca_ordered),
##     "Danish" = "Denmark",
##     "Swedish" = "Sweden",
##     "Norwegian" = "Norway",
##     "Icelandic" = "Iceland", 
##     "Belgian" = "Belgium",
##     "Dutch" = "Netherlands",
##     "German" = "Germany",
##     "French" = "France",
##     "English" = "England",
##     "Welsh" = "Wales",
##     "Scottish" = "Scotland",
##     "Orcadian" = "Orkney",
##     "Northern_Irish" = "Northern Ireland", 
##     "Irish" = "Ireland",
##     "Lithuanian" = "Lithuania",
##     "Estonian" = "Estonia",
##     "Polish" = "Poland",
##     "Russian" = "Russia",
##     "Belarusian" = "Belarus",
##     "Ukrainian" = "Ukraine"
##     )

rownames(pca_ordered) <- recode(rownames(pca_ordered),
    "Danish" = "DK",
    "Swedish" = "SE",
    "Norwegian" = "NO",
    "German" = "DE",
    "Dutch" = "NL",
    "Belgian" = "BE",
    "English" = "GB-ENG",
    "Welsh" = "GB-WLS",
    "Scottish" = "GB-SCT",
    "Orcadian" = "GB-ORK",
    "Northern_Irish" = "GB-NIR",
    "Irish" = "IE",
    "French" = "FR",
    "Polish" = "PL",
    "Russian" = "RU",
    "Lithuanian" = "LT",
    "Estonian" = "EE",
    "Italian" = "IT",
    "Finnish" = "FI",
    "Croatian" = "HR",
    "Greek" = "GR",
    "Turkish" = "TR",
    "Icelandic" = "IS"

    )


## pca color
#main_breaks <- c(0, 0.001, 0.01, 0.05, 0.5)
#main_breaks = c(0,   10,    20,   30,  40)
main_breaks = c(0, qchisq(1-0.05,df=7) %>% sqrt, qchisq(1-0.001,df=7) %>% sqrt, 10, 20, 30)
main_breaks = c(0, qchisq(1-0.05,df=7) %>% sqrt, qchisq(1-0.001,df=7) %>% sqrt, 10, 20, 30)
main_breaks_pconvert = sapply(c( 0.001, 0.01, 0.05, 0.5)   ,function(x) {qchisq(1-x,df=7) %>% sqrt}) %>% rev
main_breaks = c(round(min(pca_ordered),1), main_breaks_pconvert,30)

# Create finer breaks with 5 subdivisions between each main break
fine_breaks <- unlist(lapply(1:(length(main_breaks)-1), function(i) {
  seq(main_breaks[i], main_breaks[i+1], length.out = 6)[-6]  # 5 subdivisions, remove last to avoid duplication
}))
#fine_breaks <- c(fine_breaks, 0.5)
fine_breaks = c(fine_breaks,30)

pca_colors <- scico(length(fine_breaks),palette ='bilbao',begin=0,end=.9)  %>% rev
pca_colors <- scico(length(fine_breaks),palette ='bilbao',begin=.2,end=.95)  ## this one is ok
col_fun <- colorRamp2(
  breaks = fine_breaks,  # Your 26 break points
  colors = pca_colors    # Your 26 colors
)
## column groups
pca_col_groups <- c(rep("Scandinavia", 3), 
                   rep("Benelux", 3),
                   rep("Britian", 6),
                   rep("France",1),
                   rep("EastEurope",2),
                   rep("Baltic",2),
                   rep("Italy",1),
                   rep('Uralic',1),
                   rep('SoutheastEurope',3),
                   rep('Iceland',1)
                   )

pca_col_groups <- factor(pca_col_groups,levels=unique(pca_col_groups))

row_ha.pca = rowAnnotation(PCAsource=rep(NA,nrow(pca_ordered)),col=list(IBDsource =c("non-Local Scandinavian"= "#E9C46A", "non-Scandinavian"="#005F73", "Local Scandinavian"="#E76F51")),
                           show_legend = FALSE,show_annotation_name = FALSE,na_col = "transparent")
spacer_ha <- rowAnnotation(
    spacer = anno_empty(border = FALSE, width = unit(0, "mm"))
)
combined_left_annotation <- c( row_ha.pca,spacer_ha)


pca_lund <- Heatmap(pca_ordered,
        cluster_columns=FALSE,
        cluster_rows=FALSE,
        show_row_names = TRUE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 24),
        row_title = "PCA based MD to",  # Add this line
        row_title_side = "left",            # Position the title
        row_title_gp = gpar(fontsize = 24,fontface = "bold"),
        show_column_names=FALSE,
        column_names_side = "top",
        column_names_rot = 0,

        col = col_fun,
        na_col = 'white',
        column_split = list(annotate.df.row.lund$Age_assignment,annotate.df.row.lund$Group),
        column_gap = unit(c(1.5,1.5,4,1.5,1.5,4,1.5,1.5), "mm"),  # Note: changed from column_gap to row_gap
        column_title = NULL,
        row_split = pca_col_groups,
        row_gap = unit(1.5,'mm'),
        show_heatmap_legend=FALSE,
        ## left_annotation = combined_left_annotation
        ## top_annotation =  combined_top_annotation
        ## width = unit(ncol(plot.mat.lund)*0.8, "mm")
        )

plot(pca_lund)
group_block_anno(7:9, "empty2", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),  # This controls the BOX
                 label = "Post",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))  # This controls the TEXT

group_block_anno(4:6, "empty2", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),
                 label = "During", 
                 label_gp = gpar(fontsize = 22, fontface = "bold"))

group_block_anno(1:3, "empty2", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),
                 label = "Pre",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))

group_block_anno(1:9, "empty1", 
                 gp = gpar(fill = NA, col = NA, lwd = 2),
                 label = "Lund",
                 label_gp = gpar(fontsize = 22, fontface = "bold")) 

############# end of Lund #################


############# Trondheim  #################
plot.mat.2 <- read.table("/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Trondheim.IBD_standard.2026April28.tsv",header = T,row.names=1,sep="\t") %>% as.matrix
info <- read.table("/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Trondheim.IBD_info.2026April28.tsv",header=T,sep="\t")
table(info$Group,info$Age_assignment)
  
#### pre
###pre-nonScand
##pca
pre_nonScand_pca <- info %>% filter(Age_assignment =='Pre' & Group == 'non-Scandinavian' & Method == 'PCA') %>% pull(ID)
pre_nonScand_pca.data = plot.mat.2[,pre_nonScand_pca]
hm.pre_nonScand_pca = heatmap(pre_nonScand_pca.data)
##ibd
pre_nonScand_ibd <- info %>% filter(Age_assignment =='Pre' & Group == 'non-Scandinavian' & Method == 'IBD') %>% pull(ID)
pre_nonScand_ibd.data = plot.mat.2[,pre_nonScand_ibd]
hm.pre_nonScand_ibd = heatmap(pre_nonScand_ibd.data)

###pre-nonLocal
##pca
pre_nonLocal_pca <- info %>% filter(Age_assignment =='Pre' & Group == 'nonlocal Scandinavian'& Method=='PCA') %>% pull(ID)
pre_nonLocal_pca.data = plot.mat.2[,pre_nonLocal_pca]
hm.pre_nonLocal_pca = heatmap(pre_nonLocal_pca.data)
##ibd
pre_nonLocal_ibd <- info %>% filter(Age_assignment =='Pre' & Group == 'nonlocal Scandinavian'& Method=='IBD') %>% pull(ID)
pre_nonLocal_ibd.data = plot.mat.2[,pre_nonLocal_ibd]
hm.pre_nonLocal_ibd = heatmap(pre_nonLocal_ibd.data)

###pre-Local
##pca
pre_Local <- info %>% filter(Age_assignment =='Pre' & Group == 'local Scandinavian' ) %>% pull(ID)
pre_Local.data = plot.mat.2[,pre_Local]
hm.pre_Local = heatmap(pre_Local.data)

col_order.pre <- c(colnames(pre_nonScand_pca.data)[hm.pre_nonScand_pca$colInd],
                   colnames(pre_nonScand_ibd.data)[hm.pre_nonScand_ibd$colInd],
                   colnames(pre_nonLocal_pca.data)[hm.pre_nonLocal_pca$colInd],
                   colnames(pre_nonLocal_ibd.data)[hm.pre_nonLocal_ibd$colInd],
                   colnames(pre_Local.data)[hm.pre_Local$colInd])




#### post
###post-nonScand
##pca
post_nonScand_pca <- info %>% filter(Age_assignment =='Post' & Group == 'non-Scandinavian' & Method == 'PCA') %>% pull(ID)
# no data, skip heatmap
post_nonScand_pca.data = plot.mat.2[,post_nonScand_pca]
##hm.post_nonScand_pca = heatmap(post_nonScand_pca.data)
##ibd
post_nonScand_ibd <- info %>% filter(Age_assignment =='Post' & Group == 'non-Scandinavian' & Method == 'IBD') %>% pull(ID)
#only one sample, so skip heatmap
post_nonScand_ibd.data = plot.mat.2[,post_nonScand_ibd]
##hm.post_nonScand_ibd = heatmap(post_nonScand_ibd.data)

###post-nonLocal
##pca
post_nonLocal_pca <- info %>% filter(Age_assignment =='Post' & Group == 'nonlocal Scandinavian' & Method=='PCA') %>% pull(ID)
post_nonLocal_pca.data = plot.mat.2[,post_nonLocal_pca]
hm.post_nonLocal_pca = heatmap(post_nonLocal_pca.data)
##ibd
post_nonLocal_ibd <- info %>% filter(Age_assignment =='Post' & Group == 'nonlocal Scandinavian'& Method=='IBD') %>% pull(ID)
post_nonLocal_ibd.data = plot.mat.2[,post_nonLocal_ibd]
hm.post_nonLocal_ibd = heatmap(post_nonLocal_ibd.data)

###post-Local
##pca
post_Local <- info %>% filter(Age_assignment =='Post' & Group == 'local Scandinavian' ) %>% pull(ID)
post_Local.data = plot.mat.2[,post_Local]
hm.post_Local = heatmap(post_Local.data)

col_order.post <- c(
                    post_nonScand_ibd,
                   colnames(post_nonLocal_pca.data)[hm.post_nonLocal_pca$colInd],
                   colnames(post_nonLocal_ibd.data)[hm.post_nonLocal_ibd$colInd],
                   colnames(post_Local.data)[hm.post_Local$colInd])



col_order <- c(col_order.pre,col_order.post)
plot.mat.trondheim <- plot.mat.2[,col_order]

### make a heatmap
annotate.df.row.trond <- info[match(col_order,info$ID),]
all.equal(annotate.df.row.trond$id ,col_order)
annotate.df.row.trond$Age_assignment <- factor(annotate.df.row.trond$Age_assignment,levels=unique(annotate.df.row.trond$Age_assignment))
annotate.df.row.trond$Group <- factor(annotate.df.row.trond$Group, levels = unique(annotate.df.row.trond$Group))
annotate.df.row.trond$Pop = 'Trondheim'

annotate.df.col.trond <- data.frame(IBDsource = c(rep("non-Local Scandinavian",13),rep("Local Scandinavian",1),rep("non-Scandinavian",5)))

rownames(plot.mat.trondheim) <- gsub("_",":", rownames(plot.mat.trondheim))
rownames(plot.mat.trondheim) <- c(
  "DK", "SE:Sk", "SE:So", "SE:Gt", "SE:St", "SE:Ce", "SE:No",
  "NO:Os", "NO:So", "NO:Ce", "NO:Ea", "NO:We", "NO:No", "NO:Tr",
  "NL", "DE", "UK:Sh", "UK:Or", "IS"
)

## row annotation
row_ha.left = rowAnnotation(IBDsource=annotate.df.col.trond$IBDsource,col=list(IBDsource =c("non-Local Scandinavian"= "#E9C46A", "non-Scandinavian"="#005F73", "Local Scandinavian"="#E76F51")),
                            show_legend = FALSE,show_annotation_name = FALSE)
spacer_ha <- rowAnnotation(
    spacer = anno_empty(border = FALSE, width = unit(2, "mm"))  # Adjust width for spacing                                                                                                                                                                                                                                                                              
)
combined_left_annotation <- c(row_ha.left , spacer_ha)
## row annotation
col_ha = HeatmapAnnotation(Period=annotate.df.row.trond$Age_assignment, col=list(Period=c("Pre"="blue","During"="red","Post"="green")),na_col='white',show_legend = FALSE,show_annotation_name=FALSE)
col_ha2 =  HeatmapAnnotation(Group1=annotate.df.row.trond$Group, col=list(Group1=c("non-Scandinavian"="#005F73","nonlocal Scandinavian"="#E9C46A",'local Scandinavian'="#E76F51")),na_col='white',show_legend = FALSE,show_annotation_name = FALSE,simple_anno_size  = unit(12, "mm"))
## col_ha3 =  HeatmapAnnotation(Group2=annotate.df.row.trond$Method, col=list(Group2=c("PCA"="beige","IBD"="bisque")),na_col='white',show_legend = FALSE,show_annotation_name = FALSE,simple_anno_size  = unit(12, "mm"))
col_ha3 =  HeatmapAnnotation(Group2=annotate.df.row.trond$Method, col=list(Group2=c("PCA"="#55514D","IBD"="#A39F9A")),na_col='white',show_legend = FALSE,show_annotation_name = FALSE,simple_anno_size  = unit(12, "mm"))


ha1.trond = HeatmapAnnotation(
    empty3 = anno_empty(border = FALSE, width = unit(6, "mm")),
    gp = gpar(fontsize = 16)
)

ha2.trond = HeatmapAnnotation(
    empty4 = anno_empty(border = FALSE, width = unit(6, "mm")),
    gp = gpar(fontsize = 16)
)
spacer1 <- HeatmapAnnotation(
  empty5 = anno_empty(width = unit(22, "mm"),border=FALSE)
)

spacer2 <- HeatmapAnnotation(
  empty6 = anno_empty(width = unit(7, "mm"),border=FALSE)
)

combined_top_annotation <- c(ha1.trond, ha2.trond)
combined_bottom_annotation <- c(spacer1,col_ha2, spacer2,col_ha3)

ibd_trond <- Heatmap(plot.mat.trondheim,
        cluster_columns=FALSE,
        cluster_rows=FALSE,
        show_column_names = FALSE,
        show_row_names=T,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 24),
        ## col = pheatmap_colors,
        col = ibd_color,
        na_col = 'white',
        column_split = list(annotate.df.row.trond$Age_assignment,annotate.df.row.trond$Group),  # Note: changed from column_split to row_split
        column_gap = unit(rev(c(1.5,1.5,4,1.5,1.5)), "mm"),  # Note: changed from column_gap to row_gap
        column_title = NULL,
        left_annotation = combined_left_annotation,

        ##top_annotation = combined_top_annotation,
        bottom_annotation = combined_bottom_annotation,
        show_heatmap_legend = FALSE 
        )

plot(ibd_trond)
group_block_anno = function(group, empty_anno, gp = gpar(), 
    label = NULL, label_gp = gpar()) {

    seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
    loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
    seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
    loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))

    seekViewport("global")
    grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
        just = c("left", "bottom"), gp = gp)
    if(!is.null(label)) {
        grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5, gp = label_gp,rot  = 0)
    }
}



group_block_anno(4:6, "empty4", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),  # This controls the BOX
                 label = "Post",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))  # This controls the TEXT

group_block_anno(1:3, "empty4", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),
                 label = "Pre",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))

group_block_anno(1:6, "empty3", 
                 gp = gpar(fill = NA, col = NA, lwd = 2),
                 label = "Trondheim",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))

### trondheim PCA ###
trondheim = read.table("/maps/projects/ilab/people/pls394/plague/ancestry_version1/PCA/output/Trondheim.mahadist.europe_more3.DF7.202507.csv",sep=",",header=T,row.names=1) %>% t()  %>% as.data.frame()

pop_target <- c("Danish","Swedish","Norwegian",
                "German", "Dutch", "Belgian", # Germany and Benelux
                "English","Welsh","Scottish","Orcadian","Northern_Irish","Irish", # Britain and Ireland
                "French", # France
                "Polish","Russian", # East Europe
                "Lithuanian","Estonian", # Baltic
                "Italian", # Italy
                "Finnish", # Uralic
                "Croatian","Greek","Turkish", # Southeast Europe
                "Icelandic" # Iceland)
                )

trondheim.pca <- trondheim %>% filter( row.names(.) %in% pop_target)
pca_ordered <- trondheim.pca[match(pop_target, rownames(trondheim.pca)),match(colnames(plot.mat.trondheim),colnames(trondheim.pca)) ] %>% as.matrix()
rownames(pca_ordered) <- recode(rownames(pca_ordered),
    "Danish" = "DK",
    "Swedish" = "SE",
    "Norwegian" = "NO",
    "German" = "DE",
    "Dutch" = "NL",
    "Belgian" = "BE",
    "English" = "GB-ENG",
    "Welsh" = "GB-WLS",
    "Scottish" = "GB-SCT",
    "Orcadian" = "GB-ORK",
    "Northern_Irish" = "GB-NIR",
    "Irish" = "IE",
    "French" = "FR",
    "Polish" = "PL",
    "Russian" = "RU",
    "Lithuanian" = "LT",
    "Estonian" = "EE",
    "Italian" = "IT",
    "Finnish" = "FI",
    "Croatian" = "HR",
    "Greek" = "GR",
    "Turkish" = "TR",
    "Icelandic" = "IS"
    )

## pca color
#main_breaks <- c(0, 0.001, 0.01, 0.05, 0.5)
#main_breaks = c(0,   10,    20,   30,  40)
main_breaks = c(0,  qchisq(1-0.05,df=7) %>% sqrt,qchisq(1-0.001,df=7) %>% sqrt, 10, 20, 30)
main_breaks_pconvert = sapply(c( 0.001, 0.01, 0.05, 0.5)   ,function(x) {qchisq(1-x,df=7) %>% sqrt}) %>% rev
main_breaks = c(1, main_breaks_pconvert,30)

# Create finer breaks with 5 subdivisions between each main break
fine_breaks <- unlist(lapply(1:(length(main_breaks)-1), function(i) {
  seq(main_breaks[i], main_breaks[i+1], length.out = 6)[-6]  # 5 subdivisions, remove last to avoid duplication
}))
#fine_breaks <- c(fine_breaks, 0.5)
fine_breaks = c(fine_breaks,30)

pca_colors <- scico(length(fine_breaks),palette ='bilbao',begin=0,end=.9)  %>% rev
pca_colors <- scico(length(fine_breaks),palette ='bilbao',begin=.2,end=.95)  ## this one is ok
pca_colors <- scico(length(fine_breaks),palette ='bilbao',begin=.2,end=.95)  ## this one is ok

col_fun <- colorRamp2(
  breaks = fine_breaks,  # Your 26 break points
  colors = pca_colors    # Your 26 colors
)

## column groups
pca_col_groups <- c(rep("Scandinavia", 3), 
                   rep("Benelux", 3),
                   rep("Britian", 6),
                   rep("France",1),
                   rep("EastEurope",2),
                   rep("Baltic",2),
                   rep("Italy",1),
                   rep('Uralic',1),
                   rep('SoutheastEurope',3),
                   rep('Iceland',1)
                   )
pca_col_groups <- factor(pca_col_groups,levels=unique(pca_col_groups))
row_ha.pca = rowAnnotation(PCAsource=rep(NA,nrow(pca_ordered)),col=list(IBDsource =c("non-Local Scandinavian"= "#E9C46A", "non-Scandinavian"="#005F73", "Local Scandinavian"="#E76F51")),
                           show_legend = FALSE,show_annotation_name = FALSE,na_col = "transparent")
spacer_ha <- rowAnnotation(
    spacer = anno_empty(border = FALSE, width = unit(0, "mm"))
    )
combined_left_annotation <- c(row_ha.pca,spacer_ha)

pca_trondheim <- Heatmap(pca_ordered,
        cluster_columns=FALSE,
        cluster_rows=FALSE,
        show_row_names = TRUE,
        show_column_names=FALSE,
        row_names_side = "left",
        row_names_rot = 0,
        row_names_gp = gpar(fontsize = 24),

        col = col_fun,
        na_col = 'white',
        column_split = list(annotate.df.row.trond$Age_assignment,annotate.df.row.trond$Group),  # Note: changed from column_split to row_split
        column_gap = unit(rev(c(1.5,1.5,4,1.5,1.5)), "mm"),  # Note: changed from column_gap to row_gap 
        column_title = NULL,
        row_split = pca_col_groups,
        row_title =NULL,
        row_gap = unit(1.5,'mm'),
        show_heatmap_legend=FALSE,
        ## left_annotation = combined_left_annotation,
        ## top_annotation = combined_top_annotation  # Adjust annotations
        )

plot(pca_trondheim)
group_block_anno(4:6, "empty4", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),  # This controls the BOX
                 label = "Post",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))  # This controls the TEXT

group_block_anno(1:3, "empty4", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),
                 label = "Pre",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))

group_block_anno(1:6, "empty3", 
                 gp = gpar(fill = NA, col = NA, lwd = 2),
                 label = "Trondheim",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))
############# end of Trondheim  #################



############# Vilnius  #################
##plot.mat.2 <- read.table("/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Vilnius.IBD_standard.2026April28.tsv",header = T,row.names=1,sep="\t") %>% as.matrix
plot.mat.2 <- read.table("/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Vilnius.IBD_standard.2026May28.tsv",header = T,row.names=1,sep="\t") %>% as.matrix

##rownames(plot.mat.2)[2] <- "Estonia/Latvia"
##rownames(plot.mat.2)[4] <- "Bela/Rus"
info <- read.table("/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Vilnius.IBD_info.2026April.tsv",header=T,sep="\t")
info$Age_assignment[info$Age_assignment=='During'] = 'Post'
table(info$Group,info$Age_assignment)
  
#### pre
###pre-nonBaltic
##pca
pre_nonBaltic_pca <- info %>% filter(Age_assignment =='Pre' & Group == 'non-Baltic' & Method == 'PCA') %>% pull(ID)
pre_nonBaltic_pca.data = plot.mat.2[,pre_nonBaltic_pca]
hm.pre_nonBaltic_pca = heatmap(pre_nonBaltic_pca.data)
##ibd
pre_nonBaltic_ibd <- info %>% filter(Age_assignment =='Pre' & Group == 'non-Baltic' & Method == 'IBD') %>% pull(ID)
pre_nonBaltic_ibd.data = plot.mat.2[,pre_nonBaltic_ibd]
hm.pre_nonBaltic_ibd = heatmap(pre_nonBaltic_ibd.data)


###pre-Baltic
pre_Local <- info %>% filter(Age_assignment =='Pre' & Group == 'Baltic' ) %>% pull(ID)
pre_Local.data = plot.mat.2[,pre_Local]
hm.pre_Local = heatmap(pre_Local.data)

col_order.pre <- c(colnames(pre_nonBaltic_pca.data)[hm.pre_nonBaltic_pca$colInd],
                   colnames(pre_nonBaltic_ibd.data)[hm.pre_nonBaltic_ibd$colInd],
                   colnames(pre_Local.data)[hm.pre_Local$colInd])


#### post
###post-nonBaltic
##pca
post_nonBaltic_pca <- info %>% filter(Age_assignment =='Post' & Group == 'non-Baltic' & Method== 'PCA') %>% pull(ID)
post_nonBaltic_pca.data = plot.mat.2[,post_nonBaltic_pca]
hm.post_nonBaltic_pca = heatmap(post_nonBaltic_pca.data)
##ibd
post_nonBaltic_ibd <- info %>% filter(Age_assignment =='Post' & Group == 'non-Baltic' & Method == 'IBD') %>% pull(ID)
## no sample, so skip hm
post_nonBaltic_ibd.data = plot.mat.2[,post_nonBaltic_ibd]
##hm.post_nonBaltic_ibd = heatmap(post_nonBaltic_ibd.data)


###post-Baltic
##pca
post_Local <- info %>% filter(Age_assignment =='Post' & Group == 'Baltic' ) %>% pull(ID)
post_Local.data = plot.mat.2[,post_Local]
hm.post_Local = heatmap(post_Local.data)

col_order.post <- c(post_nonBaltic_pca,

                   colnames(post_Local.data)[hm.post_Local$colInd])


col_order <- c(col_order.pre,col_order.post)

### make heatmap
plot.mat.vilnius <- plot.mat.2[,col_order]

annotate.df.row.vilnius <- info[match(col_order,info$ID),]
all.equal(annotate.df.row.vilnius$id ,col_order)
annotate.df.row.vilnius$Age_assignment <- factor(annotate.df.row.vilnius$Age_assignment,levels=unique(annotate.df.row.vilnius$Age_assignment))
annotate.df.row.vilnius$Group <- factor(annotate.df.row.vilnius$Group, levels = unique(annotate.df.row.vilnius$Group))
annotate.df.row.vilnius$Pop = 'Vilnius'

annotate.df.col.vilnius <- data.frame(IBDsource = c(rep("Baltic",1),rep("non-Baltic",8)))

na.matrix <- matrix(NA,ncol=57,nrow=10)
plot.mat.vilnius2<- rbind(plot.mat.vilnius,na.matrix)
rownames(plot.mat.vilnius2) <- c(
    "LT",        # Lithuania
    "EE/LV",    # Estonia/Latvia
    "PL",        # Poland
    "RU/BY",        # Bela/Rus (Belarus)
    "UA",        # Ukraine
    "DK",        # Denmark
    "SE",        # Sweden
    "NO",        # Norway
    "NL",        # Netherlands
    "DE",        # Germany
    "", "", "", "", "", "", "", "", ""  # Keep empty strings as is
)

rownames(plot.mat.vilnius2) <- c(
    "Baltic",        # Lithuania
    #"EE/LV",    # Estonia/Latvia
    "PL",        # Poland
    #"RU/BY/UA",        # Bela/Rus (Belarus)
    "E Slavic",
    "Se EUR",
    "DK",        # Denmark
    "SE",        # Sweden
    "NO",        # Norway
    "NL",        # Netherlands
    "DE",        # Germany
    "", "", "", "", "", "", "", "", "" ,"" # Keep empty strings as is
)


## row annotation
row_ha.left = rowAnnotation(IBDsource=c(annotate.df.col.vilnius$IBDsource,rep('Fill',10)),col=list(IBDsource =c("non-Baltic"="#005F73", "Baltic"="#E76F51","Fill"="white")),
                               show_legend = FALSE,show_annotation_name = FALSE,
                               na_col=NULL)
spacer_ha <- rowAnnotation(
    spacer = anno_empty(border = FALSE, width = unit(2, "mm"))  # Adjust width for spacing
)
combined_left_annotation <- c(row_ha.left , spacer_ha)
## column annotation
col_ha = HeatmapAnnotation(Period=annotate.df.row.vilnius$Age_assignment, col=list(Period=c("Pre"="blue","During"="red","Post"="green")),na_col='white',show_legend = FALSE,show_annotation_name=FALSE)

## old col_ha2
col_ha2 =  HeatmapAnnotation("Classification"=annotate.df.row.vilnius$Group, col=list("Classification"=c("non-Baltic"="#005F73","nonlocal Scandinavian"="#E9C46A",'Baltic'="#E76F51")),na_col='white',show_legend = FALSE,show_annotation_name = FALSE,simple_anno_size  = unit(12, "mm"))


## new col_ha2,stripped pattern, does not work out
## grp = annotate.df.row.vilnius$Group

## classification_graphics = list(
##   "non-Baltic" = function(x, y, w, h) {
##     grid.rect(x, y, w, h, gp = gpar(fill = "#005F73", col = NA))
##   },
##   "nonlocal Scandinavian" = function(x, y, w, h) {
##     grid.rect(x, y, w, h, gp = gpar(fill = "#E9C46A", col = NA))
##   },
##   "Baltic" = function(x, y, w, h) {
##     # base fill
##     grid.rect(x, y, w, h, gp = gpar(fill = "white", col = NA))

##     # diagonal stripes
##     for (off in seq(-1, 1, by = 0.12)) {
##       grid.lines(
##         x = unit(c(off, off + 1), "npc"),
##         y = unit(c(0, 1), "npc"),
##         gp = gpar(col = "#E9C46A", lwd = 2)
##       )
##     }

##     # optional second color between the gold stripes
##     for (off in seq(-0.94, 1.06, by = 0.12)) {
##       grid.lines(
##         x = unit(c(off, off + 1), "npc"),
##         y = unit(c(0, 1), "npc"),
##         gp = gpar(col = "#E76F51", lwd = 2)
##       )
##     }
##   }
## )

## col_ha2 = HeatmapAnnotation(
##   Classification = anno_customize(grp, graphics = classification_graphics),
##   show_legend = FALSE,
##   show_annotation_name = FALSE,
##   simple_anno_size = unit(12, "mm")
## )



### vertical split
grp = annotate.df.row.vilnius$Group
classification_graphics = list(
  "non-Baltic" = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "#005F73", col = NA))
  },
  "nonlocal Scandinavian" = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "#E9C46A", col = NA))
  },
  "Baltic" = function(x, y, w, h) {
    # top 1/6
    grid.rect(
      x = x, y = y + h/4,
      width = w, height = h/2,
      gp = gpar(fill = "#E9C46A", col = NA)
    )
    # bottom half
    grid.rect(
      x = x, y = y - h/4,
      width = w, height = h/2,
      gp = gpar(fill = "#E76F51", col = NA)
    )
  }
)

classification_graphics = list(
  "non-Baltic" = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "#005F73", col = NA))
  },
  
  "nonlocal Scandinavian" = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "#E9C46A", col = NA))
  },
  
  "Baltic" = function(x, y, w, h) {
    n = 4
    part_h = h / n
    
    # Example: alternate colors (you can customize this vector)
    fills = c("#E9C46A", "#E76F51", "#E9C46A", "#E76F51", "#E9C46A", "#E76F51")
    
    for (i in seq_len(n)) {
      grid.rect(
        x = x,
        y = y - h/2 + (i - 0.5) * part_h,
        width = w,
        height = part_h,
        gp = gpar(fill = fills[i], col = NA)
      )
    }
  }
)

col_ha2 = HeatmapAnnotation(
  Classification = anno_customize(grp, graphics = classification_graphics),
  gp = gpar(col = NA), 
  border = FALSE,      
  show_legend = FALSE,
  show_annotation_name = FALSE,
  simple_anno_size = unit(12, "mm")
)

col_ha2 = HeatmapAnnotation(
  Classification = anno_customize(
    grp,
    graphics = classification_graphics
  ),
  height = unit(12, "mm"),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  border = FALSE
)


## col_ha3 =  HeatmapAnnotation(Method=annotate.df.row.vilnius$Method, col=list(Method=c("PCA"="beige","IBD"="bisque")),na_col='white',show_legend = FALSE,show_annotation_name = F,simple_anno_size  = unit(12, "mm"))
col_ha3 =  HeatmapAnnotation(Method=annotate.df.row.vilnius$Method, col=list(Method=c("PCA"="#55514D","IBD"="#A39F9A")),na_col='white',show_legend = FALSE,show_annotation_name = F,simple_anno_size  = unit(12, "mm"))

ha1 = HeatmapAnnotation(
    empty1 = anno_empty(border = FALSE, width = unit(6, "mm")),
    gp = gpar(fontsize = 16)
)

ha2 = HeatmapAnnotation(
    empty2 = anno_empty(border = FALSE, width = unit(6, "mm")),
    gp = gpar(fontsize = 16)
)
spacer1 <- HeatmapAnnotation(
  empty3 = anno_empty(width = unit(22, "mm"),border=FALSE)
)
spacer2 <- HeatmapAnnotation(
  empty4 = anno_empty(width = unit(7, "mm"),border=FALSE)
)

combined_top_annotation <- c(ha1, ha2)
combined_bottom_annotation <- c(spacer1,col_ha2, spacer2,col_ha3)

ibd_vilnius <- Heatmap(plot.mat.vilnius2,
        name='Vilnius_IBD',
        cluster_columns=FALSE,
        cluster_rows=FALSE,
        show_row_names = TRUE,
        show_column_names=FALSE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 24),

        ## col = pheatmap_colors,
        col = ibd_color,
        na_col = "white",
        column_split = list(annotate.df.row.vilnius$Age_assignment,annotate.df.row.vilnius$Group),  # Note: changed from column_split to row_split
        ## column_gap = unit(c(1.5,4,1.5,4,1.5), "mm"),  # Note: changed from column_gap to row_gap
        column_gap = unit(c(1.5,4,1.5), "mm"),  # Note: changed from column_gap to row_gap
        column_title = NULL,
        left_annotation = combined_left_annotation,  # Adjust annotations accordingly
        ## top_annotation = combined_top_annotation,
        bottom_annotation = combined_bottom_annotation,
        ##heatmap_legend_param = list(title = NULL)  ,
        show_heatmap_legend = FALSE 
        )

plot(ibd_vilnius)
group_block_anno = function(group, empty_anno, gp = gpar(), 
    label = NULL, label_gp = gpar()) {

    seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
    loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
    seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
    loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))

    seekViewport("global")
    grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
        just = c("center", "center"), gp = gp)
    if(!is.null(label)) {
        grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5, gp = label_gp,rot  = 90)
    }
}



group_block_anno(5:6, "empty2", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),  # This controls the BOX
                 label = "Post",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))  # This controls the TEXT

group_block_anno(3:4, "empty2", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),
                 label = "During", 
                 label_gp = gpar(fontsize = 22, fontface = "bold"))

group_block_anno(1:2, "empty2", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),
                 label = "Pre",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))

group_block_anno(1:6, "empty1", 
                 gp = gpar(fill = NA, col = NA, lwd = 2),
                 label = "Vilnius",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))

### vilnius PCA ###
vilnius = read.table("/maps/projects/ilab/people/pls394/plague/ancestry_version1/PCA/output/Lithuania.mahadist.europe_more3.DF7.202507.csv",sep=",",header=T,row.names=1) %>% t() %>% as.data.frame()

pop_target <- c("Danish","Swedish","Norwegian",
                "German", "Dutch", "Belgian", # Germany and Benelux
                "English","Welsh","Scottish","Orcadian","Northern_Irish","Irish", # Britain and Ireland
                "French", # France
                "Polish","Russian", # East Europe
                "Lithuanian","Estonian", # Baltic
                "Italian", # Italy
                "Finnish", # Uralic
                "Croatian","Greek","Turkish", # Southeast Europe
                "Icelandic" # Iceland)
)

vilnius.pca <- vilnius %>% filter( row.names(.) %in% pop_target)

pca_ordered <- vilnius.pca[match(pop_target, rownames(vilnius.pca)),match(colnames(plot.mat.vilnius),colnames(vilnius.pca))] %>% as.matrix()

rownames(pca_ordered) <- recode(rownames(pca_ordered),
    "Danish" = "DK",
    "Swedish" = "SE",
    "Norwegian" = "NO",
    "German" = "DE",
    "Dutch" = "NL",
    "Belgian" = "BE",
    "English" = "GB-ENG",
    "Welsh" = "GB-WLS",
    "Scottish" = "GB-SCT",
    "Orcadian" = "GB-ORK",
    "Northern_Irish" = "GB-NIR",
    "Irish" = "IE",
    "French" = "FR",
    "Polish" = "PL",
    "Russian" = "RU",
    "Lithuanian" = "LT",
    "Estonian" = "EE",
    "Italian" = "IT",
    "Finnish" = "FI",
    "Croatian" = "HR",
    "Greek" = "GR",
    "Turkish" = "TR",
    "Icelandic" = "IS"

    )

## pca color
main_breaks <- c(0, 0.001, 0.01, 0.05, 0.5)
#main_breaks = c(0,   10,    20,   30,  40)
main_breaks = c(0,  qchisq(1-0.05,df=7) %>% sqrt,qchisq(1-0.001,df=7) %>% sqrt, 10, 20, 30)
#main_breaks = c(0 ,qchisq(1-0.001,df=7) %>% sqrt, 10, 20, 30)
main_breaks_pconvert = sapply(c( 0.001, 0.01, 0.05, 0.5)   ,function(x) {qchisq(1-x,df=7) %>% sqrt}) %>% rev
main_breaks = c(1, main_breaks_pconvert,30)

# Create finer breaks with 5 subdivisions between each main break
fine_breaks <- unlist(lapply(1:(length(main_breaks)-1), function(i) {
  seq(main_breaks[i], main_breaks[i+1], length.out = 6)[-6]  # 5 subdivisions, remove last to avoid duplication
}))
#fine_breaks <- c(fine_breaks, 0.5)
fine_breaks = c(fine_breaks,30)

#pca_colors <- scico(length(fine_breaks),palette ='hawaii')
#pca_colors <- scico(length(fine_breaks),palette ='davos')  %>% rev

pca_colors <- scico(length(fine_breaks),palette ='bilbao',begin=.2,end=.95)  ## this one is ok

col_fun = colorRamp2(breaks=fine_breaks,
                     colors = pca_colors)
## column groups
pca_col_groups <- c(rep("Scandinavia", 3), 
                   rep("Benelux", 3),
                   rep("Britian", 6),
                   rep("France",1),
                   rep("EastEurope",2),
                   rep("Baltic",2),
                   rep("Italy",1),
                   rep('Uralic',1),
                   rep('SoutheastEurope',3),
                   rep('Iceland',1)
                   )

pca_col_groups <- factor(pca_col_groups,levels=unique(pca_col_groups))
row_ha.pca = rowAnnotation(PCAsource=rep(NA,nrow(pca_ordered)),col=list(IBDsource =c("non-Local Scandinavian"= "#E9C46A", "non-Scandinavian"="#005F73", "Local Scandinavian"="#E76F51")),
                           show_legend = FALSE,show_annotation_name = FALSE,na_col = "transparent")

spacer_ha <- rowAnnotation(
    spacer = anno_empty(border = FALSE, width = unit(0, "mm"))  # Adjust width for spacing
    )
combined_left_annotation <- c(row_ha.pca,spacer_ha)
pca_vilnius <- Heatmap(pca_ordered,
        cluster_columns=FALSE,
        cluster_rows=FALSE,
        show_row_names = TRUE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 24),
        show_column_names=FALSE,
        column_names_side = "top",
        column_names_rot = 0,

        col = col_fun,
        na_col = 'white',
        column_split = list(annotate.df.row.vilnius$Age_assignment,annotate.df.row.vilnius$Group),
        ## column_gap = unit(c(1.5,4,1.5,4,1.5), "mm"),  # Note: changed from column_gap to row_gap
        column_gap = unit(c(1.5,4,1.5), "mm"),  # Note: changed from column_gap to row_gap
        column_title = NULL,
        row_split = pca_col_groups,
        row_title =NULL,
        row_gap = unit(1.5,'mm'),
        show_heatmap_legend=FALSE,
        ##left_annotation = combined_left_annotation
        ## top_annotation =  combined_top_annotation
        )

plot(pca_vilnius)
group_block_anno(5:6, "empty2", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),  # This controls the BOX
                 label = "Post",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))  # This controls the TEXT

group_block_anno(3:4, "empty2", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),
                 label = "During", 
                 label_gp = gpar(fontsize = 22, fontface = "bold"))

group_block_anno(1:2, "empty2", 
                 gp = gpar(fill = NA, col = "black", lwd = 2),
                 label = "Pre",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))

group_block_anno(1:6, "empty1", 
                 gp = gpar(fill = NA, col = NA, lwd = 2),
                 label = "Vilnius",
                 label_gp = gpar(fontsize = 22, fontface = "bold"))


############# end of Vilnius #################


p1.pca = grid.grabExpr(draw(pca_lund))
p1.ibd = grid.grabExpr(draw(ibd_lund))
## p1.ibd = grid.grabExpr({
##     draw(ibd_lund)
##     decorate_heatmap_body("IBD_Lund", {
##         # Box around specific region (adjust coordinates as needed)
##         grid.rect(x = unit(0.2, "npc"), y = unit(0.3, "npc"),
##                  width = unit(0.4, "npc"), height = unit(0.4, "npc"),
##                  gp = gpar(fill = NA, col = "blue", lwd = 2))
##     })
## })
p2.pca = grid.grabExpr(draw(pca_trondheim))
p2.ibd = grid.grabExpr(draw(ibd_trond))
p3.pca = grid.grabExpr(draw(pca_vilnius))
p3.ibd = grid.grabExpr(draw(ibd_vilnius))

##bitmap('/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/PCA_IBD.heatmap.2026Mar12.bilbao.col_fun.png',w=40,h=20,res=300)
fig = cowplot::plot_grid(
             p1.pca, p2.pca, p3.pca,
             NULL, NULL, NULL,
             p1.ibd ,p2.ibd, p3.ibd,
             ncol = 3,
             align = 'hv',  # 'h' for horizontal, 'v' for vertical alignment
             rel_widths = c(ncol(plot.mat.lund),ncol(plot.mat.trondheim)*1.5, ncol(plot.mat.vilnius2)*2),
             rel_heights =c(1, 0.05,1.2)
             
             )
ggsave('/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/PCA_IBD.heatmap.2026Jun03.bilbao.col_fun.svg', fig, width = 40, height = 20,unit='in')




