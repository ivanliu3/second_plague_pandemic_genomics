library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(RColorBrewer)
library(openxlsx)
library(ComplexHeatmap)

source("/maps/projects/ilab/people/pls394/plague/ancestry_version1/script/t-test.R")
standardize_rows <- function(data) {
  # Check if input is a data frame or matrix
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input must be a data frame or matrix")
  }
  
  # Convert to matrix for calculations
  mat <- as.matrix(data)
  
  # Calculate row means and standard deviations
  row_means <- apply(mat, 1, mean, na.rm = TRUE)
  row_sds <- apply(mat, 1, sd, na.rm = TRUE)
  
  # Handle cases where SD is 0 (to avoid division by 0)
  row_sds[row_sds == 0] <- 1
  
  # Standardize the data (subtract mean, divide by SD)
  standardized <- (mat - row_means) / row_sds
  
  # Convert back to original format
  if (is.data.frame(data)) {
    standardized <- as.data.frame(standardized)
    colnames(standardized) <- colnames(data)
    rownames(standardized) <- rownames(data)
  }
  
  return(standardized)
}


masterFile <- "https://docs.google.com/spreadsheets/d/1V42OBU_UWXDCuF4P_75L8-GaSSlBh0pGvts0hVexTBY/edit?gid=350344439#gid=350344439"
master <- as.data.frame(gsheet::gsheet2tbl(masterFile))
samples.df <-master %>% filter(`Sampling Site` == 'Vilnius' & `Included in Analyses?` =='yes' )  %>% select(Identity,`Cohort Assignment`)
sample.pre <- samples.df %>% filter(`Cohort Assignment` == 'pre-BD') %>% pull(Identity)
sample.bet <- samples.df %>% filter(`Cohort Assignment` == 'during-BD') %>% pull(Identity)
sample.post <- samples.df %>% filter(`Cohort Assignment` == 'post-BD') %>% pull(Identity)

samples <- samples.df %>% pull(Identity)


#myCol <- rev(brewer.pal(4, "Set3"))
myCol <- c("#F8766D", "#7CAE00" ,"#C77CFF")
lit <- read.table("/maps/projects/ilab/people/pls394/plague/ancestry_version1/PCA/output/Lithuania.pvalue.europe_more3.DF7.202507.csv",header=T,sep=",",row.names=1)
p.thres <- 0.001

##### Europe ####
### step1: PCA-outlier under P-value of 0.001
## # "English"        "Irish"          "Northern_Irish" "Scottish"       "Welsh"          "Orcadian"
outlier.pca <- lit %>% filter(row.names(lit) %in% samples) %>%
    ##    filter(Lithuanian < p.thres & Estonian < p.thres  & Latvian <p.thres) %>% rownames() # n = 9, there are only 7 Latvian, the same outliers found if ignore Latvian
    ##filter(Lithuanian < p.thres) %>% rownames() # n = 10, there are only 7 Latvian, the same outliers found if ignore Latvian
     filter(Lithuanian < p.thres & Estonian < p.thres) %>% rownames() # n = 10, there are only 7 Latvian, the same outliers found if ignore Latvian
length(outlier.pca) # n = 9
nonoutlier.pca <- samples[! samples %in% outlier.pca] #n=48

### step2: IBD outliers
ibd_folder = '/projects/ilab/people/pls394/plague/data_decode_2024Nov/IBD/UK_biobank_2026Fig/'
t<-read.table(paste0(ibd_folder,"IBD_ALL_v_xbi_and_nonxbi_regrouped.masked.7cM.100maskedLociPerOrigCm.country_decaf_plus_changed_xbi.summary_to_send.txt"),header=T)
nonxbi <- t %>% filter(totally_excluded != TRUE ) 

locals <- c("Baltic")
nonlocals <- c("Poland","Belarus_Russia_Ukraine","Southeast_Europe","Denmark","Sweden","Norway","Netherlands","Germany")

outlier.ibd.new <- c()
#all_ibdtest <- data.frame()
for (sample in nonoutlier.pca) {
    df.sub = nonxbi %>% filter(ID == sample)
    
    df.nonBaltic = df.sub %>% filter(ID2Grp %in% nonlocals)
    #df.nonScan.top = df.nonScan[(which.max(df.nonScan$cm_mean)),]
    df.nonBaltic.top = df.nonBaltic[which.max(df.nonBaltic$shareID2Mean),]
  
    df.Baltic = df.sub %>% filter( ID2Grp %in% locals)
    df.Baltic.top = df.Baltic[which.max(df.Baltic$shareID2Mean),]

    t_output <- t.test2(
        m1 = df.nonBaltic.top$shareID2Mean,
        m2 = df.Baltic.top$shareID2Mean,
        s1 = df.nonBaltic.top$shareID2SD,
        s2 = df.Baltic.top$shareID2SD,
        n1 = df.nonBaltic.top$ID2Cnt,
        n2 = df.Baltic.top$ID2Cnt,
        m = 0,
        equal.variance = FALSE,
        alternative = "greater"
    )

    if (t_output['p-value'] < 0.05) {
        outlier.ibd.new <- append(outlier.ibd.new,sample)
  }
}
outlier.ibd.new #n=2





outlier.level1 <- c(outlier.pca,outlier.ibd.new) # n = 11
nonoutlier.level1 <- samples[!samples %in% outlier.level1]
print(outlier.level1)


nonLocal.ibd.sample <- outlier.level1
nonLocal.type = data.frame(ID=outlier.level1,Type=c(rep('PCA',length(outlier.pca)), rep('IBD',length(outlier.ibd.new))))
#write.table(outlier.level1,"/projects/ilab/people/pls394/plague/ancestry_version2/output/Cambridgeshire.outlier.finePCADF7-IBD.level1.202508.txt",col.names = F,row.names=F,quote=FALSE,sep='\t')

Local.ibd.sample<- samples.df$Identity[! samples.df$Identity  %in% outlier.level1]
Local.type = data.frame(ID=Local.ibd.sample,Type=c(rep(NA,length(Local.ibd.sample))))


### prepare matrix of standardized IBD and classification table
ancestry.vilnius.info =rbind(nonLocal.type,Local.type)
names(ancestry.vilnius.info) = c('ID','Method')
ancestry.vilnius.info$Group = c(rep("non-Baltic",nrow(nonLocal.type)),rep("Baltic", nrow(Local.type)))
ancestry.vilnius.info = ancestry.vilnius.info %>% left_join(samples.df, by=c('ID'='Identity'))


colOrder <-c("Lithuania", "Estonia_and_Latvia","Poland","Belarus_Russia" ,"Ukraine","Denmark","Sweden","Norway","Netherlands","Germany")
colOrder <- c("Baltic","Poland","Belarus_Russia_Ukraine","Southeast_Europe","Denmark","Sweden","Norway","Netherlands","Germany")

ibd.matrix =    nonxbi %>%
    filter(ID2Grp %in% colOrder) %>% 
    filter( ID %in% ancestry.vilnius.info$ID) %>%
    select(ID,shareID2Mean, ID2Grp) %>%
    pivot_wider(names_from = ID2Grp, values_from = shareID2Mean)%>%
    as.data.frame()

ibd.matrix.std = standardize_rows(ibd.matrix[,-1])
rownames(ibd.matrix.std) <- ibd.matrix$ID
ibd.matrix.std.m = as.matrix(ibd.matrix.std)[,colOrder,drop=FALSE]
ibd.matrix.std.transpos.m <- t(ibd.matrix.std.m)


##write.table(ibd.matrix.std.transpos.m,file='/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Vilnius.IBD_standard.2026April28.tsv',sep='\t',col.names = T,row.names=T,quote=FALSE)
write.table(ibd.matrix.std.transpos.m,file='/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Vilnius.IBD_standard.2026May28.tsv',sep='\t',col.names = T,row.names=T,quote=FALSE)
write.table(ancestry.vilnius.info, file='/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Vilnius.IBD_info.2026April.tsv',sep='\t',col.names = T,row.names=F,quote=FALSE)








nonxbi2 = nonxbi %>% filter(ID %in% samples.df$Identity) %>%
    filter(ID2Grp %in% c(locals,nonlocals))

max_group_by_ID <- nonxbi2 %>%
  group_by(ID) %>%
  slice_max(order_by = shareID2Mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(ID, ID2Grp, shareID2Mean)

max_group_by_ID

t2 = max_group_by_ID %>% left_join(samples.df, by=c('ID'='Identity'))
table(t2$`Cohort Assignment`,t2$ID2Grp)

tab = matrix(c(14, 33, 5, 5), nrow = 2)
rownames(tab) <- c("Group1", "Group2")
colnames(tab) <- c("Category1", "Category2")
tab
fisher.test(tab)
