library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(VennDiagram)
library(RColorBrewer)
library(openxlsx)

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
samples.df <-master %>% filter(sampling_site == 'Trondheim' & included_analysis =='yes' )  %>% select(id,Age_assignment)
sample.pre <- samples.df %>% filter(Age_assignment == 'Pre') %>% pull(id)
sample.post <- samples.df %>% filter(Age_assignment == 'Post') %>% pull(id)
samples <- samples.df %>% pull(id)

myCol <- c("#F8766D", "#7CAE00" ,"#C77CFF")
trond <- read.table("/maps/projects/ilab/people/pls394/plague/ancestry_version1/PCA/output/Trondheim.pvalue.europe_more3.DF7.202507.csv",header=T,sep=",",row.names=1)
p.thres <- 0.001


##### Europe ####
### step1: PCA-outlier under P-value of 0.001
outlier.pca <- trond  %>% filter(row.names(trond) %in% samples) %>% filter(Danish < p.thres & Swedish < p.thres & Norwegian < p.thres) %>% rownames() # n=7
nonoutlier.pca <- samples[! samples %in% outlier.pca] #n=128



### step2: IBD outliers
ibdsummary <-read.table("/maps/projects/ilab/people/pls394/plague/data_decode_2024Nov/IBD/omni_2026Fig/IBD_ALL_v_omni.masked.6cM.100maskedLociPerOrigCm.summary_to_send.txt",header=T)
ibdsummary$ID[which(ibdsummary$ID == "AO48")] <-"AØ48"
outlier.ibd.new <- c()
#all_ibdtest <- data.frame()
for (sample in nonoutlier.pca) {
  df.sub = ibdsummary %>% filter(ID == sample) %>% filter(! ID2Grp  %in% c("British","Irish")) # Icelandic, British, Irish are only used in figure for comparison
                                                                                              # in this version, Icelandic is not used for identifying non-Scand outliers, but used for identifying non-local Scand.
  df.nonScan = df.sub %>% filter(ID2Grp %in% c("Netherlands","Germany","Orkney","Shetland","Iceland"))
  #df.nonScan.top = df.nonScan[(which.max(df.nonScan$cm_mean)),]
  df.nonScan.top = df.nonScan[which.max(df.nonScan$shareID2Mean),]
  
  df.Scan = df.sub %>% filter(! ID2Grp %in% c("Netherlands","Germany","Orkney","Shetland","Iceland"))
  df.Scan.top = df.Scan[which.max(df.Scan$shareID2Mean),]

  t_output <- t.test2(
      m1 = df.nonScan.top$shareID2Mean,
      m2 = df.Scan.top$shareID2Mean,
      s1 = df.nonScan.top$shareID2SD,
      s2 = df.Scan.top$shareID2SD,
      n1 = df.nonScan.top$ID2Cnt,
      n2 = df.Scan.top$ID2Cnt,
      m = 0,
      equal.variance = FALSE,
      alternative = "greater"
      )

  if (t_output['p-value'] < 0.05) {
      outlier.ibd.new <- append(outlier.ibd.new,sample)
  }
}
outlier.ibd.new #n=10

outlier.level1 <- c(outlier.pca,outlier.ibd.new) # n = 17
print(outlier.level1)
nonScand.ibd.sample <- outlier.ibd # europe output
nonScand.type = data.frame(ID = outlier.level1,Type=c(rep('PCA',length(outlier.pca)), rep('IBD',length(outlier.ibd.new))))

#write.table(outlier.level1,"/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Trondheim.outlier.finePCADF7-IBD.nonScandinavian.2026Feb.txt",sep="\t",col.names=F,row.names=F,quote=FALSE)


##### Scand level ####
trond <- read.table("/maps/projects/ilab/people/pls394/plague/ancestry_version1/PCA/output/Trondheim.pvalue.scandice.DF15.202507.csv",header=T,sep=",",row.names=1)
p.thres <- 0.001

trond <- trond[!rownames(trond) %in% outlier.level1,] 
nonoutlier.level1 <- samples[! samples %in% outlier.level1] #n=118
trond <- trond[rownames(trond) %in% nonoutlier.level1,]  # dim = 118 15

### step1: PCA-outlier under P-value of 0.001
outlier.pca <- trond %>% filter(NO_Trondelag < p.thres ) %>% rownames() # n=18
nonoutlier.pca <- rownames(trond)[!rownames(trond) %in% outlier.pca] #n=100


### step2: IBD outliers
local <- c('NO_Trondelag')
nonlocal <- c("Denmark", "SE_Skane", "SE_South", "SE_Gotaland","NO_South" ,"SE_Stockholm", "NO_East" ,"SE_Central","NO_Central" ,"NO_Oslo" ,"NO_West" ,"SE_North" ,"NO_North") #it does not matter if Germany, Netherlands, Shetland and Orkeney are included or not

outlier.ibd <- c()
all_ibdtest <- data.frame()
for (sample in nonoutlier.pca) {
    df.sub = full %>% filter(id == sample) %>% filter(! county %in% c("British","Irish"))
    df.nonlocal = df.sub %>% filter(county %in% nonlocal)
  df.local= df.sub %>% filter(county %in% local)
    
  ind_ibdtest <- data.frame()
  for (i in seq_len(nrow(df.nonlocal))) {
      for (j in seq_len(nrow(df.local))) {
          t_output <- t.test2(
              m1 = df.nonlocal$cm_mean[i],
              m2 = df.local$cm_mean[j],
              s1 = df.nonlocal$cm_sd[i],
              s2 = df.local$cm_sd[j],
              n1 = df.nonlocal$county_n[i],
              n2 = df.local$county_n[j],
              m = 0,
              equal.variance = FALSE,
              alternative = "greater"
          )
          ind_ibdtest <- rbind(ind_ibdtest,
                               data.frame(
                                   id = sample,
                                   nonlocal_county = df.nonlocal$county[i],
                                   local_county = df.local$county[j],
                                   nonlocal_mean = df.nonlocal$cm_mean[i],
                                   local_mean = df.local$cm_mean[j],
                                   nonlocal_sd = df.nonlocal$cm_sd[i],
                                   local_sd = df.local$cm_sd[j],
                                   nonlocal_n = df.nonlocal$county_n[i],
                                   local_n = df.local$county_n[j],
                                   t = t_output[['t']],
                                   p_value = t_output[['p-value']]

                               )
                               )
      }
  }
  all_ibdtest = rbind(all_ibdtest,ind_ibdtest)

  ## summarize ind_ibdtest
  ind.summary = ind_ibdtest %>% group_by(nonlocal_county) %>%  summarise(significant = all(p_value < 0.05),  .groups = "drop") %>% as.data.frame()
  if (any(ind.summary$significant)) {
      outlier.ibd <- append(outlier.ibd, sample)
    }
  
} #loop through sample

outlier.ibd #n=21



outlier.ibd.new.level2 <- c()
for (sample in nonoutlier.pca) {
  df.sub = ibdsummary %>% filter(ID == sample) %>% filter(! ID2Grp  %in% c("British","Irish")) # Icelandic, British, Irish are only used in figure for comparison
                                                                                              # in this version, Icelandic is not used for identifying non-Scand outliers, but used for identifying non-local Scand.
  df.nonlocal = df.sub %>% filter(ID2Grp %in% nonlocal)
  df.nonlocal.top = df.nonlocal[which.max(df.nonlocal$shareID2Mean),]
  
  df.local= df.sub %>% filter(ID2Grp %in% local)
  df.local.top = df.local[which.max(df.local$shareID2Mean),]
  
  t_output <- t.test2(
      m1 = df.nonlocal.top$shareID2Mean,
      m2 = df.local.top$shareID2Mean,
      s1 = df.nonlocal.top$shareID2SD,
      s2 = df.local.top$shareID2SD,
      n1 = df.nonlocal.top$ID2Cnt,
      n2 = df.Scan.top$ID2Cnt,
      m = 0,
      equal.variance = FALSE,
      alternative = "greater"
      )

  if (t_output['p-value'] < 0.05) {
      outlier.ibd.new.level2 <- append(outlier.ibd.new.level2,sample)
  }
} 
outlier.ibd.new.level2
outlier.level2 <- c(outlier.pca,outlier.ibd.new.level2)%>% unique() # n = 39
print(outlier.level2)
nonLocal.ibd.sample <- outlier.level2 # Scand output
Local.ibd.sample<- rownames(trond)[! rownames(trond) %in% outlier.level2] # n = 79
nonLocal.type = data.frame(ID=outlier.level2,Type=c(rep('PCA',length(outlier.pca)), rep('IBD',length(outlier.ibd.new.level2))))
#write.table(outlier.level2,"/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Trondheim.outlier.finePCADF15-IBD.nonlocal-Scandinavian.2026Feb.txt",sep="\t",col.names=F,row.names=F,quote=FALSE)

Local.type = data.frame(ID=Local.ibd.sample,Type=c(rep(NA,length(Local.ibd.sample))))


### prepare matrix of standardized IBD and classification table
ancestry.trondheim.info =rbind(nonScand.type,nonLocal.type,Local.type)
names(ancestry.trondheim.info) = c('ID','Method')
ancestry.trondheim.info$Group = c(rep('non-Scandinavian',nrow(nonScand.type)), rep("nonlocal Scandinavian",nrow(nonLocal.type)),rep("local Scandinavian", nrow(Local.type)))
ancestry.trondheim.info = ancestry.trondheim.info %>% left_join(samples.df, by=c('ID'='id'))


colOrder <- c("Denmark","SE_Skane","SE_South","SE_Gotaland","SE_Stockholm","SE_Central","SE_North",
              "NO_Oslo","NO_South","NO_Central","NO_East","NO_West","NO_North",
              "NO_Trondelag",
              "Netherlands","Germany","Shetland","Orkney","Iceland"#,"British", "Irish"
              )
ibd.matrix = ibdsummary %>% filter( ID %in% ancestry.trondheim.info$ID) %>%
    select(ID,shareID2Mean,ID2Grp) %>%
    pivot_wider(names_from = ID2Grp, values_from = shareID2Mean)%>%
    as.data.frame()

ibd.matrix.std = standardize_rows(ibd.matrix[,-1])
rownames(ibd.matrix.std) <- ibd.matrix$ID
names(ibd.matrix.std)[names(ibd.matrix.std) == 'ICE'] <- 'Iceland'
ibd.matrix.std.m = as.matrix(ibd.matrix.std)[,colOrder,drop=FALSE]
ibd.matrix.std.transpos.m <- t(ibd.matrix.std.m)

write.table(ibd.matrix.std.transpos.m,file='/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Trondheim.IBD_standard.2026April28.tsv',sep='\t',col.names = T,row.names=T,quote=FALSE)
write.table(ancestry.trondheim.info, file='/maps/projects/ilab/people/pls394/plague/ancestry_final_2026Feb/output/Trondheim.IBD_info.2026April28.tsv',sep='\t',col.names = T,row.names=F,quote=FALSE)
