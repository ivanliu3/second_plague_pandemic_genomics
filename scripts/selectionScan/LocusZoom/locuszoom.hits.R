##library(EnsDb.Hsapiens.v75) #Note this is for build37
library(AnnotationHub)
ah <- AnnotationHub()

ensdb105 <- query(ah, c("EnsDb", "Hsapiens", "105"))[[1]]
library(locuszoomr)
library(dplyr)
library(data.table)
library(stringr)
#source("//maps/projects/ilab/people/pls394/plague/Association/Plink2/output/2025_preVSmodern/LocusZoom/functions/locus_zoom.R")

data(SLE_gwas_sub)  ## limited subset of data from SLE GWAS
head(SLE_gwas_sub)


### SNP1 chr1:225574884:SG.ld
## significnat in Lund
lund.preVmod = fread("/maps/projects/ilab/people/pls394/plague/Association_modern/GEMMA_output/2026Feb20/filter/Lund-danish_skane.pre_vs_present.mm.7_pcs.maf005.rmMULTI-INDEL.rmLOWINFO095.rmLOWGP.rmQ.rmHWE1e-6.modernDenmarkHWE1e-6.rmwithin5bp.tsv",header=T)
lund.preVmod$CHR <- as.integer(lund.preVmod$chr)
lund.preVmod$BP <- as.numeric(lund.preVmod$ps)
lund.preVmod = lund.preVmod %>% arrange(CHR,BP)

## modern Danish+Skane
ldfile <- "/maps/projects/ilab/people/pls394/plague/LocusZoom/package_made_2026May/CHR1_225574884_danish_skane.GZ"
ld <- fread(ldfile, header = TRUE) %>% as.data.frame()
outfile = "/maps/projects/ilab/people/pls394/plague/LocusZoom/package_made_2026Jun23/Lund.chr1.colorByLD.1M.png"
snp_id <- strsplit(ldfile, "\\.")[[1]][2]
snp_id = 'chr1:225574884:SG'
chrom<- sub("chr","",sub(".*(chr[0-9]+):.*", "\\1", ldfile) ) %>% as.integer
chrom = as.integer(1 )
position <- as.numeric(sub(".*chr[0-9]+:([0-9]+):.*", "\\1", ldfile))
position = 225574884
pos1=position-500000
pos2=position+500000


lund.subset = lund.preVmod %>% filter(CHR == chrom & BP >=pos1 & BP<=pos2)
merged_data <- lund.subset %>%
  left_join(ld, by = c("rs" = "SNP_B"))
names(merged_data)[1] = 'rsid'
names(merged_data)[ncol(merged_data)] = 'r2'

loc <- locus(data=merged_data, 
             seqname=chrom,xrange=c(pos1,pos2), 
             chrom="CHR",
             pos="BP",
             p="p_score",
             ens_db = ensdb105, #Ensembl or AnnotationHub database
             index_snp = "chr1:225574884:SG", # index snp, which calculate ld with
             LD = "r2",  # Color by R2
             )


bitmap(file=outfile,h=6,w=8,res=300)
locus_plot(loc, labels = c("index"),
           label_x = c(4, -5),
           pcutoff=NULL)
dev.off()



### SNP2 chr2:160798465:SG:1 

## modern Danish+Skane
ldfile <- "/maps/projects/ilab/people/pls394/plague/LocusZoom/package_made_2026May/CHR2_160798465_danish_skane.GZ"
ld <- fread(ldfile, header = TRUE) %>% as.data.frame()
outfile = "/maps/projects/ilab/people/pls394/plague/LocusZoom/package_made_2026Jun23/Lund.chr2.colorByLD.1M.png"

snp_id <- strsplit(ldfile, "\\.")[[1]][2]
snp_id = 'chr2:160798465:SG:1'
chrom <- sub("chr","",sub(".*(chr[0-9]+):.*", "\\1", ldfile) ) %>% as.integer
chrom = as.integer(2 )
position <- as.numeric(sub(".*chr[0-9]+:([0-9]+):.*", "\\1", ldfile))
position = 160798465
pos1=position-500000
pos2=position+500000


lund.subset = lund.preVmod %>% filter(CHR == chrom & BP >=pos1 & BP<=pos2)
merged_data <- lund.subset %>%
  left_join(ld, by = c("rs" = "SNP_B"))
names(merged_data)[1] = 'rsid'
names(merged_data)[ncol(merged_data)] = 'r2'

loc <- locus(data=merged_data, 
             seqname=chrom,
             xrange=c(pos1,pos2), 
             chrom="CHR",
             pos="BP",
             p="p_score",
             ens_db = ensdb105, #Ensembl or AnnotationHub database
             index_snp = "chr2:160798465:SG:1", # index snp, which calculate ld with
             LD = "r2",  # Color by R2
             )


bitmap(file=outfile,h=6,w=8,res=300)
locus_plot(loc, labels = c("index"),
           label_x = c(4, -5),
           pcutoff=NULL)
dev.off()

### SNP3 chr3:171104061:IG
## modern Danish+Skane
ldfile <- "/maps/projects/ilab/people/pls394/plague/LocusZoom/package_made_2026May/CHR3_171104061_danish_skane.GZ"
ld <- fread(ldfile, header = TRUE) %>% as.data.frame()
outfile = "/maps/projects/ilab/people/pls394/plague/LocusZoom/package_made_2026Jun23/Lund.chr3.colorByLD.1M.png"

snp_id <- strsplit(ldfile, "\\.")[[1]][2]
snp_id = 'chr3:171104061:IG'
chrom <- sub("chr","",sub(".*(chr[0-9]+):.*", "\\1", ldfile) ) %>% as.integer
chrom = as.integer(3)
position <- as.numeric(sub(".*chr[0-9]+:([0-9]+):.*", "\\1", ldfile))
position = 171104061
pos1=position-500000
pos2=position+500000


lund.subset = lund.preVmod %>% filter(CHR == chrom & BP >=pos1 & BP<=pos2)
merged_data <- lund.subset %>%
  left_join(ld, by = c("rs" = "SNP_B"))
names(merged_data)[1] = 'rsid'
names(merged_data)[ncol(merged_data)] = 'r2'

loc <- locus(data=merged_data, 
             seqname=chrom,
             xrange=c(pos1,pos2), 
             chrom="CHR",
             pos="BP",
             p="p_score",
             ens_db = ensdb105, #Ensembl or AnnotationHub database
             index_snp = "chr3:171104061:IG",
             LD = "r2",  # Color by R2
             )


bitmap(file=outfile,h=6,w=8,res=300)
locus_plot(loc, labels = c("index"),
           label_x = c(4, -5),
           pcutoff=NULL)
dev.off()


### SNP4 chr4:139027159:SG  
## modern Danish+Skane
ldfile <- "/maps/projects/ilab/people/pls394/plague/LocusZoom/package_made_2026May/CHR4_139027159_danish_skane.GZ"
ld <- fread(ldfile, header = TRUE) %>% as.data.frame()
outfile = "/maps/projects/ilab/people/pls394/plague/LocusZoom/package_made_2026Jun23/Lund.chr4.colorByLD.1M.png"

snp_id <- strsplit(ldfile, "\\.")[[1]][2]
snp_id = 'chr4:139027159:SG'
chrom <- sub("chr","",sub(".*(chr[0-9]+):.*", "\\1", ldfile) ) %>% as.integer
chrom = as.integer(4)
position <- as.numeric(sub(".*chr[0-9]+:([0-9]+):.*", "\\1", ldfile))
position = 139027159
pos1=position-500000
pos2=position+500000


lund.subset = lund.preVmod %>% filter(CHR == chrom & BP >=pos1 & BP<=pos2)
merged_data <- lund.subset %>%
  left_join(ld, by = c("rs" = "SNP_B"))
names(merged_data)[1] = 'rsid'
names(merged_data)[ncol(merged_data)] = 'r2'

loc <- locus(data=merged_data, 
             seqname=chrom,
             xrange=c(pos1,pos2), 
             chrom="CHR",
             pos="BP",
             p="p_score",
             ens_db = ensdb105, #Ensembl or AnnotationHub database
             index_snp = "chr4:139027159:SG",
             LD = "r2",  # Color by R2
             )


bitmap(file=outfile,h=6,w=8,res=300)
locus_plot(loc, labels = c("index"),
           label_x = c(4, -5),
           pcutoff=NULL)
dev.off()

### SNP5 chr14:22193127:SG                                                                                                                                                                                                                                                                                                                                             ## only signiciant in Lund
## modern Danish+Skane
ldfile <- "/maps/projects/ilab/people/pls394/plague/Hits_allGWAS/LD_modernpop/chr14:22193217:SG.danish_skane.ld_2026jan.ld.gz"
ld <- fread(ldfile, header = TRUE) %>% as.data.frame()
outfile = "/maps/projects/ilab/people/pls394/plague/LocusZoom/package_made_2026Jun23/Lund.chr14.colorByLD.1M.png"

snp_id <- strsplit(ldfile, "\\.")[[1]][2]
snp_id = 'chr14:22193127:SG'
chrom <- sub("chr","",sub(".*(chr[0-9]+):.*", "\\1", ldfile) ) %>% as.integer
chrom = as.integer(14)
position <- as.numeric(sub(".*chr[0-9]+:([0-9]+):.*", "\\1", ldfile))
position = 22193127
pos1=position-500000
pos2=position+500000


lund.subset = lund.preVmod %>% filter(CHR == chrom & BP >=pos1 & BP<=pos2)
merged_data <- lund.subset %>%
  left_join(ld, by = c("rs" = "SNP_B"))
names(merged_data)[1] = 'rsid'
names(merged_data)[ncol(merged_data)] = 'r2'

loc <- locus(data=merged_data, 
             seqname=chrom,
             xrange=c(pos1,pos2), 
             chrom="CHR",
             pos="BP",
             p="p_score",
             ens_db = ensdb105, #Ensembl or AnnotationHub database
             index_snp = "chr14:22193127:SG",
             LD = "r2",  # Color by R2
             )


bitmap(file=outfile,h=6,w=8,res=300)
locus_plot(loc, labels = c("index"),
           label_x = c(4, -5),
           pcutoff=NULL)
dev.off()

