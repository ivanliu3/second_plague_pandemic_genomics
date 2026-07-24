module load plink/2.0.0

gemma=/home/pls394/Software/GEMMA-0985/gemma-0.98.5-linux-static-AMD64
in_dir=/maps/projects/ilab/people/pls394/plague/Association/GEMMA/input_2026April
out_dir=/maps/projects/ilab/people/pls394/plague/Association/GEMMA/output_2026April

VARSFORGRM=${in_dir}/GRM.INFO098_SNP_noXchr_noQ_rmAccess_rmLD_rmHWE.markerlist.txt

############# Lund ##############
PHENOTYPE=${in_dir}/Lund.forGEMMA.pheno.txt

### step0: extract only Pre- and Post-
cat /maps/projects/ilab/people/pls394/plague/imputation_WWW/Lund_WWW_221011_impPhased.fam  | grep -f <(cat ${in_dir}/Lund.Pre-Post.info.tsv Lund.Pre-Post.info.tsv | cut -f1) > ${in_dir}/Lund.Pre-Post.fam 
plink2 --bfile /maps/projects/ilab/people/pls394/plague/imputation_WWW/Lund_WWW_221011_impPhased --keep-fam ${in_dir}/Lund.Pre-Post.fam	--make-bed --out ${in_dir}/Lund.forGEMMA

### step1: calcualte GRM
GRMOUT=Lund.Pre-Post.GRM

${gemma} -bfile ${in_dir}/Lund.forGEMMA \
	 -gk 2 \ # standarized
	 -snps ${VARSFORGRM} \
	 -p ${PHENOTYPE} \
	 -outdir ${in_dir} \
	 -o ${GRMOUT} \
	 -maf 0.05 \
	 -miss 1 \ # Make sure GEMMA doesn't do any additional filtering  
	 -r2 1 \ # Make sure GEMMA doesn't do any additional filtering 

### step2: run LMM
${gemma} -bfile ${in_dir}/Lund.forGEMMA \
	 -p ${PHENOTYPE} \
	 -k ${in_dir}/${GRMOUT}.sXX.txt \
	 -c ${in_dir}/Lund.forGEMMA.intcpt_sex_PC1-7.txt \
	 -lmm 3 \
	 -maf 0.0001 \
	 -miss 1 \
	 -r2 1\
	 -outdir ${out_dir} \
	 -o Lund.Pre-Post

############# end of Lund ##############

############# Trondheim ##############
PHENOTYPE=${in_dir}/Trondheim.forGEMMA.pheno.txt

### step0: extract only Pre- and Post-
## use R to prepare Trondheim.Pre-Post.fam, excluding WF531, having 124 inds left
plink2 --bfile /maps/projects/ilab/people/pls394/plague/imputation_WWW/Trondheim_WWW_221011_impPhased --keep-fam ${in_dir}/Trondheim.Pre-Post.fam   --make-bed --out ${in_dir}/Trondheim.forGEMMA

### step1: calcualte GRM
GRMOUT=Trondheim.Pre-Post.GRM

${gemma} -bfile ${in_dir}/Trondheim.forGEMMA \
	 -gk 2 \ # standarized
	 -snps ${VARSFORGRM} \
	 -p ${PHENOTYPE} \
	 -outdir ${in_dir} \
	 -o ${GRMOUT} \
	 -maf 0.05 \
	 -miss 1 \ # Make sure GEMMA doesn't do any additional filtering  
	 -r2 1 \ # Make sure GEMMA doesn't do any additional filtering 

### step2: run LMM
${gemma} -bfile ${in_dir}/Trondheim.forGEMMA \
	 -p ${PHENOTYPE} \
	 -k ${in_dir}/${GRMOUT}.sXX.txt \
	 -c ${in_dir}/Trondheim.forGEMMA.intcpt_sex_PC1-7.txt \
	 -lmm 3 \
	 -maf 0.0001 \
	 -miss 1 \
	 -r2 1 \
	 -outdir ${out_dir} \
	 -o Trondheim.Pre-Post

############# end of Trondheim ##############

##### Vilnius #####
PHENOTYPE=${in_dir}/Vilnius.forGEMMA.pheno.txt

### step0: extract only Pre- and Post-
ln -s /maps/projects/ilab/people/pls394/plague/imputation_WWW/Lit_WWW_221011_impPhased.fam Vilnius.forGEMMA.fam
ln -s /maps/projects/ilab/people/pls394/plague/imputation_WWW/Lit_WWW_221011_impPhased.bim Vilnius.forGEMMA.bim
ln -s /maps/projects/ilab/people/pls394/plague/imputation_WWW/Lit_WWW_221011_impPhased.bed Vilnius.forGEMMA.bed

### step1: calcualte GRM
GRMOUT=Vilnius.Pre-Post.GRM

${gemma} -bfile ${in_dir}/Vilnius.forGEMMA \
	 -gk 2 \ # standarized
	 -snps ${VARSFORGRM} \
	 -p ${PHENOTYPE} \
	 -outdir ${in_dir} \
	 -o ${GRMOUT} \
	 -maf 0.05 \
	 -miss 1 \ # Make sure GEMMA doesn't do any additional filtering  
	 -r2 1 \ # Make sure GEMMA doesn't do any additional filtering 

### step2: run LMM
${gemma} -bfile ${in_dir}/Vilnius.forGEMMA \
	 -p ${PHENOTYPE} \
	 -k ${in_dir}/${GRMOUT}.sXX.txt \
	 -c ${in_dir}/Vilnius.forGEMMA.intcpt_sex_PC1-7.txt \
	 -lmm 3 \
	 -maf 0.0001 \
	 -miss 1 \
	 -r2 1\
	 -outdir ${out_dir} \
	 -o Vilnius.Pre-Post
