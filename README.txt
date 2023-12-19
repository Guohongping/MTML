 Supplementary information / reproducible research files for the manuscript 
 Title: "MTML: An efficient multi-trait multi-locus GWAS method based on Cauchy combination test"
 
 Authors: Hongping Guo, Tong Li, Yao Shi and Xiao Wang.
 Code was written by Hongping Guo and Tong Li.
 In case of questions or comments please contact guohongping@hbnu.edu.cn
 
 The code was written/evaluated in R with the following software versions:
 R version 3.6.3 (2020-02-29)
 Platform: x86_64-conda-linux-gnu (64-bit)
 Running under: CentOS Linux 7 (Core)
 
 Matrix products: default
 BLAS/LAPACK: ./anaconda3/envs/R3.7/lib/libopenblasp-r0.3.15.so

 locale:
  [1] LC_CTYPE=en_US.utf8       LC_NUMERIC=C
  [3] LC_TIME=en_US.utf8        LC_COLLATE=en_US.utf8
  [5] LC_MONETARY=en_US.utf8    LC_MESSAGES=en_US.utf8
  [7] LC_PAPER=en_US.utf8       LC_NAME=C
  [9] LC_ADDRESS=C              LC_TELEPHONE=C
 [11] LC_MEASUREMENT=en_US.utf8 LC_IDENTIFICATION=C

 attached base packages:
 [1] stats     graphics  grDevices utils     datasets  methods   base

 loaded via a namespace (and not attached):
 [1] compiler_3.6.3

 R packages：
 [1] dplyr_1.1.4
 [2] MASS_7.3.54
 [3] readxl_1.3.1
 [4] mvtnorm_1.1.3
 [5] parallel_3.6.3
 [6] foreach_1.5.2
 [7] iterators_1.0.14
 [8] doParallel_1.0.17
 [9] bigmemory_4.6.1
 
 
 This folder in MTML contains the following code and data files that can be used to reproduce all process of the manuscript.
 The explanation is as follows:
 
 .MTML_code_data/data/simu_data/:
 
	    snps.csv
	    The dataset for genetic markers in Monte Carlo simulation experiments. 
	    Totally include 4,000 markers, 2000 markers from each of the first two chromosomes (chr.1 and 2). 
	    The positions range from 11,226,256 bp to 12,038,776 bp on chr.1, from 5,045,828 bp to 6,412,875 bp on chr.2.
	    
	    QTNs.csv
	    Four supposed pleiotropic QTNs (11,298,364 bp on chr.1; 5,066,968 bp on chr.2; 5,134,228 bp on chr.2; and 6,137,189 bp on chr.2). 
	    
	    simulate.R
	    R code for generating phenotype datasets using genotype data in three Monte Carlo simulation experiments. 
	    
	    simu phenotype/:
	    The dataset files for phenotypes with three different simulation settings.
	    For example "Y1_2", "1" represents the first Monte Carlo simulation experiment, and "2" represents the scenario of two traits.
	   
 .MTML_code_data/code/:

		single_trait.R
		R code for single-trait analysis by EMMAX method.

		ridge.R 
		R code for multi-locus analysis by deshrinking ridge regression.

		cauchyComb.R
		R code for aggregating the individual p-values of multiple traits by Cauchy combination test.
        
		statistic.R
		R code for calculating the statistical powers of four QTNs, average power and type I error rate in three simulation scenarios.
        
		MTML.R        
		R code for multi-trait multi-locus GWAS, which combines the three steps in the study .
        
		MTML-master.R     
		The main R script of the proposed MTML in three simulation experiments.
		The script will run successfully by sourcing the above code files (simulate.R, single_trait.R, ridge.R, cauchyComb.R, statistic.R, MTML.R). 
		In the output file, it contains the results in nine different scenarios (i.e., three Monte Carlo simulation experiments under three different traits (2,5,10)). 
		To view the intermediate results, simply add the file "write. csv" at the end of each step.
		Enter R space, input "source("./MTML_code_data/code/MTML-master.R")" 
		
 .MTML_code_data/data/real_data/:	
        
		snps_all.csv
		The genotype dataset of Arabidopsis thaliana, it totally includes 199 samples, 5 chromosomes and 216,130 SNPs.
		
		gen_deal.csv	
		Data preprocessing for the genotype data of Arabidopsis thaliana. The markers with missing samples and minor allele frequency (MAF) less than 0.01 are deleted, 
                then 215,947 markers are remained for analysis.
		
		traits_3.csv
		The three flowering-time related traits are days to flowering time under long days with vernalization (LDV), days to flowering time under short days with 
                vernalization (SDV), and days to flowering time under long days with two weeks vernalization (2W).
		
		trait_f3.csv
		Data preprocessing for the phenotype datasets (three flowering-time related traits) of Arabidopsis thaliana.
		
		real_ph_f3.3.Q
		Calculate the respective population structure matrix of three flowering-time related traits using ADMIXTURE software 
                (http://software.genetics.ucla.edu/admixture/download.html).
		
		MTML_real.R
		The main R script for Arabidopsis thaliana real data analysis using the proposed MTML method.
		Enter R space, input "source("./MTML_code_data/data/real_data/MTML_real.R")"
		
 .MTML_code_data/compare_method_code/:
        
		FASTmrEMMA/:
		
		    myFASTmrEMMA.R
			R code of FASTmrEMMA (a single-trait multi-locus GWAS method, Wen et al., 2018). 
			
			FASTmrEMMA_step.R
			The main R script of FASTmrEMMA method.
			Enter R space, input "source("./MTML_code_data/compare_method_code/FASTmrEMMA/FASTmrEMMA_step.R")"
			
			FASTmrEMMA_real.R
			The main R script for Arabidopsis thaliana real data analysis using FASTmrEMMA method.
		        Enter R space, input "source("./MTML_code_data/compare_method_code/FASTmrEMMA/FASTmrEMMA_real.R")"
			
		GEMMA_mvLMM/:
		    
			plink.map.R  plink.ped.R  plink.ph.R format_conversion.R
			Convert genotype data into binary format functions which are required in PLINK software, then save the output files.
			The output files are requried as the inputs in both single-trait single-locus GWAS method GEMMA and multi-trait single-locus GWAS method mvLMM.
			Enter R space, input "source("./MTML_code_data/compare_method_code/GEMMA_mvLMM/format_conversion.R")"
			
			GEMMA_mvLMM_step.sh
			%%The first part of the code is to convert the data into the format required by GEMMA and mvLMM methods using PLINK software. 
			  PLINK software can be downloaded and installed through http://zzz.bwh.harvard.edu/plink/. Please refer to the website for specific instructions.
			%%The second part of the code is to use GEMMA software to calculate the correlation coefficient matrix between phenotypes, and save the output file for further analysis.
			Download and install GEMMA software from https://github.com/genetics-statistics/GEMMA."
			%%The third part of the code is GEMMA (a single-trait single-locus GWAS method based on linear mixed model) analysis using GEMMA software.
			%%The forth part of the code is mvLMM (a multi-trait single-locus GWAS method based on multivariate linear mixed model) analysis using GEMMA software.
			This code is a string of shell scripts for linux systems that can be run with " source GEMMA_mvLMM_step.sh". If this generates an error, copy it segment 
                        by segment and run it.
			
			GEMMA_mvLMM_statistic.R
			Calculate the statistical powers of four QTNs, average power and type I error rate using GEMMA and mvLMM in nine different scenarios.
			
			real_format_conversion.R
			Convert Arabidopsis thaliana genotype data into binary format functions required by PLINK software.
			
			GEMMA_mvLMM_real.sh
			Perform GEMMA and mvLMM in Arabidopsis thaliana real datasets.
			
			GEMMA_mvLMM_real.R
			Show the results of GEMMA and mvLMM in Arabidopsis thaliana real datasets.

 .MTML_code_data/results_figures_tables/:
        
		figures/:
		The figures of this manuscript are produced by the academic mapping software Origin, which can be downloaded and installed at:​ 
                http://www.originsoft.cn/Single/Index/origin-download.
		    
			MTML_QTN_power.opju MTML_average_power.opju
			This file is an origin mapping file, and you can modify it at will after downloading Origin.
		    
			Figure 1 .png  Figure 1 .tif  Figure 1 .eps
		        Statistical powers of MTML, GEMMA, FASTmrEMMA and mvLMM methods for the four QTNs in three simulation scenarios.
			
			Figure 2 .png  Figure 2 .tif  Figure 2 .eps
			Average powers of MTML, GEMMA, FASTmrEMMA and mvLMM methods in three simulation scenarios (d = 5).
			
			Figure 3 .png  Figure 3 .tif  Figure 3 .eps
			Average powers of MTML in various numbers of traits (d = 2, 5, 10) in three simulation scenarios.
			
			Figure 4 .png  Figure 4 .tif  Figure 4 .eps
			Type 1 error rates (0.01%) of MTML, GEMMA, FASTmrEMMA and mvLMM methods in three simulation scenarios (d = 5).
			
			Figure 5 .png  Figure 5 .tif  Figure 5 .eps
			Type 1 error rates (0.01%) of MTML in various numbers of traits (d = 2, 5, 10) in three simulation scenarios.
			
		tables/:
		    
			QTN statistic power.xlsx
			Correspond to  Figure 1. Statistical powers of MTML, GEMMA, FASTmrEMMA and mvLMM methods for the four QTNs in three simulation scenarios.
         
		        power_type 1 error .xlsx
			Correspond to Figure 2, Figure 3, Figure 4, Figure 5.
			
			Arabidopsis genes.xlsx
			The comprehensive overview of the genes identified by MTML in three Arabidopsis thaliana flowering-time related traits (LDV, SDV and 2W).
			
			
​
     
	    
