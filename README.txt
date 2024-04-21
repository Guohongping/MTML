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
 
  Notice: (1)there is a documentation called “run_order.txt” in current directory to illustrate the file running order.
          (2)there is a documentation called “para_simu.txt” in current directory to illustrate the simulation hyperparameters.

 MTML contains the following code and data files which are required to reproduce the results of the manuscript.
 
 ./data/simu_data/:
 
	    snps.csv
	    The dataset for genetic markers in Monte Carlo simulation experiments. 
	    Totally include 4,000 markers, 2000 markers from each of the first two chromosomes (chr.1 and 2). 
	    The positions range from 11,226,256 bp to 12,038,776 bp on chr.1, from 5,045,828 bp to 6,412,875 bp on chr.2.
	    
	    QTNs.csv
	    Four supposed pleiotropic QTNs (11,298,364 bp on chr.1; 5,066,968 bp on chr.2; 5,134,228 bp on chr.2; and 6,137,189 bp on chr.2). 
	    
	    simulate.R
	    R code for generating phenotype datasets using genotype data in three Monte Carlo simulation experiments. 
	    
	    simu_phenotype/:
	    The files for phenotypes under three different simulation settings.
	    For example "Y1_2", "1" represents the first Monte Carlo simulation experiment, and "2" represents the scenario of two traits.
	   
 ./code/:

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
		Enter R space, input "source("./code/MTML-master.R")" 
		
 ./data/real_data/:	
        
                gen_deal.csv (This file is larger than 25MB, it can not be uploaded directly in Github, so we compress it in gen_deal.zip, please unzip it.)
		Data preprocessing for the genotype data of Arabidopsis thaliana. The markers with missing samples and minor allele frequency (MAF) less than 0.01 are deleted, then 215,947 markers are remained for analysis.
		
		trait_f3.csv (This file is larger than 25MB, it can not be uploaded directly in Github, so we compress it in trait_f3.zip, please unzip it.)
		Data preprocessing for the phenotype datasets (three flowering-time related traits) of Arabidopsis thaliana. The three flowering-time related traits are:
		days to flowering time under long days with vernalization (LDV), days to flowering time under short days with vernalization (SDV), and days to flowering time under long days with two weeks vernalization (2W).
		
		real_ph_f3.3.Q
		Calculate the respective population structure matrix of three flowering-time related traits using ADMIXTURE software 
                (http://software.genetics.ucla.edu/admixture/download.html).
		
		real_format_conversion.R
		Convert Arabidopsis thaliana genotype data into binary format functions required by PLINK software.
		
		GEMMA_mvLMM_real.sh
		Perform GEMMA and mvLMM in Arabidopsis thaliana real datasets.
		
		gen_info.csv
		The imformation of SNP loci and their mapped genes in Arabidopsis thaliana. (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/download.html)
		
		real_gene.R
		The main R script for Arabidopsis thaliana real data analysis using the proposed MTML, FASTmrEMMA, GEMMA and mvLMM method.
		Enter R space, input "source("./data/real_data/real_gene.R")"
		
		
 ./compare_method_code/:
        
		FASTmrEMMA/:
		
		    myFASTmrEMMA.R
			R code of FASTmrEMMA (a single-trait multi-locus GWAS method, Wen et al., 2018). 
			
			FASTmrEMMA_step.R
			The main R script of FASTmrEMMA method.
			Enter R space, input "source("./compare_method_code/FASTmrEMMA/FASTmrEMMA_step.R")"
			
			FASTmrEMMA_real.R
			The main R script for Arabidopsis thaliana real data analysis using FASTmrEMMA method.
		        Enter R space, input "source("./compare_method_code/FASTmrEMMA/FASTmrEMMA_real.R")"
			
		GEMMA_mvLMM/:
		    
                        plink.map.R   plink.ped.R
			They are two functions to convert the genotype data into two specified formats, *.ped and *.map, respectively. 
			
			plink.ph.R
			It is a funtion to convert the phenotype data into specified format (*.txt). 
			
			format_conversion.R
			It is the code to obtain the output files (*.ped and *.map), which are requried as the inputs in both GEMMA and mvLMM.
			
			Enter R space, input "source("./compare_method_code/GEMMA_mvLMM/format_conversion.R")"
			
			GEMMA_mvLMM_step.sh
			There are four steps in this file, please follow the instructions step-by-step.
			
			%%The first step: use PLINK software to convert the files (*.ped and *.map) to binary formats (*.bed, *.bim, and *.fam), which are required in GEMMA and mvLMM.
			
			Notice: PLINK (v1.90b6.10 64-bit (17 Jun 2019)) software can be downloaded from http://www.cog-genomics.org/plink/1.9/, 
			        the installation instruction can be seen in http://zzz.bwh.harvard.edu/plink.
					After installing PLINK 1.9, please copy the executable file to the current path.
					
			%%The second step: use GEMMA software to calculate the correlation coefficient matrix between phenotypes, and save the output file for further analysis.
			
			Notice: GEMMA (0.98.5 (2021-08-25)) software can be downloaded from https://github.com/genetics-statistics/GEMMA,
			        and the installation instruction can be also seen in this link.
					After installing GEMMA 0.98.5, please copy the executable file to the current path.
					
			%%The third step: perform single-trait single-locus GWAS analysis using GEMMA software.
			
			%%The forth step: perform mvLMM analysis (i.e., a multi-trait single-locus GWAS method based on multivariate linear mixed model), which can be also implemented by GEMMA software.
			
			The above four steps can be conducted via "source GEMMA_mvLMM_step.sh". If there is something wrong, please run the four steps separately.
						
			GEMMA_mvLMM_statistic.R
			Calculate the statistical powers of four QTNs, average power and type I error rate using GEMMA and mvLMM in nine different scenarios.
			
			real_format_conversion.R
			Convert Arabidopsis thaliana genotype data into binary format functions required by PLINK software.
			
			GEMMA_mvLMM_real.sh
			Perform GEMMA and mvLMM in Arabidopsis thaliana real datasets.
			
			GEMMA_mvLMM_real.R
			Show the results of GEMMA and mvLMM in Arabidopsis thaliana real datasets.

 ./results_figures_tables/:
        
		figures/:
		The figures of this manuscript are produced by the academic mapping software Origin, which can be downloaded and installed at:​ 
                http://www.originsoft.cn/Single/Index/origin-download.
		    
			MTML_QTN_power.opju MTML_average_power.opju
			This file is an origin mapping file, and you can modify it at will after downloading Origin.
		    
			Figure1.png  Figure1.tif  Figure1.eps
		        Statistical powers of MTML, GEMMA, FASTmrEMMA and mvLMM methods for the four QTNs in three simulation scenarios.
			
			Figure2.png  Figure2.tif  Figure2.eps
			Average powers of MTML, GEMMA, FASTmrEMMA and mvLMM methods in three simulation scenarios (d = 5).
			
			Figure3.png  Figure3.tif  Figure3.eps
			Average powers of MTML in various numbers of traits (d = 2, 5, 10) in three simulation scenarios.
			
			Figure4.png  Figure4.tif  Figure4.eps
			Type 1 error rates (0.01%) of MTML, GEMMA, FASTmrEMMA and mvLMM methods in three simulation scenarios (d = 5).
			
			Figure5.png  Figure5.tif  Figure5.eps
			Type 1 error rates (0.01%) of MTML in various numbers of traits (d = 2, 5, 10) in three simulation scenarios.
			
		tables/:
		    
			QTN_statistic_power.xlsx
			Correspond to  Figure1. Statistical powers of MTML, GEMMA, FASTmrEMMA and mvLMM methods for the four QTNs in three simulation scenarios.
         
		        power_type1_error.xlsx
			Correspond to Figure2, Figure3, Figure4, Figure5.
			
			Arabidopsis_genes.xlsx
			The comprehensive overview of the genes identified by MTML in three Arabidopsis thaliana flowering-time related traits (LDV, SDV and 2W).
			
			
			
​
     
	    
