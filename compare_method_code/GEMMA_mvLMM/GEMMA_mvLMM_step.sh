#!/bin/bash
#****************************************** Plink bed et al ************************************************
./plink --file ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_2 --make-bed --out ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_2
./plink --file ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_2 --make-bed --out ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_2
./plink --file ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_2 --make-bed --out ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_2
./plink --file ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_5 --make-bed --out ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_5
./plink --file ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_5 --make-bed --out ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_5 
./plink --file ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_5 --make-bed --out ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_5 
./plink --file ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_10 --make-bed --out ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_10 
./plink --file ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_10 --make-bed --out ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_10
./plink --file ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_10 --make-bed --out ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_10  
#****************************************** GEMMA calculate the correlation coefficient matrix ************************************************
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_2 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_2.txt -gk 1 -o related1_2 
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_2 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_2.txt -gk 1 -o related2_2 
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_2 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_2.txt -gk 1 -o related3_2 
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_5 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_5.txt -gk 1 -o related1_5 
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_5 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_5.txt -gk 1 -o related2_5 
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_5 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_5.txt -gk 1 -o related3_5 
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_10 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_10.txt -gk 1 -o related1_10
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_10 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_10.txt -gk 1 -o related2_10
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_10 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_10.txt -gk 1 -o related3_10
mv ./output/related1_2.cXX.txt ./output/related2_2.cXX.txt ./output/related3_2.cXX.txt ./output/related1_5.cXX.txt ./output/related2_5.cXX.txt ./output/related3_5.cXX.txt ./output/related1_10.cXX.txt ./output/related2_10.cXX.txt ./output/related3_10.cXX.txt ./MTML_code_data/compare_method_code/GEMMA_mvLMM/
# #******************************************* GEMMA ************************************************
for ((i=1;i<=10;i++))
do
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_2 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_2.txt -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related1_2.cXX.txt -lmm 1 -n $i -o gemma1_2$i
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_2 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_2.txt -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related2_2.cXX.txt -lmm 1 -n $i -o gemma2_2$i
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_2 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_2.txt -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related3_2.cXX.txt -lmm 1 -n $i -o gemma3_2$i
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_5 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_5.txt -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related1_5.cXX.txt -lmm 1 -n $i -o gemma1_5$i
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_5 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_5.txt -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related2_5.cXX.txt -lmm 1 -n $i -o gemma2_5$i
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_5 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_5.txt -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related3_5.cXX.txt -lmm 1 -n $i -o gemma3_5$i
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_10 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_10.txt -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related1_10.cXX.txt -lmm 1 -n $i -o gemma1_10$i
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_10 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_10.txt -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related2_10.cXX.txt -lmm 1 -n $i -o gemma2_10$i
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_10 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_10.txt -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related3_10.cXX.txt -lmm 1 -n $i -o gemma3_10$i
done
# #******************************************* mvLMM ************************************************
# d=2
for ((i=1;i<=10;i++))
do
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_2 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_2.txt  -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related1_2.cXX.txt  -lmm 1 -n $[2$i-1] $[2$i] -o mvlmm1_2$i
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_2 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_2.txt  -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related2_2.cXX.txt  -lmm 1 -n $[2$i-1] $[2$i] -o mvlmm2_2$i
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_2 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_2.txt  -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related3_2.cXX.txt  -lmm 1 -n $[2$i-1] $[2$i] -o mvlmm3_2$i  
done
# d=5
for ((i=1;i<=5;i++))
do
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_5 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_5.txt  -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related1_5.cXX.txt  -lmm 1 -n $[5$i-4] $[5$i-3] $[5$i-2] $[5$i-1] $[5$i] -o mvlmm1_5$i  
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_5 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_5.txt  -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related2_5.cXX.txt  -lmm 1 -n $[5$i-4] $[5$i-3] $[5$i-2] $[5$i-1] $[5$i] -o mvlmm2_5$i  
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_5 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_5.txt  -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related3_5.cXX.txt  -lmm 1 -n $[5$i-4] $[5$i-3] $[5$i-2] $[5$i-1] $[5$i] -o mvlmm3_5$i  
done
# d=10
for ((i=1;i<=2;i++))
do
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_10 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y1_10.txt  -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related1_10.cXX.txt  -lmm 1 -n $[10$i-9] $[10$i-8] $[10$i-7] $[10$i-6] $[10$i-5] $[10$i-4] $[10$i-3] $[10$i-2] $[10$i-1] $[10$i] -o mvlmm1_10$i  
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_10 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y2_10.txt  -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related2_10.cXX.txt  -lmm 1 -n $[10$i-9] $[10$i-8] $[10$i-7] $[10$i-6] $[10$i-5] $[10$i-4] $[10$i-3] $[10$i-2] $[10$i-1] $[10$i] -o mvlmm2_10$i  
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_10 -p ./MTML_code_data/compare_method_code/GEMMA_mvLMM/y3_10.txt  -k ./MTML_code_data/compare_method_code/GEMMA_mvLMM/related3_10.cXX.txt  -lmm 1 -n $[10$i-9] $[10$i-8] $[10$i-7] $[10$i-6] $[10$i-5] $[10$i-4] $[10$i-3] $[10$i-2] $[10$i-1] $[10$i] -o mvlmm3_10$i  
done
