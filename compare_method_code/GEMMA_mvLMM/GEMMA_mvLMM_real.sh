


#****************************************** Plink bed et al ************************************************
./plink --file ./MTML_code_data/data/real_data/trait_f3 --make-bed --out ./MTML_code_data/data/real_data/trait_f3

#******************************************** Gemma related matrix************************************************
./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/data/real_data/trait_f3 -p ./MTML_code_data/data/real_data/trait_f3.txt -gk 1 -o relatedf3

mv ./output/relatedf3.cXX.txt ./MTML_code_data/data/real_data/
# #******************************************* Gemma ************************************************
for ((i=1;i<=3;i++))
do
  ./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/data/real_data/trait_f3 -p ./MTML_code_data/data/real_data/trait_f3.txt  -k ./MTML_code_data/data/real_data/relatedf3.cXX.txt -lmm 1 -n $i -o gemmaf3$i
done


# #******************************************* mvLMM ************************************************

./gemma-0.98.5-linux-static-AMD64 -bfile ./MTML_code_data/data/real_data/trait_f3 -p ./MTML_code_data/data/real_data/trait_f3.txt  -k ./MTML_code_data/data/real_data/relatedf3.cXX.txt  -lmm 1 -n 1 2 3 -o mvlmmf3 

