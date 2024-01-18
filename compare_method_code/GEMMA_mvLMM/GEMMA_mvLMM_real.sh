


#****************************************** Plink bed et al ************************************************
./plink --file ./data/real_data/trait_f3 --make-bed --out ./data/real_data/trait_f3

#******************************************** Gemma ************************************************
./gemma-0.98.5-linux-static-AMD64 -bfile ./data/real_data/trait_f3 -p ./data/real_data/trait_f3.txt -gk 1 -o relatedf3
mv ./output/relatedf3.cXX.txt ./data/real_data/

for ((i=1;i<=3;i++))
do
  ./gemma-0.98.5-linux-static-AMD64 -bfile ./data/real_data/trait_f3 -p ./data/real_data/trait_f3.txt  -k ./data/real_data/relatedf3.cXX.txt -lmm 1 -n $i -o gemmaf3$i
done
