#!/bin/bash

plink --bfile 70K_200K_maf_geno_mind_v5 --autosome-num 32 --recode vcf-iid --out 70K_200K_maf_geno_mind_v5
bgzip -c 70K_200K_maf_geno_mind_v5.vcf > 70K_200K_maf_geno_mind_v5.vcf.gz
tabix -fp vcf 70K_200K_maf_geno_mind_v5.vcf.gz