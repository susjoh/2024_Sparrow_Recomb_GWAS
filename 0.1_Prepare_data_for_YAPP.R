#
# Prepare the data for Yapp.
# SEJ
#

library(ggplot2)
library(dplyr)
library(reshape2)

#~~ Update the pedigree information in the fam file

pedigree <- read.table("70k_data/Helgeland_pedigree_2023-03-17.txt", header = T)
famfile <- read.table("70k_data/70K_200K_maf_geno_mind.fam")

names(famfile)[2] <- "id"
famfile <- left_join(famfile, pedigree)
head(famfile)

table(famfile$V4[which(famfile$V4 != 0)] == famfile$dam[which(famfile$V4 != 0)])
table(famfile$V3[which(famfile$V3 != 0)] == famfile$sire[which(famfile$V3 != 0)])

famfile$dam[which(!famfile$dam %in% famfile$id)] <- 0
famfile$sire[which(!famfile$sire %in% famfile$id)] <- 0

write.table(famfile[,c("V1", "id", "sire", "dam")], "70k_data/update_parents.txt", row.names = F, col.names = F, quote = F)

#~~ Check sex - update from the 180k project and make sires and dams 1 or 2

newsex <- read.table("70k_data/update_sexes_v2.txt", header = T)[,2:3]
head(newsex)
names(newsex) <- c("id", "sex")

famfile <- left_join(famfile, newsex)
famfile$sex[which(is.na(famfile$sex))] <- 0
famfile$sex[which(famfile$id %in% famfile$sire)] <- 1
famfile$sex[which(famfile$id %in% famfile$dam)] <- 2

write.table(famfile[,c("V1", "id", "sex")], "70k_data/update_sexes_v3.txt", row.names = F, quote = F, col.names = F)

#~~ Only include autosomes in the yapp study.

bimfile <- read.table("70k_data/70K_200K_maf_geno_mind.bim")
head(bimfile)
table(bimfile$V1)
writeLines(subset(bimfile, V1 %in% 1:29)$V2, "70k_data/autosomal_SNPs.txt")

#~~ Create PLINK files...

system("plink --bfile 70k_data/70K_200K_maf_geno_mind --autosome-num 32 --update-parents 70k_data/update_parents.txt --extract 70k_data/autosomal_SNPs.txt --make-bed --out 70k_data/70K_200K_maf_geno_mind_v2")
system("plink --bfile 70k_data/70K_200K_maf_geno_mind_v2 --autosome-num 32 --update-sex 70k_data/update_sexes_v3.txt --make-bed --out 70k_data/70K_200K_maf_geno_mind_v3")

#~~ Check Mendelian errors & missing values

system("plink --bfile 70k_data/70K_200K_maf_geno_mind_v3 --autosome-num 32 --mendel --out 70k_data/70K_200K_maf_geno_mind_v3")
system("plink --bfile 70k_data/70K_200K_maf_geno_mind_v3 --autosome-num 32 --missing --out 70k_data/70K_200K_maf_geno_mind_v3")

beepr::beep()

#~~ Problem loci?

snpmiss <- read.table("70k_data/70K_200K_maf_geno_mind_v3.lmiss", header = T)
snpmend <- read.table("70k_data/70K_200K_maf_geno_mind_v3.lmendel", header = T)

snpmiss <- left_join(snpmiss, snpmend)
ggplot(snpmiss, aes(N_MISS, N)) + geom_point()
ggplot(snpmiss, aes(N)) + geom_histogram(binwidth = 10) + geom_vline(xintercept = 100)

#~~ Remove loci with more than 100 Mendelian errors

removesnps <- subset(snpmiss, N > 100)$SNP
writeLines(removesnps, "70k_data/remove_snps.txt")

#~~ Problem ids?

idmiss <- read.table("70k_data/70K_200K_maf_geno_mind_v3.imiss", header = T)
ggplot(idmiss, aes(N_MISS)) + geom_histogram(binwidth = 10) + geom_vline(xintercept = 100)

#~~ Don't remove IDs with more than 2000 Missing genos for now...

# removeids <- subset(idmiss, N_MISS > 2000)[,c("FID", "IID")]
# table(removeids$IID %in% c(famfile$dam, famfile$sire))

#~~ Update the dataset again again:

system("plink --bfile 70k_data/70K_200K_maf_geno_mind_v3 --autosome-num 32 --exclude 70k_data/remove_snps.txt --make-bed --out 70k_data/70K_200K_maf_geno_mind_v4")

system("plink --bfile 70k_data/70K_200K_maf_geno_mind_v4 --autosome-num 32 --mendel --out 70k_data/70K_200K_maf_geno_mind_v4")
system("plink --bfile 70k_data/70K_200K_maf_geno_mind_v4 --autosome-num 32 --missing --out 70k_data/70K_200K_maf_geno_mind_v4")

#~~ Problem ids again?

idmiss <- read.table("70k_data/70K_200K_maf_geno_mind_v4.imiss", header = T)
ggplot(idmiss, aes(N_MISS)) + geom_histogram(binwidth = 10) + geom_vline(xintercept = 100)

idmend <- read.table("70k_data/70K_200K_maf_geno_mind_v4.imendel", header = T)
head(idmend)

idmiss <- left_join(idmiss, idmend)
idmiss$N[which(is.na(idmiss$N))] <- 0

ggplot(idmiss, aes(N, N_MISS)) + geom_point()

#~~ examine Mendelian errors

menderr <- read.table("70k_data/70K_200K_maf_geno_mind_v4.mendel", skip = 1)
head(menderr)
menderr$ErrorSource <- "Both"
menderr$ErrorSource[which(menderr$V5 %in% c(3, 6))] <- "Father"
menderr$ErrorSource[which(menderr$V5 %in% c(4, 7))] <- "Mother"
table(menderr$ErrorSource)

x <- menderr %>% group_by(V2, ErrorSource) %>% summarise(ErrorCount = n())
head(x)
names(x)[1] <- c("Child")

x <- dcast(x, Child ~ ErrorSource)
head(x)

ggplot(x, aes(Both)) + geom_histogram(binwidth = 100) + geom_vline(xintercept = 1500)
ggplot(x, aes(Mother)) + geom_histogram(binwidth = 100) + geom_vline(xintercept = 750)
ggplot(x, aes(Father)) + geom_histogram(binwidth = 100) + geom_vline(xintercept = 750)

ggplot(x, aes(Both, Mother)) + geom_point() + geom_vline(xintercept = 1500)
ggplot(x, aes(Both, Father)) + geom_point() + geom_vline(xintercept = 1500)

#~~ It all looks fine!


#~~ make vcf

system("plink --bfile 70k_data/70K_200K_maf_geno_mind_v4 --autosome-num 32 --recode vcf-iid --out 70k_data/70K_200K_maf_geno_mind_v4")

#~~ Log onto UBUNTU partition and run:

# cd /mnt/g/My\ Drive/House\ Sparrows/20230606_Sparrow_70k_analysis/70k_data/
# bgzip -c 70K_200K_maf_geno_mind_v4.vcf > 70K_200K_maf_geno_mind_v4.vcf.gz
# tabix -fp vcf 70K_200K_maf_geno_mind_v4.vcf.gz



