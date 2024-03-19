library(GenABEL)
library(dplyr)
library(ggplot2)

load("70k_data/70K_200K_maf_geno_mind.RData")

#~~ Define IDs based on whether they were on the 200k or 70 chips

fam200 <- read.table("../../20230606_Sparrow_180k_analysis/data/Pdo_200k_n3960_21032017.fam")

sparabel200 <- sparabel70[which(idnames(sparabel70) %in% fam200$V2),]
sparabel70 <- sparabel70[which(!idnames(sparabel70) %in% fam200$V2),]

nids(sparabel200)
nids(sparabel70)

snp70 <- summary.snp.data(sparabel70@gtdata)
snp200 <- summary.snp.data(sparabel200@gtdata)

snp200$Chromosome <- NULL
snp200$Position <- NULL
snp200$Strand <- NULL

names(snp200) <- paste0(names(snp200), "_200")

snp70$SNP <- row.names(snp70)
snp200$SNP <- row.names(snp200)

snp70 <- left_join(snp70, snp200)

table(snp70$A1 == snp70$A1_200)



#~~ Compare allele frequences

# MAF
snp70$MAFdiff <- snp70$Q.2 - snp70$Q.2_200
ggplot(snp70, aes(MAFdiff)) +
  geom_histogram(binwidth = 0.01) +
  coord_cartesian(xlim = c(-0.2, 0.2))

snp70$KeepFlag <- "keep"
snp70$KeepFlag[which(snp70$MAFdiff < (mean(snp70$MAFdiff) - 3*sd(snp70$MAFdiff)))] <- "remove"
snp70$KeepFlag[which(snp70$MAFdiff > (mean(snp70$MAFdiff) + 3*sd(snp70$MAFdiff)))] <- "remove"

# Deviation from HWE
snp70$Fdiff   <- snp70$Fmax - snp70$Fmax_200
ggplot(snp70, aes(Fdiff)) +
  geom_histogram(binwidth = 0.005) +
  coord_cartesian(xlim = c(-0.2, 0.2))

snp70$KeepFlag[which(snp70$Fdiff < (mean(snp70$Fdiff) - 3*sd(snp70$Fdiff)))] <- "remove"
snp70$KeepFlag[which(snp70$Fdiff > (mean(snp70$Fdiff) + 3*sd(snp70$Fdiff)))] <- "remove"

table(snp70$KeepFlag)

# Proportion of heterozygotes
snp70$PropHet = snp70$P.12/(snp70$P.11 + snp70$P.12 + snp70$P.22)
snp70$PropHet_200 = snp70$P.12_200/(snp70$P.11_200 + snp70$P.12_200 + snp70$P.22_200)
snp70$Hetdiff <- snp70$PropHet - snp70$PropHet_200
ggplot(snp70, aes(Hetdiff)) +
  geom_histogram(binwidth = 0.005) +
  coord_cartesian(xlim = c(-0.2, 0.2))

snp70$KeepFlag[which(snp70$Hetdiff < (mean(snp70$Hetdiff) - 3*sd(snp70$Hetdiff)))] <- "remove"
snp70$KeepFlag[which(snp70$Hetdiff > (mean(snp70$Hetdiff) + 3*sd(snp70$Hetdiff)))] <- "remove"


ggplot(snp70, aes(Fdiff, Hetdiff, col = KeepFlag)) + 
  geom_point(alpha = 0.1)

ggplot(snp70, aes(Fdiff, MAFdiff, col = KeepFlag)) + 
  geom_point(alpha = 0.1)

# Difference in P for HWE test
snp70$Pdiff <- log(snp70$Pexact)-log(snp70$Pexact_200)
snp70$Pdiff <- ifelse(snp70$Pdiff < 0, snp70$Pdiff * -1, snp70$Pdiff)
snp70$Pdiff[is.infinite(snp70$Pdiff)] <- 1000

ggplot(snp70, aes(Pdiff, MAFdiff, col = KeepFlag)) + 
  geom_point(alpha = 0.1)

ggplot(snp70, aes(Pdiff, Fdiff, col = KeepFlag)) + 
  geom_point(alpha = 0.1)

ggplot(snp70, aes(-log(Pexact), -log(Pexact_200), col = KeepFlag)) + 
  geom_point(alpha = 0.1)

snps_to_remove <- subset(snp70, KeepFlag == "remove")$SNP

writeLines(snps_to_remove, "70k_data/poor_qual_SNPs.txt")

#~~ Output the new ID information

load("70k_data/70K_200K_maf_geno_mind.RData")

nsnps(sparabel70)
sparabel70 <- sparabel70[,!snpnames(sparabel70) %in% snps_to_remove]

x <- perid.summary(sparabel70)
x$id <- row.names(x)

write.table(x, "results/0.3_ID_Summary_Autosomal.txt", row.names = F, sep = "\t", quote = F)

#~~ Make new PLINK files...


system("plink --bfile 70k_data/70K_200K_maf_geno_mind_v4 --autosome-num 32 --exclude 70k_data/poor_qual_SNPs.txt --make-bed --out 70k_data/70K_200K_maf_geno_mind_v5")


#~~ make vcf

system("plink --bfile 70k_data/70K_200K_maf_geno_mind_v5 --autosome-num 32 --recode vcf-iid --out 70k_data/70K_200K_maf_geno_mind_v5")


# cd /mnt/g/My\ Drive/House\ Sparrows/20230606_Sparrow_70k_analysis/70k_data/
# bgzip -c 70K_200K_maf_geno_mind_v5.vcf > 70K_200K_maf_geno_mind_v5.vcf.gz
# tabix -fp vcf 70K_200K_maf_geno_mind_v5.vcf.gz

