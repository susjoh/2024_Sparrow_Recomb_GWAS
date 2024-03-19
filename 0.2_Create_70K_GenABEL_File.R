#
# Convert merged data to GenABEL
# SEJ

#~~ Create GenABEL map

plinkmap <- read.table("70k_data/70K_200K_maf_geno_mind.bim", stringsAsFactors = F)[,c(1, 2, 4)]
write.table(plinkmap, "70k_data/70K_200K_maf_geno_mind.genabelmap", row.names = F, col.names = F, quote = F)

#~~ Create GenABEL phenos

famfile <- read.table("70k_data/70K_200K_maf_geno_mind.fam",  stringsAsFactors = F)
genpheno <- famfile[,c(2, 5)]
head(genpheno)
genpheno$V6 <- genpheno$V5
genpheno$V5[which(genpheno$V5 == 2)] <- 0
names(genpheno) <- c("id", "sex", "oldsex")

write.table(genpheno, "70k_data/70K_200K_maf_geno_mind.pheno", row.names = F, quote = F)

#~~ Create GenABEL object

system("plink --bfile 70k_data/70K_200K_maf_geno_mind --autosome-num 33 --recode --out 70k_data/70K_200K_maf_geno_mind")

GenABEL::convert.snp.ped(pedfile = "70k_data/70K_200K_maf_geno_mind.ped",
                         mapfile = "70k_data/70K_200K_maf_geno_mind.genabelmap",
                         outfile = "70k_data/70K_200K_maf_geno_mind.genabel",
                         traits = 1, strand = "u", mapHasHeaderLine = F)
library(GenABEL)

sparabel70 <- load.gwaa.data(phenofile = "70k_data/70K_200K_maf_geno_mind.pheno",
                               genofile  = "70k_data/70K_200K_maf_geno_mind.genabel")
save(sparabel70, file = "70k_data/70K_200K_maf_geno_mind.RData")

load("70k_data/70K_200K_maf_geno_mind.RData")

table(chromosome(sparabel70))

sparabel70 <- sparabel70[,chromosome(sparabel70) %in% 1:29]

x <- perid.summary(sparabel70)
xsnp <- GenABEL::summary.snp.data(sparabel70@gtdata)

x$id <- row.names(x)
xsnp$SNP <- row.names(xsnp)

write.table(x, "results/0_ID_Summary_Autosomal.txt", row.names = F, sep = "\t", quote = F)
write.table(xsnp, "results/0_SNP_Summary_Autosomal.txt", row.names = F, sep = "\t", quote = F)

