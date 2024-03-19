library(dplyr)
library(ggplot2)

load("results/2_Full_QC_Sparrow_Yapp_Data.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Check CO Count per region  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(rectab_v2)
rectab_v2$bin <- ceiling(rectab_v2$mid_CO_pos/5e5)

x <- rectab_v2 %>% group_by(chrom, bin) %>% summarise(NCO = n())
head(x)
x$Cumu <- 1:nrow(x)

ggplot(x, aes(Cumu, NCO, col = factor(chrom %% 2), group = chrom)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set1")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Check LD                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(reshape2)

for(i in 1:32){
  
  system(paste0("plink --bfile 70k_data/70K_200K_maf_geno_mind_v5 --autosome-num 32 --chr ", i, " --r2 square --make-bed --out test"))
  
  try({
    ldtab <- read.table("test.ld")
    ldtab <- cbind(V0 = 1:nrow(ldtab), ldtab)
    
    ldtab <- melt(ldtab, id.vars = "V0")
    ldtab$variable <- as.numeric(as.character(gsub("V", "", ldtab$variable)))
    
    ldtab <- subset(ldtab, V0 > variable)
    head(ldtab)
    
    ggplot(ldtab, aes(V0, variable, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      theme_bw() +
      labs(x = "SNP Order", y = "SNP Order") +
      ggtitle(paste0("Chromosome ", i, " LD"))
    
    ggsave(paste0("figs/2_Chr", i, "_LD.png"), width = 20, height = 20)
  })
}


# SNP DENSITY

snptab <- read.table("70k_data/70K_200K_maf_geno_mind_v5.bim")
head(snptab)

distvec <- diff(snptab$V4)
distvec <- distvec[-which(distvec < 0)]

mean(distvec)
median(distvec)

# Z SNP DENSITY

snptab <- read.table("70k_data/70K_200K_maf_geno_mind.bim")
snptab <- subset(snptab, V1 == 32)
distvec <- diff(snptab$V4)

mean(distvec)
median(distvec)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Check LD Decay               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# We investigated the linkage disequilibrium (LD) decay profile for each
# chromosome for all loci occurring within 500kb windows using the flag
# --ld-window-kb 500 in PLINK v1.90. LD decayed to r2 = 0.05 at a distance of
# ~100kb, with LD decay slightly faster for micro-chromosomes (Figure S1).

system("plink --bfile 70k_data\\70K_200K_maf_geno_mind_v5 --autosome-num 33 --maf 0.05 --r2 --ld-window-r2 0 --ld-window-kb 500 --out 70k_data\\70K_200K_maf_geno_mind_v5")

ldtab <- read.table("70K_data/70K_200K_maf_geno_mind_v5.ld", header = T)
ldtab$type <- ifelse(ldtab$CHR_A %in% c(1:12, 29), "Macro-chromosomes", "Micro-chromosomes")


head(ldtab)
ldtab$Dist <- ldtab$BP_B - ldtab$BP_A
ldtab$Bin <- ceiling(ldtab$Dist/1e5)

#ggplot(ldtab[sample(nrow(ldtab),10000),], aes(Dist/1e6, R2, col = type)) +
ggplot(data = ldtab, aes(x = Dist/1e6, y = R2, col = type)) +
  #geom_point(alpha = 0.01) +
  theme_bw() +
  theme(legend.position = "top") +
  scale_colour_manual(values = c("black", "red")) +
  labs(x = "Distance between SNPs (Mb)", y = "r2", col = "") +
  stat_smooth(method = 'lm', formula = y ~ log(x))

newld <- ldtab %>%
  group_by(Bin, CHR_A) %>%
  summarise(mean_r2 = mean(R2),
            SNPcount = n())

newld$type <- ifelse(newld$CHR_A %in% c(1:12, 29), "Macro-chromosomes", "Micro-chromosomes")


ggplot(newld, aes(Bin, mean_r2, col = type)) +
  geom_point(alpha = 0.1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Look at population structure           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

system("plink --bfile 70k_data\\70K_200K_maf_geno_mind_v5 --cluster --autosome-num 29 --out 70k_data\\70K_200K_maf_geno_mind_v5")

x <- read.table("70K_data/70K_200K_maf_geno_mind_v5.cluster3")
head(x)
