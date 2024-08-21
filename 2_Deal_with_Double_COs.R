#
# Yapp sparrow sanity check DCO
# JBM & SEJ
# July 2022, now June 2023
#

# Please note: this script originally investigated the raw phasing from the vcf
# files and phases module in yapp. AFTER DISCUSSING WITH BERTRAND ON 03/12/2022:
# the output of the recombination module should be used as the phasing is not
# associated with the decision of the HMM. The HMM may call false crossovers in
# regions of error and/or larger deletions that are uncharacterised. These can
# be removed with a cut-off as we do below. 

library(dplyr)
library(tidyr)
library(ggplot2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Load and format data            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Read in yapp & Crimap results

load("results/1.3_CO_data_after_QC1.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Visualise the number of COs     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~ Number of COs per ID

ggplot(recsumm, aes(yapp_CO_count)) +
  geom_histogram(binwidth = 1, col = "grey50")

#~~ Number of COs per chromosome

ggplot(subset(recchr, n != "9"), aes(chrom, fill = n)) +
  geom_bar() +
  scale_fill_brewer(palette = "GnBu") +
  labs(x = "Chromosome", y = "Proportion of Meioses") +
  scale_y_continuous(breaks = seq(0, nrow(recsumm), nrow(recsumm)/5), labels = seq(0, 1, 0.2))

#~~ What is the distance between COs?

ggplot(rectab, aes(distance_to_next/1e6)) +
  geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = 2, linetype = "dashed") +
  theme_bw() +
  labs(x = "Distance to next CO (Mb)") +
  ggtitle("B:")
ggsave("figs/2_distance_to_next_CO.png", width = 4, height = 4)

ggplot(rectab, aes(CO_window/1e6)) +
  geom_histogram(binwidth = 0.01) +
  theme_bw() +
  labs(x = "CO window size (Mb)") +
  coord_cartesian(xlim = c(0, 5))+
  ggtitle("A:")
ggsave("figs/2_CO_Window_Size.png", width = 4, height = 4)

ggplot(rectab, aes(mid_to_mid)) +
  geom_histogram(binwidth = 5e5) +
  geom_vline(xintercept = 4e6)

ggplot(rectab, aes(distance_to_next, mid_to_mid)) +
  geom_point(alpha = 0.2) +
  geom_vline(xintercept = 2e6) +
  geom_hline(yintercept = 4e6)

ggplot(rectab, aes(distance_to_next/1e6)) +
  geom_histogram(binwidth = 0.2) +
  facet_wrap(~chrom, scales = "free") +
  geom_vline(xintercept = 2)

ggplot(rectab, aes(mid_to_mid/1e6)) +
  geom_histogram(binwidth = 0.2) +
  facet_wrap(~chrom, scales = "free") +
  geom_vline(xintercept = 4)

ggplot(subset(rectab, chrom == 2), aes(distance_to_next/1e6)) +
  geom_histogram(binwidth = 0.2) +
  geom_vline(xintercept = 2)

ggplot(subset(rectab, chrom == 2), aes(mid_to_mid/1e6)) +
  geom_histogram(binwidth = 0.2) +
  geom_vline(xintercept = 4)

ggplot(rectab, aes(log10(distance_to_next))) +
  geom_histogram(binwidth = 0.2)

ggplot(subset(rectab, chrom %in% c(1:20, 29)), aes(distance_to_next/1e6)) +
  geom_histogram(binwidth = 0.2) +
  geom_vline(xintercept = 2, col = "red")

ggplot(subset(rectab, chrom %in% c(1:20, 29)), aes(log10(distance_to_next))) +
  geom_histogram(binwidth = 0.2) +
  geom_vline(xintercept = log10(2000000))

ggplot(rectab, aes(log10(distance_to_next))) +
  geom_histogram(binwidth = 0.025) +
  facet_wrap(~chrom, scales = "free")

# plot(ecdf(rectab$distance_to_next))
# plot(ecdf(log10(rectab$distance_to_next)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Get the SNP information for each position #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# library(vcfR)
# 
# vcfsnps <- NULL
# 
# for(chr in sort(unique(rectab$chrom))){
#   x <- read.vcfR(paste0("yapp_output_v3/70K_200K_maf_geno_mind_v5_", chr,"_phased.vcf.gz"))
#   x <- data.frame(x@fix)
#   x$ORDER <- 1:nrow(x)
#   vcfsnps <- rbind(vcfsnps, x)
# }
# 
# str(vcfsnps)
# vcfsnps$POS <- as.numeric(vcfsnps$POS)
# vcfsnps$CHROM <- as.numeric(vcfsnps$CHROM)
# 
# vcfsumm <- vcfsnps %>%
#   group_by(CHROM) %>%
#   summarise(nsnps = n(),
#             chr_len = max(POS))
# 
# write.table(vcfsnps, "results/2_VCF_SNP_information_v5.txt", row.names = F, sep = "\t", quote = F)
# 
# rm(x, ped, linkmap, chr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Do QC                                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# USE DISTANCE TO NEXT

# Flag the SNPs on either side of a run of less than 2e6

cutoff <- 2e6

rectab <- arrange(rectab, meiosis, chrom, left)

rectab$shortCO <- "no"
rectab$shortCO[which(rectab$distance_to_next < cutoff)] <- "yes"
rectab$shortCO[which(rectab$shortCO  == "yes")+1] <- "yes"

rectab$n <- NULL

#~~ Dodgy CO run length

x <- table(rectab$meiosis, rectab$chrom, rectab$shortCO) %>% data.frame()
head(x)
names(x) <- c("meiosis", "chrom", "x", "shortCO_count")
x <- subset(x, x == "yes")
x$x <- NULL
x$meiosis <- as.character(x$meiosis)
x$chrom <- as.numeric(as.character(x$chrom))

rectab <- left_join(rectab, x)
rm(x)

# After thinking about this - if there are only 2 short COs on the chromosome,
# then they are related and both can be removed. If there are > 2 short COs on a
# chromosome, they need to be classified into runs. If there are e.g. 3 adjacent
# COs that are preceded or succeeded by a long CO, then those 3 COs are 1 CO.
# Need to check that 3COs are always adjacent, as it could be separated and they
# are independent....


rectab_hold <- rectab
rectab <- rectab_hold

rectab$KEY <- paste(rectab$meiosis, rectab$chrom, sep = "_")
rectab <- add_count(rectab, KEY)

#~~ Make a rectab_v2 that has all the good COs in it.

rectab_v2 <- subset(rectab, shortCO == "no")
rectab_fix <- subset(rectab, shortCO == "yes")

#write.table(rectab_fix, "test.txt", row.names = F, sep = "\t", quote = F)

ggplot(rectab_fix, aes(n)) + geom_histogram()

# If there are only 2 COs in rectab_fix, then just remove them. What about 4
# COs....?? If the distance to next in line 2 is greater than the cutoff, then
# it is two separate double COs and they can both be removed. But also if it is
# 4 sequential phase changes, they won't change the phase! So actually... all
# even numbered multiple COs be removed. In this case, all with 6 COs were also
# short DCos with no phase change.

rectab_fix <- subset(rectab_fix, !shortCO_count %in% c(2, 4, 6))

# Triple crossovers will result in a phase change. They are always sequential, so
# need to take the left of the first CO and the right of the third.

for(i in unique(rectab_fix$KEY)){
  
  x <- subset(rectab_fix, KEY == i)
  
  if(x$shortCO_count[1] == 3){
    
    x$right[1] <- x$right[3]
    x <- x[1,]
    x[,c(8:13, 15)] <- NA
    rectab_v2 <- rbind(rectab_v2, x)
    rectab_fix <- rectab_fix[-which(rectab_fix$KEY == i),]
    
  }
  rm(x)
}

#~~ what's left?

# These IDs have trouble phasing over 5 COs.

for(i in c("8730974__8730973_M2_8730561_K17_2",
           "8732254_8168605_2",
           "8M31744_F13_8N72121_J11_9",
           "8N72121_J11_8N13256_P16_9")){
  x <- subset(rectab_fix, KEY == i)
  x$right[1] <- x$right[5]
  x <- x[1,]
  x[,c(8:13, 15)] <- NA
  rectab_v2 <- rbind(rectab_v2, x)
  rectab_fix <- rectab_fix[-which(rectab_fix$KEY == i),]
}

#~~ The rest are a combination of 2 & 3 COs. I hand edited this (SEJ) as a
#little time constrained. It is only 4 chromatids.

rectab_fix <- read.table("2_rectab_DCO_hand_edit.txt", header = T)

rectab_fix <- subset(rectab_fix, !shortCO_count %in% c(2, 4, 6))

for(i in unique(rectab_fix$KEY)){
  
  x <- subset(rectab_fix, KEY == i)
  
  if(x$shortCO_count[1] == 3){
    
    x$right[1] <- x$right[3]
    x <- x[1,]
    x[,c(8:13, 15)] <- NA
    rectab_v2 <- rbind(rectab_v2, x)
    rectab_fix <- rectab_fix[-which(rectab_fix$KEY == i),]
    
  }
  rm(x)
}


#~~ RECALCULATE DISTANCES, WINDOWS etc

head(rectab_v2)
rectab_v2 <- arrange(rectab_v2, KEY, left)
rectab_v2[,c(8:13, 15)] <- NA
rectab_v2$CO_window <- rectab_v2$right - rectab_v2$left
rectab_v2$mid_CO_pos <- rectab_v2$left + (0.5*rectab_v2$CO_window)

for(i in 1:(nrow(rectab_v2)-1)){
  if(rectab_v2$KEY[i] == rectab_v2$KEY[i+1]){
    rectab_v2$distance_to_next[i] <- rectab_v2$left[i+1] - rectab_v2$right[i]
    rectab_v2$mid_to_mid[i] <- rectab_v2$mid_CO_pos[i+1] - rectab_v2$mid_CO_pos[i]
  }
}

head(rectab_v2)
rectab_v2 <- subset(rectab_v2, select = -c(mid_start, mid_stop, flag, shortCO, shortCO_count, n))

rectab_v2 <- add_count(rectab_v2, KEY)


rectab_v2$meiosis <- paste0(rectab_v2$parent, "_", rectab_v2$offspring)

#~~ all looks good

rm(rectab, rectab_hold, i, rectab_fix)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Create new response variables             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab <- rectab_v2

# Look at the distribution of COs per chromosome

recchr <- left_join(expand.grid(meiosis = unique(rectab$meiosis),
                                chrom = unique(rectab$chrom)),
                    rectab %>% group_by(meiosis, chrom) %>% count())
recchr$n[which(is.na(recchr$n))] <- 0

recchr <- arrange(recchr, meiosis, chrom)
recchr$n <- factor(recchr$n, levels = rev(sort(unique(recchr$n))))
recchr$chrom <- factor(recchr$chrom)
nmei <- length(unique(recchr$meiosis))

ggplot(recchr, aes(chrom, fill = n)) +
  geom_bar() +
  scale_fill_brewer(palette = "GnBu") +
  labs(x = "Chromosome", y = "Proportion of Meioses") +
  scale_y_continuous(breaks = seq(0, nmei, nmei/5), labels = seq(0, 1, 0.2))

ggsave("figs/2_Distribution_of_COs_Per_Chromosome_QCed.png", width = 6, height = 4)


#~~ Get the new CO count information

newsumm <- rectab_v2 %>%
  group_by(meiosis) %>%
  summarise(yapp_CO_count_QCed = n())

recsumm <- left_join(recsumm, newsumm)
head(recsumm)

ggplot(recsumm, aes(yapp_CO_count, yapp_CO_count_QCed)) +
  geom_point(alpha = 0.05)

ggplot(recsumm, aes(factor(Sex), yapp_CO_count)) +
  geom_boxplot()

ggplot(recsumm, aes(yapp_CO_count_QCed)) +
  geom_histogram(binwidth = 1, col = "grey") +
  facet_wrap(~Sex)

ggplot(recsumm, aes(yapp_CO_count)) +
  geom_histogram(binwidth = 1, col = "grey")

#~~ What about micro/macro

newsumm <- subset(rectab_v2, !chrom %in% c(21:28)) %>%
  group_by(meiosis) %>%
  summarise(yapp_CO_count_no_micro = n())
recsumm <- left_join(recsumm, newsumm)

newsumm <- subset(rectab_v2, chrom %in% c(1:12, 29)) %>%
  group_by(meiosis) %>%
  summarise(yapp_CO_count_macro = n())
recsumm <- left_join(recsumm, newsumm)

ggplot(recsumm, aes(yapp_CO_count_no_micro)) +
  geom_histogram(binwidth = 1, col = "grey")

ggplot(recsumm, aes(yapp_CO_count_macro)) +
  geom_histogram(binwidth = 1, col = "grey")

recchr_v2 <- recchr

recsumm <- data.frame(recsumm)

# ggplot(recsumm, aes(yapp_CO_count_old, yapp_CO_count_QCed)) + geom_point(alpha = 0.2)
# 
# cor.test(recsumm$yapp_CO_count_old, recsumm$yapp_CO_count_QCed)

recchr_v2$chrom <- as.numeric(as.character(recchr_v2$chrom))
recchr_v2$n <- as.numeric(as.character(recchr_v2$n))

str(recchr_v2)

#save(rectab_v2, recchr_v2, recsumm, file = "2_Full_QC_Sparrow_Yapp_Data.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Addendum - add coverage data         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#load("2_Full_QC_Sparrow_Yapp_Data.RData")

covertab <- NULL

for(i in unique(rectab_v2$chrom)){
  covertab <- rbind(covertab,
                    read.table(paste0("yapp_output_v3/70K_200K_maf_geno_mind_v5_", i , "_yapp_recomb_coverage.txt"), stringsAsFactors = F, header = T))
  
}

covertab$meiosis <- paste0(covertab$parent, "_", covertab$offspring)
head(covertab)
covertab$coverage <- covertab$right - covertab$left

str(recchr_v2)

recchr_v2 <- left_join(recchr_v2, covertab)

#~~ Get total coverage

x <- covertab %>%
  group_by(meiosis) %>%
  summarise(Total_Coverage = sum(coverage))

recsumm <- left_join(recsumm, x)

save(rectab_v2, recchr_v2, recsumm,
     file = "2_Full_QC_Sparrow_Yapp_Data.RData")
