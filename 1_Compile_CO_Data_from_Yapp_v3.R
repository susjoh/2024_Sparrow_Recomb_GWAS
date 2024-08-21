# 
# Extract COs - QC1
# SEJ
#

# This script examines yapp after its third run (yapp_output_v3).
# QC that was done here:
# * Remove CO counts of > 45
# * Remove IDs with more than 9 COs on a single chromosome
# * Remove chromosome 26
# * Remove the DCO enriched regions on chromosomes 4 and 29 (1A)
# * Remove IDs with more than 5 short DCO events (<3Mb)

library(dplyr)
library(ggplot2)
library(tidyr)
library(kinship2)


# load("prev_data/1_Cleaned_Sparrow_Yapp_Data.RData", verbose = T)
# recchr_old <- recchr
# recsumm_old <- recsumm
# rectab_old <- rectab
# 
# rm(recchr, rectab, recsumm, ped, linkmap)

#~~ Make a vector of recombinations files

x <- dir("yapp_output_v3/")
x <- x[grep("recombinations", x)]
x <- paste0("yapp_output_v3/", x)

#~~ Load pedigree

ped <- read.table("yapp_output_v3/70K_200K_maf_geno_mind_v5.fam", stringsAsFactors = F)[,c(2:5)]
names(ped) <- c("ID", "Father", "Mother", "Sex")
ped$ParentCount <- apply(ped[,c("Father", "Mother")], 1, function(x) length(which(!x == 0)))
head(ped)

fullped <- read.table("70k_data/Helgeland_pedigree_2023-03-17.txt", header = T)

#~~ Load SNP and ID statistics

perid <- read.table("results/0_ID_Summary_Autosomal.txt", header = T)
persnp <- read.table("results/0_SNP_Summary_Autosomal.txt", header = T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Extract the CO positions               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab <- NULL

for(i in x) rectab <- rbind(rectab, read.table(i, header = T))

rectab$meiosis <- paste0(rectab$parent, "_", rectab$offspring)

rectab$CO_window <- rectab$right - rectab$left
rectab$mid_CO_pos <- rectab$left + (0.5*rectab$CO_window)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Compile the distance between crossovers #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab <- arrange(rectab, meiosis, chrom, left)
rectab$distance_to_next <- NA
rectab$mid_to_mid <- NA
rectab <- add_count(rectab, meiosis, chrom)

for(i in 1:(nrow(rectab)-1)){
  if(i %in% seq(1, nrow(rectab), 10000)) message(paste("Running row", i, "of", nrow(rectab)))
  if(rectab$meiosis[i] == rectab$meiosis[i+1] & rectab$chrom[i] == rectab$chrom[i+1]){
    rectab$distance_to_next[i] <- rectab$left[i+1] - rectab$right[i]
    rectab$mid_to_mid[i] <- rectab$mid_CO_pos[i+1] - rectab$mid_CO_pos[i]
  }
}

ggplot(rectab, aes(CO_window)) + geom_histogram(binwidth = 1e4)

#~~ remove chromosome 26 because it is dodgy AF

rectab <- subset(rectab, chrom != 26)

#~~ What is the distance between COs?

ggplot(rectab, aes(distance_to_next/1e6)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw()
ggsave("figs/1.3_Distance_to_Next_CO.png", width = 4, height = 3)

ggplot(rectab, aes(mid_to_mid/1e6)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw()
ggsave("figs/1.3_Distance_from_midCO_to_midCO.png", width = 4, height = 3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Screen for regions of high levels of double COs #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# From the previous round of investigation, there are a few regions that are
# enriched for double crossovers. These are not associated with any clear
# Mendelian errors, etc. and did not appear to grossly inflate the linkage map.
# Therefore, it may be that there are some chromosomal rearrangements in these
# regions. Chromosome 26 also seems to be very problematic so we remove those
# COs here.

ggplot(rectab, aes(distance_to_next/1e6)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw()

ggplot(rectab, aes(mid_to_mid/1e6)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  geom_vline(xintercept = 3)

#~~ regions enriched for DCOs of 3Mb or less:

x <- subset(rectab, !is.na(mid_to_mid))
x <- subset(x, mid_to_mid < 3e6)

x$mid_start <- x$mid_CO_pos
x$mid_stop <- x$mid_CO_pos+x$mid_to_mid

table(x$chrom)

#~~ Look in bins of 100kb:

binwidth = 1e5

bintab <- NULL

for(i in 1:nrow(x)){
  
  bintab <- rbind(bintab,
                  data.frame(window = floor(x$mid_start[i]/binwidth):floor(x$mid_stop[i]/binwidth),
                             chrom = x$chrom[i],
                             meiosis = x$meiosis[i]))
  
}

bintab2 <- bintab %>% group_by(chrom, window) %>% summarise(count = n())

ggplot(bintab2, aes(count)) + geom_histogram(binwidth = 1)

#~~ Some 100kb windows have very high counts. Which are they?

bintab2 <- arrange(bintab2, chrom, window)
bintab2$diff <- c(1, diff(bintab2$window))
bintab2$diff[which(bintab2$diff < 0)] <- 1
bintab2$cumu <- cumsum(bintab2$diff)

ggplot(bintab2, aes(cumu, count, group = chrom, col = factor(chrom %% 2))) +
  geom_line() +
  theme(legend.position = "none")
#ggsave("figs/1.3_Regions_enriched_for_short_DCOs.png", width = 4, height = 3)

#~~ chromosomes 4 and 29 seem to be enriched. Regions with > 100 double COs:

bintab2 <- subset(bintab2, count > 100)

badregions <- bintab2 %>% group_by(chrom) %>% summarise(start = min(window * binwidth),
                                                        stop = ((1+max(window)) * binwidth))

badregions

# get rid COs in these two regions. As they are at the end of chromosomes then
# it might be less influential.

rectab$mid_start <- rectab$mid_CO_pos
rectab$mid_stop  <- rectab$mid_CO_pos + rectab$mid_to_mid
rectab$flag <- "fine"
rectab$flag[which(rectab$chrom == 4 & rectab$left < 3e6)] <- "chr4"
rectab$flag[which(rectab$chrom == 29 & rectab$right > 68e6)] <- "chr29"

rectab <- subset(rectab, flag == "fine")

#~~ now recalculate the distance to next and mid_to_mid...

rectab <- arrange(rectab, meiosis, chrom, left)
rectab$distance_to_next <- NA
rectab$mid_to_mid <- NA
rectab$n <- NULL
rectab <- add_count(rectab, meiosis, chrom)

for(i in 1:(nrow(rectab)-1)){
  if(i %in% seq(1, nrow(rectab), 10000)) message(paste("Running row", i, "of", nrow(rectab)))
  if(rectab$meiosis[i] == rectab$meiosis[i+1] & rectab$chrom[i] == rectab$chrom[i+1]){
    rectab$distance_to_next[i] <- rectab$left[i+1] - rectab$right[i]
    rectab$mid_to_mid[i] <- rectab$mid_CO_pos[i+1] - rectab$mid_CO_pos[i]
  }
}

rm(bintab, bintab2, badregions, x, binwidth)

#~~ Lets take a look again...

ggplot(rectab, aes(distance_to_next/1e6)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw()

ggplot(rectab, aes(mid_to_mid/1e6)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  geom_vline(xintercept = 3)

ggplot(subset(rectab, chrom < 20), aes(mid_to_mid/1e6)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  geom_vline(xintercept = 3)

ggplot(rectab, aes(distance_to_next/1e6)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  labs(x = "Distance to next CO (Mb)") +
  geom_vline(xintercept = 3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Compile the yapp errors                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#yapperr <- readLines("yapp_output_v3/recombinations_rho/70K_200K_maf_geno_mind_v4_yapp.log")
#yapperr <- yapperr[grep("WARNING", yapperr)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Summarise CO Count and add all relevant variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Make a vector of recombination coverage files

x <- dir("yapp_output_v3/")
x <- x[grep("coverage", x)]
x <- paste0("yapp_output_v3/", x)

covertab <- NULL
for(i in x) covertab <- rbind(covertab, read.table(i, header = T))
covertab$meiosis <- paste0(covertab$parent, "_", covertab$offspring)
covertab$coverage <- covertab$right - covertab$left
coversumm <- covertab %>% group_by(meiosis) %>% summarise(coverage = sum(coverage))

#~~ format perid

ped$yapp_kindepth <- kindepth(ped$ID, ped$Father, ped$Mother)

names(perid)[ncol(perid)] <- "ID"
ped <- left_join(ped, perid)
names(ped)[1] <- "parent"

#~~ Summarise total CO count

recsumm <- rectab %>% group_by(meiosis, parent, offspring, sex) %>% summarise(yapp_CO_count = n())
recsumm_nomicro <- rectab %>% 
  filter(chrom %in% c(1:20,29)) %>%
  group_by(meiosis, parent, offspring, sex) %>%
  summarise(yapp_CO_no_micro = n())


#~~ Merge all the information so far

recsumm <- left_join(recsumm, recsumm_nomicro)
recsumm <- left_join(recsumm, coversumm)
recsumm <- left_join(recsumm, ped)
names(perid) <- paste0(names(perid), "_offspring")
names(perid)[ncol(perid)] <- "offspring"

recsumm <- left_join(recsumm, perid)

#~~ add matriline

fullped <- read.table("70k_data/Helgeland_pedigree_2023-03-17.txt", header = T)

fullped <- fullped[,c(1:2)]
fullped$dummysire <- 0
fullped$matriline <- ifelse(is.na(fullped$dam), fullped$id, NA)
fullped$mat_kindepth <- kindepth(fullped$id, fullped$dummysire, fullped$dam)

for(i in 1:max(fullped$mat_kindepth)){
  message(paste("Running generation", i))
  for(j in 1:nrow(fullped)){
    if(fullped$mat_kindepth[j] == i){
      fullped$matriline[j] <- fullped$matriline[which(fullped$id == fullped$dam[j])]
    }
  }
}

length(unique(fullped$matriline))
head(fullped)
fullped <- subset(fullped, select = c(id, matriline, mat_kindepth))
names(fullped)[1] <- c("parent")

recsumm <- left_join(recsumm, fullped)

# #~~ Add the old recsumm as a sanity check:
# 
# x <- subset(recsumm_old, select = c(meiosis, yapp_CO_count))
# names(x)[2] <- "yapp_CO_count_old"
# recsumm <- left_join(recsumm, x)

#~~ add parent count

x <- recsumm %>% group_by(parent) %>% summarise(parent_count = n())
recsumm <- left_join(recsumm, x)

#~~ add number of short DCOs

x <- subset(rectab, mid_to_mid < 3e6)
x <- x %>% group_by(meiosis) %>% summarise(short_DCO_count = n())

recsumm <- left_join(recsumm, x)
recsumm$short_DCO_count[which(is.na(recsumm$short_DCO_count))] <- 0

#~~ How does it look?

ggplot(recsumm, aes(yapp_CO_count)) + geom_histogram(binwidth = 1)
ggplot(recsumm, aes(yapp_CO_no_micro)) + geom_histogram(binwidth = 1)

#~~ At this point just remove the badly phased meioses

recsumm <- subset(recsumm, yapp_CO_count < 45)
rectab <- subset(rectab, meiosis %in% recsumm$meiosis)

recsumm <- subset(recsumm, !meiosis %in% rectab$meiosis[which(rectab$n > 9)])
rectab <- subset(rectab, meiosis %in% recsumm$meiosis)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Check associations with different factors #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# #~~ Does it match with the old yapp CO count?
# 
# ggplot(recsumm, aes(yapp_CO_count, yapp_CO_count_old)) +
#   geom_point(alpha = 0.1) +
#   stat_smooth(method = "lm")
# 
# cor.test(recsumm$yapp_CO_count, recsumm$yapp_CO_count_old)

#~~ More COs are present when heterozygosity is higher... has this effect disappeared?

ggplot(recsumm, aes(F, yapp_CO_count)) + 
  geom_point() +
  stat_smooth()

ggplot(recsumm, aes(Het, yapp_CO_count)) + 
  geom_point(alpha = 0.1) +
  stat_smooth()

ggplot(recsumm, aes(CallPP, yapp_CO_count)) + 
  geom_point(alpha = 0.01) +
  stat_smooth()

ggplot(recsumm, aes(short_DCO_count, yapp_CO_count)) + 
  geom_point(alpha = 0.01) +
  stat_smooth()

ggplot(recsumm, aes(short_DCO_count)) + geom_histogram(binwidth = 1)

table(recsumm$short_DCO_count)

#~~ Remove IDs with > 5 short DCOs

recsumm <- subset(recsumm, short_DCO_count <= 5)
rectab <- subset(rectab, meiosis %in% recsumm$meiosis)

table(rectab$n)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Variation between chromosomes        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

recchr <- left_join(expand.grid(meiosis = unique(rectab$meiosis),
                                chrom = unique(rectab$chrom)),
                    rectab %>% group_by(meiosis, chrom) %>% count())
recchr$n[which(is.na(recchr$n))] <- 0

recchr <- arrange(recchr, meiosis, chrom)
recchr$n <- factor(recchr$n, levels = rev(sort(unique(recchr$n))))
recchr$chrom <- factor(recchr$chrom)
nmei <- length(unique(recchr$meiosis))

head(recchr)
recchr <- subset(recchr, n != "9") %>% droplevels

ggplot(recchr, aes(chrom, fill = n)) +
  geom_bar() +
  scale_fill_brewer(palette = "GnBu") +
  labs(x = "Chromosome", y = "Proportion of Meioses") +
  scale_y_continuous(breaks = seq(0, nmei, nmei/5), labels = seq(0, 1, 0.2))

ggsave("figs/1.3_Distribution_of_COs_Per_Chromosome.png", width = 6, height = 4)

#~~ AT THIS STAGE, WRITE TO FILE to deal with DCOs and final compiling separately.

write.table(recsumm, "results/1.3_recsumm_COs_per_meiosis_QC1.txt", row.names = F, quote = F, sep = "\t")
write.table(rectab, "results/1.3_rectab_COs_per_meiosis_per_chr_QC1.txt", row.names = F, quote = F, sep = "\t")
write.table(recchr, "results/1.3_recchr_COs_per_chr_QC1.txt", row.names = F, quote = F, sep = "\t")
write.table(ped, "results/1.3_pedigree_used_for_phasing_v5_QC1.txt", row.names = F, quote = F, sep = "\t")

save(recsumm, rectab, recchr, ped, file = "results/1.3_CO_data_after_QC1.RData")
