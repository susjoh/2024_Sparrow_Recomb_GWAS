# 
# Extract COs - QC1
# SEJ
#

# This script examines yapp after its second run (yapp_output_v2).

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

x <- dir("yapp_output_v2/")
x <- x[grep("recombinations", x)]
x <- paste0("yapp_output_v2/", x)

#~~ Load pedigree

ped <- read.table("yapp_output_v2/70K_200K_maf_geno_mind_v5.fam", stringsAsFactors = F)[,c(2:5)]
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

#~~ What is the distance between COs?

ggplot(rectab, aes(distance_to_next/1e6)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  facet_wrap(~chrom, scales = "free")

ggplot(rectab, aes(mid_to_mid/1e6)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  facet_wrap(~chrom, scales = "free")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Screen for regions of high levels of double COs #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# From the previous round of investigation, there are a few regions that are
# enriched for double crossovers. These are not associated with any clear
# Mendelian errors, etc. and did not appear to grossly inflate the linkage map.
# Therefore, it may be that there are some chromosomal rearrangements in these
# regions. Chromosome 26 also seems to be very problematic so we remove those
# COs here.

rectab <- subset(rectab, chrom != 26)

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

bintab2$chrom[which(bintab2$chrom == 29)] <- 1.5
bintab2 <- arrange(bintab2, chrom, window)
bintab2$diff <- c(1, diff(bintab2$window))
bintab2$diff[which(bintab2$diff < 0)] <- 1
bintab2$cumu <- cumsum(bintab2$diff)
bintab2$chrom2 <- bintab2$chrom + 1
bintab2$chrom2[which(bintab2$chrom2 == 2)] <- 1
bintab2$chrom2[which(bintab2$chrom2 == 2.5)] <- 2


chrtab <- bintab2 %>% group_by(chrom) %>% summarise(cumu = mean(cumu))
chrtab$chrom[2] <- "1A"
chrtab$chrom2 <- c(1, "1A", 2:10, "", "", 13, "", "", 17, "", "", 20, "", "", "", 24, "", "", 28)

ggplot(subset(bintab2, chrom != 26), aes(cumu, count, group = chrom, col = factor(chrom2 %% 2))) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Chromosome") +
  scale_x_continuous(breaks = chrtab$cumu, labels = chrtab$chrom2)
ggsave("figs/1.3_Regions_enriched_for_short_DCOs.png", width = 6, height = 4)

#~~ chromosomes 4 and 29 seem to be enriched. Regions with > 100 double COs:

bintab2 <- subset(bintab2, count > 100)

badregions <- bintab2 %>% group_by(chrom) %>% summarise(start = min(window * binwidth),
                                                        stop = ((1+max(window)) * binwidth))

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Compile the yapp errors                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#yapperr <- readLines("yapp_output_v2/recombinations_rho/70K_200K_maf_geno_mind_v4_yapp.log")
#yapperr <- yapperr[grep("WARNING", yapperr)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Summarise CO Count and add all relevant variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Make a vector of recombination coverage files

x <- dir("yapp_output_v2/")
x <- x[grep("coverage", x)]
x <- paste0("yapp_output_v2/", x)

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

#~~ Add the old recsumm as a sanity check:

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

#~~ More COs are present when heterozygosity is higher...

ggplot(recsumm, aes(F, yapp_CO_count)) + 
  geom_point() +
  stat_smooth()

ggplot(recsumm, aes(Het, yapp_CO_count)) + 
  geom_point(alpha = 0.1) +
  stat_smooth(method = "lm")

ggplot(recsumm, aes(CallPP, yapp_CO_count)) + 
  geom_point(alpha = 0.01) +
  stat_smooth(method = "lm")

ggplot(recsumm, aes(short_DCO_count, yapp_CO_count)) + 
  geom_point(alpha = 0.01) +
  stat_smooth()

ggplot(recsumm, aes(CallPP, Het)) + 
  geom_point(alpha = 0.01) +
  stat_smooth()

ggplot(perid, aes(CallPP_offspring)) +
  geom_histogram(binwidth = 0.001) +
  geom_vline(xintercept = 0.99)

ggplot(perid, aes(Het_offspring)) +
  geom_histogram(binwidth = 0.001) +
  geom_vline(xintercept = (mean(recsumm$Het) - (3*(sd(recsumm$Het))))) +
  geom_vline(xintercept = (mean(recsumm$Het) + (3*(sd(recsumm$Het)))))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Remove bad IDs for the next round         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Keep only IDs that have a call rate of > 0.98
#~~ Keep only IDs that have a Het of fewer than 3 sds of the mean
#~~ Remove parents from meioses that have a CO count of > 45
#~~ Remove parents from meioses that have more than 10 COs on a chromosome

badids <- rbind(subset(perid, Het_offspring < mean(recsumm$Het_offspring) - (3*(sd(recsumm$Het_offspring)))),
                subset(perid, Het_offspring > mean(recsumm$Het_offspring) + (3*(sd(recsumm$Het_offspring)))),
                subset(perid, CallPP_offspring < 0.98))$offspring %>% unique

badparents <- c(subset(recsumm, yapp_CO_count > 45)$parent)
table(badparents)

badoffspring <- c(subset(recsumm, yapp_CO_count > 45)$offspring,
                  subset(rectab, n > 10)$offspring) %>% unique

#~~ Create new famfile for run v2

new_famfile <- read.table("yapp_output_v2/70K_200K_maf_geno_mind_v5.fam")
new_famfile %>% head

new_famfile$V3[which(new_famfile$V3 %in% badids)] <- 0
new_famfile$V4[which(new_famfile$V4 %in% badids)] <- 0

new_famfile$V3[which(new_famfile$V2 %in% badoffspring)] <- 0
new_famfile$V4[which(new_famfile$V2 %in% badoffspring)] <- 0

write.table(new_famfile, "yapp_output_v3/70K_200K_maf_geno_mind_v5.fam", row.names = F, col.names = F, quote = F)

# old_famfile <- read.table("yapp_output_v2/70K_200K_maf_geno_mind_v5.fam")
# 
# table(old_famfile$V2 == new_famfile$V2)


