# 
# Extract COs - RAW
# SEJ
#

# This script examines yapp after its first run (yapp_output). At the end some
# changes are made to the family structure and then it is rerun.

library(dplyr)
library(ggplot2)
library(tidyr)
library(kinship2)


load("prev_data/1_Cleaned_Sparrow_Yapp_Data.RData", verbose = T)
recchr_old <- recchr
recsumm_old <- recsumm
rectab_old <- rectab

rm(recchr, rectab, recsumm, ped, linkmap)

#~~ Make a vector of recombinations files

x <- dir("yapp_output/")
x <- x[grep("recombinations", x)]
x <- paste0("yapp_output/", x)

#~~ Load pedigree

ped <- read.table("yapp_output/70K_200K_maf_geno_mind_v5.fam", stringsAsFactors = F)[,c(2:5)]
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
rectab$dodgy_chromosome <- NULL
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Compile the yapp errors                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#yapperr <- readLines("yapp_output/recombinations_rho/70K_200K_maf_geno_mind_v4_yapp.log")
#yapperr <- yapperr[grep("WARNING", yapperr)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Summarise CO Count                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Make a vector of recombination coverage files

x <- dir("yapp_output/")
x <- x[grep("coverage", x)]
x <- paste0("yapp_output/", x)

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


#~~ add info

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

ggplot(recsumm, aes(yapp_CO_count)) + geom_histogram(binwidth = 1)
ggplot(recsumm, aes(yapp_CO_no_micro)) + geom_histogram(binwidth = 1)

x <- recsumm %>% group_by(parent) %>% summarise(parent_count = n())
recsumm <- left_join(recsumm, x)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Check associations with different factors #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

names(recsumm)
#bad_recsumm <- rbind(subset(recsumm, yapp_CO_count > 60))

x <- subset(recsumm_old, select = c(meiosis, yapp_CO_count))
names(x)[2] <- "yapp_CO_count_old"
recsumm <- left_join(recsumm, x)

ggplot(subset(recsumm, yapp_CO_count < 50), aes(yapp_CO_count, yapp_CO_count_old)) +
  geom_point(alpha = 0.1) +
  stat_smooth(method = "lm")

cor.test(recsumm$yapp_CO_count, recsumm$yapp_CO_count_old)

#recsumm <- subset(recsumm, yapp_CO_count < 61)

ggplot(recsumm, aes(F, yapp_CO_count)) + 
  geom_point() +
  stat_smooth()

ggplot(recsumm, aes(Het, yapp_CO_count)) + 
  geom_point() +
  stat_smooth() +
  geom_vline(xintercept = 0.36)

ggplot(subset(recsumm, Het < 0.36 & yapp_CO_count < 60), aes(Het, yapp_CO_count)) + 
  geom_point() +
  stat_smooth()

badparents <- subset(recsumm, Het > 0.36)$parent %>% unique

#~~ Create new famfile for run v2

new_famfile <- read.table("yapp_output/70K_200K_maf_geno_mind_v5.fam.new")
new_famfile %>% head


recsumm <- filter(recsumm, Het < 0.36)

badoffspring <- subset(recsumm, yapp_CO_count > 60)$offspring %>% unique

#~~ Compile the IDs with t00 many COs on a chromosome

x <- subset(rectab, n > 9)
x$meiosis %>% unique
x <- subset(x, !parent %in% badparents)
x <- subset(x, !offspring %in% badoffspring)
x <- unique(subset(x, select = c(meiosis, parent, offspring)))
table(x$parent) %>% sort
table(x$offspring) %>% sort

badoffspring <- c(badoffspring, unique(x$offspring))

new_famfile$V3[which(new_famfile$V2 %in% badoffspring)] <- 0
new_famfile$V4[which(new_famfile$V2 %in% badoffspring)] <- 0

new_famfile$V3[which(new_famfile$V3 %in% badparents)] <- 0
new_famfile$V4[which(new_famfile$V4 %in% badparents)] <- 0

head(new_famfile)

write.table(new_famfile, "yapp_output_v2/70K_200K_maf_geno_mind_v5.fam", row.names = F, col.names = F, quote = F)

# old_famfile <- read.table("yapp_output/70K_200K_maf_geno_mind_v5.fam")
# 
# table(old_famfile$V2 == new_famfile$V2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. regions enriched for double COs? # JUST EXPLORATION AT PRESENT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab$meiosischr <- paste0(rectab$meiosis)

x <- subset(rectab, !is.na(mid_to_mid))
x <- subset(x, mid_to_mid < 5e6)
head(x)

table(x$chrom)

#~~ bin enrichment...

x$mid_start <- x$mid_CO_pos
x$mid_stop <- x$mid_CO_pos+x$mid_to_mid

x <- subset(x, chrom != 26)

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

bintab2 <- arrange(bintab2, chrom, window)
bintab2$diff <- c(1, diff(bintab2$window))
bintab2$diff[which(bintab2$diff < 0)] <- 1
bintab2$cumu <- cumsum(bintab2$diff)

ggplot(bintab2, aes(cumu, count, group = chrom, col = factor(chrom %% 2))) + geom_line()

#~~ regions with > 100 double COs 

bintab2 <- subset(bintab2, count > 100)

badregions <- bintab2 %>% group_by(chrom) %>% summarise(start = min(window * binwidth),
                                                        stop = ((1+max(window)) * binwidth))

snptab <- read.table("results/0_SNP_Summary_Autosomal.txt", header = T)

snptab$flag <- "fine"
snptab$Diff <- c(1, diff(snptab$Position))
snptab$Diff[which(snptab$Diff < 0)] <- 1e4
snptab$Cumu <- cumsum(snptab$Diff)

for(i in 1:nrow(badregions)){
  snptab$flag[which(snptab$Chromosome == badregions$chrom[i] & snptab$Position > badregions$start[i] & snptab$Position < badregions$stop)] <- "dodgy"
}

ggplot(snptab, aes(Q.2, Fmax, col = flag)) + geom_point()

ggplot(snptab, aes(Cumu, Fmax, col = Chromosome %% 2, shape = flag)) + geom_point()

# system("plink --bfile 70k_data/70K_200K_maf_geno_mind_v5 --mendel --autosome-num 32 --out 70k_data/70K_200K_maf_geno_mind_v5")
# 
snpmend <- read.table("70k_data/70K_200K_maf_geno_mind_v5.lmendel", header = T)

snptab <- left_join(snptab, snpmend)

ggplot(snptab, aes(Cumu, N, col = Chromosome %% 2)) + geom_point()

head(snpmend)

table(snptab$flag)

temp <- subset(rectab, chrom == 4)

ggplot(temp, aes(mid_CO_pos)) + geom_histogram(binwidth = 10000)

badregions
write.table(badregions, "results/1_High_DoubleCO_Regions.txt", row.names = F, sep = "\t", quote = F)


temp <- subset(rectab, chrom == 2)
ggplot(temp, aes(mid_CO_pos)) + geom_histogram(binwidth = 10000)

rectab$mid_start <- rectab$mid_CO_pos
rectab$mid_stop  <- rectab$mid_CO_pos + rectab$mid_to_mid
rectab$flag <- "fine"
rectab$flag[which(rectab$chrom == 4 & rectab$left < 3e6)] <- "chr4"
rectab$flag[which(rectab$chrom == 26)] <- "chr26"
rectab$flag[which(rectab$chrom == 29 & rectab$right > 68e6)] <- "chr29"
rectab$flag[which(rectab$chrom == 2 & rectab$mid_start > 64e6 & rectab$mid_start < 73e6 & rectab$mid_stop > 64e6 & rectab$mid_stop < 73e6)] <- "chr2_dco"

rectab <- arrange(rectab, meiosis, chrom, mid_CO_pos)
rectab$flag[which(rectab$flag == "chr2_dco")+1] <- "chr2_dco"

table(rectab$flag)

#~~ REMOVE the flags!!!

rectab <- subset(rectab, flag == "fine")

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




#~~ Remove the first 3 million bases from chromosome 4
#~~ Remove the last 2 million bases from chromosome 29 (1A)



