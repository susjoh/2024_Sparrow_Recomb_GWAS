#
# Intra-chromosomal shuffling
# JBM & SEJ
#
#

library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

load("2_Full_QC_Sparrow_Yapp_Data.RData")
spargff <- read.table("prev_data/sparrow_genes_from_homa.gff", stringsAsFactors = F, sep = "\t")
maptab <- read.table("70k_data/70K_200K_maf_geno_mind_v5.bim")
names(maptab) <- c("Chr", "SNP", "cM", "Position", "ref1", "ref2")

#~~ Tidy up workspace and do minor edits

recchr <- recchr_v2
rectab <- rectab_v2

rm(recchr_v2, rectab_v2)

linksumm <- maptab %>%
  group_by(Chr) %>%
  summarise(Chr_length = max(Position)) %>% data.frame()
names(linksumm)[1] <- "chrom"

rectab <- left_join(rectab, linksumm)

names(recchr)[7:8] <- c("left_coverage", "right_coverage")

# recsumm is the summarised data per id
# recchr is the number of COs and phasing coverage of each chromosome
head(recchr)
# rectab is the position of each CO
head(rectab)

table(is.na(rectab$meiosis))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Get Phase information        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rectab <- arrange(rectab, meiosis, chrom, left)

# Give each CO an order (1st, 2nd, 3rd)

rectab$CO_Order <- 1

for(i in 2:nrow(rectab)){
  if(i %in% seq(1, nrow(rectab), 10000)) message(paste("Running line", i, "of", nrow(rectab)))
  rectab$CO_Order[i] <- ifelse(rectab$KEY[i] != rectab$KEY[i-1], 1, rectab$CO_Order[i-1] + 1)
}

# Add total number of COs on the chr and then the coverage data

rectab$n <- NULL
rectab <- add_count(rectab, KEY)
rectab <- left_join(rectab, subset(recchr, select = c(meiosis, chrom, left_coverage, right_coverage)))

#~~ characterise each chromosome fragment (2COs = 3 fragments)

chunktab <- list()
counter <- 1

for(i in 1:nrow(rectab)){
  
  if(i %in% seq(1, nrow(rectab), 1000)) message(paste("Running line", i, "of", nrow(rectab)))
  x <- NULL
  
  #~~ if it is the first chunk
  
  if(rectab$CO_Order[i] == 1){
    x <- rbind(x,
               data.frame(start = 1,
                          start_true = rectab$left_coverage[i],
                          stop = rectab$mid_CO_pos[i],
                          stop_true = rectab$left[i],
                          chunk_order = 1)
    )
  }
  
  if(any(x$stop < x$start)) stop(print(i))
  
  #~~ if it is left of a mid chunk
  
  if(rectab$CO_Order[i] != rectab$n[i]){
    x <- rbind(x, data.frame(start = rectab$mid_CO_pos[i]+1,
                             start_true = rectab$right[i],
                             stop = rectab$mid_CO_pos[i+1],
                             stop_true = rectab$left[i+1],
                             chunk_order = rectab$CO_Order[i] + 1)
    )
  }
  
  if(any(x$stop < x$start)) stop(print(i))
  
  #~~ if it is the last chunk
  
  if(rectab$CO_Order[i] == rectab$n[i]){
    x <- rbind(x, data.frame(start = rectab$mid_CO_pos[i]+1,
                             start_true = rectab$right[i],
                             stop = rectab$Chr_length[i],
                             stop_true = rectab$right_coverage[i],
                             chunk_order = rectab$CO_Order[i] + 1)
    )
  }
  
  if(any(x$stop < x$start)) stop(print(i))
  
  x$KEY <- rectab$KEY[i]
  
  chunktab[[counter]] <- x
  counter <- counter + 1
  rm(x)
}

chunktab <- bind_rows(chunktab)

#~~ Do we get the expected number of chunks??
(sum(recchr$n) + nrow(subset(recchr, n != 0))) == nrow(chunktab)

# WOOP WOOP

#~~ Finally, add some phase information. NB THIS DOES NOT TRANSLATE ACROSS
# CHROMOSOMES! We will always call odd chunks 1 and even chunks 0.

head(chunktab)
chunktab$Phase <- chunktab$chunk_order %% 2

#~~ Tidy up table

x <- subset(rectab, select = c(parent, sex, offspring, chrom, meiosis, KEY, n)) %>% unique
x$chunk_count = x$n+1
x$n <- NULL

chunktab <- left_join(chunktab, x)

chunktab <- chunktab[,c("parent", "sex", "offspring", "chrom", "meiosis", "KEY",  "start",
                        "start_true", "stop", "stop_true", "chunk_order",  "chunk_count", "Phase")]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. How much of the chromosome?  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Get the chunk lengths

chunktab$chunk_length <- chunktab$stop - chunktab$start
chunktab$chunk_length_true <- chunktab$stop_true - chunktab$start_true

#~~ Get chr lengths and chunk proportions

x1 <- chunktab %>%
  group_by(KEY) %>%
  summarise(chrom_len = sum(chunk_length),
            coverage_len = sum(chunk_length_true))

x2 <- chunktab %>%
  group_by(KEY, Phase) %>%
  summarise(chrom_len_phase = sum(chunk_length),
            coverage_len_phase = sum(chunk_length_true))

x2 <- melt(x2, id.vars = c("KEY", "Phase"))
x2$variable <- paste0(x2$variable, "_", x2$Phase)
x2$Phase <- NULL

x2 <- dcast(x2,formula = KEY ~ variable)
x2 <- left_join(x2, x1)

head(x2)

#~~ Add to recchr
recchr$KEY <- paste0(recchr$meiosis, "_", recchr$chrom)
shufftab <- left_join(recchr, x2)
shufftab <- left_join(shufftab, linksumm)
head(shufftab)

na.lines <- which(is.na(shufftab$chrom_len))
shufftab$chrom_len[na.lines] <- shufftab$Chr_length[na.lines] 
shufftab$chrom_len_phase_0[na.lines] <- shufftab$chrom_len[na.lines] 
shufftab$chrom_len_phase_1[na.lines] <- 0 
shufftab$coverage_len_phase_0[na.lines] <- shufftab$coverage[na.lines] 
shufftab$coverage_len_phase_1[na.lines] <- 0 

rm(x, x1, x2, counter, i, na.lines)

head(shufftab)
shufftab$test <- shufftab$coverage_len_phase_0 + shufftab$coverage_len_phase_1

plot(test ~ coverage_len, data = shufftab)

shufftab <- subset(shufftab, select = -c(Chr_length, chrom_len_phase_1, coverage_len_phase_1, left_coverage, right_coverage))

#~~ Now add proportion of genome

genome_len <- sum(linksumm$Chr_length)
genome_len_macro <- sum(subset(linksumm, chrom %in% c(1:12, 29))$Chr_length)
genome_len_no_micro <- sum(subset(linksumm, !chrom %in% c(21:28))$Chr_length)

#~~ And do lots of calculations

shufftab$IS_full <- 2*
  (shufftab$chrom_len_phase_0/shufftab$chrom_len)*
  (1-(shufftab$chrom_len_phase_0/shufftab$chrom_len))*
  (shufftab$chrom_len/genome_len)^2

shufftab$IS_macro <- 2*
  (shufftab$chrom_len_phase_0/shufftab$chrom_len)*
  (1-(shufftab$chrom_len_phase_0/shufftab$chrom_len))*
  (shufftab$chrom_len/genome_len_macro)^2

shufftab$IS_macro[which(!shufftab$chrom %in% c(1:12, 29))] <- 0

shufftab$IS_no_micro <- 2*
  (shufftab$chrom_len_phase_0/shufftab$chrom_len)*
  (1-(shufftab$chrom_len_phase_0/shufftab$chrom_len))*
  (shufftab$chrom_len/genome_len_no_micro)^2

shufftab$IS_no_micro[which(shufftab$chrom %in% c(21:28))] <- 0

write.table(shufftab, "results/3_Full_Shuffling_Summary_by_Chr.txt", row.names = F, sep = "\t", quote = F)

head(shufftab)

shuffsumm <- shufftab %>%
  group_by(meiosis) %>%
  summarise(intra_shuff = sum(IS_full),
            intra_shuff_macro = sum(IS_macro),
            intra_shuff_no_micro = sum(IS_no_micro))

recsumm <- left_join(recsumm, shuffsumm)

rm(shufftab, shuffsumm)

ggplot(recsumm, aes(intra_shuff)) + geom_histogram()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. How many genes?              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ format the genes

spargff <- subset(spargff, V3 == "gene")[,c(1, 4, 5)]
spargff <- spargff[grep("chr", spargff$V1),]
spargff$V1 <- gsub("chr", "", spargff$V1)
spargff$V1 <- gsub("1A", "29", spargff$V1)
spargff <- subset(spargff, V1 %in% 1:29)
spargff$V1 <- as.numeric(spargff$V1)

#~~ add information to the chunks

head(chunktab)

chunktab$genecount <- NA
chunktab$genecount_true <- NA

for(i in 1:nrow(chunktab)){
  
  chunktab$genecount[i] <- length(which(spargff$V1 == chunktab$chrom[i] &
                                          spargff$V5 < chunktab$stop[i] &
                                          spargff$V4 > chunktab$start[i]))
  chunktab$genecount_true[i] <- length(which(spargff$V1 == chunktab$chrom[i] &
                                               spargff$V5 < chunktab$stop_true[i] &
                                               spargff$V4 > chunktab$start_true[i]))
  
}


#~~ Get chr lengths and chunk proportions

x1 <- chunktab %>%
  group_by(KEY) %>%
  summarise(chrom_gene = sum(genecount),
            coverage_gene = sum(genecount_true))

x2 <- chunktab %>%
  group_by(KEY, Phase) %>%
  summarise(chrom_gene_phase = sum(genecount),
            coverage_gene_phase = sum(genecount_true))

x2 <- melt(x2, id.vars = c("KEY", "Phase"))
x2$variable <- paste0(x2$variable, "_", x2$Phase)
x2$Phase <- NULL

x2 <- dcast(x2,formula = KEY ~ variable)
x2 <- left_join(x2, x1)

head(x2)

#~~ Add to recchr

totalgenes <- spargff %>%
  group_by(V1) %>% 
  summarise(total_genes = n())
names(totalgenes)[1] <- "chrom"

shufftab <- left_join(recchr, x2)
shufftab <- left_join(shufftab, totalgenes)
head(shufftab)

na.lines <- which(is.na(shufftab$chrom_gene))
shufftab$chrom_gene[na.lines] <- shufftab$total_genes[na.lines]
shufftab$chrom_gene_phase_0 [na.lines] <- shufftab$chrom_gene[na.lines] 
shufftab$chrom_gene_phase_1 [na.lines] <- 0 
shufftab$coverage_gene_phase_0[na.lines] <- shufftab$coverage_gene [na.lines] 
shufftab$coverage_gene_phase_1 [na.lines] <- 0 

rm(x, x1, x2, i, na.lines)

head(shufftab)

shufftab <- subset(shufftab, select = -c(total_genes, chrom_gene_phase_1 , coverage_gene_phase_1, left_coverage, right_coverage))

#~~ Now add proportion of genome

genome_len <- sum(totalgenes$total_genes)
genome_len_macro <- sum(subset(totalgenes, chrom %in% c(1:12, 29))$total_genes)
genome_len_no_micro <- sum(subset(totalgenes, !chrom %in% c(21:28))$total_genes)

#~~ And do lots of calculations

shufftab$IS_gene_full <- 2*
  (shufftab$chrom_gene_phase_0/shufftab$chrom_gene)*
  (1-(shufftab$chrom_gene_phase_0/shufftab$chrom_gene))*
  (shufftab$chrom_gene/genome_len)^2

shufftab$IS_gene_macro <- 2*
  (shufftab$chrom_gene_phase_0/shufftab$chrom_gene)*
  (1-(shufftab$chrom_gene_phase_0/shufftab$chrom_gene))*
  (shufftab$chrom_gene/genome_len_macro)^2

shufftab$IS_gene_macro[which(!shufftab$chrom %in% c(1:12, 29))] <- 0

shufftab$IS_gene_no_micro <- 2*
  (shufftab$chrom_gene_phase_0/shufftab$chrom_gene)*
  (1-(shufftab$chrom_gene_phase_0/shufftab$chrom_gene))*
  (shufftab$chrom_gene/genome_len_no_micro)^2

shufftab$IS_gene_no_micro[which(shufftab$chrom %in% c(21:28))] <- 0

head(shufftab)

shuffsumm <- shufftab %>%
  group_by(meiosis) %>%
  summarise(intra_shuff_gene = sum(IS_gene_full),
            intra_shuff_gene_macro = sum(IS_gene_macro),
            intra_shuff_gene_no_micro = sum(IS_gene_no_micro))

recsumm <- left_join(recsumm, shuffsumm)

rm(shufftab, shuffsumm)

ggplot(recsumm, aes(intra_shuff_gene)) + geom_histogram()
ggplot(recsumm, aes(sex, intra_shuff_gene)) + geom_boxplot()

head(chunktab)

write.table(chunktab, "results/3_Phased_Chromosomal_Chunks.txt", row.names = F, sep = "\t", quote = F)
write.table(recsumm, "results/3_Full_Recombination_Phenotypes_QCed.txt", row.names = F, sep = "\t", quote = F)
write.table(recchr, "results/3_COs_per_chromosome_QCed.txt", row.names = F, sep = "\t", quote = F)
write.table(rectab, "results/3_Individual_COs_QCed.txt", row.names = F, sep = "\t", quote = F)


# recsumm_old <- read.table("0_superseded/results_20221206/4_Full_Recombination_Phenotypes_QCed.txt", header = T, stringsAsFactors = F)
# recsumm_old <- subset(recsumm_old, select = c(meiosis, intra_shuff, intra_shuff_gene))
# names(recsumm_old)[2:3] <-c("X2", "X3")
# 
# recsumm <- left_join(recsumm, recsumm_old)

# ggplot(recsumm, aes(intra_shuff, X2)) +
#   geom_point(alpha = 0.2)
# 
# ggplot(recsumm, aes(intra_shuff_gene, X3)) +
#   geom_point(alpha = 0.2)
