library(dplyr)
library(reshape2)
library(ggplot2)
library(ggridges)
library(ggrepel)
library(gridExtra)
library(ggtext)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Figure 2                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

x <- read.table("results/4_Linkage_Map_Chromosome_Summary.txt", header = T)
x1 <- melt(x, id.vars = c("chr", "ChrLen", "NSNPs"))
for(i in c(2, 3, 5)) x1[,i] <- as.numeric(x1[,i])
x1 <- subset(x1, chr != "Total")
x1$variable <- ifelse(x1$variable == "MCM", "Males", "Females")

#~~ Make the plots

head(x1)
x1$PointLabelMacro = x1$chr
x1$PointLabelMacro[which(!x1$chr %in% c(1:9, "1A"))] <- ""
x1$PointLabelMicro = x1$chr
x1$PointLabelMicro[which(x1$chr %in% c(1:9, "1A"))] <- ""
x1$PointMacro<- "46"
x1$PointMacro[which(!x1$chr %in% c(1:9, "1A"))] <- "16"
x1$A <- "A."
x1$B <- "B."

## A

plotA <- ggplot(x1, aes(ChrLen/1e6, value, col = variable)) +
  scale_colour_brewer(palette = "Set1") +
  stat_smooth(method = "lm", linetype = 0) +
  geom_line(stat = "smooth", method = "lm", alpha = 0.7, size = 1) +
  geom_point() +
  geom_text_repel(aes(label = PointLabelMacro), box.padding = 0.2, show_guide = F) +
  labs(x = "Chromosome Length (Mb)", y = "Linkage Map Length (cM)", col = "Sex") +
  facet_wrap(~A)+
  theme_bw() +
  theme(legend.position = c(.5, 1.06),
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_text(hjust=0, size = 16))

## B

plotB <- ggplot(subset(x1, chr != 25), aes(ChrLen/1e6, value/(ChrLen/1e6), col = variable)) +
  scale_colour_brewer(palette = "Set1") +
  stat_smooth(method = "lm", linetype = 0, formula = y ~ log(x)) +
  geom_line(stat = "smooth", method = "lm", alpha = 0.7, size = 1, formula = y ~ log(x)) +
  geom_point() +
  geom_text_repel(aes(label = PointLabelMacro), box.padding = 0.2, show_guide = F) +
  labs(x = "Chromosome Length (Mb)", y = "Recombination rate (cM/Mb)", col = "Sex") +
  facet_wrap(~B)+
  theme_bw() +
  theme(legend.position = c(.5, 1.06),
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_text(hjust=0, size = 16))

# pdf("20221123_Heritability_Manuscript/MBE Submission v2/Figure_2.pdf", width = 8.5, height = 4.25)
# grid.arrange(plotA, plotB, ncol=2)
# dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Figure 3                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load("results/2_Full_QC_Sparrow_Yapp_Data.RData")

recchr_v2$chrom[which(recchr_v2$chrom == 29)] <- "1A"
recchr_v2$chrom <- factor(recchr_v2$chrom, levels = c(1, "1A", 2:15, 17:25, 27, 28))
recchr_v2$n <- factor(recchr_v2$n, levels = 8:0)

ggplot(recchr_v2, aes(chrom, fill = n)) +
  geom_bar() +
  geom_hline(yintercept = nrow(recsumm)/2, colour = "white", linetype = "dashed") +
  theme_bw() +
  scale_fill_brewer(palette = "GnBu") +
  labs(x = "Chromosome", y = "Proportion of Gametes", fill = "CO Count") +
  scale_y_continuous(breaks = seq(0, nrow(recsumm), nrow(recsumm)/5), labels = seq(0, 1, 0.2))

ggsave("20221123_Heritability_Manuscript/MBE Submission v2/Figure 3.pdf", width = 6, height = 4, dpi = 300)

#~~ Sex specific for the Supp Mat after reviewers comments.

x1 <- subset(recchr_v2, select = c(meiosis, sex)) %>% na.omit %>% unique
names(x1) <- c("meiosis", "sex2")
recchr_v2 <- left_join(recchr_v2, x1)
recchr_v2$sex2 <- ifelse(recchr_v2$sex2 == "F", "A. Females", "B. Males")

ggplot(recchr_v2, aes(chrom, fill = n)) +
  geom_bar() +
  theme_bw() +
  facet_wrap(~sex2, ncol = 1) +
  scale_fill_brewer(palette = "GnBu") +
  labs(x = "Chromosome", y = "Proportion of Gametes", fill = "CO Count") +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust=0, size = 12))
  scale_y_continuous(breaks = seq(0, nrow(recsumm), nrow(recsumm)/5), labels = seq(0, 1, 0.2))
ggsave("figs/4_Distribution_of_CO_counts_per_Chromosome_by_sex.png", width = 6, height = 6)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Figure 4                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

recsumm <- read.table("clean_data_for_John/3_Full_Recombination_Phenotypes_QCed.txt", header = T)
recsumm$sex2 <- ifelse(recsumm$sex == "M", "Male", "Female")

x1 <- subset(recsumm, select = c(sex2, yapp_CO_count_no_micro, intra_shuff_no_micro))
x1 <- melt(x1)
x1$pheno <- ifelse(x1$variable == "yapp_CO_count_no_micro", "A.", "B.")

ggplot(x1, aes(y = sex2, x = value, fill = sex2)) +
  geom_density_ridges(scale = 0.4, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.alpha = 0.2, position = position_nudge(y=-0.08), outlier.size = 1) +
  facet_wrap(~pheno, scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14, hjust = 0)) +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  labs(y = "Sex", x = expression("Autosomal Crossover Count (ACC)                          Intra-chromosomal shuffling (r"["intra"]~")"))


ggsave("20221123_Heritability_Manuscript/MBE Submission v2/Figure 4.pdf", width = 8, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Figure 5: GWAS Results                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

gwasres <- read.table("results/6_GWAS_Results.txt", header = T)
gwasres$Chromosome2 <- gwasres$Chromosome
gwasres$Chromosome2[which(gwasres$Chromosome == 29)] <- "1A"
gwasres$Chromosome2[which(gwasres$Chromosome == 30)] <- 0
gwasres$Chromosome2[which(gwasres$Chromosome == 32)] <- "Z"
gwasres$Chromosome2 <- factor(gwasres$Chromosome2, levels = c(0:1, "1A", 2:15, 17:28, "Z"))


x <- subset(gwasres, Model == "ACC" & Sex == "F")
x <- arrange(x, Chromosome2, Position)
x$Position[which(x$Chromosome2 == 0)] <- 1:length(x$Position[which(x$Chromosome2 == 0)])
x$Diff <- c(1, diff(x$Position))
x$Diff[which(x$Diff < 0)] <- 1
x$Cumu <- cumsum(x$Diff)
x <- subset(x, select = c(SNP.Name, Cumu))

bonf <- 0.05/length(unique(gwasres$SNP.Name))

gwasres$Cumu <- NULL
gwasres <- left_join(gwasres, x)

table(gwasres$Chromosome2)

gwasres$colour <- as.numeric(gwasres$Chromosome2) %% 2
gwasres$Sex2 <- ifelse(gwasres$Sex == "F", "Females", "Males")
#gwasres$Model <- ifelse(gwasres$Model == "ACC", "Autosomal Crossover Count", "Intra-chromosomal shuffling")

gwasres$Model[which(gwasres$Model == "ACC")] <- "Autosomal Crossover Count (ACC)"
gwasres$Model[which(gwasres$Model == "IntraR")] <- "Intra-chromosomal shuffling (r<sub>intra</sub>)"
gwasres$Model[which(gwasres$Model == "IntraR_cofix")] <- "Intra-chromosomal shuffling (Adjusted)"


chrinfo <- gwasres %>%
  group_by(Chromosome2) %>%
  summarise(Mid =  min(Cumu) + ((max(Cumu) - min(Cumu))/2))
chrinfo$Chr_labs <- c("0", "1", "1A", "2", "3", "4", "5", 
                      "6", "7", "8", "9", "", "11", "", "", "14", "", "", 
                      "", "", "", "", "22", "", "", "", "", "", "", 
                      "Z")

gwasres <- subset(gwasres, !Model %in% c("Intra-chromosomal shuffling (Adjusted)"))

x <- gwasres[sample(1:nrow(gwasres), 1000),]


ggplot(gwasres, aes(Cumu, -log10(Pc1df), col = factor(colour))) +
  geom_point(alpha = 0.35, size = 0.5) +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed") +
  facet_grid(Sex2 ~ Model, labeller = labeller(type = label_parsed)) +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "none", strip.text = ggtext::element_markdown()) +
  scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$Chr_labs) +
  labs(x = "Chromosomes", y = "-log10(P)")

ggsave("20221123_Heritability_Manuscript/MBE Submission v2/Figure 5.jpg", width = 7, height = 4, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Figure 6                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

x <- read.table("results/7_Chr_h2_FULL_results.txt", header = T, sep = "\t")

x$RespVar <- ifelse(x$RespVar == "ACC", "ACC", "r<sub>intra</sub>")
x$RespModel <- paste0(x$RespVar, ": ", x$Model)
x$RespSex <- paste0(x$RespVar, ": ", x$Sex)

fullslopes <- read.table("results/7_Chr_h2_regression_results.txt", header = T, sep = "\t")
fullslopes <- subset(fullslopes, Analysis == "all")

fullslopes$RespVar <- ifelse(fullslopes$RespVar == "ACC", "ACC", "r<sub>intra</sub>")
fullslopes$RespSex <- paste0(fullslopes$RespVar, ": ", fullslopes$Sex)
fullslopes$Slope <- paste0("slope = ", round(fullslopes$Estimate, 5))
fullslopes$Intercept <- paste0("intercept = ", round(fullslopes$Estimate, 5))
fullslopes$P  <- round(fullslopes$P, 3)

fullslopes$Label  <- paste0("Adj R<sup>2</sup> = ", round(fullslopes$AdjR2, 3), "<br>P ", ifelse(fullslopes$P < 0.001, "< 0.001", paste0("= ", fullslopes$P)))

ggplot() +
  geom_hline(yintercept = 0) +
  stat_smooth(data = x, aes(Length,  Estimate), method = "lm") +
  geom_errorbar(data = x, aes(Length,  Estimate, ymin = Estimate - SE, ymax = Estimate + SE), alpha = 0.5) +
  geom_point(data = x, aes(Length,  Estimate), size = 2, alpha = 0.5) +
  geom_richtext(data = subset(fullslopes, Effect == "Length"),
                aes(x = 0, y = 0.09, label = Label),
                size = 3.5, hjust = 0, colour = "red3", label.size = 0, fill = NA, label.colour = NA) +
  facet_grid(Model~RespSex, scales = "fixed") +
  coord_cartesian(ylim = c(-0.01, 0.105)) +
  theme_bw() +
  theme(legend.position = "none", strip.text = ggtext::element_markdown()) +
  labs(x = "Chromosome Length (Mb)", y = "Proportion of Phenotypic Variance")

#ggsave("20221123_Heritability_Manuscript/MBE Submission v2/Figure 6.pdf", width = 8, height = 5, dpi = 300)

head(fullslopes)

fullslopes <- subset(fullslopes, RespVar == "ACC")
x <- subset(x, RespVar == "ACC")

ggplot() +
  geom_hline(yintercept = 0) +
  stat_smooth(data = x, aes(Length,  Estimate), method = "lm") +
  geom_errorbar(data = x, aes(Length,  Estimate, ymin = Estimate - SE, ymax = Estimate + SE), alpha = 0.5) +
  geom_point(data = x, aes(Length,  Estimate), size = 2, alpha = 0.5) +
  geom_richtext(data = subset(fullslopes, Effect == "Length"),
                aes(x = 0, y = 0.09, label = Label),
                size = 3.5, hjust = 0, colour = "red3", label.size = 0, fill = NA, label.colour = NA) +
  facet_grid(Model~Sex, scales = "fixed") +
  coord_cartesian(ylim = c(-0.01, 0.105)) +
  theme_bw() +
  theme(legend.position = "none", strip.text = ggtext::element_markdown()) +
  labs(x = "Chromosome Length (Mb)", y = "Proportion of Phenotypic Variance")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# sex differences in landscape...#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

linkmap <- read.table("results/passer_domesticus_linkage_map_v6.txt", header = T)
linkmap$chr[which(linkmap$chr == 29)] <- "1A"
linkmap$chr <- factor(linkmap$chr, levels = c(1, "1A", 2:15, 17:28))
linkmap$bin <- ceiling(linkmap$Mb/1e6)

x <- linkmap %>%
  group_by(chr, bin) %>%
  summarise(Females = (max(female_cM) - min(female_cM))/(max(Mb) - min(Mb))*1e6,
            Males = (max(male_cM) - min(male_cM))/(max(Mb) - min(Mb))*1e6,
            map_dist = max(Mb) - min(Mb),
            count = n())

x$Cumu <- 1:nrow(x)

930

x <- subset(x, count > 5 & map_dist > 500000)

ggplot(x, aes(map_dist, count)) +
  geom_point()


x <- melt(x, id.vars = c("Cumu", "chr", "bin", "map_dist", "count"))
x <- subset(x, chr != 25)
x$chr2 <- factor(x$chr, levels = c(1, "1A", 2:15, 17:28))
x$sex <- ifelse(x$variable == "Females", "A. Females", "B. Males")

chrinfo <- x %>%
  group_by(chr) %>%
  summarise(Mid =  min(Cumu) + ((max(Cumu) - min(Cumu))/2))
chrinfo$Chr_labs <- c("1", "1A", "2", "3", "4", "5", 
                      "6", "7", "8", "9", "", "11", "", "", "14", "", "", 
                      "", "", "", "", "22", "", "", "", "", "")


ggplot(x, aes(Cumu, value, col = factor(as.numeric(chr2) %% 2), group = chr2)) +
  geom_line() +
  facet_wrap(~sex, ncol = 1) +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12)) +
  scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$Chr_labs) +
  labs(x = "Chromosomes", y = "Recombination Rate (cM/Mb)")

ggsave("figs/999_Recombination_landscape_by_sex.png", height = 5, width = 7)
  

