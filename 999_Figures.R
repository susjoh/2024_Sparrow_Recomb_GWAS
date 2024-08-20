library(dplyr)
library(reshape2)
library(ggplot2)

recsumm <- read.table("clean_data_for_John/3_Full_Recombination_Phenotypes_QCed.txt", header = T)

recsumm$sex2 <- ifelse(recsumm$sex == "M", "Male", "Female")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Distributions of recombination measures #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ggplot(recsumm, aes(yapp_CO_count_no_micro)) +
  geom_histogram(binwidth = 1, col = "grey") +
  theme_bw() +
  labs(x = "Autosomal Crossover Count")

#~~ Boxplot of ACC by sex

ggplot(recsumm, aes(sex2, yapp_CO_count_no_micro, fill = sex2)) +
  geom_boxplot(width = 0.6, outlier.size = 1, outlier.alpha = 0, alpha) +
  ggtitle("A.") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  labs(x = "Sex", y = "Autosomal Crossover Count (ACC)")

ggsave("figs/999_Sex_Diff_ACC.png", width = 3, height = 3)

ggplot(recsumm, aes(sex2, intra_shuff_no_micro, fill = sex2)) +
  geom_boxplot(width = 0.6, outlier.size = 1, outlier.alpha = 0.5) +
  ggtitle("B.") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  labs(x = "Sex", y = "Intra-chromosomal shuffling (r_intra)")
ggsave("figs/999_Sex_Diff_rintra.png", width = 3, height = 3)


library(ggridges)

ggplot(recsumm, aes(y = sex2, x = yapp_CO_count_no_micro, fill = sex2)) +
  geom_density_ridges(scale = 0.4, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.alpha = 0.2, position = position_nudge(y=-0.08), outlier.size = 1) +
  ggtitle("A.") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  labs(y = "Sex", x = "Autosomal Crossover Count (ACC)")
ggsave("figs/999_Sex_Diff_ACC_density.png", width = 4, height = 4)
ggsave("figs/999_Sex_Diff_ACC_density.pdf", width = 4, height = 4)


ggplot(recsumm, aes(y = sex2, x = intra_shuff_no_micro, fill = sex2)) +
  geom_density_ridges(scale = 0.4, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.alpha = 0.2, position = position_nudge(y=-0.08), outlier.size = 1) +
  ggtitle("B.") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  labs(y = "Sex", x = "Intra-chromosomal shuffling (r_intra)")
ggsave("figs/999_Sex_Diff_rintra_density.png", width = 4, height = 4)


ggplot(recsumm, aes(y = sex2, x = yapp_CO_count_no_micro, fill = sex2)) +
  geom_density_ridges(scale = 0.4, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.alpha = 0.2, position = position_nudge(y=-0.08), outlier.size = 1) +
  ggtitle("A. Crossover Count") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12)) +
  labs(y = "Sex", x = "Autosomal Crossover Count (ACC)")


ggplot(recsumm, aes(y = sex2, x = intra_shuff_no_micro, fill = sex2)) +
  geom_density_ridges(scale = 0.4, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.alpha = 0.2, position = position_nudge(y=-0.08), outlier.size = 1) +
  ggtitle("B. Intra-chromosomal shuffling") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12)) +
  labs(y = "Sex", x = "Intra-chromosomal shuffling (r_intra)")



# jittered_points = TRUE,
# position = position_points_jitter(width = 0.5, height = 0),
# point_shape = '|', point_size = 3

#~~ Histograms....

ggplot(recsumm, aes(yapp_CO_count_no_micro)) +
  geom_histogram(binwidth = 1, col = "grey") +
  ggtitle("A.") +
  theme_bw() +
  facet_wrap(~sex2, ncol = 1) +
  labs(x = "Autosomal Crossover Count (ACC)")
ggsave("figs/999_Sex_Diff_ACC_histogram.png", width = 4, height = 4)

ggplot(recsumm, aes(intra_shuff_no_micro)) +
  geom_histogram(binwidth = 0.0005, col = "grey") +
  ggtitle("B.") +
  theme_bw() +
  facet_wrap(~sex2, ncol = 1) +
  labs(x = "Intra-chromosomal shuffling (r_intra)")
ggsave("figs/999_Sex_Diff_rintra_histogram.png", width = 4, height = 4)

#~~ Summary stats

recsumm %>% group_by(sex) %>% summarise(mean_Acc = mean(yapp_CO_count_no_micro, na.rm=T),
                                        mean_rintra = mean(intra_shuff_no_micro, na.rm=T),
                                        var_Acc = var(yapp_CO_count_no_micro, na.rm=T),
                                        var_rintra = var(intra_shuff_no_micro, na.rm=T)) %>% data.frame


ggplot(recsumm, aes(sex2, yapp_CO_count_no_micro, fill = sex2)) +
  geom_violin() +
  ggtitle("A.") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  labs(x = "Sex", y = "Autosomal Crossover Count (ACC)")

ggsave("figs/999_Sex_Diff_ACC.png", width = 3, height = 3)

ggplot(recsumm, aes(sex2, intra_shuff_no_micro, fill = sex2)) +
  geom_boxplot(width = 0.6, outlier.size = 1, outlier.alpha = 0.5) +
  ggtitle("B.") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  labs(x = "Sex", y = "Intra-chromosomal shuffling (r_intra)")
ggsave("figs/999_Sex_Diff_rintra.png", width = 3, height = 3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Correlations of ACC and r_intro and r_gene #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


x <- subset(recsumm, select = c(yapp_CO_count_no_micro, intra_shuff_no_micro, intra_shuff_gene_no_micro)) %>% na.omit

ggplot(x, aes(yapp_CO_count_no_micro, intra_shuff_no_micro)) +
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Autosomal Crossover Count (ACC)", y = "Intra-chromosome shuffling (r_intra)")
ggsave("figs/999_ACC_vs_rintra_corr.png", width = 4, height = 4)

cor.test(x$yapp_CO_count_no_micro, x$intra_shuff_no_micro)


ggplot(x, aes(intra_shuff_gene_no_micro, intra_shuff_no_micro)) +
  geom_point(alpha = 0.3) +
  stat_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Intra-chromosome shuffling (r_intra)", y = "Intra-chromosome shuffling (r_gene)")
ggsave("figs/999_rgene_vs_rintra_corr.png", width = 4, height = 4)

cor.test(x$intra_shuff_gene_no_micro, x$intra_shuff_no_micro)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. GWAS Results                               #
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

gwasres$Model[which(gwasres$Model == "ACC")] <- "Autosomal Crossover Count"
gwasres$Model[which(gwasres$Model == "IntraR")] <- "Intra-chromosomal shuffling"
gwasres$Model[which(gwasres$Model == "IntraR_cofix")] <- "Intra-chromosomal shuffling (Adjusted)"


chrinfo <- gwasres %>%
  group_by(Chromosome2) %>%
  summarise(Mid =  min(Cumu) + ((max(Cumu) - min(Cumu))/2))
chrinfo$Chr_labs <- c("0", "1", "1A", "2", "3", "4", "5", 
                      "6", "7", "8", "9", "", "11", "", "", "14", "", "", 
                      "", "", "", "", "22", "", "", "", "", "", "", 
                      "Z")

x <- gwasres[sample(1:nrow(gwasres), 1000),]

ggplot(gwasres, aes(Cumu, -log10(Pc1df), col = factor(colour))) +
  geom_point(alpha = 0.35, size = 1) +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed") +
  facet_grid(Sex2 ~ Model) +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$Chr_labs) +
  labs(x = "Chromosomes", y = "-log10(P)")

ggsave("figs/999_GWAS_Results_full.png", width = 8, height = 4)


ggplot(subset(gwasres, Model %in% c("Autosomal Crossover Count", "Intra-chromosomal shuffling")), aes(Cumu, -log10(Pc1df), col = factor(colour))) +
  geom_point(alpha = 0.35, size = 1) +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed") +
  facet_grid(Sex2 ~ Model) +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$Chr_labs) +
  labs(x = "Chromosomes", y = "-log10(P)")

ggsave("figs/999_GWAS_Results.png", width = 8, height = 4)


ggplot(gwasres, aes(-log10(ExpP), -log10(Pc1df))) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed") +
  geom_point(alpha = 0.35) +
  facet_grid(Sex2 ~ Model) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Expected -log10(P)", y = "Observed -log10(P)")
ggsave("figs/999_GWAS_Results_PP_plot_full.png", width = 8, height = 4)

x <- gwasres[grep("Intra", gwasres$Model),]
x <- subset(x, select = c(Model, SNP.Name, Sex2, Pc1df))
x <- dcast(x, formula = SNP.Name ~ Model + Sex2)
head(x)

plot(log10(x[,2]) ~ log10(x[,4]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Empirical Bayes Results                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

bayesres <- read.table("results/6_GWAS_Results_EB.txt", header = T)


bayesres$Sex2 <- ifelse(bayesres$Sex == "M", "Male", "Female")
table(bayesres$Sex2)

bayesres$Chromosome2 <- bayesres$Chromosome
bayesres$Chromosome2[which(bayesres$Chromosome == 29)] <- "1A"
bayesres$Chromosome2[which(bayesres$Chromosome == 30)] <- 0
bayesres$Chromosome2[which(bayesres$Chromosome == 32)] <- "Z"
bayesres$Chromosome2 <- factor(bayesres$Chromosome2, levels = c(0:1, "1A", 2:15, 17:28, "Z"))

bonf <- 0.05/length(unique(bayesres$SNP.Name))

x <- subset(bayesres, Model == "ACC" & Sex2 == "Female")
x <- arrange(x, Chromosome2, Position)
x$Position[which(x$Chromosome2 == 0)] <- 1:length(x$Position[which(x$Chromosome2 == 0)])
x$Diff <- c(1, diff(x$Position))
x$Diff[which(x$Diff < 0)] <- 1
x$Cumu <- cumsum(x$Diff)
x <- subset(x, select = c(SNP.Name, Cumu))


bayesres$Cumu <- NULL
bayesres <- left_join(bayesres, x)


table(bayesres$Chromosome2)

bayesres$colour <- as.numeric(bayesres$Chromosome2) %% 2
bayesres$Model <- ifelse(bayesres$Model == "ACC", "Autosomal Crossover Count", "Intra-chromosomal shuffling")

chrinfo <- bayesres %>%
  group_by(Chromosome2) %>%
  summarise(Mid =  min(Cumu) + ((max(Cumu) - min(Cumu))/2))
chrinfo$Chr_labs <- c("0", "1", "1A", "2", "3", "4", "5", 
                      "6", "7", "8", "9", "", "11", "", "", "14", "", "", 
                      "", "", "", "", "22", "", "", "", "", "", "", 
                      "Z")

x <- bayesres[sample(1:nrow(bayesres), 1000),]

ggplot(bayesres, aes(Cumu, -log10(lfsr), col = factor(colour))) +
  geom_point(alpha = 0.35, size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  facet_grid(Sex2 ~ Model) +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$Chr_labs) +
  labs(x = "Chromosomes", y = "-log10(P)")

ggsave("figs/999_empirical_Bayes_lfsr_plot.png", width = 7, height = 4)

ggplot(bayesres, aes(NegativeProb)) +
  geom_histogram() +
  facet_wrap(Sex2 ~ Model)

bayesres %>% group_by(Sex2, Model) %>% summarise(NoSNPs_05 = length(which(lfsr < 0.05)),
                                                 NoSNPs_01 = length(which(lfsr < 0.01)),
                                                 NoSNPs_001 = length(which(lfsr < 0.001)),
                                                  PCsnps_05 = length(which(lfsr < 0.05))/length(unique(SNP.Name))*100,
                                                 PCsnps_01 = length(which(lfsr < 0.01))/length(unique(SNP.Name))*100,
                                                 PCsnps_0001 = length(which(lfsr < 0.001))/length(unique(SNP.Name))*100)

ggplot(bayesres, aes(-log10(lfsr), -log10(Pc1df))) +
  geom_point(alpha = 0.1) +
  geom_vline(xintercept = 2, linetype = "dashed") +
  facet_grid(Sex2 ~ Model) +
  theme_bw()

#~~ Output the non-zero SNPs

# x <- subset(bayesres, lfsr < 0.05)
# x <- subset(x, select = -c(Strand, Sex, Chromosome, Cumu, colour))
# x <- x[,c("Model", "Sex2", "Chromosome2", "Position", "A1", "A2", "Q.2", "effB", "se_effB", "chi2.1df", 
#           "Pc1df", "SNP.Name", "NegativeProb", "PositiveProb", "lfsr", "svalue", "lfdr", "qvalue", "PosteriorMean", "PosteriorSD")]
# names(x)[2:3] <- c("Sex", "Chromosome")
# x$Model <- ifelse(x$Model == "Autosomal Crossover Count", "ACC", "r_intra")
# 
# write.table(x, "results/6_Significant_SNPs_GWAS_EB.txt", row.names = F, sep = "\t", quote = F)


#~~ Heritability plot

h2tab <- read.table("results/5_Univar_Prop_of_Variance.txt", header = T, sep = "\t")

x <- data.frame(Effect = c("vm(parent, ainv)", "ide(parent)", "units!units", "vm(parent, grminv)"),
                Effect2 = c("Additive Genetic", "Identity", "Residual", "Additive Genetic"))

h2tab <- left_join(h2tab, x)
h2tab$Effect2 <- factor(h2tab$Effect2, levels = rev(c("Additive Genetic", "Identity", "Residual")))

h2tab <- subset(h2tab, GRM == "GRM" & RespVar != "intra_shuff_no_micro_cofix" & Effect2 != "Residual")

h2tab$RespVar <- ifelse(h2tab$RespVar == "yapp_CO_count_no_micro", "Crossover Count", "Intra-chr. Shuffling")

h2tab$Sex <- ifelse(h2tab$Sex == "M", "Male", "Female")

h2tab$Fill <- NA
h2tab$Fill[which(h2tab$Sex == "Male" & h2tab$Effect2 == "Additive Genetic")] <- "#1F78B4"
h2tab$Fill[which(h2tab$Sex == "Female" & h2tab$Effect2 == "Additive Genetic")] <- "#E31A1C"
h2tab$Fill[which(h2tab$Sex == "Male" & h2tab$Effect2 == "Identity")] <- "#A6CEE3"
h2tab$Fill[which(h2tab$Sex == "Female" & h2tab$Effect2 == "Identity")] <- "#FB9A99"

h2tab$Fill <- factor(h2tab$Fill, levels = c("#A6CEE3", "#1F78B4", "#FB9A99", "#E31A1C"))

ggplot(subset(h2tab, GRM == "GRM"), aes(Sex, Estimate, fill = Fill)) +
  geom_bar(stat = "identity", colour = "grey50") +
  theme_bw() +
  facet_wrap(~RespVar) +
  scale_fill_identity() +
  coord_cartesian(ylim = c(0, 0.4)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme(axis.text.x  = element_text (size = 11),
        axis.text.y  = element_text (size = 11),
        strip.text.x = element_text (size = 11),
        strip.text.y = element_text (size = 11),
        axis.title.y = element_text (size = 11),
        axis.title.x = element_text (size = 11))



ggplot(subset(h2tab, Effect2 == "Additive Genetic"), aes(Sex, Estimate, fill = Sex)) +
  geom_bar(stat = "identity", colour = "grey50") +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = 0.1) +
  theme_bw() +
  facet_wrap(~RespVar) +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  labs(y = "Heritability (h^2)") +
  theme(axis.text.x  = element_text (size = 11),
        axis.text.y  = element_text (size = 11),
        strip.text.x = element_text (size = 11),
        strip.text.y = element_text (size = 11),
        axis.title.y = element_text (size = 11),
        axis.title.x = element_text (size = 11),
        legend.position = "none")
