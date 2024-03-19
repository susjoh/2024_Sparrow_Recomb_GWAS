#
# LINKAGE MAPPING STATISTICS
# 

library(dplyr)
library(ggplot2)
library(reshape2)

linkmap <- read.table("results/passer_domesticus_linkage_map_v6.txt", header = T)

head(linkmap)
linkmap <- melt(linkmap, id.vars = c("chr", "snpid", "Mb"))
linkmap$variable <- ifelse(linkmap$variable == "male_cM", "Males", "Females")
linkmap$chr[which(linkmap$chr == 29)] <- "1A"
linkmap$chr <- factor(linkmap$chr, levels = c(1, "1A", 2:15, 17:28))


ggplot(linkmap, aes(Mb/1e6, value, col = variable)) +
  geom_point(alpha = 0.2, size = 1) +
  facet_wrap(~chr, scales = "free") +
  labs(x = "Chromosome Position (Mb)", y = "Linkage Map Position (cM)", col = "") +
  theme_bw() +
  scale_colour_brewer(palette = "Set1")
ggsave("figs/4_Linkage_Maps.png", width = 12, height = 8)

#~~ Make summary table

linkmap <- read.table("results/passer_domesticus_linkage_map_v6.txt", header = T)
linkmap$chr[which(linkmap$chr == 29)] <- "1A"
linkmap$chr <- factor(linkmap$chr, levels = c(1, "1A", 2:15, 17:28))

x <- linkmap %>% group_by(chr) %>% summarise(FCM = max(female_cM),
                                             MCM = max(male_cM),
                                             ChrLen = max(Mb),
                                             NSNPs = n()) %>% data.frame()
x$chr <- as.character(x$chr)
for(i in 2:ncol(x)) x[,i] <- as.numeric(x[,i])
x <- rbind(x, c("Total", colSums(x[,2:ncol(x)])))

write.table(x, "results/4_Linkage_Map_Chromosome_Summary.txt", row.names = F, sep = "\t", quote = F)

# unique positions?

linkmap <- linkmap[,c(1, 4, 5)] %>% unique()
nrow(linkmap)

#~~ Comparisons...

x1 <- melt(x, id.vars = c("chr", "ChrLen", "NSNPs"))
for(i in c(2, 3, 5)) x1[,i] <- as.numeric(x1[,i])
x1 <- subset(x1, chr != "Total")
x1$variable <- ifelse(x1$variable == "MCM", "Males", "Females")

ggplot(x1, aes(ChrLen/1e6, value, col = variable)) +
  geom_text(aes(label = chr), size = 3) +
  theme_bw() +
  labs(x = "Chromosome Length (Mb)", y = "Linkage Map Length (cM)", col = "Sex") +
  scale_colour_brewer(palette = "Set1") +
  stat_smooth(method = "lm") +
  theme(legend.position = "top")
ggsave("figs/4_Mb_vs_cM_chromosomes.png", width = 4, height = 4)


x <- subset(x, chr != "Total")
for(i in 2:ncol(x)) x[,i] <- as.numeric(x[,i])
x$ChrLen2 <- x$ChrLen/1e6

summary(glm(MCM ~ ChrLen2, data = x))
cor.test(x$MCM, x$ChrLen2)

summary(glm(FCM ~ ChrLen2, data = x))
cor.test(x$FCM, x$ChrLen2)
x1$RR <- x1$value/(x1$ChrLen/1e6)

ggplot(subset(x1, chr != 25), aes(ChrLen/1e6, value/(ChrLen/1e6), col = variable)) +
  geom_text(aes(label = chr), size = 3) +
  theme_bw() +
  labs(x = "Chromosome Length (Mb)", y = "Recombination rate (cM/Mb)", col = "Sex") +
  scale_colour_brewer(palette = "Set1") +
  stat_smooth(method = "lm", formula = y ~ log(x)) +
  theme(legend.position = "top")
ggsave("figs/4_Mb_vs_cMMb_chromosomes.png", width = 4, height = 4)

