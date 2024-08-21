#### Set working directory ####
setwd("~/Documents/R/Technion_NGZ/Erans_data")
### call data
source("Data_Preprocessing.R")

#### IR1-4 I decidied to combine them ####
View(data_1ir)

boxplot(Pe_ratio ~ group ,data = data_1ir)
ggplot(data = data_1ir, aes(x = Pe_ratio, y = tissue2)) + geom_boxplot()

par(mfrow = c(2, 2))

hist(data_1ir$Pe_ratio, main = "47", xlab = "Pe_ratio", col = "lightblue", prob = TRUE)
lines(density(data_1ir$Pe_ratio), col = "red", lwd = 2)

hist(data_2ir$Pe_ratio, main = "49", xlab = "Pe_ratio", col = "lightgreen", prob = TRUE,ylim = c(0, 3))
lines(density(data_2ir$Pe_ratio), col = "red", lwd = 2)

hist(data_3ir$Pe_ratio, main = "44", xlab = "Pe_ratio", col = "lightblue", prob = TRUE)
lines(density(data_3ir$Pe_ratio), col = "red", lwd = 2)

hist(data_4ir$Pe_ratio, main = "42", xlab = "Pe_ratio", col = "lightgreen", prob = TRUE)
lines(density(data_4ir$Pe_ratio), col = "red", lwd = 2)

# Reset the plotting area to default (optional)
par(mfrow = c(1, 1))

#### Combining 47 + 49 data ####
hist(data_riboz47_49$Pe_ratio, main = "47,49", xlab = "Pe_ratio",
     col = "lightpink", prob = TRUE, breaks = 20)
lines(density(data_riboz$Pe_ratio), col = "red", lwd = 2)

title <- data_riboz %>% group_by(tissue2) %>% mutate(count_tissue = n()) %>% reframe(count_tissue) %>% unique()
boxplot(data = data_riboz47_49, Pe_ratio ~ tissue2, main = title)

p <- ggplot(data_riboz47_49, aes(x = tissue2, y = Pe_ratio)) +
  geom_boxplot() +
  ggtitle("Pe_ratio Across Different Tissues") + 
  geom_signif(comparisons = list(c("L", "M"), c("L", "S"), c("M", "S")),
              map_signif_level = TRUE)


l <- data_riboz47_49 %>% filter(tissue2 == "L") %>% pull(Pe_ratio)
m <- data_riboz47_49 %>% filter(tissue2 == "M") %>% pull(Pe_ratio)
s <- data_riboz47_49 %>% filter(tissue2 == "S") %>% pull(Pe_ratio)

wilcox.test(m,s)

#### All IR boxplot ####
bact.ole.genes

par(mfrow = c(3,1))
for (ir in bact.ole.genes) {
  data.temp <- id_count_df_alldata_bact %>% filter(ID == ir)
  hist(data.temp$Pe_ratio, main = ir, xlab = "Pe_ratio", col = "lightblue", prob = TRUE, breaks = 7)
  lines(density(data.temp$Pe_ratio), col = "red", lwd = 2)
}

# filter(group == "probiotics without antibiotics group - naive",
        #tissue2 == "S", phase == "baseline")
boxplot(data = data_riboz47_49_baseline, Pe_ratio ~ day, main = "47,49")
boxplot(data = data_riboz47_baseline, Pe_ratio ~ day, main = "47")
boxplot(data = data_riboz49_baseline, Pe_ratio ~ day, main = "49")

# filter(group == "probiotics without antibiotics group - naive",
#tissue2 == "S", phase == "intervention")

par(mfrow = c(1,1))
boxplot(data = data_riboz47_49_intervention, Pe_ratio ~ day, main = "47,49")
boxplot(data = data_riboz47_intervention, Pe_ratio ~ day, main = "47")
boxplot(data = data_riboz49_intervention, Pe_ratio ~ day, main = "49")

par(mfrow = c(3,1))
ggplot(data_riboz47_49_S, aes(x = interaction(tissue2, group, phase, day), y = Pe_ratio)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Tissue, Group, Phase", y = "Pe_ratio")

ggplot(data_riboz47_49_L, aes(x = interaction(tissue2, group, phase, day), y = Pe_ratio)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Tissue, Group, Phase", y = "Pe_ratio")

ggplot(data_riboz47_49_M, aes(x = interaction(tissue2, group, phase, day), y = Pe_ratio)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Tissue, Group, Phase", y = "Pe_ratio")





