#### Set working directory ####
setwd("~/Documents/R/Technion_NGZ/Erans_data")

#### call libraries ####
library(dplyr)
library(readxl)
library(writexl)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(ggsignif)
library(tidyr)


#### Filter participants 901 and 911 ####
data <- read_xlsx("Data_filtered_probiotics_subsample.xlsx")
data <- data %>%
  filter(!(participant == 901 | participant == 911))

#write_xlsx(data, "data_probiotics_no901and911.xlsx")


#### Call data ####
data <- read_xlsx("data_probiotics_no901and911.xlsx")


#### Data prepossessing ####
id_count_df <- data %>%
  group_by(ID) %>%
  mutate(count_id = n()) %>%
  reframe(count_id) %>%
  arrange(desc(count_id)) %>%
  unique()

# Pe_ratio
id_pv_mean <- data %>%
  group_by(ID) %>%
  mutate(mean_pv = mean(Pe_ratio)) %>%
  reframe(mean_pv) %>%
  arrange(desc(mean_pv)) %>%
  unique()

# match ids count with ids pv
id_count_mean_pv <- merge(id_count_df, id_pv_mean, by = 'ID')
id_count_mean_pv_filtered <- id_count_mean_pv %>% filter(count_id >= 4, )

id_count_df_alldata <- data %>%
  group_by(ID) %>%
  mutate(count_id = n(), mean_pv = mean(Pe_ratio)) %>%
  filter(count_id >= 4, mean_pv > 0) %>%
  arrange(desc(count_id))

#### Filling missing data ####
##### complete intervention #####
nona_id_count_df_alldata <- id_count_df_alldata %>%
  mutate(phase = ifelse(grepl("during intervention", group), "intervention", phase))

##### submitter id take care #####

participant_submiter <- c()
tiss_submiter <- c()
day_submitter <- c()

# Loop through each value in the vector
for (i in id_count_df_alldata$Submitter.Id) {
  parts <- unlist(strsplit(i, "\\."))
  
  # Assign values based on the number of parts
  participant_submiter <- c(participant_submiter, parts[1])
  if (length(parts) == 3) {
    tiss_submiter <- c(tiss_submiter, parts[2])
    day_submitter <- c(day_submitter, parts[3])
  } else if (length(parts) == 2) {
    tiss_submiter <- c(tiss_submiter, NA)
    day_submitter <- c(day_submitter, parts[2])
  }
}
id_count_df_alldata$participant_submiter <- participant_submiter
#id_count_df_alldata$tissue <- tiss_submiter
#id_count_df_alldata$day_submitter <- day_submitter

id_count_df_alldata$tiss_num <- ifelse(grepl("Stool", id_count_df_alldata$tissue), 0,
                                            ifelse(grepl("lumen", id_count_df_alldata$tissue), 1,
                                                   ifelse(grepl("mucosa", id_count_df_alldata$tissue), 2, NA)))


##### Assigning the names of day, tissue #####
#id_count_df_alldata$day <- day_submitter
id_count_df_alldata$tissue2 <- ifelse(grepl(1, id_count_df_alldata$tiss_num), "L",
                                        ifelse(grepl(2, id_count_df_alldata$tiss_num), "M", "S"))

#### Bacteria ID df ####
genes.bact <- read.csv('/Users/alonkedem/Documents/R/phase_ibd/PS_IDs_multiple.csv')
genes.bact$bact_ecc <- sub(":.*", "", genes.bact$ID)
genes.bact <- genes.bact[, c("bact_ecc", "Bacteria")]
genes.bact <- genes.bact %>% unique()

# Add the eccession name as a cloumn to id_count_df_alldata
id_count_df_alldata$bact_ecc <- sub(":.*", "", id_count_df_alldata$ID)

id_count_df_alldata_bact <- left_join(id_count_df_alldata, genes.bact, by = "bact_ecc")
id_count_df_alldata <- id_count_df_alldata_bact #Bandate

#### Participants / resistant & permissive ####
participants <- c(412, 413, 415, 418, 512, 513, 515, 518, 411, 414, 416, 417, 419, 420, 511, 514, 516, 517, 519, 520)
status <- c('resistant', 'resistant', 'resistant', 'resistant', 'resistant', 'resistant', 'resistant', 'resistant', 
            'permissive', 'permissive', 'permissive', 'permissive', 'permissive', 'permissive', 'permissive', 'permissive', 
            'permissive', 'permissive', 'permissive', 'permissive')
df <- data.frame(participants, status)

id_count_df_alldata <- left_join(id_count_df_alldata, df, by = c("participant" = "participants"))

##### retrive relavent columns #####
# Create dummie variable. 
id_count_df_alldata$ID_num <- as.numeric(factor(id_count_df_alldata$ID, levels = unique(id_count_df_alldata$ID)))
id_count_df_alldata <- id_count_df_alldata[, c("ID", "Pe_ratio", "err_accession", "reads_count",
                                                "Pe_ratio_mean", "Description", "sex", 
                                                "tissue", "age", "group", "participant","status", "day", 
                                                "phase", "count_id", "tissue2", 
                                                "ID_num", "bact_ecc","Bacteria")]

###### new df name = less mistakes ######
id_count_df_alldata1 <- id_count_df_alldata
id_count_df_alldata <- id_count_df_alldata1

#### trying to address NA values (group_by group, phase, day) ####
GroupedBy_group_phase_day <- id_count_df_alldata %>% group_by(group,phase,day) %>% mutate(n=n()) %>% reframe(n) %>% unique()
id_count_df_alldata2 <- id_count_df_alldata1
id_count_df_alldata3 <- id_count_df_alldata1

id_count_df_alldata2$phase <- ifelse(test = (id_count_df_alldata2$phase == "NA" & id_count_df_alldata2$group == "placebo group - during intervention" & id_count_df_alldata2$day == 4), "intervention", id_count_df_alldata2$phase)
id_count_df_alldata2$phase <- ifelse(test = (id_count_df_alldata2$phase == "NA" & id_count_df_alldata2$group == "placebo group - naive" & id_count_df_alldata2$day == 4), "baseline", id_count_df_alldata2$phase)
id_count_df_alldata2$phase <- ifelse(test = (id_count_df_alldata2$phase == "NA" & id_count_df_alldata2$group == "probiotics without antibiotics group - during intervention" & id_count_df_alldata2$day == 4), "intervention", id_count_df_alldata2$phase)
id_count_df_alldata2$phase <- ifelse(test = (id_count_df_alldata2$phase == "NA" & id_count_df_alldata2$group == "probiotics without antibiotics group - naive" & id_count_df_alldata2$day == 4), "baseline", id_count_df_alldata2$phase)

id_count_df_alldata2$phase <- ifelse(test = (id_count_df_alldata2$phase == "NA" & id_count_df_alldata2$group == "placebo group - during intervention" & id_count_df_alldata2$day == 5), "intervention", id_count_df_alldata2$phase)
id_count_df_alldata2$phase <- ifelse(test = (id_count_df_alldata2$phase == "NA" & id_count_df_alldata2$group == "placebo group - naive" & id_count_df_alldata2$day == 5), "baseline", id_count_df_alldata2$phase)
id_count_df_alldata2$phase <- ifelse(test = (id_count_df_alldata2$phase == "NA" & id_count_df_alldata2$group == "probiotics without antibiotics group - during intervention" & id_count_df_alldata2$day == 5), "intervention", id_count_df_alldata2$phase)
id_count_df_alldata2$phase <- ifelse(test = (id_count_df_alldata2$phase == "NA" & id_count_df_alldata2$group == "probiotics without antibiotics group - naive" & id_count_df_alldata2$day == 5), "baseline", id_count_df_alldata2$phase)

id_count_df_alldata2$phase <- ifelse(test = (id_count_df_alldata2$phase == "NA" & id_count_df_alldata2$group == "placebo group - during intervention" & id_count_df_alldata2$day == 6), "intervention", id_count_df_alldata2$phase)
id_count_df_alldata2$phase <- ifelse(test = (id_count_df_alldata2$phase == "NA" & id_count_df_alldata2$group == "placebo group - naive" & id_count_df_alldata2$day == 6), "baseline", id_count_df_alldata2$phase)
id_count_df_alldata2$phase <- ifelse(test = (id_count_df_alldata2$phase == "NA" & id_count_df_alldata2$group == "probiotics without antibiotics group - during intervention" & id_count_df_alldata2$day == 6), "intervention", id_count_df_alldata2$phase)
id_count_df_alldata2$phase <- ifelse(test = (id_count_df_alldata2$phase == "NA" & id_count_df_alldata2$group == "probiotics without antibiotics group - naive" & id_count_df_alldata2$day == 6), "baseline", id_count_df_alldata2$phase)

id_count_df_alldata2$phase <- ifelse(test = (id_count_df_alldata2$phase == "NA" & id_count_df_alldata2$group == "probiotics without antibiotics group - during intervention" & id_count_df_alldata2$day == 7), "intervention", id_count_df_alldata2$phase)

id_count_df_alldata2$phase <- ifelse(test = id_count_df_alldata2$day >= 8, "intervention", id_count_df_alldata2$phase)

## ChatGPT shortest way ##
#id_count_df_alldata3$phase
#id_count_df_alldata3 <- id_count_df_alldata1 %>%
#  mutate(phase = case_when(
#    phase == "NA" & group %in% c("placebo group - during intervention", "probiotics without antibiotics group - during intervention") & day %in% 4:7 ~ "intervention",
#    phase == "NA" & group %in% c("placebo group - naive", "probiotics without antibiotics group - naive") & day %in% 4:6 ~ "baseline",
#    day >= 8 ~ "intervention",
#    TRUE ~ phase  # Keep the original phase if no condition is met
#  ))
id_count_df_alldata <- id_count_df_alldata2

#View(id_count_df_alldata)
## 
data_1ir <- id_count_df_alldata %>% filter(ID_num == 1)
data_2ir <- id_count_df_alldata %>% filter(ID_num == 2)
data_3ir <- id_count_df_alldata %>% filter(ID_num == 3)
data_4ir <- id_count_df_alldata %>% filter(ID_num == 4)

##### I decide to connect these 4 areas because if I see a tendency ####
#that affects all areas then there is more power for this explanation 

data_riboz <- id_count_df_alldata %>%
  filter(ID_num %in% c(1, 2, 3, 4))

# For each ERR somethimes three areas were found in the same timepoint, i will
# not seperate or subsuet the data because this region was found 3 times.
# Why this region was found to be invertable more then the others? Those 
# regions were found to be same sequence but invertable differently? 
# Is this region really same sequence? 

data_riboz %>% group_by(ID) %>%
  mutate(id_count = n()) %>%
  reframe(Pe_ratio_mean, id_count) %>% unique()

#MAFFT is a fast and highly accurate tool for multiple sequence alignments.
# I used vectorbuilder to allign between the sequences i downloaded

data_riboz47_49 <- id_count_df_alldata %>%
  filter(ID_num %in% c(1, 2))

#### Bacteroides oleiciplenus YIT 12058 All regions ####
bact.ole.genes <- id_count_df_alldata %>%
  filter(Bacteria == "Bacteroides oleiciplenus YIT 12058") %>%
  pull(ID) %>% unique()

data_riboz5 <- id_count_df_alldata_bact %>% filter(ID == "NZ_JH992942.1:1067-1078-1259-1270")
data_riboz6 <- id_count_df_alldata_bact %>% filter(ID == "NZ_JH992949.1:2195-2206-2387-2398")
mean(data_riboz5$Pe_ratio)
mean(data_riboz6$Pe_ratio)

#### filter data_riboz47_49 Stool to exlore features ####
data_riboz47_49_baseline <- data_riboz47_49 %>% filter(group == "probiotics without antibiotics group - naive",
                                                tissue2 == "S",
                                                phase == "baseline")
data_riboz47_baseline <- data_riboz47_49_baseline %>% filter(ID == "NZ_JH992947.1:185030-185043-185707-185720")
data_riboz49_baseline <- data_riboz47_49_baseline %>% filter(ID == "NZ_JH992949.1:3027-3041-3705-3719")

data_riboz47_49_intervention <- data_riboz47_49 %>% filter(group == "probiotics without antibiotics group - naive",
                                                tissue2 == "S",
                                                phase == "intervention")
data_riboz47_49_intervention$day <- as.numeric(data_riboz47_49_intervention$day)

data_riboz47_intervention <- data_riboz47_49_intervention %>% filter(ID == "NZ_JH992947.1:185030-185043-185707-185720")
data_riboz49_intervention <- data_riboz47_49_intervention %>% filter(ID == "NZ_JH992949.1:3027-3041-3705-3719")

data_riboz47_49_S <- data_riboz47_49 %>% filter(tissue2 == "S")
data_riboz47_49_L <- data_riboz47_49 %>% filter(tissue2 == "L")
data_riboz47_49_M <- data_riboz47_49 %>% filter(tissue2 == "M")

riboz4749_S_in_pro <- data_riboz47_49_S %>% filter(phase == "intervention", )



#### Linear Model Data processing ####
data_riboz47_492 <- data_riboz47_49 %>%
  mutate_at(vars(status, group, day, tissue, tissue2, ID,sex, participant), as.factor)

library(lmerTest)
lm1 <- lmer(Pe_ratio ~ status + group + tissue2 + sex + (1|ID), data = data_riboz47_492)
summary(lm1)

# Load the ggplot2 package
library(ggplot2)

# Create a histogram of Pe_ratio by status
ggplot(data_riboz47_492, aes(x = Pe_ratio, fill = as.factor(status))) +
  geom_histogram(binwidth = 0.1, position = "dodge", alpha = 0.7) +
  labs(x = "Pe_ratio", y = "Count", fill = "Status") +
  theme_minimal()

ggplot(data_riboz47_492, aes(x = Pe_ratio)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(x = "Pe_ratio", y = "Count") +
  facet_wrap(~ status) +
  theme_minimal()

ggplot(data_riboz47_492, aes(y = Pe_ratio, x = as.factor(sex))) + geom_boxplot()
