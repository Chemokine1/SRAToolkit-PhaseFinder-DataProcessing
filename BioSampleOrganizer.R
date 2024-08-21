#### Call libraries ####
library(stringr)
library(dplyr)
library(ggplot2)

#### generate PhaseFinder results paths ####
file_list1 <- list.files(path = "~/Documents/R/Technion_NGZ/Results_PhaFin_Eran_1_to_320th_EERs", full.names = TRUE)
file_list2 <- list.files(path = "~/Documents/R/Technion_NGZ/Results_PhaFin_Eran_320th_to_629th_ERRs", full.names = TRUE)
file_list <- c(file_list1,file_list2)

#### metadata row ####
data <- read.delim("~/Documents/R/Technion_NGZ/biosample_result_Antibiotics.txt")
data_lines <- readLines("~/Documents/R/Technion_NGZ/biosample_result_Antibiotics.txt") #bandate

##### constructing metadata & PhaseFinder results #####
# Empty lists to store the data
ids <- c()
descriptions <- c()
attributes_list <- list()

# Regular expressions to capture relevant data
id_pattern <- "BioSample: (SAMEA\\d+)"
desc_pattern <- "^Description:"

# Empty data frame to store the combined results
combined_data <- data.frame()

# counter tracking fikes
count <- 0

## phase finder combine all results ##
for (file in file_list) {
  file_name_with_extension <- basename(file)
  
  # processing only with .txt files
  if (grepl("\\.txt$", file_name_with_extension)) {
    # Remove the specific extension from the file name
    file_name <- sub("\\.ratio.txt$", "", file_name_with_extension)
    print(file_name)
    print(file)
    
    # Read the file if it is a .txt file
    data <- read.delim(file, header = TRUE)
    
    # Add a new column to store the file identifier
    data$err_accession <- file_name
    
    # Combine the data from this file with the main dataframe
    combined_data <- bind_rows(combined_data, data)
    
    # Increment the file count
    count = count + 1
    print(count)
  }
}

# Rename saving copy
combined_data1 <- combined_data

# count reads of IRs
combined_data1$counts <- combined_data1$Pe_F + combined_data1$Pe_R
combined_data1_filtered <- combined_data1 %>% filter(counts > 12)
combined_data1_filtered <- combined_data1_filtered %>% 
  group_by(ID) %>%
  mutate(Pe_ratio_mean = mean(Pe_ratio)) %>%
  filter(Pe_ratio_mean > 0.05)

hist(combined_data1_filtered$Pe_ratio)
boxplot(combined_data1_filtered$Pe_ratio_mean)

# Clean the ERRs names
combined_data1_filtered$err_accession <- gsub("out_", "", combined_data1_filtered$err_accession)
combined_data1_filtered$err_accession <- gsub("_\\d$", "", combined_data1_filtered$err_accession)

length(unique(combined_data1_filtered$err_accession)) #158 ERRs

{
  # Read the data file
  data_lines <- readLines("~/Documents/R/Technion_NGZ/biosample_result_Antibiotics.txt")
  
  # Initialize empty lists to store the data
  ids <- c()
  descriptions <- c()
  attributes_list <- list()
  
  # Regular expressions to capture relevant data
  id_pattern <- "BioSample: (SAMEA\\d+)"
  desc_pattern <- "^Description:"
  
  # Process each line
  for (i in seq_along(data_lines)) {
    if (grepl("Identifiers:", data_lines[i])) {
      # Capture the SAMEA identifier from the Identifiers line
      ids <- c(ids, str_match(data_lines[i], id_pattern)[,2])
    }
    if (grepl(desc_pattern, data_lines[i])) {
      # Capture the description from the line immediately following "Description:"
      descriptions <- c(descriptions, data_lines[i + 1])
    }
    if (grepl("Attributes:", data_lines[i])) {
      # Capture all attributes until another non-attribute line is reached
      j <- i + 1
      attributes <- list()
      while(j <= length(data_lines) && grepl("^\\s+/", data_lines[j])) {
        matches <- str_match(data_lines[j], "^\\s+/([^=]+)=\"([^\"]+)\"")
        attributes[[matches[2]]] <- matches[3]
        j <- j + 1
      }
      attributes_list[[length(attributes_list) + 1]] <- attributes
    }
  }
  
  # Create a data frame from the lists
  df <- data.frame(ID = ids, Description = descriptions, stringsAsFactors = FALSE)
  
  # Add attribute columns to the dataframe
  all_keys <- unique(unlist(lapply(attributes_list, names)))
  for (key in all_keys) {
    df[[key]] <- sapply(attributes_list, function(x) x[[key]])
  }
  
  # Replace NULL with NA in the dataframe
  df[is.null(df)] <- NA
  
  # Print the dataframe
  print(df)
  
} # create the metadata data frame #df_selected

df_selected <- df[,c("ID", "Description", "INSDC status", "Submitter Id", "scientific_name", "sex", "tissue", "#Mouse", "Feces_TimePoint", "region", "treatment", "age", "group", "participant", "day", "phase")]
View(df_selected)

# Loading the files of the SAMEA to join by the coresponding ERRs
data1 <- read.delim("ERR2SAME_Data.txt")
data1 <- data1 %>% select("run_accession","sample_accession")
data_phase_ERRsSAMEA <- left_join(combined_data1_filtered,data1, by=c("err_accession"="run_accession"))

# Clean the metadata based on the corresponding SAMEA of the Phase data
samea_needed <- unique(data_phase_ERRsSAMEA$sample_accession)
meta_data <- df_selected[df_selected$ID %in% samea_needed,]
View(meta_data)

#### Creates data to work on
aa <- inner_join(data_phase_ERRsSAMEA, meta_data, by = c("sample_accession" = "ID"))
View(aa)

aaa <- aa %>%
  group_by(ID) %>%
  mutate(IDcounts = n()) %>%
  filter(IDcounts > 30)

View(aaa)

#### EDA 
# make columns are factors
aaa$group <- as.factor(unlist(aaa$group))
aaa$age <- as.factor(unlist(aaa$age))
aaa$tissue <- as.factor(unlist(aaa$tissue))

# Create the boxplot
ggplot(aaa, aes(x = as.factor(ID), y = Pe_ratio, color = c(group))) + 
  geom_boxplot()

aaa <- aa %>% filter(ID == "NZ_LR699004.1:3222397-3222422-3222582-3222607")
#NZ_LR699004.1:3222397-3222422-3222582-3222607

boxplot(aaa$Pe_ratio)
library(ggplot2)
aaa$day <- na.omit(aaa$day)
aaa$day <- as.numeric(aaa$day)
aaa$tissue <- as.character(aaa$tissue)

ggplot(aaa, aes(x = day, y = Pe_ratio)) + geom_boxplot()

unique(aaa$day)
str(aaa)

