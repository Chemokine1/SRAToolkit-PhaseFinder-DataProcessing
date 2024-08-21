#### Set working directory ####
setwd("~/Documents/R/Technion_NGZ/Erans_data")

source("Data_Preprocessing.R")

#### EDA ####
hist(id_count_df$count_id)
hist(id_pv_mean$mean_pv)

for (id in unique(id_count_df_alldata$ID)) {
  temp.data <- id_count_df_alldata %>% filter(ID == id)
  title = paste0("name: ",id,". ID count: ",temp.data$count_id)
  boxplot(temp.data$Pe_ratio, main = title)
}
View(id_count_df_alldata)


# plot for all of the IRs
plot_list <- list()
for (id in unique(id_count_df_alldata_bact$ID)) {
  temp.data <- id_count_df_alldata_bact %>% filter(ID == id)
  title <- paste0(temp.data$Bacteria[1], "\nUnformal IR name: ", temp.data$ID_num[1], ". Count: ", temp.data$count_id)
  
  p <- ggplot(temp.data, aes(y = Pe_ratio)) +
    geom_boxplot() +
    ggtitle(title) +
    theme_minimal()
  
  plot_list[[length(plot_list) + 1]] <- p
  
  if (length(plot_list) == 8) {
    grid.arrange(grobs = plot_list, ncol = 2, nrow = 4)
    plot_list <- list()
  }
}

#### Interesting genes ####
# IDs(1,2,3,4)

data_specific_ir <- data_4ir
plot_list <- list()
for (id in unique(data_specific_ir$ID)) {
  temp.data <- data_specific_ir %>% filter(ID == id)
  title <- paste0(temp.data$Bacteria[1], "\n", temp.data$ID[1], ". Count: ", temp.data$count_id)
  
  p <- ggplot(temp.data, aes(y = Pe_ratio)) +
    geom_boxplot() +
    ggtitle(title) +
    theme_minimal()
  
  plot_list[[length(plot_list) + 1]] <- p
  
  if (length(plot_list) == 8) {
    grid.arrange(grobs = plot_list, ncol = 2, nrow = 4)
    plot_list <- list()
  }
}
plot_list

