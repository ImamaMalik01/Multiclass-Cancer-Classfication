setwd("C:/Users/HP/Desktop/R_Analysis")
file_path <- 'SKC_R.csv'
data<- read.csv(file_path)
head(data)
names(data)
duplicated(data$X)
length(data$X)
unique(data$X)
length(unique(df$X))
library('tidyverse')
df %>% distinct()
data1 <- data %>% distinct()
length(data1)
dim(data1)
data[!duplicated(data$X)]

