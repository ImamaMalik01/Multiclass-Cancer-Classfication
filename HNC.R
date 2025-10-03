intall.packages("limma")
install.packages("edgeR")
file_path <- 'SKC_R.csv'
df <- read.csv(file_path)
head(df)
names(df)
duplicated(df$X)
length(df$X)
unique(df$X)
length(unique(df$X))
library(tidyverse)
df %>% distinct()
df1 <-df %>% distinct()
length(df1)
dim(df1)
df[!duplicated(df$X), ]
df1<-df[!duplicated(df$X), ]

row.names(df1)<- df1$X
row.names(df1)
class(df1)
dim(df1)
head(df1)
df1<- df1[,-116]
head(df1)
tail(df1)
library(edgeR)
library(limma)
summary(df1)
sum_size <- colSums(df1[, -1])
print(sum_size)
expt<-rep(c("T","N"),c(9,43))
expt
expt<- factor(expt, levels = c("T","N")) # ensure the order
expt

rawcounts <- df1
rawcounts[rawcounts<0] <- 0
df2<-rawcounts[rawcounts<0] <- 0
sum_size
names(df2)
class(df2)
class(rawcounts)
rawcounts < 0
rawcounts[rawcounts<0] <- 0
rawcounts
rawcounts < 0
count(rawcounts < 0)
rawcounts[rawcounts<0] <- 0
sum_size <- colSums(rawcounts)
class(rawcounts)
write.csv(rawcounts, "rawcounts_SKC.csv")
rawcounts <- rawcounts[-1]
sum_size <- colSums(rawcounts)
d<- DGEList(counts=as.matrix(rawcounts), lib.size=sum_size, group=expt)
class(d)
names(d)
d$counts
d$samples
d<- calcNormFactors(d)
d$samples
d<- estimateCommonDisp(d)
nbt<- exactTest(d,pair=c("N","T"))
names(nbt)
nbt$table
# The p.values have not been adjusted for multiple hypothesis testing; This can be done with the topTags function, similar to the topTable limma function. It will perform FDR (False Discovery Rate correction)
nbt.corrected<- topTags(nbt,n=Inf)
nbt.corrected #see deg for top 5 genes
#### This result can be output by
write.csv(nbt.corrected$table,file="NBT_SKC.csv")
# How many genes are significant at FDR p < 0.05?
sum(nbt.corrected$table$FDR<0.05)
nbttable<- nbt.corrected$table
deg<- nbttable[nbttable$FDR<0.05,]
head(deg)
dim(deg)
deg$ID<-rownames(deg)
# separate up-regulated and down-regulated genes
deg_up<- deg[(deg$FDR<0.05&deg$logFC>0),]
deg_up<- deg_up[order(deg_up$logFC, decreasing=T),]
dim(deg_up)
head(deg_up)
deg$ID<-rownames(deg)
head(deg)
# separate up-regulated and down-regulated genes
deg_up<- deg[(deg$FDR<0.05&deg$logFC>0),]
deg_up<- deg_up[order(deg_up$logFC, decreasing=T),]
dim(deg_up)
head(deg_up)
deg_down<- deg[(deg$FDR<0.05&deg$logFC<0),]
deg_down<- deg_down[order(deg_down$logFC),]
dim(deg_down)
head(deg_down)
#Output your deg
write.csv(deg_up,file="UPRegulated_SKC.csv")
write.csv(deg_down,file="DownRegulated_SKC.csv")
# combine deg with raw counts
rawcounts$ID<-rownames(rawcounts)
head(rawcounts)
combined<-merge(rawcounts, deg, by='ID', all.ID=T)
head(combined)
deg_up<- deg_up[order(deg_up$logFC, decreasing=T),]
combined<-combined[order(combined$logFC, decreasing=T),]
head(combined)
write.csv(combined, "Combined_SKC.csv")
library(dplyr)
library(pheatmap)
library(gplots)

# Read the CSV file
combined <- read.csv("combined_CRC.CSV")

# Convert columns to appropriate data types
numeric_combined <- combined %>%
  mutate(across(where(is.numeric), as.numeric),
         across(where(is.character), as.character))

# Extract gene IDs and numeric values
gene_ids <- numeric_combined$ID
numeric_values <- as.matrix(numeric_combined[, -1])

# Select a subset of rows and columns (adjust indices as needed)
subset_data <- numeric_values[1:20, 2:30]
subset_gene_ids <- gene_ids[1:20]

# Create a heatmap using the heatmap.2 function
heatmap.2(
  subset_data,
  Rowv = FALSE,  # Disable row clustering
  Colv = FALSE,  # Disable column clustering
  dendrogram = "none",  # Do not show dendrograms
  trace = "none",  # Do not show trace lines
  margins = c(5, 10),  # Adjust margins
  labRow = subset_gene_ids,  # Use gene IDs as row labels
  main = "Heatmap of CRC",
  key = TRUE,  # Show color key
  keysize = 1.5,  # Adjust size of the color key
  density.info = "none"  # Do not show density plot
)

head(combined)
dim(combined)
install.packages("dplyr")
install.packages("viridis")
library(viridis)
library(gplots)
# Install and load the 'profvis' package
#install.packages("profvis")
library(profvis)
# Install and load the 'dplyr' package
# install.packages("dplyr")
library(dplyr)
library(pheatmap)

# Assuming you have combined, gene_ids, and numeric_values as defined in your code

# Convert character columns to numeric
numeric_combined <- combined %>%
  mutate(across(where(is.numeric), as.numeric),
         across(where(is.character), as.character))

# Extract gene IDs and numeric values
gene_ids <- numeric_combined$ID
numeric_values <- as.matrix(numeric_combined[, -1])  # Exclude the 'ID' column

# Create a data frame for row annotations
annotation_df <- data.frame(combined[,-1])

# Define custom color palette
background_color <- "#f0f0f0"  # Light grey
expression_color <- "red"

# Create a heatmap using pheatmap
pheatmap(numeric_values,
         annotation_row = annotation_df,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c(background_color, expression_color))(50),
         breaks = seq(min(numeric_values), max(numeric_values), length = 51),
         main = "Heatmap of CRC")

# Save the heatmap as an image file (e.g., PNG, PDF, JPEG)
ggsave("heatmap_output_custom_colors.png", width = 800, height = 600, units = "px", dpi = 300)
library(gplots)

# Read data from your file (replace 'your_data_file.csv' with your actual file name)
your_data <- read.csv("combined_CRC.CSV", header = TRUE, row.names = 1)

# Extract gene names for row labels
row_labels <- rownames(your_data)

# Specify the colors for heatmap
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Create the heatmap
heatmap.2(as.matrix(your_data),
          Colv = FALSE,                  # Do not cluster columns
          Rowv = FALSE,                  # Do not cluster rows
          dendrogram = "none",           # Do not show dendrograms
          col = my_palette,              # Set the color palette
          trace = "none",                # Do not show trace lines
          margins = c(10, 6),            # Set margins for row and column names
          labRow = row_labels,           # Use gene names as row labels
          key = TRUE,                    # Show color key
          keysize = 1.5,                 # Set size of the color key
          cexRow = 0.8,                  # Adjust size of row labels
          cexCol = 0.8,                  # Adjust size of column labels
          col_labels = row_labels        # Use gene names as column labels
)
numeric_combined <- combined %>%
  mutate(across(where(is.numeric), as.numeric),  # Convert numeric columns to numeric
         across(where(is.character), as.character))  # Convert character columns to character

# Extract gene IDs and numeric values
gene_ids <- numeric_combined$ID
numeric_values <- as.matrix(numeric_combined[, -1])  # Exclude the 'ID' column

# Check the structure of numeric_values
str(numeric_values)

# Create a heatmap using the numeric matrix
heatmap.2(numeric_values, Rowv = NA, Colv = NA,
        labRow = gene_ids, main = "Heatmap of CRC")

# Save the heatmap as an image file (e.g., PNG, PDF, JPEG)
# Choose the appropriate file type based on your needs
# For example, save as PNG:
png("heatmap_output.png", width=800, height=600, units="px", res=300)
# Or save as PDF:
# pdf("heatmap_output.pdf", width=8, height=6)
# Or save as JPEG:
# jpeg("heatmap_output.jpg", width=800, height=600, units="px", res=300)

# Close the device after saving
dev.off()
write.csv(combined, "DEG.csv")
write.csv(combined, "DEG_CRC.csv")
install.packages("pheatmap")
pheatmap(
  data_matrix,
  color = colorRampPalette(c("darkblue", "white", "darkred"))(100),  # Define the color palette
  cluster_rows = TRUE,  # Cluster rows
  cluster_cols = TRUE,  # Cluster columns
  annotation_row = gene_ids,  # Display gene IDs as row annotations
  annotation_col = colnames(data_matrix),  # Display sample names as column annotations
  main = "Heatmap Title",  # Add a title
  fontsize_row = 6,  # Adjust font size for row annotations
  fontsize_col = 6,  # Adjust font size for column annotations
  cellwidth = 15,  # Adjust cell width
  cellheight = 8,  # Adjust cell height
  show_colnames = TRUE,  # Show column names
  show_rownames = TRUE,  # Show row names
  border_color = NA,  # Remove border color
  gaps_col = 1,  # Add a small gap between columns
  gaps_row = 1   # Add a small gap between rows
)
