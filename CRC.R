#raw.data<- read.csv("CRC_R.csv",header=T, row.names=1)
#class(raw.data)# Because it's a data frame, we can ask what the dimensions are:
#dim(raw.data)  # how many rows and columns we have
#head(raw.data) # Look at the first few (6) rows of the data
# Replace 'your_file.csv' with the actual file path
file_path <- 'CRC_R.csv'

# Read the CSV file into a data frame
df <- read.csv(file_path)
head(df)
# Fill missing values with 0
df[is.na(df)] <- 0

# Save the modified data frame back to the CSV file
write.csv(df, file = file_path, row.names =1)

cat("Missing values filled with 0 and file saved successfully.\n")
class(df1)
dim(df1)
head(df1)
tail(df1)
library(edgeR)
library(limma)
library(affy)
library(affycoretools)
summary(df1)
sum_size <- colSums(rawcounts)
print(sum_size)
expt<-rep(c("T","N"),c(71,43))
expt
expt<- factor(expt, levels = c("T","N")) # ensure the order
expt
library(edgeR)
rawcounts <- df1
 rawcounts[rawcounts<0] <- 0
d<- DGEList(counts=as.matrix(rawcounts), lib.size=sum_size, group=expt)
class(d)
names(d)
d$counts
d$samples
d<- calcNormFactors(d)
d$samples
plotMDS.DGEList(d,main="MDS plot")
d<- estimateCommonDisp(d)
nbt<- exactTest(d,pair=c("N","T"))
names(nbt)
nbt$table
# The p.values have not been adjusted for multiple hypothesis testing; This can be done with the topTags function, similar to the topTable limma function. It will perform FDR (False Discovery Rate correction)
nbt.corrected<- topTags(nbt,n=Inf)
nbt.corrected[1:5,] #see deg for top 5 genes
#### This result can be output by
write.csv(nbt.corrected$table,file="NBT.csv")
# How many genes are significant at FDR p < 0.05?
sum(nbt.corrected$table$FDR<0.05)
nbttable<- nbt.corrected$table
deg<- nbttable[nbttable$FDR<0.05,]
head(deg)
dim(deg)   
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
write.csv(deg_up,file="UPRegulated.csv")
write.csv(deg_down,file="DownRegulated.csv")
# combine deg with raw counts
rawcounts$ID<-rownames(rawcounts)
head(rawcounts)
combined<-merge(rawcounts, deg, by='ID', all.ID=T)
head(combined)
deg_up<- deg_up[order(deg_up$logFC, decreasing=T),]
combined<-combined[order(combined$logFC, decreasing=T),]
head(combined)
heatmap(as.matrix(combined), Colv=NA, Rowv=NA)
write.csv(combined, "DEG.csv")