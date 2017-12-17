#=======================================================================================================================================
# These are the R commands used to analyze the RNAseq data (the prepDE.py output from the last step of the shell pipeline).
# For details, check out Final.Rnw (or the PDF)
#=======================================================================================================================================

#=======================================================================================================================================

#Installation of the bioconductor packages needed for this part (only done once, if they're absent):

source("http://www.bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("goseq")
biocLite("org.Hs.eg.db")
#=======================================================================================================================================
# To start with a clean working environment, you can execute the line below. Some people get really angry seeing this line in a script for whatever reason. 

rm(list = ls()) > ls() 

#=======================================================================================================================================

#Loading required libraries
library(car)
library(dplyr)
library(DESeq2)
library(edgeR)
library(goseq)
library(org.Hs.eg.db)
library(rtracklayer)

#Starting prepDE counts files are in this directory:

directory="~/5430_RNAseq/prepDE/"
setwd(directory)

# Importing the csv file gene_count_matrix 

countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
head(countData)
dfcountData <- data.frame(countData)


#=======================================================================================================================================
# Plotting replicates
#=======================================================================================================================================

plot(log2(dfcountData[,1]), log2(dfcountData[,2]), 
     pch=20, 
     col='gray60',
     cex=0.2, 
     main="Log2 scatterplot of the E2 replicates", 
     xlab="Log2(E2 rep 1 counts)", 
     ylab="Log2(E2 rep 2 counts)")

# When adding the fit line and finding the correlation, there is a problem for genes whose counts are zero. To fix this, I added 0.1 to the counts number for those commands (does not alter the countData matrix).

model1 <- lm(log2(dfcountData[,1]+0.1) ~ log2(dfcountData[,2]+0.1))
abline(model1, col="red")
model1_cor <- cor(log2(dfcountData[,1]+0.1), log2(dfcountData[,2]+0.1))
model1_cor
text(14, 1, paste("R squared =", round(model1_cor, digits=4)))

# Same for the untreated replicates 1 and 2:

plot(log2(dfcountData[,3]), log2(dfcountData[,4]), 
     pch=20, 
     col='gray60',
     cex=0.2,
     main="Log2 scatterplot of the untreated replicates", 
     xlab="Log2(untr_rep 1 counts)",
     ylab="Log2(untr_rep 2 counts)")

model2 <- lm(log2(dfcountData[,3]+0.1) ~ log2(dfcountData[,4]+0.1))
abline(model2, col="red")

model2_cor <- cor(log2(dfcountData[,3]+0.1), log2(dfcountData[,4]+0.1))
model2_cor
text(14, 1, paste("R squared =", round(model2_cor, digits=4)))


#=======================================================================================================================================
# Differential Gene Expression analysis 
#=======================================================================================================================================


# Get the name of the columns for the samples

sampleNames <- colnames(countData)

# Take a look at how countData object is displayed

head(countData)

# Creating groups (all samples that are untreated, all samples that are
# treated) to be used for DGEList object later. The group syntax is such
# that you name the label for the first group and then how many columns
# are grouped, then same for 2nd group. Must be in that order that you 
# see in head(countData). Here, first two columns are the treated ones,
# next 2 columns are teh untreated ones.

groups=c(rep("treated",2), rep('untreated',2))

# DGEList function creates the DGEList object for the edgeR. We call it cds here

cds <- DGEList(countData, group=groups)

names(cds)   

# this displays the components of the cds object
# the cds object right now has the dataset with the counts and also some
# information about the groups of samples that you have

head(cds$counts)
head(cds$samples)

# Filtering poorly expressed genes
cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size,
                                                       dim(cds)) > 1) >=3, ]
dim(cds)
cds <- calcNormFactors(cds)
cds$samples

#Effective library sizes after normalization

cds$samples$lib.size * cds$samples$norm.factors
cds <- estimateCommonDisp(cds)
cds <- estimateTagwiseDisp(cds, prior.df = 10)

names(cds) # now cds has much more info than before

# View common dispersion
cds$common.dispersion

# View tagwise dispersion values for several datapoints
head(cds$tagwise.dispersion)

tag_variance <- sqrt(as.numeric(cds$tagwise.dispersion))

# View tagwise variance values for several datapoints
head(tag_variance)
common_variance <- sqrt(as.numeric(cds$common.dispersion))

# View overall variance
common_variance

# Using tagwise dispersion model to call significantly changed genes
de.tgw <- exactTest(cds, dispersion = cds$tagwise.dispersion, 
                    pair = c("untreated", "treated"))
head(de.tgw$table)
resultsTbl.tgw <- topTags(de.tgw, n=nrow(de.tgw$table))$table
head(resultsTbl.tgw)

#Creating new table with significantly changed genes - here entries with
# FDR values above 0.05 are being discarded
de.genes.tbl.tgw <- resultsTbl.tgw[ resultsTbl.tgw$FDR <= 0.05, ]

# This displays how many genes have FDR <= 0.05
dim(de.genes.tbl.tgw)

#=======================================================================================================================================
# PCA plots
#=======================================================================================================================================

condition <- c( "treated", "treated", "untreated","untreated")

colData <- data.frame(condition)

row.names(colData) <- sampleNames

colData

all(rownames(sampleNames) %in% colnames(countData))

countData <- countData[, sampleNames]

all(rownames(colData) == colnames(countData))

dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData,
                              design = ~ condition)
dds<- DESeq(dds)
res <- results(dds, alpha=0.05)
rld <- rlogTransformation(dds, blind=TRUE)
plotPCA(rld, intgroup=c('condition'))

#=======================================================================================================================================
# MA/MDS plots 
#=======================================================================================================================================

# The MA plot using plot function

plot(resultsTbl.tgw$logCPM, resultsTbl.tgw$logFC,
     pch=20, 
     col='grey60',
     cex=0.5,main = "Tagwise dispersion MA plot using plot() and points()",
     xlab="log2(total CPM)", 
     ylab="log2(Fold Change)")

# adding horizontal lines for y = -1 and 1
abline(h=-1, lty=5)
abline(h=1, lty=5)

points(de.genes.tbl.tgw$logCPM, de.genes.tbl.tgw$logFC, 
       pch=20, 
       col='red',
       cex=0.75)

# The MA plot using EdgeR's plotSmear function
  
de.genes.tgw <- rownames(resultsTbl.tgw)[resultsTbl.tgw$FDR <= 0.05]

plotSmear(cds, de.tags=de.genes.tgw, 
          main="Tagwise dispersion MA plot using plotSmear", 
          pair=c("untreated", "treated"), 
          cex=0.35, 
          xlab="log CPM", 
          ylab="Log fold-change")

abline(h=-1, lty=2)
abline(h=1, lty=2)

write.table(resultsTbl.tgw , file="all_genes_tgw.txt", col.names =T, row.names=T)

resultsTbl.tgw <- read.table("all_genes_tgw.txt")

# Filter based on FDR

de.genes.tgw <- resultsTbl.tgw[ resultsTbl.tgw$FDR <= 0.05,]

dim(de.genes.tgw)

E2_up <- de.genes.tgw [de.genes.tgw $logFC >0,]
E2_down <- de.genes.tgw [de.genes.tgw $logFC <0,]

dim(E2_up)
dim(E2_down)

Genes.all <- rownames(resultsTbl.tgw)
E2_up_names <- rownames(E2_up)
E2_down_names <- rownames(E2_down)

# The above gene names still have the transcript variant attached after the "." separator. In order to obtain the gene name, they need to be split:

Genes.all_split <- sapply(Genes.all, function(x) unlist(strsplit(x, split="[.]"))[1])

E2up_split <-sapply(E2_up_names, function(x) unlist(strsplit(x, split="[.]"))[1])

E2down_split <-sapply(E2_down_names, function(x) unlist(strsplit(x, split="[.]"))[1])

write.table(as.data.frame(Genes.all_split), file="Genes_all_split.txt", row.names=TRUE, col.names=TRUE)

write.table(as.data.frame(E2up_split), file="E2_upregulated_split.txt", row.names=TRUE, col.names=TRUE)

write.table(as.data.frame(E2down_split), file="E2_downregulated_split.txt", row.names=TRUE, col.names=TRUE)

head(Genes.all_split)


#=======================================================================================================================================
# GO analysis 
#=======================================================================================================================================

  
up.vector=as.integer(Genes.all_split%in%E2up_split)

down.vector=as.integer(Genes.all_split%in%E2down_split)

head(up.vector)

# Add names from all genes to vector

names(up.vector)<-Genes.all_split

names(down.vector)<-Genes.all_split

head(up.vector)

pwf_up=nullp(up.vector,"hg19","ensGene")


pwf_down=nullp(down.vector,"hg19","ensGene")

head(pwf_up)

GO.wall_up=goseq(pwf_up,"hg19","ensGene", method = "Wallenius")

GO.wall_down=goseq(pwf_down,"hg19","ensGene", method = "Wallenius")

head(GO.wall_up)

head(GO.wall_down)


enriched.GO.BP_up=GO.wall_up$category[p.adjust(GO.wall_up$over_represented_pvalue, method="BH")<0.05]

#The FDR of 0.01 would not give many results here, so I had to use 0.05 and worse.

enriched.GO.BP_down=GO.wall_down$category[p.adjust(GO.wall_down$over_represented_pvalue, method="BH")<0.1]

head(enriched.GO.BP_up)

head(enriched.GO.BP_down)

library(GO.db)

for(go in enriched.GO.BP_up[1:5])
  {
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
    }

#=======================================================================================================================================
# Comparison with ChIP
#=======================================================================================================================================

#loading the BED file into R
E2_all_closest <- read.table("ER_E2_closest_all.bed", header=FALSE, sep="\t", stringsAsFactors = FALSE)

head(E2_all_closest)

#loading the lists of names in R

all_gene_names <- read.table("Gene_names_all.txt", header=TRUE, stringsAsFactors = FALSE)
upreg_names <- read.table("Upregulated_names_sorted.txt", header=TRUE, stringsAsFactors = FALSE)
downreg_names <- read.table("Downregulated_names_sorted.txt", header=TRUE, stringsAsFactors = FALSE)

# To subset the E2_all_closest, I first selected the rows on the condition that V9 (the ENSG name) is within the gene list names, then picked those rows from E2_all_closest

selected_rows <- (E2_all_closest$V9 %in% upreg_names$E2upsplit)

Upreg_closest <- E2_all_closest[selected_rows,]

#this is the dataframe with only upregulated genes

selected_rows <- (E2_all_closest$V9 %in% downreg_names$E2downsplit) 

Downreg_closest <- E2_all_closest[selected_rows,] 

selected_rows <- (E2_all_closest$V9 %in% all_gene_names$Genesallsplit)

Allgenes_closest <- E2_all_closest[selected_rows,]

# To obtain the non-regulated table, I then removed the regulated genes rows from the E2_all_closest dataset 

Regulated_closest <- union(Upreg_closest, Downreg_closest)

Nonregulated_closest <- setdiff(Allgenes_closest, Regulated_closest)

up_closest_cdf <- ecdf(abs(Upreg_closest[,12]))
down_closest_cdf <- ecdf(abs(Downreg_closest[,12]))
non_closest_cdf <- ecdf(abs(Nonregulated_closest[,12]))

plot(up_closest_cdf, 
     xlim=c(0,50000), 
     ylim=c(0.1, 0.9), 
     col="forestgreen", 
     main = "CDF plot of ER E2 distance to TSS", 
     xlab="Distance (bps)", 
     ylab = "Fraction", 
     las=1)

par(new=TRUE)

plot(down_closest_cdf, 
     xlim=c(0,50000), 
     ylim=c(0.1, 0.9),
     col="red",
     main = "CDF plot of ER E2 distance to TSS", 
     xlab="Distance (bps)",
     ylab = "Fraction",
     las=1)

par(new=TRUE)

plot(non_closest_cdf, 
     xlim=c(0,50000), 
     ylim=c(0.1, 0.9), 
     col="black", 
     main = "CDF plot of ER E2 distance to TSS", 
     xlab="Distance (bps)", 
     ylab = "Fraction", 
     las=1)                       

abline(h=0.5, lty=2)
abline(h=0.2, lty=2)

legend( x="topleft", 
        legend=c("Upregulated genes", "Downregulated genes", "Non-regulated genes"), 
        col=c("forestgreen","red", "black"),
        lwd=1,
        lty=c(1,1,1), 
        pch=c(NA,NA) )

Upreg_closest$effect<-"Upregulated"
Downreg_closest$effect<-"Downregulated"
Nonregulated_closest$effect<-"Nonregulated"

Genes_with_effects<- rbind(Upreg_closest, Downreg_closest, Nonregulated_closest)

boxplot(abs(Genes_with_effects$V12) ~ Genes_with_effects$effect, 
        outline=F, 
        border=c("red", "grey40", "forestgreen"),
        main="Distribution of distances between ChIP peaks and gene TSS",
        lty=1)
