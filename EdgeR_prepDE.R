#=======================================================================================================================================
# These are the R commands used to analyze the RNAseq data (the prepDE.py output from the last step of the shell pipeline)
#=======================================================================================================================================

#=======================================================================================================================================

#Installation of the bioconductor packages needed for this part (only done once, if they're absent):

source("http://www.bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("goseq")

#=======================================================================================================================================
# To start with a clean working environment, you can execute the line below. Some people get really angry seeing this line in a script
# for whatever reason. 

rm(list = ls()) > ls() 

#=======================================================================================================================================

#======================================================================================================================================
# Differential gene expression analysis using EdgeR, starting with the CSV tables from prepDE.
#=======================================================================================================================================

#Loading EdgeR library

library(edgeR)

#Starting prepDE counts files are in this directory:

directory="~/5430_RNAseq/prepDE/"
setwd(directory)

# Importing the csv file gene_count_matrix 

countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id")) 

# Get the name of the columns for the samples

sampleNames <- colnames(countData)

# Take a look at how countData object is displayed

head(countData)

# Creating groups (all samples that are untreated, all samples that are treated) to be used for DGEList object later.
# The group syntax is such as you name the label for the first group and then how many columns are grouped, then same for 2nd group.
# Must be in that order that you see in head(countData). Here, first two columns are the treated ones, next 2 columns are teh untreated ones.

groups=c(rep("treated",2), rep('untreated',2))

# DGEList function creates the DGEList object for the edgeR. We call it cds here

cds <- DGEList(countData, group=groups)

names(cds)   # this displays the components of the cds object

#the cds object right now has the dataset with the counts and also some information about the groups of samples that you have

head(cds$counts)
head(cds$samples)

#Filtering poorly expressed genes

cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >=3, ]
dim(cds)

cds <- calcNormFactors(cds)

cds$samples

#Effective library sizes after normalization

cds$samples$lib.size * cds$samples$norm.factors
cds <- estimateCommonDisp(cds)

names(cds) # now cds has much more stuff than before

cds$common.dispersion
variance <- sqrt(as.numeric(cds$common.dispersion))
variance

cds<-estimateTagwiseDisp(cds, prior.df=10)
plotMDS(cds, method="bcv", col=as.numeric(cds$samples$group))
legend("bottomright", as.character(unique(cds$samples$group)), col=1:3, pch=20, inset=0)       


# Using tagwise dispersion model to call significantly changed genes

de.tgw <- exactTest(cds, dispersion = cds$tagwise.dispersion, pair = c("untreated", "treated"))

names(de.tgw)

de.tgw$comparison
head(de.tgw$table)

resultsTbl.tgw <- topTags(de.tgw, n=nrow(de.tgw$table))$table

#Creating new table with significantly changed genes:

de.genes.tbl.tgw <- resultsTbl.tgw[ resultsTbl.tgw$FDR <= 0.05, ]
dim(de.genes.tbl.tgw)



#MA plot

plot(resultsTbl.tgw$logCPM, resultsTbl.tgw$logFC, pch=20, col='grey60', cex=0.5, main = "Tagwise dispersion", xlab="log(total CPM)", ylab="log(Fold Change)")
points(de.genes.tbl.tgw$logCPM, de.genes.tbl.tgw$logFC, pch=20, col='red', cex=0.75)
