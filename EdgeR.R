#Installation:

source("http://www.bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("goseq")


#Loading EdgeR library

library(edgeR)

#Starting HTseq counts files are in this directory:

directory="~/5430_RNAseq/HTseq_counts/"
setwd(directory)

#Creating the files and conditions for the table:

sampleFiles <- c("untr_rep1.counts", "untr_rep2.counts", "E2_rep1.counts", "E2_rep2.counts")
sampleConditions <- c("untreated", "untreated", "treated", "treated")

#Building the experiment table

counts=data.frame(matrix(ncol=4, nrow=20366))  #creates empty df to be filled with the information from all 4 samples I have

#filling in the table:

for (ii in sampleFiles)
{
  print(ii)
  temp<- read.table(paste(directory, ii, sep=""), header=F, stringsAsFactors = F)
  counts[,which(sampleFiles==ii)] <- temp[,2]
}

#adding labels to the table
rownames(counts) <- temp[,1]
colnames(counts) <- sampleFiles

head(counts)

# Creating groups (all samples that are untreated, all samples that are treated) to be used for DGEList()

groups=c(rep("untreated",2), rep('treated',2))


# DGEList function creates the object for the edgeR

cds <- DGEList(counts, group=groups)

names(cds)
head(cds$counts)


#Filtering poorly expressed genes

cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >=3, ]
dim(cds)

cds <- calcNormFactors(cds)
cds$samples

#Effective library sizes after normalization

cds$samples$lib.size * cds$samples$norm.factors
cds <- estimateCommonDisp(cds)
names(cds)
cds$common.dispersion

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
