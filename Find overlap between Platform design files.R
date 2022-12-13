#============================================
# Find overlap between Platform design files:
#============================================

source("https://bioconductor.org/biocLite.R")
biocLite("ChIPpeakAnno")
library(GenomicRanges)
library(Biostrings)
library(ChIPpeakAnno)
library(HelloRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

# Load GRanges files for each platform
load(file="Platform_Design_GRanges.RData")

# Load the BSgenome package for homo sapiens.
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19


#-----------------------------------
# Make platforms GRangesList object
#----------------------------------

platforms_design.grl <- GRangesList(Agilent_SureSelect.gr, Roche_EpiGiant.gr, Illumina_EPIC.gr, RRBS.gr.sub)
names(platforms_design.grl) <- c("Agilent_SureSelect", "Roche_EpiGiant", "Illumina_EPIC", "RRBS_sub")

# Get target design size of each platform
sum(width(reduce(platforms_design.grl)))

# Get sequences from GRangesList object:  
grl.seq <- getSeq(genome, platforms_design.grl) # The return values are DNAStringSetList objects.

# Count number of occurances for CG in a given region of list object:
lapply(grl.seq, function(x) sum(vcountPattern("CG", x)))


# Find overlaps between several GRanges objects:

# returns an object of overlappingPeaks, which contains there elements: venn_cnt, peaklist (a list of overlapping peaks or unique peaks), and overlappingPeaks (a list of data frame consists of the annotation of all the overlapping peaks).
ol <- findOverlapsOfPeaks(Agilent_SureSelect.gr, Roche_EpiGiant.gr,Illumina_EPIC.gr, RRBS.gr.sub)



# Visualize the number of common and specific peaks with Venn diagram
makeVennDiagram(ol)

# Or in one shot
res <- makeVennDiagram(Peaks=list(Agilent_SureSelect.gr, Roche_EpiGiant.gr, Illumina_EPIC.gr, RRBS.gr.sub), ignore.strand=TRUE,
                       NameOfPeaks=c("Agilent", "Roche", "Illumina", "RRBS"),
                       main="Venn Diagram for platform coverage design",
                       fill=c(1,2,3,4))

###### 4-way diagram using annotated feature instead of chromosome ranges
,scaled=FALSE, euler.d=FALSE

# A pie chart is used to demonstrate the overlap features of the common peaks.
pie1(table(ol$overlappingPeaks[["Agilent_SureSelect.gr///Roche_EpiGiant.gr"]]$overlapFeature)) 




# After finding the overlapping peaks, you can use annotatePeakInBatch to annotate the overlapping peaks with the genomic features in the AnnotationData within certain distance away specified by maxgap, which is 5kb in the following example.


# But if you have them as a single GRanges stratified by type metadata, you could make it GRangesList and plot it like this:
  
# peaks1$type<-"TF1"
# peaks2$type<-"TF2"
# peaks3$type<-"TF3"
# gr <- c(peaks1, peaks2, peaks3) # like your data
# 
# grl <- splitAsList(gr, gr$type)
# res <- makeVennDiagram(Peaks=grl, NameOfPeaks=names(grl))




# If you want a nice and colorful Venn diagram, then check https://github.com/js229/Vennerable/blob/master/README.md to install Venerable (I could not get it from Rforge) and do it like this:
biocLite(c("RBGL","graph"))
install.packages("devtools"); library(devtools);
install_github("js229/Vennerable"); library(Vennerable);
vignette("Venn")

library(Vennerable)
venn_cnt2venn <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts") - 1
  SetNames=colnames(venn_cnt)[1:n]
  Weight=venn_cnt[,"Counts"]
  names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  Venn(SetNames=SetNames, Weight=Weight)
}
v <- venn_cnt2venn(res$vennCounts)
plot(v)

venn_cnt <- res$vennCounts
n <- which(colnames(venn_cnt)=="Counts") - 1
SetNames=colnames(venn_cnt)[1:n]
Weight=venn_cnt[,"Counts"]
names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
v <- Venn(SetNames=SetNames, Weight=Weight)
plot(v)


# Following works for getting the exact intersects between all the ranges.
Reduce(intersect, list(gr, gr1, gr2))