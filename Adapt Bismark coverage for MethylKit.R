# Adapt coverage files for import by MethylKit


# read Bismark covrage file
df <- fread.gzipped(paste0(path2Nugen_methcall, "/", "NuGEN_ZYMO-UM_R1_val_1_trimmed_bismark_bt2_pe_n6dupsRemoved_NameSorted.bismark.cov.gz"), data.table=FALSE)
head(df)

# set m:
min.cov <- 3
assembly <- "hg19"
context <- "CpG"
sample.id <- "ZYMO"

# remove low coverage stuff
df=df[(df[,5]+df[,6]) >= min.cov ,]



# make the object (arrange columns of df), put it in a list
result[[1]] = new("methylRaw",data.frame(chr=df[,1],start=df[,2],end=df[,3],
                                        strand="*",coverage=(df[,5]+df[,6]),
                                        numCs=df[,5],numTs=df[,6]),
                 sample.id=sample.id[[1]],
                 assembly=assembly,context=context,resolution="base")
                 


# You can read the methylation ratio file by using "pipeline" argument in read() function. 
# You need to provide a list of column numbers corresponding to chr,start,end,strand,coverage and ratio of methylation. 
# Actually, you can read any generic methylation ratio or percentage file using this method. 
# The file needs to have the location information (chr,start,end and strand), coverage information and methylation percentage or ratio information.

library(methylKit)

obj=read("bsmap_output.txt",sample.id="test",assembly="mm9",header=TRUE, context="CpG", resolution="base",
         pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                       coverage.col=6,strand.col=3,freqC.col=5 )
)

obj

obj=read("bismark.txt",sample.id="test",assembly="hg19",header=TRUE, 
         context="CpG", resolution="base",
         pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                       coverage.col=6,strand.col=3,freqC.col=5 )
)


#----------------------------------------------------------------------

# Simon's script:
library(devtools)
install_github("al2na/methylKit",ref="development",build_vignettes=FALSE)
library(methylKit)

myobj=read("/Data/Simon/data/hiseq1/samSort/Agilent_Hela_1.sorted_CpG_prefix2.methylKit",
           pipeline=list(fraction=TRUE,chr.col=2,start.col=3,end.col=3,coverage.col=5,strand.col=4,freqC.col=6),
           sample.id="agilentHela1",
           assembly="hg18",
           treatment=1)

library("graphics")
getMethylationStats(myobj, plot = T, both.strands = F)

getCoverageStats(myobj, plot = T, both.strands = F)

dev.print(pdf, '/Data/Simon/methylKit/filename.pdf')

# Used MethylDackel to generate comatible files:
