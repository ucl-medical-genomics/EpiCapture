

# Import metric on target capture efficiencies for 3 hybridization platforms:
hyb_capture.df <- read.delim("clipboard", header = TRUE, stringsAsFactors = TRUE )

dim(hyb_capture.df)
str(hyb_capture.df)

save(hyb_capture.df,file = "S:/EpiCapture/RESULTS/3_TargetCaptureMetrics/hyb_capture.RData")
load("S:/EpiCapture/RESULTS/3_TargetCaptureMetrics/hyb_capture.RData")

# Arrange for 3 plots - 1 per row/ 3 per row
par(mfrow=c(1,1))

# plot for each Platform
names(hyb_capture.df)
Agilent="Agilent \n SureSelect MethylSeq"
Roche= "Roche \n SeqCap EpiGiant"
Illumina="Illumina \n TruSeq Methyl Cature EPIC"

head(hyb_capture.df[,c(1:4,21:28)])
pct_cov.df <- hyb_capture.df[,c(1:4,21:28)]

head(pct_cov.df)
library(dplyr)

x <- c(1,2,10,20,30,40,50,100)

#install.packages("wesanderson")  # Install
#library(wesanderson)


#-------------------
# Agilent:
pct_cov.df %>% filter(Platform=="Agilent") -> y

rownames(y) <- y[,3]
colnames(y)
y <- y[,5:12]

# set color pallete:

names(wes_palettes)
col = wes_palette("Zissou1", 16, type = "continuous")

# plot the first line:
plot(x,y[1,], type="l", lty=1, xaxt='n', col= col[1], ylab = "Fraction of capture targets >= Depth", xlab = "Depth of coverege", main = Agilent)

# plot x-axis
depth <- c( "1x", "2x","10x","20x", "30x","40x","50x", "100x")
axis(side = 1, at = x, labels = depth, tcl = -0.2)
abline(v=x, col="gray", lty=2)
abline(h=0.8, col="red", lty=1)
abline(h=0.6, col="blue", lty=1)

# Add a second line
for (i in 2:16){
  lines(x, y[i,], type = "l", lty=1, col= col[i])
}

# Add a legend to the plot
legend("topright", legend=rownames(y),
       col=col, lty = 1, cex=0.7)
#----------------
# Roche:

pct_cov.df %>% filter(Platform=="Roche") -> y

rownames(y) <- y[,3]
colnames(y)
y <- y[,5:12]

# set color pallete:
colors <- RColorBrewer::brewer.pal(16,)

names(wes_palettes)
col = wes_palette("Zissou1", 16, type = "continuous")

# plot the first line:
plot(x,y[1,], type="l", lty=1, xaxt='n', col= col[1], ylab = "Fraction of targets", xlab = "Depth of coverege", main = Roche)

# plot x-axis
depth <- c( "1x", "2x","10x","20x", "30x","40x","50x", "100x")
axis(side = 1, at = x, labels = depth, tcl = -0.2)
abline(v=x, col="gray", lty=2)
abline(h=0.8, col="red", lty=1)
abline(h=0.6, col="blue", lty=1)

# Add a second line
for (i in 2:16){
  lines(x, y[i,], type = "l", lty=1, col= col[i])
}

# Add a legend to the plot
legend("topright", legend=rownames(y),
       col=col, lty = 1, cex=0.7)


#------------
# Illumina: 
pct_cov.df %>% filter(Platform=="Illumina") -> y

rownames(y) <- y[,3]
colnames(y)
y <- y[,5:12]

# set color pallete:
colors <- RColorBrewer::brewer.pal(16,)

names(wes_palettes)
col = wes_palette("Zissou1", 16, type = "continuous")

# plot the first line:
plot(x,y[1,], type="l", lty=1, xaxt='n', col= col[1], ylab = "Fraction of targets", xlab = "Depth of coverege", main = Illumina)

# plot x-axis
depth <- c( "1x", "2x","10x","20x", "30x","40x","50x", "100x")
axis(side = 1, at = x, labels = depth, tcl = -0.2)
abline(v=x, col="gray", lty=2)
abline(h=0.8, col="red", lty=1)
abline(h=0.6, col="blue", lty=1)

# Add a second line
for (i in 2:16){
  lines(x, y[i,], type = "l", lty=1, col= col[i])
}

# Add a legend to the plot
legend("topright", legend=rownames(y),
       col=col, lty = 1, cex=0.7)


# Plot boxplot/barplot with error bars:


