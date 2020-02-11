#' ---
#' title: "Class 5 Winter 20"
#' author: "Kiana Munoz"
#' date: "January 21st, 2020"
#' ---

# Class 5
# Data Visualization and graphs in R

plot(1:10, col = "blue", typ = "o")

# Need to import/read input data file first
baby <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)

# Plot weight vs age
plot(baby$Age, baby$Weight, type = "o", 
     pch = 15, cex = 1.5, lwd = 2, 
     ylim = c(2,10), xlab = "Age", 
     ylab = "Weight (kg)", main = "Baby Weight Vs Age",
     col = "blue")

# Barplot data with 'feature_counts' text
mouse <- read.table("bimm143_05_rstats/feature_counts.txt", header = TRUE, sep = "\t")
par(mar = c(5, 11, 2, 1))
barplot(mouse$Count, horiz = TRUE, col = "darkgreen", names.arg = mouse$Feature, las = 1, ylab = "")

# Section 3A
color <- read.table("bimm143_05_rstats/male_female_counts.txt", header= TRUE, sep = "\t")
barplot(color$Count, names.arg = color$Sample, col = c("blue", "white"), las = 2)
