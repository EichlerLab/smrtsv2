library(ggplot2)
library(gridExtra)

args <- commandArgs(trailingOnly=TRUE)
input.file <- args[1]
output.file <- args[2]

data <- read.table(input.file, header=TRUE)
aligned.data <- data[data$alignment_status == "mapped",]

median.subread.length <- median(data$subread_length)
sd.subread.length <- sd(data$subread_length)

median.aligned.length <- median(aligned.data$aligned_length)
sd.aligned.length <- sd(aligned.data$aligned_length)

p1 <- ggplot(data, aes(x=subread_length)) + geom_histogram(binwidth=500, fill="darkgreen") + geom_vline(xintercept=median.subread.length, colour="black") + xlab("Subread length (bp)") + ylab("Frequency") + theme_bw() + annotate("text", x = median.subread.length + 2 * sd.subread.length, y = 1500000, label = paste("Median:", median.subread.length))

p2 <- ggplot(aligned.data, aes(x=aligned_length)) + geom_histogram(binwidth=500, fill="darkgreen") + geom_vline(xintercept=median.aligned.length, colour="black") + xlab("Aligned length (bp)") + ylab("Frequency") + theme_bw() + annotate("text", x = median.aligned.length + 2 * sd.aligned.length, y = 1500000, label = paste("Median:", median.aligned.length))

pdf(output.file, width=12, height=5)
grid.arrange(arrangeGrob(p1, p2, ncol=2))
dev.off()

