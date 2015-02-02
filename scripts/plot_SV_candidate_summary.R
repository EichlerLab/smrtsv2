library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
input.file <- args[1]
output.file <- args[2]
standard.deviations <- 4

candidates <- read.table(input.file, header=TRUE)

x.limit <- mean(candidates$mean_length) + standard.deviations * sd(candidates$mean_length)

pdf(output.file, width=8, height=6)
ggplot(candidates, aes(x=mean_length, fill=event_type)) + geom_histogram(binwidth=50) + theme_bw() + xlab("Mean SV candidate length (bp)") + ylab("Number of SV candidates") + xlim(c(0, x.limit))
dev.off()
