library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
input.file <- args[1]
lengths.output.file <- args[2]
support.output.file <- args[3]
standard.deviations <- 4

candidates <- read.table(input.file, header=TRUE)

x.limit <- mean(candidates$mean_length) + standard.deviations * sd(candidates$mean_length)
pdf(lengths.output.file, width=8, height=6)
ggplot(candidates, aes(x=mean_length, fill=event_type)) + geom_histogram(binwidth=50) + theme_bw() + xlab("Mean SV candidate length (bp)") + ylab("Number of SV candidates") + xlim(c(0, x.limit))
dev.off()

x.limit <- max(candidates$support)
pdf(support.output.file, width=8, height=6)
ggplot(candidates, aes(x=support, fill=event_type)) + geom_histogram(binwidth=1) + theme_bw() + xlab("SV candidate support") + ylab("Number of SV candidates") + xlim(c(0, x.limit))
dev.off()
