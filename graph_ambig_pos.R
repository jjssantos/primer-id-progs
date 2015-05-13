#!/usr/bin/env Rscript
# Example Input:
# Position	count	Gap
# 1	10	F
# 2	19	F
# 3	1	F
# 4	20	F
# 5	1	F
# 6	1	F
# 7	11	F
# 8	5	F
# 9	7	F
# 
# graph_ambig_pos.R <input>
# Makes a barplot of ambiguous nucleotide counts by position

args <- commandArgs(TRUE)
pdfname <- sub("[.][^.]*$", ".pdf", args[1], perl=TRUE)
#print(pdfname)
data <- read.delim(args[1])
#print(head(data))
library(ggplot2)
max <- subset(data, data$count == max(data$count) )
#print(max)
pdf(pdfname, width=12, height=8)
ggplot(data, aes(x=Position, y=count, color=Gap)) + geom_bar(stat = "identity") + annotate("text", x=max$Position, y=max$count * 1.02, label=paste("(", max$Position, ", ", max$count, ")", sep="")) + labs(title="Ambiguous Nucleotides Found in Consensus Reads", y="Count", x="Position") + scale_color_manual(values = c("black", "red"))
dev.off()
