#!/usr/bin/env Rscript

# Example input file
# name	position	coverageDepth	numNonConsensus
# 4	4	3265	5
# 5	5	3265	0
# 6	6	3265	1

# graph_aaf.R <input>
# Makes a scatterplot of alternate allele frequency by position
# Takes output of add_consensus_columns_to_frequency_tables.pl or ( add_consensus_columns_to_frequency_tables.pl + merge_tally_overlapping_regions.pl ) 
args <- commandArgs(TRUE)
pdfname <- sub("[.][^.]*$", ".pdf", args[1], perl=TRUE)
#print(pdfname)
data <- read.delim(args[1])
#data$freqNonConsensus = data$numNonConsensus/data$coverageDepth
data$freqNonConsensus = data$numNonConsensus/data$unambigCoverageDepth
meanfreq <- mean(data$freqNonConsensus)
medianfreq <- median(data$freqNonConsensus)
medianfreq.round <- signif(medianfreq, digits=3)
print(paste("Replacing zero with median frequency:", medianfreq.round))
# How do I get the average background??  medianfreq is closer to the average background.
for (i in row.names(data)){if(data[i,"freqNonConsensus"] == 0){ data[i,"freqNonConsensus"] = medianfreq } }		# Replace zero values with median frequency

# Make sure we have a 'position' column.
if (!('position' %in% colnames(data))){
	#print("no position");
	if ('aminoAcidPosition' %in% colnames(data)){
		data$position = data$aminoAcidPosition
	} else if ('nucleotidePosition' %in% colnames(data)){
		data$position = data$nucleotidePosition;
	} else {
		print("no position column")
		quit(save = "no", status = 1, runLast = FALSE)
	}
}


#print(head(data))
library(ggplot2)

pdf(pdfname, width=12, height=8)
#p <- ggplot(data, aes(x=position, y=freqNonConsensus)) + geom_point(stat = "identity") + theme_bw() + coord_trans(y="log") + scale_y_continuous(limits=c(0.00001,1), breaks=c(0.00001,0.0001,0.001,0.01,0.1,1), labels=c("0.00001","0.0001","0.001","0.01","0.1","1"), minor_breaks = c(0.00005,0.0005,0.005,0.05,0.5)) + labs(x="Position", y="Alternate Allele Frequency") + theme(axis.text.x = element_text(size=9, hjust=1, vjust = 0.5, angle=90 )) 
#p + geom_hline(yintercept=medianfreq) + annotate("text", x=max(data$position) * 1.05, y=medianfreq * 1.2, label=paste("Median:", medianfreq.round), size=2) 
if ('merged' %in% colnames(data)){
	p <- ggplot(data, aes(x=position, y=freqNonConsensus, color=merged))
} else {
	p <- ggplot(data, aes(x=position, y=freqNonConsensus))
} 
p <- p + geom_point(stat = "identity") + theme_bw() + coord_trans(y="log") + scale_y_continuous(limits=c(0.00001,1), breaks=c(0.00001,0.0001,0.001,0.01,0.1,1), labels=c("0.00001","0.0001","0.001","0.01","0.1","1"), minor_breaks = c(0.00005,0.0005,0.005,0.05,0.5)) + labs(x="Position", y="Alternate Allele Frequency") + theme(axis.text.x = element_text(size=9, hjust=1, vjust = 0.5, angle=90 )) + scale_colour_manual(values = c("TRUE" = "red","FALSE" = "black"))
p + geom_hline(yintercept=medianfreq) + annotate("text", x=max(data$position) * 1.05, y=medianfreq * 1.2, label=paste("Median:", medianfreq.round), size=2) 
dev.off()
