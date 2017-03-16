args=(commandArgs(TRUE))
expFile=args[1]
exp=read.table(expFile)
print(paste("mean expression:", round(mean(exp$V2),2)))
print(paste("median expression:", round(median(exp$V2),2)))
print(paste("maximam expression:", round(max(exp$V2),0)))
print(paste("isoforms with exp > 10 TPM:", length(exp$V2[exp$V2>10])))
print(paste("isoforms with exp > 100 TPM:", length(exp$V2[exp$V2>100])))
print(paste("isoforms with exp > 1000 TPM:", length(exp$V2[exp$V2>1000])))
print(paste("isoforms with exp > 10000 TPM:", length(exp$V2[exp$V2>10000])))
pdf("isoformTPM-hist-plot.pdf")
hist(log(exp$V2,10),xlab="Log10 expression",ylab="Isoform count",main="Histogram of log TPM expression", breaks=25)
dev.off()

