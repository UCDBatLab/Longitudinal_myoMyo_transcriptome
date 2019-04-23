library(DESeq2)

countData=read.table("6_vs_7+",header=TRUE,row.names=1)
colData=DataFrame(condition=factor(c(rep("A",3),rep("B",4))))
dds=DESeqDataSetFromMatrix(countData,colData,formula(~condition))
dds=DESeq(dds)
res=results(dds)
resOrdered=res[order(res$padj),]
write.table(resOrdered,"DESeq2_6_vs_7+")

