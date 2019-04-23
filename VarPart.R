library(variancePartition)
library(limma)
library(edgeR)

Matrix=read.table("ptcd_cohort_count",header=TRUE,row.names=1)
info=read.table("info",header=TRUE,row.names=1)
gExpr=DGEList(counts=Matrix)
gExpr=calcNormFactors(gExpr)
design=model.matrix(~ Age, info)
vobjGenes=voom(gExpr,design)
form=~(1|Individual)+(1|Site)+(1|Year)+(1|Batch)+Age
varPart=fitExtractVarPartModel(vobjGenes,form,info)
write.table(round(varPart,2),"varPart")

table=read.table("varPart",header=TRUE)
library(reshape2)
melt=melt(table,vars="Gene")

ggplot(melt,aes(x=variable,y=value,fill=variable))+geom_boxplot(aes(color=variable))+background+labs(x="",y="Variance explained %")+theme(legend.title=element_blank())+theme(axis.text.x=element_text(angle=45,hjust=1))+scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")+stat_summary(geom="crossbar",color="black",width=0.75,fun.data=function(x){return(c(y=median(x),ymin=median(x),ymax=median(x)))})+theme(legend.position="none")

