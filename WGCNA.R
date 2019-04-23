#! /usr/bin/env Rscript

library(WGCNA)
options(stringsAsFactors=FALSE)
enableWGCNAThreads()

#Data Input

data=read.table("cohort_deg_gene_tmm",header=TRUE,row.names=1)
datExpr0=as.data.frame(t(data))
gsg=goodSamplesGenes(datExpr0,verbose=3)
gsg$allOK

#Cluster the samples

sampleTree=hclust(dist(datExpr0),method="average")
pdf("Sample_clustering.pdf",width=12,height=9)
plot(sampleTree,main="Sample clustering to detect outliers",sub="",xlab="",cex.lab=1.5,cex.axis=1.5,cex.main=2)
dev.off()

clust=cutreeStatic(sampleTree,cutHeight=40000000,minSize=10)
table(clust)
keepSamples=(clust==1)
datExpr=datExpr0[keepSamples,]
nGenes=ncol(datExpr)
nSamples=nrow(datExpr)

#Loading trait data

datTraits=read.table("sample_trait",header=TRUE,row.names=1)

#Softpower test

powers=c(c(1:10),seq(from=12,to=20,by=2))
sft=pickSoftThreshold(datExpr,powerVector=powers,verbose=5)
pdf("Soft_threshold_power.pdf",width=5,height=5)
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold Power",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main=paste("Scale independence"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=0.9,col="red")
dev.off()

#Co-expression similarity and adjacency

softPower=4
adjacency=adjacency(datExpr,power=softPower)
TOM=TOMsimilarity(adjacency)
dissTOM=1-TOM

#Cluster using TOM

geneTree=hclust(as.dist(dissTOM),method="average")
minModuleSize=200
dynamicMods=cutreeDynamic(dendro=geneTree,distM=dissTOM,deepSplit=2,pamRespectsDendro=FALSE,minClusterSize=minModuleSize)
table(dynamicMods)
dynamicColors=labels2colors(dynamicMods)
pdf("Module_cluster.pdf",width=8,height=6)
plotDendroAndColors(geneTree,dynamicColors,"Dynamic Tree Cut",dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05,main="Gene dendrogram and module colors")
dev.off()

#Merge modules

MEList=moduleEigengenes(datExpr,colors=dynamicColors)
MEs=MEList$eigengenes
MEDiss=1-cor(MEs)
METree=hclust(as.dist(MEDiss),method="average")
pdf("merge_module.pdf")
plot(METree)
dev.off()
MEDissThres=Threshold
merge=mergeCloseModules(datExpr,dynamicColors,cutHeight=MEDissThres,verbose=3)
mergedColors=merge$colors
mergedMEs=merge$newMEs
pdf("geneDendro.pdf",width=12,height=9)
plotDendroAndColors(geneTree,cbind(dynamicColors,mergedColors),c("Dynamic Tree Cut","Merged Dynamic"), dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
dev.off()

#Rename merged variables

moduleColors=mergedColors
colorOrder=c("grey",standardColors(50))
moduleLabels=match(moduleColors,colorOrder)-1
MEs=mergedMEs

#Extract eigengenes

hubs=chooseTopHubInEachModule(datExpr,moduleColors)

#Quantify module-trait associations

MEs0=moduleEigengenes(datExpr,moduleColors)$eigengenes
MEs=orderMEs(MEs0)
moduleTraitCor=cor(MEs,datTraits,use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)

pdf("Module_trait_correlation.pdf",width=10,height=6)
textMatrix=paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
labeledHeatmap(Matrix=moduleTraitCor,xLabels=names(datTraits),yLabels=names(MEs),ySymbols=names(MEs),colorLabels=FALSE,colors=blueWhiteRed(50),textMatrix=textMatrix,setStdMargins=FALSE,cex.text=0.5,zlim=c(-1,1),main=paste("Module_trait relationship"))
dev.off()

#Module trait correlation

age=as.data.frame(datTraits$Age)
names(age)='age'
GS.age=as.numeric(cor(datExpr,age,use="p"))
GS.ageColor=numbers2colors(GS.age,signed=T)
pdf("module_trait.pdf",width=12,height=9)
plotDendroAndColors(geneTree,cbind(dynamicColors,GS.ageColor),c("Modules","Age"),dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
dev.off()

#Extract genes in each module

allLLIDs=names(datExpr)
module="green"
modGenes=(moduleColors==module)
modLLIDs=allLLIDs[modGenes]
fileName=paste(module,sep="")
write.table(as.data.frame(modLLIDs),file=fileName,row.names=FALSE,col.names=FALSE)



