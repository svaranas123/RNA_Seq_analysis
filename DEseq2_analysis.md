### DEseq2 analysis
```
library(DESeq2)

colData = read.csv("colData.csv", header = TRUE, row.names = 1)
countData = read.csv('count_data.csv', header = TRUE, row.names=1)
colnames(countData)= paste0(colData$phenotype)
rownames(colData) = paste0(colData$phenotype)

dds = DESeqDataSetFromMatrix(countData = countData,
                             colData = colData,
                             design = ~ stress)
                             
dim(dds)  ## before filtering
dds = dds[rowSums(counts(dds))>1, ]
dim(dds)
dds = DESeq(dds)
res = results(dds)
res
summary(res)
rld<-rlog(dds, blind = FALSE)
head(assay(rld), 3)
```
### For Volcano plot with labels

```
res =results(dds, contrast=c("stress", "Memory", "Activated"))
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"

volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))
with(subset(res, padj<.005 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
    with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
     with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
    if (labelsig) {
        require(calibrate)
      with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
    legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(res, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()
```
### For heatmap for 20 differentially regulated genes
```
library(pheatmap)
topVarGenes=head(order(rowVars(assay(rld)), decreasing = TRUE), 20)
mat=assay(rld)[topVarGenes, ]
mat=mat-rowMeans(mat)
df=as.data.frame(colData(rld)[,c("stress","phenotype")])
pheatmap(mat,annotation_col = df)
```
### For plotting individual gene
```
geneCounts=plotCounts(dds, gene="Gzmb", intgroup = c("stress"), returnData = TRUE)
ggplot(geneCounts, aes(x=stress, y=count, fill=stress))+scale_y_log10()+geom_dotplot(binaxis="y", stackdir = "center")
```
### For KEGG pathway and GO analysis
```
res =results(dds, contrast=c("stress", "Activated", "Resting"))
res$entrez<-mapIds(org.Mm.eg.db, keys=row.names(res), column="ENTREZID", keytype = "SYMBOL", multiVals = "first")
source("https://bioconductor.org/biocLite.R")
biocLite("pathview")
biocLite("gage")
biocLite("gageData")
library(pathview)
library(gage)
library(gageData)
library(dplyr)
data("kegg.sets.mm")
data("sigmet.idx.mm")
kegg.sets.mm= kegg.sets.mm[sigmet.idx.mm]
head(kegg.sets.mm, 3)
foldchanges=res$log2FoldChange
names(foldchanges) =res$entrez
head(foldchanges)
keggres=gage(foldchanges, gsets = kegg.sets.mm, same.dir = TRUE)
lapply(keggres, head)

keggrespathways=data.frame(id=rownames(keggres$greater), keggres$greater) %>% tbl_df() %>% filter(row_number()<=5) %>% .$id %>% as.character()
keggrespathways

data(go.sets.mm)
data("go.subs.mm")
gobsets=go.sets.mm[go.subs.mm$BP]
gopres=gage(foldchanges, gsets = gobsets, same.dir = TRUE)
lapply(gopres, head)
```
### For PCA plots
```
rld_pca <- function (rld, intgroup = "stress", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=10, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}
png("qc-pca.png", 1000, 1000, pointsize=15)
rld_pca(rld, colors=mycols, intgroup="stress", xlim=c(-95, 55))
dev.off()
```
