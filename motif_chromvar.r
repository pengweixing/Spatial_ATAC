#################################################
#  File Name:motif_chromvar.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Fri 24 Sep 2021 04:07:02 PM UTC
#################################################

library(chromVAR)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
library(gridExtra)

args=commandArgs(T)

peaks <- getPeaks(args[1])
peak2=resize(sort(peaks),500,fix="center")
#data("homer_pwms")
#data("mouse_pwms_v1")
fragment_counts <- getCounts(args[2],  
                             peak2,
                             paired =  TRUE,
                             by_rg = TRUE,
                             format = "bam",
                             colData = DataFrame(celltype = "brain"))

			  
if (args[3]=="mm"){
fragment_counts <- addGCBias(fragment_counts, genome = "BSgenome.Mmusculus.UCSC.mm10")
jaspar_motifs <- getJasparMotifs("Mus musculus")
counts_filtered <- filterSamples(fragment_counts, min_depth = 4000,
                                 min_in_peaks = 0.01, shiny = FALSE)
counts_filtered <- filterPeaks(sort(counts_filtered))
motif_ix <- matchMotifs(jaspar_motifs, counts_filtered, genome = BSgenome.Mmusculus.UCSC.mm10, out = "scores")
} else if (args[3]=="hs"){
fragment_counts <- addGCBias(fragment_counts, genome = "BSgenome.Hsapiens.UCSC.hg38")
jaspar_motifs <- getJasparMotifs("Homo sapiens")
counts_filtered <- filterSamples(fragment_counts, min_depth = 4000,
                                 min_in_peaks = 0.01, shiny = FALSE)
counts_filtered <- filterPeaks(sort(counts_filtered))
motif_ix <- matchMotifs(jaspar_motifs, counts_filtered, genome =BSgenome.Hsapiens.UCSC.hg38, out = "scores")
} else {
print("Please specify the reference genome")
}


bg <- getBackgroundPeaks(object = counts_filtered)
dev <- computeDeviations(object = counts_filtered,
                         annotations = motif_ix,
                         background_peaks = bg)
deviation=as.data.frame(deviations(dev))
zscore=as.data.frame(deviationScores(dev))
variability <- computeVariability(dev)

name_list = strsplit(rownames(deviation),"_")
name_list2=c()
for(i in 1:length(name_list)){
 name_list2 = c(name_list2,name_list[[i]][2])
   }
rownames(deviation)=name_list2
write.table(deviation,file=paste0(args[4],"_deviation.txt"),quote=F,sep="\t")
tsne_results <- deviationsTsne(dev, threshold = 1, perplexity = 10)
tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation_name = variability$name[1],  sample_column = "celltype",  shiny = FALSE)
pdf(file=paste0(args[4],"_deviationsTsne.pdf"),width=10,height=5)
grid.arrange(tsne_plots[[1]],tsne_plots[[2]],nrow = 1)
dev.off()
save.image(paste0(args[4],'.Rdata'))
