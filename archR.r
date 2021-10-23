#################################################
#  File Name:archR.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Thu 09 Sep 2021 09:36:53 AM UTC
#################################################


library(ArchR)
args=commandArgs(T)
addArchRThreads(threads = 20)
addArchRGenome(args[1])
sample <- read.table(args[2])
list1 <- paste(sample$V1,sample$V2,sep="/")
names(list1)=sample$V1
ArrowFiles <- createArrowFiles(inputFiles = list1, sampleNames = names(list1),minTSS = 0,minFrags = 0, addTileMat = TRUE,addGeneScoreMat = TRUE)
projHeme1 <- ArchRProject( ArrowFiles = ArrowFiles,outputDirectory = "project",copyArrows = TRUE )
new_list = list1
tss_data = data.frame()
for(each in names(new_list)){
idxSample <- BiocGenerics::which(projHeme1$Sample %in% each)
cellsSample <- projHeme1$cellNames[idxSample]
temp <- plotTSSEnrichment(ArchRProj = projHeme1[cellsSample, ])
tss_data = rbind(tss_data,temp$data)
assign(each,temp)
}
pal <- paletteDiscrete(values = unique(tss_data$group))
p_tss <- ggplot(tss_data, aes(x, v, color = group)) + geom_line(size = 1) + theme_ArchR()+ 
  xlab("Distance From Center (bp)") +  ylab("Normalized Insertion Profile")+ 
  scale_y_continuous(limits = c(0, max(tss_data$v)*1.05, expand = c(0, 0))) + 
  scale_x_continuous(limits = c(min(tss_data$x),   max(tss_data$x)), expand = c(0, 0))+ 
  scale_color_manual(values = pal) 
ggsave("all_tss_enrichent.pdf",plot=p_tss,width = 6,height = 6)

frag_data = data.frame()
for(each in names(new_list)){
idxSample <- BiocGenerics::which(projHeme1$Sample %in% each)
cellsSample <- projHeme1$cellNames[idxSample]
temp <- plotFragmentSizes(ArchRProj = projHeme1[cellsSample, ])
frag_data = rbind(frag_data,temp$data)
assign(each,temp)
}
pal <- paletteDiscrete(values = unique(frag_data$group))
p_frag <- ggplot(frag_data, aes(fragmentSize, fragmentPercent, color = group)) + geom_line(size = 1) + theme_ArchR()+ 
  xlab("ATAC-seq Fragment Size (bp)") +  ylab("Percentage of Fragments")+ 
  scale_y_continuous(limits = c(0, max(frag_data$fragmentPercent)*1.05, expand = c(0, 0))) + 
  scale_x_continuous(limits = c(min(frag_data$fragmentSize),   max(frag_data$fragmentSize)), expand = c(0, 0))+ 
  scale_color_manual(values = pal) 
ggsave("all_frag_length_dist.pdf",plot=p_frag,width = 6,height = 6)

