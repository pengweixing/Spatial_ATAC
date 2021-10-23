#################################################
#  File Name:plot_qc.r
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Tue 14 Sep 2021 03:56:50 PM UTC
#################################################

library(ggplot2)
library(ggpubr)
args = commandArgs(T)
data = read.table(args[1],header = T)
p=ggscatter(data,x="fragments_filtered_with_blacklist",y="FRiP",size = 1,xlab="The number of unique fragments",ylab="FRiP")
out= gsub(".txt","",args[1])
ggsave(paste0(out,".pdf"),plot=p,width=6,height=4)
