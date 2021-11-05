#################################################
#  File Name:plot_gene_matrix_to_spatial.py 
#  Author: xingpengwei
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se
#  Created Time: Fri 17 Sep 2021 11:00:06 AM UTC
#################################################


import pandas as pd
import argparse
import sys
from draw_spatial import SpatialPlot


parser = argparse.ArgumentParser(usage="python motif_plot.py")
parser.add_argument('-m',"--matrix",help="the deviation/gene sore matrix", required=True)
parser.add_argument('-n',"--name_list",help="the name list of motif/gene within the matrix", required=True)
parser.add_argument('-b',"--barcode",help="the barcode txt", required=True)
parser.add_argument('-cm',"--colormap",help="the colormap for spatial plot", default='jet')
parser.add_argument('-cl',"--color",help="the color for tissue projection", default='red')
parser.add_argument('-im',"--image",help="the tissue section picture .png", required=False)
parser.add_argument('-a',"--alpha",help="the transparency of each pixel",default=1, type=float)
parser.add_argument('-dpi',"--dpi",help="the dpi of output figure", type=float,default=300)
parser.add_argument('-fm',"--format",help="the format of output figure (pdf,png,eps,svg)",default='pdf')
parser.add_argument('-ht',"--height",help="the height of figure", type=int,default=10)
parser.add_argument('-wt',"--width",help="the width of figure",type=int, default=12)
#parser.add_argument('-o',"--output",help="the output name", required=True)
args = parser.parse_args()
ht = args.height
wt = args.width
data = pd.read_table(args.matrix) ####deviation
barcode = pd.read_table(args.barcode,header=None)  ###barcode
name = args.name_list
tissue_image = args.image
alpha = args.alpha
colormap = args.colormap
color = args.color
oformat = args.format
dpi = args.dpi

try:
    f = open(name,'r')
except FileNotFoundError:
    message = 'Your input {genename} is not a gene list, try to use it as a gene name.'
    print(message.format(genename=name))
    all_names = [name]
else:
    print("Your input is a gene list")
    all_names = [line.strip() for line in f.readlines()]

barcodeA = barcode.iloc[:,0]
barcodeB = barcode.iloc[:,2]
barcodeA = barcodeA.dropna()
barcodeB = barcodeB.dropna()

all_motif = data.index
all_motif2 = [line.replace("::","_").replace("(","_").replace(")","") for line in all_motif]
data.index=all_motif2

for each_gene in all_names:
    data2 = pd.DataFrame(data = None,index=barcodeA, columns =barcodeB)
    if not each_gene in list(data.index):
        message = "{gene} does not exit in the matrix, then skip it to next gene"
        print(message.format(gene=each_gene))
        continue
    gene = data.loc[each_gene]
    for each,value in zip(gene.index,gene):
        ba = each.split('_')[0]
        bb = each.split('_')[1]
        if ba in list(barcodeA):
            if bb in list(barcodeB):
                data2.loc[ba,bb]=value
    if tissue_image:
         each_plot = SpatialPlot(data = data2,oformat = oformat,dpi=dpi,height=ht,width=wt,alpha = alpha, offset_r = 0.75, offset_l=0.25)
         each_plot.project_GeneMatrix_to_tissue(data=data2,output=each_gene,image=tissue_image,color=color)
    else:
        each_plot = SpatialPlot(data = data2,height=ht,dpi=dpi,oformat = oformat,width=wt,alpha = alpha,colormap = colormap)
        each_plot.project_GeneMatrix_to_spatial(data = data2,output = each_gene)
    print("%s has been plotted" % each_gene)
   