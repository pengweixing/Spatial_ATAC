#################################################
#  File Name:plot_cluster_spatial.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com,pengwei.xing@igp.uu.se
#  Created Time: Fri 17 Sep 2021 11:00:06 AM UTC
#################################################
"""
 cluster_matrix  --cluster_matrix
 BC1     BC2     cluster
 Y18     X11     C3
 Y18     X14     C1
 Y18     X13     C2

 barcode txt --barcode
 The first line and the third line are required and it should be in order.
 Y1      AACGTGAT        X1      AGAGTCAA
 Y2      AAACATCG        X2      ACTATGCA
 Y3      ATGCCTAA        X3      ACGTATCA

 color code list --color_code_list
 colour    group
 #272E6A   C2
 #89288F   C4
 #208A42   C3
 #D51F26   C1

 optional
--image 
tissue.png 

"""

import pandas as pd
from draw_spatial import SpatialPlot
import argparse
import sys
import numpy as np

parser = argparse.ArgumentParser(usage="python plot_cluster_spatial.py")
parser.add_argument('-m',"--cluster_matrix",help="the cluster for each index (required)", required=True)
parser.add_argument('-b',"--barcode",help="the barcode txt (required)", required=True)
parser.add_argument('-im',"--image",help="the tissue section picture .png (optional)", required=False)
parser.add_argument('-wc',"--which_cluster",help="which cluster do you want to plot on the tissue picture,e.g. C1 or C1,C3", default='All')
parser.add_argument('-c',"--color_code_list",help="the color code for each group", required=True)
parser.add_argument('-o',"--output",help="the output name (required)", required=True)
parser.add_argument('-fm',"--format",help="the format of output figure (pdf,png,eps,svg)",default='pdf')
parser.add_argument('-dpi',"--dpi",help="the dpi of output figure", type=float,default=300)
parser.add_argument('-a',"--alpha",help="the transparency of each pixel, only works on tissue image",default=1,type=float)
parser.add_argument('-ht',"--height",help="the height of figure", type=int,default=10)
parser.add_argument('-wt',"--width",help="the width of figure",type=int, default=10)
args = parser.parse_args()


data = pd.read_table(args.cluster_matrix)
barcode = pd.read_table(args.barcode,header=None)  ###barcode
color = pd.read_table(args.color_code_list,header=0)  ###barcode
output = args.output
ht = args.height
wt = args.width
which_cluster = args.which_cluster
tissue_image = args.image
alpha = args.alpha
oformat = args.format
dpi = args.dpi

if which_cluster == 'All':
    all_cluster = list(color['group'])
elif ',' in which_cluster:
    all_cluster = which_cluster.split(',')
else:
    all_cluster = which_cluster

barcodeA = barcode.iloc[:,0]
barcodeB = barcode.iloc[:,2]
barcodeA = barcodeA.dropna()
barcodeB = barcodeB.dropna()

data2 = pd.DataFrame(data = None,index=barcodeA, columns =barcodeB)
for i in range(data.shape[0]):
    BC_X = data["BC2"][i]
    BC_Y = data["BC1"][i]
    each_cl = data["cluster"][i]
    data2.loc[BC_Y,BC_X] = each_cl
if tissue_image:
    plot = SpatialPlot(data = data2,oformat = oformat,height=ht,width=wt,dpi=dpi,alpha = alpha,offset_r = 0.75, offset_l=0.25)
    plot.project_cluster_to_tissue(data=data2,color_list=color,output=output,image=tissue_image,which_cluster=all_cluster)
else:
    plot = SpatialPlot(data = data2,height=ht,width=wt,dpi=dpi,oformat = oformat)
    plot.project_cluster_to_spatial(data2,color,output,all_cluster)
