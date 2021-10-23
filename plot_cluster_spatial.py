#################################################
#  File Name:plot_cluster_spatial.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Fri 17 Sep 2021 11:00:06 AM UTC
#################################################

# cluster_matrix <- 
# BC1     BC2     cluster
# Y18     X11     C3
# Y18     X14     C1
# Y18     X13     C2

# barcode txt
# Y1      AACGTGAT        X1      AGAGTCAA
# Y2      AAACATCG        X2      ACTATGCA
## Y3      ATGCCTAA        X3      ACGTATCA

# color code list
# colour    group
# #272E6A   C2
# #89288F   C4
# #208A42   C3
# #D51F26   C1

import pandas as pd
import draw_spatial
import argparse
import sys
import numpy as np

parser = argparse.ArgumentParser(usage="python plot_cluster_spatial.py")
parser.add_argument('-m',"--cluster_matrix",help="the cluster for each index", required=True)
parser.add_argument('-b',"--barcode",help="the barcode txt", required=True)
parser.add_argument('-c',"--color_code_list",help="the color code for each group", required=True)
#parser.add_argument('-o',"--output",help="the output name", required=True)
args = parser.parse_args()

data = pd.read_table(args.cluster_matrix)
barcode = pd.read_table(args.barcode,header=None)  ###barcode
color = pd.read_table(args.color_code_list,header=0)  ###barcode

barcodeA = barcode.iloc[:,0]
barcodeB = barcode.iloc[:,2]
barcodeA = barcodeA.dropna()
barcodeB = barcodeB.dropna()

data2 = pd.DataFrame(data = -100,index=barcodeA, columns =barcodeB)
for i in range(data.shape[0]):
    BC_X = data["BC2"][i]
    BC_Y = data["BC1"][i]
    each_cl = data["cluster"][i]
    data2.loc[BC_Y,BC_X] = each_cl

draw_spatial.draw_spatial_cluster(data2,color)
