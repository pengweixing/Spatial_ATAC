#################################################
#  File Name:motif.plot.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Fri 17 Sep 2021 11:00:06 AM UTC
#################################################


import pandas as pd
import draw_spatial
import argparse
import sys

parser = argparse.ArgumentParser(usage="python motif_plot.py")
parser.add_argument('-d',"--deviation",help="the deviation matrix", required=True)
#parser.add_argument('-m',"--motif",help="the name of motif within the deviation matrix", required=True)
parser.add_argument('-b',"--barcode",help="the barcode txt", required=True)
parser.add_argument('-min',"--vmin",help="min value after normlized", default=-1,type=int)
#parser.add_argument('-o',"--output",help="the output name", required=True)
args = parser.parse_args()

data = pd.read_table(args.deviation) ####deviation
barcode = pd.read_table(args.barcode,header=None)  ###barcode
vmin = args.vmin

barcodeA = barcode.iloc[:,0]
barcodeB = barcode.iloc[:,2]
barcodeA = barcodeA.dropna()
barcodeB = barcodeB.dropna()

all_motif = data.index
all_motif2 = [line.replace("::","_").replace("(","_").replace(")","") for line in all_motif]
data.index=all_motif2

for each_motif in all_motif2:
    data2 = pd.DataFrame(data = -100,index=barcodeA, columns =barcodeB)
    motif = data.loc[each_motif]
    for each,value in zip(motif.index,motif):
        ba = each.split('_')[0]
        bb = each.split('_')[1]
        if ba in list(barcodeA):
            if bb in list(barcodeB):
                data2.loc[ba,bb]=value
    draw_spatial.draw_spatial(data2,each_motif,vmin)
#data2.to_csv("motif.txt",sep="\t")
