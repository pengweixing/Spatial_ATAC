#################################################
#  File Name:debarcode.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Wed 10 Mar 2021 11:03:39 PM UTC
#################################################

import sys
import gzip
import argparse
import collections

def fargv():
    parser = argparse.ArgumentParser(usage="python debarcode.py -r1 R1.fastq.gz -r2 R2.fastq.gz -b barcode.list -o output_name")
    parser.add_argument('-r1',"--R1",help="the R1 of file ", required=True)
    parser.add_argument('-r2',"--R2",help="the R2 of file ", required=True)
    parser.add_argument('-pA',"--pos_barcodeA",help=" the position of barcode A,default is [0-8)", default="0-8")
    parser.add_argument('-pB',"--pos_barcodeB",help="the position of barcode B, default is [38-46)", default="38-46")
    parser.add_argument('-l',"--pos_of_R1",help="the position of genomic DNA in R1, default is 95", default="95",type=int)
    parser.add_argument('-b',"--barcode",help="the barcode list ", required=True)
    parser.add_argument('-o',"--output_name",help="the name of output ", required=True)
    args = parser.parse_args()
    return args

def generate_10xbarcode(spatial_barcode):
    f_10x = open('/disk1/pengweixing/software/cellranger-atac-2.0.0/barcodes/737K-cratac-v1.txt','r')
    f_output = open('spatial_barcode.map.10x.txt','w')
    all_barcodes = spatial_barcode
    all_10x = [line.strip() for line in f_10x.readlines()]
    all_len = len(all_barcodes)
    barcode_dict = {}
    b_index = 1
    for i in range(all_len):
        barcodeA = all_barcodes[i]
        barcodeA = barcodeA.split()
        for j in range(all_len):
            barcodeB = all_barcodes[j]
            barcodeB = barcodeB.split()
            print(barcodeA[0],barcodeA[1],barcodeB[2],barcodeB[3],all_10x[b_index],sep="\t",file=f_output)
            barcode_dict[barcodeA[0]+"_"+barcodeB[2]] = all_10x[b_index]
            b_index += 1
    f_10x.close()
    f_output.close()
    return barcode_dict

class barcode_space:
    pass

def gener_bar_list(all_barcode):
    mybarcode_space = barcode_space()
    mybarcode_space.BarcodeA_id = []
    mybarcode_space.BarcodeA_seq = []
    mybarcode_space.BarcodeB_id = []
    mybarcode_space.BarcodeB_seq = []
    for line in all_barcode:
        line1 = line.split()
        mybarcode_space.BarcodeA_id.append(line1[0])
        mybarcode_space.BarcodeA_seq.append(line1[1])
        mybarcode_space.BarcodeB_id.append(line1[2])
        mybarcode_space.BarcodeB_seq.append(line1[3])
    return mybarcode_space

def split_file(all_barcode_list,R1_rds,R2_rds,pA,pB,pos_R1,bar_10x,outname):
    mynames = collections.OrderedDict()
    mybarcode_space = all_barcode_list
    BarcodeA_id = mybarcode_space.BarcodeA_id
    BarcodeA_seq = mybarcode_space.BarcodeA_seq
    BarcodeB_id = mybarcode_space.BarcodeB_id
    BarcodeB_seq = mybarcode_space.BarcodeB_seq
    count=1
    pA_s = int(pA.split('-')[0])
    pA_e = int(pA.split('-')[1])
    pB_s = int(pB.split('-')[0])
    pB_e = int(pB.split('-')[1])
    undecoded_R1 = gzip.open("undecoded_R1.fastq.gz","wb")
    undecoded_R2 = gzip.open("undecoded_R2.fastq.gz","wb")
    cache_output = 200000
    index_output = 0
    while 1:
        if index_output % cache_output == 0 and index_output > 0:
            output(mynames,BarcodeA_id,BarcodeB_id,bar_10x,outname)
            mynames = collections.OrderedDict()
            index_output = 0
        else:
            index_output += 1
            p1_line = R1_rds.readline()
            p2_line = R2_rds.readline()
            if not p1_line:
                break
            if count ==1:
                seqhead1 = p1_line
                seqhead2 = p2_line
            elif count ==2:
                seq1 = p1_line.rstrip()
                seq1_barA = seq1[pA_s:pA_e].decode()
                seq1_barB = seq1[pB_s:pB_e].decode()
                seq1_genome = seq1[pos_R1:]
                seq2 = p2_line.rstrip()
            elif count ==3:
                qualhead1 = p1_line
                qualhead2 = p2_line
            elif count ==4:
                qual1 = p1_line.rstrip()
                qual1_genome = qual1[pos_R1:]
                qual2 = p2_line.rstrip()

                for A in range(len(BarcodeA_id)):
                    each_barcodeA = BarcodeA_seq[A]
                    each_barcodeA_id = BarcodeA_id[A]
                    if each_barcodeA == seq1_barA:
                        index_A = A
                        break
                else:
                    undecoded_R1.write(seqhead1);undecoded_R1.write(seq1+b"\n")
                    undecoded_R1.write(qualhead1);undecoded_R1.write(qual1+b"\n")
                    undecoded_R2.write(seqhead2);undecoded_R2.write(seq2[0:55]+b"\n")
                    undecoded_R2.write(qualhead2);undecoded_R2.write(qual2[0:55]+b"\n")
                    count = count + 1
                    if count == 5:
                        count = 1
                    else:
                        count = count
                    continue
                 
                for B in range(len(BarcodeB_id)):

                    each_barcodeB = BarcodeB_seq[B]
                    each_barcodeB_id = BarcodeB_id[B]
                    if each_barcodeB == seq1_barB:
                        
                        index_B = B
                        break
                else:
                    undecoded_R1.write(seqhead1);undecoded_R1.write(seq1+b"\n")
                    undecoded_R1.write(qualhead1);undecoded_R1.write(qual1+b"\n")
                    undecoded_R2.write(seqhead2);undecoded_R2.write(seq2[0:55]+b"\n")
                    undecoded_R2.write(qualhead2);undecoded_R2.write(qual2[0:55]+b"\n")
                    count = count + 1
                    if count == 5:
                        count = 1
                    else:
                        count = count
                    continue

                temp_id = BarcodeA_id[index_A] + "_" + BarcodeB_id[index_B]
                if temp_id in mynames:
                    mynames[temp_id].append([seqhead1,seq1_genome,qualhead1,qual1_genome,seqhead2,seq2[0:55],qualhead2,qual2[0:55]])
                else:
                    mynames[temp_id] = [[seqhead1,seq1_genome,qualhead1,qual1_genome,seqhead2,seq2[0:55],qualhead2,qual2[0:55]]]

            count = count + 1
            if count == 5:
                count = 1
            else:
                count = count
    if len(mynames) >0:
        output(mynames,BarcodeA_id,BarcodeB_id,bar_10x,outname)        
    
def output(mynames,BarcodeA_id,BarcodeB_id,bar_10x,outname):
    allcombination = [i+"_"+j for i in BarcodeA_id for j in BarcodeB_id]
    f_I = gzip.open(outname+'_S1_L001_I1_001.fastq.gz','ab')
    f1 = gzip.open(outname+'_S1_L001_R1_001.fastq.gz','ab')
    f2 = gzip.open(outname+'_S1_L001_R2_001.fastq.gz','ab')
    f3 = gzip.open(outname+'_S1_L001_R3_001.fastq.gz','ab')
    for each in allcombination:
        if each in mynames:
            all_lines = mynames[each]
            barcode_10x = bar_10x[each]
            for each_block in all_lines:
                R1_block = each_block[0:4]
                R2_block = each_block[4:8]
                f_I.write(R1_block[0]);f_I.write(b'AAAGCATA\n')
                f_I.write(b'+\n');f_I.write(b'FFFFFFFF\n')
                f1.write(R1_block[0]);f1.write(R1_block[1]+b'\n')
                f1.write(R1_block[2]);f1.write(R1_block[3]+b'\n')
                f2.write(R1_block[0])
                f2.write(barcode_10x.encode()+b'\n')
                f2.write(b'+\n');f2.write(b'FFFFFFFFFFFFFFFF\n')
                f3.write(R2_block[0]);f3.write(R2_block[1]+b'\n')
                f3.write(R2_block[2]);f3.write(R2_block[3]+b'\n')
    f_I.close()
    f1.close()
    f2.close()
    f3.close()

def main(kwargs):
    args = kwargs
    R1_input = args.R1
    R2_input = args.R2
    pA=args.pos_barcodeA
    pB=args.pos_barcodeB
    pos_R1 = args.pos_of_R1
    barcode_list = args.barcode
    outname = args.output_name
    f_barcode = open(barcode_list,'r')
    all_barcode = [line.strip() for line in f_barcode.readlines()]
    bar_10x = generate_10xbarcode(all_barcode)
    all_barcode_list = gener_bar_list(all_barcode)
    R1_rds = gzip.open(R1_input,'rb')
    R2_rds = gzip.open(R2_input,'rb')
    split_file(all_barcode_list,R1_rds,R2_rds,pA,pB,pos_R1,bar_10x,outname)

if __name__ == "__main__":
    kwargs = fargv()
    main(kwargs)

