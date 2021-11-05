#################################################
#  File Name:downsample.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Fri 22 Oct 2021 11:12:17 AM UTC
#################################################

import gzip
import argparse
import numpy as np
import random
import os

parser = argparse.ArgumentParser(usage="python ")
parser.add_argument('-n',"--number",help="the number of reads reserved", default=200000,type=int)
parser.add_argument('-r1',"--reads1",help="the R1 with fastq.gz format", required=True)
parser.add_argument('-r2',"--reads2",help="the R2 with fastq.gz format", required=True)

args = parser.parse_args()
number = args.number
file_R1 = args.reads1
file_R2 = args.reads2

f_R1 = gzip.open(file_R1,'rb')
f_R2 = gzip.open(file_R2,'rb')

fastq_block_R1 = []
fastq_block_R2 = []

def output(file_R1,file_R2,fastq_block_R1_array,fastq_block_R2_array):
    if not os.path.exists('down_sample/'):
        os.makedirs('down_sample')
    f_R1_o = gzip.open('down_sample/'+file_R1,'wb')
    f_R2_o = gzip.open('down_sample/'+file_R2,'wb')
    for reads1,reads2 in zip(fastq_block_R1_array,fastq_block_R2_array):
        f_R1_o.write(reads1[0])
        f_R1_o.write(reads1[1])
        f_R1_o.write(reads1[2])
        f_R1_o.write(reads1[3])

        f_R2_o.write(reads2[0])
        f_R2_o.write(reads2[1])
        f_R2_o.write(reads2[2])
        f_R2_o.write(reads2[3])
    f_R1_o.close()
    f_R2_o.close()

for line in f_R1:
    each_block = [line,next(f_R1),next(f_R1),next(f_R1)]
    fastq_block_R1.append(each_block)
for line in f_R2:
    each_block = [line,next(f_R2),next(f_R2),next(f_R2)]
    fastq_block_R2.append(each_block)
if len(fastq_block_R1) > number:
    all_number = range(1,len(fastq_block_R1))
    selected = random.sample(all_number, number)
    fastq_block_R1_array = np.array(fastq_block_R1)
    fastq_block_R2_array = np.array(fastq_block_R2)
    fastq_block_R1_array1 = fastq_block_R1_array[selected]
    fastq_block_R2_array1 = fastq_block_R2_array[selected]
    output(file_R1,file_R2,fastq_block_R1_array1,fastq_block_R2_array1)
else:
    output(file_R1,file_R2,fastq_block_R1,fastq_block_R2)

