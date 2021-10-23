#################################################
#  File Name:get_qc_each_spatial.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Mon 13 Sep 2021 08:28:48 PM UTC
#################################################

import sys
import os
import gzip
import numpy as np

sample_list = open(sys.argv[1],'r')## dir.list
all_sample = [each.strip() for each in sample_list]

for each_sampe in all_sample:
    all_cell = [cell.strip() for cell in os.popen("ls %s/*_R1_001.trim.fastq.gz" % each_sampe)]
    all_name = [each.split('/')[1].replace('_L001_R1_001.trim.fastq.gz','') for each in all_cell]
    f = open(each_sampe+".spatial_QC.txt",'w') 
    header = "Name\tTotal_reads\tMapping_rate\tMapped_pair_reads\tChrM_rate\tDup_rate\tfragments\tfragments_filtered_with_blacklist\tfrag_size_median\tfrag_size_mean\tFRiP\n"
    f.write(header)
    for each_cell,each_name in zip(all_cell,all_name):
        with open(each_sampe+"/01.Mapping/"+each_name+".sort.bam.flagstat") as f_mprate:
            all_line = [line.strip() for line in f_mprate.readlines()]
            map_rate = all_line[4].split()[4][1:]
            total_reads = all_line[0].split()[0]
        with open(each_sampe+"/01.Mapping/"+each_name+".sort.bam.flagstat") as f_mprate:
            all_line = [line.strip() for line in f_mprate.readlines()]
            map_pair = all_line[8].split()[0]
        with open(each_sampe+"/01.Mapping/"+each_name+".sort.chrM.txt") as f_chrM:
            all_line = [line.strip() for line in f_chrM.readlines()]
            chrM_rate = "%.2f%%" % (int(len(all_line))/int(total_reads)/2*100)
        with open(each_sampe+"/02.QC/"+each_name+".Picard.log") as f_dup:
            all_line = [line.strip() for line in f_dup.readlines()]
            dup_num = [line.split() for line in all_line if 'Marking' in line][0][5]
            dup_rate = "%.2f%%" % (int(dup_num)/int(total_reads)/2*100)
        with open(each_sampe+"/03.fragments/"+each_name+".fragments") as f_frag:
            all_frag = [line.strip() for line in f_frag.readlines()]
            frag_num = len(all_frag)
        with open(each_sampe+"/03.fragments/"+each_name+".filterBL.fragments") as f_fragBL:
            all_fragBL = [line.strip().split() for line in f_fragBL.readlines()]
            frag_BL = len(all_fragBL)
            if frag_BL > 0:
                all_fragBL_np = np.array(all_fragBL)
                frag_size = all_fragBL_np[:,2].astype(int)-all_fragBL_np[:,1].astype(int)
                frag_size_median = np.median(frag_size)
                frag_size_mean ="%.1f" % np.mean(frag_size)
            else:
                frag_size_median = 0
                frag_size_mean = 0
        with open(each_sampe+"/04.FRiP/"+each_name+".within.peaks.bed") as f_frip:
            all_frip = [line.strip() for line in f_frip.readlines()]
            if frag_BL > 0:
                frip_num = "%.5f" % (len(all_frip)/frag_BL)
            else:
                frip_num = 0
        print(each_name,total_reads,map_rate,map_pair,chrM_rate,dup_rate,frag_num,frag_BL,frag_size_median,frag_size_mean,frip_num,sep="\t",file=f)
