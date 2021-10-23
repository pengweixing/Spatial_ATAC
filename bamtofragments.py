#################################################
#  File Name:bamtofragments.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Thu 02 Sep 2021 02:38:24 PM UTC
#################################################

import sys
import pysam
import argparse
import collections

def process(bam_input,MapQ):
    f_bam_input = pysam.AlignmentFile(bam_input,mode='rb',check_header=True)
    fragDict = collections.OrderedDict()
    for each in f_bam_input:
        if each.mapping_quality > MapQ:
            if each.is_proper_pair:
                if each.template_length > 0:
                    if not each.is_reverse:
                        #### shift ###
                        fragStart = each.reference_start + 5
                        fragend = fragStart + each.template_length -9
                        theChr = each.reference_name
                        if (theChr,fragStart,fragend) in fragDict.keys():
                            fragDict[(theChr,fragStart,fragend)] = fragDict[(theChr,fragStart,fragend)] + 1
                        else:
                            fragDict[(theChr,fragStart,fragend)] = 1
    return fragDict

def main():
    parser = argparse.ArgumentParser(usage="python ")
    parser.add_argument('-i',"--input",help="the input of bam ", required=True)
    parser.add_argument('-o',"--output",help="the name of output ", default='Unnamed')
    parser.add_argument('-b',"--barcode",help="the barcode of this bam file", required=True)
    parser.add_argument('-mq',"--mapquality",help="the map quality value [0-60]",type=int,default=10)
    args = parser.parse_args()
    bam_input = args.input
    bam_output = args.output
    barcode = args.barcode
    MapQ = args.mapquality
    fragDict =  process(bam_input,MapQ)
    with open(bam_output,'w') as f:
        for key,value in fragDict.items():
            print(key[0],key[1],key[2],barcode,value,sep="\t",file=f)

if __name__ == "__main__":
    main() 
