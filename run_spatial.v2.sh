#################################################
#  File Name:run.sh
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Thu 11 Mar 2021 12:23:25 PM UTC
#################################################

ls *_L001_R1_001.fastq.gz > sample_list

## debarcoding ####################################
echo `date`
echo "start to spatial debarcoding....."
cat /dev/null > dir.list
cat /dev/null > Spatial_qc.stat
echo -e "name\tcell_num\ttotal_reads\tundecoded\tundecoded_rate" >> Spatial_qc.stat
for line in `cat sample_list`
do
R1=$line
R2=`echo $line|sed s#_R1#_R2#g`
out=`echo $R1|sed s/_L001_R1_001.fastq.gz//g`
if [ ! -d $out ];then
mkdir $out
fi
cd $out
echo -e "/disk1/pengweixing/software/anaconda3/bin/python /disk1/pengweixing/pipeline/spatial_atac/debarcode.py -r1 ../$R1 -r2 ../$R2 -b ../barcode.txt" > devarcode.sh
echo -e "echo \"$out has been done\"" >> devarcode.sh
sh devarcode.sh 1>log.o 2>log.e &
cd ..
ls -d $out >> dir.list
done
wait

echo `date`
echo "debarcoding processing has been finished"
echo "start to generate the QC for spatial data...."
##### qc stat #########################################
for line in `cat sample_list`
do
R1=$line
out=`echo $R1|sed s/_L001_R1_001.fastq.gz//g`
num_well=`ls ./$out/*_R1_001.fastq.gz |wc -l`
num_undecoded=`zcat ./$out/undecoded_R1.fastq.gz |grep ^@ |wc -l`
all_num=`zcat $R1|grep ^@ |wc -l`
rate=`echo "scale=4;$num_undecoded/$all_num" |bc`
echo -e "$out\t$num_well\t$all_num\t$num_undecoded\t$rate" >> Spatial_qc.stat.temp
done
echo `date`
echo "QC for spatial data has been generated"


#echo "start to Mapping...."
##### Mapping ####
refrence_genome=$1
if [ ARG"$2" == ARG ];
then processor=10
else
processor=$2
fi

echo $refrence_genome
if [[ "$refrence_genome" == "mm" ]]
then
ref_index=/disk1/pengweixing/database/mm10/bowtie2/mm10
ref_size=/disk1/pengweixing/database/mm10/mm10.chrom.sizes
ref_anno=/disk1/pengweixing/database/mm10/mm10.tss.bed
blacklist=/disk1/pengweixing/database/blacklist/mm10-blacklist.v2.merge.mitochon.bed
elif [[ "$refrence_genome" == "hs" ]] 
then
ref_index=/disk1/pengweixing/database/hg38/index/hg38.fa
ref_size=/disk1/pengweixing/database/hg38/hg38.chrom.sizes
ref_anno=/disk1/pengweixing/database/hg38/hg38.tss.bed
blacklist=/disk1/pengweixing/database/blacklist/hg38-blacklist.v2.merge.mitochon.bed
fi
for line in `cat dir.list`
do
cd $line
echo "###start to mapping" > run_map.sh
for R1 in `ls *_L001_R1_001.fastq.gz`
do
out=`echo $R1|sed s/_L001_R1_001.fastq.gz//g`
R2=`echo $R1|sed s#_R1#_R2#g`
R2_trim=`echo $R2 |sed s#fastq.gz#trim.fastq.gz#g`
R1_trim=`echo $R1 |sed s#fastq.gz#trim.fastq.gz#g`
echo -e "if [ ! -d "01.Mapping" ];then\nmkdir 01.Mapping\nfi" >> run_map.sh
echo -e "if [ ! -d "02.QC" ];then\nmkdir 02.QC\nfi" >> run_map.sh
echo -e "if [ ! -d "03.fragments" ];then\nmkdir 03.fragments\nfi" >> run_map.sh
echo -e "/home/xingqichen/SOFTWARE/Code/00_bin/pyadapter_trim.py -a $R1 -b $R2" >> run_map.sh
echo -e "bowtie2 -p $processor --very-sensitive   -x $ref_index -1 $R1_trim -2 $R2_trim -S ./01.Mapping/$out.sam --rg-id $out" >> run_map.sh
echo -e "awk ' \$3!=\"chrM\" ' ./01.Mapping/$out.sam |samtools view -S -b -f 0x2 -q 10 - |samtools sort -  -o ./01.Mapping/$out.pe.q10.sort.bam" >> run_map.sh
echo -e "samtools index ./01.Mapping/$out.pe.q10.sort.bam" >> run_map.sh
echo -e "samtools flagstat ./01.Mapping/$out.pe.q10.sort.bam > ./01.Mapping/$out.sort.bam.flagstat" >> run_map.sh
echo -e "samtools view -Sb ./01.Mapping/$out.sam > ./01.Mapping/$out.bam" >> run_map.sh
echo -e "samtools sort -@ 10 ./01.Mapping/$out.bam > ./01.Mapping/$out.sort.bam" >> run_map.sh
echo -e "samtools index ./01.Mapping/$out.sort.bam" >> run_map.sh
echo -e "samtools flagstat ./01.Mapping/$out.sort.bam > ./01.Mapping/$out.sort.bam.flagstat" >> run_map.sh
echo -e "samtools view ./01.Mapping/$out.sort.bam chrM > ./01.Mapping/$out.sort.chrM.txt" >> run_map.sh
echo -e "rm ./01.Mapping/$out.bam" >> run_map.sh
echo -e "rm ./01.Mapping/$out.sam" >> run_map.sh
echo -e "java -Xmx4G -jar /home/xingqichen/SOFTWARE/picard-tools-1.119/MarkDuplicates.jar INPUT=./01.Mapping/$out.pe.q10.sort.bam OUTPUT=./01.Mapping/$out.pe.q10.sort.rmdup.bam METRICS_FILE=./01.Mapping/$out.Picard_Metrics_unfiltered_bam.txt  VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true &> ./02.QC/$out.Picard.log"  >> run_map.sh
echo -e "samtools index ./01.Mapping/$out.pe.q10.sort.rmdup.bam" >> run_map.sh
echo -e "/disk1/pengweixing/software/anaconda3/bin/python /disk1/pengweixing/pipeline/spatial_atac/bamtofragments.py -i ./01.Mapping/$out.pe.q10.sort.rmdup.bam -o ./03.fragments/$out.fragments -b $out" >> run_map.sh
echo -e "bedtools intersect -a ./03.fragments/$out.fragments -b $blacklist -v > ./03.fragments/$out.filterBL.fragments" >> run_map.sh
done
bash run_map.sh 2>log.o 1>log.e &
cd ..
done
wait
echo `date`
echo "Mapping has been finished"

######## run bulk #########
echo `date`
echo "Peak calling for bulk ATAC...."

if [ ! -d "Merge_bulk" ];then
mkdir Merge_bulk
fi

for line in `cat dir.list`
do
zcat $line/*_L001_R1_001.trim.fastq.gz |gzip > ./Merge_bulk/${line}_L001_R1_001.trim.fastq.gz
zcat $line/*_L001_R2_001.trim.fastq.gz |gzip > ./Merge_bulk/${line}_L001_R2_001.trim.fastq.gz
done
if [[ "$refrence_genome" == "mm" ]]
then
cp  /disk1/pengweixing/pipeline/spatial_atac/run_bulk.txt ./Merge_bulk/
cd ./Merge_bulk/
nohup bash run_bulk.txt  mm &
cd ..
elif [[ "$refrence_genome" == "hs" ]]
then
cp /disk1/pengweixing/pipeline/spatial_atac/run_bulk.txt ./Merge_bulk/
cd ./Merge_bulk/
nohup bash run_bulk.txt hs &
cd ..
fi
wait
echo `date`
echo "Peak calling for bulk ATAC has been done"

######## FRiP ##########
echo `date`
echo "calculating the FRiP...."
for line in `cat dir.list`
do
cd $line
if [ ! -d "04.FRiP" ];then
mkdir 04.FRiP 
fi
for frag_list in `ls ./03.fragments/*filterBL.fragments`
do
out_name=`echo $frag_list |cut -d '/' -f 3 |sed s#.filterBL.fragments##g`
bedtools intersect -a $frag_list -b ../Merge_bulk/Peak_calling/${line}_L001.filterBL.bed > ./04.FRiP/$out_name.within.peaks.bed
done
cd ..
done
echo `date`
echo "FRiP done"

######  QC stat #####
echo -e "final_reads\tpeak_number" > temp.frag
for line in `cat dir.list`
do
cd $line
cat ./03.fragments/*.filterBL.fragments > $line.filterBL.all.fragments
if [ -s $line.filterBL.all.fragments  ]; then
cat $line.filterBL.all.fragments |sort -k1,1 -k2,2n > $line.filterBL.all.sort.fragments
bgzip -c $line.filterBL.all.sort.fragments > $line.filterBL.all.sort.fragments.gz
fi
#num_peak=`cat ./04.Merge_bulk/${line}_peaks.filterBL.bed |wc -l`
#final_reads=`samtools flagstat ./04.Merge_bulk/${line}.bulk.pe.q10.rmdup.bam |head -1 |cut -d ' ' -f 1`
#echo -e "$final_reads\t$num_peak"  >> ../temp.frag 
cd ..
done
#paste Spatial_qc.stat.temp temp.frag  > Spatial_qc.stat
ls */*.filterBL.all.sort.fragments.gz |sed s#/#\\t#g > fragments.list
#rm Spatial_qc.stat.temp
#rm temp.frag

echo `date`
echo "start to do archR analysis..."
/disk1/pengweixing/software/anaconda3/bin/python /disk1/pengweixing/pipeline/spatial_atac/get_qc_each_spatial.py dir.list



#############  ArchR analysis #########
export R_LIBS="/disk1/pengweixing/software/R-4.1.1/library"
if [[ "$refrence_genome" == "mm" ]]
then
/disk1/pengweixing/software/R-4.1.1/bin/Rscript /disk1/pengweixing/pipeline/spatial_atac/archR.r mm10 fragments.list
elif [[ "$refrence_genome" == "hs" ]]
then
/disk1/pengweixing/software/R-4.1.1/bin/Rscript /disk1/pengweixing/pipeline/spatial_atac/archR.r hg38 fragments.list
fi
echo `date`
echo "All done!"
for line in `cat dir.list`
do
/disk1/pengweixing/software/R-4.1.1/bin/Rscript /disk1/pengweixing/pipeline/spatial_atac/plot_qc.r $line.spatial_QC.txt
done
rm Rplots.pdf
rm *arrow
