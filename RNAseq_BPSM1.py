#!/bin/bash
cp -r /localdisk/data/BPSM/ICA1/* .
mv fastq/* .
mv Tcongo_genome/* .

#1 quality analysis
# use fastqc installed package to check sequence quality
find *fq.gz | parallel -j 48 fastqc #run 48 jobs in parallel

#2 quality summary
# decompress fastqc result
find *.zip | parallel -j 48 unzip 
# fastqc result summary
echo "sample_file" > fastfilename.txt
echo "PASS" > fastfile-pass.txt
echo "FAIL" > fastfile_fail.txt
echo "WARN" > fastfile-warn.txt
echo "Failed_item" > fast_failname.txt
# count the number of passes fails and warns for each sequence file
for fastfile in  *fastqc
do
echo ${fastfile} >>fastfilename.txt
grep -c PASS ${fastfile}/summary.txt >> fastfile-pass.txt
grep -c FAIL ${fastfile}/summary.txt >> fastfile_fail.txt
grep -c WARN ${fastfile}/summary.txt >> fastfile-warn.txt
if test $(grep -c FAIL ${fastfile}/summary.txt) -ne 0   # grab failed item
then
fail=$(grep FAIL ${fastfile}/summary.txt | cut -f 2)
echo ${fail} >> fast_failname.txt
else
 echo "NA" >> fast_failname.txt
 fi
done
paste fastfilename.txt fastfile-pass.txt fastfile-warn.txt fastfile_fail.txt fast_failname.txt > fastfile_result.txt
# move result to new directory: all_result_folder/quality_result 
mkdir all_result_folder
mkdir all_result_folder/quality_result
mv *fastqc.html all_result_folder/quality_result/.
mv *fastqc all_result_folder/quality_result/.
mv fastfile_result.txt all_result_folder/quality_result/.
echo -e "quality analysis finished.\n"

#3 alignment
find *fq.gz |parallel -j 48 gzip -d # uncompress fq.gz sequence files
gzip -d *fasta.gz # uncompress genome file
# align the read pairs to the Trypanosoma congolense genome using the installed programme bowtie2
hisat2-build -p 16 TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta genome # create index
for file in *1.fq
do
# -p number of threads, -x genome index, -s output sam
hisat2 -x genome -1 ${file} -2 ${file:0:9}2.fq -S ${file:0:8}.sam & 
# do parallel
done
wait
# wait till above finished
echo -e "sam files generated.\n"
for file in *sam
do
#convert the output "sam" to indexed "bam" format with samtools
samtools view -b -S ${file} > ${file:0:8}.bam & 
done
wait
echo -e "bam files generated.\n"
for file in *.bam
do
# sort the lines in bam files into the same order
samtools sort ${file:0:8}.bam -o ${file:0:8}_sort_bam & 
done
wait
echo -e "bam files sorted.\n"
echo -e "alignment finished.\n"

#4 count aligned reads
# count the number of reads that align to the regions of the genome that code for genes (assume all genes have no introns)
# use bedtools intersect to count number of aligned reads
for file in *sort_bam
do
bedtools bamtobed -i ${file} > ${file:0:13}.bed  &  #convert the bam files to bed format
done
wait
# for each feature in the “a” file, the number of overlapping features in the “b” file
for file in *sort_bam
do
bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b ${file:0:13}.bed -c > ${file:0:8}_aligned.bed &
# -c report the number of hits in b
done
# chromosome    chr_start       chr_end     gene_name     gene_description     number_of_aligned_reads
wait
echo -e "count aligned reads finished.\n"

#5 calculate the statistical mean of expression levels of each gene for each group
# generate files that contain sample names that are in the same group (same type, same time and same treatment)
for type in WT Clone1 Clone2
do
for treat in Uninduced Induced
do
for time in {0,24,48}
do
grep ${type} Tco.fqfiles | grep ${treat} | grep -w ${time} | sort -k3 > ${type}_${time}_${treat}.txt
done
done
done
echo -e "categories generated.\n"
# *aligned.bed files are already in order and also in same length, so sum up each line of files in the same group
for type in WT Clone1 Clone2
do
for treat in Uninduced Induced
do
for time in {0,24,48}
do
# make expression table in that category
count=0
cut -f 1 ${type}_${time}_${treat}.txt > category_file # contain file_names in that category
cut -f 4,5 Tco-5053_aligned.bed > ${type}_${time}_${treat}_expression.txt # columns of gene names, they are the same in all bed files
while read -r line  # go through each file_name in that category
do 
let count=count+1 # count how many files in this category (how many expression columns)
cut -f 6 Tco-${line:3:6}_aligned.bed > expression # the column of the expression of this file
paste ${type}_${time}_${treat}_expression.txt expression > pasting_file # paste gene name and expressions from each file together
cat pasting_file > ${type}_${time}_${treat}_expression.txt 
done < category_file
echo ${count} samples in ${type}_${time}_${treat}_expression.txt
# calculate average the expression of each gene in that category
awk -F "\t" 'BEGIN{sum=0;ave=0;}{
        if('$count'==3)
        {sum=$3+$4+$5;ave=sum/'$count';print ave}
        else {if('$count'==4){sum=$3+$4+$5+$6;ave=sum/'$count';print ave}
        else{print "error:one file have more than 4 expressions!!"}}
        }' ${type}_${time}_${treat}_expression.txt > expression_average
cut -f 4,5 Tco-5053_aligned.bed > head_gene
paste head_gene expression_average > ${type}_${time}_${treat}_expression_average.txt 
# table: gene_name gene_description expression_level
done
done
rm -f ${type}_0_Induced_expression.txt
rm -f ${type}_0_Induced_expression_average.txt 
done
echo -e "average calculation finished.\n"

#6 generate "fold change" data for the "group-wise" comparisons
# calculate the fold change of Uninduced and Induced groups that are in the same type and time category
for type in WT Clone1 Clone2
do
for time in {24,48}
do
    ./FoldChange.sh  ${type}_${time}_Uninduced_expression_average.txt ${type}_${time}_Induced_expression_average.txt
done
done
# move results into all_result_folder
mv *aligned.bed all_result_folder/.
mv *expression.txt all_result_folder/.
mv *expression_average.txt all_result_folder/.
mv *.FoldChange.txt all_result_folder/.
rm -f -r Tcongo_genome fastq
for file in *
do
if test ${file} == "script"  || test ${file} == "FoldChange.sh" # space are required between ==
then
echo ${file}
else
rm -f ${file}
fi
done
echo -e " To do other kinds of "group-wise" comparisons, simply type in:\n./FoldChange.sh <original_file> <subsequent_file>"
echo -e "File names that can be inputted are end with “_expression_average.txt”."
echo -e "For example:\n cd all_result_folder \n ./FoldChange.sh WT_0_Uninduced_expression_average.txt  WT_48_Uninduced_expression_average.txt"
cp FoldChange.sh all_result_folder/.
cp script all_result_folder/.
cd all_result_folder

#!/bin/bash
# for all *_expression_average.txt files, generate fold change by inputing an original file(X) and a subsequent file(Y), 
#using formular: foldechange = (Y-X)/X.
# ./FoldChange.sh X Y
        echo "Calculating Fold Change of " ${2//_expression_average.txt/} " with respect to " ${1//_expression_average.txt/} " using formular of (Y-X)/X "
        cut -f 1,2 $1 > head_gene
        cut -f 3 $1 > fileX
        cut -f 3 $2 > fileY
        paste fileX fileY > file         # only the two expression columns
        # caluculate fold change
        awk -F '\t' 'BEGIN{fold=0;}{
                if($1!=0){fold=($2-$1)/$1;print fold}
                else{print "Infinity"}
                }' file > foldchanges
        paste head_gene foldchanges > FoldChanges.txt
        echo -e "Gene_ID\tGene_Name\tFold_Change" > ${1//_expression_average.txt/}_to_${2//_expression_average.txt/}.FoldChange.txt
        # fold-change sorted in decresing order
        sort -t$'\t' -k3,3nr FoldChanges.txt >> ${1//_expression_average.txt/}_to_${2//_expression_average.txt/}.FoldChange.txt
        echo The result is in ${1//_expression_average.txt/}_to_${2//_expression_average.txt/}.FoldChange.txt
        rm -f fileX fileY file head_gene foldchange FoldChanges.txt
# Running Example: ./FoldChange.sh WT_0_Uninduced_expression_average.txt WT_48_Uninduced_expression_average.txt


