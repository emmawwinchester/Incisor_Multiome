#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4 
#SBATCH --mem-per-cpu=16G
#SBATCH -p general
#SBATCH --qos=general
#SBATCH --mail-user=cotney@uchc.edu
#SBATCH --mail-type=end
#SBATCH --job-name=hic2juicer
#SBATCH --export=ALL
#SBATCH -o %j.out
#SBATCH -e %j.err

#SORT YOUR DAMN CHROMSIZES FILES!!!!!!!!!!!!!!!!!!!!!!! THE BED AND THE TXT FILE OR IT WILL NOT RUN!!!!!!!!!!!!!!!!!!!!!!!!!!
#DO NOT RUN JUICEBOX IN THE CONDA ENVIRONMENT, IT WILL NOT RUN PROPERLY. 
module load juicer
export ANALYSISDIR=`pwd`
export SAMPLE=AWIL107
export JOB_NAME=${SAMPLE}_Juicer_convert
export BAMDIR=/home/FCAM/awilderman/ANALYSIS/HiC/Human/CS17_03-29-21/HiCPro/Results/bowtie_results/bwt2/AWIL107
export GENOME=hg19_nh #hg19, mm10, etc.
export JUICERTOOLS=/isg/shared/apps/juicer/1.8.9/SLURM/scripts/juicer_tools.jar #<path to juicer_tools.jar>  
export LC_COLLATE=C
export BEDGRAPHS=/home/FCAM/ewentworth/cotney/analysis/abcenhancer/human_cf/bedgraphs
export NARROWPEAKS=/home/FCAM/jcotney/ANALYSIS/ABC_Enhancer_Gene_Prediciton/enhancer_lists
export HIC=/home/FCAM/ewentworth/cotney/analysis/abcenhancer/human_cf/hg19/HiC-out/

export ABC=/home/FCAM/ewentworth/cotney/analysis/abcenhancer/tooth
export SEGMENTATION=/home/FCAM/ewentworth/cotney/rawdata/chromatinsegmentations/human_25state
export bigwigs=/linuxshare/users/jcotney/ANALYSIS/ChIP/ChromImpute/CHROMIMPUTE_APPLY
export TOOLS=/home/FCAM/jcotney/TOOLS/ABC-Enhancer-Gene-Prediction
export CHROMSIZES=/home/FCAM/ewentworth/cotney/genome/hg19/hg19.chrom.sizes
export CHROMSIZESBED=/home/FCAM/ewentworth/cotney/genome/hg19/hg19.chrom.sizes.bed
export TARAGENES=/home/FCAM/tyankee/ANALYSIS/RNA-seq/Cranio_facial/CS/Gencode-v25.bed
export GENES=/home/FCAM/ewentworth/cotney/genome/hg19/ArchR_gene_coordinates_hg19.tsv
export UBIQUITOUS=/home/FCAM/ewentworth/cotney/analysis/abcenhancer/human_cf/ubiquitous_genes_symbol.txt
#export HIC=/home/FCAM/awilderman/ANALYSIS/HiC/Human/CS17_03-29-21/HiCPro/Results/hic_results/data/AWIL107/AWIL107_allValidPairs
export HIC=$ABC/../human_cf/HiC-out/
export RNA=/home/FCAM/tyankee/ANALYSIS/RNA-seq/Cranio_facial/CS/CF_all_TPM.tsv

source /home/FCAM/jcotney/.bashrc_miniconda3 
conda activate abc_enhancer_gene
module load bedtools
module load kent-tools


#!/bin/bash
module load samtools/1.9
module load htslib/1.10.2
module load juicer/1.8.9
module load juicebox/1.9.8
module load kent-tools
conda activate archr
R
load('~/cotney/one2one_ortho.Rdata')
list<-list.files(pattern='_TPM.tsv')
list<-list[grep('hg38', list, invert=TRUE)]
for (filename in list){
tab<-read.table(filename)
tab<-tab[grep('Rik|LINC|Linc|Loc', tab$V1, ignore.case=TRUE, invert=TRUE),]
tab<-tab[which(tab$V1 %in% genesV2$MGI.symbol),]
standardmarkers<-genesV2[which(genesV2$MGI.symbol %in% tab$V1),]
tab2<-merge(tab, standardmarkers, by.x='V1', by.y='MGI.symbol')
tab2<-tab2[,c(5,2)]
write.table(tab2, file=paste('hg19-', filename, sep=''), sep='\t', quote=F, col.names=F, row.names=F)
}
list<-list.files(pattern='.bed')
list<-list[grep('Tooth|padded|unmapped|hg38|hg19', list, invert=TRUE)]
for (filename in list){
tab<-read.table(filename)
tab$V7<-gsub('\\..*', '', tab$V7)
tab$V8<-gsub('\\..*', '', tab$V8)
tab$V9<-gsub('\\..*', '', tab$V9)
tab$V10<-gsub('\\..*', '', tab$V10)
write.table(tab, file=filename, sep='\t', quote=F, col.names=F, row.names=F)
}
q()
n
conda activate abc_enhancer_gene

egrep -v -hi 'Rik|LINC|LOC' $GENES | awk '{print $1, $2-2000, $2+2000, $4, $5, $6}' | sed 's/ /\t/g' > padded_TSS.bed


for sample in $(ls *.bed | egrep -v 'hg38|enh|ABC|unmapped|temp')
do
export NAME=`echo $sample | sed 's/.bed//g'`
liftOver $sample -minMatch=0.25 ~/cotney/genome/mm10/mm10ToHg19.over.chain.gz hg19-$sample unmapped.bed
egrep -hi -v 'Un|random|M' hg19-$sample > temp.txt
bedClip temp.txt $CHROMSIZES hg19-$sample
done


for sample in $(ls *bw | egrep -v 'hg')
do
export NAME=`echo $sample | sed 's/.bw//g'`
echo -e "#/bin/bash\n#SBATCH --job-name=liftover\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -c 4\n#SBATCH --partition=general\n#SBATCH --qos=general\n#SBATCH --mail-type=END\n#SBATCH --mem=150G\n#SBATCH --mail-user=wentworth@uchc.edu\n#SBATCH -o %j.out\n#SBATCH -e %j.err\nmodule load bedtools\nmodule load kent-tools\nbigWigToBedGraph "$sample $NAME".bgr\nawk '{print \$0, \$3-\$2}' "$NAME".bgr | sed 's/ /\t/g' > "$NAME".temp.txt\nliftOver "$NAME".temp.txt -minMatch=0.25 ~/cotney/genome/mm10/mm10ToHg19.over.chain.gz hg19-"$NAME".bgr unmapped.bed\negrep -hi -v 'Un|random|chrm' hg19-"$NAME".bgr | awk '{if (\$5 <= 1.5*(\$3-\$2)){print \$1, \$2, \$3, \$4}}' > "$NAME".temp.txt\nbedClip "$NAME".temp.txt "$CHROMSIZES" "$NAME".temp.bgr\nbedtools sort -i "$NAME".temp.bgr | bedtools merge -d -1 -c 4 -o sum -i - > "$NAME".temp1.bbb\nbedGraphToBigWig "$NAME".temp1.bbb "$CHROMSIZES" hg19-"$NAME".bw" | sed 's/\/bin/!\/bin/g' > "$NAME".sh
sbatch "$NAME".sh
done

bedClip hg19_blacklist.txt $CHROMSIZES temp.txt
mv temp.txt hg19_blacklist

#######################make predictions
conda activate abc_enhancer_gene
for sample in $(ls *.bed | egrep 'hg19' | egrep -v 'ooth|un|txt' )
do
export NAME=`echo $sample | sed 's/.bed//g'`
##############MAKE SURE THE OUTPUT IS SORTED
#egrep -v 'Un|random|chrM' "$NAME".bgr > temp2.txt
#bedClip temp2.txt $CHROMSIZES temp.txt
#mv temp.txt "$NAME".bgr
#egrep -hi -v 'Un|random|chrm' "$NAME".bgr | awk '{if ($5 <= 1.5*($3-$2)){print $1, $2, $3, $4}}' > "$NAME".temp.txt
#bedClip "$NAME".temp.txt "$CHROMSIZES" "$NAME".temp.bgr
#bedtools sort -i "$NAME".temp.bgr | bedtools merge -d -1 -c 4 -o sum -i - > "$NAME".temp1.bbb
#bedGraphToBigWig "$NAME".temp1.bbb "$CHROMSIZES" "$NAME".bw
python $TOOLS/src/makeCandidateRegions.py --narrowPeak $sample --bam "$NAME".bw --outDir $ABC/candidateregions --chrom_sizes $CHROMSIZES --ignoreSummits --nStrongestPeaks 150000 --regions_blocklist hg19_blacklist.txt
python $TOOLS/src/run.neighborhoods.py --candidate_enhancer_regions candidateregions/"$NAME".bed.candidateRegions.bed --outdir neighborhoods --genes $GENES --expression_table "$NAME"_general_celltypes_TPM.tsv --chrom_sizes $CHROMSIZES --default_accessibility_feature ATAC --ATAC "$NAME".bw --ubiquitously_expressed_genes $UBIQUITOUS --cellType $NAME 
python $TOOLS/src/predict.py --enhancers neighborhoods/EnhancerList.txt --chrom_sizes $CHROMSIZES --genes neighborhoods/GeneList.txt --threshold 0.005 --make_all_putative --HiCdir $HIC --hic_resolution 5000 --scale_hic_using_powerlaw --cellType $NAME --outdir $ABC/DHS_Filtered_Predictions_wHiC/ --chromosomes chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chrX
cat $ABC/DHS_Filtered_Predictions_wHiC/EnhancerPredictions.bedpe | awk '{if ($2 < $6) {print $1"\t"$2"\t"$6"\t"$7"\t0\t"$8"\tface\t#7A67EE\t"$1"\t"$2"\t"$3"\t"$7"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$10} else {print $4"\t"$6"\t"$2"\t"$7"\t0\t"$8"\tface\t#7A67EE\t"$4"\t"$5"\t"$6"\t"$7"\t"$9"\t"$1"\t"$2"\t"$3"\t"$7"\t"$10}}' | egrep -v "track" | sort -k1,1 -k2,2n | awk '{if ($6 >= 0.01) print $0}' > $ABC/"$NAME"_faceEnhancerPredictions.05.interact
rm -r neighborhoods/ candidateregions/ DHS_Filtered_Predictions_wHiC/
done



conda activate archr
Rscript HighestInteractions.Rscript


for sample in *keep*
do
sed 's/|/_/g' $sample
done



####################make loadable on mouse genome
rm hg19_loadtracks.txt
for sample in hg19*05.interact
do
export NAME=`echo $sample | sed 's/_faceEnhancerPredictions.05.interact//g'`
sed 's/|/_/g' $sample | egrep -f "$NAME"-keep_interactions.txt - > temp.txt
mv temp.txt $sample
bedToBigBed -type=bed5+13 -as=interact.as -extraIndex=name -unc $sample $CHROMSIZES /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/"$NAME".05.bigInteract
echo -e "track type=bigInteract name=\""$NAME" ABC\" description=\""$NAME" ABC Loops of face Enhancers and Target Genes\" useScore=on visibility=full maxHeightPixels=50 bigDataUrl=http://graveleylab.cam.uchc.edu/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/"$NAME".05.bigInteract" >> hg19_loadtracks.txt
scp -r "$NAME".bw /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/.
echo -e "track type=bigWig name=\""$NAME" bigwig\" description=\""$NAME" bigwig\" visibility=dense bigDataUrl=http://graveleylab.cam.uchc.edu/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/"$NAME".bw" >> hg19_loadtracks.txt
chmod -R 755 /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth
done

####################make loadable on human genome
for sample in hg19*.bed hg19*enhancers.bed 
do
export NAME=`echo $sample | sed 's/.bed//g'`
bedtools sort -i $sample | cut -f1,2,3 | egrep -v -hi 'Un|random' | sed 's/ /\t/g' > temp1.txt
bedToBigBed -type=bed3 -unc temp1.txt $CHROMSIZES /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/$NAME.bigBed
echo -e "track type=bigBed name=\""$NAME"\" description=\""$NAME"\" useScore=on visibility=dense bigDataUrl=http://graveleylab.cam.uchc.edu/WebData/cotney/hubs/ChIP/ewentworth/mousetooth/"$NAME".bigBed" >> hg19_loadtracks.txt
chmod -R 755 /tgc/TGCore_User_Data/WebData/cotney/hubs/ChIP/ewentworth/mousetooth
done
#in R


