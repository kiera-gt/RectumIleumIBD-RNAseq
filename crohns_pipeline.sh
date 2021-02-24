#!/bin/bash
scriptname=$0   # $0 is the name of the program
readLength=151
kmersize=35
trim=no
batch=""
vars=~/data/scripts/analysis_vars_ibd.sh
###
#help function
HELP () {
	echo -e "\n\tUSAGE: $scriptname -b BATCH [-t (yes|no)] [-v VARIABLES_FILE]\n"
	echo -e "\tMANDATORY OPTIONS:\n"
	echo -e "\t\t-b BATCH\t\tName of the sequencing run or project"
	echo -e "\tADDITIONAL OPTIONS:\n"
	echo -e "\t\t-v VARIABLES_FILE\tLocation of file containing variable information for analysis scripts"
	echo -e "\t\t\t\t\tDEFAULT: ~/data/scripts/analysis_vars_ibd.sh\n"						
	echo -e "\t\t-t TRIM\t\t\tWhether the FASTQ files to be aligned have been trimmed"
	echo -e "\t\t\t\t\tDEFAULT: no\n"
	echo -e "\t\t-l READ_LENGTH\t\tRead Length. Must be an integer."
	echo -e "\t\t\t\t\tDEFAULT: 151\n"
	exit 0
}
###
#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  HELP
fi
###
#Get Options from command line
while getopts :t:l:b:v:h opt
do
	case $opt in
	t)trim=$OPTARG;;
	b)batch=$OPTARG;;
	l)readLength=$OPTARG;;
	v)vars=$OPTARG;;
    h) HELP; exit;;
    \?) echo "ERROR: Invalid option: -$OPTARG" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
    :) echo "ERROR: Option -$OPTARG requires an argument" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
    *) echo "ERROR: Unimplemented option: -$OPTARG" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
  esac
done
shift $(($OPTIND -1))
###
#check to make sure variables file exists
if [ ! -f $vars ] 
then
    echo -e "ERROR: File $vars DOES NOT exist" >&2
	echo -e "Edit this scripts default or use the -v option to point the script to your file that contains variables for analysis scripts" >&2
    exit 1
fi
###
source $vars
###
#check whether or not to trim, set fastq_dir appropriately
fastq_type=""
if [ $trim == 'yes' ]
then 
	fastq_type=fastq_files/trimmed
	echo "Will align TRIMMED fastq files"
elif [ $trim == 'no' ]
then
	fastq_type=fastq_files
	echo "Will align UNTRIMMED fastq files"
else
	echo -e "ERROR: \"$trim\" is not a valid argument for [-t] option. Valid arguments are \"yes\" or \"no\""
	exit 1
fi
###
#ensure directory containing samples folder and samples.txt exists
if [ ! -d $basedir/$batch ] 
then
    echo -e "ERROR: Directory $basedir/$batch DOES NOT exist" >&2
	echo -e "Check that option arguments -b is correct" >&2
    exit 1
fi
###
#ensure samples.txt exists
if [ ! -f $basedir/$batch/samples.txt ] 
then
    echo -e "ERROR: File $basedir/$batch/samples.txt DOES NOT exist" >&2
	echo -e "Create the standard samples.txt file for this batch" >&2
    exit 1
fi
###
#set batchdir for ease of use later in script
batchdir=$basedir/$batch
###
#check if output folders already exist, create if not
if [ ! -d $batchdir/bamfiles ]
then
	mkdir $batchdir/bamfiles
fi
###
if [ ! -d $batchdir/star ]
then
	mkdir $batchdir/star
fi
###
if [ ! -d $batchdir/starlogs ]
then
	mkdir $batchdir/starlogs
fi
###
if [ ! -d $batchdir/starsj ]
then
	mkdir $batchdir/starsj
fi
###
if [ ! -d $batchdir/psi ]
then
	mkdir $batchdir/psi
fi
###
if [ ! -d $batchdir/genecounts ]
then
	mkdir $batchdir/genecounts
fi
###
if [ ! -d $batchdir/qorts ]
then
	mkdir $batchdir/qorts
fi
###
if [ ! -d $batchdir/variants ]
then
	mkdir $batchdir/variants
fi
###
if [ ! -d $batchdir/variants/temp ]
then
	mkdir $batchdir/variants/temp
fi
###
if [ ! -d $batchdir/variants/gvcfs ]
then
	mkdir $batchdir/variants/gvcfs
fi
###

###
#star 2pass function
star_2pass () {
	local sample=$1
	local dir=$2
	local readgroup=$3
	local tissue=$4
	if [ ! -d $dir/star/$sample ]
	then
		mkdir $dir/star/$sample #STAR needs outFileNamePrefix to already exist before running
	fi
	local step=star
	local pbsfile=$dir/star/$sample.$step.pbs
	read1=""
	read2=""
	counter=1
	rgline=""
	echo -e "#PBS -N $sample.$step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for star run

	for file in $dir/$fastq_type/$sample*_R1_001.fastq.gz
	do
		if [ "$counter" -eq 1 ]
		then
			read1="${file}"
			counter=$((counter+1))
		else
			read1="${read1},${file}"
		fi
			#filesub=`echo -e "$(basename $file)" | sed 's/\(_R1_001.fastq.gz\)//g'`
	done
	counter=1
	for file in $dir/$fastq_type/$sample*_R2_001.fastq.gz
	do
		if [ "$counter" -eq 1 ]
		then
			read2="${file}"
			lane=`echo "${file: -17}" | sed 's/_R2_001.fastq.gz//g'`
			rgline="ID:${readgroup}.${lane} PL:illumina PU:${readgroup}.${lane} LB:${tissue} SM:${sample}"
			counter=$((counter+1))
		else
			read2="${read2},${file}"
			lane=`echo "${file: -17}" | sed 's/_R2_001.fastq.gz//g'`
			rgline="${rgline} , ID:${readgroup}.${lane} PL:illumina PU:${readgroup}.${lane} LB:${tissue} SM:${sample}"
			counter=$((counter+1))
		fi
	done
	#counter=$((counter-1))
	#if [ "$counter" -ne 1 ]
	#then
	#	for i in $(seq 2 $counter)
	#	do
	#		rgline="${rgline} , ID:${readgroup}.${i} PL:illumina PU:${readgroup}.${i} LB:${seqtype} SM:${sample}"
	#	done
	#fi
	echo "$STAR --runThreadN 16 --genomeDir $genomeDir --sjdbGTFfile $annotations --readFilesIn $read1 $read2 --readFilesCommand zcat --outFileNamePrefix $dir/star/$sample/ --sjdbOverhang 100 --outSAMtype BAM SortedByCoordinate --outSAMattrRGline $rgline --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMultimapNmax 1 --outSJfilterCountUniqueMin 5 5 5 5 --outSJfilterCountTotalMin 5 5 5 5 --outFilterType BySJout --outSJfilterReads Unique --alignInsertionFlush Right --quantMode GeneCounts --twopassMode Basic" >> $pbsfile #add commands for star run
	echo -e "mv $dir/star/$sample/Aligned.sortedByCoord.out.bam $dir/bamfiles/$sample.bam" >> $pbsfile
	echo -e "mv $dir/star/$sample/Log.final.out $dir/starlogs/$sample.log.final.out" >> $pbsfile
	echo -e "mv $dir/star/$sample/SJ.out.tab $dir/starsj/$sample.sj.out.tab" >> $pbsfile
	echo -e "mv $dir/star/$sample/ReadsPerGene.out.tab $dir/genecounts/$sample.genecounts.out.tab" >> $pbsfile
	echo -e "samtools index $dir/bamfiles/$sample.bam" >> $pbsfile
	echo -e "cd $dir/psi/" >> $pbsfile
	echo -e "qsub $dir/psi/$sample.psi.pbs" >> $pbsfile
}
###
#psi function
psi () {
	local sample=$1
	local dir=$2
	local step=psi
	local pbsfile=$dir/psi/$sample.$step.pbs
	echo -e "#PBS -N $sample.$step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "$psi_script -s $sample -b $dir -l $readLength -v $vars" >> $pbsfile
	#echo -e "sed -i \"1i #\t#\t$sample\t$sample\t$sample\" $dir/psi/$sample.exonic_parts.psi" >> $pbsfile
}
###
qorts () {
	local sample=$1
	local dir=$2
	local readgroup=$3
	local step=qorts
	local pbsfile=$dir/qorts/$sample.$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step

	for file in $dir/$fastq_type/$sample*_R1_001.fastq.gz
	do
		read="${file}"
		lane=`echo "${file: -17}" | sed 's/_R1_001.fastq.gz//g'`
		if [ ! -d $dir/qorts/$sample.$readgroup.$lane ]
		then
			mkdir $dir/qorts/$sample.$readgroup.$lane
		fi
		echo -e "java -jar $qorts QC --verbose --stranded --readGroup ${readgroup}.${lane} $dir/bamfiles/$sample.bam $annotations $dir/qorts/${sample}.${readgroup}.${lane}/" >> $pbsfile
	done
}
###
mark_duplicates () {
	local sample=$1
	local dir=$2
	local step=markdups
	local pbsfile=$dir/variants/$sample.$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "java -jar $picard MarkDuplicates I=$dir/bamfiles/$sample.bam O=$dir/variants/temp/$sample.dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=true M=$dir/variants/$sample.markdups_metrics.txt" >> $pbsfile
	echo -e "cd $dir/variants/" >> $pbsfile #move into variant_calling dir so pbs outfile goes to here
	echo -e "qsub $dir/variants/$sample.splitn.pbs\n" >> $pbsfile
}
###
split_reads () {
	local sample=$1
	local dir=$2
	local step=splitn
	local pbsfile=$dir/variants/$sample.$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "java -jar $gatk -T SplitNCigarReads -R $wholeGenomeFasta -I $dir/variants/temp/$sample.dedup.bam -o $dir/variants/temp/$sample.dedup.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS" >> $pbsfile
	echo -e "cd $dir/variants/" >> $pbsfile #move into variant_calling dir so pbs outfile goes to here
	echo -e "qsub $dir/variants/$sample.gvcf.pbs\n" >> $pbsfile
}
###
gvcf () {
	local sample=$1
	local dir=$2
	local step=gvcf
	local pbsfile=$dir/variants/$sample.$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "java -jar $gatk -T HaplotypeCaller -R $wholeGenomeFasta -I $dir/variants/temp/$sample.dedup.split.bam --dbsnp $dbsnp -kmerSize $kmersize -dontUseSoftClippedBases -stand_call_conf 20.0 --emitRefConfidence GVCF -L $intervals -o $dir/variants/gvcfs/$sample.plink.snps.indels.g.vcf" >> $pbsfile
	echo -e "$dir/variants/gvcfs/$sample.plink.snps.indels.g.vcf" >> $dir/variants/jointcall/plink.gvcfs.list
}
###
counter=1
while read sampleid fastqname readgrp tissuetype
do
	if [ $sampleid == '11469AR' ]
	then
		counter=2
	fi
	if [ "$counter" -eq 2 ]
	then

		#star_2pass $sampleid $batchdir $readgrp $tissuetype
		#psi $sampleid $batchdir 
		#qorts $sampleid $batchdir $readgrp
		mark_duplicates $sampleid $batchdir
		split_reads $sampleid $batchdir
		gvcf $sampleid $batchdir
		###
		#cd $batchdir/star/
		#qsub $batchdir/star/$sampleid.star.pbs
		cd $batchdir/variants/
		qsub $batchdir/variants/$sampleid.markdups.pbs
	else
		continue
	fi
done <$batchdir/samples.txt
###
###