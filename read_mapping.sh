#! /bin/bash

################
# Requirements #
################

# TODO

###########
# Options #
###########

# default values
FASTQ_FILE_2=""
SORT_BOOL=0

while getopts "c:g:h:i:o:s:t" param
do
	case $param in
		# BAM compression level
		c) BAM_COMPRESSION_LEVEL=$OPTARG;;
		# Path to ungzipped fastq file 1
		g) FASTQ_FILE_1=$OPTARG;;
		# Path to ungzipped fastq file 2
		h) FASTQ_FILE_2=$OPTARG;;
		# Path to genome indexes directory
		i) GENOME_INDEXES_DIR=$OPTARG;;
		# Path to output filename prefix
		o) OUTPUT_FILENAME_PREFIX=$OPTARG;;
		# Sample ID
		s) SAMPLE_ID=$OPTARG;;
		# Sort BAM file 
		t) SORT_BOOL=1;;
		:) echo "Option -$OPTARG expects an argument. Exit." ; exit;;
		\?) echo "Unvalid option -$OPTARG. Exit." ; exit;;
	esac
done



function create_directory {
	DIR_TO_CREATE=$1
	if [ ! -d ${DIR_TO_CREATE} ]; then
		echo "create directory: ${DIR_TO_CREATE}"
		mkdir -p ${DIR_TO_CREATE}
	fi
}

function run_cmd {
	local cmd=$@
	
	echo "${cmd}"
	eval time ${cmd}
	local error_code=$?
	if [ ${error_code} -ne 0 ]; then
		local message=${cmd%%*}
		echo "Error while ${message} execution."
		exit ${error_code}
	fi
}

##############
# Parameters #
##############

# STAR parameters
NB_THREADS=1

if [ ! -z ${FASTQ_FILE_2} ]; then
	FASTQ_FILES_PARAM="${FASTQ_FILE_1} ${FASTQ_FILE_2}"
else
	FASTQ_FILES_PARAM="${FASTQ_FILE_1}"
fi

if [ ${SORT_BOOL} -eq 1 ]; then
	SORT_PARAM="SortedByCoordinate"
	SORT_EXTENSION=".sortedByCoord"
else
	SORT_PARAM="Unsorted"
	SORT_EXTENSION=""
fi

if [ ! -z ${BAM_COMPRESSION_LEVEL} ]; then
	BAM_COMPRESSION_OPTION=" --outBAMcompression ${BAM_COMPRESSION_LEVEL}"
else
	BAM_COMPRESSION_OPTION=""
fi

# STAR read mapping
CMD="STAR --runThreadN ${NB_THREADS} --genomeDir ${GENOME_INDEXES_DIR} --readFilesIn ${FASTQ_FILES_PARAM} --outFileNamePrefix ${OUTPUT_FILENAME_PREFIX} --outSAMtype BAM ${SORT_PARAM}${BAM_COMPRESSION_OPTION}"
run_cmd ${CMD}

# remove STAR tmp directory
STAR_TMP_DIR=${OUTPUT_FILENAME_PREFIX}_STARtmp
if [ -d ${STAR_TMP_DIR} ]; then
	echo "remove tmp directory: ${STAR_TMP_DIR}"
	CMD="rm -Rf ${STAR_TMP_DIR}"
	run_cmd ${CMD}
fi

# Read mapping statistics
PATH_TO_BAM_FILE=${OUTPUT_FILENAME_PREFIX}Aligned${SORT_EXTENSION}.out.bam
if [ -e ${PATH_TO_BAM_FILE} ]; then
	CMD="samtools flagstat ${PATH_TO_BAM_FILE} > ${PATH_TO_BAM_FILE%.bam}.flagstat.txt"
	run_cmd ${CMD}
	
	CMD="samtools stats ${PATH_TO_BAM_FILE} > ${PATH_TO_BAM_FILE%.bam}.stats.txt"
	run_cmd ${CMD}
else
	echo "ERROR: the BAM file ${PATH_TO_BAM_FILE} does not exist. Exit."
	exit 1
fi
