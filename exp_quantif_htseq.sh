#! /bin/bash

################
# Requirements #
################

# TODO

###########
# Options #
###########

# default values
ANNOTATION="NCBI_RefSeq"
REMOVE_BAM=0
HTSEQ_COUNT_MODE="intersection-nonempty"
OUTPUT_ANNOTATED_BAM_FILE_BOOL=0
PAIRED_END=0
SORT_BOOL=0

FORMAT="bam"

# get options
while getopts "a:b:cg:m:o:prst:" param
do
	case $param in
		# Annotation
		a) ANNOTATION=$OPTARG;;
		# Path to BAM file
		b) BAM_FILE=$OPTARG;;
		# Cleaning: remove BAM files
		c) REMOVE_BAM=1;;
		# Path to genome annotation GTF file
		g) GTF_FILE=$OPTARG;;
		# htseq-cound mode to handle reads overlapping more than one feature
		m) HTSEQ_COUNT_MODE=$OPTARG;;
		# Path to output directory for gene quantification
		o) OUTPUT_DIR=$OPTARG;;
		# Output an annotated BAM file
		p) OUTPUT_ANNOTATED_BAM_FILE_BOOL=1;;
		# Paired-end reads
		r) PAIRED_END=1;;
		# Sort BAM file 
		s) SORT_BOOL=1;;
		# Strandedness during cDNA synthesis
		t) STRANDED=$OPTARG;;
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


# create output directory if needed
create_directory ${OUTPUT_DIR}
# set tmp directory variable
TMP_DIR=${OUTPUT_DIR}/tmp_exp_quantif_htseq


# sort BAM file if needed
if [ ${SORT_BOOL} -eq 1 ]; then
	FILENAME=$(basename ${BAM_FILE})
	BAM_FILE_TO_USE=${OUTPUT_DIR}/${FILENAME%.bam}.sort.bam
	CMD="samtools sort ${BAM_FILE} -O BAM -o ${BAM_FILE_TO_USE}"
	run_cmd ${CMD}
else
	BAM_FILE_TO_USE=${BAM_FILE}
fi

# gunzip GTF file if needed 
if [[ ${GTF_FILE} == *.gz ]]; then
	# create tmp directory if needed
	create_directory ${TMP_DIR}

	# gunzip
	echo "gunzip GTF file"
	FILENAME=$(basename ${GTF_FILE})
	if [[ ${ANNOTATION} == "Ensembl" ]]; then
		echo "gunzip -c ${GTF_FILE} | gunzip -c > ${TMP_DIR}/${FILENAME%.*}"
		gunzip -c ${GTF_FILE} | gunzip -c > ${TMP_DIR}/${FILENAME%.*}
	else
		echo "gunzip -c ${GTF_FILE} > ${TMP_DIR}/${FILENAME%.*}"
		gunzip -c ${GTF_FILE} > ${TMP_DIR}/${FILENAME%.*}
	fi
	GTF_FILE=${TMP_DIR}/${FILENAME%.*}
fi

# htseq-count
FILENAME=$(basename ${BAM_FILE_TO_USE})
CMD="samtools index ${BAM_FILE_TO_USE}"
run_cmd ${CMD}

# set output basename
OUTPUT_BASENAME=${FILENAME%.bam}

if [ ${OUTPUT_ANNOTATED_BAM_FILE_BOOL} -eq 1 ]; then
	O_OPTION="-o ${OUTPUT_DIR}/${OUTPUT_BASENAME}.annotated.bam "
	P_OPTION="-p ${FORMAT} "
else
	O_OPTION=""
	P_OPTION=""
fi

if [ ${PAIRED_END} -eq 1 ]; then
	R_OPTION="-r pos "
else
	R_OPTION=""
fi

if [ ! -z ${STRANDED} ]; then
	S_OPTION="-s ${STRANDED} "
else
	S_OPTION=""
fi


CMD="htseq-count ${BAM_FILE_TO_USE} ${GTF_FILE} -f ${FORMAT} -m ${HTSEQ_COUNT_MODE} ${O_OPTION}${P_OPTION}${R_OPTION}${S_OPTION}> ${OUTPUT_DIR}/${OUTPUT_BASENAME}.counts.txt"
run_cmd ${CMD}

# remove input sorted BAM file
if [ ${REMOVE_BAM} -eq 1 ]; then
	CMD="rm -f ${BAM_FILE_TO_USE}"
	run_cmd ${CMD}
fi

# annotated BAM file
if [ ${OUTPUT_ANNOTATED_BAM_FILE_BOOL} -eq 1 ]; then
	ANNOTATED_BAM_BASENAME=${OUTPUT_BASENAME}.annotated
	ANNOTATED_BAM_FILE=${OUTPUT_DIR}/${ANNOTATED_BAM_BASENAME}.bam
	if [ -e ${ANNOTATED_BAM_FILE} ]; then
		# only keep additional tag in the annotated output BAM file
		if [ ${PAIRED_END} -eq 1 ]; then
			CMD="samtools view ${ANNOTATED_BAM_FILE} | cut -f 1,2,3,4,5,7,8,16 | gzip -c > ${OUTPUT_DIR}/${ANNOTATED_BAM_BASENAME}.tsv.gz"
		else
			CMD="samtools view ${ANNOTATED_BAM_FILE} | cut -f 1,2,3,4,5,16 | gzip -c > ${OUTPUT_DIR}/${ANNOTATED_BAM_BASENAME}.tsv.gz"
		fi
		run_cmd ${CMD}
		# remove annotated BAM file
		if [ ${REMOVE_BAM} -eq 1 ]; then
			CMD="rm -f ${ANNOTATED_BAM_FILE}"
			run_cmd ${CMD}
		fi
	else
		echo "WARNING: the annotated BAM file ${ANNOTATED_BAM_FILE} does not exist."
	fi
fi


# remove tmp directory if needed
if [ -d ${TMP_DIR} ]; then
	echo "remove tmp directory: ${TMP_DIR}"
	rm -Rf ${TMP_DIR}
fi

