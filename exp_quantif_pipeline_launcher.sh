#! /bin/bash

################
# Requirements #
################

# TODO

###########
# Options #
###########

# default values
OUTPUT_ANNOTATED_BAM_FILE_BOOL=0
DRY_RUN=0
REMOVE_BAM=0
PAIRED_END=0
SORT_BOOL=0

while getopts "a:bc:def:g:i:l:m:o:pr:s:tu:" param
do
	case $param in
		# Annotation
		a) ANNOTATION=$OPTARG;;
		# Output an annotated BAM file
		b) OUTPUT_ANNOTATED_BAM_FILE_BOOL=1;;
		# BAM compression level
		c) BAM_COMPRESSION_LEVEL=$OPTARG;;
		# Dry run
		d) DRY_RUN=1;;
		# Cleaning: remove BAM files
		e) REMOVE_BAM=1;;
		# Path to genome annotation GTF file
		f) PATH_TO_GTF_FILE=$OPTARG;;
		# Path to genome indexes directory
		g) GENOME_INDEXES_DIR=$OPTARG;;
		# Path to file storing SRA sample IDs, one ID per line
		i) SRA_ACC_LIST_FILE=$OPTARG;;
		# Minimum read length
		l) MIN_LENGTH=$OPTARG;;
		# htseq-cound mode to handle reads overlapping more than one feature
		m) HTSEQ_COUNT_MODE=$OPTARG;;
		# Path to output directory for read mapping
		o) OUTPUT_DIR=$OPTARG;;
		# Paired-end reads
		p) PAIRED_END=1;;
		# Path to raw data directory: stores fastq files in subdirectories, one subdirectory per sample
		r) RAW_DATA_DIR=$OPTARG;;
		# Path to src directory
		s) SRC_DIR=$OPTARG;;
		# Sort BAM file 
		t) SORT_BOOL=1;;
		# Strandedness during cDNA synthesis
		u) STRANDED=$OPTARG;;
		:) echo "Option -$OPTARG expects an argument. Exit." ; exit;;
		\?) echo "Unvalid option -$OPTARG. Exit." ; exit;;
	esac
done



function time_cmd {
	local cmd=$@
	
	echo "${cmd}"
	eval time ${cmd}
}

if [ ! -z ${ANNOTATION} ]; then
	A_OPTION="-a ${ANNOTATION} "
else
	A_OPTION=""
fi

if [ ${OUTPUT_ANNOTATED_BAM_FILE_BOOL} -eq 1 ]; then
	B_OPTION="-b "
else
	B_OPTION=""
fi

if [ ${REMOVE_BAM} -eq 1 ]; then
	D_OPTION="-d "
else
	D_OPTION=""
fi

if [ ! -z ${BAM_COMPRESSION_LEVEL} ]; then
	C_OPTION="-c ${BAM_COMPRESSION_LEVEL} "
else
	C_OPTION=""
fi

if [ ! -z ${MIN_LENGTH} ]; then
	L_OPTION="-l ${MIN_LENGTH} "
else
	L_OPTION=""
fi

if [ ! -z ${HTSEQ_COUNT_MODE} ]; then
	M_OPTION="-m ${HTSEQ_COUNT_MODE} "
else
	M_OPTION=""
fi

if [ ${PAIRED_END} -eq 1 ]; then
	P_OPTION="-p "
else
	P_OPTION=""
fi

if [ ${SORT_BOOL} -eq 1 ]; then
	T_OPTION=" -t"
else
	T_OPTION=""
fi

if [ ! -z ${STRANDED} ]; then
	U_OPTION=" -u ${STRANDED}"
else
	U_OPTION=""
fi


GENOME_INDEXES_NAME=$(basename ${GENOME_INDEXES_DIR})
while IFS= read -r line;
do
	SAMPLE_ID=$line
	SAMPLE_OUTPUT_DIR=${OUTPUT_DIR}/${SAMPLE_ID}/${GENOME_INDEXES_NAME}
	CMD="${SRC_DIR}/exp_quantif_pipeline.sh ${A_OPTION}${B_OPTION}${C_OPTION}${D_OPTION}-f ${PATH_TO_GTF_FILE} -g ${GENOME_INDEXES_DIR} -i ${SAMPLE_ID} ${L_OPTION}${M_OPTION}-o ${SAMPLE_OUTPUT_DIR} ${P_OPTION}-r ${RAW_DATA_DIR} -s ${SRC_DIR}${T_OPTION}${U_OPTION}"
	echo ${CMD}
	if [ ${DRY_RUN} -eq 0 ]; then
		eval time ${CMD}
		# if error during exp_quantif_pipeline.sh, then move on to the next sample
		if [ $? -ne 0 ]; then
			message=${CMD%%*}
			echo "Error while ${message} execution."
			continue
		fi
	fi
done < ${SRA_ACC_LIST_FILE}

