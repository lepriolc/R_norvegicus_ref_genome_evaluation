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
OUTPUT_ANNOTATED_BAM_FILE_BOOL=0
REMOVE_BAM=0
HTSEQ_COUNT_MODE="intersection-nonempty"
PAIRED_END=0
SORT_BOOL=0

while getopts "a:bc:df:g:i:l:m:o:pr:s:tu:" param
do
	case $param in
		# Annotation
		a) ANNOTATION=$OPTARG;;
		# Output an annotated BAM file
		b) OUTPUT_ANNOTATED_BAM_FILE_BOOL=1;;
		# BAM compression level
		c) BAM_COMPRESSION_LEVEL=$OPTARG;;
		# Cleaning: remove BAM files
		d) REMOVE_BAM=1;;
		# Path to genome annotation GTF file
		f) PATH_TO_GTF_FILE=$OPTARG;;
		# Path to genome indexes directory
		g) GENOME_INDEXES_DIR=$OPTARG;;
		# Sample ID
		i) SAMPLE_ID=$OPTARG;;
		# Minimum read length
		l) MIN_LENGTH=$OPTARG;;
		# htseq-cound mode to handle reads overlapping more than one feature
		m) HTSEQ_COUNT_MODE=$OPTARG;;
		# Path to output directory for expression quantification pipeline
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

function ungz {
	local PATH_FILE_TO_UNGZ=$1
	local PATH_TO_TMP_FILE=$2
	
	TMP_DIR=$(dirname ${PATH_TO_TMP_FILE})
	if [ -e ${PATH_FILE_TO_UNGZ} ]; then
		if (expr match "${PATH_FILE_TO_UNGZ}" '.*\(.gz$\)'); then
			# create tmp directory if needed
			create_directory ${TMP_DIR}			
			# gunzip
			CMD="gunzip -c ${PATH_FILE_TO_UNGZ} > ${PATH_TO_TMP_FILE}"
			run_cmd ${CMD}
		else
			echo "ERROR: the file ${PATH_FILE_TO_UNGZ} does not have the '.gz' extension. Exit."
			exit 1
		fi
	else
		echo "ERROR: the file ${PATH_FILE_TO_UNGZ} does not exist. Exit."
		exit 1
	fi
}

function check_cmd_exit_code {
	local cmd=$@
	
	echo "${cmd}"
	eval ${cmd}
	local exit_code=$?
	if [ ${exit_code} -ne 0 ]; then
		local message=${cmd%%*}
		echo "Error while ${message} execution, exit code: ${exit_code}"
		exit ${exit_code}
	fi
}

function check_file_exists {
	local PATH_TO_FILE=$1
	
	if [ ! -e ${PATH_TO_FILE} ]; then
		echo "WARNING: the file ${PATH_TO_FILE} does not exist. Return."
		exit 0 # exit code 0 to exit this particular script and eventually move on to the next iteration of a calling script
	fi
}


#################
# Set variables #
#################

# create output directory if needed
create_directory ${OUTPUT_DIR}

# Path to Fastq files
FASTQ_FILE_1=${RAW_DATA_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_1.fastq.gz
check_file_exists ${FASTQ_FILE_1}
FASTQ_FILES_PARAM=${FASTQ_FILE_1}
if [ ${PAIRED_END} -eq 1 ]; then
	FASTQ_FILE_2=${RAW_DATA_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_2.fastq.gz
	check_file_exists ${FASTQ_FILE_2}
	FASTQ_FILES_PARAM="${FASTQ_FILE_1} ${FASTQ_FILE_2}"
fi

# Path to tmp directory
TMP_DIR=${OUTPUT_DIR}/tmp_exp_quantif_pipeline


############
# Trimming #
############
if [ ! -z ${MIN_LENGTH} ]; then
	LENGTH_OPTION="--length ${MIN_LENGTH} "
else
	LENGTH_OPTION=""
fi

if [ ${PAIRED_END} -eq 1 ]; then
	PAIRED_OPTION="--paired "
else
	PAIRED_OPTION=""
fi

TRIM_OUTPUT_DIR=${OUTPUT_DIR}/trim_galore
create_directory ${TRIM_OUTPUT_DIR}
CMD="trim_galore --fastqc ${LENGTH_OPTION}-o ${TRIM_OUTPUT_DIR} ${PAIRED_OPTION}${FASTQ_FILES_PARAM}"
run_cmd ${CMD}

if [ ${PAIRED_END} -eq 1 ]; then
	FILENAME=$(basename ${FASTQ_FILE_1})
	FASTQ_FILE_1_TRIM=${OUTPUT_DIR}/trim_galore/${FILENAME%%.*}_val_1.fq.gz
	FILENAME=$(basename ${FASTQ_FILE_2})
	FASTQ_FILE_2_TRIM=${OUTPUT_DIR}/trim_galore/${FILENAME%%.*}_val_2.fq.gz
else
	FILENAME=$(basename ${FASTQ_FILE_1})
	FASTQ_FILE_1_TRIM=${OUTPUT_DIR}/trim_galore/${FILENAME%%.*}_trimmed.fq.gz
fi
	
######################
# gunzip Fastq files #
######################
# gunzip Fastq files if needed
check_file_exists ${FASTQ_FILE_1_TRIM}
echo "gunzip Fastq file 1"
FILENAME=$(basename ${FASTQ_FILE_1_TRIM})
PATH_TO_UNGZ_FILE=${TMP_DIR}/${FILENAME%.*}
ungz ${FASTQ_FILE_1_TRIM} ${PATH_TO_UNGZ_FILE}
FASTQ_FILE_1=${PATH_TO_UNGZ_FILE}

if [ ${PAIRED_END} -eq 1 ]; then
	check_file_exists ${FASTQ_FILE_2_TRIM}
	echo "gunzip Fastq file 2"
	FILENAME=$(basename ${FASTQ_FILE_2_TRIM})
	PATH_TO_UNGZ_FILE=${TMP_DIR}/${FILENAME%.*}
	ungz ${FASTQ_FILE_2_TRIM} ${PATH_TO_UNGZ_FILE}
	FASTQ_FILE_2=${PATH_TO_UNGZ_FILE}
fi

################
# Read mapping #
################
if [ ! -z ${BAM_COMPRESSION_LEVEL} ]; then
	C_OPTION="-c ${BAM_COMPRESSION_LEVEL} "
else
	C_OPTION=""
fi

if [ ${PAIRED_END} -eq 1 ]; then
	H_OPTION="-h ${FASTQ_FILE_2} "
else
	H_OPTION=""
fi

if [ ${SORT_BOOL} -eq 1 ]; then
	T_OPTION=" -t"
	SORT_EXTENSION=".sortedByCoord"
else
	T_OPTION=""
	SORT_EXTENSION=""
fi

GENOME_INDEXES_NAME=$(basename ${GENOME_INDEXES_DIR})
OUTPUT_FILENAME_PREFIX=${OUTPUT_DIR}/${SAMPLE_ID}_${GENOME_INDEXES_NAME}_

CMD="${SRC_DIR}/read_mapping.sh ${C_OPTION}-g ${FASTQ_FILE_1} ${H_OPTION}-i ${GENOME_INDEXES_DIR} -o ${OUTPUT_FILENAME_PREFIX} -s ${SAMPLE_ID}${T_OPTION}"
run_cmd ${CMD}

# remove intermediary files
## tmp directory
if [ -d ${TMP_DIR} ]; then
	echo "remove tmp directory: ${TMP_DIR}"
	rm -Rf ${TMP_DIR}
fi
## trimmed fastq files
echo "remove trimmed fastq file 1: ${FASTQ_FILE_1_TRIM}"
CMD="rm -f ${FASTQ_FILE_1_TRIM}"
run_cmd ${CMD}
if [ ${PAIRED_END} -eq 1 ]; then
	echo "remove trimmed fastq file 2: ${FASTQ_FILE_2_TRIM}"
	CMD="rm -f ${FASTQ_FILE_2_TRIM}"
	run_cmd ${CMD}
fi

PATH_TO_BAM_FILE=${OUTPUT_FILENAME_PREFIX}Aligned${SORT_EXTENSION}.out.bam

##############
# Read count #
##############
check_file_exists ${PATH_TO_BAM_FILE}

if [ ${REMOVE_BAM} -eq 1 ]; then
	C_OPTION="-c "
else
	C_OPTION=""
fi

if [ ${OUTPUT_ANNOTATED_BAM_FILE_BOOL} -eq 1 ]; then
	P_OPTION=" -p"
else
	P_OPTION=""
fi

if [ ${PAIRED_END} -eq 1 ]; then
	R_OPTION=" -r"
else
	R_OPTION=""
fi

if [ ${SORT_BOOL} -eq 1 ]; then
	# the BAM file has already been sorted during the read mapping step
	S_OPTION=""
else
	# the BAM file must be sorted by coordinates to run htseq-count
	S_OPTION=" -s"
fi

if [ ! -z ${STRANDED} ]; then
	T_OPTION=" -t ${STRANDED}"
else
	T_OPTION=""
fi

CMD="${SRC_DIR}/exp_quantif_htseq.sh -a ${ANNOTATION} -b ${PATH_TO_BAM_FILE} ${C_OPTION}-g ${PATH_TO_GTF_FILE} -m ${HTSEQ_COUNT_MODE} -o ${OUTPUT_DIR}${P_OPTION}${R_OPTION}${S_OPTION}${T_OPTION}"
run_cmd ${CMD}	

