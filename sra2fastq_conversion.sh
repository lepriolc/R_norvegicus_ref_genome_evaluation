#! /bin/sh

################
# Requirements #
################

# TODO

###########
# Options #
###########

while getopts "i:o:" param
do
	case $param in
		# Path to SRA accession list file (one SRA accession per line)
		i) SRA_ACC_LIST_FILE=$OPTARG;;
		# Path to output directory
		o) OUTPUT_DIR=$OPTARG;;
		:) echo "Option -$OPTARG expects an argument. Exit." ; exit;;
		\?) echo "Unvalid option -$OPTARG. Exit." ; exit;;
	esac
done


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


while read line;
do
	SRA_ACC=$line
	SRA_OUTPUT_DIR=${OUTPUT_DIR}/${SRA_ACC}
	SRA_FILE=${SRA_OUTPUT_DIR}/${SRA_ACC}.sra
	if [ -e ${SRA_FILE} ]; then
		# Fastq conversion
		echo "Fastq conversion of ${SRA_FILE}"
		CMD="fastq-dump -I --split-files --gzip -O ${SRA_OUTPUT_DIR} ${SRA_FILE}"
		run_cmd ${CMD}
		# remove SRA file
		CMD="rm -f ${SRA_FILE}"
		run_cmd ${CMD}
	else
		echo "SRA file does not exist: ${SRA_FILE}"
	fi
done < ${SRA_ACC_LIST_FILE}
