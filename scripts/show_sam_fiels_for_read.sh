#!/bin/bash

# Check if QNAME argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <BAM_FILE> <QNAME>"
    exit 1
fi
if [ -z "$2" ]; then
    echo "Usage: $0 <BAM_FILE> <QNAME>"
    exit 1
fi
# Execute the command with QNAME argument
samtools view "$1" | grep "$2" | awk '{print "\033[34mQname:\033[0m", $1, "\n\033[34mflag:\033[0m", $2, "\n\033[34mrname:\033[0m", $3, "\n\033[34mpos:\033[0m", $4, "\n\033[34mmapq:\033[0m", $5, "\n\033[34mcigar:\033[0m", $6, "\n\033[34mrnext:\033[0m", $7, "\n\033[34mpnext:\033[0m", $8, "\n\033[34mtlen:\033[0m", $9, "\n\033[34mseq:\033[0m", $10, "\n\033[34mqual:\033[0m", $11}'
