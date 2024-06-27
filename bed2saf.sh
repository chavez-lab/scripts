#!/bin/bash

# Usage: bed2saf.sh bed.bed
awk 'OFS="\t" {print $1":"$2+1"-"$3, $1, $2+1, $3, "."}' $1
