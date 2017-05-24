#!/bin/bash

export project_root=$(pwd)
export reference=${project_root}/reference
export reads=${project_root}/reads
export maps=${project_root}/maps
export calls=${project_root}/calls
export src=$project_root/src
export logs=$project_root/log
export split=${project_root}/split
export merges=${project_root}/merges
export alignments=${project_root}/alignments

mkdir -p $reads $maps $calls $logs $merges $split

export picard="java -jar $src/picard.jar"
export gatk="java -jar $src/GenomeAnalysisTK.jar"

source $src/functions.sh
