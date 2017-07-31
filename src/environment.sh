#!/bin/bash

export project_root=$(pwd)
export reference=${project_root}/reference
export reads=${project_root}/reads
export maps=${project_root}/maps
export calls=${project_root}/calls
export src=$project_root/src
export log=$project_root/log
export split=${project_root}/split
export merges=${project_root}/merges
export alignments=${project_root}/alignments

mkdir -p $reference $reads $maps $calls $log $merges $split $alignments

export picard="${src}/jdk1.8.0_144/bin/java -jar $src/picard.jar"
#export gatk="java -Xmx16G -Djava.io.tmpdir=/tmp -jar $src/GenomeAnalysisTK.jar"
export gatklaunch="$src/gatk-4.beta.3-SNAPSHOT/gatk-launch"

source $src/functions.sh
