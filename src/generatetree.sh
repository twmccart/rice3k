#!/bin/bash
set -euo pipefail

fasta=$1
cultivar_to_names="s/IRIS_313-8326/JC1/g
s/IRIS_313-11737/CHUNDI/g
s/IRIS_313-10603/UCP122/g
s/IRIS_313-11802/JIECAOZHAN/g
s/CX357/PSBRC80/g
s/IRIS_313-11643/NCS603B/g
s/IRIS_313-11723/DISSIGBE/g
s/IRIS_313-11794/TSIPALA/g
s/IRIS_313-8493/MELABOH/g
s/IRIS_313-11812/ASFALA/g
s/IRIS_313-10177/DAGANGZHAN/g
s/IRIS_313-11671/DERAWA/g
s/IRIS_313-8658/CP231/g
s/IRIS_313-8669/CALROSE76/g
s/IRIS_313-11800/FANGENG6/g
s/IRIS_313-8665/M7/g
s/IRIS_313-11736/MALAGKIT/g
s/IRIS_313-15910/CYPRESS/g
s/IRIS_313-11790/VARIRANGAHY/g
s/IRIS_313-11924/NAMJAM/g"

colorize="s/(JC1[^']*)'/\1\'[\&!color=#ffdf00\]/g
s/(CHUNDI[^']*)'/\1\'[\&!color=#af00d7\]/g
s/(UCP122[^']*)'/\1\'[\&!color=#af00d7\]/g
s/(JIECAOZHAN[^']*)'/\1\'[\&!color=#0000af\]/g
s/(PSBRC80[^']*)'/\1\'[\&!color=#0000af\]/g
s/(NCS603B[^']*)'/\1\'[\&!color=#0000af\]/g
s/(DISSIGBE[^']*)'/\1\'[\&!color=#005fff\]/g
s/(TSIPALA[^']*)'/\1\'[\&!color=#005fff\]/g
s/(MELABOH[^']*)'/\1\'[\&!color=#00d7ff\]/g
s/(ASFALA[^']*)'/\1\'[\&!color=#00d7ff\]/g
s/(DAGANGZHAN[^']*)'/\1\'[\&!color=#00d7ff\]/g
s/(DERAWA[^']*)'/\1\'[\&!color=#00005f\]/g
s/(CP231[^']*)'/\1\'[\&!color=#008700\]/g
s/(CALROSE76[^']*)'/\1\'[\&!color=#008700\]/g
s/(FANGENG6[^']*)'/\1\'[\&!color=#ff8700\]/g
s/(M7[^']*)'/\1\'[\&!color=#ff8700\]/g
s/(MALAGKIT[^']*)'/\1\'[\&!color=#800000\]/g
s/(CYPRESS[^']*)'/\1\'[\&!color=#800000\]/g
s/(VARIRANGAHY[^']*)'/\1\'[\&!color=#800000\]/g
s/(NAMJAM[^']*)'/\1\'[\&!color=#ff0000\]/g"


< $fasta sed -e "$cultivar_to_names" | sed 's/ <unknown description>//g' > ${fasta}.treeable
rm -f RAxML_${fasta%%.vcf.fasta}.ERROR
rm -f RAxML_${fasta%%.vcf.fasta}.log
rm -f RAxML*${fasta%%.vcf.fasta}*
## Model ASC_GTRGAMMA must be used to correct for the fact that we're only using SNPs

# -f d means rapid hill-climbing algorithm
#raxmlHPC -f d -m ASC_GTRGAMMA --asc-corr=lewis -n ${fasta%%.vcf.fasta} -p 12345 -s ${fasta}.treeable 2>&1 > RAxML_${fasta%%.vcf.fasta}.log || touch RAxML_${fasta%%.vcf.fasta}.ERROR

# -f a means bootstrap analysis and bestTree in one run.
raxmlHPC -f a -m ASC_GTRGAMMA --asc-corr=lewis -n ${fasta%%.vcf.fasta} -N 100 -p 12345 -s ${fasta}.treeable -x 12345 2>&1 > RAxML_${fasta%%.vcf.fasta}.log || touch RAxML_${fasta%%.vcf.fasta}.ERROR

# This can go a lot faster with openmpi
#mpirun -n 6 raxmlHPC-MPI -f a -m ASC_GTRGAMMA --asc-corr=lewis -n ${fasta%%.vcf.fasta} -N 100 -p 12345 -s ${fasta}.treeable -x 12345 2>&1 > RAxML_${fasta%%.vcf.fasta}.log || touch RAxML_${fasta%%.vcf.fasta}.ERROR


#for file in RAxML_bestTree.${fasta%%.vcf.fasta}; do
for file in RAxML_bestTree.${fasta%%.vcf.fasta} RAxML_bipartitions.${fasta%%.vcf.fasta}; do
	# Put single quotes around each taxon
	tree=$(< $file sed -E "s/([^(),']+):/'\1':/g")
	# Count the number of taxa
	ntax=$(echo $tree | grep -o -E "'[^(),']+'" | wc -l)
	# Print out a file in the format FigTree uses
	cat << EOF > ${file}.figtree
#NEXUS
begin taxa;
	dimensions ntax=${ntax};
	taxlabels
	$(echo $tree | grep -o -E "'[^(),']+'")
;
end;

begin trees;
	tree tree_1 = [&R] $tree
end;

begin figtree;
	set appearance.backgroundColorAttribute="Default";
	set appearance.backgroundColour=#ffffff;
	set appearance.branchColorAttribute="User selection";
	set appearance.branchColorGradient=false;
	set appearance.branchLineWidth=1.0;
	set appearance.branchMinLineWidth=0.0;
	set appearance.branchWidthAttribute="Fixed";
	set appearance.foregroundColour=#000000;
	set appearance.hilightingGradient=false;
	set appearance.selectionColour=#2d3680;
	set branchLabels.colorAttribute="User selection";
	set branchLabels.displayAttribute="Branch times";
	set branchLabels.fontName="sansserif";
	set branchLabels.fontSize=8;
	set branchLabels.fontStyle=0;
	set branchLabels.isShown=false;
	set branchLabels.significantDigits=4;
	set layout.expansion=0;
	set layout.layoutType="RADIAL";
	set layout.zoom=0;
	set legend.attribute="bootstrap";
	set legend.fontSize=10.0;
	set legend.isShown=false;
	set legend.significantDigits=4;
	set nodeBars.barWidth=4.0;
	set nodeBars.displayAttribute=null;
	set nodeBars.isShown=false;
	set nodeLabels.colorAttribute="User selection";
	set nodeLabels.displayAttribute="bootstrap";
	set nodeLabels.fontName="sansserif";
	set nodeLabels.fontSize=8;
	set nodeLabels.fontStyle=0;
	set nodeLabels.isShown=false;
	set nodeLabels.significantDigits=4;
	set nodeShape.colourAttribute="User selection";
	set nodeShape.isShown=false;
	set nodeShape.minSize=10.0;
	set nodeShape.scaleType=Width;
	set nodeShape.shapeType=Circle;
	set nodeShape.size=4.0;
	set nodeShape.sizeAttribute="Fixed";
	set polarLayout.alignTipLabels=false;
	set polarLayout.angularRange=0;
	set polarLayout.rootAngle=0;
	set polarLayout.rootLength=100;
	set polarLayout.showRoot=true;
	set radialLayout.spread=0.0;
	set rectilinearLayout.alignTipLabels=false;
	set rectilinearLayout.curvature=0;
	set rectilinearLayout.rootLength=100;
	set scale.offsetAge=0.0;
	set scale.rootAge=1.0;
	set scale.scaleFactor=1.0;
	set scale.scaleRoot=false;
	set scaleAxis.automaticScale=true;
	set scaleAxis.fontSize=8.0;
	set scaleAxis.isShown=false;
	set scaleAxis.lineWidth=1.0;
	set scaleAxis.majorTicks=1.0;
	set scaleAxis.origin=0.0;
	set scaleAxis.reverseAxis=false;
	set scaleAxis.showGrid=true;
	set scaleBar.automaticScale=true;
	set scaleBar.fontSize=10.0;
	set scaleBar.isShown=true;
	set scaleBar.lineWidth=1.0;
	set scaleBar.scaleRange=0.0;
	set tipLabels.colorAttribute="User selection";
	set tipLabels.displayAttribute="Names";
	set tipLabels.fontName="sansserif";
	set tipLabels.fontSize=8;
	set tipLabels.fontStyle=0;
	set tipLabels.isShown=true;
	set tipLabels.significantDigits=4;
	set trees.order=false;
	set trees.orderType="increasing";
	set trees.rooting=false;
	set trees.rootingType="User Selection";
	set trees.transform=false;
	set trees.transformType="cladogram";
end;

EOF
	< ${file}.figtree sed -r -e "$colorize" > ${file}.figtree.colorized
    done
