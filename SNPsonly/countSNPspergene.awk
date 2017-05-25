#!/bin/awk -f

' BEGIN {
	for (i=0; i<=32914; i++) {
		genes["gene" i]=0
	}
}
/^chr/ && /gene[0123456789]+/ {
	match($8, /gene[0123456789]+/, array);
	if (array[0] in genes) {
		genes[array[0]]++
	}
	else {
		genes[array[0]]=1
	}
}
END {
	for (gene in genes) {
		print gene, genes[gene];
		sum+=genes[gene]
	};
#	print "#average of 32915 genes:", sum/32915
}'
