#!/bin/awk -f

BEGIN {
	for (i=0; i<=32914; i++) {
		promoters["promoter" i]=0
	}
}
/^chr/ && /promoter[0123456789]+/ {
	match($8, /promoter[0123456789]+/, array);
	if (array[0] in promoters) {
		promoters[array[0]]++
	}
	else {
		promoters[array[0]]=1
	}
}
END {
	command="sort -V" ;
	for (promoter in promoters) {
		print promoter, promoters[promoter] | command;
		sum+=promoters[promoter]
	};
	print "#average of 32915 promoters:", sum/32915
}
