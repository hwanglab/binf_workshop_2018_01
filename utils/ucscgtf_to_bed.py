#!/usr/bin/env python

import fileinput
import sys,re
'''
1       hg19_refGene    start_codon     67000042        67000044        0.000000        +       .       gene_id "NM_001350218"; transcript_id "NM_001350218"; 
1       hg19_refGene    CDS     67000042        67000051        0.000000        +       0       gene_id "NM_001350218"; transcript_id "NM_001350218"; 
1       hg19_refGene    exon    66999639        67000051        0.000000        +       .       gene_id "NM_001350218"; transcript_id "NM_001350218"; 
'''
for line in fileinput.input():
	i = line.split('\t')
	mObj = re.search(r"transcript_id\s\"(\w+)\"",i[8])
	if not mObj:
		raise RuntimeError('check gtf file stdin')
	ele = [j.strip() for j in [i[0], i[3], i[4], mObj.group(1),'0',i[6]]]
	sys.stdout.write("%s\n"%'\t'.join(ele))