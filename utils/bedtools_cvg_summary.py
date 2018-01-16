#!/usr/bin/env python

import os, gzip, argparse

def open2(fn,open_mode):
	fp = None
	if open_mode == 'r':
		if not os.path.exists(fn):
			return fp
		if fn.endswith('.gz'):
			try:
				fp=gzip.open(fn,'r')
			except IOError:
				try:
					fp=gzip.open(fn,'rb')
				except IOError:
					raise IOError('cannot open the file[%s] to read'%fn)
		else:
			fp=open(fn,'r')
	elif open_mode == 'w':
		if fn.endswith('.gz'):
			try:
				fp=gzip.open(fn,'w')
			except IOError:
				try:
					fp=gzip.open(fn,'wb')
				except IOError:
					raise IOError('cannot open the file[%s] to write'%fn)
		else:
			fp=open(fn,'w')
	elif open_mode == 'a':
		if fn.endswith('.gz'):
			try:
				fp=gzip.open(file,'a')
			except IOError:
				try:
					fp=gzip.open(file,'ab')
				except IOError:
					raise IOError('cannot open the file[%s] to append'%file)
		else:
			fp=open(fn,'a')
	else:
		raise IOError('check file[%s] and open_mode[%s]'%(fn,open_mode))
	return fp

def bedHist2cvgReport(bedhist_exon,eCvgFn,rd=5):

	if not os.path.exists(bedhist_exon):
		raise IOError('check if [%s] exists'%bedhist_exon)

	if bedhist_exon.endswith('gz'):
		cmd = "zgrep -v ^all %s" % bedhist_exon
	else:
		cmd  = "grep -v ^all %s" % bedhist_exon
	tmp_fn = "%s.tmp"%bedhist_exon
	cmd += " | sort -k4,4V -k1,1V -k2,2n -k3,3n > %s"%(tmp_fn)
	os.system(cmd)

	fp=open2(tmp_fn,'r')
	fp2=open2(eCvgFn,'w')

	fp2.write('#bedtools_hist_file\t%s\n'%bedhist_exon)
	fp2.write('#%s\ttotal_bp\tcoverage_%dx[%%]\tmean_depth\n' % ("transcript_id",rd))

	visit1=True
	pChr,pStbp,pEdbp,pAnnot=['.',0,0,'.']
	P = 100.
	eCvg = 0.
	cReads = 0.
	total_bp = 1

	for i in fp:
		if i.startswith('all'):break
		chrm,stBp,edBp,annot,readx,bpCvd,eLen = parse_bedtools_exon_line(i)

		if annot!=pAnnot:
			if visit1:
				visit1=False
			else: #wrapup prevRecord
				fp2.write('%s\t%g\t%g\t%g\n'%(pAnnot, total_bp, P*eCvg/total_bp,1.*cReads/total_bp))
			eCvg = 0. #reset bp cvg count
			cReads = 0
			total_bp = 0
			pAnnot = annot

		if readx >= rd:
			eCvg += bpCvd

		#cumulative reads
		cReads+=(bpCvd*readx)

		total_bp += eLen


	#wraup the final instance
	fp2.write('%s\t%g\t%g\t%g\n' % (pAnnot, total_bp, P * eCvg / total_bp, 1. * cReads / total_bp))
	fp.close()
	fp2.close()
	#os.unlink(tmp_fn)

def parse_bedtools_exon_line(bedtoolHistExonLine):

	chrm,stBp,edBp,annot,_,_,readx,bpCvd,eLen,cvdRate\
		=bedtoolHistExonLine.rstrip().split('\t')
	
	[stBp,edBp,readx,bpCvd,eLen] = map(int,[stBp,edBp,readx,bpCvd,eLen])
	return [chrm,stBp,edBp,annot,readx,bpCvd,eLen]


def main():
	parser = argparse.ArgumentParser(description="bedtools hist to exon coverage summary")
	parser.add_argument('-i', '--bedtools_hist', dest='bedtools_hist', required=True, help='bedtools coverage histogram output file')
	parser.add_argument('-o', '--out_file', dest='out_file', required=False, default=None, help='summary output file')

	parser.add_argument('-d', '--min_depth', action='store', dest='min_depth', required=False, default=5, type=int,
											help='minimum read depth of coverage to consider covering [5]')

	args = parser.parse_args()
	bedHist2cvgReport(args.bedtools_hist,args.out_file,rd=args.min_depth)

if __name__ == '__main__':
	main()