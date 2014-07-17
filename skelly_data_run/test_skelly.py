import pysam
import itertools

byfile = pysam.Samfile('by.sam')
rmfile = pysam.Samfile('rm.sam')


i = 0
for by, rm in itertools.izip(byfile, rmfile):
	if len(by.seq) != len(rm.seq):
		print ('at {} : {} != {} : {} {}'.format(i, len(by.seq), len(rm.seq), by.qname, rm.qname))
	i+=1
