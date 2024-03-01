#2023.8.24 Ruojia Wu: sort reads by Unique Sample Barcodes, count sequence types in each USB, and write each USB to a file

from collections import defaultdict
import pysam

from collections import Counter

import difflib 

from argparse import ArgumentParser
import os 
from math import ceil

def main():

	# parser

	parser = ArgumentParser()

	parser.add_argument("--BamName", action="store", dest="input_bam_name", help="Input Bam file name", required=True)

	o = parser.parse_args()

	caller(o.input_bam_name)


def caller(input_bam_name):

	usb_len = 8
	same_seq_min = 3 #a calling should have more (>=) reads than this
	usb_family_min = 10000 #each USB needs more (>=) reads than this to be printed

	usb_list = {}


	#bamname = 'aligntest.bam'

	bam = pysam.AlignmentFile(input_bam_name, 'rb')

# sort reads by usb
	readcount = 0
	for read1, read2 in read_pair_generator(bam):
		readcount = readcount+1
		if readcount%100000 == 0:
			print(str(readcount)+' Reads processed')

		read_usb = read1.query_name[-usb_len:]
		read_nt = read1.query_sequence + 'xx' + read2.query_sequence
		read_cigar = read1.cigarstring + 'xx' + read2.cigarstring
		read_info = read_cigar + 'xx' + read_nt

		if read_usb not in usb_list:
			usb_list[read_usb] = list([read_info])
		else:
			usb_list[read_usb].append(read_info)


	print('Writing...')
	# list contents for each usb; write each usb to a file (remove small usb families)
	for curr_usb in usb_list:

		#print(len(usb_list[curr_usb]))

		# only keep usbs with enough size (real samples, not mutated usbs)
		cur_usb_family = len(usb_list[curr_usb])
		#print(cur_usb_family)

		if cur_usb_family >= usb_family_min:
			usb_list3 = Counter(usb_list[curr_usb])

			final_outfile_name = input_bam_name[:-4]+curr_usb+'_ResultSummary.txt'
			no_indel_count = 0

			with open(final_outfile_name,'a') as ff:
				ff.write('CIGAR1'+'\t'+'CIGAR2'+'\t'+'ReadCount'+'\t'+'Seq1'+'\t'+'Seq2'+'\n')
				for curmutinfo, curcount in usb_list3.most_common():
					if curcount >= same_seq_min:
						if curmutinfo[:12]=='242Mxx242Mxx': #no indel
							no_indel_count = no_indel_count + curcount
						else:
							xxpos1 = curmutinfo.find('xx')
							xxpos2 = curmutinfo.find('xx',xxpos1+1)
							xxpos3 = curmutinfo.find('xx',xxpos2+1)
							ff.write(curmutinfo[:xxpos1]+'\t'+curmutinfo[xxpos1+2:xxpos2]+'\t'+str(curcount)+'\t'+curmutinfo[xxpos2+2:xxpos3]+'\t'+curmutinfo[xxpos3+2:]+'\n')
							#ff.write(curmutinfo +'\t'+str(curcount)+'\n')
				ff.write('242M'+'\t'+'242M'+'\t'+str(no_indel_count)+'\n')

	print('Done')

 	

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]






main()