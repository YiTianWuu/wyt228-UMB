# Find Unique Sample Barcodes (4nt) in R1 adn R2, merge and put in R1 and R2 names, remove reads w/o FLT3 primers

from argparse import ArgumentParser

def main():

	# parser

	parser = ArgumentParser()

	parser.add_argument("--FastqR1", action="store", dest="input_name1", help="Input FASTQ R1 file", required=True)

	parser.add_argument("--FastqR2",  action="store",  dest="input_name2", help="Input FASTQ R2 file", required=True)

	parser.add_argument("--OutputName1",  action="store",  dest="output_name1", help="Output trimmed FASTQ file R1", required=True)
	
	parser.add_argument("--OutputName2",  action="store",  dest="output_name2", help="Output trimmed FASTQ file R2", required=True)

	o = parser.parse_args()

	clean_fastq_wPrimer(o.input_name1, o.input_name2, o.output_name1, o.output_name2)


def reverse_complement(seq):
	complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c'}
	bases = list(seq)
	bases = [complement[base] for base in bases]
	return ''.join(bases[::-1])


def clean_fastq_wPrimer(input_name1, input_name2, output_name1, output_name2):
	ITDfP1 = 'TTCCTAACTGACTCATCATTTCATCTCT'
	ITDrP1 = 'ATCTTTGTTGCTGTCCTTCCACTATA'
	from itertools import izip
	with open(input_name1, "r") as r1, open(input_name2, "r") as r2, open(output_name1, "w") as w1, open(output_name2, "w") as w2:
		counter = 0
		last = False
		bunchmax = 1000000
		bunchsize = 0
		bunchnum = 0
		bunch1 = []
		bunch2 = []
		#minlen = 60
		#cutlen = 60
		for l1, l2 in izip(r1, r2):
			if counter%4 == 0: #name
				#create two read elements
				read1 = ['','','','']
				read2 = ['','','','']
				read1[0] = l1
				read2[0] = l2
			elif counter%4 == 1:#seq
				#seq1 = l1[:]
				#umi_seq1 = l2[4:19] #UMI
				#seq2 = l2[19:]

				seq1 = l1[8:]
				usb_seq1 = l1[4:8]+l2[4:8] #unique sample barcodes left + right
				seq2 = l2[8:]

				to_write = 1
				if 'N' in seq1 or 'N' in usb_seq1 or 'N' in seq2:
					to_write = 0
				if (ITDfP1 not in seq1) or (ITDrP1 not in seq2):
					to_write = 0
				if to_write == 1:
					read1[1] = seq1
					spacepos = read1[0].find(' ')
					read1[0] = read1[0][:spacepos]+'USB'+usb_seq1+'\n'
					read2[1] = seq2
					spacepos2 = read2[0].find(' ')
					read2[0] = read2[0][:spacepos2]+'USB'+usb_seq1+'\n'
			elif counter%4 == 2:#+
				read1[2] = l1
				read2[2] = l2
			elif counter%4 == 3:#Qscore
				if to_write == 1:
					#read1[3] = l1[:]
					#read2[3] = l2[19:]
					read1[3] = l1[8:]
					read2[3] = l2[8:]
					bunch1.extend(read1)
					bunch2.extend(read2)
					bunchsize += 1

			if bunchsize > bunchmax:
				w1.writelines(bunch1)
				w2.writelines(bunch2)
				bunchsize = 0 
				bunchnum += 1
				bunch1 = []
				bunch2 = []
			counter += 1 
		w1.writelines(bunch1)
		w2.writelines(bunch2)
		print counter/4
		print bunchnum*bunchmax+bunchsize
		print ''


main()
