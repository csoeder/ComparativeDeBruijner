# shamelessly ripped off of Pall Melsted
# https://pmelsted.wordpress.com/2013/11/23/naive-python-implementation-of-a-de-bruijn-graph/

import collections
import sys
from Bio import Seq, SeqIO, SeqRecord

# functions to deal with k-mers
def twin(km):
	#	return reverse complement of input sequence km
	return Seq.reverse_complement(km)

def kmers(seq,k):
	# return a generator which outputs the k-mers of input seq
	for i in xrange(len(seq)-k+1):
		yield seq[i:i+k]
 
def fw(km):
	#	for a given k-mer, generate the 4 possible successor k-mers
	for x in 'ACTG':
		yield km[1:]+x

def bw(km):
	#	for a given k-mer, generate the 4 possible predecessor transition k-mers
	for x in 'ACTG':
		yield x + km[:-1]



