# shamelessly ripped off of Pall Melsted
# 	https://pmelsted.wordpress.com/2013/11/23/naive-python-implementation-of-a-de-bruijn-graph/
# visualization ref:
# 	https://plot.ly/ipython-notebooks/network-graphs/
import collections
import sys
from Bio import Seq, SeqIO, SeqRecord
import matplotlib.pyplot as plt
sys.path.append('/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/')
import networkx as nx


# functions to deal with k-mers


def twin(km):
	# 	return reverse complement of input sequence km
	return Seq.reverse_complement(km)


def kmers(seq, k):
	# 	return a generator which outputs the k-mers of input seq
	for i in xrange(len(seq) - k + 1):
		yield seq[i:i+k]


def fw(km):
	# 	for a given k-mer, generate the 4 possible successor k-mers
	for x in 'ACTG':
		yield km[1:]+x


def bw(km):
	# 	for a given k-mer, generate the 4 possible predecessor transition k-mers
	for x in 'ACTG':
		yield x + km[:-1]


def build(fn, k=31, limit=0, format='q'):
	# 	a function to build a k-mer dictionary from a FASTA filelist, fn
	d = collections.defaultdict(int)

	for f in fn:
		reads = SeqIO.parse(f, 'fast' + format)
		for read in reads:
			seq_s = str(read.seq)
			seq_l = seq_s.split('N')
			for seq in seq_l:
				for km in kmers(seq, k):
					d[km] += 1
				seq = twin(seq)
				for km in kmers(seq, k):
					d[km] += 1
	# 	noise trimming
	d1 = [x for x in d if d[x] <= limit]
	for x in d1:
		del d[x]
	return d


def graph_scraper(seq, k=11):
	murrs = collections.defaultdict(int)
	transitions = collections.defaultdict(int)
	seq_split = seq.split('N')
	for subseq in seq_split:
		jenny = kmers(subseq, k)

		new_kmer = jenny.next()
		murrs[new_kmer] += 1
		old_kmer = new_kmer

		for new_kmer in jenny:
			murrs[new_kmer] += 1
			transitions[(old_kmer, new_kmer)] += 1
			old_kmer = new_kmer

		jenny = kmers(twin(subseq), k)

		new_kmer = jenny.next()
		murrs[new_kmer] += 1
		old_kmer = new_kmer

		for new_kmer in jenny:
			murrs[new_kmer] += 1
			transitions[(old_kmer, new_kmer)] += 1
			old_kmer = new_kmer
	return murrs, transitions


nodes, edges = graph_scraper("ATGCCGTA", k=4)
G = nx.Graph()
for node in nodes:
	G.add_node(node)
for edge in edges:
	G.add_edge(edge[0], edge[1])

nx.draw(G)
plt.show()










