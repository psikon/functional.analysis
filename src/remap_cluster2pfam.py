#!/usr/bin/env python

# Author: Youri Lammers
# Contact: youri.lammers@naturalis.nl / youri.lammers@gmail.com

# Filter the cluster output from CD-hit based on
# the minimum number of sequences per cluster

# command: CD-hit_filter.py -f [cluster .fasta file] -c [.clstr cluster file] -m [min # seqs per cluster]

# the output of the CD-hit_filter.py script is a fasta file [input cluster .fasta file with the min # seqs in the file name]
# with the OTU sequences for the clusters with more than or equal number of sequences to the user specified minimum.

# import modules used by the script
import argparse, os, itertools
import re

# set argument parser
parser = argparse.ArgumentParser(description = 'Filter the output from CD-hit based on the minimum number of read per cluster.\nThe filtered output fasta file produced has the same name as the input file with _min_[minimum size from -c argument].fasta attachted to the name.')

parser.add_argument('-f', '--fasta', metavar='.fasta file', dest='fasta', type=str,
					help='The .fasta file containing the clusters produced by CD-hit.')
parser.add_argument('-p', '--pfam', metavar='.pfam file', dest='pfam', type=str,
					help='The .pfam file containing the annotated pfams from hmmsearch in table output.')
parser.add_argument('-c', '--cluster', metavar='.clstr file', dest='cluster', type=str,
					help='The .clstr file producec by CD-hit that contains the cluster information.')
parser.add_argument('-m', '--minimum', metavar='minimum size', dest='minimum', type=int,
					help='The minimum cluster size.')

args = parser.parse_args()

def read_clstr():
		# parse through the .clstr file and create a dictionary
		# with the sequences per cluster
		
		# open the cluster file and set the output dictionary
		cluster_file, cluster_dic = open(args.cluster), {}
		
		# parse through the cluster file and store the cluster name + sequences in the dictionary
		cluster_groups = (x[1] for x in itertools.groupby(cluster_file, key=lambda line: line[0] == '>'))
		for cluster in cluster_groups:
			name = cluster.next().strip()
			iteration = cluster_groups.next()
			seqs = []
			for seq in iteration:
				if '*' in seq:
					rep = seq.split('>')[1].split('...')[0]
				seqs.append(seq.split('>')[1].split('...')[0])
			#seqs = [seq.split('>')[1].split('...')[0] for seq in iteration]
			cluster_dic[rep] = seqs
		
		# return the cluster dictionary
		return cluster_dic

def read_pfam():

	# open the pfam file and set the pfam dictionary
	pfam_file, pfam_dict = open(args.pfam,'r'), {}
	
	for line in pfam_file:
		if '#' not in line:
			pfam_dict[re.split(r'\ +', line)[0]] = re.split(r'\ +', line)
	return pfam_dict


def create_pfam_abundance(cluster_dic, pfam_dic):
	with open("data/newout.txt", 'wb') as out:
		for key in cluster_dic:
			if key in pfam_dic:
				for id in cluster_dic[key]:
					l = pfam_dic[key][1:len(pfam_dic[key])]
					pfam = '\t'.join([str(x) for x in l])
					out.write(id + '\t' + pfam)
	
def main():
		
		# obtain a dictionary with the clusters and sequences
		cluster_dic = read_clstr()
		# read pfam file 
		pfam_dic = read_pfam()
		
		create_pfam_abundance(cluster_dic, pfam_dic)
		

if __name__ == '__main__':
	main()