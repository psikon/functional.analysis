#!/usr/bin/env python

# Author: Philipp Sehnert
# Contact: psikon86[a]gmail.com

# script to remap from representive sequences, that have been annotated with hmmsearch, 
# back to all sequences clustered together during cd-hit run.

# command: remap_cluster2pfam.py -p [pfam .hmmsearch file] -c [.clstr cluster file] -o [output file]

# the output is a hmmsearch file in table format containing all sequences of a cluster with the 
# same pfam dataset like the representive sequence

# import modules used by the script
import argparse, os, itertools, re

# set argument parser
parser = argparse.ArgumentParser(description = 'remap from representive sequences, annotated by hmmsearch,\n to sequences also occuring in the cluster.')

parser.add_argument('-p', '--pfam', metavar='.pfam file', dest='pfam', type=str,
					help='The .pfam file containing the annotated pfams from hmmsearch in table output.')
parser.add_argument('-c', '--cluster', metavar='.clstr file', dest='cluster', type=str,
					help='The .clstr file produced by CD-hit that contains the cluster information.')
parser.add_argument('-o', '--output', metavar='output file', dest='output', type=str,
					help='output file')

args = parser.parse_args()

'''parse through the .clstr file and create a dictionary
   with the sequences per cluster'''
def read_clstr():
	
		# open the cluster file and set the output dictionary
		cluster_file, cluster_dic = open(args.cluster), {}
		
		# create an iteratable object for every cluster
		cluster_groups = (x[1] for x in itertools.groupby(cluster_file, 
			key=lambda line: line[0] == '>'))
		# parse through the cluster and create a dictionary with 
		# representive sequence as key and other sequences as values
		for cluster in cluster_groups:
			# get name of the cluster
			name = cluster.next().strip()
			# store the actual iiteration
			iteration = cluster_groups.next()
			# init value array
			seqs = []
			# iterate over all sequences in a cluster
			for seq in iteration:
				# get the representive sequence
				if '*' in seq:
					rep = seq.split('>')[1].split('...')[0]
				# store the non representive contents of cluster
				seqs.append(seq.split('>')[1].split('...')[0])
			# create new entry in dictionary
			cluster_dic[rep] = seqs
		
		# return the cluster dictionary
		return cluster_dic
'''parse the content of a hmmsearch file in table output and store it in a dictionary'''
def read_pfam():

	# open the pfam file and init the pfam dictionary
	pfam_file, pfam_dict = open(args.pfam,'r'), {}
	# parse the pfam file line by line
	for line in pfam_file:
		# remove commentated lines
		if '#' not in line:
			# create a dictionary with target sequence as key and the rest of line as value
			pfam_dict[re.split(r'\ +', line)[0]] = re.split(r'\ +', line)
	# return the pfam dictionary
	return pfam_dict

'''compare the keys in the pfam dictionary with the keys in the cluster dictionary and
   create for every matched key new lines from the pfam dictionary with new target sequences'''
def create_pfam_abundance(cluster_dic, pfam_dic, output):
	# open the output file
	with open(output, 'wb') as out:
		# iterate over the keys in cluster dictionary
		for key in cluster_dic:
			# matched key
			if key in pfam_dic:
				# iterate over the values in cluster dictionary
				for id in cluster_dic[key]:
					# get the line of the pfam dictionary
					l = pfam_dic[key][1:len(pfam_dic[key])]
					# create a tab seperated string of the line
					pfam = '\t'.join([str(x) for x in l])
					# append the new created line to the output file
					# id = new target sequence, pfam = original pfam line
					out.write(id + '\t' + pfam)
	
def main():
		
		# obtain a dictionary with the clusters and sequences
		cluster_dic = read_clstr()
		# obtain a dictionary of all pfam lines
		pfam_dic = read_pfam()
		# compare the two dictionaries and create a new pfam file with more target sequences
		create_pfam_abundance(cluster_dic, pfam_dic, args.output)
		

if __name__ == '__main__':
	main()