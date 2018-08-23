#!/usr/bin/env python

"""
exitron-pipeline.py => Exitron identification, quantification and comparison pipeline
"""

from natsort import natsorted
from collections import Counter
import argparse
import subprocess
import os

def yield_junctions(f, file_format="STAR"):

	for line in open(f):
		try:

			cols = line.rstrip().split('\t')

			if file_format == "TopHat2":
				yield({
					'chr': cols[0],
					'start': cols[1],
					'end': cols[2],
					'junction_id': cols[3],
					'uniq_reads': int(cols[4]),
					'strand': cols[5],
					'rgb': cols[8],
					'block_size': cols[9],
					'blocks': cols[10],
					'junc_start': int(cols[1]) + int(cols[10].split(',')[0]),
					'junc_end': int(cols[2]) - int(cols[10].split(',')[1]),
					'misc': cols[11]
					})

			elif file_format == "STAR":
				yield({
					'chr': cols[0],
					'junc_start': int(cols[1]),
					'junc_end': int(cols[2]),
					'strand': { '0': 'undefined', '1': '+', '2': '-' }[cols[3]],
					'intron_motif': { '0': 'non-canonical', '1': 'GT-AG', '2': 'CT-AC', '3': 'GC-AG', '4': 'CT-GC', '5': 'AT-AC', '6': 'GT-AT' }[cols[4]],
					'is_annotated': cols[5] == '1',
					'uniq_reads': int(cols[6]),
					'multimap_reads': int(cols[7]),
					'max_overhang': int(cols[8])
					})

		except IndexError: pass

def get_shared_junctions(group_junctions, min_support=3):

	# Count the number of times a junction has been found
	jN = Counter({})
	for sample in group_junctions:
		for j in sample:
			if sample[j] >= 0:
				jN[j] += 1

	# Return those junctions with >= min_support occurences
	return set([ x for x in jN if jN[x] >= min_support ])

def identify_exitrons(work_dir, gtf_file, samples, junction_format="STAR", junction_filename="SJ.out.tab", min_support=3):
	""" Intersect the shared junctions with min_support occurences per 
	samples group with the CDS GTF annotation. """

	# Parse the samples per group
	groups = {}
	for line in open(samples):
		group_id, path = line.rstrip().split('\t')[:2]
		try: groups[group_id].append( { "{}:{}-{}:{}".format(j["chr"], j["junc_start"], j["junc_end"], j["strand"]): j["uniq_reads"] for j in yield_junctions(path+junction_filename, junction_format) } )
		except KeyError: groups[group_id] = [ { "{}:{}-{}:{}".format(j["chr"], j["junc_start"], j["junc_end"], j["strand"]): j["uniq_reads"] for j in yield_junctions(path+junction_filename, junction_format) } ]

	# Get all junctions with at least min_support occurences per group
	all_junctions = set.union(*[ get_shared_junctions(groups[x], min_support) for x in groups ])

	# Output all selected junctions in bed format
	with open("{}all_junctions.N{}.bed".format(work_dir, min_support), 'w') as fout:
		for j in natsorted(all_junctions):
			c, coord, strand = j.split(':')
			fout.write( "{}\t{}\t{}\tJUNCTION\t1000\t{}\n".format(c, coord.split('-')[0], coord.split('-')[1], strand) )

	# Intersect the junctions with the provided CDS GTF annotation
	# Only strand-specific (-s) and full-length matches (-f 1) are taken
	subprocess.call("bedtools intersect -s -f 1 -wa -wb -a {0}all_junctions.N{1}.bed -b {2} | awk -F\"\\t\" '{{ OFS=\"\t\"; print $1, $2, $3, $6, $10, $11, $NF }}' > {0}junctions_GTF_map.tmp".format(work_dir, min_support, gtf_file), shell=True)

	# Parse intersection
	seen = set()
	infoOut, bedOut = open(work_dir+"exitrons.info", 'w'), open(work_dir+"exitrons.bed", 'w')
	infoOut.write( "#exitron_id\ttranscript_id\tgene_id\tgene_name\tEI_length_in_nt\tEIx3\n" )

	for line in open(work_dir+"junctions_GTF_map.tmp"):
		j_chr, j_start, j_end, j_strand, cds_start, cds_end, attributes = line.rstrip().split('\t')
		attr = {}
		for a in attributes.split(';'):
			if len(a):
				attr_name, attr_value = list(filter(None, a.split(' ')))
				attr[attr_name.strip()] = attr_value.replace('\"', '')

		EI_length = abs(int(j_end)-int(j_start)) + 1
		EIx3 = "yes" if EI_length % 3 == 0 else "no"
		exitron_id = "{}:{}-{}:{}".format(j_chr, j_start, j_end, j_strand)

		if not exitron_id in seen:
			infoOut.write( "{}\t{}\t{}\t{}\t{}\t{}\n".format(exitron_id, attr["transcript_id"], attr["gene_id"], attr["gene_name"], EI_length, EIx3) )
			bedOut.write( "{}\t{}\t{}\t{};{}\t1000\t{}\n".format(j_chr, j_start, j_end, attr["gene_name"], exitron_id, j_strand) )
			seen.add(exitron_id)

	infoOut.close(), bedOut.close()

	subprocess.call("rm -f "+work_dir+"junctions_GTF_map.tmp", shell=True)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-v', '--version', action='version', version='0.1.1')
	parser.add_argument('-w', '--work-dir', default="./", help="Output working directory.")

	subparsers = parser.add_subparsers(dest='command', help="Sub-command help.")

	# Identify exitrons
	parser_a = subparsers.add_parser('identify-exitrons', help="Map all junction with minimal --min-support occurences to the provided GTF annotation.")
	parser_a.add_argument('-g', '--gtf', required=True, help="CDS annotation gtf file.")
	parser_a.add_argument('--samples', required=True, help="Tab-separated file containing the group id (e.g. wt), absolute path to the folder containing the mapping bam and junction files, and sample id (e.g. wt-1).")
	parser_a.add_argument('--junction-format', default="STAR", choices=["STAR", "TopHat2"], help="Junction file format.")
	parser_a.add_argument('--junction-filename', default="SJ.out.tab", help="Junction filename.")
	parser_a.add_argument('--min-support', type=int, default=3, help="Minimum number of replicates a junctions needs to occur in per group.")

	args = parser.parse_args()

	# Create work-dir if not exists
	work_dir = args.work_dir if args.work_dir.endswith('/') else args.work_dir + '/'
	if not os.path.exists(work_dir): os.makedirs(work_dir)

	if args.command == "identify-exitrons":
		identify_exitrons(work_dir, args.gtf, args.samples, args.junction_format, args.junction_filename, args.min_support)