#!/usr/bin/env python

"""
exitron-pipeline.py => Exitron identification, quantification and comparison pipeline
"""

from natsort import natsorted
import argparse
import subprocess
import time
import os

def log_settings(work_dir, args, write_mode='a'):

	args_order = ['version', 'work_dir', 'command', 'gtf', 'samples', 'junction_format', 'junction_filename', 'min_count', 'min_coverage', 'min_total_coverage', 'bam', 'file-handle', 'genome-fasta', 'bam-filename', 'NPROC']

	with open(work_dir+"Log.out", write_mode) as fout:
		fout.write(time.asctime( time.localtime(time.time()) )+"\n")
		for arg in args_order:
			try: fout.write("--{}\t{}\n".format(arg.replace('_', '-'), getattr(args, arg)))
			except AttributeError: pass

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
					'depth': int(cols[4]),
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
					'max_overhang': int(cols[8]),
					'depth': int(cols[6]) + int(cols[7])
					})

		except IndexError: pass

def get_shared_junctions(group_junctions, min_count=3, min_coverage=1, min_total_coverage=1):

	# Count the number of times a junction has been found
	jN = {}
	for sample in group_junctions:
		for j in sample:
			if sample[j] >= min_coverage:
				try: 
					jN[j]['N'] += 1
					jN[j]['total_cov'] += sample[j]
				except KeyError:
					jN[j] = { 'N': 1, 'total_cov': sample[j] }

	# Return those junctions with >= min_count occurences
	return set([ x for x in jN if jN[x]['N'] >= min_count and jN[x]['total_cov'] >= min_total_coverage ])

def identify_exitrons(work_dir, args):
	""" Intersect the shared junctions with min_count occurences per 
	samples group with the CDS GTF annotation. """

	# Parse the samples per group
	groups = {}
	for line in open(args.samples):
		group_id, path = line.rstrip().split('\t')[:2]
		junction_file = path+args.junction_filename
		try: groups[group_id].append( { "{}:{}-{}:{}".format(j["chr"], j["junc_start"], j["junc_end"], j["strand"]): j["depth"] for j in yield_junctions(junction_file, args.junction_format) } )
		except KeyError: groups[group_id] = [ { "{}:{}-{}:{}".format(j["chr"], j["junc_start"], j["junc_end"], j["strand"]): j["depth"] for j in yield_junctions(junction_file, args.junction_format) } ]

	# Get all junctions with at least min_count occurences per group
	all_junctions = set.union(*[ get_shared_junctions(groups[x], args.min_count, args.min_coverage, args.min_total_coverage) for x in groups ])

	# Output all selected junctions in bed format
	with open(work_dir+"all_junctions.bed", 'w') as fout:
		for j in natsorted(all_junctions):
			c, coord, strand = j.split(':')
			fout.write( "{}\t{}\t{}\tJUNCTION\t1000\t{}\n".format(c, coord.split('-')[0], coord.split('-')[1], strand) )

	# Intersect the junctions with the provided CDS GTF annotation
	# Only strand-specific (-s) and full-length matches (-f 1) are taken
	subprocess.call("bedtools intersect -s -f 1 -wa -wb -a {0}all_junctions.bed -b {1} | awk -F\"\\t\" '{{ OFS=\"\t\"; print $1, $2, $3, $6, $10, $11, $NF }}' > {0}junctions_GTF_map.tmp".format(work_dir, args.gtf), shell=True)

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

def prepare_bam_files(work_dir, bam_file, handle, genome_fasta, NPROC=4):
	""" Extract the unique reads and unique exonic reads from the 
	supplied bam file, index and output to the working directory. """

	# Extract unique junctions reads (ujr)
	uniq_reads_bam = "{}ujr.{}.bam".format(work_dir, handle) 
	cmd = "samtools view -@ %d %s | grep -w \"NH:i:1\" | perl -n -e '@line=split(/\\t/,$_); if ($line[5]=~/N/){ print \"$_\"; }' > %s" % (NPROC, bam_file, uniq_reads_bam.replace('.bam', '.sam'))
	subprocess.call(cmd, shell=True)

	# Convert sam to bam
	subprocess.call("samtools view -@ {0} -bT {1} {2} | samtools sort -o {3} -".format(NPROC, genome_fasta, uniq_reads_bam.replace('.bam', '.sam'), uniq_reads_bam), shell=True)

	# Index bam
	subprocess.call("{} index {}".format(samtools_path, uniq_reads_bam), shell=True)

	# Remove sam
	subprocess.call("rm -f {}".format(uniq_reads_bam.replace('.bam', '.sam')), shell=True)

	# Extract unique exonic reads (uer)
	uniq_exon_reads_bam = '{}uer.{}.bam'.format(work_dir, handle)
	cmd = "samtools view -@ %d %s | grep -w \"NH:i:1\" | perl -n -e '@line=split(/\\t/,$_); if ($line[5]!~/N/){ print \"$_\"; }' > %s" % (NPROC, bam_file, uniq_exon_reads_bam.replace('.bam', '.sam'))
	subprocess.call(cmd, shell=True)

	# Convert sam to bam
	subprocess.call("samtools view -@ {0} -bT {1} {2} | samtools sort -o {3} -".format(NPROC, genome_fasta, uniq_exon_reads_bam.replace(".bam", ".sam"), uniq_exon_reads_bam), shell=True)

	# Index bam
	subprocess.call("{} index {}".format(samtools_path, uniq_exon_reads_bam), shell=True)

	# Remove sam
	subprocess.call("rm -f {}".format(uniq_exon_reads_bam.replace('.bam', '.sam')), shell=True)

def prepare_multi_bam(work_dir, args):
	""" Loop through the samples in the file and 
	prepare the bam files. """

	for line in open(args.samples):
		group_id, path, sample_id = line.rstrip().split('\t')[:3]
		prepare_bam_files(work_dir, path+args.bam_filename, sample_id, args.genome_fasta, args.NPROC)

if __name__ == '__main__':

	def stay_positive(val):
		if int(val) < 1:
			raise argparse.ArgumentTypeError("--NPROC must be at least 1.")
		return val

	version = "0.2.0"
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-v', '--version', action='version', version=version, default=version)
	parser.add_argument('-w', '--work-dir', default="./", help="Output working directory.")

	subparsers = parser.add_subparsers(dest='command', help="Sub-command help.")

	# Identify exitrons
	parser_a = subparsers.add_parser('identify-exitrons', help="Map junctions to the provided GTF annotation.")
	parser_a.add_argument('-g', '--gtf', required=True, help="CDS annotation gtf file.")
	parser_a.add_argument('--samples', required=True, help="Tab-separated file containing the group id (e.g. wt), absolute path to the folder containing the mapping bam and junction files, and sample id (e.g. wt-1).")
	parser_a.add_argument('--junction-format', default="STAR", choices=["STAR", "TopHat2"], help="Junction file format.")
	parser_a.add_argument('--junction-filename', default="SJ.out.tab", help="Junction filename.")
	parser_a.add_argument('--min-count', type=int, default=3, help="Minimum number of replicates a junctions needs to occur in per group.")
	parser_a.add_argument('--min-coverage', type=int, default=1, help="Minimum junction coverage in a single replicate.")
	parser_a.add_argument('--min-total-coverage', type=int, default=1, help="Minimum junction coverage taken over all replicates of a group.")

	# Prepare bam files for PSI calculation
	parser_b = subparsers.add_parser('prepare-bam', help="Extract the unique reads required for the PSI calculation.")
	parser_b.add_argument('-b', '--bam', required=True, help="BAM alignment file.")
	parser_b.add_argument('-f', '--file-handle', required=True, help="Unique file handle. The output files will be [work_dir]/uniq_reads.[file-handle].bam and [work_dir]/uniq_exonic_reads.[file-handle].bam.")
	parser_b.add_argument('-g', '--genome-fasta', required=True, help="Genome fasta corresponding to the one used for the BAM file creation.")
	parser_b.add_argument('--NPROC', type=stay_positive, default=4, help="Number of processes to use when running samtools.")

	parser_c = subparsers.add_parser('prepare-multi-bam', help="Prepare BAM files for all samples in the files.txt file.")
	parser_c.add_argument('--samples', required=True, help="Tab-separated file containing the group id (e.g. wt), absolute path to the folder containing the mapping bam and junction files, and sample id (e.g. wt-1).")
	parser_c.add_argument('-g', '--genome-fasta', required=True, help="Genome fasta corresponding to the one used for the BAM file creation.")
	parser_c.add_argument('--bam-filename', default="Aligned.sortedByCoord.out.bam", help="BAM filename (not path). Assumes that all BAM files have the same name, but different paths.")
	parser_c.add_argument('--NPROC', type=stay_positive, default=4, help="Number of processes to use when running samtools.")	

	args = parser.parse_args()

	# Create work-dir if not exists
	work_dir = args.work_dir if args.work_dir.endswith('/') else args.work_dir + '/'
	if not os.path.exists(work_dir): os.makedirs(work_dir)

	if args.command == "identify-exitrons":
		identify_exitrons(work_dir, args)
		log_settings(work_dir, args, 'w')

	elif args.command == "prepare-bam":
		prepare_bam_files(work_dir, args.bam, args.file_handle, args.genome_fasta, args.NPROC)
		log_settings(work_dir, args, 'a')

	elif args.command == "prepare-multi-bam":
		prepare_multi_bam(work_dir, args)
		log_settings(work_dir, args, 'a')