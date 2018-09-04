#!/usr/bin/env python

"""
exitron-pipeline.py => Exitron identification, quantification and comparison pipeline
"""

from __future__ import division
from natsort import natsorted
import numpy as np
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
					'depth': int(cols[6]) # Take depth as unique reads only
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
	infoOut, bedOut, eceOut = open(work_dir+"exitrons.info", 'w'), open(work_dir+"exitrons.bed", 'w'), open(work_dir+"exitron-containing-exons.bed", 'w')
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
		exitron_id = "{}:{}-{}:{}".format(j_chr, j_start, int(j_end)+1, j_strand)

		if not exitron_id in seen:
			infoOut.write( "{}\t{}\t{}\t{}\t{}\t{}\n".format(exitron_id, attr["transcript_id"], attr["gene_id"], attr["gene_name"], EI_length, EIx3) )
			bedOut.write( "{}\t{}\t{}\t{};{}\t1000\t{}\n".format(j_chr, j_start, int(j_end)+1, attr["gene_name"], exitron_id, j_strand) )
			eceOut.write( "{}\t{}\t{}\t{};{}\t1000\t{}\n".format(j_chr, cds_start, cds_end, attr["gene_name"], exitron_id, j_strand) )
			seen.add(exitron_id)

	infoOut.close(), bedOut.close(), eceOut.close()

	subprocess.call("rm -f "+work_dir+"junctions_GTF_map.tmp", shell=True)

def prepare_bam_files(work_dir, bam_file, file_handle, genome_fasta, NPROC=4):
	""" Extract the unique reads and unique exonic reads from the 
	supplied bam file, index and output to the working directory. """

	# Extract unique junctions reads (ujr)
	uniq_reads_bam = "{}ujr.{}.bam".format(work_dir, file_handle) 
	cmd = "samtools view -@ %d %s | grep -w \"NH:i:1\" | perl -n -e '@line=split(/\\t/,$_); if ($line[5]=~/N/){ print \"$_\"; }' > %s" % (NPROC, bam_file, uniq_reads_bam.replace('.bam', '.sam'))
	subprocess.call(cmd, shell=True)

	# Convert sam to bam
	subprocess.call("samtools view -@ {0} -bT {1} {2} | samtools sort -o {3} -".format(NPROC, genome_fasta, uniq_reads_bam.replace('.bam', '.sam'), uniq_reads_bam), shell=True)

	# Index bam
	subprocess.call("samtools index {}".format(uniq_reads_bam), shell=True)

	# Remove sam
	subprocess.call("rm -f {}".format(uniq_reads_bam.replace('.bam', '.sam')), shell=True)

	# Extract unique exonic reads (uer)
	uniq_exon_reads_bam = '{}uer.{}.bam'.format(work_dir, file_handle)
	cmd = "samtools view -@ %d %s | grep -w \"NH:i:1\" | perl -n -e '@line=split(/\\t/,$_); if ($line[5]!~/N/){ print \"$_\"; }' > %s" % (NPROC, bam_file, uniq_exon_reads_bam.replace('.bam', '.sam'))
	subprocess.call(cmd, shell=True)

	# Convert sam to bam
	subprocess.call("samtools view -@ {0} -bT {1} {2} | samtools sort -o {3} -".format(NPROC, genome_fasta, uniq_exon_reads_bam.replace(".bam", ".sam"), uniq_exon_reads_bam), shell=True)

	# Index bam
	subprocess.call("samtools index {}".format(uniq_exon_reads_bam), shell=True)

	# Remove sam
	subprocess.call("rm -f {}".format(uniq_exon_reads_bam.replace('.bam', '.sam')), shell=True)

def prepare_multi_bam(work_dir, args):
	""" Loop through the samples in the file and 
	prepare the bam files. """

	for line in open(args.samples):
		group_id, path, sample_id = line.rstrip().split('\t')[:3]
		prepare_bam_files(work_dir, path+args.bam_filename, sample_id, args.genome_fasta, args.NPROC)

def quality_score(A, B, C, D):
	""" Score the exitron based on the number of reads """

	from scipy.stats import binom_test

	if D >= 100 or A >= 100 and C >= 100: RS = "SOK"
	elif D >= 20 or any([ A >= 20 and C >= 15, A >= 15 and C >= 20 ]): RS = "OK"
	elif D >= 15 or any([ A >= 15 and C >= 10, A >= 10 and C >= 15 ]): RS = "LOW"
	elif D >= 10 or any([ A >= 10 and C >= 5, A >= 5 and C >= 10 ]): RS = "VLOW"
	else: RS = "N"

	''' 
		Filter EI events based on (median) read coverage and read balance. 
		The goal of the binomial test (balance) was to exclude events in which there is a 
		high imbalance in read counts among the two exon-intron junctions and the intron 
		body sequence. Such imbalances can arise from neighboring alternative 5' and/or 3'
		splice sites or overlapping genes, confound PSI estimates, and lead to the false 
		detection of EI.

		p-value (binomial{ 	M = min(A, B, C),
							N = min(A, B, C) + max(A, B, C),
							P = 1/3.5,
							alternative = lower })

		M = number of successes
		N = number of trials
		P = probability of success

		scipy.stats.binom_test(M, N, 1/3.5, alternative="less") 
	'''

	if ( np.median([ A, B, C ]) + D > 10 ) and \
	binom_test(
		min(A, B, C),
		min(A, B, C) + max(A, B, C),
		1/3.5,
		alternative="less"
	) >= 0.05: balance = "OKAY"
	else: balance = "NOT_OKAY"

	return RS, balance

def calculate_PSI(work_dir, exitron_info, ujr_bam, uer_bam, file_handle):
	''' Calculate exitron PSI values, based on the coverage of the 
	unique exonic reads. 
						
						A		   B		 C
					   ----		  ----		----
	EEEEEEEEEEEEEEEEEEEEEXXXXXXXXXXXXXXXXXXXXXEEEEEEEEEEEEEEEEEEEEE
					   --	                  --
						 \                   / 
						  \        D        / 							
                           -----------------
                           	
    E = exon, X = exitron
    A = Reads aligning from -10 to +10 around exitron 5'SS
    B = Reads aligning from -10 to +10 around exitron middle point
    C = Reads aligning from -10 to +10 around exitron 3'SS
    D = Reads supporting exitron splicing
    
	PSI = ( mean(A, B, C) / mean(A, B, C) + D) * 100 '''

	import re

	# Get junction support
	rc, info = {}, {}
	track_seen, ABCOut = set(), open(work_dir+"tmp_coverageBed_input.bed", 'w')
	for line in open(exitron_info):
		if not line.startswith('#'):
			exitron_id, t_id, gene_id, gene_name, EI_len, EIx3 = line.rstrip().split('\t')
			c, coord, strand = exitron_id.split(':')
			s, e = map(int, coord.split('-'))
			info[exitron_id] = { 't_id': t_id, 'gene_id': gene_id, 'gene_name': gene_name, 'EI_length': EI_len, 'EIx3': EIx3 }

			if exitron_id not in track_seen: 

				adjusted_e = e-1

				### Unique junction support ###
				subprocess.call("samtools view {} {}:{}-{} > {}tmp.sam".format(ujr_bam, c, s, adjusted_e, work_dir), shell=True)
				N = "{}N".format(e-s) # Get the required read/junction gap signature

				uniq_junction_count = 0
				for aln in open(work_dir+"tmp.sam"):
					qname = aln.split('\t')[0]
					pos = int(aln.split('\t')[3])
					cigar = aln.split('\t')[5]
					try: start = (pos + int(re.search('^([0-9]+)M', cigar).group(1)))
					except AttributeError: start = None

					# Check if the junction is at the correct position 
					# and if the junction size is correct
					if N in cigar and start == int(s): uniq_junction_count += 1

				rc[exitron_id] = { 'A': 0, 'B': 0, 'C': 0, 'D': uniq_junction_count }

				### Bed file generation for A, B, C read support ###
				middle_point = int(s + ((adjusted_e-s)/2))
				ABCOut.write( "{}\t{}\t{}\t{}_A\n".format(c, s-10, s+10, exitron_id) )
				ABCOut.write( "{}\t{}\t{}\t{}_B\n".format(c, middle_point-10, middle_point+10, exitron_id) )
				ABCOut.write( "{}\t{}\t{}\t{}_C\n".format(c, e-10, e+10, exitron_id) )

				# Track seen
				track_seen.add(exitron_id)

	ABCOut.close()

	### Get the unique A, B, C read support ###
	subprocess.call("samtools view -H {} | grep -P \"@SQ\tSN:\" | sed 's/@SQ\tSN://' | sed 's/\tSN://' | sed 's/\tLN:/\t/' > {}tmp_genome.txt".format(uer_bam, work_dir), shell=True)
	subprocess.call("sort -k1,1V -k2,2n {0}tmp_coverageBed_input.bed > {0}tmp_coverageBed_input.sorted.bed".format(work_dir), shell=True)
	subprocess.call("bedtools coverage -sorted -counts -g {0}tmp_genome.txt -a {0}tmp_coverageBed_input.sorted.bed -b {1} > {0}tmp_coverageBed_output.bed".format(work_dir, uer_bam), shell=True)
	for line in open(work_dir+"tmp_coverageBed_output.bed"):
		c, start, end, locus, coverage = line.rstrip().split('\t')
		exitron_id, letter = locus.split('_')
		rc[exitron_id][letter] = int(coverage)

	### Calculate PSI ###
	with open("{}{}.psi".format(work_dir, file_handle), 'w') as fout:
		fout.write( "exitron_id\ttranscript_id\tgene_id\tgene_name\tEI_length\tEIx3\tA\tB\tC\tD\tPSI\tCOVERAGE_SCORE\tREAD_BALANCE\n" )
		for EI in natsorted(rc):
			A, B, C, D = [ rc[EI][x] for x in ['A', 'B', 'C', 'D'] ]
			try: PSI = (np.mean([A, B, C]) / (np.mean([A, B, C]) + D)) * 100
			except ZeroDivisionError: PSI = 'nan'

			# Get the read score and read balance labels
			rs, rb = quality_score(A, B, C, D)
			fout.write( "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(EI, info[EI]["t_id"], info[EI]["gene_id"], info[EI]["gene_name"], info[EI]["EI_length"], info[EI]["EIx3"], A, B, C, D, PSI, rs, rb) )

	# Clean-up
	subprocess.call("rm -f {}tmp*".format(work_dir), shell=True)

def calculate_multi_PSI(work_dir, args):
	""" Loop through the samples in the files file and 
	calculate the PSI for all exitrons for each sample. """

	for line in open(args.samples):
		group_id, path, sample_id = line.rstrip().split('\t')[:3]
		ujr = "{}ujr.{}.bam".format(args.bam_dir, sample_id)
		uer = "{}uer.{}.bam".format(args.bam_dir, sample_id)
		calculate_PSI(work_dir, args.exitrons_info, ujr, uer, sample_id)

def monte_carlo_permutation_test(x, y, nmc=10000):

	import math

	# Calculate the unique number of combinations
	# without replacement
	n = len(x) + len(y)
	r = len(x)
	N = math.factorial(n) / math.factorial(r) / math.factorial(n-r)

	# Count the number of times the difference of the mean
	# of the randomly permutated samples is bigger than the
	# observed difference of the mean for x and y
	k, diff = 0.0, abs(np.nanmean(x) - np.nanmean(y))
	xy = np.concatenate([x, y])

	if N > nmc: # Do nmc random permutations
		for i in xrange(nmc):
			xy = np.concatenate([x, y])
			np.random.shuffle(xy)
			k += diff < abs(np.nanmean(xy[:r]) - np.nanmean(xy[r:]))
		return k / nmc

	else: # Check all permutations
		from itertools import combinations

		xy = np.concatenate([x, y])
		xy_indx = xrange(n)
		for comb in combinations(xy_indx, r):
			new_x = [ xy[i] for i in comb ]
			new_y = [ xy[i] for i in set(xy_indx).difference(set(comb)) ]
			k += diff < abs(np.nanmean(new_x) - np.nanmean(new_y))
		return k / N

def paired_permutation_test(x, y, max_2k=10000):

	"""
		https://arxiv.org/pdf/1603.00214.pdf -> Haven't read it, but is about permuting incomplete paired data
		http://axon.cs.byu.edu/Dan/478/assignments/permutation_test.php

		1) Obtain a set of k pairs of estimates {(a1, b1), (a2, b2), ..., (ak, bk)} for (M1, M2)
		2) Calculate the average difference, diff = sum([ ai-bi for i in k ]) / k
		3) Let n = 0
		4) For each possible permutation of swapping ai and bi (swapping labels)
			a) Calculate the average difference, new_diff = sum([ ai-bi for i in k ]) / k
			b) if abs(new_diff) >= abs(diff) => n += 1
		5) Report p = n/2**k
	"""

	from itertools import product

	# Calculate the differences and average difference
	k = len(x)
	diff = [ x[i]-y[i] for i in xrange(k) ]
	avg_diff = np.nanmean(diff)

	# Go through the permutations
	n = 0.0
	for i, signs in enumerate(product([1, -1], repeat=k)):
		if i < max_2k:
			new_avg_diff = np.nanmean([ diff[j]*signs[j] for j in xrange(k) ])
			if abs(new_avg_diff) >= abs(avg_diff): n += 1
		else: break
	return n / (2**k)

def CI_difference_between_means(x, y, confidence=0.95):
    """ Calculate the confidence interval for the difference
    of two means. 

    M1 - M2 +/- (tCL)(SM1-M2)
    """

    import math
    import scipy.stats as stats

    # http://onlinestatbook.com/2/estimation/difference_means.html
    M1, M2 = np.mean(x), np.mean(y)
    s1, s2 = np.var(x, ddof=1), np.var(y, ddof=1)
    diff = M2-M1

    # Degrees of freedom and finding t
    df = (len(x)-1) + (len(y)-1)
    t = stats.t.ppf((1+confidence)/2, df)

    if len(x) == len(y): # Equal sample sizes

        # Compute the estimate of the standard error
        # of the difference between means
        MSE = (s1+s2) / 2
        SM1M2 = math.sqrt( (2*MSE) / len(x) )

    else: # Unequal sample sizes

        # Compute the sum of squares error (SSE)
        SSE = sum([ (X-M1)**2 for X in x ]) + sum([ (X-M2)**2 for X in y ])
        MSE = SSE/df

        harmonic_mean = 2 / (1/len(x) + 1/len(y))
        SM1M2 = math.sqrt( (2*MSE) / harmonic_mean )

    return diff, diff - t*SM1M2, diff + t*SM1M2

def parse_PSI(f, readScore = ['VLOW', 'LOW', 'OK', 'SOK'], readBalance = ['OKAY']):
	""" Parse the PSI values. If an exitron's read score and 
	read balance are not in the readScore or readBalance arrays, 
	the PSI value will be parsed as nan. """

	psi = {}
	for line in open(f):
		if not line.startswith('exitron_id'):
			cols = line.rstrip().split('\t')
			psi[cols[0]] = {
				'transcript_id': cols[1],
				'gene_id': cols[2],
				'gene_name': cols[3],
				'EI_length': cols[4],
				'EIx3': cols[5],
				'PSI': float(cols[10]) if cols[10] != "NA" else float('nan'),
				'QS': all([ cols[11] in readScore, cols[12] in readBalance ])
			}
	return psi

def compare(work_dir, args):
	""" Compare the exitrons of two groups """

	import numpy as np

	"""
		Methods: 
			1) Delta: Simply calculate the difference between the two groups
			2) Mann-Whitney U test
			2) Paired permutation test
			3a) Wilcoxon signed rank test
			3b) Paired wilcoxon signed rank test
			...
	"""
	
	# Get the PSI values for each sample
	psi_per_sample, info = {}, {}
	for line in open(args.samples):
		group_id, path, sample_id = line.rstrip().split('\t')[:3]
		try: psi_per_sample[group_id][sample_id] = parse_PSI(args.psi_dir+sample_id+".psi")
		except KeyError: 
			psi_per_sample[group_id] = { sample_id: parse_PSI(args.psi_dir+sample_id+".psi") }
			info = parse_PSI(args.psi_dir+sample_id+".psi")

	if args.paired:
		pairs = {}
		for line in open(args.samples):
			group_id, path, sample_id, pair_id = line.rstrip().split('\t')[:4]
			try: pairs[pair_id][group_id] = sample_id
			except KeyError: pairs[pair_id] = { group_id: sample_id }

	# Get the PSI values per exitron, split up
	# for the different groups
	psi_vals, QS = {}, {}
	groups = [args.reference, args.test]
	some_sample = [ x for x in psi_per_sample[groups[0]] ][0]
	for ei in psi_per_sample[args.reference][some_sample]:

		if args.paired:
			psi_vals[ei] = { g: [] for g in groups }
			for p in pairs:
				for g in pairs[p]:
					psi_vals[ei][pairs[p][g]].append(psi_per_samples[g][pairs[p][g]][ei]["PSI"])

		else: 
			psi_vals[ei] = { g: [ psi_per_sample[g][r][ei]["PSI"] for r in psi_per_sample[g] ] for g in groups }
		
		QS[ei] = "A{};B{}".format(*[ sum([ psi_per_sample[g][r][ei]["QS"] for r in psi_per_sample[g] ]) for g in groups ])

	# TBD: Significance test selection
	results = {}
	for ei in psi_vals:
		results[ei] = { 'pvalue': paired_permutation_test(psi_vals[ei][args.reference], psi_vals[ei][args.test]) if args.paired else monte_carlo_permutation_test(psi_vals[ei][args.reference], psi_vals[ei][args.test]) }

	# FDR

	# For now I'm just going to output the sample means,
	# the difference of the means + the 95% CI for the
	# difference between the means
	with open("{}{}_{}.diff".format(work_dir, args.reference, args.test), 'w') as fout:
		fout.write("#exitron_id\ttranscript_id\tgene_id\tgene_name\tEI_length_in_nt\tEIx3\t{}_mean\t{}_mean\tdiff\tpvalue\n".format(args.reference, args.test))
		for ei in natsorted(psi_vals):
			ref_psi = psi_vals[ei][args.reference]
			test_psi = psi_vals[ei][args.test]
			
			fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.9f}\t{}\n".format(
				ei,
				info[ei]["transcript_id"],
				info[ei]["gene_id"],
				info[ei]["gene_name"],
				info[ei]["EI_length"],
				info[ei]["EIx3"],
				np.nanmean(ref_psi),
				np.nanmean(test_psi),
				np.nanmean(test_psi)-np.nanmean(ref_psi),
				results[ei]['pvalue'], 
				QS[ei])
			)

if __name__ == '__main__':

	def stay_positive(val):
		if int(val) < 1:
			raise argparse.ArgumentTypeError("--NPROC must be at least 1.")
		return val

	version = "0.2.6"
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

	parser_c = subparsers.add_parser('prepare-multi-bam', help="Prepare BAM files for all samples in the samples.txt file.")
	parser_c.add_argument('--samples', required=True, help="Tab-separated file containing the group id (e.g. wt), absolute path to the folder containing the mapping bam and junction files, and sample id (e.g. wt-1).")
	parser_c.add_argument('-g', '--genome-fasta', required=True, help="Genome fasta corresponding to the one used for the BAM file creation.")
	parser_c.add_argument('--bam-filename', default="Aligned.sortedByCoord.out.bam", help="BAM filename (not path). Assumes that all BAM files have the same name, but different paths.")
	parser_c.add_argument('--NPROC', type=stay_positive, default=4, help="Number of processes to use when running samtools.")	

	parser_d = subparsers.add_parser('calculate-PSI', help="Calculate the exitron PSI.")
	parser_d.add_argument('--exitrons-info', required=True, help="exitrons.info file (from identify-exitrons).")
	parser_d.add_argument('--ujr', required=True, help="Unique junction reads BAM file (from prepare-bam).")
	parser_d.add_argument('--uer', required=True, help="Unique exonic reads BAM file (from prepare-bam).")
	parser_d.add_argument('-f', '--file-handle', required=True, help="Unique file handle. The output file will be [work_dir]/[file-handle].psi.")

	parser_e = subparsers.add_parser('calculate-multi-PSI', help="Calculate PSI for all samples in the samples.txt file.")
	parser_e.add_argument('--bam-dir', required=True, help="Path to the prepared bam files (from prepare-bam).")
	parser_e.add_argument('--exitrons-info', required=True, help="exitrons.info file (from identify-exitrons).")
	parser_e.add_argument('--samples', required=True, help="Tab-separated file containing the group id (e.g. wt), absolute path to the folder containing the mapping bam and junction files, and sample id (e.g. wt-1).")

	parser_f = subparsers.add_parser('compare', help="Compare exitrons of two groups.")
	parser_f.add_argument('--samples', required=True, help="Tab-separated file containing the group id (e.g. wt), absolute path to the folder containing the mapping bam and junction files, and sample id (e.g. wt-1).")
	parser_f.add_argument('--psi-dir', required=True, help="PSI files directory.")
	parser_f.add_argument('--reference', required=True, help="Group id (as listed in the samples file) of the test group.")
	parser_f.add_argument('--test', required=True, help="Group id (as listed in the samples file) of the test group.")
	parser_f.add_argument('--paired', action='store_true', help="Treat data as paired data. Requires the 4th column in the samples file to be the pair_id.")

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

	elif args.command == "calculate-PSI":
		if all([ os.path.exists(args.exitrons_info), os.path.exists(args.ujr), os.path.exists(args.uer) ]):
			calculate_PSI(work_dir, args.exitrons_info, args.ujr, args.uer)
			log_settings(work_dir, args, 'a')
		else:
			sys.stderr.write("One (or multiple) input files could not be found.")

	elif args.command == "calculate-multi-PSI":
		calculate_multi_PSI(work_dir, args)
		log_settings(work_dir, args, 'a')

	elif args.command == "compare":
		compare(work_dir, args)