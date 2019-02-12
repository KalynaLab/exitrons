#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
exitron-pipeline.py => Exitron identification, quantification and comparison pipeline
"""

from __future__ import division
from natsort import natsorted
import numpy as np
import scipy.stats
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
		if not line.startswith('#'):
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

	# Extract unique reads
	uniq_bam = "{}uniq.{}.bam".format(work_dir, file_handle)
	cmd = "samtools view -@ {0} {1} | grep -w \"NH:i:1\" | samtools view -@ {0} -bT {2} - | samtools sort -o {3} -".format(NPROC, bam_file, genome_fasta, uniq_bam)
	subprocess.call(cmd, shell=True)

	# Index bam
	subprocess.call("samtools index {}".format(uniq_bam), shell=True)

def prepare_multi_bam(work_dir, args):
	""" Loop through the samples in the file and 
	prepare the bam files. """

	for line in open(args.samples):
		if not line.startswith('#'):
			group_id, path, sample_id = line.rstrip().split('\t')[:3]
			print sample_id
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

def mean_CI(x, confidence=0.95):
	n = len(x)
	m, sem = np.mean(x), scipy.stats.sem(x)
	h = sem * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
	return m, m-h, m+h

def calculate_PSI(work_dir, exitron_info, uniq_bam, file_handle):

	import re

	"""
		http://bioinformatics.cvr.ac.uk/blog/tag/cigar-string/
		CIGAR
		D => Deletion; the nucleotide is present in the reference but not in the read
		I => Insertion; the nucleotide is present in the read, but not in the reference
		M => Match; can be either an alignment match or mismatch. The nucleotide is present in the reference
		N => Skipped region; a region of nucleotides is not present in the read
		P => Padding; padded area in the read and not in the reference
		S => Soft clipping; the clipped nucleotides are present in the read

		Case1 (M/X/=):
			start at the specified mapping position, set counter to 1
			Add 1 to both the counts of the bases from that position and the counter.
			Move to the next position.
			Repeat this process till counter is the same as the number associated with the operator.
		Case2 (N/D):
			Move the specified mapping position by the number associated with the operator.
		Case3 (I/S/H/P):
			Do nothing
	"""

	# Extract the reads overlapping with the exitron +/- 10nt region
	rc, info = {}, {}
	for line in open(exitron_info):
		if not line.startswith('#'):
			exitron_id, t_id, gene_id, gene_name, EI_len, EIx3 = line.rstrip().split('\t')
			info[exitron_id] = { 't_id': t_id, 'gene_id': gene_id, 'gene_name': gene_name, 'EI_length': EI_len, 'EIx3': EIx3 }

			c, coord, strand = exitron_id.split(':')
			s, e = map(int, coord.split('-'))
			adjusted_e = e-1
			m = s + ((adjusted_e-s)/2)

			# The 20 nt ranges that should be fully covered by a read in
			# order to be counted.
			ARange = set(xrange(s-10, s+10))
			BRange = set(xrange(int(m)-10, int(m)+10))
			CRange = set(xrange(e-10, e+10))

			# Get the required read/junction gap signature
			N = "{}N".format(e-s) 

			# Extract the reads from the bam file
			subprocess.call("samtools view {} {}:{}-{} > {}tmp.sam".format(uniq_bam, c, s-10, adjusted_e+10, work_dir), shell=True)

			A, B, C, D = 0, 0, 0, 0
			EIFreq = { str(x): 0 for x in xrange(s, e) } # EI per position coverage
			for aln in open(work_dir+"tmp.sam"):
				pos = int(aln.split('\t')[3])
				cigar = aln.split('\t')[5]

				# Go through the read alignment
				current_pos = pos
				for m in re.finditer('(\d+)(\w)', cigar):
					alnLen, alnOperator = m.groups()
					if alnOperator in ['M', 'X', '=']:

						alnRange = set(xrange(current_pos, (current_pos + int(alnLen))))
						current_pos += int(alnLen)

						# Check if a read fully covers the A, B, and/or C regions
						if ARange <= alnRange: A += 1
						if BRange <= alnRange: B += 1
						if CRange <= alnRange: C += 1

						# Track the exitron per base coverage
						for x in alnRange:
							if str(x) in EIFreq:
								EIFreq[str(x)] += 1

					elif alnOperator in ['N', 'D']:
						if current_pos == s and str(alnLen)+alnOperator == N:
							D += 1
						current_pos += int(alnLen)

					else:
						pass

			# Store exitron coverage
			rc[exitron_id] = { 'A': A, 'B': B, 'C': C, 'D': D, 'cov': [ EIFreq[x] for x in EIFreq] }

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

	# Calculate PSI and some coverage stats/metrics
	with open("{}{}.psi".format(work_dir, file_handle), 'w') as fout:
		fout.write( "exitron_id\ttranscript_id\tgene_id\tgene_name\tEI_length\tEIx3\tA\tB\tC\tD\tPSI\tCOVERAGE_SCORE\tREAD_BALANCE\t5th/median/95th\tLOW_COV\n" )
		for EI in natsorted(rc):
			A, B, C, D = [ rc[EI][x] for x in ['A', 'B', 'C', 'D'] ]
			try: PSI = (np.mean([A, B, C]) / (np.mean([A, B, C]) + D)) * 100
			except ZeroDivisionError: PSI = 'nan'

			# Get the coverage metrics
			rs, rb = quality_score(A, B, C, D)
			#m, ci_low, ci_high = mean_CI(rc[EI]['cov'])
			#CI = "{}≤{}≤{}".format(int(ci_low), int(m), int(ci_high))
			ninetyfive = "{}/{}/{}".format(int(np.percentile(rc[EI]['cov'], 5)), int(np.median(rc[EI]['cov'])), int(np.percentile(rc[EI]['cov'], 95)))
			fout.write( "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(EI, info[EI]["t_id"], info[EI]["gene_id"], info[EI]["gene_name"], info[EI]["EI_length"], info[EI]["EIx3"], A, B, C, D, PSI, rs, rb, ninetyfive, min(rc[EI]['cov'])) )

	# Clean-up
	subprocess.call("rm -f {}tmp*".format(work_dir), shell=True)

def calculate_multi_PSI(work_dir, args):
	""" Loop through the samples in the files file and 
	calculate the PSI for all exitrons for each sample. """

	for line in open(args.samples):
		if not line.startswith('#'):
			group_id, path, sample_id = line.rstrip().split('\t')[:3]
			uniq = "{}uniq.{}.bam".format(args.bam_dir, sample_id)
			calculate_PSI(work_dir, args.exitrons_info, uniq, sample_id)

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

def paired_permutation_test(x, y, max_2k=100000):

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
	n = 0
	for i, signs in enumerate(product([1, -1], repeat=k)):
		if i < max_2k:
			new_avg_diff = np.nanmean([ diff[j]*signs[j] for j in xrange(k) ])
			if abs(new_avg_diff) >= abs(avg_diff): n += 1
		else: break
	return n / (2**k) if (2**k) < max_2k else n / max_2k

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

def parse_samples(samples_file, control_name, test_name, paired=True):
	
	tmp, pairs, control, test = { control_name: {}, test_name: {} }, {}, [], []
	fields = ["group_id", "path", "sample_id", "pair_id", "gender", "tumor_stage", "etnicity"]
	for line in open(samples_file):
		if not line.startswith('#'):
			cols = line.rstrip().split('\t')
			d = { fields[i]: cols[i] if i in xrange(len(cols)) else None for i in xrange(len(fields)) }
			tmp[d["group_id"]][d["sample_id"]] = d

			if paired:
				try: pairs[d["pair_id"]][d["group_id"]] = tmp[d["group_id"]][d["sample_id"]]
				except KeyError: pairs[d["pair_id"]] = { d["group_id"]: tmp[d["group_id"]][d["sample_id"]] }

	if paired:
		for p in pairs:
			control.append(pairs[p][control_name])
			test.append(pairs[p][test_name])

	else:
		control = [ tmp[control_name][x] for x in tmp[control_name] ]
		test = [ tmp[test_name][x] for x in tmp[test_name] ]

	return control, test

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

def parse_gene_TPM(f):
	""" Parse the gene TPM values as floats """

	tpm = {}
	for line in open(f):
		try: 
			gene_name, gene_length, eff_length, TPM, numReads = line.rstrip().split('\t')
			tpm[gene_name] = float(TPM)
		except ValueError: pass
	return tpm

"""def compare(work_dir, args):
	# Compare the exitrons of two groups

	import numpy as np

	#Methods: 
	#	1) Delta: Simply calculate the difference between the two groups
	#	2) Mann-Whitney U test
	#	2) Paired permutation test
	#	3a) Wilcoxon signed rank test
	#	3b) Paired wilcoxon signed rank test
	#	...
	
	# Get the PSI values for each sample
	psi_per_sample, info = {}, {}
	for line in open(args.samples):
		if not line.startswith('#'):
			group_id, path, sample_id = line.rstrip().split('\t')[:3]
			try: psi_per_sample[group_id][sample_id] = parse_PSI(args.psi_dir+sample_id+".psi")
			except KeyError: 
				psi_per_sample[group_id] = { sample_id: parse_PSI(args.psi_dir+sample_id+".psi") }
				info = parse_PSI(args.psi_dir+sample_id+".psi")

	if args.paired:
		pairs = {}
		for line in open(args.samples):
			if not line.startswith('#'):
				group_id, path, sample_id, pair_id = line.rstrip().split('\t')[:4]
				try: pairs[pair_id][group_id] = sample_id
				except KeyError: pairs[pair_id] = { group_id: sample_id }
		for p in pairs:
			print p, pairs[p]

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
					psi_vals[ei][g].append(psi_per_sample[g][pairs[p][g]][ei]["PSI"])

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
"""

def compare(work_dir, args):

	# Parse the samples files
	reference, test = parse_samples(args.samples, args.reference, args.test, args.paired)
	refNames = np.array([ r["sample_id"] for r in reference ])
	tstNames = np.array([ t["sample_id"] for t in test ])
	#print refNames
	#print tstNames

	# Get the PSI values in the same order
	refPSI = [ parse_PSI(args.psi_dir+r["sample_id"]+".psi") for r in reference ]
	tstPSI = [ parse_PSI(args.psi_dir+t["sample_id"]+".psi") for t in test ]

	# Parse the gene TPM
	if args.expr_filter:
		refTPM = [ parse_gene_TPM(r["path"]+args.gene_TPM_file) for r in reference ]
		tstTPM = [ parse_gene_TPM(t["path"]+args.gene_TPM_file) for t in test ]

	# Compare samples
	results = {}
	for ei in refPSI[0]:

		rPSI = np.array([ r[ei]["PSI"] for r in refPSI ])
		tPSI = np.array([ t[ei]["PSI"] for t in tstPSI ])

		# Select PSI values based on the gene TPM cut-off
		if args.expr_filter:

			eiGene = refPSI[0][ei]["gene_id"]
			#print ei, eiGene
			#print rPSI, tPSI

			rTPM = np.array([ r[eiGene] for r in refTPM ])
			tTPM = np.array([ t[eiGene] for t in tstTPM ])

			# Select only those pairs for which the TPM cut-off
			# is exceeded for both the reference and test
			if args.paired:
				fltr = map( lambda x, y: (x >= args.min_TPM) & (y >= args.min_TPM), rTPM, tTPM )
				rPSI = rPSI[fltr]
				tPSI = tPSI[fltr]
				usedSamples = '{}'.format(','.join([ x.replace(args.reference, '') for x in refNames[fltr] ]))
				nSamples = sum(fltr)
			
			# Filter reference and test individually
			else:
				rFltr = [ x >= args.min_TPM for x in rTPM ]
				tFltr = [ x >= args.min_TPM for x in tTPM ]
				rPSI = rPSI[rFltr]
				tPSI = tPSI[tFltr]
				usedSamples = '{};{}'.format(','.join(refNames[rFltr]), ','.join(tstNames[tFltr]))
				nSamples = '{};{}'.format(sum(rFltr), sum(tFltr))

			# Significance testing
			results[ei] = {
				'nSamples': nSamples,
				'usedSamples': usedSamples,
				'p-value': paired_permutation_test(rPSI, tPSI) if args.paired else monte_carlo_permutation_test(rPSI, tPSI),
				'meanRefPSI': np.mean(rPSI),
				'meanTstPSI': np.mean(tPSI),
				'dPSI': np.mean(tPSI) - np.mean(rPSI)
			}
			results[ei]['p-value'] = results[ei]['p-value'] if results[ei]['nSamples'] else float('nan')

	fResults = open("{}{}_{}.{}.diff".format(work_dir, args.reference, args.test, args.file_handle), 'w')
	fResults.write("#exitron_id\ttranscript_id\tgene_id\tgene_name\tEI_length_in_nt\tEIx3\t{}_mean\t{}_mean\tdiff\tpvalue\n".format(args.reference, args.test))

	fFilter = open("{}{}_{}.{}.fltr".format(work_dir, args.reference, args.test, args.file_handle), 'w')
	fFilter.write('#exitron_id\tN\tsamples\n')

	for ei in natsorted(results):
		fResults.write("{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.9f}\t{}\n".format(
			ei,
			refPSI[0][ei]["transcript_id"],
			refPSI[0][ei]["gene_id"],
			refPSI[0][ei]["gene_name"],
			refPSI[0][ei]["EI_length"],
			refPSI[0][ei]["EIx3"],
			results[ei]["meanRefPSI"],
			results[ei]["meanTstPSI"],
			results[ei]["dPSI"],
			results[ei]["p-value"],
			results[ei]["nSamples"])
		)
		fFilter.write("{}\t{}\t{}\n".format(ei, results[ei]["nSamples"], results[ei]["usedSamples"]))

	fFilter.close(), fResults.close()

if __name__ == '__main__':

	def stay_positive(val):
		if int(val) < 1:
			raise argparse.ArgumentTypeError("--NPROC must be at least 1.")
		return int(val)

	version = "0.4.0"
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
	parser_a.add_argument('--min-count', type=int, default=3, help="Minimum number of replicates a junction needs to occur in per group.")
	parser_a.add_argument('--min-coverage', type=int, default=1, help="Minimum junction coverage in a single replicate.")
	parser_a.add_argument('--min-total-coverage', type=int, default=1, help="Minimum junction coverage taken over all replicates of a group.")

	# Prepare bam files for PSI calculation
	parser_b = subparsers.add_parser('prepare-bam', help="Extract the unique reads required for the PSI calculation.")
	parser_b.add_argument('-b', '--bam', required=True, help="BAM alignment file.")
	parser_b.add_argument('-f', '--file-handle', required=True, help="Unique file handle. The output files will be [work_dir]/ujr.[file-handle].bam and [work_dir]/uer.[file-handle].bam.")
	parser_b.add_argument('-g', '--genome-fasta', required=True, help="Genome fasta corresponding to the one used for the BAM file creation.")
	parser_b.add_argument('--NPROC', type=stay_positive, default=4, help="Number of processes to use when running samtools.")

	parser_c = subparsers.add_parser('prepare-multi-bam', help="Prepare BAM files for all samples in the samples.txt file.")
	parser_c.add_argument('--samples', required=True, help="Tab-separated file containing the group id (e.g. wt), absolute path to the folder containing the mapping bam and junction files, and sample id (e.g. wt-1).")
	parser_c.add_argument('-g', '--genome-fasta', required=True, help="Genome fasta corresponding to the one used for the BAM file creation.")
	parser_c.add_argument('--bam-filename', default="Aligned.sortedByCoord.out.bam", help="BAM filename (not path). Assumes that all BAM files have the same name, but different paths.")
	parser_c.add_argument('--NPROC', type=stay_positive, default=4, help="Number of processes to use when running samtools.")	

	parser_d = subparsers.add_parser('calculate-PSI', help="Calculate the exitron PSI.")
	parser_d.add_argument('--exitrons-info', required=True, help="exitrons.info file (from identify-exitrons).")
	parser_d.add_argument('--uniq', required=True, help="Unique reads BAM file (from prepare-bam).")

	parser_e = subparsers.add_parser('calculate-multi-PSI', help="Calculate PSI for all samples in the samples.txt file.")
	parser_e.add_argument('--bam-dir', required=True, help="Path to the prepared bam files (from prepare-bam).")
	parser_e.add_argument('--exitrons-info', required=True, help="exitrons.info file (from identify-exitrons).")
	parser_e.add_argument('--samples', required=True, help="Tab-separated file containing the group id (e.g. wt), absolute path to the folder containing the mapping bam and junction files, and sample id (e.g. wt-1).")

	parser_f = subparsers.add_parser('compare', help="Compare exitrons of two groups.")
	parser_f.add_argument('--samples', required=True, help="Tab-separated file containing the group id (e.g. wt), absolute path to the folder containing the mapping bam and junction files, and sample id (e.g. wt-1).")
	parser_f.add_argument('--psi-dir', required=True, help="PSI files directory.")
	parser_f.add_argument('--file-handle', default="all", help="Additional file handle to be appended.")
	parser_f.add_argument('--reference', required=True, help="Group id (as listed in the samples file) of the test group.")
	parser_f.add_argument('--test', required=True, help="Group id (as listed in the samples file) of the test group.")
	parser_f.add_argument('--paired', action='store_true', help="Treat data as paired data. Requires the 4th column in the samples file to be the pair_id.")
	parser_f.add_argument('--min-TPM', type=float, default=1, help="Minimum gene TPM cut-off.")
	parser_f.add_argument('--expr-filter', action='store_true', help="Filter based on minimum gene TPM. It's assumed that the gene expression file is located in the path directory specified in the samples file.")
	parser_f.add_argument('--gene-TPM-file', default="quant.genes.sf", help="Name of the Salmon per gene TPM quantification file.")

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
		if all([ os.path.exists(args.exitrons_info), os.path.exists(args.uniq) ]):
			calculate_PSI(work_dir, args.exitrons_info, args.uniq, args.file_handle)
			log_settings(work_dir, args, 'a')
		else:
			sys.stderr.write("One (or multiple) input files could not be found.")

	elif args.command == "calculate-multi-PSI":
		calculate_multi_PSI(work_dir, args)
		log_settings(work_dir, args, 'a')

	elif args.command == "compare":
		compare(work_dir, args)