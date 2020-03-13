#!/usr/bin/env python

""" Extract basic features (splice sites, GC%, splice site scores)
from all exitrons and also detect disordered regions and protein
domains for the EIx3's. """

import os
import uuid
import argparse
import subprocess
from natsort import natsorted
from math import log

iupred_cmd = "iupred2a"
vsl2b_cmd = "java -jar ~/.local/bin/VSL2.jar"

def yield_fasta(f):
	''' Simple fasta parser that yield's the identifier and sequence of each record '''

	class SeqRecord:
		def __init__(self, seq_id, seq):
			self.id = seq_id
			self.seq = seq

	seq = ''
	for line in open(f):
		if line.startswith('>'):
			if len(seq):
				yield(SeqRecord(identifier, seq))
				seq = ''
			identifier = line[1:].rstrip()
		else: seq += line.rstrip()
	yield(SeqRecord(identifier, seq))

def exitron2fasta(work_dir, exitron_bed, genome_fasta, output_file, upstream_nt=0, downstream_nt=0):
	""" Extract the exitron fasta sequence, with going -x nt and +y nt
	in the neighboring exons. """

	# Create a temporary bed file w/ adjustment for
	# the upstream_nt and downstream_nt parameters
	tmpBed = uuid.uuid4().hex + ".bed"
	while os.path.isfile(tmpBed):
		tmpBed = uuid.uuid4().hex + ".bed"

	with open(tmpBed, 'w') as fout:
		for line in open(exitron_bed):
			c, s, e, bed_id, score, strand = line.rstrip().split('\t')
			fout.write( "{}\t{}\t{}\t{}\t{}\t{}\n".format(c, (int(s)-1)-upstream_nt, int(e)+downstream_nt, bed_id, score, strand) )

	# Extract the fasta sequences
	subprocess.call("bedtools getfasta -s -name -fi {0} -bed {1} | sed 's/\(::\|(\).*//g' > {2}{3}".format(genome_fasta, tmpBed, work_dir, output_file), shell=True)

	# Clean-up temporary bed file
	os.remove(tmpBed)

### Get exitron relative protein sequence ###
def get_relative_exitron_position(work_dir, exitrons_info, genome_fasta, ccds_gtf, reference_proteins):

	# Get the fasta sequence of the EIx3's
	ei = {}
	with open(work_dir+"exitrons_EIx3.bed", 'w') as fout:
		for line in open(exitrons_info):
			if not line.startswith('#'):
				exitron_id, t_id, g_id, gene_name, EI_len, EIx3 = line.rstrip().split('\t')
				c, coord, strand = exitron_id.split(':')
				s, e = map(int, coord.split('-'))

				if EIx3 == "yes":
					fout.write( "{}\t{}\t{}\t{}\t1000\t{}\n".format(c, s, e, exitron_id, strand) )
					ei[exitron_id] = t_id

	exitron2fasta(work_dir, work_dir+"exitrons_EIx3.bed", genome_fasta, "exitrons_EIx3.fa")

	# Get the reference transcript protein sequences
	ref_trs = set([ ei[x] for x in ei ])
	ref_prot = {}
	for record in yield_fasta(reference_proteins):
		try: t_id = record.id.split('|')[1]
		except IndexError: t_id = record.id.split(' ')[0]
		if t_id in ref_trs:
			ref_prot[t_id] = record.seq

	# Parse the exitron reference transcripts coordinates
	gtf = {}
	for line in open(ccds_gtf):
		if not line.startswith('#'):
			c, source, feature, start, end, score, strand, frame, attributes = line.rstrip().split('\t')
			attr = {}
			for a in attributes.split(';'):
				if len(a):
					attr_name, attr_value = list(filter(None, a.split(' ')))
					attr[attr_name.strip()] = attr_value.replace('\"', '')

			if feature == "CDS" and attr["transcript_id"] in ref_trs:
				try: gtf[attr["transcript_id"]].append([ int(start), int(end) ])
				except KeyError: gtf[attr["transcript_id"]] = [[ int(start), int(end) ]]

	# Get the relative exitron position on
	# the reference isoforms
	refProtOut = open(work_dir+"EIx3_ref_proteins.fa", 'w')
	eiPosOut = open(work_dir+"EIx3_ref_positions.txt", 'w')

	for record in yield_fasta(work_dir+"exitrons_EIx3.fa"):

		#t_id = ei[record.id]
		ccds = gtf[ei[record.id]]
		ei_len = len(record.seq)

		# Get the EIx3 AA position
		c, coord, strand = record.id.split(':')
		ei_start, ei_end = map(int, coord.split('-'))

		track = 0
		if strand == '+':
			for x, y in ccds:
				if x < ei_start < y and x < ei_end < y:
					leftover = (ei_start-x) % 3
					ei_from = track + (((ei_start-x)-leftover) / 3) + 1
					ei_to = ei_from + (ei_len / 3) - 1
				else:
					leftover = (y-x) % 3
					track += (((y-x)-leftover) / 3) + 1

		elif strand == '-':
			for x, y in natsorted(ccds, reverse=True):
				if x < ei_start < y and x < ei_end < y:
					leftover = (y-ei_end) % 3
					ei_from = track + (((y-ei_end)-leftover) / 3) + 1
					ei_to = ei_from + (ei_len / 3) - 1
				else:
					leftover = (y-x) % 3
					track += (((y-x)-leftover) / 3) + 1

		try:
			refProtOut.write( ">{};{}\n{}\n".format(ei[record.id], record.id, ref_prot[ei[record.id]]) )
			eiPosOut.write( "{}\t{:.0f}\t{:.0f}\n".format(record.id, ei_from, ei_to) )
		except KeyError: pass

	eiPosOut.close()
	refProtOut.close()

### Get protein domains ###
def get_protein_domains(work_dir, pfam_db, reference_proteins, exitron_reference_pos, max_e_value=0.001):

	import re
	from collections import Counter

	subprocess.call("hmmsearch --tblout {0}/exitron_hmmsearch_Pfam.txt -E 1e-2 --cpu 4 {1} {2} > {0}/exitron_hmmsearch_Pfam_alignments.txt".format(work_dir, pfam_db, reference_proteins), shell=True)

	# Parse exitron reference positions
	EI_pos = {}
	for line in open(exitron_reference_pos):
		exitron_id, ei_from, ei_to = line.rstrip().split('\t')
		EI_pos[exitron_id] = { "ei_from": int(ei_from), "ei_to": int(ei_to) }

	# Parse exitron overlapping domains
	parse, lineCount = False, 0
	domains, domainCount, domainDesc, onePerGene = {}, Counter(), {}, {}
	domainRE = re.compile("Query:\s+(.*?)\s+\[M=(\d+)\]")
	for line in open(work_dir+"exitron_hmmsearch_Pfam_alignments.txt"):

		# Get the domain name and size
		if line.startswith("Query:"):
			domain, domain_size = re.search(domainRE, line.rstrip()).groups()

		# Get the domain description
		if line.startswith("Description:"):
			desc = line.rstrip().split("Description: ")[1]
			domainDesc[domain] = desc

		# Get the matching exitron and match position
		if line.startswith(">>"):
			ref_t, exitron_id = line.rstrip().split(">> ")[1].split(';')
			parse = True

		if parse and line.rstrip() == "":
			parse, lineCount = False, 0

		elif parse and lineCount >= 3:
			domain_match = re.sub('\s+', ';', line.strip()).split(';')
			e_value, hmm_from, hmm_to, hmm_code, ali_from, ali_to = domain_match[5:11]
			env_from, env_to, env_code = domain_match[12:15]

			if float(e_value) <= max_e_value:

				# See if the exitron overlaps with the detected domain
				ei_from = EI_pos[exitron_id]["ei_from"]
				ei_to = EI_pos[exitron_id]["ei_to"]
				ei_range = set(range(ei_from, ei_to+1))
				ali_range = set(range(int(ali_from), int(ali_to)+1))

				# Domain/exitron overlap label and percentage
				domain_ei_overlap = ali_range.intersection(ei_range)
				if ei_range.issubset(ali_range):
					overlap_label = "exitron_inside_domain"
					domain_ei_perc = 100
				elif ali_range.issubset(ei_range):
					overlap_label = "domain_inside_exitron"
					domain_ei_perc = 100
				else:
					overlap_label = "part_domain_exitron"
					domain_ei_perc = (len(domain_ei_overlap) / int(domain_size)) * 100

				if len(domain_ei_overlap) and any([ overlap_label == "exitron_inside_domain", domain_ei_perc >= 20 ]):
					output_line = "{}\t{}\t{}\t{}\t{}-{}\t{}\t{}-{}\n".format(ref_t, exitron_id, domain, e_value, env_from, env_to, overlap_label, ei_from, ei_to)
					try: domains[exitron_id].append(output_line)
					except KeyError: domains[exitron_id] = [output_line]

					try: onePerGene[ref_t].add(domain)
					except KeyError: onePerGene[ref_t] = set([ domain ])

					domainCount[domain] += 1

			lineCount += 1

		elif parse:
			lineCount += 1

	# Output the Pfam domain matches
	with open(work_dir+"exitron_Pfam_domains.txt", 'w') as fout:
		fout.write("Reference transcript\tExitron ID\tPfam domain\tE-value\tDomain location protein\tExitron location domain\tExitron location protein\n")
		for exitron_id in natsorted(domains):
			for line in domains[exitron_id]:
				fout.write(line)

	# Sum up onePerGene domain counts
	opg = Counter()
	for g in onePerGene:
		for x in onePerGene[g]:
			opg[x] += 1

	# Output the Pfam domain counts
	with open(work_dir+"exitron_Pfam_domain_counts.txt", 'w') as fout:
		fout.write("Pfam domain\tPfam domain counts\tOne domain per gene counts\tPfam domain description\n")
		for d in domainCount.most_common():
			fout.write( "{}\t{}\t{}\t{}\n".format(d[0], d[1], opg[d[0]], domainDesc[d[0]]) )

	subprocess.call("rm {}exitron_hmmsearch*".format(work_dir), shell=True)

### Get disordered regions ###
def parse_iupred(f, ei_from, ei_to, col=2):

	# Parse disordered regions (at least 20 consecutive AA
	# with values above 0.5)

	disordered, region, start_pos = [], [], None
	lineCount = 0
	for line in open(f):
		if not line.startswith('#'):
			lineCount += 1
			pos, AA = line.split('\t')[:2]
			score = float(line.rstrip().split('\t')[col])

			if score > 0.5:
				region.append(AA)
				if start_pos == None: start_pos = int(pos)

			elif score <= 0.5 and len(region) >= 20:
				disordered.append([start_pos, int(pos)-1])
				region, start_pos = [], None

			elif score <= 0.5:
				region, start_pos = [], None

	if len(region) >= 20:
		disordered.append([start_pos, int(pos)])

	# Calculate the disordered percentage
	sumDisAA = sum([ abs(y-x) for x, y in disordered ])
	percDis = (sumDisAA/lineCount) * 100

	return disordered, percDis, any([ len(set(range(ei_from, ei_to)).intersection(set(range(x, y)))) >= 20 for x, y in disordered ])

def get_disordered_regions(work_dir, reference_proteins, exitron_reference_pos):

	# Parse exitron reference positions
	EI_pos = {}
	for line in open(exitron_reference_pos):
		exitron_id, ei_from, ei_to = line.rstrip().split('\t')
		EI_pos[exitron_id] = { "ei_from": int(ei_from), "ei_to": int(ei_to) }

	# Iterate over the exitron fasta sequences
	iupred, anchor, vsl2b = {}, {}, {}
	for record in yield_fasta(reference_proteins):

		ei_from, ei_to = [ EI_pos[record.id.split(';')[1]][x] for x in ['ei_from', 'ei_to'] ]

		# Write individual sequence to file, since
		# IUPred2A only takes one sequence at a time
		with open(work_dir+'tmp.fa', 'w') as fout:
			fout.write( '{}\n'.format(record.seq) )

		### IUPred2A & ANCHOR2 ###
		# Run IUPred2A in short mode
		subprocess.call("{0} -a {1}tmp.fa short > {1}tmp_iupred_output.txt".format(iupred_cmd, work_dir), shell=True)

		# IUPred2A
		disRegions, disPerc, eiOverlap = parse_iupred(work_dir+"tmp_iupred_output.txt", ei_from, ei_to)
		if eiOverlap:
			iupred[record.id] = '{}\t{:.2f}%'.format(','.join([ '{}-{}'.format(x, y) for x, y in disRegions]), disPerc)

		# ANCHOR2
		bindRegions, bindPerc, eiOverlap = parse_iupred(work_dir+"tmp_iupred_output.txt", ei_from, ei_to, 3)
		if eiOverlap:
			anchor[record.id] = '{}\t{:.2f}%'.format(','.join([ '{}-{}'.format(x, y) for x, y in bindRegions]), bindPerc)

		### VSL2B ###
		# Run VSL in VSL2B mode
		subprocess.call("{0} -s:{1}tmp.fa > {1}tmp_vsl2b_output.txt".format(vsl2b_cmd, work_dir), shell=True)

		# Parse disordered regions (at least 20
		# consecutive AA with values above 0.5)
		disRegions, region, start_pos, parse = [], [], None, False
		for line in open(work_dir+"tmp_vsl2b_output.txt"):
			if line.startswith('---'):
				parse = True

			elif parse and line.startswith('==='):
				parse = False

				# End of file check
				if len(region) >= 20:
					disRegions.append([start_pos, int(pos)])

			elif parse:

				pos, AA, score, disorder = line.rstrip().split('\t')
				if float(score) > 0.5:
					region.append(AA)
					if start_pos == None: start_pos = int(pos)

				elif float(score) <= 0.5 and len(region) >= 20:
					disRegions.append([start_pos, int(pos)-1])
					region, start_pos = [], None

				elif float(score) <= 0.5:
					region, start_pos = [], None

		sumDisAA = sum([ abs(y-x) for x, y in disRegions ])
		percDis = (sumDisAA/len(record.seq)) * 100
		if any([ len(set(range(ei_from, ei_to)).intersection(set(range(x, y)))) >= 20 for x, y in disRegions ]):
			vsl2b[record.id] = '{}\t{:.2f}%'.format(','.join([ '{}-{}'.format(x, y) for x, y in disRegions ]), percDis)

	# Remove IUPred2A and VSL2B temporary files
	subprocess.call("rm {}tmp*".format(work_dir), shell=True)

	# Output all disordered regions
	with open(work_dir+"exitron_disordered_regions.txt", 'w') as fout:
		fout.write("Reference Transcript\tExitron ID\tExitron location protein\tIUPred2A disordered regions\tIUPred2A disordered percentage\tANCHOR2 binding regions\tANCHOR2 binding percentage\tVSL2B disordered regions\tVSL2B disordered percentage\n")
		allIDs = set().union(*[iupred, anchor, vsl2b])
		for full_id in natsorted(allIDs):
			if any([ full_id in iupred, full_id in anchor, full_id in vsl2b ]):
				ref_t, exitron_id = full_id.split(';')
				fout.write( "{}\t{}\t{}-{}\t{}\t{}\t{}\n".format(ref_t, exitron_id, EI_pos[exitron_id]["ei_from"], EI_pos[exitron_id]["ei_to"],
					iupred[full_id] if full_id in iupred else "none_found\tNA",
					anchor[full_id] if full_id in anchor else "none_found\tNA",
					vsl2b[full_id] if full_id in vsl2b else "none_found\tNA") )

### Get splice site scores ###
def read_PWM_from_file(PWM_file, region=None):

	# Parse only a specific region of the PWM
	if region:
		region = [ str(x+1) if x >= 0 else str(x) for x in [ i for i in range(*map(int, region.split(','))) ] ]

	PWM, parse, bases = [], True, []
	for line in open(PWM_file):
		if not line.startswith('#'):

			if line in ['\n', '\r\n']: parse = False
			elif parse:
				cols = line.rstrip().split('\t')
				try:
					total = sum([ float(x) for x in cols[1:] ])

					# Only parse position specified in region
					if region:
						if cols[0] in region:
							PWM.append( { bases[i-1]: float(cols[i])/total for i in range(1, len(cols)) } )

					# No region specified, so parse everything
					else:
						PWM.append( { bases[i-1]: float(cols[i])/total for i in range(1, len(cols)) } )

				except ValueError:
					bases = cols[1:]

	min_seq = ''.join([ min(p, key=p.get) for p in PWM ])
	max_seq = ''.join([ max(p, key=p.get) for p in PWM ])

	return PWM, min_seq, max_seq

def score_splice_site(seq, PWM):

	return sum([ log((PWM[i][base]/0.25)+0.0001, 2) if base in PWM[i].keys() else log(0.0001, 2) for i, base in enumerate(seq)  ])

def update_min_max(seq, PWM):

	score = score_splice_site(seq, PWM['pwm'])
	PWM['min'] = score if score < PWM['min'] else PWM['min']
	PWM['max'] = score if score > PWM['max'] else PWM['max']

	return PWM

def rescale_score(PWM_score, ss_min, ss_max):
    ''' Scale the PWM LOD score to 0-100

        ((b-a)*(PWM_score-min)/(max-min))+a

        If the PWM score is positive -> a = 50, b = 100, min = 0
        If the PWM score is negative -> a = 0, b = 50, max = 0
    '''

    if PWM_score > 0: norm_score = ((50*PWM_score)/ss_max)+50
    else: norm_score = -(50*(PWM_score-ss_min)/ss_min)

    return norm_score

def score_sites(work_dir, bkgd_fa, test_fa, PWM5_file, PWM3_file, five_prime_region, three_prime_region, file_handle):

	# POTENTIAL TO-DO: MAYBE IT WOULD BE BETTER TO SCALE AGAINST THE MINIMUM
	# AND MAXIMAL SCORES OF THE PWM, INSTEAD OF AGAINST A BACKGROUND SET

	# Read PWM, take the splice site regions
	# into account
	#PWM5 = { 'pwm': read_PWM_from_file(PWM5_file, five_prime_region), 'min': 0, 'max': 0 }
	#PWM3 = { 'pwm': read_PWM_from_file(PWM3_file, three_prime_region), 'min': 0, 'max': 0 }

	p5, min5, max5 = read_PWM_from_file(PWM5_file, five_prime_region)
	PWM5 = { 'pwm': p5, 'min': score_splice_site(min5, p5), 'max': score_splice_site(max5, p5) }
	p3, min3, max3 = read_PWM_from_file(PWM3_file, three_prime_region)
	PWM3 = { 'pwm': p3, 'min': score_splice_site(min3, p3), 'max': score_splice_site(max3, p3) }

	five_size = sum([ abs(x) for x in map(int, five_prime_region.split(',')) ])
	three_size = sum([ abs(x) for x in map(int, three_prime_region.split(',')) ])

	# Score the background splice sites
	#for record in yield_fasta(bkgd_fa):
	#	PWM5 = update_min_max(record.seq[:five_size].upper(), PWM5)
	#	PWM3 = update_min_max(record.seq[-three_size:].upper(), PWM3)

	# Score the test splice sites
	test, order = {}, []
	for record in yield_fasta(test_fa):

		five_seq = record.seq[:five_size].upper()
		three_seq = record.seq[-three_size:].upper()

		test[record.id] = {
			'5seq': five_seq,
			'3seq': three_seq,
			'5score': score_splice_site(five_seq, PWM5['pwm']),
			'3score': score_splice_site(three_seq, PWM3['pwm'])
		}
		#PWM5 = update_min_max(five_seq, PWM5)
		#PWM3 = update_min_max(three_seq, PWM3)

		order.append(record.id)

	# Output
	with open(work_dir+file_handle+".txt", 'w') as fout:
		for x in order:
			fout.write( "{}\t{}\t{}\t{:.3f}\t{:.3f}\n".format(x,
				test[x]['5seq'],
				test[x]['3seq'],
				rescale_score(test[x]['5score'], PWM5['min'], PWM5['max']),
				rescale_score(test[x]['3score'], PWM3['min'], PWM3['max']))
			)

### Get other stats ###
def get_misc(work_dir, exitron_fasta, get_length=False, get_GC=False, get_dinu=False):

	if any([ get_length, get_GC, get_dinu ]):

		with open(work_dir+"misc.txt", 'w') as fout:

			for record in yield_fasta(exitron_fasta):
				output = [record.id]

				# Get sequence length
				if get_length: output.append(str(len(record.seq)))

				# Get GC content
				if get_GC: output.append('{:.3f}'.format(sum([ record.seq.count(char) for char in ['G', 'C'] ]) / len(record.seq) * 100))

				# Get dinucleotides, assuming that the fasta
				# sequence is an intronic sequence
				if get_dinu: output.append('{}-{}'.format(record.seq[:2], record.seq[-2:]))

				fout.write( '\t'.join(output)+'\n' )

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-w', '--work-dir', default="./features/", help="Output directory.")

	subparsers = parser.add_subparsers(dest="command", help="sub-command help")

	parser_a = subparsers.add_parser("get-sequences", help="Get the exitron DNA and protein sequences.")
	parser_a.add_argument('-b', '--bed', required=True, help="Exitron bed file.")
	parser_a.add_argument('-g', '--genome-fasta', required=True, help="Genome fasta.")
	parser_a.add_argument('-o', '--output-file', required=True, help="Output fasta file.")
	parser_a.add_argument('--upstream-nt', type=int, default=0, help="Number of nucleotides to read into the upstream (5'SS) exon.")
	parser_a.add_argument('--downstream-nt', type=int, default=0, help="Number of nucleotides to read into the downstream (3'SS) exon.")

	parser_b = subparsers.add_parser("prepare", help="Get the EIx3 sequences, reference protein sequences and EIx3 relative positions.")
	parser_b.add_argument('--exitrons-info', required=True, help="exitrons.info file (from identify-exitrons).")
	parser_b.add_argument('-g', '--genome-fasta', required=True, help="Genome fasta.")
	parser_b.add_argument('-c', '--ccds', required=True, help="CCDS GTF file.")
	parser_b.add_argument('-r', '--reference-proteins', required=True, help="Reference transcript translations.")

	parser_c = subparsers.add_parser("get-domains", help="Get Pfam protein domains in the exitrons.")
	parser_c.add_argument('-p', '--pfam', required=True, help="Pfam hmm database.")
	parser_c.add_argument('-f', '--fasta', required=True, help="Exitron-containing reference transcript protein sequences.")
	parser_c.add_argument('--exitron-pos', required=True, help="Exitron position on reference transcript protein sequences.")
	parser_c.add_argument('-e', '--e-value', default="0.001", type=float, help="Max e-value cut-off.")

	parser_d = subparsers.add_parser("get-disordered", help="Identify disordered regions with IUPred2A and VSL2B.")
	parser_d.add_argument('-f', '--fasta', required=True, help="Exitron-containing reference transcript protein sequences.")
	parser_d.add_argument('--exitron-pos', required=True, help="Exitron position on reference transcript protein sequences.")

	parser_e = subparsers.add_parser("score-sites", help="Score splice sites based on supplied PWMs.")
	parser_e.add_argument('-b', '--bkgd', required=True, help="Background fasta sequences for score scaling.")
	parser_e.add_argument('-t', '--test', required=True, help="Fasta sequences to score.")
	parser_e.add_argument('--PWM5', required=True, help="5' splice site PWM file.")
	parser_e.add_argument('--PWM3', required=True, help="3' splice site PWM file.")
	parser_e.add_argument('-f', '--file-handle', default="all", help="File handle for the output file.")
	parser_e.add_argument('--five-prime-region', default="-3,+6", help="Specify the number of exonic bases (-x) and intronic bases (+y) for the 5' splice site.")
	parser_e.add_argument('--three-prime-region', default="-6,+3", help="Specify the number of intronic base (-x) and exonic bases (+y) for the 3' splice site.")
	group = parser_e.add_mutually_exclusive_group()
	group.add_argument('--ath', action='store_true', help="Default preset for Arabidopsis thaliana: 5'ss => -3,+10, 3'ss => -14,+3.")
	group.add_argument('--hsa', action='store_true', help="Default preset for Homo sapiens: 5'ss => -3,+6, 3'ss => -6,+3.")

	parser_f = subparsers.add_parser("misc", help="Get misc stats like GC content and dinucleotides.")
	parser_f.add_argument('-f', '--fasta', required=True, help="Exitron fasta sequences.")
	parser_f.add_argument('--length', action='store_true', help="Output the sequence length.")
	parser_f.add_argument('--gc', action='store_true', help="Output the sequence GC content.")
	parser_f.add_argument('--dinu', action='store_true', help="Assumes the fasta sequence is an intron sequence. Output the dinucleotides.")

	args = parser.parse_args()

	work_dir = args.work_dir if args.work_dir.endswith('/') else args.work_dir+'/'
	if not os.path.exists(work_dir): os.makedirs(work_dir)

	if args.command == "get-sequences":
		exitron2fasta(work_dir, args.bed, args.genome_fasta, args.output_file, args.upstream_nt, args.downstream_nt)

	elif args.command == "prepare":
		get_relative_exitron_position(work_dir, args.exitrons_info, args.genome_fasta, args.ccds, args.reference_proteins)

	elif args.command == "get-domains":
		get_protein_domains(work_dir, args.pfam, args.fasta, args.exitron_pos, args.e_value)

	elif args.command == "get-disordered":
		get_disordered_regions(work_dir, args.fasta, args.exitron_pos)

	elif args.command == "score-sites":

		# Set presets
		if args.ath:
			args.five_prime_region = "-3,+10"
			args.three_prime_region = "-14,+3"
		elif args.hsa:
			args.five_prime_region = "-3,+6"
			args.three_prime_region = "-6,+3"

		score_sites(work_dir, args.bkgd, args.test, args.PWM5, args.PWM3, args.five_prime_region, args.three_prime_region, args.file_handle)

	elif args.command == "misc":
		get_misc(work_dir, args.fasta, args.length, args.gc, args.dinu)

# Update log:
# * 11-12-2019: fixed mistake with relative positioning of negative strand exitrons
