# imports
import os

# VARS FOR SWEEP (that all aligner x pc combo uses)
genome = ['600bp.fa', '1pct.fa'] 
covs = [1, 2, 3, 4, 5]
# pvals = [0.05, 0.01, 0.001]
peak_nums = [1, 2, 3]
peak_tallness = [10, 100, 1000]
peak_broadness = [10, 100, 1000]
#aligner = ['bowtie2', 'bwa-mem']
#peak_callers = ['macs2', 'epic2']

# VARS FOR SWEEP (specific to aligner x pc combo)
# stores both genome name, genome.fa, and genome index file in a tuple; for bowtie+bwa mem
bowtie_genome_and_idx = [ 
	('600bp', '/genomes/600bp.fa''/path/to/600bp_genome_index'), 
	('1pct', '/genome/1pct.fa','/path/to/1pct_c_elegans_index')
	]
bowtie_genome_and_idx = [ 
	('600bp', '/genomes/600bp.fa''/path/to/600bp_genome_index'), 
	('1pct', '/genome/1pct.fa','/path/to/1pct_c_elegans_index')
	]
nomodel = [True, False] # only for macs2


# CONSTANT VARS
RESULTS = '/results/'

# chipseq reads
READS = os.path.join(RESULTS, 'reads/')

# reads aligned
ALIGN = os.path.join(RESULTS, 'align/')

# peaks
PEAKS = os.path.join(RESULTS, 'peaks/')

"""
Rules
"""

# note: have one rule each for bowtie2 and bwa-mem alignment
	# no matter how complex the DAG is, snakemake should be able to handle it
	# as long as the output file paths are unique
	
rule all:
	# final name includes all vars for sweep
	# concern: we would have to include parameters specific to specifc aligner x pc combos
		# solution: could just create separate expand() for each specific aligner x pc combo

	# for bowtie2 x macs2
		# {nomodel} is unique to macs2
	expand('bowtie2_macs2_g_{gen}_cov_{cov}'
			'_pnum_{pnum}_ptall_{ptall}'
			'_pbroad_{pbroad}_nomodel_{nomodel}.peak', 
			cov=covs, pval=pvals, peak_num=peak_nums)

	# bowtie2 x epic2
		# no unique vars
	expand('bowtie2_epic2_g_{gen}_cov_{cov}'
			'_pnum_{pnum}_ptall_{ptall}'
			'_pbroad_{pbroad}.peak',
			cov=covs, pval=pvals, peak_num=peak_nums)

	# bwa x macs2
		# {nomodel} is unique to macs2
	expand('bwa_macs2_g_{gen}_cov_{cov}'
			'_pnum_{pnum}_ptall_{ptall}'
			'_pbroad_{pbroad}_nomodel_{nomodel}.peak', 
			cov=covs, pval=pvals, peak_num=peak_nums)

	# bwa x epic2
		# no unique vars
	expand('bwa_epic2_g_{gen}_cov_{cov}'
		'_pnum_{pnum}_ptall_{ptall}'
		'_pbroad_{pbroad}.peak',
		cov=covs, pval=pvals, peak_num=peak_nums)
	
rule chipseq:
	params:
		exp_con = ['exp', 'con']
	output: 
		os.path.join(EXP_READS, \ # might need to use expand() here?
		'g_{gen}_cov_{cov}'
		'_pnum_{pnum}_ptall_{ptall}'
		'_pbroad_{pbroad}_chipseq_reads.fa')

		os.path.join(CONTROL_READS, \
		'g_{gen}_cov_{cov}'
		'_pnum_{pnum}_ptall_{ptall}'
		'_pbroad_{pbroad}_chipseq_reads.fa')
	shell:
		python3 chip_seq.py --genome {gen} --cov {cov} # FINISH THE ARGUMENTS
		> chipseq.output[1] # NOT SURE IF THIS IS CORRECT SYNTAX BUT 
		# IT SHOULD BE SAME PATH AS FIRST LISTED OUTPUT PATH
		python3 chip_seq.py --genome {gen} --cov {cov} # FINISH THE ARGUMENTS
		> chipseq.output[2] # IT SHOULD BE SAME PATH AS SECOND LISTED OUTPUT PATH

rule align_bowtie2:
	input:
		chipseq.output[1]
		chipseq.output[2]
	params:
		exp_con = ['exp', 'con']
	output:
		# NOTE: most peakcallers will not use default SAM or BAM, so exclude
		
		# sorted bam
		os.path.join(ALIGN, \
			'{params.exp_con}_' # Q. is there a way to store these universal parameter sweep variables
			'g_{gen}_'				# in a var somewhere so I don't have to write it out everytime?
			'cov_{cov}_'
			'pnum_{pnum}_'
			'ptall_{ptall}_'
			'pbroad_{pbroad}'
			'.sorted.bam')

		# sorted bai
		os.path.join(ALIGN, \
			'{params.exp_con}_'
			'g_{gen}_'
			'cov_{cov}_'
			'pnum_{pnum}_'
			'ptall_{ptall}_'
			'pbroad_{pbroad}'
			'.sorted.bai')

	shell:
		# something like this
		# first need to create genome index
		samtools {gen} > {os.path.join(ALIGN, '{gen}_index')} # something like this

		# For EXP
		# then align using created genome index
		bowtie2 align -genome_index {os.path.join(ALIGN, '{gen}_index')} -reads_fasta chipseq.output[1] | \
		bowtie2 ... | \ # stdout (the SAM) pipe into BAM 
		bowtie2 ... > # convert to sorted.bam save the file in ALIGNMENT

	

rule align_bwa:








	