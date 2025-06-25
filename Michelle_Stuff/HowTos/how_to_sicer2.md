</h2>How to use Sicer2</h2><br>

You do need your python in your yml file to run 3.10 as Sicer2 does not support python 3.13

1. Convert BAM file to BED file<br>
    - make sure both samtools and bedtools are in yml dependencies<br>
    - have both your experiment and control bam files in the pwd
    - run both bam files: `bedtools bamtobed -i <name.bam> > <name.bed>`

2. Install Sicer2
    - Numpy and scipy are needed to run Sicer 2: `conda install -y numpy scipy`<br>
    - Install sicer: `pip install SICER2`<br>

3. Adding your own genome sizes to Sicer2
    - locate your `GenomeData.py`
    - I found mine in miniconda3\envs\sicer2\lib\python3.10\site-packages\sicer\lib
    - add your genome in: `<genome_name>_chroms = ['chr1', 'chr2',...]`
    - set your chromosome sizes: `<genome_name>_chrom_lengths = {'chr1': <size>, 'chr2': <size>,...}`
    - at the bottom, at species_chroms and species_chrom_lengths, define your genome name: `"random1": random1_chroms` and `"random1": random1_chrom_lengths`
    - save your changes

4. Run Sicer2
    - example:
        - `sicer -t random_genome_1_cov_100_output.bed -c unbiased_control_cov_100_random_genome_1_output.bed --species random1 --o random_1_sicer2`
    - -t = treatment/experimental bed file
    - -c = control bed file
    - -s/--species = default species or your own
    - --o = creates an output directory

Other possible parameters:<br>
    - -w/--window_size = window size. Default: 200bp<br>
    - -f/--fragment_size = amount of shift from the beginning of the read to the center of the DNA fragment represented by the read. Default: 150bp<br>
    - -egf/--effective_genome_fraction = effective_genome_fraction. Default: 0.74<br>
    - -fdr/--false_discovery_rate. Default: 0.01<br>
    - -g/--gap_size = minimum length of a "gap" such that neighboring window is an "island." Default: 600bp<br>
    - -e/--e_value: needed when no control library is provided. Default: 1000<br>
