# How `bt2_2_homer.py` works

This also works soft links to directories.

There are 3 required arguments and 1 optional argument.

In this order,
- path to genome (in .fa or .fa.gz)
- path to directory containing all experimental reads
- path to directory containing all control reads

The optional argument is `--output`. This indicates where all the outputs for this pipeline will go. The default is set to the current directory.

After running the script with its arguments, 2 additional directories will be created.
- `index`: index directory is required for bowtie2 to run
- `results`: all the results of the peakcaller will appear here as `.txt` files