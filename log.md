# Log for this project
## Question
How well do different peak-finding programs actually work when given proper Chip-seq data?

How do different peak-finding programs handle small sample sizes?
> At what sample size do different peak-finding programs become accurate at distingushing randomness to actual peaks?

## Progress
### Feb 10
Kicked off project after meeting with Dr. Korf. Supposed to be simple program.

### Feb 14
Created repo `https://github.com/RussRussBus/Creating_Controls.git` and talked to Dr. Korf to see if I am on track. Created random version of the program but have to create the enumerate version of it. Need to also make outputs of randogen more clear.

### Feb 24
Create for loop that accepts inputs: number of peaks, size of the genome, and the height of the peak. It then creates a list that has every possible combination of positions of the peaks in the genome.

### Feb 28
After further testing and inspection the for loop only finds some of the combinations. Plan to fix now to work as intended.

Got the program to emunmate APC for genome size 5 and 3 peaks. Need to make it work with other combination of genome sizes and number of peaks.

Program works with all number of peaks except 1 and 2. Will probably have to write separate code for these cases.

Wrote a function to handle situations when number of peaks is 1 or 2. Also made the program into a function. Allows either function to be called on depending on number of peaks.

Need to beautify and comment out code. Also need to add height and multipliers.

### Mar 2
I retract previous statement, code does not work with number of peaks greater than 3.

Restarting the portion of enumerate every possible peak position in the genome. Did not like how messy and confusing code was getting.

No real progress today just going in circles trying to get all unique combinations.

Choose function good for knowing the total number of combinations.

### Mar 25
Forgot to mark down notes after meeting with Dr.Korf. He talked about steps for the project.

Step 1: Creating random genome

Step 2: Create multiplier deterministically

Step 3: Create peaks

As of today, a first attempt was made at writing a script for the step 1 and step 2. They are finished in the sense that they produce the intended output. However, they are no where close to being using in any other experiments.

The genome is just randomly generated nucleotides. There is no biological situations taken into consideration when in the genome generation.

The multiplier is just a value that has been assigned to regions of the genome, and this conveyed by having an ASCII value below each nucleotide. The ASCII character "A" is equivalent to decimal value "0", and the decimal values are in log base 2.

The way peaks are created needs to be discussed because they can be generate randomly, based on the nucleotide sequence, deterministically, or some other way.