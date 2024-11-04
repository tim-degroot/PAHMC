# PAH-MC
Monte Carlo Simulation on deuterated PAHs

Required Python Libraries:
numpy
os
re
sys
random
multiprocessing
copy
getopt


How to run the program:

python main.py [-o <output>] [-l <log>] <inputfile> <ratedef-file> <rate-directory> <#cores>

The two statements in brackets are optional for if you want to specify a specific filename for either your output file or your log file, if not it standardly takes the name of your input file with .out or .log as an extension.

The number of cores specifies the number of parallel processes run.

If you are going to run this program on a supercomputer/computer cluster that uses SLURM for the job management, run the following:

sbatch mc_slurm.sh <inputfile> <ratedef-file> <rate-directory>

Number of processes and the optional arguments can be edited in the mc_slurm.sh file. Take mind to keep the number of cores specified for the program the same as the number of cores used in total (number of nodes times the number of tasks per node, both specified in the SBATCH statements).

This script has been successfully used on the Snellius computer cluster/supercomputer by Surf.


This program was written by Emma Carels for a master project. If any errors or questions arise, feel free to contact me at emma@carels.info.

