# PAH-MC

Monte Carlo Simulation on deuterated PAHs.

Associated research which include usage:

1. [Carels, E. J. Isotopic Effects Revealed Upon Photolysis of PAHs. (Universiteit van Amsterdam, 2023).](https://scripties.uba.uva.nl/search?id=record_53847)

## Installation

Install the prerequisites using Mamba or another environment manager:

```python
mamba create -n PAH-MC numpy
```

## Usage

How to run the program:

```bash
usage: PAHMC [-h] [-o OUTPUT] [-l LOG] inputfile cores

Process input file and rate definition file.

positional arguments:
  inputfile            Input YAML file
  cores                Number of parallel processes to run

options:
  -h, --help           show this help message and exit
  -o, --output OUTPUT  Output file
  -l, --log LOG        Log file
```

If no Output or Log files are specified the program uses the name of `inputfile` with the `.out` and `.log` extensions.

If you are going to run this program on a supercomputer/computer cluster that uses SLURM for the job management create a script and run it:

```bash
#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH --tasks-per-node 36

module load 2022
module load Python/3.10.4-GCCcore-11.3.0
module load Anaconda3/2022.05

python src/main.py $1 36
```

```bash
sbatch mc_slurm.sh <inputfile> 
```

This script has been successfully used on the Snellius computer cluster/supercomputer by Surf.

## How it works

![Code2flow diagram](out.png)

## Credits and acknowledgements

This program was written by Emma Carels for a master project. It was updated by Tim de Groot for a subsequent bachelor project.
