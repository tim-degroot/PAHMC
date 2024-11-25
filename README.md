# PAH-MC

Perform Monte Carlo simulation of scrambling and photodissociation reactions on PAHs.

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
usage: PAHMC [-h] [-c CORES] [-o OUTPUT] [-l LOG] [-d] inputfile

Perform Monte Carlo simulation of scrambling and photodissociation reactions on PAHs.

positional arguments:
  inputfile            Input YAML file

options:
  -h, --help           show this help message and exit
  -c, --cores CORES    Number of parallel processes to run
  -o, --output OUTPUT  Output file
  -l, --log LOG        Log file
  -d, --debug          Enable debugging
```

If no Output or Log files are specified the program uses the name of `inputfile` with the `.out` and `.log` extensions.

If you are going to run this program on a supercomputer/computer cluster that uses SLURM for the job management create a script and run it:

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 128
#SBATCH --partition=rome

module load 2024
module load Python/3.12.3-GCCcore-13.3.0
module load Anaconda3/2024.06-1

python ../PAHMC/src/main.py $1
```

```bash
sbatch mc_slurm.sh <inputfile> 
```

This script has been successfully used on the Snellius computer cluster/supercomputer by Surf.

## How it works

![Code2flow diagram](out.png)

## Credits and acknowledgements

This program was written by Emma Carels for a master project. It was updated by Tim de Groot for a subsequent bachelor project.
