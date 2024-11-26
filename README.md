# PAH-MC

Perform Monte Carlo simulation of scrambling and photodissociation reactions on PAHs.

Associated research:

1. [Carels, E. J. Isotopic Effects Revealed Upon Photolysis of PAHs. (Universiteit van Amsterdam, 2023).](https://scripties.uba.uva.nl/search?id=record_53847)

## Installation

Clone the repository.

Install the prerequisites using Mamba or another environment manager:

```python
mamba create -n PAHMC numpy, matplotlib, pyyaml
```

If installing on a computer cluster you can skip the environment and instead load the modules using [SLURM](https://slurm.schedmd.com) or a similar workload manager (see [Running the program](#running-the-program)).

## Usage

### Preparing the `.yaml` input file

The input file uses [YAML](https://yaml.org). While the documentation explains the file structure in a certain order this order is not necessary for the program to work.

We will prepare a simulation of Anthracene in this example to explain the input file structure.

Provide a name for the simulation (optional, currently not used for anything):

```yaml
Name: Anthracene
```

The edge consists of groups of substituents that are linked by links. The edge loops around to itself.

Define the initial edge substituents (0 for none, H, D, HH or HD). Only one aliphatic site (HD/HH) supported. Label them according to nomeclature using edge numbers.

```yaml
Initial edge: (H,H,H,H) (HD) (H,H,H,H) (H)
Initial edge numbers: (5,6,7,8) (9) (1,2,3,4) (10)
```

Define the initial tertiary links (0 for None). The links are added inbetween the edge groups starting after the first edge. Starting with aliphatic tertiary is currently not supported. Label them according to nomenclature using link numbers.

```yaml
Initial links: (0) (0) (0) (0)
Initial link numbers: (8a) (9a) (4a) (10a)
```

Define any cross-links (optional). Cross-links can be any reaction that is not between neighbours on the edge structure. These are defined in pairs.

```yaml
Initial cross-links: (8a,10a) (9a,4a)
```

Define the simulation energy in eV:

```yaml
Energy: 4.66
```

Define the number of iterations:

```yaml
Iterations: 2000
```

Define the maximum simulation time (s):

```yaml
Maximum time: 0.02
```

Define error handling (w for warming, q for quit).

```yaml
Error handling: w
```

Define the rates by pointing to the path where the RRKM rate files (prepared with [rrkm.py](https://github.com/tim-degroot/PAHMC/tree/main/RRKM)) are stored and matching all reactions (`HAtoB` and `DAtoB` for every neighbour and cross-link pair). Symmetrical reactions can point to the same rate file.

```yaml
Rates path: /RRKM/
Rates:
  H9to9a,H9to8a,H10to4a,H10to10a: H9to9a.txt
  D9to9a,D9to8a,D10to4a,D10to10a: D9to9a.txt
  H9ato9,H8ato9,H4ato10,H10ato10: H9ato9.txt
  D9ato9,D8ato9,D4ato10,D10ato10: D9ato9.txt
  ...
  H9diss,H10diss: diss_rate_H9.txt
  D9diss,D10diss: diss_rate_D9.txt
  ...
```

### Running the program

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

If no Output or Log files are specified the program uses the name of `inputfile` with different extensions (see [Understanding the output](#understanding-the-output)).

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

This script has been successfully used on the Snellius computer cluster by Surf using the above configuration..

### Understanding the output

The filename of the output is taken from the yaml (`ABC.yaml` gives `filename = ABC`) or can be set manually using the `--OUTPUT` parameter (see [Running the program](#running-the-program)).

The program generates several files. They will be discussed in alphabetical order (the same way your system  usually sorts them).

`{filename}_data.log` contains the following:

```md
MC#       	Diss atom 	Diss pos  	Diss time 	# hops    	# D hops  	# cross hops	HH time   	HD time   	DD time   
MC0       	None      	None      	2.0004e-06	4556      	348       	21          	1.5672e-06	4.331e-07 	0         
```

Where `MC#` is the iteration number; `Diss atom` is the type of atom that dissociated, either H(ydrogen) or D(euterium); `Diss pos` is the edge number where the dissociation occured; `Diss time` is the time the simulation ran before dissociation occured, if no dissociation happened but the simulation was ended because of the Maximum time this value is the time at that point; `# hops` is the total number of hops that happened in the simulation; `# D hops` is the number of hops of a deuterium atom; `# cross hops` is the number of hops across a cross-link; `HH/HD/DD time` is the time spent with a HH/HD/DD aliphatic site respectively.

`{filename}_key_hops.out` contains the following:

```md
MC#	D1to2	D1to9a	D2to1	D2to3	D9ato1	D9ato4a	D9ato9	D9to9a	H1to2	H1to9a	H2to1	H2to3	H9ato1	H9ato4a	H9ato9	H9to9a
MC0	147  	1     	147  	4    	2     	0      	23    	24    	1757 	40    	1757 	142  	39    	21     	226   	226   
```

Where the amount of occurences of every unique hop is stored. Symmetrical ones are summed under the label of the first in the list.

`{filename}_diss_time.png` shows a histogram of the dissociation times per iteration.

`{filename}_hops.png` shows a histogram of the number of hops per iteration.

`{filename}.log` contains the programs logger.

`{filename}.out`

If `--debug` is enabled `{filename}_{iteration}_mol_structures.log` contains the complete structure history of the molecular structure of the first iteration (iteration 0) is stored.

## How it works

![Code2flow diagram](out.png)

## Credits and acknowledgements

This program was created by Emma Carels for a master project. It was updated and is currently maintained by [Tim de Groot](tim.degroot@student.uva.nl) for a subsequent bachelor project.

Supervision for both projects was done by [dr. Alessandra Candian](a.candian2@uva.nl) and dr. ir. Annemieke Petrignani

All research was performed at the [Anton Pannekoek Institute for Astronomy](https://api.uva.nl) & [Van 't Hoff Institute for Molecular Sciences](https://hims.uva.nl) at the University of Amsterdam.
