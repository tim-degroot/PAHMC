# RRKM

The `rrkm.py` script used Rice-Ramsperger-Kassel-Marcus (RRKM) theory to evaluate the reaction rates for scrambling and dissociation rates of deuterium and hydrogen in Polycyclic Aromatic Hydrocarbons (PAHs) and to prepare input files for use with the accompanying PAHMC model for Monte Carlo simulations.

## Installation

Make sure the prerequisites are installed: `numpy, matplotlib`

Test to see if `rrkm.py` is functioning. The current version works with MULTWELL 2023.2.1. Depending on changes to the densum output file you might have to update the `start_line` calculation.

The Python script can be made accessible by adding it to the bin and changing permission.

```bash
cp rrkm.py /usr/local/bin/rrkm.py
chmod 755 /usr/local/bin/rrkm.py
```

## Usage

```docs
Usage:
    rrkm.py [options] <reactant> <TS> <product>
    rrkm -h | --help
Options:
    -q  Quiet output.
    -p  Plot rate.
    -l  Loose transition state -> needed for diss of H/D (not an actual TS but gradual thing).
    -f  Scaling factor [default: 0.97].
    -t  Temperature [default: 1000].
    -g  Energy grain in cm-1 -> 100=0.012eV (CORO); 400=0.048eV (OVA) [default: 100].
    -m  Maximum internal energy to evaluate (in cm-1) [default: 242000].
    -o  Output file [default rrkm.txt].
    -r  Reverse output file.
```

Example:

```bash
rrkm.py -o /RRKM/Phenanthrene/C14H10-D9to10.txt /DFT/Phenanthrene/C14H10-D9.log /DFT/Phenanthrene/C14H10-TS-D9toD10.log /DFT/Phenanthrene/C14H10-D10.log
```
