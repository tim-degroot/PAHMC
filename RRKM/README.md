# RRKM

The `rrkm.py` script uses the `densum` program included in the `MULTIWELL`[^1][^2][^3] suite to calculate reaction rates in accordance with Rice-Ramsperger-Kassel-Marcus (RRKM) Theory and prepare reaction rate input files for use with the PAHMC model for Monte Carlo simulations.

[^1]: R. Barker, T. L. Nguyen, J. F. Stanton, C. Aieta, M. Ceotto, F. Gabas, T. J. D. Kumar, C. G. L. Li, L. L. Lohr, A. Maranzana, N. F. Ortiz, J. M. Preses, J. M. Simmie, J. A. Sonk, and P. J. Stimac; MultiWell-2023-2.1 Software Suite; J. R. Barker, University of Michigan, Ann Arbor, Michigan, USA, 2023; [MultiWell Site](https://multiwell.engin.umich.edu/).
[^2]: John R. Barker, Int. J. Chem. Kinetics, 33, 232-45 (2001).
[^3]: John R. Barker, Int. J. Chem. Kinetics, 41, 748-763 (2009).

## Installation

Make sure the prerequisites are installed: `numpy, matplotlib`

Test to see if `rrkm.py` is functioning. The current version works with MULTWELL 2023.2.1. Depending on changes to the densum output file you might have to update the `start_line` calculation.

The Python script can be made accessible by making it executable and adding it to the path.

```bash
chmod +x rrkm.py
export PATH="/path/to/PAHMC/RRKM:$PATH"
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
rrkm.py -o /RRKM/D10to10a.txt -r /RRKM/D10ato10.txt /DFT/C14H10-D10.log /DFT/C14H10-TS-D10toD10a.log /DFT/C14H10-D10a.log
```
