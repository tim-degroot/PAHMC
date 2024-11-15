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
usage: RRKM [options] reactant TS product

Calculate RRKM rate based on Gaussian16 frequency calculations using densum

positional arguments:
  reactant    Reactant .log file
  TS          Transition State .log file
  product     Product .log file

options:
  -h, --help  show this help message and exit
  -q          Quiet mode (default: False)
  -p          Plot mode (default: False)
  -l          Loose transition state, i.e. dissociation (default: False)
  -f F        Scaling factor (default: 0.97)
  -t T        Temperature (default: 1000)
  -g G        Energy grain in cm-1 (default: 100)
  -m M        Maximum internal energy to evaluate (default: 242000)
  -o O        Output file (default: rrkm.txt)
  -r R        Reverse output file (default: None)
```

Example:

```bash
rrkm.py -o /RRKM/D10to10a.txt -r /RRKM/D10ato10.txt /DFT/C14H10-D10.log /DFT/C14H10-TS-D10toD10a.log /DFT/C14H10-D10a.log
```
