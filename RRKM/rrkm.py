#! /usr/bin/env python3

import os, sys, tempfile, re, getopt
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import argparse

cwd = os.getcwd()

### Scaling factor
sf = 0.97
### Conv. eV to cm-1
conv = 8065.73
### Plack constant (cm-1*s)
h = 4.135667516e-15 * conv

### Ensamble temperature (K)
T = 1000.0
### Plack constant (J*cm)
h_ent = 1.98644568e-23
### Boltzmann constant (J/K)
kb = 1.38064852e-23
### Gas constant (cal/K/mol)
R = 8.3144598


parser = argparse.ArgumentParser(
    prog="RRKM",
    usage="RRKM [options] reactant TS product",
    description="Calculate RRKM rate based on Gaussian16 frequency calculations using densum",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("reactant", type=str, help="Reactant .log file")
parser.add_argument("TS", type=str, help="Transition State .log file")
parser.add_argument("product", type=str, help="Product .log file")
parser.add_argument("-q", action="store_true", help="Quiet mode")
parser.add_argument("-p", action="store_true", help="Plot mode")
parser.add_argument(
    "-l", action="store_true", help="Loose transition state, i.e. dissociation"
)
parser.add_argument("-f", type=float, default=0.97, help="Scaling factor")
parser.add_argument("-t", type=int, default=1000, help="Temperature")
parser.add_argument("-g", type=int, default=100, help="Energy grain in cm-1")
parser.add_argument(
    "-m", type=int, default=242000, help="Maximum internal energy to evaluate"
)
parser.add_argument("-o", type=str, default="rrkm.txt", help="Output file")
parser.add_argument("-r", type=str, help="Reverse output file")
args = parser.parse_args()

flq = args.q
flp = args.p
flloose = args.l
sf = args.f
temp = args.t
Egrain = args.g
Emax2 = args.m
outfile = args.o
outfile_rev = args.r if args.r else ""


E0 = []
energ = []
dos = []
sos = []
Stemp = []

functionals = []

tmpdir = tempfile.TemporaryDirectory()
structures = [args.reactant, args.TS, args.product]

for struct in structures:
    if struct[-4:] != ".log" or not os.path.isfile(struct):
        print(f"File {struct} is not valid!")
        sys.exit(2)

    fh_freqfile = tempfile.NamedTemporaryFile(mode="w+", dir=tmpdir.name)
    freqfile = fh_freqfile.name
    fh_densfile = tempfile.NamedTemporaryFile(mode="w+", dir=tmpdir.name, delete=False)
    densfile = fh_densfile.name

    ### Read frequencies from Gaussian file
    if not flq:
        print(f"Getting frequencies from {struct}...")
    subprocess.call(["grep", "Frequencies", struct], stdout=fh_freqfile)

    ### Get ground state energies
    grelec = 0.0
    ZPVE = 0.0
    funbset = ""
    for line in open(struct):
        if "#" in line and funbset == "":
            funbset = re.findall(r".+?\s", line)[1]
            funbset = funbset[:-1]
        if "Done" in line:
            grelec = float(re.findall(r"[-+]?\d+[\.]?\d*[eE]?[-+]?\d*\b", line)[0])
        if "Kcal" in line:
            ZPVE = float(re.findall(r"[-+]?\d+[\.]?\d*[eE]?[-+]?\d*\b", line)[0])
    grelec *= 27.2107  ## Conv. to eV
    ZPVE *= 0.0433634  ## Conv. to eV
    tot_ener = grelec + ZPVE

    functionals.append(funbset)
    E0.append(tot_ener)

    ### Create freq. array
    freqs = np.loadtxt(freqfile, usecols=(2, 3, 4))
    fh_freqfile.close()
    freqs = freqs.flatten()
    freqs *= sf
    ai = freqs.size

    ### Remove imaginary frequencies
    freqs = freqs[freqs > 0]
    af = freqs.size
    if not flq:
        print(f"Freqs. read: {ai}")
        print(f"Imaginary freqs.: {ai - af}")

    theta = h_ent * freqs / kb

    Q = np.prod(1 / (1 - np.exp(-theta / T)))
    U = R * np.sum(theta / (np.exp(theta / T) - 1))

    S = R * np.log(Q) + U / T
    Stemp.append(S)

    ### Title for densum file
    title = "dummy"
    ### SET line1 for densum input file: TITLE
    line1 = f"{title} sf={sf:.2f} {funbset}"
    ### SET line2 for densum input file: filename for multiwell
    line2 = title
    ### SET line3 for densum input file
    IWR = 1  # 0=Exact Count; 1=Whitten-Rabinovich
    ### DoF, IWR, 'HAR'(harmonic)/'OBS'(fundamental), units ('AMUA')
    line3 = [af, IWR, "'HAR'", "'AMUA'"]
    ### SET line4 for densum input file
    Imax1 = 10000  # number of array elements in first segment of double array
    Isize = 500000  # total size of double array (number of elements)
    ### Egrain, Imax1, Isize, Emax2
    line4 = [Egrain, Imax1, Isize, Emax2]

    ### Create densum file
    print(line1, file=fh_densfile)
    print(line2, file=fh_densfile)
    print(*line3, sep=",", file=fh_densfile)
    print(*line4, sep=",", file=fh_densfile)
    for i, freq in enumerate(freqs):
        print(i + 1, "vib", freq, 0.0, 1, sep=",", file=fh_densfile)
    ### Final empty line is required for densum
    print("", file=fh_densfile)
    fh_densfile.close()
    if not flq:
        print(" Densum file written.")

    ### Run densum
    os.chdir(tmpdir.name)
    fh_densfile.close()

    subprocess.call(["cp", densfile, os.path.join(tmpdir.name, "densum.dat")])
    if not flq:
        print(" Running densum...")
    subprocess.call(["densum", densfile], stdout=tempfile.TemporaryFile())
    if not flq:
        print(" Done.")
        print("")

    fh = open("densum.out", "r")

    start_line = 65 + af

    en = []
    dens = []
    suma = []

    for i, line in enumerate(fh):
        if i >= start_line - 1:
            data = line.split()
            if not flq:
                print(data)
            en.append(float(data[1]))

            check = re.compile(r"E\+")
            if not check.search(data[2]):
                data[2] = data[2].replace("+", "E+")
            dens.append(float(data[2]))
            if not check.search(data[3]):
                data[3] = data[3].replace("+", "E+")
            suma.append(float(data[3]))

    sos.append(np.array(suma))
    dos.append(np.array(dens))
    energ = np.array(en)

    fh.close()
    os.chdir(cwd)

if not flq:
    if functionals[0] == functionals[1] and functionals[1] == functionals[2]:
        print("Functionals/basis set are the same for all structures.")
    else:
        print("WARNING: Functional/basis set are not the same")

barrier = E0[2] - E0[0] - 13.6 if flloose else E0[1] - E0[0]
# barrier = 3.18
barrier *= conv
delta = E0[2] - E0[0] - 13.6 if flloose else E0[2] - E0[0]
# delta = 2.36

DS = Stemp[1] - Stemp[0]

rate = np.zeros_like(energ)

rate[energ > barrier] = (
    sos[1][energ < np.max(energ) - barrier] / dos[0][energ > barrier]
)
rate /= h

barrier /= conv

fh_rrkmfile = open(outfile, "w")

print(
    "RRKM rate, sf={:.3f}, S*={:.2f}, E0={:.2f} and Delta={:.2f}".format(
        sf, DS, barrier, delta
    ),
    file=fh_rrkmfile,
)
print("{:>10s} {:>10s}".format("En (eV)", "k (s-1)"), file=fh_rrkmfile)

for E, k in zip(energ / conv, rate):
    print("{:10.5g} {:10.5g}".format(E, k), file=fh_rrkmfile)

if not flq:
    print("")
    print("RRKM file written.")
    print("sf={:.3f}; Delta={:.2f} eV".format(sf, delta))
    print("E0={:.2f} eV; S*={:.2f}".format(barrier, DS))
    print("Number of points: " + str(len(energ)))

if flp:
    plt.plot(energ / conv, rate, label="Forward")

fh_rrkmfile.close()

if outfile_rev:
    barrier -= delta  # E0[0]-E0[2]-13.6 if flloose else E0[1]-E0[2]
    barrier *= conv
    delta = -delta  # E0[0]-E0[2]-13.6 if flloose else E0[0]-E0[2]

    DS = Stemp[1] - Stemp[2]

    rate = np.zeros_like(energ)

    rate[energ > barrier] = (
        sos[1][energ < np.max(energ) - barrier] / dos[2][energ > barrier]
    )
    rate /= h

    barrier /= conv

    fh_rrkmfile_rev = open(outfile_rev, "w")

    print(
        "RRKM rate, sf={:.3f}, S*={:.2f}, E0={:.2f} and Delta={:.2f}".format(
            sf, DS, barrier, delta
        ),
        file=fh_rrkmfile_rev,
    )
    print("{:>10s} {:>10s}".format("En (eV)", "k (s-1)"), file=fh_rrkmfile_rev)

    for E, k in zip(energ / conv, rate):
        print("{:10.5g} {:10.5g}".format(E, k), file=fh_rrkmfile_rev)

    if not flq:
        print("")
        print("RRKM (rev.) file written.")
        print("sf={:.3f}; Delta={:.2f} eV".format(sf, delta))
        print("E0={:.2f} eV; S*={:.2f} J/K/mol".format(barrier, DS))
        print("Number of points: " + str(len(energ)))

    if flp:
        plt.plot(energ / conv, rate, label="Reverse")

    fh_rrkmfile_rev.close()

if flp:
    plt.yscale("log")

    plt.xlabel(r"$E_\mathrm{int}$ (eV)")
    plt.ylabel("$k$ (s$^{-1}$)")

    if outfile_rev:
        plt.legend()

    plt.show()
