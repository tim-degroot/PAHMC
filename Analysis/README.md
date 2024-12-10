# Analysis

An analysis toolkit has been developped to quickly analyse and visualise the data generated using the `RRRKM` program and `PAHMC` model.

This documentation is a major work in progress and cannot be relied on at this moment.

## Usage

Start by importing the analysis module and iinitialising some parameters that will be shared by all plots. You might have to add the path to the module to the path as well.

```python
sys.path.append('/path/to/PAHMC/Analysis/')
import analysis
settings = analysis.Settings(figsize=(5, 4), dpi=150)
```

### RRKM

Initialise a dataset using the `analysis.RRKM(folder)` class which should point to a folder containg all the `RRRKM` output `.txt` files.

```python
molecule = analysis.RRKM(folder="path/to/folder/")
```

#### RRKM `plot`

The `plot` function generates a plot for a specified subset of data files. It allows for optional filtering, setting axis limits, and customizing the plot title.

**Parameters:**

- `subset` (str): The key to identify the subset of data files to plot.
- `filter_function` (callable, optional): A function to filter the filenames in the subset. Defaults to `None`.
- `title` (str, optional): The title of the plot. Defaults to an empty string.
- `ylim` (list, optional): A list specifying the y-axis limits. Defaults to an empty list.
- `xlim` (list, optional): A list specifying the x-axis limits. Defaults to `[0.0, 5.2]`.

**Example usage:**

```python
molecule.plot(subset="H_diss", filter_function=lambda x: 'a' not in x, title="Hydrogen Dissociation", ylim=[0, 1e5], xlim=[0.0, 5.2])
```

#### RRKM `plot_pairs`

The `plot_pairs` function generates plots for pairs of data files, which are defined as the same reaction for another isotope (i.e. H1to2.txt and D1to2.txt). It allows for setting axis limits and customizing the plot title.

- `subset` (str): The key to identify the subset of data files to plot.
- `title` (str, optional): The title of the plot. Defaults to an empty string.
- `ylim` (list, optional): A list specifying the y-axis limits. Defaults to an empty list.
- `xlim` (list, optional): A list specifying the x-axis limits. Defaults to `[0.0, 5.2]`.

**Example usage:**

```python
molecule.plot_pairs(subset="H_pairs", title="Hydrogen Pairs", ylim=[0, 1e5], xlim=[0.0, 5.2])
```

### Monte Carlo

Initialise a dataset using the `analysis.MonteCarlo` class which should point to a folder containg all the `PAHMC` output and log files.

```python
anthracene = analysis.MonteCarlo(folder="../../Data/Anthracene/Monte_Carlo/", simulations=["1-Anthracene", "2-Anthracene"], numbering="(5,6,7,8) (9,10) (1,2,3,4)", name="D$^+$Anthracene-H$_{10}$")
```

Functions:

- `molecule.data_tables()` prints several calculated data points
- `molecule.mass_spectra(mass)` plots a simulated mass spectra
- `molecule.all_plots("^D")` plots all available plots with regex filter `"^D"`
