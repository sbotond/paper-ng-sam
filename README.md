# Pipeline for simulating NG-SAM experiments

This repository contains the simulation pipeline used in the manuscript:

Botond Sipos, Tim Massingham, Adrian M. Stütz, Nick Goldman: *An improved protocol for sequencing of repetitive genomic regions and structural variations using mutagenesis and Next Generation Sequencing*

## Files

* **bin** - scripts:
    * calibrate_mut - script calculating the branch length scaling factor
    * pcr_coal.R - R script simulating PCR amplifications using [pcrcoal](https://github.com/sbotond/pcrcoal) and dilutions by sampling from [Poisson distributions](http://en.wikipedia.org/wiki/Poisson_distribution)
    * sim_exp - simulate a single NG-SAM experiment with the specified target sequence and parameters
    * run_seq_sim - simulate NG-SAM experiments on different target sequences
    * run_dil_sim - simulate NG-SAM experiments with a range of dilution factors
    * plot_dil_res - plot the results of seq_sim
    * plot_seq_res - plot the results of dil_sim
* **dat** - data files:
    * bl_scaler.txt - the branch length scaling factor and the corresponding Hamming distance as calculated by the "calibrate_mut" script
    * dmel_eater.fas - the coding sequence of the *Drosophila melanogaster* *eater* gene
    * eater_root.fas - the target sequence used in the second simulation setup (run_dil_sim)
    * MH22.fas - sequence used to calibrate the branch length scaling factor from <a href="http://www.ncbi.nlm.nih.gov/pubmed/8568899">Zaccolo et al.</a>
    * mutation_model.tab - the mutation spectrum observed by <a href="http://www.ncbi.nlm.nih.gov/pubmed/8568899">Zaccolo et al.</a>, used to build the mutation model
    * s_8_4x.runfile - the [simNGS](http://www.ebi.ac.uk/goldman-srv/simNGS/) runfile used in the simulations
* **lib/*.py** - python classes used by the scripts under bin/
* **seq_sim** - output directory for the first simulation setup
* **dil_sim** - output directory for the second simulation setup
* **Makefile** - makefile containing utility targets
* **simulations.mk** - makefile containing simulation targets and parameters
* **reports** - plots:
    * calibration_report.pdf - diagnostic plots from the "calibrate_mut" script
    * seq_sim.pdf - the results of the first simulation setup
    * dil_sim.pdf - the results of the second simulation setup

## Requirements

The simulation pipeline runs in a standard UNIX environment and uses the Platform LSF workload manager to distribute simulations between multiple compute nodes. It also requires the following software to be installed:

* [R](http://www.r-project.org/) (>= 2.14.1) with the [pcrcoal](http://cran.r-project.org/web/packages/pcrcoal) package installed.
* [python](http://www.python.org/) (>= 2.7.1) with the following non-standard packages:
    * [Biopython](http://pypi.python.org/pypi/biopython/) (>= 1.59)
    * [DendroPy](http://pypi.python.org/pypi/DendroPy/) (>= 3.11.0)
    * [numpy](http://pypi.python.org/pypi/numpy/) (>= 1.6.1)
    * [matplotlib](http://pypi.python.org/pypi/matplotlib/) (>= 1.1.0)
* [exonerate](http://www.ebi.ac.uk/~guy/exonerate/) (>= 2.12.3)
* [muscle](http://www.drive5.com/muscle/) (>= 3.8.31)
* [simNGS](http://www.ebi.ac.uk/goldman-srv/simNGS/) (1.5.1)
* [velvet](https://github.com/dzerbino/velvet) (latest version)

## Running simulations

The simulation parameters of interest are stored in the makefile simulations.mk.
The simulations and the plotting scripts are launched by calling the following make targets:

* **seq_sim** - Submit the jobs for the first simulation setup. The results and random target sequences are saved under "seq_sim".
* **plot_seq_res** - process the output of seq_sim
* **dil_sim** - Submit the jobs for the first simulation setup. The results are saved under "dil_sim".
* **plot_seq_res** - process the output of dil_sim

Other useful make targets:

* **t** - test the simulation framework
* **calibration** - recalculate the branch length scaling factor and save it in dat/bl_scaler.txt
