# Integrated NMR Structure Determination and MD Simulation Workflow

This repository contains the computational workflow used for NMR structure determination, molecular dynamics (MD) simulations, and validation of structural ensembles through back-calculation of experimental observables.

The workflow integrates:

- XPLOR-NIH structure calculations
- Molecular dynamics simulations
- Chemical shift (CS) back-calculation
- Dipolar coupling (DC/RDC) analysis
- Structural refinement and validation
- Deposited structural ensembles

---

## Overview

This project provides a reproducible pipeline for:

1. Preparing NMR-derived restraints
2. Calculating structures using XPLOR-NIH
3. Refining structures with MD simulations
4. Back-calculating experimental observables
5. Comparing simulations with experimental NMR data
6. Preparing deposited structural ensembles

---

## Methods

### Structure Determination
Structures were calculated using:

- XPLOR-NIH
- Simulated annealing protocols
- Experimental restraints including:
  - NOEs
  - HBond restraints
  - Dihedral restraints
  - RDC/DC restraints

### Molecular Dynamics Simulations
MD simulations were performed using:

- NAMD
- Explicit solvent and membrane environments

### Back-Calculation and Validation
Experimental observables were validated using:

- Chemical shift back-calculation
- Dipolar coupling/RDC calculations
- Ensemble averaging
- Structural clustering and analysis

---

## Repository Structure

```bash
.
├── XPLOR/                 # XPLOR-NIH structure calculation scripts
├── MD/                    # Molecular dynamics simulation inputs
├── Analysis/              # Trajectory analysis scripts
├── CS_BackCalculation/    # Chemical shift back-calculation
├── RDC_BackCalculation/   # RDC/DC analysis
├── Structures/            # Final structural ensembles
├── Figures/               # Manuscript figures
└── Scripts/               # script for automatic generation of psf file from sequence
```

---

## Software Requirements

### Structure Calculation
- XPLOR-NIH
- VMD

### Molecular Dynamics
- NAMD 
- CUDA-enabled GPUs (recommended)

### Analysis
- cpptraj
- Python ≥ 3.9
- NumPy
- pandas
- matplotlib
- MDAnalysis (optional)

---

## Example Workflow

### 1. Run XPLOR-NIH Structure Calculation

```bash
python run_xplor.py
```

### 2. Run MD Simulation

```bash
pmemd.cuda -O -i mdin -p system.parm7 -c equil.rst7 \
-o md.out -r md.rst7 -x md.nc
```

### 3. Perform Trajectory Analysis

```bash
cpptraj -i analysis.in
```

### 4. Back-Calculate Experimental Observables

```bash
python backcalculate_cs.py
python backcalculate_rdc.py
```

---

## Data Availability

Final structural ensembles and associated data have been deposited in public repositories (protein data bank, PDB).

---

## Citation

If you use this workflow or repository, please cite:

```text
[Authors], [Title], [Journal], [Year]
```

---

