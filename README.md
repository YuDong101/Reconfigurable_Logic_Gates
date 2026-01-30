# Reconfigurable Logic Gates and Latches in Noise-Driven Bistable Systems

This repository contains the MATLAB source code for simulating the behavior of noise-driven bistable systems. The code is designed to investigate the operational stability and reliability of reconfigurable logic gates and latches within these systems using the Mean First-Passage Time (MFPT) approach.

**Reference:** *Reconfigurable Logic Gates and Latches in Noise-Driven Bistable Systems: A Mean First-Passage Time Approach* by Zi-yi Chen, Long Jiang, Lu-lu Lu, Qi-ming Pei, and Dong Yu.

---

## 1. File Structure

Ensure the following MATLAB files are in the same directory. The workflow relies on the interaction between the main script and the calculation modules.

* **`main.m`**: The main execution script. It sets up the simulation parameters (such as coupling strengths, noise levels), loads initial conditions, runs the ensemble simulations, and aggregates results for analysis.
* **`trial_time_series_.m`**: Contains the code for generating time series data for the trial simulations. It handles the generation of individual oscillator data.
* **`trial_tot.m`**: The core engine for the trial simulations. It aggregates the trial data and computes the overall synchronization dynamics.
* **`find_xm_only.m`**: A utility function used to find and track the synchronization phase of the system across trials.

> **Note:** The main entry point is `main.m`. Do not rename these files, as the script calls them directly.

---

## 2. Usage Instructions

### 2.1 Prerequisites
* MATLAB (Recommended version: R2018b or later).
* No specific toolboxes are required beyond standard MATLAB functions.

### 2.2 Configuration
In `main.m`, you can modify the system parameters to reproduce different results from the manuscript. **The default configuration provided in this repo is set to simulate the bistable system and analyze noise-driven logic gates.**

* **Coupling Parameters**:
    * `K1`: Pairwise coupling strength (default `0.5`).
    * `K2`: Higher-order coupling strength (default `8`).
* **Noise Parameters**:
    * `D1`: Pairwise interaction noise intensity. The script generates a logarithmic range (`Dmin` to `Dmax`).
    * `D2`: Higher-order interaction noise intensity (default `0`).
* **Simulation Settings**:
    * `N_nod`: Number of oscillators (default `10000`).
    * `xIni`: Initial conditions loaded from data files.

### 2.3 Running the Simulation
1. Open `main.m` in the MATLAB editor.
2. Ensure all related files (`trial_time_series_.m`, `trial_tot.m`, `find_xm_only.m`) are in the same directory as `main.m`.
3. Run the script by clicking the **Run** button or typing `main` in the Command Window.

**What the program does:**
* Simulates a bistable system with noise-driven dynamics.
* Analyzes synchronization phases and operational reliability of logic gates.
* Validates the theoretical criterion for determining reliable parameter ranges in noise-driven systems.

### 2.4 Reproducing Other Results
**Important:** All numerical simulation results can be reproduced by adjusting the parameters in `main.m`.

---
