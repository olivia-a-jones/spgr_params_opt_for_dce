# spgr_params_opt_for_dce

# Interactive SPGR signal analysis for optimising parameters for a brain DCE-MRI scan

This repository contains an interactive tool to visualise how modifying scan parameters (repitition time, flip angle) change the signal and T1-sensitivity of a SPGR scan, with a view to selecting optimised parameters for detecting subtle blood-brain barrier leakage. T10 (aka target tissue) and concentration (amount of BBB leakage) can also be changed to assess how the signal will change in the context of BBB leakage. 

## Interactive Plots
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/olivia-a-jones/spgr_params_opt_for_dce/main?labpath=spgr_interactive_plot.ipynb)

Click the Binder badge above to launch the visualization without downloading code.

## Features
- Interactive sliders for TR, T1₀, and maximum contrast agent concentration
- Comprehensive visualization including:
  - Signal and R1 sensitivity vs flip angle
  - Signal vs concentration for different flip angles
  - Optimal angles vs TR for brain tissues
  - Heatmaps for Ernst angle and optimal R1 sensitivity angle

## Files
- `spgr_interactive.py`: SPGR signal equation calculations and plotting
- `spgr_interactive_plot.ipynb`: Jupyter notebook with interactive plots
