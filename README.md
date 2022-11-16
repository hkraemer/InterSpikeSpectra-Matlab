# InterSpikeSpectra for Matlab

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hkraemer.github.io/InterSpikeSpectra-Matlab/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hkraemer.github.io/InterSpikeSpectra-Matlab/dev)
[![Build Status](https://github.com/hkraemer/InterSpikeSpectra-Matlab/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hkraemer/InterSpikeSpectra.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/hkraemer/InterSpikeSpectra-Matlab/branch/main/graph/badge.svg)](https://codecov.io/gh/hkraemer/InterSpikeSpectra-Matlab)

A MATLAB toolbox for obtaining inter-spike spectra from signals. It is recommended to analyze "spiky" signals with 
this method. As the authors showed in the corresponding paper (Kraemer et al. 2022, Spike spectra for recurrences) 
this method can yield reasonable results and insights into $\tau$-recurrence rate signals obtained from 
recurrence plots of the underlying signal.

# Installation
Simply double-click `InterSpikeSpectra-Matlab.mltbx` for installing the toolbox.

# Functionality
The main function to call is `inter_spike_spectrum()`, simply type 
```matlab
help inter_spike_spectrum
```
into you MATLAB-IDE and read the documenation including an example how to use it.