# GasTranSim.jl
[![Build Status](https://github.com/kaarthiksundar/GasTranSim.jl/workflows/CI/badge.svg?branch=master)](https://github.com/kaarthiksundar/GasTranSim.jl/actions?query=workflow%3ACI) 
[![codecov](https://codecov.io/gh/kaarthiksundar/GasTranSim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/kaarthiksundar/GasTranSim.jl)
[![Documentation](https://github.com/lanl-ansi/PowerModels.jl/workflows/Documentation/badge.svg)](https://kaarthiksundar.github.io/GasTranSim.jl/dev/)
[![GitHub commits since tagged version](https://img.shields.io/github/commits-since/kaarthiksundar/GasTranSim.jl/v0.1.0.svg)](https://github.com/kaarthiksundar/GasTranSim.jl/commits/master)


This Julia package implements an explicit staggered grid method to perform forward simulation of the full transient gas flow dynamics in a pipeline network. 
``GasTranSim.jl`` is a registered Julia package and can be installed using

```julia 
using Pkg
Pkg.add("GasTransSim")
```

For usage, please refer to the documentation or the ``examples/`` and ``test/`` directories. The json schemas for the data with the validator is provided in the ``schemas/`` folder. 

## Citation
If you find the package useful in your work, we kindly request that you cite the following paper ([arxiv link](https://arxiv.org/abs/1803.00418)): 

```bibtex
@article{GyryaZlotnik2019,
  title={An explicit staggered-grid method for numerical simulation of large-scale natural gas pipeline networks},
  author={Gyrya, Vitaliy and Zlotnik, Anatoly},
  journal={Applied Mathematical Modelling},
  volume={65},
  pages={34--51},
  year={2019},
  publisher={Elsevier}
}
```
