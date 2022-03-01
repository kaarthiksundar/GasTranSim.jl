# GasTranSim.jl

This Julia package implements an explicit staggered grid method to perform forward simulation of the full transient gas flow dynamics in a pipeline network. 
``GasTranSim.jl`` is not a registered Julia package. Hence installation of the package should be done as follows:

```julia 
using Pkg
Pkg.add("https://github.com/kaarthiksundar/GasTranSim.jl.git")
```

For usage, please refer to the documentation or the ``examples/`` and ``test/`` directories. The json schemas for the data with the validator is provided in the ``schemas/`` folder. 