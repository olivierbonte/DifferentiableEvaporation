# DifferentiableEvaporation

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://olivierbonte.github.io/DifferentiableEvaporation/dev/)

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> DifferentiableEvaporation

It is authored by Olivier Bonte.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are not included in the git-history and need to be downloaded independently (for more info on how to do so, the reader is referred to `scripts/data_download/README.md`). This is done easily with git:
   ```
   git clone https://github.com/olivierbonte/DifferentiableEvaporation
   ```
1. Open a Julia console, navigate to the folder DifferentiableEvaporation, and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> using DrWatson
   julia> @quickactivate "DifferentiableEvaporation" 
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths. To avoid dependency conflicts, it is recommend to use the Julia version given in the `compat` section of `Project.toml`. 

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "DifferentiableEvaporation"
```
which auto-activate the project and enable local path handling from DrWatson.

Besides Julia, also Python is used in this project. The required dependencies to run the python code is given in `environment.yml`. Install via [(Mini)conda](https://docs.anaconda.com/miniconda/) or [Mamba](https://mamba.readthedocs.io/en/latest/) (in this case, replace `conda` by `mamba` in the command below) in the command line interface:
```
conda env create -f environment.yml
conda activate DifferentiableEvaporation
```

## Data processing
For more information on the data processing, follow the instructions in the folders in the order given below:
1. Data download information, see [here](scripts/data_download/README.md)
2. Data processing information, see [here](scripts/data_processing/README.md)