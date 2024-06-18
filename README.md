# DifferentiableEvaporation

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> DifferentiableEvaporation

It is authored by Olivier Bonte.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently. This is done easily with git:
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
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "DifferentiableEvaporation"
```
which auto-activate the project and enable local path handling from DrWatson.
