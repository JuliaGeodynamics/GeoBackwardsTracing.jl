# GeoBackwardsTracing

:warning:**Warning**:warning: This package is still under development and the API is not stable yet. 

In geodynamic simulations, it is often useful to perform postprocessing steps to compute, for example, *P-T* evolution or finite strain tensors at specific points at the end of a simulation. Whereas this can be done during a simulation (by adding passive tracers), it is often better for visualisations to do this at the end of a simulation and trace points backwards in time.

This package starts with a specific timestep and advects passive tracers backwards in time, using timesteps saved on disk. In a next step, tracers can again be advected forward and parameters such as finite strain can be computed. It works along with the [GeophysicalModelGenerator](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl) `CartData` file format and was developed to work with LaMEM simulations and is mostly a package build on [JustPIC](https://github.com/JuliaGeodynamics/JustPIC.jl). It should be possible to adapt it to other numerical codes with minor changes, provided they use the `CartData` format.

