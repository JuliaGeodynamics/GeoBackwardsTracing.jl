module GeoBackwardsTracing

using LaMEM, GeophysicalModelGenerator
using JustPIC, StaticArrays

# set_backend("CUDA_Float64_3") # need to restart session if this changes
set_backend("Threads_Float64_3D") # need to restart session if this changes

using CellArrays, ParallelStencil
@init_parallel_stencil(Threads, Float64, 3)

const TA = backend == "CUDA" ? JustPIC.CUDA.CuArray : Array


include("init.jl")  # particle initialisation
include("utils.jl")





end # module GeoBackwardsTracing
