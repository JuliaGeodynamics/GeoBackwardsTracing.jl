using Test, GeoBackwardsTracing

@test 1==1
include("LaMEM_RisingSphere.jl")

# perform a simple LaMEM simulation
run_LaMEM(ntime=5)

# 
FileName = "output"
DirName  = pwd()

# Read all timesteps:
Timesteps, Filenames, time = ReadAllTimesteps(FileName, DirName)

t_start = time[end]

# Read the 1D coordinate vectors of the full grid
grid = GetGrid(FileName, DirName, Timesteps[end])

grid_p = grid   # where are particles located (can be subset of grid)?


@testset "Backwards advection" begin

    particles, p_tag = init_particles(grid_p);  # init particles @ center of cell

    # advect backwards and record when a strain of 0.01 is reached
    particles, p_tag, p_Strain, p_tcrit = advect_backwards(particles, p_tag, FileName, DirName, t_start; 
                                                            fac=0.1, Code="LaMEM", t_end=0.0, strain_max = 0.01);

    @test p_tag[100][1] == 100.0
    @test p_Strain[100][1] ≈ -0.029795379278230134
    @test p_tcrit[100][1] ≈ 0.18976066870597108


    particles, p_tag = init_particles(grid_p);  # init particles @ center of cell

    # advect backwards and record when a strain of 0.01 is reached
    particles, p_tag, p_Strain, p_tcrit = advect_backwards(particles, p_tag, FileName, DirName, t_start; 
    fac=0.1, Code="LaMEM", t_end=0.0, strain_max = 0.01, only_initial_velocity=true);

    @test p_tag[100][1] == 100.0
    @test p_Strain[100][1] ≈ -0.030579771444418327
    @test p_tcrit[100][1] ≈ 0.18976066870597108

end
