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

particles, p_tag = init_particles(grid_p);  # init particles @ center of cell


particles, p_tag, p_Strain = advect_backwards(particles, p_tag, FileName, DirName, t_start; fac=0.1, Code="LaMEM", t_end=0.0)



#=
# Read timesteps back to julia

data, t =   Read_LaMEM_timestep(FileName, 0)
xvi     =   xv, yv, zv = data.x.val[:,1,1], data.y.val[1,:,1], data.z.val[1,1,:]
dxi     =   dx, dy, dz = xv[2] - xv[1], yv[2] - yv[1], zv[2] - zv[1]
ni      =   length.(xvi)
grid    =   xvi;

nxcell, max_xcell, min_xcell = 24, 24, 24
particles = init_particles(
    nxcell, max_xcell, min_xcell, xvi..., dxi..., ni.-1
);

# Define fields to be tracked on particles
p_tag, p_J2, p_Strain = init_cell_arrays(particles, Val(3));

# Add tags
add_particle_tags!(p_tag, particles)

# Define cell array with tensors that will hold the velocity gradient
#celldims = (3,3)
#Cell = SMatrix{celldims..., Float64, prod(celldims)}
#VelocityGradient = CPUCellArray{Cell}(undef, size(particles.coords[1])...);

# make backup of original particles
particles0 = deepcopy(particles);

dt_LaMEM = diff(time)
dt_LaMEM = vcat(dt_LaMEM, dt_LaMEM[end])

# velocities in LaMEM are in cm/yr & coordinates and time in km and Myrs
cmYr_to_kmMyr = 1e-5*1e6
SecMyr = 3600*24*365.25*1e6

data, t = Read_LaMEM_timestep(FileName, 0)

t = 0
fac = 0.1

pxv = particles.coords[1].data;
pyv = particles.coords[2].data;
pzv = particles.coords[3].data;
idxv = particles.index.data;

I = findfirst(p_tag.data .== 10000.0)

p0 = [(pxv[I][1], pyv[I][1], pzv[I][1])]

while t>-0.5
    global t, fac, grid, p_Strain, p_tag, p_J2, VelocityGradient

    # Velocity is given @ vertexes in cm/year
    Vx = data.fields.velocity[1];
    Vy = data.fields.velocity[2];
    Vz = data.fields.velocity[3];
    V  = Float64.(Vx.*cmYr_to_kmMyr), Float64.(Vy.*cmYr_to_kmMyr), Float64.(Vz.*cmYr_to_kmMyr);
    J2 = data.fields.j2_strain_rate
    dt = -fac*min(dx / maximum(abs.(V[1])),  dy / maximum(abs.(V[2])),  dz / maximum(abs.(V[3]))) # time step

    advection_RK!(
        particles, V, grid, grid, grid, dt, 2 / 3
    )

    # Interpolate strainrate invariant to particles
    grid2particle!(p_J2, xvi, data.fields.j2_strain_rate, particles.coords)

    # update velocity gradient tensor

    # Update accumulated strain on particles
    p_Strain .= p_Strain + p_J2.*dt.*SecMyr

    # redistribute in memory
    shuffle_particles!(particles, xvi, (p_tag, p_J2, p_Strain))

    t += dt
    @show t
end

pxv = particles.coords[1].data;
pyv = particles.coords[2].data;
pzv = particles.coords[3].data;
I = findfirst(p_tag.data .== 10000.0)
p = [(pxv[I][1], pyv[I][1], pzv[I][1])]

=#