# Advection routines
export advect_backwards

"""
    particles, p_tag, p_Strain, p_tcrit = advect_backwards(particles::Particles, p_tag, FileName, DirName="", t_start=nothing; fac=0.1, Code="LaMEM", t_end=0.0, verbose=false, strain_max = 1)

Advect tracers backwards with the flow, starting with time `t_start` until `t_end`. Also computes the time, `p_tcrit`, at which particles reach the strain `strain_max`. 

Input:
- `particles`: particles array
- `p_tag`: array with the tag (number) of every particle
- `FileName`: filename of the dataset
- `DirName`: directory of the dataset. If not specified we take the current directory
- `t_start`: starting time. If not specified we take the last timestep 
- `fac`: timestep CFL factor 
- `Code`: Code to be used
- `t_end`: time at which we stop advecting
- `verbose`: display info if `true`
- `strain_max`: maximum strain

Output:
- `particles`: advected particles @ end
- `p_tag`: advected tags
- `p_Strain`: finite strain at time = `t_end`
- `p_tcrit`: time at which the particle exceeded the strain `strain_max`. This can be used during forward advection as starting point.

"""
function advect_backwards(particles::Particles, p_tag, FileName, DirName="", t_start=nothing; fac=0.1, Code="LaMEM", t_end=0.0, verbose=false, strain_max = 1)

    # Read all timesteps:
    Timesteps, Filenames, time = ReadAllTimesteps(FileName, DirName)
    dt_numerical = minimum(diff(time))

    if isnothing(t_start)
        t_start = time[end]
    end

    # Define fields to be tracked on particles
    p_J2, p_Strain, p_Strain_new, p_tcrit = init_cell_arrays(particles, Val(4));

    SecMyr = 3600*24*365.25*1e6

    # Start velocity
    fields  = ("velocity","j2_strain_rate")
    grid, V, J2   = ReadTimestep(FileName, DirName, t_start, fields=fields; Code=Code)
    t       = t_start

    while t>t_end

        # Velocity is given @ vertexes in cm/year in LaMEM
        dx,dy,dz    =   minimum.(diff.(grid))   # in km
        dt          =   -fac*min(dx / maximum(abs.(V[1])),  dy / maximum(abs.(V[2])),  dz / maximum(abs.(V[3]))) # time step
        dt          =   minimum([dt_numerical, dt])
        if (t+dt)<t_end
            dt = t_end - t
        end

        # Advect backwards with 2nd order Runga Kutta
        advection_RK!( particles, V, grid, grid, grid, dt, 2 / 3 )
    
        # Interpolate strainrate invariant to particles
        grid2particle!(p_J2, grid, J2, particles.coords)
    
        # Update accumulated strain on particles (note that it is negative)
        p_Strain .= p_Strain_new
        
        p_Strain_new .= p_Strain + p_J2.*dt.*SecMyr

        # check whether we exceed a critical strain & record the time when that happens
        for I in CartesianIndices(p_Strain_new)
            N = length(p_Strain_new[I])
            p_Strain_new_cell = p_Strain_new[I]
            p_Strain_cell   = p_Strain[I]
            p_tcrit_cell = p_tcrit[I]
            cell = fill(NaN,N)

      
            for k=1:N
                if !isnan(p_Strain_new_cell[k])
                    cell[k] = p_tcrit_cell[k]
                    if ((p_Strain_new_cell[k]<= -strain_max) &  (p_Strain_cell[k] > -strain_max))
                        # estimate exact time that this happens
                        Δ_Strain = p_Strain_cell[k]-p_Strain_new_cell[k]
                        f = (p_Strain_cell[k]+strain_max)/Δ_Strain
                        cell[k] = t + dt*f
                    end
                    p_tcrit[I] = SVector{N}(cell)
                end
            end
        end

        
        # redistribute in memory
        shuffle_particles!(particles, grid, (p_tag, p_J2, p_Strain, p_Strain_new, p_tcrit))
        
        if verbose
            @show t, dt
        end
        t += dt
    end

    return particles, p_tag, p_Strain, p_tcrit
end