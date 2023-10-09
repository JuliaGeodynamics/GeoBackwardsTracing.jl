# Advection routines
export advect_backwards

"""

Advect tracers backwards with the flow, starting with time `t_start` until `t_end`.
We 
"""
function advect_backwards(particles::Particles, p_tag, FileName, DirName, t_start; fac=0.1, Code="LaMEM", t_end=0.0, verbose=false)

    # Read all timesteps:
    Timesteps, Filenames, time = ReadAllTimesteps(FileName, DirName)
    dt_numerical = minimum(diff(time))

    # Define fields to be tracked on particles
    p_J2, p_Strain = init_cell_arrays(particles, Val(2));

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
        p_Strain .= p_Strain + p_J2.*dt.*SecMyr
    
        # redistribute in memory
        shuffle_particles!(particles, grid, (p_tag, p_J2, p_Strain))
        
        if verbose
            @show t, dt
        end
        t += dt
    end

    return particles, p_tag, p_Strain
end