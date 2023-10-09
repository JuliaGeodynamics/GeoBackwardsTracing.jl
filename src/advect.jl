# Advection routines
export advect_backwards

"""

Advect tracers backwards with the flow, starting with time `t_start` until `t_end`.
We 
"""
function advect_backwards(particles::Particles, p_tag, FileName, DirName, t_start; fac=0.1, Code="LaMEM", t_end=0.0)

    # Read all timesteps:
    Timesteps, Filenames, time = ReadAllTimesteps(FileName, DirName)

    cmYr_to_kmMyr = 1e-5*1e6
    SecMyr = 3600*24*365.25*1e6

    # Define fields to be tracked on particles
    p_J2, p_Strain = init_cell_arrays(particles, Val(2));

    # Start velocity
    fields  = ("velocity","j2_strain_rate")
    data    = ReadVelocity(FileName, DirName, t_start, fields=fields; Code=Code)
    t       = t_start

    while t>t_end

        # Velocity is given @ vertexes in cm/year
        Vx, Vy, Vz  =   data.fields.velocity[1], data.fields.velocity[2], data.fields.velocity[3];
        V           =   Float64.(Vx.*cmYr_to_kmMyr), Float64.(Vy.*cmYr_to_kmMyr), Float64.(Vz.*cmYr_to_kmMyr);
        J2          =   data.fields.j2_strain_rate
        grid        =   GetGrid(data)
        dx,dy,dz    =   minimum.(diff.(grid))
        dt          =   -fac*min(dx / maximum(abs.(V[1])),  dy / maximum(abs.(V[2])),  dz / maximum(abs.(V[3]))) # time step
        if (t+dt)<t_end
            dt = t_end - t
        end

        # Advect backwards with 2nd order Runga Kutta
        advection_RK!( particles, V, grid, grid, grid, dt, 2 / 3 )
    
        # Interpolate strainrate invariant to particles
        grid2particle!(p_J2, grid, data.fields.j2_strain_rate, particles.coords)
    
        # Update accumulated strain on particles
        p_Strain .= p_Strain + p_J2.*dt.*SecMyr
    
        # redistribute in memory
        shuffle_particles!(particles, grid, (p_tag, p_J2, p_Strain))
    
        t += dt
        @show t
    end

    return particles, p_tag, p_Strain
end