
export init_particles

"""
    particles = init_particles(grid, nxcell, max_xcell, min_xcell; random=true)

This initializes particles on the `grid` (consisting of 1D vectors)
"""
function init_particles(grid::Tuple, nxcell::Int64=1, max_xcell::Int64=10, min_xcell::Int64=1; random=false)
    x,y,z = grid
    ni    = length.(grid) .- 1

    ncells = prod(ni)
    np = max_xcell * ncells
    px, py, pz = ntuple(_ -> @fill(NaN, ni..., celldims=(max_xcell,)) , Val(3))

    inject = @fill(false, ni..., eltype=Bool)
    index = @fill(false, ni..., celldims=(max_xcell,), eltype=Bool) 
    
    @parallel_indices (i, j, k) function fill_coords_index(px, py, pz, index)
        # lower-left corner of the cell
        x0, y0, z0 = (x[i+1]-x[i])/2, (y[j+1]+y[j])/2, (z[k+1]+z[k])/2
        dx, dy, dz = x[i+1]-x[i], y[j+1]-y[j], z[k+1]-z[k]
        
        # fill index array
        for l in 1:nxcell
            @cell    px[l, i, j, k] = x0 + dx * rand(-0.49:1e-5:0.49) * Float64(random)
            @cell    py[l, i, j, k] = y0 + dy * rand(-0.49:1e-5:0.49) * Float64(random)
            @cell    pz[l, i, j, k] = z0 + dz * rand(-0.49:1e-5:0.49) * Float64(random)
            @cell index[l, i, j, k] = true
        end
        return nothing
    end

    @parallel (JustPIC.@idx ni) fill_coords_index(px, py, pz, index)    

    particles = Particles(
        (px, py, pz), index, inject, nxcell, max_xcell, min_xcell, np, ni
    )
    
    # add tags to particles
    p_tag,  = init_cell_arrays(particles, Val(1));

    add_particle_tags!(p_tag, particles)

    return particles, p_tag
end


"""
    add_particle_tags!(particle_tag, particles)
Add tags (numbers) to particles
"""
function add_particle_tags!(particle_tag, particles)

    counter = 0
    for i in eachindex(particles.index.data)
        !(particles.index.data[i]) && continue 
        counter += 1
        particle_tag.data[i] = counter
    end

    return nothing
end