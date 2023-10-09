
export ReadAllTimesteps, GetGrid, ReadTimestep

"""
    Timesteps, Filenames, time = ReadAllTimesteps(FileName="output", DirName=pwd(); Code="LaMEM")

Reads all timesteps of filename `FileName` in directory `DirName`. By default, this is done for LaMEM; other codes can be added here
"""
function ReadAllTimesteps(FileName="output", DirName=pwd(); Code="LaMEM")
    if Code=="LaMEM"
        Timesteps, Filenames, time = Read_LaMEM_simulation(FileName,DirName)
    else
        error("other methods not yet implemented")
    end

    return Timesteps, Filenames, time
end

"""
    grid = GetGrid(FileName, DirName, starting_timestep; Code="LaMEM")

Retrieves the 1D coordinate vectors that describe the cartesian grid of the simulation
"""
function GetGrid(FileName, DirName, starting_timestep; Code="LaMEM")

    if Code=="LaMEM"
        data, t =   Read_LaMEM_timestep(FileName, starting_timestep, DirName)
        grid    =   GetGrid(data)
    else
        error("other methods not yet implemented")
    end

    return grid
end

function GetGrid(data::CartData)

    grid    =   data.x.val[:,1,1], data.y.val[1,:,1], data.z.val[1,1,:]

    return grid
end



"""

    grid, data_fields... = ReadTimestep(FileName, DirName, time, fields=("velocity","j2_strain_rate"); Code="LaMEM")

Interpolate the data at time `time` from a series of timesteps (read from disk). Tuples of fields are returned and if the time is inbetween 2 timesteps the data fields are interpolated
"""
function ReadTimestep(FileName, DirName, time; fields=("velocity","j2_strain_rate"), Code="LaMEM")
    
    Timesteps, Filenames, time_vec = ReadAllTimesteps(FileName, DirName; Code=Code)

    if time<minimum(time_vec)
        error("time less than minimum time")
    elseif time>maximum(time_vec)
        error("time more than maximum time")
    end

    # interpolate the timesteps
    id0 = findlast( (time_vec .- time).<= 0)
    id1 = findfirst((time_vec .- time).>= 0)
    
    dt  = time_vec[id1] - time_vec[id0]
    if dt>0
        fac = (time-time_vec[id0])/dt
    else
        fac = 1.0
    end 


    if Code=="LaMEM"
        # read 2 surrounding timesteps
        data0, _ = Read_LaMEM_timestep(FileName, Timesteps[id0], DirName, fields=fields)
        data1, _ = Read_LaMEM_timestep(FileName, Timesteps[id1], DirName, fields=fields)

        # interpolate & create a new field
        data_interp = CartData(data0.x.val, data0.y.val, data0.z.val, (z=data0.z.val,))
        names = keys(data0.fields)
        for i=1:length(data0.fields)
            d0 = data0.fields[i]
            d1 = data1.fields[i]
            d_average = fac.*d0 .+ (1.0-fac).*d1

            # scale
            cmYr_to_kmMyr = 1e-5*1e6
            if (fields[i] == "velocity")
                # scale to km/Myr
                d_average = d_average.*cmYr_to_kmMyr
            end

            data_interp = AddField(data_interp,String(names[i]),d_average) 
        end
        data_interp = RemoveField(data_interp,"z") 

    else
        error("other codes not yet implemented")
    end
    grid    =   GetGrid(data_interp)
    
    return grid, data_interp.fields...
end