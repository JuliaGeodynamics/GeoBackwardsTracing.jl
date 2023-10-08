
export ReadAllTimesteps, GetGrid

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
        grid    =   data.x.val[:,1,1], data.y.val[1,:,1], data.z.val[1,1,:]
    else
        error("other methods not yet implemented")
    end

    return grid
end
