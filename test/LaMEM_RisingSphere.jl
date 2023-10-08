using LaMEM, GeophysicalModelGenerator

# Create some LaMEM output datasets
function run_LaMEM(;ntime=5)
    # Define grids for Interpolations
    n = 17

    # Setup a 3D LaMEM model
    model  = Model(Grid(nel=(n-1,n-1,n-1), x=[-50,50], y=[-50,50], z=[-40,10]))
    matrix = Phase(ID=0,Name="matrix",eta=1e20,rho=3000);
    sphere = Phase(ID=1,Name="sphere",eta=1e23,rho=3200)
    add_phase!(model, sphere, matrix)

    AddSphere!(model,cen=(0.0,0.0,-20), radius=(10, ))

    # Run LaMEM
    model.Time.nstep_max=ntime          # do 5 timesteps
    model.Output.out_strain_rate=1      # output strainrate tensor
    model.Output.out_vel_gr_tensor=1    # velocity gradient
    model.Output.out_j2_strain_rate=1   # strainrate invariant
    
    run_lamem(model,2)

    return nothing

end
