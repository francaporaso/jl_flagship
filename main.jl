using Printf
using Distributed

const NCORES = 16

addprocs(NCORES)
@everywhere include("radial_profile.jl")


"""
Calcula el perfil stackeado de las lentes seleccionadas.
"""
function radial_profile(RMIN, RMAX, NBINS, 
                        Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max,
                        flag, lensname, tracname,
                        ncores)


    ## split lens catalogue
    Lsplit = array_split(L, ncores)
    tot_num = length(Lsplit)
    
    mass   = zeros(NBINS)
    NHalos = zeros(NBINS)
    MassBall  = 0.0
    HalosBall = 0.0

    for (l,L_l) in enumerate(Lsplit)
        println("Vuelta $l de $tot_num")
        num = size(L_l)[1]
        if num==1
            entrada = @views [RMIN, RMAX, NBINS,
                              L_l[:,2], L_l[:,6], L_l[:,7], L_l[:,8]]

            resmap = partial_profile(entrada...)
        else    
            resmap = pmap(partial_profile, 
                       fill(RMIN,num), fill(RMAX,num), fill(NBINS,num), L_l[:,2], L_l[:,6], L_l[:,7], L_l[:,8])
        end
        
        for res in resmap
            mass      += res[1]
            NHalos    += res[2]
            MassBall  += res[3]
            HalosBall += res[4]
        end
        
    end
    println("Paralelizado terminado!")

    ### ----------------------------------------------------- BALL
    println("Calculando stacking...")
    
    DR = (RMAX-RMIN)/NBINS

    MeanDen = MassBall/(4pi/3 * (2RMAX)^3)
    MeanHalos = HalosBall/(4pi/3 * (2RMAX)^3)
    masscum = cumsum(mass)
    NHalosCum = cumsum(NHalos)
    
    Delta    = zeros(NBINS)
    DeltaCum = zeros(NBINS)
    DenHalos = zeros(NBINS)
    DenHalosCum = zeros(NBINS)
    for k in 1:NBINS
        Vol = (4pi/3) * ((k*DR + RMIN)^3 - ((k-1.0)*DR + RMIN)^3)
        Delta[k] = (mass[k]/Vol)/MeanDen - 1.0
        DenHalos[k] = (NHalos[k]/Vol)/MeanHalos

        Vol = (4pi/3) * (k*DR + RMIN)^3
        DeltaCum[k] = ((masscum[k])/Vol)/MeanDen - 1.0
        DenHalosCum[k] = ((NHalosCum[k])/Vol)/MeanHalos
    end
    
    ### ----------------------------------------------------- SHELL
    # Delta = zeros(NBINS)
    # for k in 1:NBINS
    #     V = (4pi/3) * ((k*DR + RMIN)^3 - ((k-1.0)*DR + RMIN)^3)

    #     Delta[k] = (mass[k] * VolShell)/(V * MassShell)
    # end

    println("Done!")

    println("......................")
    println("Saving...")    

    open("pru_stack_par.csv", "w") do io 
        writedlm(io, [Delta DeltaCum DenHalos DenHalosCum], ',')
    end

    # open("pru_individual.csv", "w") do io 
    #     writedlm(io, [Rho RhoCum], ',')
    # end

    # open("pru_meanden.csv", "w") do io 
    #     writedlm(io, MeanDen, ',')
    # end

    println("Done!")
    println("......................")
end


# ------------------------------------------------ #
#                                             main #
@everywhere begin  
    RMIN, RMAX, NBINS = 0.05, 5.0, 100
    Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 10.0, 12.0, 0.2, 0.25, -1.0, -0.8, -1.0, 100.0, 2
    lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
    tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice_halos_cut.fits"
    # lensname = "/mnt/simulations/MICE/voids_MICE.dat"
    # tracname = "/home/fcaporaso/cats/MICE/mice_halos_centralesF.fits"
end

## reading cats
println("......................")
println("Lenses with:")

println("Rvmin......: $Rv_min Mpc")
println("Rvmax......: $Rv_max Mpc")
println("zmin.......: $z_min")
println("zmax.......: $z_max")
println("rho1min....: $rho1_min")
println("rho1max....: $rho1_max")
println("rho2min....: $rho2_min")
println("rho2max....: $rho2_max")

@everywhere begin
    
    L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag=flag, lensname=lensname)
    S = traccat_load(z_min, z_max, tracname=tracname)
    
    nvoids = size(L)[1]
    if nvoids < ncores
        ncores = nvoids
    end
end

println("......................")
println("Calculating profile for:")

println("Nvoids.....: $nvoids")
println("RMIN.......: $RMIN")
println("RMAX.......: $RMAX")
println("NBINS......: $NBINS")
println("NCORES.....: $ncores")

t = @elapsed begin
    radial_profile(RMIN, RMAX, NBINS, Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag, lensname, tracname, ncores)
    # calculate_all_indiv(RMIN, RMAX, NBINS, Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag, lensname, tracname)
end

t /= 60.0
@printf("Ended in %.2f minutes \n", t)
