function main(resmap)
    masa = zeros(NBINS)
    nhalos = zeros(NBINS)
    massball = 0.0
    halosball = 0.0

    for res in resmap
        masa      += res[1]
        nhalos    += res[2]
        massball  += res[3]
        halosball += res[4]
    end

    masscum = cumsum(masa)
    nhaloscum = cumsum(nhalos)

    DR = (RMAX-RMIN)/NBINS

    MeanDen = massball/(4pi/3 * (2RMAX)^3)
    MeanHalos = halosball/(4pi/3 * (2RMAX)^3)

    Delta    = zeros(NBINS)
    DeltaCum = zeros(NBINS)
    DenHalos = zeros(NBINS)
    DenHalosCum = zeros(NBINS)
    for k in 1:NBINS
        Vol = (4pi/3) * ((k*DR + RMIN)^3 - ((k-1.0)*DR + RMIN)^3)
        Delta[k] = (masa[k]/Vol)/MeanDen - 1.0
        DenHalos[k] = (nhalos[k]/Vol)/MeanHalos

        Vol = (4pi/3) * (k*DR + RMIN)^3
        DeltaCum[k] = ((masscum[k])/Vol)/MeanDen - 1.0
        DenHalosCum[k] = ((nhaloscum[k])/Vol)/MeanHalos
    end


    println("Saving...")    

    open("pru_stack_par.csv", "w") do io 
        writedlm(io, [Delta DeltaCum DenHalos DenHalosCum], ',')
    end
end

## ------------------------------------------ MAIN -------------

using Distributed

NCORES = 8

addprocs(NCORES-1)
nc = nprocs()
println("NCORES: $nc")

@everywhere begin 
    include("radial_profile.jl")
    using DelimitedFiles, FITSIO, DataFrames
    using Statistics
    # using Base.Threads
    #using ProgressMeter
    
    include("tools.jl")
end

@everywhere begin 
    RMIN, RMAX, NBINS = 0.05, 5.0, 100
    Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 10.0, 12.0, 0.2, 0.25, -1.0, -0.9, -1.0, 100.0, 2
    # lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
    # tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice_halos_cut.fits"
    lensname = "/mnt/simulations/MICE/voids_MICE.dat"
    tracname = "/home/fcaporaso/cats/MICE/mice_halos_centralesF.fits"
end

println("Cargando catálogos")
@everywhere begin
    S = traccat_load(z_min, z_max, tracname=tracname, lmhalo_min=1)
    L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag=flag, lensname=lensname)
    nvoids = size(L)[1]
    i = 1:nvoids
end

println("NVOIDS: $nvoids")
println("Corriendo en paralelo....")

t = @elapsed begin
    resmap = pmap(partial_profile, fill(RMIN,nvoids), fill(RMAX,nvoids), fill(NBINS,nvoids), L[i,2], L[i,6], L[i,7], L[i,8], batch_size=NCORES)    
end

t /= 60.0
println("Hecho!")
println("Termiando en $t mins")

println("Calculando stacking...")

main(resmap)

println("Done!")