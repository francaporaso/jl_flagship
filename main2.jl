using Distributed

NCORES = 4

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
    Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 10.0, 12.0, 0.2, 0.25, -1.0f0, -0.9f0, -1.0, 100.0f0, 2
    lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
    tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice_halos_cut.fits"
end


println("Cargando catálogos")
@everywhere begin
    S = traccat_load(z_min, z_max, tracname=tracname)
    L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag=flag, lensname=lensname)
    nvoids = size(L)[1]
    i = 1:nvoids
end

println("NVOIDS: $nvoids")
println("Corriendo en paralelo....")

resmap = pmap(partial_profile, fill(RMIN,nvoids), fill(RMAX,nvoids), fill(NBINS,nvoids), L[i,2], L[i,6], L[i,7], L[i,8])

println("Hecho!")

println("Calculando stacking...")

mass_p   = zeros(NBINS, nvoids)
NHalos_p = zeros(NBINS, nvoids)
MassBall_p  = zeros(nvoids)
HalosBall_p = zeros(nvoids)

for l in 1:nvoids
    for j in eachindex(resmap)
        mass_p[:,l]    += resmap[j][1]
        NHalos_p[:,l]  += resmap[j][2]
        MassBall_p[l]  += resmap[j][3]
        HalosBall_p[l] += resmap[j][4]
    end
end

mass = vec(sum(mass_p, dims=2))
masscum = cumsum(mass)
NHalos = vec(sum(NHalos_p, dims=2))
NHalosCum = cumsum(NHalos)
MassBall = sum(MassBall_p)
HalosBall = sum(HalosBall_p)

DR = (RMAX-RMIN)/NBINS

MeanDen = MassBall/(4pi/3 * (2RMAX)^3)
MeanHalos = HalosBall/(4pi/3 * (2RMAX)^3)

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

println("Saving...")    

open("pru_stack.csv", "w") do io 
    writedlm(io, [Delta DeltaCum DenHalos DenHalosCum], ',')
end

println("Done!")