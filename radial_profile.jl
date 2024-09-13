using DelimitedFiles
using DataFrames
using FITSIO
using Statistics

include("cosmology.jl")

"""
Loads the lenses catalog
"""
function lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max; 
                      flag=2.0, lensname="/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat")
    ## L[1] = id
    ## L[2] = rv
    ## L[5] = z
    ## L[6] = xv
    ## L[7] = yv
    ## L[8] = zv
    ## L[9] = rho1
    ## L[10] = rho2
    ## L[12] = flag
    L = DataFrame(readdlm(lensname), :auto)

    m_rv   = @. (L[!,2] >= Rv_min) && (L[!,2] <= Rv_max)
    m_z    = @. (L[!,5] >= z_min) && (L[!,5] <= z_max)
    m_rho  = @. (L[!,9] >= rho1_min) && (L[!,9] <= rho1_max) && (L[!,10] >= rho2_min) && (L[!,10] <= rho2_max)
    m_flag = @. (L[!,12] >= flag) 

    mask = @. m_rv && m_z && m_rho && m_flag

    return L[mask,:]
end

"""
Dado un cto de centros (xv, yv, zv) con su correspondiente radio rv
encuentra los trazadores alrededor del centro hasta (1+DR)RMAX*rv
"""
function get_tracers(RMAX::Float64, DR::Float64,
                    rv::Vector{Float64}, xv::Vector{Float64}, yv::Vector{Float64}, zv::Vector{Float64}; 
                    tracname="/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits")

    tcat = Matrix(DataFrame(FITS(tracname)[2])) ## FITS lee .fits; DataFrame transforma en tabla
    ## tcat[1] = id
    ## tcat[2] = flagcentral
    ## tcat[3] = lmhalo
    ## tcat[4] = xhalo
    ## tcat[5] = yhalo
    ## tcat[6] = zhalo

    nvoids = length(rv)
    trac_lists = Vector{Matrix{Float64}}(undef, nvoids)

    ### Máscara en una bola con centro (xv,yv,zv) y radio (1+DR)RMAX*rv
    for v in 1:nvoids
        distance = @views @. sqrt((tcat[:,4] - xv[v])^2 + (tcat[:,5] - yv[v])^2 + (tcat[:,6] - zv[v])^2)
        mask = distance .<= ((1.0+DR)*RMAX*rv[v])

        tr = hcat(tcat[mask,3:6], distance[mask])
        trac_lists[v] = tr[sortperm(tr[:,end]),:] ## ordena de menor a mayor distance
    end

    return trac_lists
end

"""
Calcula el perfil de 1 void dados los trazadores tcat al rededor de (xv,yv,zv)
"""
function partial_profile(tcat::Matrix{Float64}, 
                        RMIN, RMAX, DR,
                        rv, z, # redshift
                        xv, yv, zv)
    
    # MeanDen = mean_density(z, 70, 0.25, 0.75)
    MeanDen = critical_density(z, 70, 0.25, 0.75)

    NBINS = length(RMIN:DR:RMAX)
    rmin = rv*RMIN
    rmax = rv*RMAX
    dr   = rv*DR

    NTrac = zeros(NBINS)
    mass = zeros(NBINS)

    rad = RMIN

    for i in 1:NBINS
        tr = @views tcat[(tcat[:,end] .< (rad+dr)) .&& (tcat[:,end] .>= rad), :]
        mass[i] = sum(10.0 .^ tr[:,1])
        volume = (4.0 *pi/ 3.0 )*((rad+dr)^3 - rad^3)
        
        NTrac[i] = (length(tr[:,1])/volume)/MEAN_NTRAC
        rad += dr
    end

    vol = 4pi/3*[(ri+dr)^3 - ri^3 for ri in rmin:dr:rmax]
    volcum = 4pi/3*[(ri+dr)^3 for ri in rmin:dr:rmax]
    Delta = (mass./vol)/MeanDen
    DeltaCum = (cumsum(mass)./volcum)/MeanDen

    return Delta, DeltaCum, NTrac
end

"""
Calcula todos los perfiles de las lentes seleccionadas
"""
function radial_profile(RMIN, RMAX, dr, 
                        Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max,
                        flag, lensname::String, tracname::String)
    ## reading cats
    println("Reading lenses...")
    L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag=flag, lensname=lensname)
    println("Reading tracers...")
    tracers_lists = get_tracers(RMAX, dr, L[!,2], L[!,6], L[!,7], L[!,8], tracname=tracname)
    
    nvoids = nrow(L)
    NBINS = length(RMIN:DR:RMAX)
    println("nvoids: $nvoids")

    println("Calculating profile...")
    
    Delta   = Matrix{Float64}(undef, nvoids, NBINS)
    DeltaCum = Matrix{Float64}(undef, nvoids, NBINS)
    NTracers = Matrix{Float64}(undef, nvoids, NBINS)

    for i in 1:nvoids
        Delta[i,:], DeltaCum[i,:], NTracers[i,:] = partial_profile(tracers_lists[i], RMIN, RMAX, DR, L[i,2], L[i,5], L[i,6], L[i,7], L[i,8])
    end
    
    Delta_stack = sum(Delta, dims=1)'/nvoids
    Delta_std  = std(Delta, dims=1)'/nvoids
    
    DeltaCum_stack = sum(DeltaCum, dims=1)'/nvoids
    DeltaCum_std  = std(DeltaCum, dims=1)'/nvoids

    NTrac_stack = sum(NTracers, dims=1)'/nvoids
    NTrac_std = std(NTracers, dims=1)'/nvoids

    println("Saving...")

    open("pru_stack.csv", "w") do io 
        writedlm(io, [Delta_stack Delta_std DeltaCum_stack DeltaCum_std NTrac_stack NTrac_std], ',')
    end

    open("pru_individual.csv", "w") do io 
        writedlm(io, [Delta DeltaCum NTracers], ',')
    end

    println("End!")
end


# ------------------------------------------------ #
#                                             main #
const RMIN, RMAX, DR = 0.01, 5., 0.05
const Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 10., 15., 0.2, 0.3, -1., -0.9, -1.0, 100.0, 2
# const lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
# const tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits"
const lensname = "/mnt/simulations/MICE/voids_MICE.dat"
const tracname = "/home/fcaporaso/cats/MICE/micecat2_halos_full.fits"

const NTRACS = length(read(FITS(tracname)[2], "unique_halo_id"))
const LBOX = 3072 #Mpc/h box of mice
const MEAN_NTRAC = NTRACS/LBOX^3
# const MEANDENSITY = mean_den(0.0, 70, 0.3, 0.7)

radial_profile(RMIN, RMAX, DR, Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag, lensname, tracname)

