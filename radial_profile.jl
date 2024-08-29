using DelimitedFiles
using DataFrames ## USAR! con df[!, col] el ! usa @views, es más rapido que @views A[:,col] con A matrix
using FITSIO
using Plots
using Statistics

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

    m_rv   = @. (L[!,2] >= Rv_min) & (L[!,2] <= Rv_max)
    m_z    = @. (L[!,5] >= z_min) & (L[!,5] <= z_max)
    m_rho  = @. (L[!,9] >= rho1_min) & (L[!,9] <= rho1_max) & (L[!,10] >= rho2_min) & (L[!,10] <= rho2_max)
    m_flag = @. (L[!,12] >= flag) 

    mask = @. m_rv & m_z & m_rho & m_flag

    return L[mask,:]
end

"""
Dado un cto de centros (xv, yv, zv) con su correspondiente radio rv
encuentra los trazadores al rededor del centro hasta 1.01 RMAX*rv
"""
function get_tracers(RMAX::Float64, rv::Vector{Float64}, xv::Vector{Float64}, yv::Vector{Float64}, zv::Vector{Float64}; 
                    tracname="/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits")

    tcat = Matrix(DataFrame(FITS(tracname)[2])) ## FITS lee .ftis; DataFrame transforma en tabla
    ## tcat[1] = id
    ## tcat[2] = flagcentral
    ## tcat[3] = lmhalo
    ## tcat[4] = xhalo
    ## tcat[5] = yhalo
    ## tcat[6] = zhalo

    nvoids = length(rv)
    trac_lists = Vector{Matrix{Float64}}(undef, nvoids)

    for v in 1:nvoids
        distance = @views @. sqrt((tcat[:,4] - xv[v])^2 + (tcat[:,5] - yv[v])^2 + (tcat[:,6] - zv[v])^2)
        mask = distance .<= ((1.0+dr)*RMAX*rv[v])

        tr = hcat(tcat[mask,:], distance[mask])
        trac_lists[v] = tr[sortperm(tr[:,end]),:]
    end

    return trac_lists
end

function partial_profile(tcat::Matrix{Float64}, 
                        RMIN, RMAX, dr,
                        rv, xv, yv, zv)
    
    MeanDen   = 1.0
    MeanNTrac = 1.0
    NBINS = length(RMIN:dr:RMAX)
    
    RMIN *= rv
    RMAX *= rv
    dr   *= rv
    #sq_distance = @. (tcat[:,4] - xv)^2 + (tcat[:,5] - yv)^2 + (tcat[:,6] - zv)^2

    Delta = zeros(Float64, NBINS)
    DeltaCum = zeros(Float64, NBINS)
    NTrac = zeros(Float64, NBINS)

    for m in 1:size(tcat,1)
        if (tcat[m,end] >= RMIN) && (tcat[m,end] <= RMAX)
            ibin = ceil(Int, (tcat[m,end]-RMIN)/dr)
            Delta[ibin] += 10.0 ^ tcat[m, 3]
            NTrac[ibin] += 1.0
        end
    end

    for i in 1:NBINS
        ### cuidado con esto... lo hace Polaco pero no estoy seguro por qué
        if NTrac[i] < 3.0
            Delta[i] = 0.0
            NTrac[i] = 0.0
        end

        DeltaCum[i] += Delta[i]
    end

    for k in 1:NBINS
        Ri = k*dr + RMIN
        Rm = (k+0.5)*dr + RMIN
        Rs = (k+1.0)*dr + RMIN

        volumen = (4.0/3.0)*pi*(Rs^3 - Ri^3)
        # Delta[k] = Delta[k]/volumen/MeanDen - 1.0
        Delta[k] = Delta[k]/volumen

        volcum = (4.0/3.0)*pi*(Rs^3)
        # DeltaCum[k] = DeltaCum[k]/volumen/MeanDen - 1.0
        DeltaCum[k] = DeltaCum[k]/volcum
    end

    # for i in 1:NBINS
    #     #tr = @views tcat[(sq_distance .>= rad^2) .& (sq_distance .< (rad+dr)^2), :]
    #     tr = @views tcat[(tcat[:,end] .< (rad+dr)^2) .&& (tcat[:,end] .>= rad^2), :]
    #     mass   = sum(10.0 .^ tr[:,3])
    #     volume = (4.0 *pi/ 3.0 )*((rad+dr)^3 - rad^3)
        
    #     Delta[i]   = (mass/volume)/MeanDenTrac - 1.0
    #     DenTrac[i] = length(tr[:,1])#/volume

    #     rad += dr
    # end

    return Delta, DeltaCum, NTrac
end


function radial_profile(RMIN, RMAX, dr, 
                        Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max,
                        flag, lensname::String, tracname::String)
    ## reading cats
    println("Reading lenses...")
    L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag=flag, lensname=lensname)
    println("Reading tracers...")
    tracers_lists = get_tracers(RMAX, L[!,2], L[!,6], L[!,7], L[!,8], tracname=tracname)
    
    Nvoids = nrows(L[!,1])
    NBINS = length(RMIN:dr:RMAX)
    println("Nvoids: $Nvoids")

    println("Calculating profile...")
    A = partial_profile.(tracers_lists, RMIN, RMAX, dr, L[!,2], L[!,6], L[!,7], L[!,8])
    
    Delta   = Matrix{Float64}(undef, Nvoids, NBINS)
    DenTrac = Matrix{Float64}(undef, Nvoids, NBINS)

    for i in 1:Nvoids
        Delta[i,:] = A[i][1]
        DenTrac[i,:] = A[i][2]
    end

    Delta_mean = mean(Delta, dims=1)'
    Delta_std  = std(Delta, dims=1)'
    
    Dentrac_mean = mean(DenTrac, dims=1)'
    Dentrac_std  = std(DenTrac, dims=1)'

    println("Saving...")

    writedlm("pru.csv", Delta, ',')

    println("End!")
end

## -- main --
function main()
    RMIN, RMAX, dr = 0.0, 5., 0.05
    Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 10., 12., 0.2, 0.25, -1., -0.9, 0.2, 100, 2
    lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
    tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits"
    # lensname = "/mnt/simulations/MICE/voids_MICE.dat"
    # tracname = "/home/fcaporaso/cats/MICE/micecat2_halos_full.fits"


    radial_profile(RMIN, RMAX, dr, Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag, lensname, tracname)
end

#main() 
