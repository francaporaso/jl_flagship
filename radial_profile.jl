using DelimitedFiles
using DataFrames
using FITSIO
using Statistics
using Plots


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
encuentra los trazadores alrededor del centro hasta (1+dr)RMAX*rv
"""
function get_tracers(RMAX::Float64, dr::Float64,
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

    ### Máscara en una bola con centro (xv,yv,zv) y radio (1+dr)RMAX*rv
    for v in 1:nvoids
        distance = @views @. sqrt((tcat[:,4] - xv[v])^2 + (tcat[:,5] - yv[v])^2 + (tcat[:,6] - zv[v])^2)
        mask = distance .<= ((1.0+dr)*RMAX*rv[v])

        tr = hcat(tcat[mask,:], distance[mask])
        trac_lists[v] = tr[sortperm(tr[:,end]),:] ## ordena de menor a mayor distance
    end

    return trac_lists
end

"""
Calcula el perfil de 1 void dados los trazadores tcat al rededor de (xv,yv,zv)
"""
function partial_profile(tcat::Matrix{Float64}, 
                        RMIN, RMAX, dr,
                        rv, xv, yv, zv)
    
    MeanDen   = 1.0
    MeanNTrac = NTRACS/LBOX
    NBINS = length(RMIN:dr:RMAX)
    
    RMIN *= rv
    RMAX *= rv
    dr   *= rv

    Delta = zeros(NBINS)
    # DeltaCum = zeros(NBINS)
    NTrac = zeros(NBINS)

    for m in 1:size(tcat,1)
        if (tcat[m,end] >= RMIN) && (tcat[m,end] <= RMAX)
            ibin = ceil(Int, (tcat[m,end]-RMIN)/dr)
            Delta[ibin] += 10.0 ^ tcat[m, 3]
            NTrac[ibin] += 1.0
        end
    end

    # for i in 1:NBINS
    #     ### cuidado con esto... lo hace Polaco pero no estoy seguro por qué
    #     if NTrac[i] < 3.0
    #         Delta[i] = 0.0
    #         NTrac[i] = 0.0
    #     end
    # end

    DeltaCum = cumsum(Delta)

    for k in 0:NBINS-1
        Ri = k*dr + RMIN
        Rm = (k+0.5)*dr + RMIN
        Rs = (k+1.0)*dr + RMIN

        volumen = (4.0/3.0)*pi*(Rs^3 - Ri^3)
        Delta[k+1] = Delta[k+1]/volumen/MeanDen
        #Delta[k] /= volumen
        NTrac[k+1] = NTrac[k+1]/volumen/MeanNTrac 

        volcum = (4.0/3.0)*pi*(Rs^3)
        DeltaCum[k+1] = DeltaCum[k+1]/volumen/MeanDen
        #DeltaCum[k] /= volcum
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
    tracers_lists = get_tracers(RMAX, dr, L[!,2], L[!,6], L[!,7], L[!,8], tracname=tracname)
    
    nvoids = nrow(L)
    NBINS = length(RMIN:dr:RMAX)
    println("nvoids: $nvoids")

    println("Calculating profile...")
    
    Delta   = Matrix{Float64}(undef, nvoids, NBINS)
    DeltaCum = Matrix{Float64}(undef, nvoids, NBINS)
    NTracers = Matrix{Float64}(undef, nvoids, NBINS)

    for i in 1:nvoids
        Delta[i,:], DeltaCum[i,:], NTracers[i,:] = partial_profile(tracers_lists[i], RMIN, RMAX, dr, L[i,2], L[i,6], L[i,7], L[i,8])
    end

    # display(Delta)
    # display(DeltaCum)

    begin
        fig = plot(legend=false, layout=3)
        for i in 1:nvoids
            plot!(fig[1], RMIN:dr:RMAX, Delta[i,:], c=:blue)
            plot!(fig[2], RMIN:dr:RMAX, DeltaCum[i,:], c=:red)
            plot!(fig[3], RMIN:dr:RMAX, NTracers[i,:], c=:green)   
        end
        # display(fig)
        savefig(fig, "pru.png")
    end
    
    # Delta_mean = mean(Delta, dims=1)'
    # Delta_std  = std(Delta, dims=1)'
    
    # Dentrac_mean = mean(DenTrac, dims=1)'
    # Dentrac_std  = std(DenTrac, dims=1)'

    # println("Saving...")

    # writedlm("pru.csv", Delta, ',')

    println("End!")
end


# ------------------------------------------------ #
#                                             main #
RMIN, RMAX, dr = 0.01, 5., 0.05
const Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 10., 12., 0.2, 0.25, -1., -0.9, 0.2, 100.0, 2
const lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
const tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits"
#const lensname = "/mnt/simulations/MICE/voids_MICE.dat"
#const tracname = "/home/fcaporaso/cats/MICE/micecat2_halos_full.fits"

const NTRACS = length(read(FITS(tracname)[2], "unique_gal_id"))
const LBOX = 3072 #Mpc/h

## ̄ρₘ = 3 (H₀)² Ωₘ₀ / 8 π G

#radial_profile(RMIN, RMAX, dr, Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag, lensname, tracname)

