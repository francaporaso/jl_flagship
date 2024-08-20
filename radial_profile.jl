using DelimitedFiles
using DataFrames
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
    L = readdlm(lensname)

    m_rv   = @views @. (L[:,2] >= Rv_min) & (L[:,2] <= Rv_max)
    m_z    = @views @. (L[:,5] >= z_min) & (L[:,5] <= z_max)
    m_rho  = @views @. (L[:,9] >= rho1_min) & (L[:,9] <= rho1_max) & (L[:,10] >= rho2_min) & (L[:,10] <= rho2_max)
    m_flag = @views @. (L[:,12] >= flag) 

    mask = @. m_rv & m_z & m_rho & m_flag

    return Matrix{Float64}(L[mask,:])
end

function get_tracers(RMAX, rv::Vector{Float64}, xv::Vector{Float64}, yv::Vector{Float64}, zv::Vector{Float64}; 
                    tracname="/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits")

    tcat = Matrix( DataFrame( FITS(tracname)[2] ) ) ## FITS lee .ftis; DataFrame transforma en tabla; Matrix transforma en Array{T,2}
    ## tcat[4] = xhalo
    ## tcat[5] = yhalo
    ## tcat[6] = zhalo
    ## tcat[3] = lmhalo

    nvoids = length(rv)
    trac_lists = Vector{Matrix{Float64}}(undef, nvoids)

    for v in 1:nvoids
        sq_distance = @views @. (tcat[:,4] - xv[v])^2 + (tcat[:,5] - yv[v])^2 + (tcat[:,6] - zv[v])^2
        mask = sq_distance .<= (RMAX*rv[v])^2

        trac_lists[v] = tcat[mask,:]
    end

    return trac_lists
end

function partial_profile(tcat::Matrix{Float64}, 
                        RMIN, RMAX, dr,
                        rv, xv, yv, zv)
    
    MeanDenTrac = 1.0
    NBINS = round(Int32,(RMAX-RMIN)/dr)
    
    RMIN *= rv
    RMAX *= rv
    dr   *= rv
    sq_distance = @. (tcat[:,4] - xv)^2 + (tcat[:,5] - yv)^2 + (tcat[:,6] - zv)^2

    Delta   = zeros(Float64, NBINS)
    DenTrac = zeros(Float64, NBINS)
    rad = RMIN

    for i in 1:NBINS
        tr = @views tcat[(sq_distance .>= rad^2) .& (sq_distance .< (rad+dr)^2), :]
        mass   = sum(10.0 .^ tr[:,3])
        volume = (4.0 *pi/ 3.0 )*((rad+dr)^3 - rad^3)
        
        Delta[i]   = (mass/volume)/MeanDenTrac - 1.0
        DenTrac[i] = length(tr[:,1])/volume

        rad += dr
    end

    return Delta, DenTrac
end


function radial_profile(RMIN, RMAX, dr, 
                        Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max,
                        flag, lensname::String, tracname::String)
    ## reading cats
    println("Reading catalogs...")

    L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag=flag, lensname=lensname)
    tracers_lists = get_tracers(RMAX, L[:,2], L[:,6], L[:,7], L[:,8], tracname=tracname)
    
    Nvoids = length(L[:,1])
    NBINS = round(Int32, (RMAX- RMIN)/dr)
    println("Nvoids: $Nvoids")

    println("Calculating profile...")
    A = partial_profile.(tracers_lists, RMIN, RMAX, dr, L[:,2], L[:,6], L[:,7], L[:,8])
    
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

    println("Plotting...")
    r = range(RMIN, RMAX, NBINS)
    p1 = plot(legend=false)
    p2 = plot(legend=false)
    # for i in 1:Nvoids
    #     plot!(p1, r, Delta[i,:], lc=:blue, alpha=0.5)
    #     plot!(p2, r, DenTrac[i,:], lc=:red, alpha=0.5) 
    # end

    plot!(p1, r, Delta_mean, yerr=Delta_std, lc=:darkblue)
    plot!(p2, r, Dentrac_mean, yerr=Dentrac_std, lc=:darkred)

    png(p1, "prof1_pru")
    png(p2, "prof2_pru")
    println("End!")
end

## -- main --
function main()
    RMIN, RMAX, dr = 0.01, 5., 0.1
    Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 5., 15., 0.2, 0.3, -1., -0.9, -1., 100., 2
    lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
    tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits"
    # lensname = "/mnt/simulations/MICE/voids_MICE.dat"
    # tracname = "/home/fcaporaso/cats/MICE/micecat2_halos_full.fits"
    radial_profile(RMIN, RMAX, dr, Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag, lensname, tracname)
end

main()
