using DelimitedFiles
using DataFrames
using FITSIO
using Plots

function load_lenscat(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag,
                      filename)
    
    L = DataFrame(readdlm(filename), ["id", "rv", "ra", "dec", "z", "xv", "yv", "zv", "rho1", "rho2", "p", "flag"])
    m_rv  = (L.rv .>= Rv_min) .& (L.rv .< Rv_max)
    m_z   = (L.z .>= z_min) .& (L.z .< z_max)
    m_rho = (L.rho1 .>= rho1_min) .& (L.rho1 .< rho1_max) .& (L.rho1 .>= rho2_min) .& (L.rho1 .< rho2_max)
    m_fl  = (L.flag .>= 2)
    mask = m_rv .& m_z .& m_rho .& m_fl

    return L[mask,:]
end

function get_tracers(tcat, RMAX,
                     rv, xv, yv, zv)
    
    sq_distance = (tcat.xhalo .- xv).^2 .+ (tcat.yhalo .- yv).^2 .+ (tcat.zhalo .- zv).^2
    tcat[!, :sq_d] = sq_distance
    sort!(tcat, :sq_d)
    # tcat[tcat.sq_d .<= (RMAX*rv)^2, :] ### no haría falta, porque el loop para
    # el perfil simplemente itera hasta que llega al maximo permitido...
    
    return tcat
end

function partial_profile(tcat, 
                         RMIN, RMAX, dr,
                         rv, xv, yv, zv)
    
    MeanDenTrac = 1.
    NBINS = round(Int32,(RMAX-RMIN)/dr)
    
    RMIN *= rv
    RMAX *= rv
    dr *= rv
    
    Delta   = zeros(Float32, NBINS)
    DenTrac = zeros(Float32, NBINS)
    rad = RMIN

    for i in 1:NBINS
        tr = tcat[(tcat.sq_d .>= rad) .&& (tcat.sq_d .< rad+dr), :]
        mass = sum(10.0 .^ tr.lmhalo)
        volume = (4. *pi/ 3. )*((rad+dr)^3 - rad^3)
        Delta[i] = (mass/volume)/MeanDenTrac - 1.0
        DenTrac[i] = nrow(tr)/volume
        rad += dr
    end

    return Delta, DenTrac
end


function radial_profile(RMIN, RMAX, dr, 
                        Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag,
                        lensname, tracname)
    ## reading cats
    println("Reading catalogs...")

    L = load_lenscat(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag, lensname)
    f = FITS(tracname)
    tracers = DataFrame(f[2])
    
    Nvoids = nrow(L)
    NBINS = round(Int32, (RMAX- RMIN)/dr)
    println("Nvoids: $Nvoids")

    println("Calculating profile...")
    #tcat = get_tracers.(tracers, RMAX, L.rv, L.xv, L.yv, L.zv)
    #Delta, DenTrac = partial_profile.(tcat, RMIN, RMAX, dr, L.rv, L.xv, L.yv, L.zv)
    
    Delta   = zeros(Float32, Nvoids, NBINS)
    DenTrac = zeros(Float32, Nvoids, NBINS)
    ## for void in voids: get_tracers, calculate profile, add to vector.
    for i in 1:Nvoids
        tcat = get_tracers(tracers, RMAX, L.rv[i], L.xv[i], L.yv[i], L.zv[i])
        Delta[i,:], DenTrac[i,:] = partial_profile(tcat, RMIN, RMAX, dr, L.rv[i], L.xv[i], L.yv[i], L.zv[i])
    end

    println("Plotting")
    r = range(RMIN, RMAX, NBINS)
    p = plot(legend=false)
    for i in 1:Nvoids
        plot!(p, r, Delta[i,:], lc=:blue)
        plot!(p, r, DenTrac[i,:], lc=:red)
    end

    png(p, "prof_pru")
    println("End!")
end

## -- main --
RMIN, RMAX, dr = 0.01, 5., 0.1
Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 10., 12., 0.2, 0.25, -1., 1., -1., 100., 2
#lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
#tracname = "../FlagShip/tests/testcat_uniform.fits"
lensname = "/mnt/simulations/MICE/voids_MICE.dat"
tracname = "/home/fcaporaso/cats/MICE/micecat2_halos_full.fits"
radial_profile(RMIN, RMAX, dr, Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag, filename)

