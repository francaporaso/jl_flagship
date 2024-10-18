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
    L = readdlm(lensname, Float32)

    m_rv   = @. (L[:,2] >= Rv_min) && (L[:,2] <= Rv_max)
    m_z    = @. (L[:,5] >= z_min) && (L[:,5] <= z_max)
    m_rho  = @. (L[:,9] >= rho1_min) && (L[:,9] <= rho1_max) && (L[:,10] >= rho2_min) && (L[:,10] <= rho2_max)
    m_flag = @. (L[:,12] >= flag) 

    mask = @views @. m_rv && m_z && m_rho && m_flag

    return L[mask,:]
end

"""
Loads the tracers catalog
"""
function traccat_load(z_min, z_max; 
                      lmhalo_min = nothing,              
                      tracname="/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits")

    f = FITS(tracname)[2]
    S = Matrix{Float32}([read(f, "xhalo") read(f, "yhalo") read(f, "zhalo") read(f, "lmhalo")]) #read(f, "flag_central")])

    if isnothing(lmhalo_min)
        return S
    end

    ## halos más grandes q 10 DM particles
    mp = 2.93e10 #Msun/h
    m_logm = @. S[:,end] > log10(10*mp)
    
    return S[m_logm, :]

    ### Filtros en flag_central y redshift:
    # S = Matrix{Float32}([read(f, "z_cgal") read(f, "xhalo") read(f, "yhalo") read(f, "zhalo") read(f, "lmhalo")]) #read(f, "flag_central")])
    # m_z = @. (S[:,1] >= (z_min-0.2)) && (S[:,1] <= (z_max+0.2))
    # m_flag = @. (S[:,end] == zero(Float32)) ## halos centrales
    # mask = @views @. m_z && m_flag
    # return S[mask,2:5]    

    # return S[m_z,2:end]
end

"""
Calcula el perfil de 1 void dados los halos S
"""
function individual_profile(S::Matrix{Float32}, 
                         RMIN, RMAX, NBINS,
                         rv, z, xv, yv, zv;
                         w=false, id=nothing)
    
    ### tcat[:,1] = logm
    ### tcat[:,2] = comovil_dist from center (xv,yv,zv) in units of void radius [rv]
    tcat = get_halos(S, RMIN, RMAX, rv, xv, yv, zv)
    MeanDen, MeanNTrac = mean_density_comovilshell(S, RMAX, rv, xv, yv, zv)

    NHalos   = zeros(NBINS)
    mass     = zeros(NBINS)
    Delta    = zeros(NBINS)
    DeltaCum = zeros(NBINS)
 
    ### calculamos el bin al que corresponde cada particula y sumando en el array 
    ### que corresponde (masa o halo)
    DR = (RMAX-RMIN)/NBINS

    for t in 1:size(tcat)[1]
        if (tcat[t,2] >= RMIN) && (tcat[t,2] <= RMAX)   
            ibin = ceil(Int32, (tcat[t,2] - RMIN)/DR)
            NHalos[ibin] += 1.0
            mass[ibin]  += 10.0 ^ tcat[t,1]
        end
    end
    
    mass_cum = cumsum(mass)
    NHalosCum = cumsum(NHalos)
    # MeanDen = mean_density_ball(tcat[:,1], rv, RMAX)
    # return mass, mass_cum, NHalos, NHalosCum, MeanDen

    for k in 0:NBINS-1
        Ri = (k*DR + RMIN)*rv
        # Rm = ((k+0.5)*DR + RMIN)*rv
        Rs = ((k+1.0)*DR + RMIN)*rv

        vol = (4pi/3) * (Rs^3 - Ri^3)
        Delta[k+1] = mass[k+1]/vol/MeanDen - 1.0
        NHalos[k+1] = NHalos[k+1]/vol/MeanNTrac - 1.0
        
        vol = (4pi/3) * (Rs^3)
        DeltaCum[k+1] = mass_cum[k+1]/vol/MeanDen - 1.0
        NHalosCum[k+1] = NHalosCum[k+1]/vol/MeanNTrac - 1.0
    end

    if w
        println("Saving in void_$id.dat")
        open("voids/void_$id.dat", "w") do io
            writedlm(io, [Delta DeltaCum NHalos NHalosCum])
        end
        return nothing
    end

    return Delta, DeltaCum, NHalos, NHalosCum, MeanDen
end

""" 
Perfil parcial, masa en el void y masa en una bola de 2RMAX de 1 void.
Para ser usada con stacking únicamente
"""
function partial_profile(RMIN, RMAX, NBINS,
                         rv, xv, yv, zv)
    
    ### tcat[:,1] = logm
    ### tcat[:,2] = comovil_dist from center (xv,yv,zv) in units of void radius [rv]
    tcat = get_halos(S, RMIN, RMAX, rv, xv, yv, zv)
    MassBall, HalosBall = mass_ball(S, 2RMAX, rv, xv, yv, zv)
    # MassShell, VolShell = mass_comovingshell(S, RMAX, rv, xv, yv, zv)

    NHalos = zeros(NBINS)
    mass   = zeros(NBINS)
 
    ### calculamos el bin al que corresponde cada particula y sumando en el array 
    ### que corresponde (masa o halo)
    DR = (RMAX-RMIN)/NBINS

    for t in 1:size(tcat)[1]
        if (tcat[t,2] >= RMIN) && (tcat[t,2] <= RMAX)   
            ibin = ceil(Int32, (tcat[t,2] - RMIN)/DR)
            NHalos[ibin] += 1.0
            mass[ibin]   += 10.0 ^ tcat[t,1]
        end
    end
    
    return mass, NHalos, MassBall, HalosBall
    # return mass, NHalos, MassShell, VolShell
end
