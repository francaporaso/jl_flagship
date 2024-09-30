using DelimitedFiles
using DataFrames
using FITSIO
using Statistics
using Base.Threads
# using Distributed

"""
Mean density of the universe in [Msun h^2 / Mpc^3], calculated like:
``\\rho_m (z) = \\rho_m (0) a^(-3) ``
            `` = \\rho_m(0) (1+z)^3 ``
            `` = 3 H_0 ^2 \\Omega_{m,0} (1+z)^3 / 8 π G ``
"""
function mean_density_universe(z; H0=100.0, Om0=0.25, Ode0=0.75)
    G = 4.3009173e-9 # (km/s)² Mpc / M_sun
    return 3.0*H0^2*Om0*(1+z)^3/(8.0*pi*G)
end

"""
Mean density in a ball centered in rv and radius RMAX*rv, in [Msun h^2 / Mpc^3]
"""
function mean_density_ball(logm, rv, RMAX)
    mass = sum(10.0 .^ logm)
    vol = (4pi/3) * (RMAX*rv)^3
    return mass/vol
end

"""
Total mass in a ball centered in rv and radius RMAX, in [Msun / h]
"""
function mass_ball(S, RMAX, rv, xv, yv, zv)
    halos = get_halos(S, 0.0, RMAX, rv, xv, yv, zv)
    mass = sum(10.0 .^ halos[:,1])
    return mass, size(halos)[1]
end

"""
Mean density in a comoving shell centered in
the observer between √(xᵥ+yᵥ+zᵥ)-RMAX*rv and √(xᵥ+yᵥ+zᵥ)-RMAX*rv , in [Msun h^2 / Mpc^3]
"""
function mean_density_comovilshell(S, RMAX,
                                   rv, xv, yv, zv)

    # cosmo = cosmology(h=1, OmegaM=0.25, Tcmb=0.0)
    χ_min = sqrt(xv^2 + yv^2 + zv^2) - RMAX*rv
    χ_max = sqrt(xv^2 + yv^2 + zv^2) + RMAX*rv

    m1 = @. sqrt(S[:,2]^2 + S[:,3]^2 + S[:,4]^2) > χ_min
    m2 = @. sqrt(S[:,2]^2 + S[:,3]^2 + S[:,4]^2) < χ_max
    logm = S[m1 .&& m2, end]

    vol = (1/8)*(4pi/3)*(χ_max^3 - χ_min^3)
    mass = sum(10.0 .^ logm)

    return mass/vol, length(logm)/vol
end

"""
Total mass in a comoving shell and its volume centered in
the observer between √(xᵥ+yᵥ+zᵥ)-RMAX*rv and √(xᵥ+yᵥ+zᵥ)-RMAX*rv.
Mass in [Msun / h] and volume in [Rv^3]
"""
function mass_comovingshell(S, RMAX,
                           rv, xv, yv, zv)

    # cosmo = cosmology(h=1, OmegaM=0.25, Tcmb=0.0)
    χ_min = sqrt(xv^2 + yv^2 + zv^2) - RMAX*rv
    χ_max = sqrt(xv^2 + yv^2 + zv^2) + RMAX*rv

    m1 = @. sqrt(S[:,2]^2 + S[:,3]^2 + S[:,4]^2) > χ_min
    m2 = @. sqrt(S[:,2]^2 + S[:,3]^2 + S[:,4]^2) < χ_max
    logm = S[m1 .&& m2, end]

    mass = sum(10.0 .^ logm)
    vol = 1/8 * (4pi/3)*(χ_max^3 - χ_min^3)/rv^3
    return mass, vol
end

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
    L = DataFrame(readdlm(lensname, Float32), :auto)

    m_rv   = @. (L[!,2] >= Rv_min) && (L[!,2] <= Rv_max)
    m_z    = @. (L[!,5] >= z_min) && (L[!,5] <= z_max)
    m_rho  = @. (L[!,9] >= rho1_min) && (L[!,9] <= rho1_max) && (L[!,10] >= rho2_min) && (L[!,10] <= rho2_max)
    m_flag = @. (L[!,12] >= flag) 

    mask = @views @. m_rv && m_z && m_rho && m_flag

    return L[mask,:]
end

"""
Loads the tracers catalog
"""
function traccat_load(z_min, z_max; 
                      tracname="/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits")
    ## S.unique_halo_id
    ## S.z_cgal (true redshift)
    ## S.xhalo
    ## S.yhalo
    ## S.xhalo
    ## S.zhalo
    ## S.lmhalo
    ## S.flag_central
    ## S.cgal (comoving distance)
    f = FITS(tracname)[2]
    S = Matrix{Float32}([read(f, "z_cgal") read(f, "xhalo") read(f, "yhalo") read(f, "zhalo") read(f, "lmhalo")]) #read(f, "flag_central")])

    m_z = @. (S[:,1] >= z_min) && (S[:,1] <= z_max)
    
    ### para el catalogo q tenemos ya están filtrados, este paso es al pedo
    # m_flag = @. (S[:,end] == zero(Float32)) ## halos centrales
    # mask = @views @. m_z && m_flag
    # return S[mask,2:5]    

    return S[m_z,2:end]
end


"""
Dado un sólo centro (xv,yv,zv) y su radio rv, encuentra los halos
al rededor del centro entre RMIN*rv hasta RMAX*rv
"""
function get_halos(S::Matrix{Float32},
                   RMIN, RMAX,
                   rv, xv, yv, zv)


    halos_list = Matrix{Float64}(undef,0,2)

    ### Máscara en una bola con centro (xv,yv,zv) y radio (1+2DR)RMAX*rv
    distance = @views @. sqrt((S[:,1] - xv)^2 + (S[:,2] - yv)^2 + (S[:,3] - zv)^2)/rv
    m1 = distance .< RMAX
    m2 = distance .> RMIN

    halos_list = vcat(halos_list, hcat(S[m1 .&& m2, end], distance[m1 .&& m2]))

    return halos_list
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
function partial_profile(S::Matrix{Float32}, 
                         RMIN, RMAX, NBINS,
                         rv, z, xv, yv, zv)
    
    ### tcat[:,1] = logm
    ### tcat[:,2] = comovil_dist from center (xv,yv,zv) in units of void radius [rv]
    tcat = get_halos(S, RMIN, RMAX, rv, xv, yv, zv)
    MassBall, HalosBall = mass_ball(S, 2RMAX, rv, xv, yv, zv)

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
end

"""
Calcula el perfil stackeado de las lentes seleccionadas.
"""
function radial_profile(RMIN, RMAX, NBINS, 
                        Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max,
                        flag, lensname, tracname)
    ## reading cats
    println("......................")
    println("Reading lenses with:")

    println("Rvmin......: $Rv_min Mpc")
    println("Rvmax......: $Rv_max Mpc")
    println("zmin.......: $z_min")
    println("zmax.......: $z_max")
    println("rho1min....: $rho1_min")
    println("rho1max....: $rho1_max")
    println("rho2min....: $rho2_min")
    println("rho2max....: $rho2_max")

    L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag=flag, lensname=lensname)

    nvoids = nrow(L)
    println("Nvoids.....: $nvoids")
    println("Done!")

    println("......................")
    println("Reading halos...")
    S = traccat_load(z_min, z_max, tracname=tracname)
    println("Done!")
    
    println("......................")
    println("Calculating profile for:")

    println("RMIN.......: $RMIN")
    println("RMAX.......: $RMAX")
    println("NBINS......: $NBINS")

    DR = (RMAX-RMIN)/NBINS
    
    mass   = zeros(NBINS)
    NHalos = zeros(NBINS)
    MassBall  = 0.0
    HalosBall = 0.0

    for i in 1:nvoids
        res = partial_profile(S, RMIN, RMAX, NBINS, L[i,2], L[i,5], L[i,6], L[i,7], L[i,8])

        mass   += res[1]
        NHalos += res[2]
        MassBall  += res[3]
        HalosBall += res[4]
    end
    
    Vol = zeros(NBINS)
    VolCum = zeros(NBINS)
    for k in 0:NBINS-1
        Vol[k+1] = (4pi/3) * (((k+1.0)*DR + RMIN)^3 - (k*DR + RMIN)^3)
        VolCum[k+1] = (4pi/3) * ((k+1.0)*DR + RMIN)^3
    end

    ## TODO no funca...
    Delta = (mass./Vol) / (MassBall/(4pi/3 * (2RMAX)^3)) .- 1 
    DeltaCum = (cumsum(mass)./VolCum) / (MassBall/(4pi/3 * (2RMAX)^3)) .- 1 
    DenHalos = (NHalos./Vol) / (MassBall/(4pi/3 * (2RMAX)^3))
    DenHalosCum = (cumsum(mass)./VolCum) / (HalosBall/(4pi/3 * (2RMAX)^3))
    
    println("Done!")

    println("......................")
    println("Saving...")

    open("pru_stack.csv", "w") do io 
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

function calculate_all_indiv(RMIN, RMAX, NBINS, 
                             Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max,
                             flag, lensname, tracname)
    println("Calculando todos los perfiles individuales para:")

    println("Rvmin......: $Rv_min Mpc")
    println("Rvmax......: $Rv_max Mpc")
    println("zmin.......: $z_min")
    println("zmax.......: $z_max")
    println("rho1min....: $rho1_min")
    println("rho1max....: $rho1_max")
    println("rho2min....: $rho2_min")
    println("rho2max....: $rho2_max")

    L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag=flag, lensname=lensname)

    nvoids = nrow(L)
    println("Nvoids.....: $nvoids")
    println("Done!")

    println("......................")
    println("Reading halos...")
    S = traccat_load(z_min, z_max, tracname=tracname)
    println("Done!")

    println("......................")
    println("Calculating profile for:")

    println("RMIN.......: $RMIN")
    println("RMAX.......: $RMAX")
    println("NBINS......: $NBINS")

    @threads for i in 1:nvoids
        id = Int(L[i,1])
        if !isfile("void_$id.dat")
            individual_profile(S, RMIN, RMAX, NBINS, L[i,2], L[i,5], L[i,6], L[i,7], L[i,8], w=true, id=Int64(L[i,1]))
        else
            continue
        end
    end

end