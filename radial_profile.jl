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
    G = 4.30091727e-9 # (km/s)² Mpc / M_sun
    return 3.0*H0^2*Om0*(1+z)^3/(8.0*pi*G)
end

"""
Mean density in a ball centered in rv and radius RMAX*rv, mean den in [Msun h^2 / Mpc^3]
"""
function mean_density_ball(logm, rv, RMAX)
    mass = sum(10.0 .^ logm)
    vol = (4pi/3) * (RMAX*rv)^3
    return mass/vol
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
    L = DataFrame(readdlm(lensname), :auto)

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
    S = Matrix{Float32}([read(f, "z_cgal") read(f, "xhalo") read(f, "yhalo") read(f, "zhalo") read(f, "lmhalo") read(f, "flag_central")])

    m_z    = @. (S[:,1] >= z_min) && (S[:,1] <= z_max)
    m_flag = @. (S[:,end] == zero(Float32)) ## halos centrales

    mask = @views @. m_z && m_flag

    return S[mask,2:5]
end


"""
Dado un sólo centro (xv,yv,zv) y su radio rv, encuentra los halos
al rededor del centro entre RMIN*rv hasta RMAX*rv
"""
function get_halos(S::Matrix{Float32},
                   RMIN::Float64, RMAX::Float64, NBINS::Int64,
                   rv::Float64, xv::Float64, yv::Float64, zv::Float64)


    halos_list = Matrix{Float64}(undef,0,2)

    ### Máscara en una bola con centro (xv,yv,zv) y radio (1+2DR)RMAX*rv
    distance = @views @. sqrt((S[:,1] - xv)^2 + (S[:,2] - yv)^2 + (S[:,3] - zv)^2)/rv
    m1 = distance .< RMAX
    m2 = distance .> RMIN

    halos_list = vcat(halos_list, hcat(S[m1 .&& m2, end], distance[m1 .&& m2]))

    return halos_list
end

# """
# Dado un cto de centros ([xv], [yv], [zv]) con sus radios [rv]
# encuentra los trazadores alrededor de cada centro hasta RMAX*[rv]
# """
# function get_halos(RMAX::Float64, NBINS::Int64,
#                     rv::Vector{Float64}, xv::Vector{Float64}, yv::Vector{Float64}, zv::Vector{Float64}; 
#                     tracname="/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits",
#                     sorted=false)

#     tcat = Matrix(DataFrame(FITS(tracname)[2])) ## FITS lee .fits; DataFrame transforma en tabla
#     ## tcat[1] = id
#     ## tcat[2] = flagcentral
#     ## tcat[3] = lmhalo
#     ## tcat[4] = xhalo
#     ## tcat[5] = yhalo
#     ## tcat[6] = zhalo

#     nvoids = length(rv)
#     trac_list = Vector{Matrix{Float64}}(undef,nvoids)

#     ### Máscara en una bola con centro (xv,yv,zv) y radio (1+NBINS)RMAX*rv
#     @threads for v in 1:nvoids
#         distance = @views @. sqrt((tcat[:,4] - xv[v])^2 + (tcat[:,5] - yv[v])^2 + (tcat[:,6] - zv[v])^2)
#         mask = distance .<= (RMAX*rv[v])

#         trac_list[v] = hcat(tcat[mask,3], distance[mask]/rv[v])
#     end

#     if sorted
#         return @views trac_list[sortperm(trac_list[:,end]), :]
#     end

#     return trac_list
# end

"""
Calcula el perfil de 1 void dados los trazadores tcat
"""
function partial_profile(S::Matrix{Float32}, 
                         RMIN::Float64, RMAX::Float64, NBINS::Int64,
                         rv::Float64, z::Float64, xv::Float64, yv::Float64, zv::Float64)
    
    ### tcat[:,1] = logm
    ### tcat[:,2] = comovil_dist from center (xv,yv,zv) in units of void radius [rv]
    tcat = get_halos(S, RMIN, RMAX, NBINS, rv, xv, yv, zv)

    NHalos   = zeros(NBINS)
    mass     = zeros(NBINS)
    Delta    = zeros(NBINS)
    DeltaCum = zeros(NBINS)
 
    ### calculamos el bin al que corresponde cada particula y sumando en el array 
    ### que corresponde (masa o halo)
    DR = (RMAX-RMIN)/NBINS

    for t in 1:size(tcat)[1]
        if (tcat[t,2] >= RMIN) && (tcat[t,2] <= RMAX)
            ibin = ceil(Int64, (tcat[t,2] - RMIN)/DR)
            NHalos[ibin] += 1.0
            mass[ibin]  += 10.0 ^ tcat[t,1]
        end
    end
    
    mass_cum = cumsum(mass)
    NHalosCum = cumsum(NHalos)
    MeanDen = mean_density_ball(tcat[:,1], rv, RMAX)

    for k in 0:NBINS-1
        Ri = (k*DR + RMIN)*rv
        # Rm = ((k+0.5)*DR + RMIN)*rv
        Rs = ((k+1.0)*DR + RMIN)*rv

        vol = 4pi/3 * (Rs^3 - Ri^3)
        Delta[k+1] = mass[k+1]/vol/MeanDen #- 1.0
        NHalos[k+1] = NHalos[k+1]/vol/MEAN_NTRAC #- 1.0

        vol = 4pi/3 * (Rs^3)
        DeltaCum[k+1] = mass_cum[k+1]/vol/MeanDen #- 1.0
        NHalosCum[k+1] = NHalosCum[k+1]/vol/MEAN_NTRAC #- 1.0
    end

    # rad = RMIN
    # for i in 1:NBINS
    #     tr = @views tcat[(tcat[:,end] .< (rad+DR)) .&& (tcat[:,end] .>= rad), :]
    #     mass[i] = sum(10.0 .^ tr[:,1])
    #     NTrac[i] = length(tr[:,1])

    #     rad += DR
    # end

    # vol = 4pi/3*[(ri+DR)^3 - ri^3 for ri in range(RMIN,RMAX,NBINS)]
    # volcum = 4pi/3*[(ri+DR)^3 - RMIN^3 for ri in range(RMIN,RMAX,NBINS)]
    # Delta = (mass./vol)/MeanDen
    # DeltaCum = (cumsum(mass)./volcum)/MeanDen
    # NTracCum = cumsum(NTrac)./volcum/MEAN_NTRAC
    # NTrac ./= vol/MEAN_NTRAC

    return Delta, DeltaCum, NHalos, NHalosCum
end

"""
Calcula todos los perfiles de las lentes seleccionadas
"""
function radial_profile(RMIN, RMAX, NBINS, 
                        Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max,
                        flag, lensname::String, tracname::String)
    ## reading cats
    println("......................")
    println("Reading lenses...")
    L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag=flag, lensname=lensname)

    nvoids = nrow(L)
    println("Nvoids.....: $nvoids")
    println("Done!")

    println("......................")
    println("Reading halos...")
    S = traccat_load(z_min, z_max, tracname=tracname)
    println("Done!")
    
    println("......................")
    println("Calculating profile...")

    println("RMIN.......: $RMIN")
    println("RMAX.......: $RMAX")
    println("NBINS......: $NBINS")

    println("......................")
    println("Rvmin......: $Rv_min Mpc")
    println("Rvmax......: $Rv_max Mpc")
    println("zmin.......: $z_min")
    println("zmax.......: $z_max")
    println("rho1min....: $rho1_min")
    println("rho1max....: $rho1_max")
    println("rho2min....: $rho2_min")
    println("rho2max....: $rho2_max")

    Delta   = Matrix{Float64}(undef, nvoids, NBINS)
    DeltaCum = Matrix{Float64}(undef, nvoids, NBINS)
    Nhalos = Matrix{Float64}(undef, nvoids, NBINS)

    for i in 1:nvoids
        Delta[i,:], DeltaCum[i,:], Nhalos[i,:] = partial_profile(S,RMIN, RMAX, NBINS, L[i,2], L[i,5], L[i,6], L[i,7], L[i,8])
    end
    
    Delta_stack = sum(Delta, dims=1)'/nvoids
    Delta_std  = std(Delta, dims=1)'/nvoids
    
    DeltaCum_stack = sum(DeltaCum, dims=1)'/nvoids
    DeltaCum_std  = std(DeltaCum, dims=1)'/nvoids

    NTrac_stack = sum(Nhalos, dims=1)'/nvoids
    NTrac_std = std(Nhalos, dims=1)'/nvoids
    println("Done!")

    println("......................")
    println("Saving...")

    open("pru_stack.csv", "w") do io 
        writedlm(io, [Delta_stack Delta_std DeltaCum_stack DeltaCum_std NTrac_stack NTrac_std], ',')
    end

    open("pru_individual.csv", "w") do io 
        writedlm(io, [Delta DeltaCum Nhalos], ',')
    end

    println("Done!")
    println("......................")
end

"""
Calcula todos los perfiles de las lentes seleccionadas
"""
function test_profile(RMIN, RMAX, NBINS, 
                        Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max,
                        flag, lensname::String, tracname::String)
    ## reading cats
    println("......................")
    println("Reading lenses...")
    L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag=flag, lensname=lensname)

    nvoids = nrow(L)
    println("Nvoids.....: $nvoids")
    println("Done!")

    println("......................")
    println("Reading halos...")
    S = traccat_load(z_min, z_max, tracname=tracname)
    println("Done!")
    
    println("......................")
    println("Calculating profile...")

    println("RMIN.......: $RMIN")
    println("RMAX.......: $RMAX")
    println("NBINS......: $NBINS")

    println("......................")
    println("Rvmin......: $Rv_min Mpc")
    println("Rvmax......: $Rv_max Mpc")
    println("zmin.......: $z_min")
    println("zmax.......: $z_max")
    println("rho1min....: $rho1_min")
    println("rho1max....: $rho1_max")
    println("rho2min....: $rho2_min")
    println("rho2max....: $rho2_max")

    Delta   = zeros(NBINS)
    DeltaCum = zeros(NBINS)
    NHalos = zeros(NBINS)
    NHalosCum = zeros(NBINS)

    Delta, DeltaCum, NHalos, NHalosCum .= partial_profile(S, RMIN, RMAX, NBINS, L[1,2], L[1,5], L[1,6], L[1,7], L[1,8])
    
    println("Done!")

    println("......................")
    println("Saving...")

    open("test_profile.csv", "w") do io 
        writedlm(io, [Delta DeltaCum NHalos NHalosCum], ',')
    end

    println("Done!")
    println("......................")
end