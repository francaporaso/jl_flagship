using Distributed
using Printf

NCORES = 16
addprocs(NCORES)

RMIN, RMAX, NBINS = 0.0, 5.0, 50
Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 6.0, 9.622, 0.2, 0.4, -1.0, -0.8, -1.0, 100.0, 2
filename = @sprintf "radialprof_stack_R_%.0f_%.0f_z%.1f_%.1f.csv" Rv_min Rv_max z_min z_max
# lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
# tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice_halos_cut.fits"
# lensname = "/mnt/simulations/MICE/voids_MICE.dat"
# tracname = "/home/fcaporaso/cats/MICE/mice_halos_centralesF.fits"

@everywhere begin 
    using FITSIO, DelimitedFiles
end

@everywhere begin
    """
    Divides a array in n/step parts. 
    If n%step = 0 then each part has the same number of elements,
    if not the last element in the resulting array has the rest
    """
    function array_split(ary, step)
        n = size(ary)[1]
        lbins = round(Int64, n/step)
        if n % step != 0
            lbins += 1
        end

        sub_ary = Vector{Matrix{Float64}}(undef, lbins)
        x = 1
        for i in 1:lbins
            try
                sub_ary[i] = ary[x:x+step-1,:]
                x += step
            catch e
                sub_ary[i] = ary[x:end,:]
            end
        end

        return sub_ary
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
    Dado un sólo centro (xv,yv,zv) y su radio rv, encuentra los halos
    al rededor del centro entre RMIN*rv hasta RMAX*rv
    """
    function get_halos(S::Matrix{Float32},
                    RMIN, RMAX,
                    rv, xv, yv, zv)

        ### Máscara en una bola con centro (xv,yv,zv) y radio (1+2DR)RMAX*rv
        distance = @views @. sqrt((S[:,1] - xv)^2 + (S[:,2] - yv)^2 + (S[:,3] - zv)^2)/rv
        m1 = distance .<= RMAX
        m2 = distance .>= RMIN

        halos_list = hcat(S[m1 .&& m2, end], distance[m1 .&& m2])

        return halos_list
    end

    """
    Loads the lenses catalog
    """
    function lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max; 
                        flag=2.0, lensname="/mnt/simulations/MICE/voids_MICE.dat")
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

        m_rv   = @views @. (L[:,2] >= Rv_min) && (L[:,2] <= Rv_max)
        m_z    = @views @. (L[:,5] >= z_min) && (L[:,5] <= z_max)
        m_rho  = @views @. (L[:,9] >= rho1_min) && (L[:,9] <= rho1_max) && (L[:,10] >= rho2_min) && (L[:,10] <= rho2_max)
        m_flag = @views @. (L[:,12] >= flag) 

        mask = @views @. m_rv && m_z && m_rho && m_flag

        return L[mask,:]
    end

    """
    Loads the tracers catalog
    """
    function traccat_load(; lmhalo_min = nothing,              
                        tracname="/home/fcaporaso/cats/MICE/mice_halos_centralesF.fits")

        f = FITS(tracname)[2]
        S = Matrix{Float32}([read(f, "xhalo") read(f, "yhalo") read(f, "zhalo") read(f, "lmhalo")]) #read(f, "flag_central")])

        if !isnothing(lmhalo_min)
            return S
        end

        ## halos más grandes q 10 DM particles
        mp = 2.93e10 #Msun/h
        m_logm = @views @. S[:,end] > log10(10*mp)
        
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
    Perfil parcial, masa en el void y masa en una bola de 2RMAX de 1 void.
    Para ser usada con stacking únicamente
    """
    function partial_profile(RMIN, RMAX, NBINS,
                            rv, xv, yv, zv;
                            tracname="/home/fcaporaso/cats/MICE/mice_halos_centralesF.fits")
                            
        #"/home/franco/FAMAF/Lensing/cats/MICE/mice_halos_cut.fits"

        ### tcat[:,1] = logm
        ### tcat[:,2] = comovil_dist from center (xv,yv,zv) in units of void radius [rv]
        S = traccat_load(lmhalo_min=1, tracname=tracname)
        tcat = get_halos(S, RMIN, RMAX, rv, xv, yv, zv)
        MassBall, HalosBall = mass_ball(S, 5RMAX, rv, xv, yv, zv)
        # MassShell, VolShell = mass_comovingshell(S, RMAX, rv, xv, yv, zv)

        NHalos = zeros(NBINS)
        mass   = zeros(NBINS)
    
        ### calculamos el bin al que corresponde cada particula y sumando en el array 
        ### que corresponde (masa o halo)
        DR = (RMAX-RMIN)/NBINS

        for t in 1:size(tcat)[1]
            ibin = ceil(Int32, (tcat[t,2] - RMIN)/DR)
            NHalos[ibin] += 1.0
            mass[ibin]   += 10.0 ^ tcat[t,1]
        end
        
        return mass, NHalos, MassBall, HalosBall

    end

end ## fin de @everywhere begin

"""
Paralelizado de partial_profile
"""
function paralellization(partial, NCORES,
                         RMIN, RMAX, NBINS, 
                         Rv_min, Rv_max, z_min, z_max, rho2_min, rho2_max; 
                         lensname="/mnt/simulations/MICE/voids_MICE.dat")

    L = lenscat_load(Rv_min, Rv_max, z_min, z_max, -1.0, -0.8, rho2_min, rho2_max, lensname=lensname)
    nvoids = size(L)[1]
    println("NVOIDS: .... $nvoids")

    resmap = pmap(partial, fill(RMIN,nvoids), fill(RMAX,nvoids), fill(NBINS,nvoids), L[:,2], L[:,6], L[:,7], L[:,8], batch_size=NCORES)

    return resmap
end

"""
Calcula el stacking de un cto de perfiles de masa resmap
"""
function stacking(resmap,
                  RMIN, RMAX, NBINS;
                  filename = "pru_stack_par.csv")

    mass  = zeros(NBINS)
    halos = zeros(NBINS)
    massball  = 0.0
    halosball = 0.0

    for i in 1:nvoids
        mass  += resmap[i][1]
        halos += resmap[i][2]
        massball  += resmap[i][3]
        halosball += resmap[i][4]
    end
        
    meandenball   = massball/(4pi/3 * (5RMAX)^3)
    meanhalosball = halosball/(4pi/3 * (5RMAX)^3)

    DR = (RMAX-RMIN)/NBINS
    
    vol = zeros(NBINS)
    volcum = zeros(NBINS)
    for k in 1:NBINS
        vol[k] = (4*pi/3)*(((k+1.0)*DR + RMIN)^3 - (k*DR + RMIN)^3)    
        volcum[k] = (4*pi/3)*((k+1.0)*DR + RMIN)^3
    end
    
    Delta    = (mass./vol)/meandenball .- 1
    DeltaCum = (cumsum(mass)./volcum)/meandenball .- 1
    DeltaHalos    = (halos./vol)/meanhalosball .- 1
    DeltaHalosCum = (cumsum(halos)./volcum)/meanhalosball .- 1

    println("Saving in: ", filename)    

    open(filename, "w") do io 
        writedlm(io, [Delta DeltaCum DeltaHalos DeltaHalosCum], ',')
    end
end

stacking(
        paralellization(
            partial_profile, NCORES,
            RMIN, RMAX, NBINS,
            Rv_min, Rv_max, z_min, z_max, rho2_min, rho2_max
        ), 
        RMIN, RMAX, NBINS,
        filename=filename
)