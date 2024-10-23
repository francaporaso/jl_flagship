using Distributed
using Printf
using ProgressMeter

# NCORES = Int32(10)
# addprocs(NCORES)
NCORES = nprocs()

RMIN, RMAX, NBINS = 0.0f0, 5.0f0, Int32(50)
Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 6.0f0, 9.622f0, 0.2f0, 0.4f0, -1.0f0, -0.8f0, -1.0f0, 0.0f0, 2.0f0
filename = @sprintf "radialprof_stack_R_%.0f_%.0f_z%.1f_%.1f.csv" Rv_min Rv_max z_min z_max
# lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
# tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice_halos_cut.fits"
lensname = "/mnt/simulations/MICE/voids_MICE.dat"
tracname = "/home/fcaporaso/cats/MICE/mice_halos_centralesF.fits"

@everywhere begin 
    using FITSIO, DelimitedFiles
end

@everywhere begin
    
    function get_halos(RMIN, RMAX, 
                       rv, xv, yv, zv;
                       tracname="/home/fcaporaso/cats/MICE/mice_halos_centralesF.fits")

        mp = 2.93f10  # Msun/h
        
        f = FITS(tracname)[2]
        xhalo = read(f, "xhalo")
        yhalo = read(f, "yhalo")
        zhalo = read(f, "zhalo")
        lmhalo = read(f, "lmhalo")
    
        distance = @views @. sqrt((xhalo - xv)^2 + (yhalo - yv)^2 + (zhalo - zv)^2) / rv
    
        mask = (distance .<= RMAX) .&& (distance .>= RMIN) .&& (lmhalo .> log10(10 * mp))
        
        lmhalo = lmhalo[mask]
        distance = distance[mask]

        massball = sum(10.0f0 .^ lmhalo)
        halosball = length(distance)
    
        return lmhalo, distance, massball, halosball
    end

    """ 
    Perfil parcial, masa en el void y masa en una bola de 2RMAX de 1 void.
    Para ser usada con stacking únicamente
    """
    function partial_profile(RMIN, RMAX, NBINS,
                            rv, xv, yv, zv;
                            tracname="/home/fcaporaso/cats/MICE/mice_halos_centralesF.fits")
                            #"/home/franco/FAMAF/Lensing/cats/MICE/mice_halos_cut.fits")

        ### tcat[:,1] = logm
        ### tcat[:,2] = comovil_dist from center (xv,yv,zv) in units of void radius [rv]
        
        logm, distance, MassBall, HalosBall = get_halos(0.0f0, 5RMAX, rv, xv, yv, zv)

        mask = @. (distance <= RMAX) && (distance >= RMIN)
        logm, distance = logm[mask], distance[mask]

        @assert length(logm) == length(distance)

        NHalos = zeros(NBINS)
        mass   = zeros(NBINS)
    
        ### calculamos el bin al que corresponde cada particula y sumando en el array 
        ### que corresponde (masa o halo)
        DR = (RMAX-RMIN)/NBINS

        for (lm, d) in zip(logm,distance)
            ibin = ceil(Int32, (d - RMIN)/DR)
            NHalos[ibin] += 1.0
            mass[ibin]   += 10.0 ^ lm
        end
        
        return mass, NHalos, MassBall, HalosBall
    end

end ## fin de @everywhere begin

"""
Loads the lenses catalog
"""
function lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max; 
                    flag=2.0f0, lensname="/mnt/simulations/MICE/voids_MICE.dat")
    ## L[1] = id
    ## L[2] = rv
    ## L[5] = z
    ## L[6] = xv
    ## L[7] = yv
    ## L[8] = zv
    ## L[9] = rho1
    ## L[10] = rho2
    ## L[12] = flag
    L::Matrix{Float32} = Matrix(readdlm(lensname, Float32))

    mask = @. (L[:,2] >= Rv_min) && (L[:,2] <= Rv_max) && (L[:,5] >= z_min) && (L[:,5] <= z_max) && (L[:,9] >= rho1_min) && (L[:,9] <= rho1_max) && 
              (L[:,10] >= rho2_min) && (L[:,10] <= rho2_max) && (L[:,12] >= flag) 

    return L[mask,:]
end

"""
Calcula el stacking.
"""
function stacking(NCORES,
                  RMIN, RMAX, NBINS,
                  Rv_min, Rv_max, z_min, z_max, rho2_min, rho2_max;
                  lensname="/mnt/simulations/MICE/voids_MICE.dat",
                  filename = "pru_stack_par.csv")

    L = lenscat_load(Rv_min, Rv_max, z_min, z_max, -1.0, -0.8, rho2_min, rho2_max, lensname=lensname)
    nvoids = size(L)[1]

    println("NVOIDS: .... $nvoids")
    println("Coriendo en paralelo")
    println("NCORES: .... $NCORES")

    ### TODO
    ### probablemente el bug de memoria desaparece si se implementa el Lsplit...

    ### ----------- intento de solución
    mass  = zeros(NBINS)
    halos = zeros(NBINS)
    massball  = 0.0
    halosball = 0.0

    i = 1
    @showprogress for j in NCORES:NCORES:nvoids
        resmap = pmap(
            partial_profile,
            fill(RMIN,NCORES), fill(RMAX,NCORES), fill(NBINS,NCORES), L[i:j,2], L[i:j,6], L[i:j,7], L[i:j,8],
            batch_size=NCORES,
            )
        
        for res in resmap
            mass  += res[1]
            halos += res[2]
            massball  += res[3]
            halosball += res[4]
        end

        i+=NCORES
    end
    ### -------------------------------

    #resmap = pmap(partial_profile, fill(RMIN,nvoids), fill(RMAX,nvoids), fill(NBINS,nvoids), L[:,2], L[:,6], L[:,7], L[:,8], batch_size=NCORES)

    println("Done!")
    println("Stacking...")


        
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

    println("DONE!")

    return 0
end

stacking(
        NCORES,
        RMIN, RMAX, NBINS,
        Rv_min, Rv_max, z_min, z_max, rho2_min, rho2_max,
        lensname=lensname,
        filename=filename,
)
