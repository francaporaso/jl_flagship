using FITSIO, DelimitedFiles
using Printf

RMIN, RMAX, NBINS = 0.0f0, 5.0f0, Int32(50)
Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 10.0f0, 12.0f0, 0.2f0, 0.3f0, -1.0f0, -0.8f0, -1.0f0, 100.0f0, 2.0f0

lensname = "/mnt/simulations/MICE/voids_MICE.dat"
tracname = "/home/fcaporaso/cats/MICE/mice_halos_centralesF.fits"

mp = 2.93f10  # Msun/h

f = FITS(tracname)[2]
xhalo = read(f, "xhalo")
yhalo = read(f, "yhalo")
zhalo = read(f, "zhalo")
lmhalo = read(f, "lmhalo")

mask_particles = (lmhalo .> log10(10 * mp))
xhalo = xhalo[mask_particles]
yhalo = yhalo[mask_particles]
zhalo = zhalo[mask_particles]
lmhalo = lmhalo[mask_particles]

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
    
    mask = @. (L[:,2] >= Rv_min) && (L[:,2] < Rv_max) && (L[:,5] >= z_min) && (L[:,5] < z_max) && (L[:,9] >= rho1_min) && (L[:,9] < rho1_max) && 
    (L[:,10] >= rho2_min) && (L[:,10] < rho2_max) && (L[:,12] >= flag) 
    
    return L[mask,:]
end

function get_halos(RMIN, RMAX, 
                    rv, xv, yv, zv;
                    tracname="/home/fcaporaso/cats/MICE/mice_halos_centralesF.fits")

    distance = @views @. sqrt((xhalo - xv)^2 + (yhalo - yv)^2 + (zhalo - zv)^2) / rv
    
    mask_ball = (distance .< 5RMAX) .&& (distance .>= 0.0)
    mask_prof = (distance .< RMAX) .&& (distance .>= RMIN)
    
    massball = sum(10.0f0 .^ lmhalo[mask_ball])
    halosball = length(distance[mask_ball])

    return lmhalo[mask_prof], distance[mask_prof], massball, halosball
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
    
    logm, distance, MassBall, HalosBall = get_halos(RMIN, RMAX, rv, xv, yv, zv)

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

N= 10
L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max)
for i in 1:N
    res = partial_profile(RMIN,RMAX,NBINS,L[i,2], L[i,6], L[i,7], L[i,8])
    name = @sprintf "void_jl_%d.csv" L[i,1]
    writedlm(name, [res[1] res[2] fill(res[3]) fill(res[4])], ',')
end