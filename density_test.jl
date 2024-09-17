include("radial_profile.jl")

const RMIN, RMAX, NBINS = 0.5, 4., 100
const Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 10., 12., 0.2, 0.25, -1., -0.9, -1.0, 100.0, 2
const lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
const tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits"
# const lensname = "/mnt/simulations/MICE/voids_MICE.dat"
# const tracname = "/home/fcaporaso/cats/MICE/micecat2_halos_full.fits"

const NTRACS = length(read(FITS(tracname)[2], "unique_gal_id"))
# const NTRACS = length(read(FITS(tracname)[2], "unique_halo_id"))
# const NTRACS = 4096^3 # según Fosalba et al. 2015
const LBOX = 3072 #Mpc/h box of mice 
const MEAN_NTRAC = NTRACS/LBOX^3

function mean_density_box(logm, rv, RMAX, DR)
    mass = sum(10.0 .^ logm)
    vol = 4pi/3 * ((1+2DR)*RMAX*rv)^3
    return mass/vol
end

### ----------------- main
L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max)
nvoids = nrow(L)
tr = get_tracers(RMAX, NBINS, L[!,2], L[!,6], L[!,7], L[!,8])

@assert typeof(tr) == Vector{Matrix{Float64}}

md_list = zeros(nvoids)
for i in 1:nvoids
    md[i] = mean_density_box(tr[i][:,1], L[i,2], RMAX, dr)
end

open("den_box_test.csv", "w") do io
    writedlm(io,  md, ',')
end

println("END!")