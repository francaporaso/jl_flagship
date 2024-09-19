using Printf
include("radial_profile.jl")

# ------------------------------------------------ #
#                                             main #
const RMIN, RMAX, NBINS = 0.05f0, 5.0f0, Int32(100)
const Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 10.0f0, 12.0f0, 0.2f0, 0.25f0, -1.0f0, -0.9f0, -1.0f0, 100.0f0, 2
# const lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
# const tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits"
const lensname = "/mnt/simulations/MICE/voids_MICE.dat"
const tracname = "/home/fcaporaso/cats/MICE/mice_halos_centralesF.fits"

# const NTRACS = length(read(FITS(tracname)[2], "unique_gal_id"))
# const NTRACS = length(read(FITS(tracname)[2], "unique_halo_id"))
# # const NTRACS = 4096^3 # según Fosalba et al. 2015
# const LBOX = 3072 #Mpc/h box of mice 
# const MEAN_NTRAC = NTRACS/LBOX^3
# const MEANDENSITY = mean_den(0.0, 70, 0.3, 0.7)

t = @elapsed begin
    # println("NTRACS.....: $NTRACS")
    # println("LBOX.......: $LBOX Mpc/h" )
    # println("MEAN_NTRAC.: $MEAN_NTRAC h/Mpc")

    # addprocs(4)

    radial_profile(RMIN, RMAX, NBINS, Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag, lensname, tracname)
end

t /= 60.0
@printf("Ended in %.2f minutes \n", t)
