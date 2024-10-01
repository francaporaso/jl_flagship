using Printf
include("radial_profile.jl")

# ------------------------------------------------ #
#                                             main #
const RMIN, RMAX, NBINS = 0.05f0, 5.0f0, Int32(100)
const Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 10.0, 12.0, 0.2, 0.25, -1.0f0, -0.9f0, -1.0, 100.0f0, 2
# const lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
# const tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice_halos_cut.fits"
const lensname = "/mnt/simulations/MICE/voids_MICE.dat"
const tracname = "/home/fcaporaso/cats/MICE/mice_halos_centralesF.fits"

t = @elapsed begin
    radial_profile(RMIN, RMAX, NBINS, Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag, lensname, tracname)
    # calculate_all_indiv(RMIN, RMAX, NBINS, Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag, lensname, tracname)
end

t /= 60.0
@printf("Ended in %.2f minutes \n", t)
