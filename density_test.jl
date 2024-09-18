include("radial_profile.jl")

function mean_density_box(logm, rv, RMAX)
    mass = sum(10.0 .^ logm)
    vol = (4pi/3) * (RMAX*rv)^3
    return mass/vol
end

function mean_density_comovilshell(logm, rv, RMAX)
    mass = sum(10.0 .^ logm)
    vol = (4pi/3) * (RMAX*rv)^3
    return mass/vol
end


### ----------------- main
# DR = (RMAX)/NBINS

const pp = true

const RMIN, RMAX, NBINS = 0.05, 5., 100
const Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 10., 12., 0.2, 0.25, -1., -0.9, -1.0, 100.0, 2
const lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
const tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits"

# leyendo cat de lentes
L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max)
nvoids = nrow(L)
# leyendo halos al rededro de xv, yv, zv hasta RMAX
tr = get_tracers(RMAX, NBINS, L[!,2], L[!,6], L[!,7], L[!,8])

# calculando la densidad media de cada void en caja
md = zeros(nvoids)
for i in 1:nvoids
    md[i] = mean_density_box(tr[i][:,1], L[i,2], RMAX)
end

open("den_box_test.csv", "w") do io
    writedlm(io,  md, ',')
end


# densidad media teorica
md_universe = mean_density.(L[!,5], H0, Om0, Ode0)

println("END!")

using Plots
if pp
    p = plot(legend=:outertopright)
    # plot!(p, L[!, 1], md, label="Mean density box")
    # plot!(p, L[!, 1], md_universe, label="Mean density universe")
    plot!(p, L[!,1], md./md_universe, label="\\rho_{box} / \\rho_{uni} ")
    plot!(yscale=:log10)
    display(p)
end
