# using Base.Threads
using Cosmology
using Unitful
include("radial_profile.jl")

function mean_density_box(logm, rv, RMAX)
    mass = sum(10.0 .^ logm)
    vol = (4pi/3) * (RMAX*rv)^3
    return mass/vol
end

function get_tracers_z(z_min, z_max, xv, yv, zv; 
    tracname="/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits")

    tcat = Matrix(DataFrame(FITS(tracname)[2])) ## FITS lee .fits; DataFrame transforma en tabla
    ## tcat[1] = id
    ## tcat[2] = flagcentral
    ## tcat[3] = lmhalo
    ## tcat[4] = xhalo
    ## tcat[5] = yhalo
    ## tcat[6] = zhalo

    trac_list = Array{Float64}(undef, 0)

    ### Máscara en una cáscara esférica con centro (xv,yv,zv) y radios comoving_radial_dist(z_min),comoving_radial_dist(z_min) 
    cosmo = cosmology(h=1, OmegaM=0.25, Tcmb=0.0)
    distance = @views @. sqrt((tcat[:,4] - xv)^2 + (tcat[:,5] - yv)^2 + (tcat[:,6] - zv)^2)
    m1 = distance .<= ustrip(comoving_radial_dist(cosmo, z_max))
    m2 = distance .>= ustrip(comoving_radial_dist(cosmo, z_min))

    trac_list = vcat(tcat[m1 .&& m2, 3])

    return trac_list
end

function mean_density_comovilshell(z_min, z_max;
                                   tracname="/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits")
    
    cosmo = cosmology(h=1, OmegaM=0.25, Tcmb=0.0)
    vol = (1/8)*(4pi/3)*(ustrip(comoving_radial_dist(cosmo, z_max))^3 - ustrip(comoving_radial_dist(cosmo, z_min))^3)
    halos = get_tracers_z(z_min, z_max, 0.0, 0.0, 0.0, tracname=tracname)
    mass = sum(10.0 .^ halos[:])

    return mass/vol
end

function test_box()
    println("testing box")

    RMIN, RMAX, NBINS = 0.05, 5., 100
    Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 10., 12., 0.2, 0.25, -1., -0.9, -1.0, 100.0, 2
    # lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
    # tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits"
    lensname = "/mnt/simulations/MICE/voids_MICE.dat"
    tracname = "/home/fcaporaso/cats/MICE/micecat2_halos_full.fits"

    # leyendo cat de lentes
    println("leyendo voids")
    L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, lensname=lensname)
    nvoids = nrow(L)
    println("$nvoids")
    # leyendo halos al rededro de xv, yv, zv hasta RMAX
    println("leyendo halos")
    tr = get_tracers(RMAX, NBINS, L[!,2], L[!,6], L[!,7], L[!,8], tracname=tracname)

    # densidad media de cada void en bola al rededor del centro hasta RMAX
    println("calculando densidades en la bola")
    ρ = zeros(nvoids)
    @threads for i in 1:nvoids
        ρ[i] = mean_density_box(tr[i][:,1], L[i,2], RMAX)
    end

    # densidad media teorica
    println("calculado den media teorica")
    ρ_universe = mean_density.(L[!,5], H0, Om0, Ode0)

    # if ρ ./ ρ_universe .<= 10.0
    #     println("DAN IGUAL!, TEST APROBADO")
    # end

    println("guardando")
    open("den_box_test.csv", "w") do io
        writedlm(io,  [ρ ρ_universe], ',')
    end

    println("TUKI!")
end

function test_comoving_shell()
    println("testing comoving shell")

    RMIN, RMAX, NBINS = 0.05, 5., 100
    Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, flag = 10., 12., 0.2, 0.25, -1., -0.9, -1.0, 100.0, 2
    # lensname = "/home/franco/FAMAF/Lensing/cats/MICE/voids_MICE.dat"
    # tracname = "/home/franco/FAMAF/Lensing/cats/MICE/mice-halos-cut.fits"
    lensname = "/mnt/simulations/MICE/voids_MICE.dat"
    tracname = "/home/fcaporaso/cats/MICE/micecat2_halos_full.fits"

    # leyendo cat de lentes
    println("leyendo voids")
    L = lenscat_load(Rv_min, Rv_max, z_min, z_max, rho1_min, rho1_max, rho2_min, rho2_max, lensname=lensname)
    nvoids = nrow(L)
    println("$nvoids")

    # densidad media en cascaron esferico comovil de la simu entre z_min, z_max
    println("calculando densidades en la bola")
    ρ = mean_density_comovilshell(z_min, z_max, tracname=tracname)

    # densidad media teorica
    println("calculado den media teorica")
    ρ_universe = mean_density.(L[!,5])

    # if ρ ./ ρ_universe .<= 10.0
    #     println("DAN IGUAL!, TEST APROBADO")
    # end

    println("guardando")
    open("den_comovingshell_test.csv", "w") do io
        writedlm(io,  [fill(ρ,nvoids) ρ_universe], ',')
    end

    println("TUKI!")
end

function plot_result(voidid, ρ, ρ_universe)
    p = plot(legend=:outertopright)
    # plot!(p, L[!, 1], md, label="Mean density box")
    # plot!(p, L[!, 1], md_universe, label="Mean density universe")
    plot!(p, voidid, ρ ./ ρ_universe, label="\\rho_{box} / \\rho_{uni} ")
    plot!(yscale=:log10)
    display(p)
end

function plot_result!(p, voidid, ρ, ρ_universe)
    # p = plot(legend=:outertopright)
    # plot!(p, L[!, 1], md, label="Mean density box")
    # plot!(p, L[!, 1], md_universe, label="Mean density universe")
    plot!(p, voidid, ρ ./ ρ_universe, label="\\rho_{box} / \\rho_{uni} ")
    plot!(yscale=:log10)
    display(p)
end


### ---------------------------------------------- main
t = @elapsed begin
    # test_box_shell()
    test_comoving_shell()
end
println("terminado en $t sec")

# pp = true

# if pp
#     using Plots
#     dat = readdlm("den_box_test.csv", ',')
#     plot_result(range(length=size(dat)[1]), dat[:,1], dat[:,2])
# end
