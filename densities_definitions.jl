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
function mean_density_ball(S, RMAX, rv, xv, yv, zv)
    halos = get_halos(S, 0.0, RMAX, rv, xv, yv, zv)
    mass = sum(10.0 .^ halos[:,1])    
    vol = (4pi/3) * (RMAX*rv)^3
    return mass/vol
end

"""
Mean density in a comoving shell centered in
the observer between √(xᵥ+yᵥ+zᵥ)-RMAX*rv and √(xᵥ+yᵥ+zᵥ)-RMAX*rv , in [Msun h^2 / Mpc^3]
"""
function mean_density_comovingshell(S, RMAX,
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