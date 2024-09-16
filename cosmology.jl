"""
Hubble parameter squared (E function)
"""
H2(z, H0, Om0, Ode0) = H0^2 * (Om0*(1+z)^3 + Ode0)

"""
critical density at redshift z
"""
function critical_density(z, H0, Om0, Ode0)
    G = 4.30091727e-9 # (km/s)² Mpc / M_sun
    return 3.0*H2(z,H0,Om0,Ode0)/(8.0*pi*G)
end


"""
``\\rho_m (z) = \\rho_m (0) a^(-3) ``
            `` = \\rho_m(0) (1+z)^3 ``
            `` = 3 H_0 ^2 \\Omega_{m,0} (1+z)^3 / 8 π G ``
"""
function mean_density(z, H0, Om0, Ode0)
    return critical_density(0.0,H0,Om0,Ode0)*Om0*(1+z)^3
end
