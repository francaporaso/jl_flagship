"""
Hubble parameter squared (E function)
"""
H2(z, H0, Om0, Ode0) = H0^2 * (Om0*(1+z)^3 + Ode0)

"""
critical density at redshift z
"""
function critical_den(z, H0, Om0, Ode0)
    G = 4.30091727e-9 # (km/s)² Mpc / M_sun
    return 3.0*H2(z,H0,Om0,Ode0)/(8.0*pi*G)
end

"""
ρₘ = 3 (H₀)² Ωₘ₀ / 8 π G * (1+z)^3
where:
- [H0] = km/s/Mpc
- [ρₘ] = M_sun/Mpc^3
"""
function mean_den(z, H0, Om0, Ode0)
    return critical_den(z,H0,Om0,Ode0)*Om0*(1+z)^3
end
