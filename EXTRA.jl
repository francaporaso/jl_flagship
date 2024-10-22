"""
Calcula el perfil de 1 void dados los halos S
"""
function individual_profile(S::Matrix{Float32}, 
                         RMIN, RMAX, NBINS,
                         rv, z, xv, yv, zv;
                         w=false, id=nothing)
    
    ### tcat[:,1] = logm
    ### tcat[:,2] = comovil_dist from center (xv,yv,zv) in units of void radius [rv]
    tcat = get_halos(S, RMIN, RMAX, rv, xv, yv, zv)
    MeanDen, MeanNTrac = mean_density_comovilshell(S, RMAX, rv, xv, yv, zv)

    NHalos   = zeros(NBINS)
    mass     = zeros(NBINS)
    Delta    = zeros(NBINS)
    DeltaCum = zeros(NBINS)
 
    ### calculamos el bin al que corresponde cada particula y sumando en el array 
    ### que corresponde (masa o halo)
    DR = (RMAX-RMIN)/NBINS

    for t in 1:size(tcat)[1]
        if (tcat[t,2] >= RMIN) && (tcat[t,2] <= RMAX)   
            ibin = ceil(Int32, (tcat[t,2] - RMIN)/DR)
            NHalos[ibin] += 1.0
            mass[ibin]  += 10.0 ^ tcat[t,1]
        end
    end
    
    mass_cum = cumsum(mass)
    NHalosCum = cumsum(NHalos)
    # MeanDen = mean_density_ball(tcat[:,1], rv, RMAX)
    # return mass, mass_cum, NHalos, NHalosCum, MeanDen

    for k in 0:NBINS-1
        Ri = (k*DR + RMIN)*rv
        # Rm = ((k+0.5)*DR + RMIN)*rv
        Rs = ((k+1.0)*DR + RMIN)*rv

        vol = (4pi/3) * (Rs^3 - Ri^3)
        Delta[k+1] = mass[k+1]/vol/MeanDen - 1.0
        NHalos[k+1] = NHalos[k+1]/vol/MeanNTrac - 1.0
        
        vol = (4pi/3) * (Rs^3)
        DeltaCum[k+1] = mass_cum[k+1]/vol/MeanDen - 1.0
        NHalosCum[k+1] = NHalosCum[k+1]/vol/MeanNTrac - 1.0
    end

    if w
        println("Saving in void_$id.dat")
        open("voids/void_$id.dat", "w") do io
            writedlm(io, [Delta DeltaCum NHalos NHalosCum])
        end
        return nothing
    end

    return Delta, DeltaCum, NHalos, NHalosCum, MeanDen
end

"""
Divides a array in n/step parts. 
If n%step = 0 then each part has the same number of elements,
if not the last element in the resulting array has the rest
"""
function array_split(ary, step)
    n = size(ary)[1]
    lbins = round(Int64, n/step)
    if n % step != 0
        lbins += 1
    end

    sub_ary = Vector{Matrix{Float64}}(undef, lbins)
    x = 1
    for i in 1:lbins
        try
            sub_ary[i] = ary[x:x+step-1,:]
            x += step
        catch e
            sub_ary[i] = ary[x:end,:]
        end
    end

    return sub_ary
end 