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

"""
Total mass in a ball centered in rv and radius RMAX, in [Msun / h]
"""
function mass_ball(S, RMAX, rv, xv, yv, zv)
    halos = get_halos(S, 0.0, RMAX, rv, xv, yv, zv)
    mass = sum(10.0 .^ halos[:,1])
    return mass, size(halos)[1]
end

"""
Dado un sólo centro (xv,yv,zv) y su radio rv, encuentra los halos
al rededor del centro entre RMIN*rv hasta RMAX*rv
"""
function get_halos(S::Matrix{Float32},
                   RMIN, RMAX,
                   rv, xv, yv, zv)

    ### Máscara en una bola con centro (xv,yv,zv) y radio (1+2DR)RMAX*rv
    distance = @views @. sqrt((S[:,1] - xv)^2 + (S[:,2] - yv)^2 + (S[:,3] - zv)^2)/rv
    m1 = distance .< 1.1RMAX
    m2 = distance .> 0.9RMIN

    halos_list = hcat(S[m1 .&& m2, end], distance[m1 .&& m2])

    return halos_list
end