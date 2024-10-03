function array_split(ary, step)
    n = size(ary)[1]
    lbins = round(Int64, n/step)
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