module Tyndall
using Jlitran
using DataFrames
export planck, dry_adiabat, isothermal, flux, compute_Hs, dry_adiabat_with_strat

# constants
p_s = 1.013e5 # p_s == surface pressure; ps == pressures
R = 287
c_p = 1006
σ = 5.67e-8
h = 6.626e-34
c = 3.e8
k = 1.38e-23
g = 9.8

# settings
no_of_bands = 200
ν_max = 3000
bands = [ν_max*(i-1/no_of_bands):ν_max*i for i = 1/no_of_bands:1/no_of_bands:1]
lnκs = -8:6

# cached values
Hs = readcsv("$(Pkg.dir())/Tyndall/data/Hs.csv")

# ν in cm-1
planck(ν,T) = 100 * 2h * c^2 * (100ν)^3 / expm1(c*h*100ν/(k*T))

# temperature profiles
isothermal(T_s) = p -> T_s
# dry adiabat
# dry_adiabat(T_s) = p -> T_s * (p/p_s)^(R/c_p)
dry_adiabat_with_strat(T_s) = p -> p > 1.e4 ? T_s * (p/p_s)^(R/c_p) : T_s * (1.e4/p_s)^(R/c_p)

flux(T, qs) = ps -> sum([flux(band_i, T, qs)(ps) for band_i = 2:length(bands)])
flux(band_i, Τ, qs) = function(ps)
    band = bands[band_i]
    ν = mean(band)

    layers = length(ps) - 1
    fluxes = Array(Float64,(layers,2))
    
    # calculate transmissions
    ts = Array(Float64, (layers+1,layers+1))
    t = transmission(band_i, qs)
    for i = 1:layers+1, j = 1:i
        ts[i,j] = t(ps[i],ps[j])
    end
    
    Bs = [π * planck(ν,.5*(Τ(ps[i])+Τ(ps[i+1]))) for i = 1:layers]
    I_s = π * planck(ν,Τ(ps[end]))
    for i = 1:layers
        fluxes[i,1] = sum([Bs[i] * (ts[i+1,j+1] - ts[i+1,j]) for j=1:i]) # down
        fluxes[i,2] = I_s * t(ps[i], ps[layers+1]) + sum([Bs[i] * (ts[j,i] - ts[j+1,i]) for j=i:layers]) # up
    end
    
    (band[end]-band[1]) * fluxes
end

transmission(band_i, qs) = function(p_1, p_2)
    if p_1 == p_2
        1.
    else
        Hs = [H(band_i, :co2, lnκ_i) for lnκ_i=1:length(lnκs)]
        sum([exp(-(exp(lnκ)+exp(lnκ+1) * qs[:co2] * abs(p_1^2 - p_2^2)) / (4p_s * g)) * (Hs[lnκ+10] - Hs[lnκ+9]) for lnκ=-8:5])
    end
end

H(band_i, gas, lnκ_i) = Hs[band_i, lnκ_i]

function compute_Hs()
    νs = 25:.001:3000
    co2 = DataFrame(
        νs = νs,
        κs = compute_κs(:co2, νs, 260, 1.e4)
    )

    Hs = Array(Float64,(length(bands), length(lnκs)))
    
    for i = 1:length(bands)
        band = bands[i]
        slice = co2[band[1] .<= co2["νs"] .<= band[end],:]
        for j = 1:length(lnκs)
            lnκ = lnκs[j]
            # println(length(slice[slice["κs"] .< exp(lnκ),:]["κs"])/length(slice[:,:]["κs"]))
            Hs[i,j] = (band[end] - band[1]) * length(slice[slice["κs"] .< exp(lnκ),:]["κs"])/length(slice[:,:]["κs"])
        end
    end
    
    writecsv("data/Hs.csv", Hs)
end

end # module
