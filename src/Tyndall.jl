module Tyndall
using Jultran
using Sundials
using JSON
export planck, dry_adiabat, isothermal, flux, compute_κ_grid, dry_adiabat_with_strat


# constants
p_s = 1.013e5 # p_s == surface pressure; ps == pressures
R = 287
c_p = 1006
σ = 5.67e-8
h = 6.626075540e-34    #Planck's constant
c = 2.99792458e8       #Speed of light
k =1.38065812e-23      #Boltzman thermodynamic constant
g = 9.8
cos_θ = .5

# settings

ν_min = 0
ν_max = 2500
no_of_ν_bands = 2 # 80
no_of_gs = 2 # 15
gs = 1:no_of_gs
ν_step = (ν_max - ν_min) / no_of_ν_bands
νs = ν_min:ν_step:ν_max
ν_bands = [ν:(ν+ν_step) for ν in νs[1:end-1]]
lnκs = -8:6

no_of_ps = 2 # 10
p_step = 1.013e5/(no_of_ps - 1)
ps = 0:p_step:1.013e5

T_min = 150
T_max = 340
no_of_Ts = 2 # 10
T_step = (T_max-T_min)/(no_of_Ts-1)
Ts = 150:T_step:340

# κ_grid = JSON.parse(open("$(Pkg.dir())/Tyndall/data/κ_grid.json"))

# ν in cm-1
planck(ν,T) = 100 * 2h * c^2 * (100ν)^3 / expm1(c*h*100ν/(k*T))

I_plus(ν_band::Range, T, q) = function(p, I_plus, I_plus_dot)
    ν = mean(ν_band)
    I_plus_dot[1] = q * κ(ν, g, p, T) * (I_plus - π * planck(ν,T(p))) / (g * cos_θ)
end

function κ(ν, g, p, T)
    p_weight = (p % p_step) / p_step
    p_i = fld(p,p_step) 
    T_weight = (T % T_step) / T_step
    T_i = fld(T,T_step) 
    
    sum(T_weight * (κ_grid[p_i:p_i+1,T_i:T,ν,g] .* [p_weight,1-p_weight]) + (1 - T_weight) * (κ_grid[p_i:p_i+1,T_i+1,ν,g] .* [p_weight,1-p_weight]))
end

function OLR(T,q; by_ν=false)
    if by_ν
        map(ν_bands) do ν_band
            sum(map(1:no_of_gs) do g
                Sundials.ode(I_plus(ν_band,g,T,q), [T(p_s)], [p_s, 0])
            end)
        end
    else
        sum(OLR(T,q,by_ν=true))
    end
end

# temperature profiles
isothermal(T_s) = p -> T_s
# dry adiabat
dry_adiabat(T_s) = p -> T_s * (p/p_s)^(R/c_p)
dry_adiabat_with_strat(T_s) = p -> p > 1.e4 ? T_s * (p/p_s)^(R/c_p) : T_s * (1.e4/p_s)^(R/c_p)

function compute_κ_grid()
    κ_grid = Array(Float64, (no_of_ps, no_of_Ts, no_of_ν_bands, no_of_gs))
    
    for (p_i,p)=enumerate(ps), (T_i,T)=enumerate(Ts)
        for (ν_i,ν_band)=enumerate(ν_bands)
            νs_in_band = ν_band[1]:.001:ν_band[end]
            println("p is $p")
            println("T is $T")
            println("ν_band is $(ν_band[1])")
            println("ν_band is $(ν_band[end])")
            zeros(Float64, length(νs_in_band))
            # cs = compute_κs(:co2, νs_in_band, T, p)
            # sorted_κs = sort(cs) # compute_κs from Jultran
            # g_size = length(νs_in_band) / no_of_gs
            # for g=gs
            #     κ_grid[p_i,T_i,ν_i,g] = mean(sorted_κs[(1:g_size) + floor(g_size * (g-1))])
            # end
        end
    end
    
    # json(open(Pkg.dir() + "/Tyndall/data/κ_grid.json", "w+"), κ_grid)
end

end # module
