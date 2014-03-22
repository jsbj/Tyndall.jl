module Tyndall
using Jultran
using Sundials
using JSON
export planck, dry_adiabat, isothermal, OLR, compute_κ_grid, dry_adiabat_with_strat


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

ν_min = 10
ν_max = 2500
no_of_ν_bands = 80
ν_step = (ν_max - ν_min) / no_of_ν_bands
νs = ν_min:ν_step:ν_max
ν_bands = [ν:(ν+ν_step) for ν in νs[1:end-1]]
no_of_γs = 15
γs = 1:no_of_γs

no_of_ps = 4
p_step = 1.013e5/(no_of_ps - 1)
ps = 0:p_step:1.013e5

T_min = 150
T_max = 340
no_of_Ts = 4
T_step = (T_max-T_min)/(no_of_Ts-1)
Ts = 150:T_step:340

κ_grid = reshape(JSON.parse(open("$(Pkg.dir())/Tyndall/data/κ_grid.json")), (no_of_ps, no_of_Ts, no_of_ν_bands, no_of_γs))

# ν in cm-1
planck(ν,T) = 100 * 2h * c^2 * (100ν)^3 / expm1(c*h*100ν/(k*T))

I_plus(ν_i, γ, T, q) = function(p, I_plus, I_plus_dot)
    ν = mean(ν_bands[ν_i])
    I_plus_dot[1] = q * κ(ν_i, γ, p, T(p)) * (I_plus[1] - π * planck(ν,T(p))) / (g * cos_θ)
end

function κ(ν_i, γ, p, T)
    if p >= 0
        p_weight = (p % p_step) / p_step
        p_i = fld(p,p_step) + 1
        if p_weight < 1.e-8
            p_weight = 1
            p_i -= 1
        end

        T_weight = (T % T_step) / T_step
        T_i = fld(T-T_min,T_step) + 1
        sum(T_weight * (κ_grid[p_i:p_i+1,T_i,ν_i,γ] .* [p_weight,1-p_weight]) + (1 - T_weight) * (κ_grid[p_i:p_i+1,T_i+1,ν_i,γ] .* [p_weight,1-p_weight]))
    else
        0
    end
end

function OLR(T,q; by_ν=false)
    if by_ν
        map(1:no_of_ν_bands) do ν_i
            sum(map(1:no_of_γs) do γ
                Sundials.ode(I_plus(ν_i,γ,T,q), [π * planck(mean(ν_bands[ν_i]),T(p_s))], [p_s, 0])[2]
            end)/no_of_γs
        end
    else
        sum(OLR(T,q,by_ν=true)) * ν_step
    end
end

# temperature profiles
isothermal(T_s) = p -> T_s
# dry adiabat
dry_adiabat(T_s) = p -> T_s * (p/p_s)^(R/c_p)
# dry_adiabat_with_strat(T_s) = p -> p > 1.e4 ? T_s * (p/p_s)^(R/c_p) : T_s * (1.e4/p_s)^(R/c_p)
dry_adiabat_with_strat(T_s) = p -> p > p_s * (160.0/T_s) ^(c_p/R) ? T_s * (p/p_s)^(R/c_p) : 160.0

function compute_κ_grid()
    κ_grid = Array(Float64, (no_of_ps, no_of_Ts, no_of_ν_bands, no_of_γs))
    
    for (p_i,p)=enumerate(ps), (T_i,T)=enumerate(Ts)
        println("p is $p")
        println("T is $T")
        ν_i = 0
        for ν_band=ν_bands
            println("ν_i is $(ν_i)")
            ν_i += 1
            νs_in_band = ν_band[1]:.01:ν_band[end]
            sorted_κs = sort(compute_κs(:co2, νs_in_band, T, p)) # compute_κs from Jultran
            γ_size = length(νs_in_band) / no_of_γs
            for γ=γs
                κ_grid[p_i,T_i,ν_i,γ] = mean(filter(x->!isnan(x),sorted_κs[(1:γ_size) + floor(γ_size * (γ-1))]))
            end
        end
    end
    
    open("$(Pkg.dir())/Tyndall/data/κ_grid.json", "w+") do f
        JSON.print(f, reshape(κ_grid,(no_of_ps * no_of_Ts * no_of_ν_bands * no_of_γs)))
    end
    
    κ_grid
end

end # module
