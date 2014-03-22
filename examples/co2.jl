using Tyndall
using PyPlot

ν_min = 10
ν_max = 2500
no_of_ν_bands = 80
ν_step = (ν_max - ν_min) / no_of_ν_bands
νs = ν_min:ν_step:ν_max
ν_bands = [ν:(ν+ν_step) for ν in νs[1:end-1]]


# T = isothermal(300)
# output = OLR(p -> min(dry_adiabat(283.15)(p), 160),.0001 * (44/29),by_ν=true)
# noutput = [π * planck(mean(ν_band), 300) for ν_band in ν_bands]
# 
