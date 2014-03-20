module TestTyndall
using Tyndall
using PyCall
using FactCheck

# facts("Planck function") do
#     @fact abs(planck(500,250)-0.09) < 0.01 => true
#     @fact abs(planck(500,350)-0.225) < 0.01 => true
# end
# 
# facts("dry_adiabat creates a dry adiabat temperature profile with a specified surface temperature") do
#     T = dry_adiabat(300)
#     p_s = 1.013e5
#     @fact T(p_s) => 300
#     @fact abs(T(p_s/2)-.82 * 300) < 1 => true
# end

facts("produces fluxes for a given atmospheric condition") do
    # context("isothermal") do
    #     ps = [1.013e5i for i=0:1/10:1]
    #     T = isothermal(300)
    #     fluxes = flux(T, {:co2 => 400})(ps)
    #     @fact abs(fluxes[1,2]-5.67e-8 * 300^4) < 4.0 => true
    # end
    
    context("dry adiabat") do
        ps = [1.013e5i for i=0:1/101:1][2:end]
        T = dry_adiabat_with_strat(300)
        co2 = 400
        
        # call the Python (which calls the Fortran...)
        @pyimport climt
        cam3 = climt.make_radiation(
            scheme="cam3", co2=co2,
            p=ps[1:end-1]/100, ps=ps[end]/100, T=[T(p) for p=ps[1:end-1]], 
            Ts=T(ps[end]), q=[1.e-10 for i=ps[1:end-1]]
        )
        
        # call the Julia
        fluxes = flux(T, [:co2 => co2])(ps)
        println(fluxes[1,2])
        println(getindex(cam3, "lwuflx")[1])
        @fact abs(fluxes[:up][1] - getindex(cam3, "lwuflx")[1]) < 0.2 => true
    end
end
end