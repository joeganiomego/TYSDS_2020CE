"""
Historic and Instant Preston Curves
Approximating LEB sole function of Tech Level.
Angus Maddison 1990 GDPPC being equal to TL = technical level (*), 
LEB = life expectancy at birth.
Code to be run on Julia (1.4.1, julialang.org).

Turin,Italy 20/June/2020 by Joe GANIO-MEGO, gganio.geo@yahoo.com, +39 3495472346.

(*) technical level as proportional to Angus Maddison 1990 GDPPC values (TLAM1990GDPPC)

"""

using Dierckx
using LsqFit
#using DifferentialEquations
#using DiffEqParamEstim

#using LeastSquaresOptim
#using Optim

using Plots


# Common data, Julia data style 

counter = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

year = [-2000.0, -1500.0, -1000.0, -500.0, 0.0, 500.0, 1000.0, 1300.0, 1500.0, 1600.0, 1700.0, 1820.0, 1870.0, 1913.0, 1950.0, 1973.0, 1998.0, 2010.0, 2020.0]

LEB = [26.0, 26.1, 26.2, 26.3, 26.4, 26.5, 26.6, 26.7, 26.8, 26.9, 27.0, 29.1, 29.7, 34.1, 45.73, 58.37, 65.74, 69.93, 72.78]

Popu = [87623823, 111239718, 141220440, 179281404, 225810000, 246560000, 267330000, 352860000, 438428000, 556148000, 603490000, 1041707970, 1275732067, 1792924702, 2527959894, 3922793352, 5921366350, 6933000000, 7794798739]

Popu2 = [77826427, 98193399, 124002434, 164868119, 239920070, 273322571, 335291843, 434021668, 480506023, 550878347, 643618799, 1005302902, 1327318152, 1827661021, 2579434878, 3859732004, 5912061676, 6933000000, 7794798739]

GDPPC = [127, 143, 151, 239, 467, 461, 453, 508, 566, 596, 615, 666, 870, 1524, 2111, 4083, 5709, 7814, 8694]

##################
# Info
#
# Technarian (Age) = year[12:19] , year[1820:2020]
# Extended_Technarian (Age) = year[11:19] , year[1700:2020]
# Popu_Technarian (Age) = year[9:19] , year[1500:2020]
# Pre_Technarian (Age) = year[9:15] , year[1500:1950]
# Extended_Pre_Technarian (Age) = year[8:15] , year[1300:1950]
####################
#
# Common data Normalized to values at CE=-2000.0

LEB_Norm=LEB/LEB[1]
Popu_Norm=Popu/Popu[1]
Popu2_Norm=Popu2/Popu2[1]
GDPPC_Norm=GDPPC/GDPPC[1]

#############################
#Spline building with dierckx.  1 smooth spline and 1 linear spline.
# Spline1D(x, y; w=ones(length(x)), k=3, bc="nearest", s=0.0)
# Spline1D(x, y, xknots; w=ones(length(x)), k=3, bc="nearest")
# derivative(spl, x; nu=1)

#Smooth Spline
spl_Technarian_p = Spline1D(year[9:19], Popu_Norm[9:19] ; k=3,bc="nearest", s=2.15)
spl_Technarian_g = Spline1D(year[9:19], GDPPC_Norm[9:19] ; k=5, bc="nearest", s=30.50)
#spl_Technarian_l = Spline1D(year[8:19], LEB_Norm[8:19]; k=5, bc="nearest", s=0.021)
spl_Technarian_lb = Spline1D(year[8:19], LEB_Norm[8:19]; k=3, bc="nearest", s=0.023)

#Linear Spline
spl_Pre_Technarian_p = Spline1D(year[9:19], Popu_Norm[9:19]; k=1, bc="nearest" , s=0.0)
spl_Pre_Technarian_g = Spline1D(year[9:19], GDPPC_Norm[9:19]; k=1, bc="nearest" , s=0.0)
#spl_Pre_Technarian_l = Spline1D(year[8:19], LEB_Norm[8:19]; k=1, bc="nearest", s=0.0)
spl_Pre_Technarian_lb = Spline1D(year[8:19], LEB_Norm[8:19]; k=1, bc="nearest", s=0.0)

#######################################################################################
#
# Exponential functions that match to world Hystoric Preston Curve are identified here, with least squares method.
# LsqFit matching to Scatter data

x = GDPPC_Norm[10:19]
y= LEB_Norm[10:19]

@.model1(t, p) = (p[1] *(t)^(p[2]))+p[3]
 p01 = [0.5, 0.30,-1.0]
 fit1 = curve_fit(model1, x, y, p01)
 param1 = fit1.param

@.model3(t, p) = (p[1] *(t)^(p[2]))+p[3]
w = (LEB_Norm[10:19]/LEB_Norm[10])
 p03 = [0.5,0.5, -1.0]
 fit3 = curve_fit(model3, x, y, (w), p03)
 param3 = fit3.param

########################################################################################
#Best Retained Match = model5
# LEBNorm =(1.27814 *(GDPPCNorm)^(0.252439))-0.885932
# function LEBNorm(GDPPCNorm)
#    return (1.27814 *(GDPPCNorm)^(0.252439))-0.885932
#end

@.model5(t, p) = (p[1] *(t)^(p[2]))+p[3]
@.w = (LEB_Norm[10:19]/LEB_Norm[10])^0.5
w5 = w
 p05 = [0.5,0.5,0.5]
 fit5 = curve_fit(model5, x, y, (w), p05)
 param5 = fit5.param

#######################################################################################
# LsqFit matching to Dierck Splines data

sd = collect(year[10]:10.0:year[19])
sx =  evaluate(spl_Technarian_g, sd)
sy = evaluate(spl_Technarian_lb, sd)

@.model2(tt, pp) = (pp[1] *(tt)^(pp[2]))+pp[3]
w = (sd/year[10])
 p02 = [0.5, 0.4,-1.0]
 fit2 = curve_fit(model2, sx, sy, w, p02)
 param2 = fit2.param

@.model4(tt, pp) = (pp[1] *(tt)^(pp[2]))+pp[3]
 p04 = [0.5, 0.4,-1.0]
 fit4 = curve_fit(model4, sx, sy, p04)
 param4 = fit4.param

#######################################################################################
#######################################################################################
# Color Plot of LEB function of Tech Level, all normalized

sxx = collect(year[9]:5.0:year[19])

plt1 = plot(layout=(1,1),title="LEB function of Tech Level ,Normalized")

scatter!(plt1[1],(GDPPC_Norm[9:19]),(LEB_Norm[9:19]),label="TechLevel_LEB_Norm Real")

plot!(plt1[1],(evaluate(spl_Technarian_g, sxx),evaluate(spl_Technarian_lb, sxx)),linewidth=1,label="TechLevel_LEB Smooth Spline")
plot!(plt1[1],(evaluate(spl_Pre_Technarian_g, sxx),evaluate(spl_Pre_Technarian_lb, sxx)),linewidth=1,label="TechLevel_LEB Linear Spline")

sxg = collect(GDPPC_Norm[9]:1.0:GDPPC_Norm[19])

plot!(plt1[1],(sxg),(model1(sxg, param1)),linewidth=1,label="model1")
plot!(plt1[1],(sxg),(model2(sxg, param2)),linewidth=1,label="model2")
plot!(plt1[1],(sxg),(model3(sxg, param3)),linewidth=1,label="model3")
plot!(plt1[1],(sxg),(model4(sxg, param4)),linewidth=1,label="model4")
plot!(plt1[1],(sxg),(model5(sxg, param5)),linewidth=4,label="model5, (retained model)",
legend=:bottomright,xaxis="[TechLevel = GDPPC_Norm]",yaxis="[LEB_Norm]",show = true)

#read(stdin, Char)
savefig("Preston_Historic_Instant_01_1")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("Preston_Historic_Instant_01_1big")

#######################################################################################
# BW Plot of LEB function of Tech Level, all normalized

plt1 = plot(layout=(1,1),size=(600,400),title="LEB function of Tech Level ,Normalized")

scatter!(plt1[1],(GDPPC_Norm[9:19]),(LEB_Norm[9:19]),label="TechLevel_LEB_Norm Real",color=:black)

plot!(plt1[1],(evaluate(spl_Technarian_g, sxx),evaluate(spl_Technarian_lb, sxx)),linewidth=1,label="TechLevel_LEB Smooth Spline", color=:black,linestyle=:dot)
plot!(plt1[1],(evaluate(spl_Pre_Technarian_g, sxx),evaluate(spl_Pre_Technarian_lb, sxx)),linewidth=1,label="TechLevel_LEB Linear Spline", color=:black,linestyle=:dot)

plot!(plt1[1],(sxg),(model1(sxg, param1)),linewidth=1,label="model1", color=:black,linestyle=:dash)
plot!(plt1[1],(sxg),(model2(sxg, param2)),linewidth=1,label="model2", color=:black,linestyle=:dash)
plot!(plt1[1],(sxg),(model3(sxg, param3)),linewidth=1,label="model3", color=:black,linestyle=:dash)
plot!(plt1[1],(sxg),(model4(sxg, param4)),linewidth=1,label="model4", color=:black,linestyle=:dash)
plot!(plt1[1],(sxg),(model5(sxg, param5)),linewidth=4,label="model5, (retained model)",
legend=:bottomright,xaxis="[TechLevel = GDPPC_Norm]",yaxis="[LEB_Norm]",color=:black,linestyle=:solid,show = true)

#read(stdin, Char)
savefig("Preston_Historic_Instant_01bw_1")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("Preston_Historic_Instant_01bw_1big")


###################################################################
# Instant Preston Curve Functions

#https://en.wikipedia.org/wiki/Preston_curve
# CE 2005 Intl Dollars
# 1.5 deflator roughly estimated from here:https://www.in2013dollars.com/  additional info here https://data.worldbank.org/
function PrestonLEB_a(r)
    return @.(6.6354*(log(r*1.5))+10.754)
end

#https://blog.euromonitor.com/economic-growth-and-life-expectancy-do-wealthier-countries-live-longer/
# CE 2012 Intl Dollars
# 1.6 deflator roughly estimated from here:https://www.in2013dollars.com/  additional info here https://data.worldbank.org/
function PrestonLEB_b(r)
    return @.(6.0406*(log(r*1.6))+16.132)
end

#https://blog.huawei.com/2020/09/04/ict-sustainable-development-through-a-preston-curve-lens/
# CE 2009 Intl Dollars
# 1.55 deflator roughly estimated from here:https://www.in2013dollars.com/  additional info here https://data.worldbank.org/
function PrestonLEB_c(r)
    return @.(7.02*(log(r*1.55))+6.9)
end

#######################################################################
# Instant Preston Approximated Power 0.15 Function
# Manual approximation of all above functions with a power 0.15 function

function PrestonLEB_Approx(r)
    return @.((1.27814 *(r/700)^(0.15))+0.95)*LEB[1]
end


#######################################################################
# Historic Preston Curve Function

function LEBmodel5(r)
    return @.((1.27814 *(r/GDPPC[1])^(0.252439))-0.885932)*LEB[1]
end

######################################################################
# Color Plot of Preston Curves, Historic Preston Curve and Instant Preston Curve

plt2 = plot(layout=(1,1),title="LEB function of Tech Level, Preston Curves")

scatter!(plt2[1],(GDPPC[9:19]),(LEB[9:19]),label="TechLevel_LEB Real")

psxg = collect(GDPPC[9]:100.0:15000)

plot!(plt2[1],(psxg),LEBmodel5(psxg),linewidth=4,label="LEBmodel5, Historic Preston")

plot!(plt2[1],(psxg),PrestonLEB_c(psxg),linewidth=1,label="Instant Preston c, CE 2009")
plot!(plt2[1],(psxg),PrestonLEB_b(psxg),linewidth=1,label="Instant Preston b, CE 2012")
plot!(plt2[1],(psxg),PrestonLEB_a(psxg),linewidth=1,label="Instant Preston a, CE 2005")
plot!(plt2[1],(psxg),PrestonLEB_Approx(psxg),linewidth=1,label="Instant Preston 0.15 Approx")

plot!(plt2[1],[7814.0,7814.0],[50.0,80.0],linewidth=2,label="CE 2010",
legend=:bottomright,xaxis="[Angus Maddison GDPPC_CE1990IntlD]",yaxis="[LEB years]",xticks = 0:2000:15000,show = true)

#read(stdin, Char)
savefig("Preston_Historic_Instant_01_2")

plot!(plt2[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("Preston_Historic_Instant_01_2big")


######################################################################
# BW Plot of Preston Curves, Historic Preston Curve and Instant Preston Curve

plt2 = plot(layout=(1,1),size=(600,400),title="LEB function of Tech Level, Preston Curves")

scatter!(plt2[1],(GDPPC[9:19]),(LEB[9:19]),label="TechLevel_LEB Real",color=:black)

psxg = collect(GDPPC[9]:100.0:15000)

plot!(plt2[1],(psxg),LEBmodel5(psxg),linewidth=4,label="LEBmodel5, Historic Preston",color=:black,linestyle=:solid)

plot!(plt2[1],(psxg),PrestonLEB_c(psxg),linewidth=1,label="Instant Preston c, CE 2009",color=:black,linestyle=:dot)
plot!(plt2[1],(psxg),PrestonLEB_b(psxg),linewidth=1,label="Instant Preston b, CE 2012",color=:black,linestyle=:dash)
plot!(plt2[1],(psxg),PrestonLEB_a(psxg),linewidth=1,label="Instant Preston a, CE 2005",color=:black,linestyle=:dashdot)
plot!(plt2[1],(psxg),PrestonLEB_Approx(psxg),linewidth=1,label="Instant Preston 0.15 Approx",color=:black,linestyle=:solid)

plot!(plt2[1],[7814.0,7814.0],[50.0,80.0],linewidth=2,label="CE 2010",color=:black,linestyle=:dashdotdot,
legend=:bottomright,xaxis="[Angus Maddison GDPPC_CE1990IntlD]",yaxis="[LEB years]",xticks = 0:2000:15000,show = true)

#read(stdin, Char)
savefig("Preston_Historic_Instant_01bw_2")

plot!(plt2[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("Preston_Historic_Instant_01bw_2big")

##########################################################
# Printing parameters of LsqFit

println(param1 ,param2, param3, param4, param5)

println("   ")

println("Best Retained Match = model5  LEBNorm =(1.27814 *(GDPPCNorm)^(0.252439))-0.885932")  
println(param5)


#######################################################################
# WRITE # WRITE # WRITE


outfile = "Preston_Hist_Inst_Results.txt"

f = open(outfile, "w") 

println(f,"   ")
println(f,"LEBNorm =(p[1] *(GDPPCNorm)^(p[2]))-p[3]")  
println(f,"param1 = ", param1)
println(f,"param2 = ", param2)
println(f,"   ")

println(f,"param3 = ", param3)
println(f,"param4 = ", param4)
println(f,"   ")

println(f,"param5 = ", param5)
println(f,"   ")
println(f,"Best Retained Match = model5  LEBNorm =(1.27814 *(GDPPCNorm)^(0.252439))-0.885932") 
println(f,"   ")

close(f)

##########################

#end
