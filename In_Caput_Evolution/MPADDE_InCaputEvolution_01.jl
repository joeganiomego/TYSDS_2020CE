"""
Early allometric model of homo sapiens species.
Based on data from CE -2000 to CE 2020 predicts human population , technical level (*), life expectancy at birth, healthy life expectancy (**) and age of training end (**) from the year CE 2020 to CE 6000.
Code run on Julia (1.5.3, julialang.org).

Turin,Italy 20/June/2020 by Giuseppe GANIO-MEGO, gganio.geo@yahoo.com, +39 3495472346.

(*) technical level as proportional to Angus Maddison 1990 GDPPC values (TLAM1990GDPPC)
(**) preliminary attempt.

"""
###############################################à
# Lsqfit used to identify core "technarian jump" S-curve parameters.
# The Technarian Age Jump is fuelled by
# IN CAPUT EVOLUTION
# Manual Parameter Adjustment used to identify additional (speculative) parameters to extend validity of equations outside of core 
# Technarian Jump
# Note:
# About 14 allometric functions candidates have been tested and veryfied.  Only 6 of them have been retained as worthy to be applied.

using Dierckx
using LsqFit
using DifferentialEquations
#using DiffEqParamEstim

#using DataFrames
using CSV

using Plots


# Common input data

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
"""
interpolate scattered data of the Technarian Age Jump (TAJ) to get derivative input data.
Model is based on the assumption that during the TAJ all movement of the species is intended to 
maximaze dpopu/dt and dTL/dt

"""
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
# Main Allometric Function is identified here. 
# several function candidates have been tested.  Only 2 performed as well as beeing worthy to be retained.
# since the IN CAPUT EVOLUTION mechanism was based on deivatives the allometry functions are identified by
# LsqFit matching to derived splines from Scatter data, limited to the Technarian age, where the current S-Curve 
# jump took place involving also LEB.

"""
maximization of dpopu/dt and dTL/dt happens with polynomial beneficial effect from popu and TL.
While the negative effect is due to the increasing difficulty being either exponetial or double exponetial
(which covers also factoria , stirling).
The beneficial effect of longer LEB is covered by direct relationship between TL and LEB (Preston Curve)

#https://en.wikipedia.org/wiki/Time_complexity

"""

sd = collect(year[10]:5.0:year[19])
popudata = derivative(spl_Technarian_p, sd; nu=1)
gdppcdata = derivative(spl_Technarian_g, sd; nu=1)
#lebdata = derivative(spl_Technarian_lb, sd; nu=1)  # to be used in case of allometric models avoiding GDPPC data but using LEB instead.
xdata = [ evaluate(spl_Technarian_p, sd)  evaluate(spl_Technarian_g, sd)   evaluate(spl_Technarian_lb, sd) ]  


# Exponential Style Allometry Equation 
# popu derivative with polynomial+exponential technical complexity
@. multimodelp1(x, p) = p[1] *x[:, 1]^(p[2])*(x[:, 2]^(p[3]))*2^(x[:, 2]*p[4])
 pp01 = [0.001, 1.5,0.5,-0.01]
 fitp1 = curve_fit(multimodelp1, xdata, popudata, pp01)
 paramp1 = fitp1.param
 sep1 = stderror(fitp1)

# Stirling Style Allometry Equation (it is a double exponential and covers factorial style equations)
# popu derivative with polynomial+stirling technical complexity
@. multimodelp2(x, p) = p[1] *x[:, 1]^(p[2])*(x[:, 2]^(p[3]))*(x[:, 2]^(x[:, 2]*p[4]))
 pp02 = [0.001, 1.5,0.5,-0.01]
 fitp2 = curve_fit(multimodelp2, xdata, popudata, pp02)
 paramp2 = fitp2.param
 sep2 = stderror(fitp2)


# Exponential Style Allometry Equation 
# GDPPC derivative with polynomial+exponential technical complexity
@. multimodelg1(x, p) = p[1] *x[:, 1]^(p[2])*x[:, 2]^(p[3])*2^(x[:, 2]*p[4])
 pg01 = [0.001, 1.5,0.5,-0.01]
 fitg1 = curve_fit(multimodelg1, xdata, gdppcdata, pg01)
 paramg1 = fitg1.param
 seg1 = stderror(fitg1)
 
 
# Stirling Style Allometry Equation (it is a double exponential and covers factorial style equations)
# GDPPC derivative with polynomial+stirling technical complexity
@. multimodelg2(x, p) = p[1] *x[:, 1]^(p[2])*(x[:, 2]^(p[3]))*(x[:, 2]^(x[:, 2]*p[4]))
 pg02 = [0.001, 1.5,0.5,-0.01]
 fitg2 = curve_fit(multimodelg2, xdata, gdppcdata, pg02)
 paramg2 = fitg2.param
 seg2 = stderror(fitg2)


###################################
# Parameter contribution to derivative
# paramp1[1] = general multiplier for popu derivative with polynomial+exponential technical complexity
# paramp1[2] =*x[1]^(paramp1[2])* exponent for popu factor of popu derivative with polynomial+exponential technical complexity
# paramp1[3] = (x[2]^(paramp1[3])) exponent for GDPPC factor of popu derivative with polynomial+exponential technical complexity
# paramp1[4] = *2^(x[2]*(paramp1[4])) exponent multiplier for GDPPC exponent of popu derivative with polynomial+exponential technical complexity
# paramg1[1] = general multiplier for GDPPC derivative with polynomial+exponential technical complexity
# paramg1[2] =*x[1]^(paramp1[2])* exponent for popu factor of GDPPC derivative with polynomial+exponential technical complexity
# paramg1[3] = (x[2]^(paramp1[3])) exponent for GDPPC factor of GDPPC derivative with polynomial+exponential technical complexity
# paramg1[4] = *2^(x[2]*(paramp1[4])) exponent multiplier for GDPPC exponent of GDPPC derivative with polynomial+exponential technical complexity

"""

@. multimodelg1paramp11(x, p) = p[p11] 
@. multimodelg1paramp12(x, p) = x^(p[p12])
@. multimodelg1paramp13(x, p) = x^(p[p13])
@. multimodelg1paramp14(x, p) = 2^(x*p[p14])

@. multimodelg1paramp21(x, p) = p[p21] 
@. multimodelg1paramp22(x, p) = x^(p[p22])
@. multimodelg1paramp23(x, p) = x^(p[p23])
@. multimodelg1paramp24(x, p) = (x^(x*p[p24]))

@. multimodelg1paramg11(x, p) = p[g11] 
@. multimodelg1paramg12(x, p) = x^(p[g12])
@. multimodelg1paramg13(x, p) = x^(p[g13])
@. multimodelg1paramg14(x, p) = 2^(x*p[g14])

@. multimodelg1paramg21(x, p) = p[g21] 
@. multimodelg1paramg22(x, p) = x^(p[g22])
@. multimodelg1paramg23(x, p) = x^(p[g23])
@. multimodelg1paramg24(x, p) = (x^(x*p[g24]))

"""

##########################################################
# ODE section
# ODE from smooth Dierckxs Splines
# Note:
# Here are identidied the core  allometric parametrs for the 2 kind of proposed functions. About 14 allometric functions candidates have been tested and veryfied.  Only 6 of them have been retained as worthy to be applied.
# the following (model1 and 2) are the simplest ones.
#
#x[1]=Popu (normalized)
#x[2]=GDPPC (normalized)


function PopuGDPPCnorm1(dx,x,p,t)
# popu derivative with polynomial+exponential technical complexity
dx[1] = paramp1[1]*x[1]^(paramp1[2])*(x[2]^(paramp1[3]))*2^(x[2]*(paramp1[4]))
# GDPPC derivative with polynomial+exponential technical complexity
dx[2] = paramg1[1]*x[1]^(paramg1[2])*(x[2]^(paramg1[3]))*2^(x[2]*(paramg1[4]))
end

function PopuGDPPCnorm2(dx,x,p,t)
# popu derivative with polynomial+stirling technical complexity
dx[1] = paramp2[1]*x[1]^(paramp2[2])*(x[2]^(paramp2[3]))*(x[2]^(x[2]*(paramp2[4])))
# GDPPC derivative with polynomial+stirling technical complexity
dx[2] = paramg2[1]*x[1]^(paramg2[2])*(x[2]^(paramg2[3]))*(x[2]^(x[2]*(paramg2[4])))
end

x0 = [Popu_Norm[19] ,  GDPPC_Norm[19]]
tspan = (year[19],4000.0)
tspanb = (year[19],0.0)

# Model1
PopuGDPPC1 = ODEProblem(PopuGDPPCnorm1,x0,tspan)
PopuGDPPC1b = ODEProblem(PopuGDPPCnorm1,x0,tspanb)
PopuGDPPCsol1 = solve(PopuGDPPC1)
PopuGDPPCsol1b = solve(PopuGDPPC1b)

# Model2
PopuGDPPC2 = ODEProblem(PopuGDPPCnorm2,x0,tspan)
PopuGDPPC2b = ODEProblem(PopuGDPPCnorm2,x0,tspanb)
PopuGDPPCsol2 = solve(PopuGDPPC2)
PopuGDPPCsol2b = solve(PopuGDPPC2b)

######################################
#  HI PRECISION
# https://docs.sciml.ai/stable/basics/faq/#Numerical-Error-1
#no need to use high precision since it changes the results of 0.000001.
#good hi precision
#PopuGDPPCsol1hp = solve(PopuGDPPC1,Tsit5(),abstol=1e-12,reltol=1e-12)
#PopuGDPPCsol1bhp = solve(PopuGDPPC1b,Tsit5(),abstol=1e-12,reltol=1e-12)
#good hi precision
#PopuGDPPCsol2hp = solve(PopuGDPPC2,Tsit5(),abstol=1e-12,reltol=1e-12)
#PopuGDPPCsol2bhp = solve(PopuGDPPC2b,Tsit5(),abstol=1e-12,reltol=1e-12)

#################################
# ODE with Manual Parameter Estimation
# Manual Parameter Adjustment used to identify additional (speculative) parameters to extend validity of equations outside of core 
# Technarian Jump
# In this case (model 3 and 4) it is supposed that at low crowding (low x[1]) there is low focus on increasing the capacity
# to accept higher x[1].
# At higher tech level (GDDPC = TL = x[2]) it is supposed that it becames more and more difficlut to icrease x[2].
"""
It is supposed that the previously identified ODEs would simulate species evolution only during the Technarian Age jump.  Outside
of the TAJ there is for sure a reduction in the GDPPC dynamic and unknown situation for the Population factor.
Indeed the original ODEs under predict at lower years , such as 0.0 CE year.
Adjusting the new added equation parameters (manually) shows that outside the TAJ both Popu and GDPPC factors need to be reduced.

"""
#
#x[1]=Popu (normalized)
#x[2]=GDPPC (normalized)

function PopuGDPPCnorm3(dx,x,p,t)
# popu derivative with polynomial+exponential technical complexity + Popu density correction
dx[1] = paramp1[1]*x[1]^(paramp1[2])*(x[2]^(paramp1[3]))*2^(x[2]*(paramp1[4]))+(-0.0508)*(1/(1+x[1]))
# GDPPC derivative with polynomial+exponential technical complexity + GDPPC level correction
dx[2] = paramg1[1]*x[1]^(paramg1[2])*(x[2]^(paramg1[3]))*2^(x[2]*(paramg1[4]))+(-0.02575)*(x[2]/(7.581+x[2]))
end

function PopuGDPPCnorm4(dx,x,p,t)
# popu derivative with polynomial+stirling technical complexity + Popu density correction
dx[1] = paramp2[1]*x[1]^(paramp2[2])*(x[2]^(paramp2[3]))*(x[2]^(x[2]*(paramp2[4])))+(-0.0585)*(1/(1+x[1]))
# GDPPC derivative with polynomial+stirling technical complexity + GDPPC level correction
dx[2] = paramg2[1]*x[1]^(paramg2[2])*(x[2]^(paramg2[3]))*(x[2]^(x[2]*(paramg2[4])))+(-0.0258)*(x[2]/(6.92+x[2]))
end



x0 = [Popu_Norm[19] ,  GDPPC_Norm[19]]
tspan = (year[19],6000.0)
tspanb = (year[19],-2000.0)


# Note: Solver BS3() chosen because only alg that delivered smooth lines.
# However all algoritms had the same overall results.

# Model3, the callback is done to get the derivative values of solution.
PopuGDPPC3 = ODEProblem(PopuGDPPCnorm3,x0,tspan)
PopuGDPPC3b = ODEProblem(PopuGDPPCnorm3,x0,tspanb)
saved_values3 = SavedValues(Float64, AbstractArray{Float64})
cb3 = SavingCallback((u,t,integrator)->integrator(t,Val{1}), saved_values3)
saved_values3b = SavedValues(Float64, AbstractArray{Float64})
cb3b = SavingCallback((u,t,integrator)->integrator(t,Val{1}), saved_values3b)
alg = BS3()
PopuGDPPCsol3 = solve(PopuGDPPC3,alg, callback=cb3)
dP3dt = [saved_values3.saveval[i][1] for i in 2:length(PopuGDPPCsol3)]
dG3dt = [saved_values3.saveval[i][2] for i in 2:length(PopuGDPPCsol3)]
PopuGDPPCsol3b = solve(PopuGDPPC3b,alg, callback=cb3b)
dP3bdt = [saved_values3b.saveval[i][1] for i in 2:length(PopuGDPPCsol3b)]
dG3bdt = [saved_values3b.saveval[i][2] for i in 2:length(PopuGDPPCsol3b)]


# Model4, the callback is done to get the derivative values of solution.
PopuGDPPC4 = ODEProblem(PopuGDPPCnorm4,x0,tspan)
PopuGDPPC4b = ODEProblem(PopuGDPPCnorm4,x0,tspanb)
saved_values4 = SavedValues(Float64, AbstractArray{Float64})
cb4 = SavingCallback((u,t,integrator)->integrator(t,Val{1}), saved_values4)
saved_values4b = SavedValues(Float64, AbstractArray{Float64})
cb4b = SavingCallback((u,t,integrator)->integrator(t,Val{1}), saved_values4b)
alg = BS3()
PopuGDPPCsol4 = solve(PopuGDPPC4,alg, callback=cb4)
dP4dt = [saved_values4.saveval[i][1] for i in 2:length(PopuGDPPCsol4)]
dG4dt = [saved_values4.saveval[i][2] for i in 2:length(PopuGDPPCsol4)]
PopuGDPPCsol4b = solve(PopuGDPPC4b,alg, callback=cb4b)
dP4bdt = [saved_values4b.saveval[i][1] for i in 2:length(PopuGDPPCsol4b)]
dG4bdt = [saved_values4b.saveval[i][2] for i in 2:length(PopuGDPPCsol4b)]



##################################################
# DDE with Manual Parameter Estimation
# In this case (model 5 and 6) it is supposed that at low crowding (low x[1]) there is low focus on increasing habitat capacity
# to accept higher x[1].
# At higher tech level (GDDPC = TL = x[2]) it is supposed that it getsvmore and more difficlut to icrease x[2].
# Additionnaly in the DDE set of equations for x[2] there is sentiment factor.  The more it grows the bigger will be the effort
# (or efficiency?) to get x[2] growing.
"""
As before but here the GDPPC growth factor is also reduced by the lack of dynamicsm.

"""
#
#x[1]=Popu (normalized)
#x[2]=GDPPC (normalized)

# tau = lag = 3 years : Delay time estimated for a perceived situation (level of x[1,2]) to trickle into a reaction.
 
tau = 3
lags = [tau]

function PopuGDPPCnorm5(dx,x,h,p,t)
# popu derivative with polynomial+exponential technical complexity + Popu density correction
dx[1] = paramp1[1]*x[1]^(paramp1[2])*(x[2]^(paramp1[3]))*2^(x[2]*(paramp1[4]))+(-0.055)*(1/(1+h(p, t-tau)[1]))
# GDPPC derivative with polynomial+exponential technical complexity + GDPPC level and sentiment correction
dx[2] = paramg1[1]*x[1]^(paramg1[2])*(x[2]^(paramg1[3]))*2^(x[2]*(paramg1[4]))+(-0.024)*(x[2]/(5.37+x[2]))*
(
((0.06)*h(p, t-tau)[2])/
(abs(x[2]-h(p, t-tau)[2])+((0.06)*h(p, t-tau)[2]))
)
end


function PopuGDPPCnorm6(dx,x,h,p,t)
# popu derivative with polynomial+stirling technical complexity + Popu density correction
dx[1] = paramp2[1]*x[1]^(paramp2[2])*(x[2]^(paramp2[3]))*(x[2]^(x[2]*(paramp2[4])))+(-0.0552)*(1/(1+h(p, t-tau)[1]))
# GDPPC derivative with polynomial+stirling technical complexity + GDPPC level and sentiment correction
dx[2] = paramg2[1]*x[1]^(paramg2[2])*(x[2]^(paramg2[3]))*(x[2]^(x[2]*(paramg2[4])))+(-0.0235)*(x[2]/(5.159+x[2]))*
(
((0.0602)*h(p, t-tau)[2])/
(abs(x[2]-h(p, t-tau)[2])+((0.0602)*h(p, t-tau)[2]))
)
end

h(p, t) = [Popu_Norm[19]-(tau/30) ,  GDPPC_Norm[19]-(tau/30)]
x0 = [Popu_Norm[19] ,  GDPPC_Norm[19]]
tspan = (year[19],6000.0)
tspanb = (year[19],-2000.0)

# Note: Solver MethodOfSteps(BS3()) chosen because only alg that delivered smooth lines.
# However all algoritms had the same overall results.

# Model5, the callback is done to get the derivative values of solution.
PopuGDPPC5 = DDEProblem(PopuGDPPCnorm5,x0,h,tspan,constant_lags = lags)
PopuGDPPC5b = DDEProblem(PopuGDPPCnorm5,x0,h,tspanb,constant_lags = -lags)
saved_values5 = SavedValues(Float64, AbstractArray{Float64})
cb5 = SavingCallback((u,t,integrator)->integrator(t,Val{1}), saved_values5)
saved_values5b = SavedValues(Float64, AbstractArray{Float64})
cb5b = SavingCallback((u,t,integrator)->integrator(t,Val{1}), saved_values5b)
alg = MethodOfSteps(BS3())
PopuGDPPCsol5 = solve(PopuGDPPC5,alg,callback=cb5)
dP5dt = [saved_values5.saveval[i][1] for i in 2:length(PopuGDPPCsol5)]
dG5dt = [saved_values5.saveval[i][2] for i in 2:length(PopuGDPPCsol5)]
PopuGDPPCsol5b = solve(PopuGDPPC5b,alg, callback=cb5b)
dP5bdt = [saved_values5b.saveval[i][1] for i in 2:length(PopuGDPPCsol5b)]
dG5bdt = [saved_values5b.saveval[i][2] for i in 2:length(PopuGDPPCsol5b)]

# Model6, the callback is done to get the derivative values of solution.
PopuGDPPC6 = DDEProblem(PopuGDPPCnorm6,x0,h,tspan,constant_lags = lags)
PopuGDPPC6b = DDEProblem(PopuGDPPCnorm6,x0,h,tspanb,constant_lags = -lags)
saved_values6 = SavedValues(Float64, AbstractArray{Float64})
cb6 = SavingCallback((u,t,integrator)->integrator(t,Val{1}), saved_values6)
saved_values6b = SavedValues(Float64, AbstractArray{Float64})
cb6b = SavingCallback((u,t,integrator)->integrator(t,Val{1}), saved_values6b)
alg = MethodOfSteps(BS3())
PopuGDPPCsol6 = solve(PopuGDPPC6,alg,callback=cb6)
dP6dt = [saved_values6.saveval[i][1] for i in 2:length(PopuGDPPCsol6)]
dG6dt = [saved_values6.saveval[i][2] for i in 2:length(PopuGDPPCsol6)]
PopuGDPPCsol6b = solve(PopuGDPPC6b,alg, callback=cb6b)
dP6bdt = [saved_values6b.saveval[i][1] for i in 2:length(PopuGDPPCsol6b)]
dG6bdt = [saved_values6b.saveval[i][2] for i in 2:length(PopuGDPPCsol6b)]


###################################
# Bulding the %Growth Values
#


dP5dtAug=vcat(dP5dt,dP5dt[end])
dG5dtAug=vcat(dG5dt,dG5dt[end])
dP5bdtAug=vcat(dP5bdt,dP5bdt[end])
dG5bdtAug=vcat(dG5bdt,dG5bdt[end])

dP6dtAug=vcat(dP6dt,dP6dt[end])
dG6dtAug=vcat(dG6dt,dG6dt[end])
dP6bdtAug=vcat(dP6bdt,dP6bdt[end])
dG6bdtAug=vcat(dG6bdt,dG6bdt[end])


GrowthP5=@.dP5dtAug/PopuGDPPCsol5[1,:]*100
GrowthTL5=@.dG5dtAug/PopuGDPPCsol5[2,:]*100
GrowthP5b=@.dP5bdtAug/PopuGDPPCsol5b[1,:]*100
GrowthTL5b=@.dG5bdtAug/PopuGDPPCsol5b[2,:]*100

GrowthP6=@.dP6dtAug/PopuGDPPCsol6[1,:]*100
GrowthTL6=@.dG6dtAug/PopuGDPPCsol6[2,:]*100
GrowthP6b=@.dP6bdtAug/PopuGDPPCsol6b[1,:]*100
GrowthTL6b=@.dG6bdtAug/PopuGDPPCsol6b[2,:]*100


GrowthGDP5=@.(GrowthP5+GrowthTL5)
GrowthGDP5b=@.(GrowthP5b+GrowthTL5b)
GrowthGDP6=@.(GrowthP6+GrowthTL6)
GrowthGDP6b=@.(GrowthP6b+GrowthTL6b)


##################################

##########################
# Automatic Parameter Fitting. 
# Not done, since manual parameter fitting proved to yield a stable enough solutions.
# Therefore automatic parameter fitting would have been cding effort for a small difference in 
# the final results.



#############################################
# Historic Preston Curve
# Prestonization of HLE and ETA.
# LEB as single function of GDPPC, see file Preston_Historic_Instant_01.jl

# Best Retained Match = Model5
# LEBNorm =(1.27814 *(GDPPCNorm)^(0.252439))-0.885932
 function LEBNorm(GDPPCNorm)
   return (1.27814 *(GDPPCNorm)^(0.252439))-0.885932
end

# hypotesis of HLE as single function of LEB
# Healthy life expectancy
 function HLE(LEB)
   return ((LEB)^(0.25*0.5))*37
end

# hypotesis of End of Training AGE as single function of LEB
# ETA
 function ETA(LEB)
   return ((LEB)^(0.25*0.5))*12
end


###############################################################################à
#  OUTPUT  # OUTPUT # OUTPUT 
####################################################
# Check Gauges Section
#  Used to do Manual Parameter Adjustement

##########################################
#denormalize ODE Simple

GDPPCNorm1= PopuGDPPCsol1[2,end]
GDPPCNorm2= PopuGDPPCsol2[2,end]

LebSol1 = LEBNorm(GDPPCNorm1)*LEB[1]
LebSol2 = LEBNorm(GDPPCNorm2)*LEB[1]

PopuFinal1 = PopuGDPPCsol1[1,end]*Popu[1]
PopuFinal2 = PopuGDPPCsol2[1,end]*Popu[1]


println("LebSol1 = " , LebSol1 , " years of age")
println("LebSol2 = " , LebSol2 , " years of age")

println("   ")

println("PopuFinal1 = ", PopuFinal1 , " people")
println("PopuFinal2 = " , PopuFinal2 , " people")

println("   ")

##########################################
#denormalize ODE Adjusted parameter

GDPPCNorm3= PopuGDPPCsol3[2,end]
GDPPCNorm4= PopuGDPPCsol4[2,end]

LebSol3 = LEBNorm(GDPPCNorm3)*LEB[1]
LebSol4 = LEBNorm(GDPPCNorm4)*LEB[1]

PopuFinal3 = PopuGDPPCsol3[1,end]*Popu[1]
PopuFinal4 = PopuGDPPCsol4[1,end]*Popu[1]


println("LebSol3 = " , LebSol3 , " years of age")
println("LebSol4 = " , LebSol4 , " years of age")

println("   ")

println("PopuFinal3 = ", PopuFinal3 , " people")
println("PopuFinal4 = " , PopuFinal4 , " people")

println("   ")

##########################################
#denormalize DDE Adjusted parameter

GDPPCNorm5= PopuGDPPCsol5[2,end]
GDPPCNorm6= PopuGDPPCsol6[2,end]

LebSol5 = LEBNorm(GDPPCNorm5)*LEB[1]
LebSol6 = LEBNorm(GDPPCNorm6)*LEB[1]

PopuFinal5 = PopuGDPPCsol5[1,end]*Popu[1]
PopuFinal6 = PopuGDPPCsol6[1,end]*Popu[1]


println("LebSol5 = " , LebSol5 , " years of age")
println("LebSol6 = " , LebSol6 , " years of age")

println("   ")

println("PopuFinal5 = ", PopuFinal5 , " people")
println("PopuFinal6 = " , PopuFinal6 , " people")

println("   ")

########################################
#denormalize High Precision
# hi precision
#GDPPCNorm1hp= PopuGDPPCsol1hp[2,end]
#GDPPCNorm2hp= PopuGDPPCsol2hp[2,end]
#LebSol1hp = LEBNorm(GDPPCNorm1hp)*LEB[1]
#LebSol2hp = LEBNorm(GDPPCNorm2hp)*LEB[1]
#PopuFinal1hp = PopuGDPPCsol1hp[1,end]*Popu[1]
#PopuFinal2hp = PopuGDPPCsol2hp[1,end]*Popu[1]
#println("LebSol1hp = " , LebSol1hp , " years of age")
#println("LebSol2hp = " , LebSol2hp , " years of age")
#println("   ")
#println("PopuFinal1hp = ", PopuFinal1hp , " people")
#println("PopuFinal2hp = " , PopuFinal2hp , " people")
#println("   ")


##################################
# Comparing measured data to ODE, DDE
# Comparing Gauges


POPUgaugeYear0_1= PopuGDPPCsol1b(0.0)[1]/Popu_Norm[5]
POPUgaugeYear1600_1= PopuGDPPCsol1b(1600.0)[1]/Popu_Norm[10]
POPUgaugeYear1820_1= PopuGDPPCsol1b(1820.0)[1]/Popu_Norm[12]
POPUgaugeYear1950_1= PopuGDPPCsol1b(1950.0)[1]/Popu_Norm[15]
POPUgaugeYear2020_1= PopuGDPPCsol1b(2020.0)[1]/Popu_Norm[19]

GDPPCgaugeYear0_1= PopuGDPPCsol1b(0.0)[2]/GDPPC_Norm[5]
GDPPCgaugeYear1600_1= PopuGDPPCsol1b(1600.0)[2]/GDPPC_Norm[10]
GDPPCgaugeYear1820_1= PopuGDPPCsol1b(1820.0)[2]/GDPPC_Norm[12]
GDPPCgaugeYear1950_1= PopuGDPPCsol1b(1950.0)[2]/GDPPC_Norm[15]
GDPPCgaugeYear2020_1= PopuGDPPCsol1b(2020.0)[2]/GDPPC_Norm[19]


POPUgaugeYear0_2= PopuGDPPCsol2b(0.0)[1]/Popu_Norm[5]
POPUgaugeYear1600_2= PopuGDPPCsol2b(1600.0)[1]/Popu_Norm[10]
POPUgaugeYear1820_2= PopuGDPPCsol2b(1820.0)[1]/Popu_Norm[12]
POPUgaugeYear1950_2= PopuGDPPCsol2b(1950.0)[1]/Popu_Norm[15]
POPUgaugeYear2020_2= PopuGDPPCsol2b(2020.0)[1]/Popu_Norm[19]

GDPPCgaugeYear0_2= PopuGDPPCsol2b(0.0)[2]/GDPPC_Norm[5]
GDPPCgaugeYear1600_2= PopuGDPPCsol2b(1600.0)[2]/GDPPC_Norm[10]
GDPPCgaugeYear1820_2= PopuGDPPCsol2b(1820.0)[2]/GDPPC_Norm[12]
GDPPCgaugeYear1950_2= PopuGDPPCsol2b(1950.0)[2]/GDPPC_Norm[15]
GDPPCgaugeYear2020_2= PopuGDPPCsol2b(2020.0)[2]/GDPPC_Norm[19]


println("   ")
println("POPUgaugeYear0_1 = " , POPUgaugeYear0_1 , " ")
println("POPUgaugeYear1600_1 = " , POPUgaugeYear1600_1 , " ")
println("POPUgaugeYear1820_1 = " , POPUgaugeYear1820_1 , " ")
println("POPUgaugeYear1950_1 = " , POPUgaugeYear1950_1 , " ")
println("POPUgaugeYear2020_1 = " , POPUgaugeYear2020_1 , " ")
println("   ")
println("GDPPCgaugeYear0_1 = " , GDPPCgaugeYear0_1 , " ")
println("GDPPCgaugeYear1600_1 = " , GDPPCgaugeYear1600_1 , " ")
println("GDPPCgaugeYear1820_1 = " , GDPPCgaugeYear1820_1 , " ")
println("GDPPCgaugeYear1950_1 = " , GDPPCgaugeYear1950_1 , " ")
println("GDPPCgaugeYear2020_1 = " , GDPPCgaugeYear2020_1 , " ")
println("   ")

println("   ")
println("POPUgaugeYear0_3 = " , PopuGDPPCsol3b(0.0)[1]/Popu_Norm[5] , " ")
println("POPUgaugeYear1600_3 = " , PopuGDPPCsol3b(1600.0)[1]/Popu_Norm[10] , " ")
println("POPUgaugeYear1820_3 = " , PopuGDPPCsol3b(1820.0)[1]/Popu_Norm[12] , " ")
println("POPUgaugeYear1950_3 = " , PopuGDPPCsol3b(1950.0)[1]/Popu_Norm[15] , " ")
println("POPUgaugeYear2020_3 = " , PopuGDPPCsol3b(2020.0)[1]/Popu_Norm[19] , " ")
println("   ")
println("GDPPCgaugeYear0_3 = " , PopuGDPPCsol3b(0.0)[2]/GDPPC_Norm[5] , " ")
println("GDPPCgaugeYear1600_3 = " , PopuGDPPCsol3b(1600.0)[2]/GDPPC_Norm[10] , " ")
println("GDPPCgaugeYear1820_3 = " , PopuGDPPCsol3b(1820.0)[2]/GDPPC_Norm[12] , " ")
println("GDPPCgaugeYear1950_3 = " , PopuGDPPCsol3b(1950.0)[2]/GDPPC_Norm[15] , " ")
println("GDPPCgaugeYear2020_3 = " , PopuGDPPCsol3b(2020.0)[2]/GDPPC_Norm[19] , " ")
println("   ")


println("   ")
println("POPUgaugeYear0_5 = " , PopuGDPPCsol5b(0.0)[1]/Popu_Norm[5] , " ")
println("POPUgaugeYear1600_5 = " , PopuGDPPCsol5b(1600.0)[1]/Popu_Norm[10] , " ")
println("POPUgaugeYear1820_5 = " , PopuGDPPCsol5b(1820.0)[1]/Popu_Norm[12] , " ")
println("POPUgaugeYear1950_5 = " , PopuGDPPCsol5b(1950.0)[1]/Popu_Norm[15] , " ")
println("POPUgaugeYear2020_5 = " , PopuGDPPCsol5b(2020.0)[1]/Popu_Norm[19] , " ")
println("   ")
println("GDPPCgaugeYear0_5 = " , PopuGDPPCsol5b(0.0)[2]/GDPPC_Norm[5] , " ")
println("GDPPCgaugeYear1600_5 = " , PopuGDPPCsol5b(1600.0)[2]/GDPPC_Norm[10] , " ")
println("GDPPCgaugeYear1820_5 = " , PopuGDPPCsol5b(1820.0)[2]/GDPPC_Norm[12] , " ")
println("GDPPCgaugeYear1950_5 = " , PopuGDPPCsol5b(1950.0)[2]/GDPPC_Norm[15] , " ")
println("GDPPCgaugeYear2020_5 = " , PopuGDPPCsol5b(2020.0)[2]/GDPPC_Norm[19] , " ")
println("   ")


println("   ")
println("POPUgaugeYear0_2 = " , POPUgaugeYear0_2 , " ")
println("POPUgaugeYear1600_2 = " , POPUgaugeYear1600_2 , " ")
println("POPUgaugeYear1820_2 = " , POPUgaugeYear1820_2 , " ")
println("POPUgaugeYear1950_2 = " , POPUgaugeYear1950_2 , " ")
println("POPUgaugeYear2020_2 = " , POPUgaugeYear2020_2 , " ")
println("   ")
println("GDPPCgaugeYear0_2 = " , GDPPCgaugeYear0_2 , " ")
println("GDPPCgaugeYear1600_2 = " , GDPPCgaugeYear1600_2 , " ")
println("GDPPCgaugeYear1820_2 = " , GDPPCgaugeYear1820_2 , " ")
println("GDPPCgaugeYear1950_2 = " , GDPPCgaugeYear1950_2 , " ")
println("GDPPCgaugeYear2020_2 = " , GDPPCgaugeYear2020_2 , " ")
println("   ")


println("   ")
println("POPUgaugeYear0_4 = " , PopuGDPPCsol4b(0.0)[1]/Popu_Norm[5] , " ")
println("POPUgaugeYear1600_4 = " , PopuGDPPCsol4b(1600.0)[1]/Popu_Norm[10] , " ")
println("POPUgaugeYear1820_4 = " , PopuGDPPCsol4b(1820.0)[1]/Popu_Norm[12] , " ")
println("POPUgaugeYear1950_4 = " , PopuGDPPCsol4b(1950.0)[1]/Popu_Norm[15] , " ")
println("POPUgaugeYear2020_4 = " , PopuGDPPCsol4b(2020.0)[1]/Popu_Norm[19] , " ")
println("   ")
println("GDPPCgaugeYear0_4 = " , PopuGDPPCsol4b(0.0)[2]/GDPPC_Norm[5] , " ")
println("GDPPCgaugeYear1600_4 = " , PopuGDPPCsol4b(1600.0)[2]/GDPPC_Norm[10] , " ")
println("GDPPCgaugeYear1820_4 = " , PopuGDPPCsol4b(1820.0)[2]/GDPPC_Norm[12] , " ")
println("GDPPCgaugeYear1950_4 = " , PopuGDPPCsol4b(1950.0)[2]/GDPPC_Norm[15] , " ")
println("GDPPCgaugeYear2020_4 = " , PopuGDPPCsol4b(2020.0)[2]/GDPPC_Norm[19] , " ")
println("   ")


println("   ")
println("POPUgaugeYear0_6 = " , PopuGDPPCsol6b(0.0)[1]/Popu_Norm[5] , " ")
println("POPUgaugeYear1600_6 = " , PopuGDPPCsol6b(1600.0)[1]/Popu_Norm[10] , " ")
println("POPUgaugeYear1820_6 = " , PopuGDPPCsol6b(1820.0)[1]/Popu_Norm[12] , " ")
println("POPUgaugeYear1950_6 = " , PopuGDPPCsol6b(1950.0)[1]/Popu_Norm[15] , " ")
println("POPUgaugeYear2020_6 = " , PopuGDPPCsol6b(2020.0)[1]/Popu_Norm[19] , " ")
println("   ")
println("GDPPCgaugeYear0_6 = " , PopuGDPPCsol6b(0.0)[2]/GDPPC_Norm[5] , " ")
println("GDPPCgaugeYear1600_6 = " , PopuGDPPCsol6b(1600.0)[2]/GDPPC_Norm[10] , " ")
println("GDPPCgaugeYear1820_6 = " , PopuGDPPCsol6b(1820.0)[2]/GDPPC_Norm[12] , " ")
println("GDPPCgaugeYear1950_6 = " , PopuGDPPCsol6b(1950.0)[2]/GDPPC_Norm[15] , " ")
println("GDPPCgaugeYear2020_6 = " , PopuGDPPCsol6b(2020.0)[2]/GDPPC_Norm[19] , " ")
println("   ")



#
################################################# 
# Listing of equations parameters

println("polynomial+exponential technical complexity")  
println("paramp1 = ", paramp1)
println("paramg1 = ", paramg1)
println("polynomial+stirling technical complexity")  
println("paramp2 = ", paramp2)
println("paramg2 = ", paramg2)

# Listing of LSQ fit standard Errors
println("   ")

println("polynomial+exponential stderror")  
println("sep1 = ",sep1) 
println("seg1 = ",seg1) 
println("polynomial+stirling stderror")  
println("sep2 = ",sep2)  
println("seg2 = ",seg2)


#polynomial+exponential technical complexity
#[0.000809814, 0.331911, 2.04555, -0.0651806][0.00068161, 0.53622, 1.64428, -0.0429536]
#polynomial+stirling technical complexity
#[0.00111146, 0.27198, 1.86948, -0.0083154][0.000837945, 0.515811, 1.50476, -0.00541452]


#######################################################################
# WRITE # WRITE # WRITE
#################################
# write .csv of ODE, DDE


outfile = "InCaputEvolutionResults.txt"

f = open(outfile, "w") 

println(f,"   ")
println(f,"polynomial+exponential technical complexity")  
println(f,"paramp1 = ", paramp1)
println(f,"paramg1 = ", paramg1)
println(f,"   ")
println(f,"polynomial+stirling technical complexity")  
println(f,"paramp2 = ", paramp2)
println(f,"paramg2 = ", paramg2)
println(f,"   ")
close(f)


##########################

CSV.write("PopuGDPPCsol1.csv",PopuGDPPCsol1)
CSV.write("PopuGDPPCsol1b.csv",PopuGDPPCsol1b)
CSV.write("PopuGDPPCsol2.csv",PopuGDPPCsol2)
CSV.write("PopuGDPPCsol2b.csv",PopuGDPPCsol2b)
CSV.write("PopuGDPPCsol3.csv",PopuGDPPCsol3)
CSV.write("PopuGDPPCsol3b.csv",PopuGDPPCsol3b)
CSV.write("PopuGDPPCsol4.csv",PopuGDPPCsol4)
CSV.write("PopuGDPPCsol4b.csv",PopuGDPPCsol4b)
CSV.write("PopuGDPPCsol5.csv",PopuGDPPCsol5)
CSV.write("PopuGDPPCsol5b.csv",PopuGDPPCsol5b)
CSV.write("PopuGDPPCsol6.csv",PopuGDPPCsol6)
CSV.write("PopuGDPPCsol6b.csv",PopuGDPPCsol6b)


##############################################################
# PLOT SECTION # PLOT SECTION # PLOT SECTION # PLOT SECTION
##############################################################
# plot Population Normalized
# THE POPULATION Normalized


plt1 = plot(layout=(1,1),title="World Population Normalized",size=(600,400))

scatter!(plt1[1],year[1:19],(Popu_Norm[1:19]),color=:black,label="Historic Data")

#############################
# Plot ODE from smooth Splines
# Population

plot!(plt1[1],PopuGDPPCsol1.t,PopuGDPPCsol1[1,:],linewidth=1,color=:yellow,label="Model 1")
plot!(plt1[1],PopuGDPPCsol1b.t,PopuGDPPCsol1b[1,:],linewidth=1,color=:yellow,primary=false )

plot!(plt1[1],PopuGDPPCsol2.t,PopuGDPPCsol2[1,:],linewidth=1,color=:brown,label="Model 2")
plot!(plt1[1],PopuGDPPCsol2b.t,PopuGDPPCsol2b[1,:],linewidth=1,color=:brown,primary=false )


###################
# Plot ODE from manually estimated parameters simple
# Population

plot!(plt1[1],PopuGDPPCsol3.t,PopuGDPPCsol3[1,:],linewidth=1,color=:red,label="Model 3")
plot!(plt1[1],PopuGDPPCsol3b.t,PopuGDPPCsol3b[1,:],linewidth=1,color=:red,primary=false )

plot!(plt1[1],PopuGDPPCsol4.t,PopuGDPPCsol4[1,:],linewidth=1,color=:orange,label="Model 4")
plot!(plt1[1],PopuGDPPCsol4b.t,PopuGDPPCsol4b[1,:],linewidth=1,color=:orange,primary=false )

###################
# Plot DDE from manually estimated parameters
# Popu

plot!(plt1[1],PopuGDPPCsol5.t,PopuGDPPCsol5[1,:],linewidth=1,color=:green,label="Model 5")
plot!(plt1[1],PopuGDPPCsol5b.t,PopuGDPPCsol5b[1,:],linewidth=1,color=:green,primary=false )

plot!(plt1[1],PopuGDPPCsol6.t,PopuGDPPCsol6[1,:],linewidth=1,color=:blue,label="Model 6")
plot!(plt1[1],PopuGDPPCsol6b.t,PopuGDPPCsol6b[1,:],linewidth=1,color=:blue,primary=false ,
xaxis="[CE years]",yaxis="[World Population Normalized]",xticks = -2000:(floor(6000/10)):6000,yticks = 0.0:20:160,
#xlims = (0.5,1.5),ylims = (0.5,1.5),
legend=:bottomright,grid= false,show = true)


#read(stdin, Char)
savefig("Population_Norm_01_1")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("Population_Norm_01_1big")

# plot2 focus
plot!(plt1[1],size=(1000,600),
xaxis="[CE years]",yaxis="[World Population Normalized]",
xticks = -2000:(floor(6000/240)):6000,
#yticks = 0.0:20:120,
xlims = (1700,2300),ylims = (0,180),
show = true)

#read(stdin, Char)
savefig("Population_Norm_01_1focus")



##############################################################
# plot Population Real
# THE POPULATION Real


plt2 = plot(layout=(1,1),title="World Population",size=(600,400))

scatter!(plt2[1],year[1:19],(Popu_Norm[1:19])*Popu[1],color=:black,label="Historic Data")

#############################
# Plot ODE from smooth Splines
# Population

plot!(plt2[1],PopuGDPPCsol1.t,PopuGDPPCsol1[1,:]*Popu[1],linewidth=1,color=:yellow,label="Model 1")
plot!(plt2[1],PopuGDPPCsol1b.t,PopuGDPPCsol1b[1,:]*Popu[1],linewidth=1,color=:yellow,primary=false)

plot!(plt2[1],PopuGDPPCsol2.t,PopuGDPPCsol2[1,:]*Popu[1],linewidth=1,color=:brown,label="Model 2")
plot!(plt2[1],PopuGDPPCsol2b.t,PopuGDPPCsol2b[1,:]*Popu[1],linewidth=1,color=:brown,primary=false)


###################
# Plot ODE from manually estimated parameters simple
# Population

plot!(plt2[1],PopuGDPPCsol3.t,PopuGDPPCsol3[1,:]*Popu[1],linewidth=1,color=:red,label="Model 3")
plot!(plt2[1],PopuGDPPCsol3b.t,PopuGDPPCsol3b[1,:]*Popu[1],linewidth=1,color=:red,primary=false)

plot!(plt2[1],PopuGDPPCsol4.t,PopuGDPPCsol4[1,:]*Popu[1],linewidth=1,color=:orange,label="Model 4")
plot!(plt2[1],PopuGDPPCsol4b.t,PopuGDPPCsol4b[1,:]*Popu[1],linewidth=1,color=:orange,primary=false)

###################
# Plot DDE from manually estimated parameters
# Popu

plot!(plt2[1],PopuGDPPCsol5.t,PopuGDPPCsol5[1,:]*Popu[1],linewidth=1,color=:green,label="Model 5")
plot!(plt2[1],PopuGDPPCsol5b.t,PopuGDPPCsol5b[1,:]*Popu[1],linewidth=1,color=:green,primary=false)

plot!(plt2[1],PopuGDPPCsol6.t,PopuGDPPCsol6[1,:]*Popu[1],linewidth=1,color=:blue,label="Model 6")
plot!(plt2[1],PopuGDPPCsol6b.t,PopuGDPPCsol6b[1,:]*Popu[1],linewidth=1,color=:blue,primary=false,
xaxis="[CE years]",yaxis="[World Population]",xticks = -2000:(floor(6000/10)):6000,
#yticks = 0.0:20:120,
#xlims = (0.5,1.5),ylims = (0.5,1.5),
legend=:bottomright,grid= false,show = true)


#read(stdin, Char)
savefig("Population_Real_01_1")

plot!(plt2[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("Population_Real_01_1big")

# plot focus
plot!(plt2[1],size=(1000,600),
xaxis="[CE years]",yaxis="[World Population]",
xticks = -2000:(floor(6000/240)):6000,
#yticks = 0.0:20:120,
xlims = (1700,2300),
#ylims = (25,120),
show = true)

#read(stdin, Char)
savefig("Population_Real_01_1focus")


######################################################

######################################################
# plot Technical Level
# THE TECH LEVEL

plt3 = plot(layout=(1,1),title="Technical Level",size=(600,400))

scatter!(plt3[1],year[1:19],(GDPPC_Norm[1:19]),color=:black,label="Historic Data")

#############################
# Plot ODE from smooth Splines
# Technical Level

plot!(plt3[1],PopuGDPPCsol1.t,PopuGDPPCsol1[2,:],linewidth=1,color=:yellow,label="Model 1")
plot!(plt3[1],PopuGDPPCsol1b.t,PopuGDPPCsol1b[2,:],linewidth=1,color=:yellow,primary=false)

plot!(plt3[1],PopuGDPPCsol2.t,PopuGDPPCsol2[2,:],linewidth=1,color=:brown,label="Model 2")
plot!(plt3[1],PopuGDPPCsol2b.t,PopuGDPPCsol2b[2,:],linewidth=1,color=:brown,primary=false)


#############################
# Plot ODE from manually estimated parameters simple
# GDPPC

plot!(plt3[1],PopuGDPPCsol3.t,PopuGDPPCsol3[2,:],linewidth=1,color=:red,label="Model 3")
plot!(plt3[1],PopuGDPPCsol3b.t,PopuGDPPCsol3b[2,:],linewidth=1,color=:red,primary=false)

plot!(plt3[1],PopuGDPPCsol4.t,PopuGDPPCsol4[2,:],linewidth=1,color=:orange,label="Model 4")
plot!(plt3[1],PopuGDPPCsol4b.t,PopuGDPPCsol4b[2,:],linewidth=1,color=:orange,primary=false)


#############################
# Plot DDE from manually estimated parameters
# GDPPC

plot!(plt3[1],PopuGDPPCsol5.t,PopuGDPPCsol5[2,:],linewidth=1,color=:green,label="Model 5")
plot!(plt3[1],PopuGDPPCsol5b.t,PopuGDPPCsol5b[2,:],linewidth=1,color=:green,primary=false)

plot!(plt3[1],PopuGDPPCsol6.t,PopuGDPPCsol6[2,:],linewidth=1,color=:blue,label="Model 6")
plot!(plt3[1],PopuGDPPCsol6b.t,PopuGDPPCsol6b[2,:],linewidth=1,color=:blue,primary=false,
xaxis="[CE years]",yaxis="[TechLevel]",xticks = -2000:(floor(6000/10)):6000,
#yticks = 0.0:20:120,
#xlims = (0.5,1.5),ylims = (0.5,1.5),
legend=:bottomright,grid= false,show = true)


#read(stdin, Char)
savefig("TechLevel_01_1")

plot!(plt3[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("TechLevel_01_1big")

# plot focus
plot!(plt3[1],size=(1000,600),
xaxis="[CE years]",yaxis="[TechLevel]",
xticks = -2000:(floor(6000/180)):6000,
#yticks = 0.0:20:120,
xlims = (1800,2800),ylims = (0,270),
show = true)

#read(stdin, Char)
savefig("TechLevel_01_1focus")



#################################################################
#################################################################
# Plot LEB , HLE, ETA data
# THE LEB

plt4 = plot(layout=(1,1),title="LEB, HLE, ETA",size=(600,400))

# Scatter input data
scatter!(plt4[1],year[1:19],(LEB[1:19]),label="Historic LEB")

# model5 and 6 LEB data
Leb5=@. LEBNorm(PopuGDPPCsol5[2,:])*LEB[1]
Leb5b=@. LEBNorm(PopuGDPPCsol5b[2,:])*LEB[1]
Leb6=@. LEBNorm(PopuGDPPCsol6[2,:])*LEB[1]
Leb6b=@. LEBNorm(PopuGDPPCsol6b[2,:])*LEB[1]

Hle5=@. HLE(Leb5)
Hle5b=@. HLE(Leb5b)
Hle6=@. HLE(Leb6)
Hle6b=@. HLE(Leb6b)

Eta5=@. ETA(Leb5)
Eta5b=@. ETA(Leb5b)
Eta6=@. ETA(Leb6)
Eta6b=@. ETA(Leb6b)

plot!(plt4[1],PopuGDPPCsol5.t,Leb5,linewidth=1,color=:red, label="LEB model5")
plot!(plt4[1],PopuGDPPCsol5b.t,Leb5b,linewidth=1,color=:red,primary=false)
plot!(plt4[1],PopuGDPPCsol5.t,Hle5,linewidth=1,color=:orange, label="HLE model5")
plot!(plt4[1],PopuGDPPCsol5b.t,Hle5b,linewidth=1,color=:orange,primary=false)
plot!(plt4[1],PopuGDPPCsol5.t,Eta5,linewidth=1,color=:pink, label="ETA model5")
plot!(plt4[1],PopuGDPPCsol5b.t,Eta5b,linewidth=1,color=:pink,primary=false)

plot!(plt4[1],PopuGDPPCsol6.t,Leb6,linewidth=1,color=:green,label="LEB model6")
plot!(plt4[1],PopuGDPPCsol6b.t,Leb6b,linewidth=1,color=:green,primary=false)
plot!(plt4[1],PopuGDPPCsol6.t,Hle6,linewidth=1,color=:blue,label="HLE model6")
plot!(plt4[1],PopuGDPPCsol6b.t,Hle6b,linewidth=1,color=:blue,primary=false)
plot!(plt4[1],PopuGDPPCsol6.t,Eta6,linewidth=1,color=:cyan,label="ETA model6")
plot!(plt4[1],PopuGDPPCsol6b.t,Eta6b,linewidth=1,color=:cyan,primary=false,
xaxis="[CE years]",yaxis="LEB[years]",xticks = -2000:(floor(6000/10)):6000,yticks = 0.0:10:120,
#xlims = (0.5,1.5),ylims = (0.5,1.5),
ylims = (0,130),
legend=:topleft,grid= false,show = true)


#read(stdin, Char)
savefig("LEB_HLE_ETA_01_1")

plot!(plt4[1],title="LEB(life expectancy at birth), HLE(healthy life expectancy), ETA(end of training age)",size=(1000,600),show = true)

#read(stdin, Char)
savefig("LEB_HLE_ETA_01_1big")

# plot focus
plot!(plt4[1],size=(1000,600),
xaxis="[CE years]",yaxis="LEB[years]",
xticks = -2000:(floor(6000/240)):6000,yticks = 0.0:10:120,
xlims = (1700,2300),ylims = (0,130),
show = true)

#read(stdin, Char)
savefig("LEB_HLE_ETA_01_1focus")

#######################################################
#######################################################
# Plot of Technarian Age Jump Functional Response
# THE LINE
#

plt5 = plot(layout=(1,1),title="In Caput Evolution Functional Response",size=(600,400))

plot!(plt5[1],(PopuGDPPCsol5[2,:]),(PopuGDPPCsol5[1,:]),linewidth=2,color=:orange,label="Model 5 Future" )
plot!(plt5[1],(PopuGDPPCsol5b[2,:]),(PopuGDPPCsol5b[1,:]),linewidth=5,color=:red,label="Model 5 Past" )

plot!(plt5[1],(PopuGDPPCsol6[2,:]),(PopuGDPPCsol6[1,:]),linewidth=2,color=:cyan,label="Model 6 Future")
plot!(plt5[1],(PopuGDPPCsol6b[2,:]),(PopuGDPPCsol6b[1,:]),linewidth=4,color=:blue,label="Model 6 Past" ,
xaxis="Technical Level",yaxis="Population_Normalized",
#xticks = -2000:(floor(6000/20)):6000,yticks = 0.0:20:120,xlims = (0.5,1.5),ylims = (0.5,1.5),
legend=:bottomright,show = true)


#read(stdin, Char)
savefig("FunctionalResponse_01_1")

plot!(plt5[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("FunctionalResponse_01_1big")




#########################################
# Plot of Derivatives of Technarian Age
# Plot derivatives from Dierckxs interpolations
# and derivatives from identified (model 1 and2) models via LsqFit
# THE INTERPOLATIONS


sx = collect(year[9]:5.0:year[19])
sxx = collect(year[9]:5.0:year[19])

plt6 = plot(layout=(1,1),title="Population Derivative and LsqFit",size=(600,400))

plot!(plt6[1],sx,(derivative(spl_Technarian_p, sx; nu=1)),linewidth=1,color=:red, label="Historic Population Spline")
plot!(plt6[1],sxx,(derivative(spl_Pre_Technarian_p, sxx; nu=1)),linewidth=1,color=:orange,label="Historic Population Straight")

plot!(plt6[1],(sd),(multimodelp1(xdata, paramp1)),linewidth=1,color=:green, label="Population LsqFit Exponential (model1)")
plot!(plt6[1],(sd),(multimodelp2(xdata, paramp2)),linewidth=1,color=:blue, label="Population LsqFit Stirling (model2)",
xaxis="[CE years]",yaxis="[Derivative of Population]",xticks = year[9]:(floor(year[19]/40)):year[19],
# yticks = 0.0:20:120, xlims = (0.5,1.5),ylims = (0.5,1.5),
legend=:topleft,grid= false,show = true)


#read(stdin, Char)
savefig("PopulationDerivatives_01_1")

plot!(plt6[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("PopulationDerivatives_01_1big")

#######

plt7 = plot(layout=(1,1),title="Technical Level Derivative and LsqFit",size=(600,400))

plot!(plt7[1],sx,(derivative(spl_Technarian_g, sx; nu=1)),linewidth=1,color=:red, label="Historic Tech Level Spline")
plot!(plt7[1],sxx,(derivative(spl_Pre_Technarian_g, sxx; nu=1)),linewidth=1,color=:orange, label="Historic Tech Level Straight")

plot!(plt7[1],(sd),(multimodelg1(xdata, paramg1)),linewidth=1,color=:green, label="Tech Level LsqFit Exponential (model1)")
plot!(plt7[1],(sd),(multimodelg2(xdata, paramg2)),linewidth=1,color=:blue, label="Tech Level LsqFit Stirling (model2)",
xaxis="[CE years]",yaxis="[Derivative of Tech Level]",xticks = year[9]:(floor(year[19]/40)):year[19],
# yticks = 0.0:20:120, xlims = (0.5,1.5),ylims = (0.5,1.5),
legend=:topleft,grid= false,show = true)

#read(stdin, Char)
savefig("TLDerivatives_01_1")

plot!(plt7[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("TLDerivatives_01_1big")


#############################################################################################
#############################################################################################
# World Population % Growth Rates (from 1950 to 2099) 
# source
# https://ourworldindata.org/grapher/population-growth-rates?tab=chart&time=earliest..latest


RealUNpopuGRyears=
[1950.0, 1951.0, 1952.0, 1953.0, 1954.0, 1955.0, 1956.0, 1957.0, 1958.0, 1959.0, 1960.0, 1961.0, 1962.0, 1963.0, 1964.0, 1965.0, 1966.0, 1967.0, 1968.0, 1969.0, 1970.0, 1971.0, 1972.0, 1973.0, 1974.0, 1975.0, 1976.0, 1977.0, 1978.0, 1979.0, 1980.0, 1981.0, 1982.0, 1983.0, 1984.0, 1985.0, 1986.0, 1987.0, 1988.0, 1989.0, 1990.0, 1991.0, 1992.0, 1993.0, 1994.0, 1995.0, 1996.0, 1997.0, 1998.0, 1999.0, 2000.0, 2001.0, 2002.0, 2003.0, 2004.0, 2005.0, 2006.0, 2007.0, 2008.0, 2009.0, 2010.0, 2011.0, 2012.0, 2013.0, 2014.0, 2015.0, 2016.0, 2017.0, 2018.0, 2019.0, 2020.0]


RealUNpopuGrowthRate=
[1.859, 1.828,	1.779,	1.755,	1.751,	1.76,	1.776,	1.796,	1.816,	1.833,	1.849,	1.867,	1.891,	1.923,	1.963,	2.005,	2.042,	2.065,	2.068,	2.053,	2.028,	2.003,	1.974,	1.937,	1.895,	1.849,	1.806,	1.772,	1.754,	1.752,	1.755,	1.756,	1.76,	1.773,	1.791,	1.812,	1.83,	1.829,	1.797,	1.736,	1.66,	1.586,	1.522,	1.471,	1.436,	1.408,	1.381,	1.354,	1.33,	1.309,	1.289,	1.272,	1.258,	1.249,	1.244,	1.241,	1.238,	1.235,	1.229,	1.22,	1.211,	1.2,	1.189,	1.175,	1.16,	1.143,	1.126,	1.107,	1.085,	1.061,	1.036]


EstimatesUNpopuGRyears=
[2020.0, 2021.0, 2022.0, 2023.0, 2024.0, 2025.0, 2026.0, 2027.0, 2028.0, 2029.0, 2030.0, 2031.0, 2032.0, 2033.0, 2034.0, 2035.0, 2036.0, 2037.0, 2038.0, 2039.0, 2040.0, 2041.0, 2042.0, 2043.0, 2044.0, 2045.0, 2046.0, 2047.0, 2048.0, 2049.0, 2050.0, 2051.0, 2052.0, 2053.0, 2054.0, 2055.0, 2056.0, 2057.0, 2058.0, 2059.0, 2060.0, 2061.0, 2062.0, 2063.0, 2064.0, 2065.0, 2066.0, 2067.0, 2068.0, 2069.0, 2070.0, 2071.0, 2072.0, 2073.0, 2074.0, 2075.0, 2076.0, 2077.0, 2078.0, 2079.0, 2080.0, 2081.0, 2082.0, 2083.0, 2084.0, 2085.0, 2086.0, 2087.0, 2088.0, 2089.0, 2090.0, 2091.0, 2092.0, 2093.0, 2094.0, 2095.0, 2096.0, 2097.0, 2098.0, 2099.0]


EstimatesUNpopuGrowthRate=
[1.036, 1.011,	0.986,	0.963,	0.941,	0.921,	0.9,	0.88,	0.86,	0.841,	0.823,	0.805,	0.787,	0.769,	0.751,	0.732,	0.714,	0.697,	0.68,	0.663,	0.646,	0.63,	0.614,	0.598,	0.582,	0.566,	0.55,	0.535,	0.519,	0.504,	0.489,	0.475,	0.46,	0.446,	0.432,	0.418,	0.404,	0.391,	0.378,	0.366,	0.354,	0.342,	0.331,	0.319,	0.308,	0.298,	0.287,	0.277,	0.267,	0.257,	0.247,	0.238,	0.229,	0.22,	0.211,	0.202,	0.194,	0.186,	0.178,	0.17,	0.162,	0.154,	0.147,	0.139,	0.133,	0.126,	0.12,	0.113,	0.107,	0.1,	0.094,	0.087,	0.081,	0.074,	0.068,	0.061,	0.054,	0.047,	0.04,	0.032]

#######################################################################################################
# World Population % GDP and GDPPC Growth Rates (from 1961 to 2019) 
# source
# https://data.worldbank.org/indicator/NY.GDP.PCAP.KD.ZG

WorldBankGDPyears=
[1961.0, 1962.0, 1963.0, 1964.0, 1965.0, 1966.0, 1967.0, 1968.0, 1969.0, 1970.0, 1971.0, 1972.0, 1973.0, 1974.0, 1975.0, 1976.0, 1977.0, 1978.0, 1979.0, 1980.0, 1981.0, 1982.0, 1983.0, 1984.0, 1985.0, 1986.0, 1987.0, 1988.0, 1989.0, 1990.0, 1991.0, 1992.0, 1993.0, 1994.0, 1995.0, 1996.0, 1997.0, 1998.0, 1999.0, 2000.0, 2001.0, 2002.0, 2003.0, 2004.0, 2005.0, 2006.0, 2007.0, 2008.0, 2009.0, 2010.0, 2011.0, 2012.0, 2013.0, 2014.0, 2015.0, 2016.0, 2017.0, 2018.0, 2019.0]


WorldBankGDPrate=
[4.2991830819569, 5.55413660359476, 5.35067756620416, 6.71355681792285, 5.5196444061129, 5.76849813773291, 4.48595098085787, 6.31352836560291, 6.11362814312773, 3.8598737570305, 4.34152325743152, 5.72382344767152, 6.50525271058334, 1.99594553474572, 0.603101336836744, 5.26935520826883, 3.93167714024989, 3.89217831329481, 4.12722123711964, 1.9030934393484, 1.92751400751163, 0.423811737894013, 2.41106597373202, 4.50486538530438, 3.71063583582453, 3.39450985129828, 3.70789618012796, 4.61908187329017, 3.67466134322667, 2.90949442488176, 1.41886837530257, 1.76386888858855, 1.53112826799003, 3.00120765116208, 3.03658854203159, 3.38452035675527, 3.46610155225279, 2.55732309815349, 3.24920337802101, 4.38807708017994, 1.96078124801235, 2.18237967746562, 2.96593220546156, 4.40840084348217, 3.91611492851847, 4.37819663780303, 4.32213629244085, 1.85232660162453, -1.67364137627817, 4.30301651760969, 3.13777412881937, 2.51850732689265, 2.66597893243782, 2.86109839141314, 2.87394930113432, 2.60578989614638, 3.29862797064233, 2.97677626566207, 2.34337753648894]


WorldBankGDPPCrate=
[2.90591241168217, 3.76501360417755, 3.20086023820465, 4.56684068975464, 3.39497989431041, 3.58460738072065, 2.39081269049237, 4.19597018647401, 3.91712796385897, 1.73323522850197, 2.1884299942819, 3.61909817538975, 4.45114203288671, 0.049072154868242, -1.23883237147832, 3.42245682947211, 2.14465536658395, 2.10746596225211, 2.3243356186575, 0.151880539658308, 0.16033657749901, -1.35194742634485, 0.616109826392702, 2.71073406526348, 1.92832180565061, 1.59600546252598, 1.89093030165705, 2.79909505928784, 1.90278206801136, 1.15381135883254, -0.243543196000459, 0.190783817588056, -0.031227064266886, 1.45714001322779, 1.50455676588375, 1.90357450182471, 2.01019031476275, 1.14971790624895, 1.87210295250017, 3.0249800679717, 0.654691714662164, 0.894156589481781, 1.68322075784278, 3.1150673615065, 2.63631553952453, 3.0959985646374, 3.04872122475383, 0.603977265755205, -2.86043758288393, 3.06290993826943, 1.94498353156864, 1.31907413857594, 1.46492184309635, 1.66157813382308, 1.68601335354921, 1.42647550355652, 2.13215782453886, 1.85204659911035, 1.25468001220281]


#####################################################################################################

######################################################
# Plot Derivative Loop of Functional Response
# THE LOOP
#

plt8 = plot(layout=(1,1),title="Derivative Loop of Functional Response",size=(600,400))

plot!(plt8[1],((derivative(spl_Technarian_g, sx; nu=1)),(derivative(spl_Technarian_p, sx; nu=1))),linewidth=1,
color=:red, label="Historic Spline")
plot!(plt8[1],((derivative(spl_Pre_Technarian_g, sxx; nu=1)),(derivative(spl_Pre_Technarian_p, sxx; nu=1))),linewidth=1,
color=:orange, label="Historic Straight")

plot!(plt8[1],(dG3dt,dP3dt),linewidth=1,color=:lightblue, label="Model 3")
plot!(plt8[1],(dG3bdt,dP3bdt),linewidth=3,color=:lightblue,primary=false)
plot!(plt8[1],(dG4dt,dP4dt),linewidth=1,color=:lightgreen, label="Model 4")
plot!(plt8[1],(dG4bdt,dP4bdt),linewidth=3,color=:lightgreen,primary=false)
plot!(plt8[1],(dG5dt,dP5dt),linewidth=1,color=:blue, label="Model 5")
plot!(plt8[1],(dG5bdt,dP5bdt),linewidth=4,color=:blue,primary=false)
plot!(plt8[1],(dG6dt,dP6dt),linewidth=1,color=:green, label="Model 6")
plot!(plt8[1],(dG6bdt,dP6bdt),linewidth=3,color=:green,primary=false,
xaxis="Technical Level Derivative",yaxis="Population_Normalized Derivative",
#xticks = -2000:(floor(6000/20)):6000,yticks = 0.0:20:120,xlims = (0.5,1.5),ylims = (0.5,1.5),
legend=:bottomright,show = true)

#read(stdin, Char)
savefig("DerivativesLoop_01_1")

plot!(plt8[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("DerivativesLoop_01_1big")

######################################################
# Plot %Growth Loop of Functional Response
# THE LOOP
#

plt9 = plot(layout=(1,1),title="%Growth Loop of Functional Response",size=(600,400))

plot!(plt9[1],(GrowthTL5,GrowthP5),linewidth=1,color=:blue, label="Model 5")
plot!(plt9[1],(GrowthTL5b,GrowthP5b),linewidth=3,color=:blue,primary=false)
plot!(plt9[1],(GrowthTL6,GrowthP6),linewidth=1,color=:green, label="Model 6")
plot!(plt9[1],(GrowthTL6b,GrowthP6b),linewidth=3,color=:green,primary=false,
xaxis="Technical Level % Growth",yaxis="Population % Growth",
#xticks = -2000:(floor(6000/20)):6000,yticks = 0.0:20:120,xlims = (0.5,1.5),ylims = (0.5,1.5),
legend=:bottomright,show = true)

#read(stdin, Char)
savefig("GrowthLoop_01_1")

plot!(plt9[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("GrowthLoop_01_1big")


############################
# Plot Derivative Runs
# THE SWITCH
#
#

plt10 = plot(layout=(1,1),title="The Technarian Age Switch",size=(600,400))

plot!(plt10[1],(sx),(derivative(spl_Technarian_p, sx; nu=1)),linewidth=1,color=:red, label="Population Spline")
plot!(plt10[1],sx,(derivative(spl_Technarian_g, sx; nu=1)),linewidth=1,color=:orange, label="Technical Level Spline")

plot!(plt10[1],saved_values5.t[2:end],(dP5dt),linewidth=1,color=:blue, label="Model5 Population Derivative")
plot!(plt10[1],saved_values5b.t[2:end],(dP5bdt),linewidth=4,color=:blue,primary=false)
plot!(plt10[1],saved_values6.t[2:end],(dP6dt),linewidth=1,color=:green, label="Model6 Population Derivative")
plot!(plt10[1],saved_values6b.t[2:end],(dP6bdt),linewidth=3,color=:green,primary=false)

plot!(plt10[1],saved_values5.t[2:end],(dG5dt),linewidth=1,color=:lightblue, label="Model5 Tech Level Derivative")
plot!(plt10[1],saved_values5b.t[2:end],(dG5bdt),linewidth=4,color=:lightblue,primary=false)
plot!(plt10[1],saved_values6.t[2:end],(dG6dt),linewidth=1,color=:lightgreen, label="Model6 Tech Level Derivative")
plot!(plt10[1],saved_values6b.t[2:end],(dG6bdt),linewidth=3,color=:lightgreen,primary=false,

xaxis="[CE years]",yaxis="[derivatives]",
xticks = -2000:(floor(6000/10)):6000,
# yticks = 0.0:20:120, xlims = (0.5,1.5),ylims = (0.5,1.5),
grid= false,show = true)

#read(stdin, Char)
savefig("TheSwitch_01_1")

plot!(plt10[1],size=(1000,600),
xticks = -2000:(floor(6000/10)):6000,
show = true)

#read(stdin, Char)
savefig("TheSwitch_01_1big")

# plot focus
plot!(plt10[1],size=(1000,600),
xticks = -2000:(floor(6000/240)):6000, yticks = 0.0:0.05:2,
xlims = (1700,2300),ylims = (0,1.2),
grid= true,show = true)

#read(stdin, Char)
savefig("TheSwitch_01_1focus")

# plot superfocus
plot!(plt10[1],size=(1000,600),
xticks = -2000:(floor(6000/840)):6000, yticks = 0.0:0.05:2,
xlims = (1950,2100),ylims = (0,1.2),
grid= true,show = true)

#read(stdin, Char)
savefig("TheSwitch_01_1superfocus")

############################
# Plot %Growth Runs
# THE SWITCH 
#
# Population and GDPPC=TL % Growth

plt11 = plot(layout=(1,1),title="The Technarian Age Switch %Growth",size=(600,400))


plot!(plt11[1],saved_values5.t[1:end],GrowthP5,linewidth=1,color=:blue, label="Model5 Population %Growth")
plot!(plt11[1],saved_values5b.t[1:end],GrowthP5b,linewidth=4,color=:blue,primary=false)
plot!(plt11[1],saved_values6.t[1:end],GrowthP6,linewidth=1,color=:green, label="Model6 Population %Growth")
plot!(plt11[1],saved_values6b.t[1:end],GrowthP6b,linewidth=3,color=:green,primary=false)
plot!(plt11[1],saved_values5.t[1:end],GrowthTL5,linewidth=1,color=:red, label="Model5 Tech Level %Growth")
plot!(plt11[1],saved_values5b.t[1:end],GrowthTL5b,linewidth=4,color=:red,primary=false)
plot!(plt11[1],saved_values6.t[1:end],GrowthTL6,linewidth=1,color=:violet, label="Model6 Tech Level %Growth")
plot!(plt11[1],saved_values6b.t[1:end],GrowthTL6b,linewidth=3,color=:violet,primary=false,

xaxis="[CE years]",yaxis="[% Growth]",
xticks = -2000:(floor(6000/10)):6000,
# yticks = 0.0:20:120, xlims = (0.5,1.5),ylims = (0.5,1.5),
grid= false,show = true)

#read(stdin, Char)
savefig("TheSwitchGrowth_01_1")

plot!(plt11[1],size=(1000,600),
xticks = -2000:(floor(6000/20)):6000,
show = true)

#read(stdin, Char)
savefig("TheSwitchGrowth_01_1big")

# plot focus
plot!(plt11[1],size=(1000,600),
xticks = -2000:(floor(6000/240)):6000, yticks = 0.0:0.1:5,
xlims = (1700,2300),ylims = (0,2.0),
grid= true,show = true)

#read(stdin, Char)
savefig("TheSwitchGrowth_01_1focus")

# plot superfocus 

plot!(plt11[1],RealUNpopuGRyears,RealUNpopuGrowthRate,linewidth=3,color=:black, label="RealUNpopuGrowthRate")
plot!(plt11[1],EstimatesUNpopuGRyears,EstimatesUNpopuGrowthRate,linewidth=1,color=:yellow,label="EstimatesUNpopuGrowthRate")

plot!(plt11[1],WorldBankGDPyears,WorldBankGDPPCrate,linewidth=1,color=:brown, label="WorldBankGDPPCrate")


plot!(plt11[1],size=(1000,600),
xticks = -2000:(floor(6000/840)):6000, yticks = 0.0:0.2:10,
xlims = (1940,2110),ylims = (0,4.0),
grid= true,show = true)

#read(stdin, Char)
savefig("TheSwitchGrowth_01_1superfocus")


########################################
# Plot GDP%Growth Runs
# THE SWITCH 
#
# GDP % Growth
#

plt12 = plot(layout=(1,1),title="The Technarian Age Switch GDP %Growth (AM 1990Intld)",size=(600,400))

plot!(plt12[1],saved_values5.t[1:end],GrowthGDP5,linewidth=1,color=:orange, label="Model5 Tech Level GDP %Growth")
plot!(plt12[1],saved_values5b.t[1:end],GrowthGDP5b,linewidth=4,color=:orange,primary=false)
plot!(plt12[1],saved_values6.t[1:end],GrowthGDP6,linewidth=1,color=:red, label="Model6 Tech Level GDP %Growth")
plot!(plt12[1],saved_values6b.t[1:end],GrowthGDP6b,linewidth=3,color=:red,primary=false,

xaxis="[CE years]",yaxis="[GDP % Growth]",
xticks = -2000:(floor(6000/10)):6000,
# yticks = 0.0:20:120, xlims = (0.5,1.5),ylims = (0.5,1.5),
grid= false,show = true)

#read(stdin, Char)
savefig("TheSwitchGDPGrowth_01_1")

plot!(plt12[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("TheSwitchGDPGrowth_01_1big")

# plot focus
plot!(plt12[1],size=(1000,600),
xticks = -2000:(floor(6000/240)):6000,yticks = 0.0:0.10:5,
xlims = (1700,2300),ylims = (0,3.8),
grid= true,show = true)

#read(stdin, Char)
savefig("TheSwitchGDPGrowth_01_1focus")

# plot superfocus
plot!(plt12[1],WorldBankGDPyears,WorldBankGDPrate,linewidth=1,color=:black, label="WorldBankGDPrate")

plot!(plt12[1],size=(1000,600),
xticks = -2000:(floor(6000/840)):6000, yticks = -3.0:0.25:10,
xlims = (1950,2100),ylims = (-2.0,7.0),
grid= true,show = true)

#read(stdin, Char)
savefig("TheSwitchGDPGrowth_01_1superfocus")

# plot superfocus2
plot!(plt12[1],WorldBankGDPyears,WorldBankGDPrate,linewidth=1,color=:black, label="WorldBankGDPrate")

plot!(plt12[1],size=(1000,600),
xticks = -2000:(floor(6000/840)):6000, yticks = 0.0:0.1:10,
xlims = (1950,2100),ylims = (0.0,3.5),
grid= true,show = true)

#read(stdin, Char)
savefig("TheSwitchGDPGrowth_01_2superfocus")


########################################
# Plot Population%Growth Runs
# THE SWITCH 
#
# Population  % Growth

plt13 = plot(layout=(1,1),title="The Technarian Age Switch Population %Growth",size=(600,400))

plot!(plt13[1],saved_values5.t[1:end],GrowthP5,linewidth=1,color=:blue, label="Model5 Population %Growth")
plot!(plt13[1],saved_values5b.t[1:end],GrowthP5b,linewidth=4,color=:blue,primary=false)
plot!(plt13[1],saved_values6.t[1:end],GrowthP6,linewidth=1,color=:green, label="Model6 Population %Growth")
plot!(plt13[1],saved_values6b.t[1:end],GrowthP6b,linewidth=3,color=:green,primary=false,

xaxis="[CE years]",yaxis="[Population % Growth]",
xticks = -2000:(floor(6000/10)):6000,
# yticks = 0.0:20:120, xlims = (0.5,1.5),ylims = (0.5,1.5),
grid= false,show = true)

#read(stdin, Char)
savefig("TheSwitchPopulationGrowth_01_1")

plot!(plt13[1],size=(1000,600),
xticks = -2000:(floor(6000/20)):6000,
show = true)

#read(stdin, Char)
savefig("TheSwitchPopulationGrowth_01_1big")

# plot focus
plot!(plt13[1],size=(1000,600),
xticks = -2000:(floor(6000/240)):6000, yticks = 0.0:0.1:5,
xlims = (1700,2300),ylims = (0,1.7),
grid= true,show = true)

#read(stdin, Char)
savefig("TheSwitchPopulationGrowth_01_1focus")

# plot superfocus 

plot!(plt13[1],RealUNpopuGRyears,RealUNpopuGrowthRate,linewidth=3,color=:black, label="RealUNpopuGrowthRate")
plot!(plt13[1],EstimatesUNpopuGRyears,EstimatesUNpopuGrowthRate,linewidth=1,color=:orange,label="EstimatesUNpopuGrowthRate")


plot!(plt13[1],size=(1000,600),
xticks = -2000:(floor(6000/840)):6000, yticks = 0.0:0.1:5.0,
xlims = (1940,2110),ylims = (0,2.2),
grid= true,show = true)

#read(stdin, Char)
savefig("TheSwitchPopulationGrowth_01_1superfocus")

########################################
# Plot GDPPC%Growth Runs
# THE SWITCH 
#  Tech LEVEL
#
# GDPPC % Growth  = Tech Level
#

plt14 = plot(layout=(1,1),title="The Technarian Age Switch TECH LEVEL %Growth",size=(600,400))

plot!(plt14[1],saved_values5.t[1:end],GrowthTL5,linewidth=1,color=:orange, label="Model5 Tech Level %Growth")
plot!(plt14[1],saved_values5b.t[1:end],GrowthTL5b,linewidth=4,color=:orange,primary=false)
plot!(plt14[1],saved_values6.t[1:end],GrowthTL6,linewidth=1,color=:violet, label="Model6 Tech Level %Growth")
plot!(plt14[1],saved_values6b.t[1:end],GrowthTL6b,linewidth=3,color=:violet,primary=false,

xaxis="[CE years]",yaxis="[TECH LEVEL % Growth]",
xticks = -2000:(floor(6000/10)):6000,
# yticks = 0.0:20:120, xlims = (0.5,1.5),ylims = (0.5,1.5),
grid= false,show = true)

#read(stdin, Char)
savefig("TheSwitchTECHLEVELGrowth_01_1")

plot!(plt14[1],size=(1000,600),
xticks = -2000:(floor(6000/20)):6000,
show = true)

#read(stdin, Char)
savefig("TheSwitchTECHLEVELGrowth_01_1big")

# plot focus
plot!(plt14[1],size=(1000,600),
xticks = -2000:(floor(6000/240)):6000, yticks = 0.0:0.1:5,
xlims = (1700,2300),ylims = (0,2.0),
grid= true,show = true)

#read(stdin, Char)
savefig("TheSwitchTECHLEVELGrowth_01_1focus")

# plot superfocus 


plot!(plt14[1],WorldBankGDPyears,WorldBankGDPPCrate,linewidth=1,color=:brown, label="WorldBankGDPPCrate")


plot!(plt14[1],size=(1000,600),
xticks = -2000:(floor(6000/840)):6000, yticks = 0.0:0.2:10,
xlims = (1940,2110),ylims = (0,4.0),
grid= true,show = true)

#read(stdin, Char)
savefig("TheSwitchTECHLEVELGrowth_01_1superfocus")


#### end ####

