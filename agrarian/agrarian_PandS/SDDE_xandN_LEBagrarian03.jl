#=
/*
03 Version.  Conditions squewed toward degenerated Elite State and over reproduction.
Version 03 has bb=0.33, and gets all of the sudden close to instability.
*/
/* Stochastic DDE LEB Agrarian Bioenergy (x) and Population (N) version   (N= total Population)*/
N = S+P, State (S) and People (N).  S = keju*N, keju as state civil servant selection exams 
/* build 1.0 */
=#
""" 
/*
01 Version.  Nominal Conditions
Code run on Julia (1.5.3, julialang.org).
*/

/*
Conditions to calculate parameters:

LEBagrarianParameters_PES_01.ods :
Spreadsheet used to identify Nominal Conditions Parameters for LEBVenatorian

aa= 0.25 by choice
(Nominal)bb= 0.24 by choice

Strong Input/Output Targets:
at NNo=1 with PNo= 0.9 and ENo=0.1 then dPNot=0 and dENot=0

All Other parameters set with common sense, trying to match a Produce/People with Defenders/Exploiters/Elite 
early agrarian society.  

cshs  : 0.5; copulation system of homo sapiens
(Nominal)Freducr : 0.9 ; fertility reduction .  0.9 instead 1.0 because with advanced agrarian society some artificial fertility reduction can be expected.
nobhs: 1.0 ; number of offsping per birth for humans. set to 1.0.
/*Freducp : 1.0 ; fertility reduction for P (people) .  Considered 0 reduction at paleo conditions
/*Freduce : 0.6 ; fertility reduction for E (Elite) .  Considered 0 reduction at paleo conditions
AoMhs : 18.0 ; Age of menarche for humans
IBIhs  :  3.0 ;  Inter birth interval for humans, includes correction for 50% miscarriges and death 
(Includes pregnancy(9 month) +LAM (14 months)+ refecundization probabylity time (6 months, 4 moths for modern humans)). 
Srp : xxr (RRo,NNo)/(xxr (RRo,NNo)+1); Perinatal and juvenile survival rate for P (People)
Sre : xxn (RRo,NNo)/(xxn (RRo,NNo)+1); Perinatal and juvenile survival rate for E (Elite)

(Nominal)HHa : 1.0  ; habitat container size
(Nominal)Klhs : 30 ; constant linking x to LEB for humans.
(Nominal)KtoX : 2.47; constant linking E and P numbers to the x output.
(Nominal)keju : 0.1 ;   keju*N = SN , people that works for the state.
(Nominal)fleecing : 0.00001 ; constant to account for a continous loss of bioenergy.
(Nominal)xx0 : 0.20 ; minimal value of bioenergy at which bioenergy output goes to 0.

*/
"""
#=
Note: E^0.25*P to account for economy of scale factor
Social mobility, cannibalism ignored
=#

#using Dierckx
#using LsqFit
using DifferentialEquations
using StochasticDelayDiffEq
#using DataFrames
using CSV
using Plots

a=2.0 	# axis max */
b=1.3 	# axis max */
dur=500  # durat of t */

#= Nominal Version Parameters
nobhs =1.0 
cshs  =0.5
#Freducp = 1.0 
#Freduce = 0.5
Freducn = 0.9
IBIhs  = 3.0 
AoMhs =18.0  

aa=0.25
bb=0.24
HHa=1.0
#EOr=5.0
KLhs=30.0
#Kev=150.0
KtoX=2.47
#KXtoS=10.98
keju = 0.1

=#

# Close to implosion Version Parameters

nobhs =1.0 
cshs  =0.5
Freducn = 0.95
IBIhs  = 3.0 
AoMhs =18.0  

aa=0.25
bb=0.33
HHa=1.08
KLhs=31.0
KtoX=3.0

keju = 0.22
fleecing = 0.0915
xx0= 0.25

#--------------------------------------------------
# x, N isoclines and quiver plot preparation

XXo(xxn,NNo)= ((xxn-xx0)^(1+bb))*((KtoX*((keju*NNo)^bb)*((1-keju)*NNo)*(HHa/(xxn*NNo))))*(1/(((xxn/1.8)^20)+1))

Srn(xxn,NNo)=xxn/(xxn+1)

LEBn(xxn,NNo) =(KLhs*(xxn)^aa)

dxxnt(xxn,NNo) = 0.5*((XXo(xxn,NNo)/(NNo))-xxn)-fleecing
#fleecing of bioenergy

dNNot(xxn,NNo)= (NNo*Freducn*cshs*nobhs*Srn(xxn,NNo)*(LEBn(xxn,NNo)-AoMhs)/IBIhs/LEBn(xxn,NNo) - NNo/LEBn(xxn,NNo))

#####################################

# plot Quiver and Isoclines

# quiver chart
meshgrid(xP, yS) = (repeat(xP, outer=length(yS)), repeat(yS, inner=length(xP)))
xP, yS = meshgrid(0.3:0.1:a, 0.1:0.1:b)
uu = @.dxxnt(xP,yS)
vv = @.dNNot(xP,yS)

# isoclines color ###############
plt1 = plot(layout=(1,1),show = true,title="keju_xxn_NNo_LEBagrarian, SN=0.22",size=(600,400))
quiver!(xP,yS,quiver=(0.03*uu,0.03*vv), color=:lightgray)
contour!(0.3:0.001:a, 0.1:0.001:b, dxxnt, levels=[0], color=:cyan, colorbar = false, w=3)
contour!(0.3:0.001:a, 0.1:0.001:b, dNNot, levels=[0], color=:orange,colorbar = false, w=3)

plot!(plt1[1],
xlims = (0.3,a),ylims = (0.2,(b-0.01)),xaxis="[xxn]",yaxis="[N/No]",xticks = 0:0.1:a,yticks = 0:0.1:(b-0.01),grid= false, show = true)

#read(stdin, Char)
savefig("SDDE_xandN_LEBagrarian03_1")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_xandN_LEBagrarian03_1big")

# isoclines BW ###############
plt1 = plot(layout=(1,1),show = true,title="keju_xxn_NNo_LEBagrarian, SN=0.22",size=(600,400))
quiver!(xP,yS,quiver=(0.03*uu,0.03*vv), color=:black)
contour!(0.3:0.001:a, 0.1:0.001:b, dxxnt, levels=[0], color=:black, colorbar = false, w=4)
contour!(0.3:0.001:a, 0.1:0.001:b, dNNot, levels=[0], color=:black, colorbar = false, w=3)

plot!(plt1[1],
xlims = (0.3,a),ylims = (0.2,(b-0.01)),xaxis="[xxn]",yaxis="[N/No]",xticks = 0:0.1:a,yticks = 0:0.1:(b-0.01),grid= false, show = true)

#read(stdin, Char)
savefig("SDDE_xandN_LEBagrarian03_bw_1")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_xandN_LEBagrarian03_bw_1big")

#-----------------------------------------------

 function keju_LEBagrarian_f(du,u,p,t)
nobhs,cshs,Freducn,IBIhs,AoMhs,aa,bb,HHa,KLhs,KtoX,keju,fleecing,xx0 = p
#u[1]=xxn
#u[2]=NNo
XXo= ((u[1]-xx0)^(1+bb))*((KtoX*((keju*u[2])^bb)*((1-keju)*u[2])*(HHa/(u[1]*u[2]))))*(1/(((u[1]/1.8)^20)+1))
Srn=u[1]/(u[1]+1)
LEBn=(KLhs*(u[1])^aa)

du[1]= 0.5*((XXo/(u[2]))-u[1])-fleecing
du[2]= (u[2]*Freducn*cshs*nobhs*Srn*(LEBn-AoMhs)/IBIhs/LEBn - u[2]/LEBn)
end

function keju_LEBagrarian_g(du,u,p,t)
nobhs,cshs,Freducn,IBIhs,AoMhs,aa,bb,HHa,KLhs,KtoX,keju,fleecing,xx0 = p
#u[1]=xxn
#u[2]=NNo
du[1]= 0.1*0.024*(u[1]+u[2])
du[2]= 0.1*0.006*(u[1]+u[2])
end


 function keju_LEBagrarian_SDDE_f(du,u,h,pdde,t)
nobhs,cshs,Freducn,IBIhs,AoMhs,aa,bb,HHa,KLhs,KtoX,keju,fleecing,xx0,lagrepro,lagdeath = pdde
#u[1]=xxn
#u[2]=NNo
histrepro = h(pdde, t-lagrepro)
histdeath = h(pdde, t-lagdeath)

XXo= ((histrepro[1]-xx0)^(1+bb))*((KtoX*((keju*histrepro[2])^bb)*((1-keju)*histrepro[2])*(HHa/(u[1]*u[2]))))*(1/(((u[1]/1.8)^20)+1))
Srn=u[1]/(u[1]+1)
LEBn=(KLhs*(u[1])^aa)

du[1]= 0.5*((XXo/(u[2]))-u[1])-fleecing
du[2]= (histrepro[2]*Freducn*cshs*nobhs*Srn*(LEBn-AoMhs)/IBIhs/LEBn - histdeath[2]/LEBn)
end

function keju_LEBagrarian_SDDE_g(du,u,h,pdde,t)
nobhs,cshs,Freducn,IBIhs,AoMhs,aa,bb,HHa,KLhs,KtoX,keju,fleecing,xx0,lagrepro,lagdeath = pdde
#u[1]=xxn
#u[2]=NNo
du[1]= 0.1*0.024*(u[1]+u[2])
du[2]= 0.1*0.006*(u[1]+u[2])
end

lagrepro=20.0
#average age of reproduction for venatorian homo sapiens
lagdeath=40.0
#average age of death for venatorian homo sapiens

pdde = [nobhs,cshs,Freducn,IBIhs,AoMhs,aa,bb,HHa,KLhs,KtoX,keju,fleecing,xx0,lagrepro,lagdeath]

p = [nobhs,cshs,Freducn,IBIhs,AoMhs,aa,bb,HHa,KLhs,KtoX,keju,fleecing,xx0]

u0 = [0.97, 0.85]; #Initial populations
tspan = (0.0,dur);
ODEprob = ODEProblem(keju_LEBagrarian_f,u0,tspan,p);
sol = solve(ODEprob);

SDEprob = SDEProblem(keju_LEBagrarian_f,keju_LEBagrarian_g,u0,tspan,p);
SDEsol = solve(SDEprob);

h(t, pdde) = [0.97,0.85]
DDEprob = DDEProblem(keju_LEBagrarian_SDDE_f, u0, h, tspan,pdde,constant_lags = [lagrepro,lagdeath])
DDEsol = solve(DDEprob)

SDDEprob = SDDEProblem(keju_LEBagrarian_SDDE_f,keju_LEBagrarian_SDDE_g, u0, h, tspan,pdde;constant_lags = [lagrepro,lagdeath])
SDDEsol = solve(SDDEprob,RKMil())

#################################
# write .csv of ODE, SDE, DDE, SDDE

CSV.write("xandN_agrarian03_ODE.csv",sol)
CSV.write("xandN_agrarian03_SDE.csv",SDEsol)
CSV.write("xandN_agrarian03_DDE.csv",DDEsol)
CSV.write("xandN_agrarian03_SDDE.csv",SDDEsol)


#####################################
# plot ODE, SDE, DDE, SDDE
#quiver, isoclines and trajectory 
#color

plt1 = plot(layout=(1,1),show = true,title="keju_xxn_NNo_LEBagrarian, SN=0.22",size=(600,400))
quiver!(xP,yS,quiver=(0.03*uu,0.03*vv), color=:lightgray)
contour!(0.3:0.001:a, 0.1:0.001:b, dxxnt, levels=[0], color=:cyan, colorbar = false, w=3)
contour!(0.3:0.001:a, 0.1:0.001:b, dNNot, levels=[0], color=:orange,colorbar = false, w=3)

plot!(plt1[1],sol,vars=(1,2),label="ODE", w=1)
plot!(plt1[1],SDEsol,vars=(1,2),label="SDE", w=2)
plot!(plt1[1],DDEsol,vars=(1,2),label="DDE", w=3)
plot!(plt1[1],SDDEsol,vars=(1,2),label="SDDE",w=1,
xlims = (0.3,a),ylims = (0.2,(b-0.01)),xaxis="[xxn]",yaxis="[N/No]",xticks = 0:0.1:a,yticks = 0:0.1:(b-0.01),grid= false, show = true)

#read(stdin, Char)
savefig("SDDE_xandN_LEBagrarian03_2")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_xandN_LEBagrarian03_2big")

#####################################
# plot ODE, SDE, DDE, SDDE
#quiver, isoclines and trajectory 
#BW

plt1 = plot(layout=(1,1),show = true,title="keju_xxn_NNo_LEBagrarian, SN=0.22",size=(600,400))
quiver!(xP,yS,quiver=(0.03*uu,0.03*vv), color=:black)
contour!(0.3:0.001:a, 0.1:0.001:b, dxxnt, levels=[0], color=:black, colorbar = false, w=4)
contour!(0.3:0.001:a, 0.1:0.001:b, dNNot, levels=[0], color=:black, colorbar = false, w=3)

plot!(plt1[1],sol,vars=(1,2),label="ODE",w=1, color=:black,linestyle=:dot)
plot!(plt1[1],SDEsol,vars=(1,2),label="SDE", w=2, color=:black,linestyle=:dashdot)
plot!(plt1[1],DDEsol,vars=(1,2),label="DDE", w=2, color=:black,linestyle=:dash)
plot!(plt1[1],SDDEsol,vars=(1,2),label="SDDE", w=1, color=:black,linestyle=:solid,
xlims = (0.3,a),ylims = (0.2,(b-0.01)),xaxis="[xxn]",yaxis="[N/No]",xticks = 0:0.1:a,yticks = 0:0.1:(b-0.01),grid= false, show = true)

#read(stdin, Char)
savefig("SDDE_xandN_LEBagrarian03_bw_2")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_xandN_LEBagrarian03_bw_2big")

##############################################
# plot ODE, SDE, DDE, SDDE
#Plot trajectories vs time
#color

plt2 = plot(layout=(1,1),title="keju_xxn_NNo_LEBagrarian, SN=0.22",size=(600,400))
plot!(plt2[1],sol,vars=(0,1),label="ODE xxn", w=1)
plot!(plt2[1],sol,vars=(0,2), label="ODE N/No", w=1)

plot!(plt2[1],SDEsol,vars=(0,1),label="SDE xxn",w=2)
plot!(plt2[1],SDEsol,vars=(0,2), label="SDE N/No",w=2)

plot!(plt2[1],DDEsol,vars=(0,1),label="DDE xxn",w=3)
plot!(plt2[1],DDEsol,vars=(0,2), label="DDE N/No",w=3)

plot!(plt2[1],SDDEsol,vars=(0,1),label="SDDE xxn",w=1)
plot!(plt2[1],SDDEsol,vars=(0,2), label="SDDE N/No",w=1,
xaxis="Time[Years]",yaxis="[xxn and N/No]",xticks = 0:(floor(dur/20)):dur,yticks = 0.5:0.1:1.5,grid= false,ylims = (0.5,1.5),show = true)

#read(stdin, Char)
savefig("SDDE_xandN_LEBagrarian03_3")

plot!(plt2[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_xandN_LEBagrarian03_3big")

##############################################
# plot ODE, SDE, DDE, SDDE
#Plot trajectories vs time
#BW

plt2 = plot(layout=(1,1),title="keju_xxn_NNo_LEBagrarian, SN=0.22",size=(600,400))
plot!(plt2[1],sol,vars=(0,1),label="ODE xxn", w=1,color=:black,linestyle=:dot)
plot!(plt2[1],sol,vars=(0,2), label="ODE N/No", w=3,color=:black,linestyle=:dot)

plot!(plt2[1],SDEsol,vars=(0,1),label="SDE xxn", w=1,color=:black,linestyle=:dashdot)
plot!(plt2[1],SDEsol,vars=(0,2), label="SDE N/No", w=3,color=:black,linestyle=:dashdot)

plot!(plt2[1],DDEsol,vars=(0,1),label="DDE xxn", w=1,color=:black,linestyle=:dash)
plot!(plt2[1],DDEsol,vars=(0,2), label="DDE N/No", w=3,color=:black,linestyle=:dash)

plot!(plt2[1],SDDEsol,vars=(0,1),label="SDDE xxn", w=1,color=:black,linestyle=:solid)
plot!(plt2[1],SDDEsol,vars=(0,2), label="SDDE N/No", w=3,color=:black,linestyle=:solid,
xaxis="Time[Years]",yaxis="[xxn and N/No]",xticks = 0:(floor(dur/20)):dur,yticks = 0.5:0.1:1.5,grid= false,ylims = (0.5,1.5),show = true)

#read(stdin, Char)
savefig("SDDE_xandN_LEBagrarian03_bw_3")

plot!(plt2[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_xandN_LEBagrarian03_bw_3big")

#end


