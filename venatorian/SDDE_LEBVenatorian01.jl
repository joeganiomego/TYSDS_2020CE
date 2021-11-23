#=
/*
01 Version.  Nominal Conditions +20% Stochastic Random on average dRRot and dNNot
with delay of reproducing cohort and dying cohort.
*/
/* Stochastic DDE LEB Venatorian */
/* build 1.0 */
/* all data normalized with supposed equilibrium at R/Ro=1 and N/No=1 (or very close) */

=#
"""
/*
01 Version.  Nominal Conditions
Code run on Julia (1.5.3, julialang.org).
*/

/*
Conditions to calculate parameters:

Parameters estimation done at file : LEBvenatorianParameters01.ods

aa= 0.25 by choice

Strong Input/Output Targets:
at RRo=1 and NNo=1 then dRRot=0   (to have isocline cross at 1,1)
at RRo=1 (and NNo=1) then dNNot=0   (to have isocline cross at 1,1)
at RRo=2.2 and NNo=0 then dRRot=0   (to have carry capacity Kc = 2.2)( from: https://www1.maine.gov/ifw/pdfs/species_planning/mammals/moose/assessment.pdf )
at RRo=1 then xxr=1.0    (by choice) (nice number)
at RRo=1 (and NNo=1) then xxn=1.0 (by choice) (nice number)
at RRo=1 and NNo=0 then LEBr=15   (to have 15 year LEB herbivore (Averaged Large+Mixed+Birds))
at RRo=1 then LEBn=30 therefore KLn : 30.0 (paleolithic conditions HS LEB)
csr  : 0.5; copulation system
Freducr : 1.0 ; fertility reduction .  Considered 0 reduction at Herbivore conditions
nobr: 1.0 ; number of offsping per birth. set to 1.0.
csn  : 0.5; copulation system
Freducn : 1.0 ; fertility reduction .  Considered 0 reduction at paleo conditions
nobn: 1.0 ; number of offsping per birth. set to 1.0.
AoMn : 18.0 ; Age of menarche N (paleo HS)
IBIn  :  3.0 ;  Inter birth interval (paleo HS), includes correction for 50% miscarriges and death 
(Includes pregnancy(9 month) +LAM (14 months)+ refecundization probabylity time (6 months, 4 moths for modern humans)). 
Srr : xxr (RRo,NNo)/(xxr (RRo,NNo)+1); Perinatal and juvenile survival rate R (Herbivore) 
Srn : xxn (RRo,NNo)/(xxn (RRo,NNo)+1); Perinatal and juvenile survival rate N (paleo HS) (https://ourworldindata.org/child-mortality/)

(at Nominal conditions : stable equilibrium at optimal competitivity condition)
therefore:
At nominal conditions, at RRo=1 and NNo=1 crosspoint to have positive jacobian determinat and negative jacobian trace.
At nominal conditions dRRot=0 isocline to have max point at RRo<1.0 but close to it.

Weaker Input/Output Targets:
at RRo=1 and NNo=0 then dRRot about 0.08 ; (to have doubling population in 6 years.  Herbivore style)

Calculated (Respecting Plausibility):
gg: 1.0;  ( http://www.zo.utexas.edu/courses/bio301/chapters/Chapter15/Chapter15.html )
hh: 0.5; (with linear relation of grazing intensity vs dry mass.  https://www.nature.com/articles/srep16434  )
Ar : hh+1; =1.5;
AoMr : 3.5 ; Age of menarche R (Large + Mixed + Birds Herbivore) (calculated)
IBIr  :  1.2 ;  Inter birth interval (Herbivore), includes correction for 50% miscarriges and death (Calculated) 
Khunt: 13.4;
KRtox :(1+gg)/Khunt;
Khunt*(0.5/Srr) : Khunt correction to account for juvenile prays.

Failed desiderata;
RRo=2 and NNo=1 then dNNot about 0.04   (to have doubling population in 10 years.  HS style)
*/
"""

#using Dierckx
#using LsqFit
using DifferentialEquations
using StochasticDelayDiffEq
#using DataFrames
using CSV
using Plots

a=2.0 	# axis max */
b=2.0 	# axis max */
dur=200  # durat of t */

nobr =1.0 
csr  =0.5 
Freducr = 1.0 
IBIr  = 1.2 
AoMr =3.5  

nobn =1.0
csn  = 0.5 
Freducn =1.0  
IBIn  =3.0  
AoMn =18.0  

aa=0.25
gg=1.0
hh=0.5
Ar=1.5
KLr=15.0
Khunt=13.4

KRtox =(1+gg)/Khunt
KLn =30.0;


#--------------------------------------------------

# isoclines and quiver plot preparation

xxr(RRo,NNo)=Ar-hh*RRo 
xxn(RRo,NNo)=KRtox*Khunt*RRo/(RRo+gg)

Srr(RRo,NNo)=xxr(RRo,NNo)/(xxr(RRo,NNo)+1)
Srn(RRo,NNo)=xxn(RRo,NNo)/(xxn(RRo,NNo)+1)

LEBr(RRo,NNo) =(KLr*(xxr(RRo,NNo))^aa)-(Khunt*(0.5/Srr(RRo,NNo))*(RRo/(RRo+gg))*(NNo/RRo))
LEBn(RRo,NNo) =(KLn*(xxn(RRo,NNo))^aa)

dRRot(RRo,NNo)= RRo*Freducr*csr*nobr*Srr(RRo,NNo)*(LEBr(RRo,NNo)-AoMr)/IBIr/LEBr(RRo,NNo) - RRo/LEBr(RRo,NNo)  
dNNot(RRo,NNo)= NNo*Freducn*csn*nobn*Srn(RRo,NNo)*(LEBn(RRo,NNo)-AoMn)/IBIn/LEBn(RRo,NNo) - NNo/LEBn(RRo,NNo)  

#####################################

# plot Quiver and Isoclines

# quiver chart
meshgrid(xR, yN) = (repeat(xR, outer=length(yN)), repeat(yN, inner=length(xR)))
xR, yN = meshgrid(0:0.1:a, 0:0.1:b)
uu = @.dRRot(xR,yN)
vv = @.dNNot(xR,yN)

# isoclines color ###############
plt1 = plot(layout=(1,1),size=(600,400),title="LEBVenatorian")
quiver!(xR,yN,quiver=(0.5*uu,0.5*vv), color=:lightgray)
contour!(0:0.01:a, 0:0.01:b, dRRot, levels=[0], color=:cyan, colorbar = false, w=3)
contour!(0:0.01:a, 0:0.01:b, dNNot, levels=[0], color=:orange,colorbar = false, w=3)

plot!(plt1[1],
xlims = (0,a),ylims = (0,(b-0.5)),xaxis="[R/Ro]",yaxis="[N/No]",xticks = 0:0.25:a,yticks = 0:0.25:(b-0.5),grid= false, show = true)

#read(stdin, Char)
savefig("SDDE_LEBVenatorian01_1")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_LEBVenatorian01_1big")

# isoclines BW ########################
plt1 = plot(layout=(1,1),size=(600,400),title="LEBVenatorian")
quiver!(xR,yN,quiver=(0.5*uu,0.5*vv), color=:black)
contour!(0:0.01:a, 0:0.01:b, dRRot, levels=[0], color=:black, colorbar = false, w=4)
contour!(0:0.01:a, 0:0.01:b, dNNot, levels=[0], color=:black, colorbar = false, w=3)

plot!(plt1[1],
xlims = (0,a),ylims = (0,(b-0.5)),xaxis="[R/Ro]",yaxis="[N/No]",xticks = 0:0.25:a,yticks = 0:0.25:(b-0.5),grid= false, show = true)

#read(stdin, Char)
savefig("SDDE_LEBVenatorian01_bw_1")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_LEBVenatorian01_bw_1big")

#-----------------------------------------------

 function LEBVenatorian_f(du,u,p,t)
nobr,csr,Freducr,IBIr,AoMr,nobn,csn,Freducn,IBIn,AoMn,aa,gg,hh,Ar,KLr,Khunt,KRtox,KLn = p
#u[1]=RRo
#u[2]=NNo

xxr=Ar-hh*u[1] 
xxn=KRtox*Khunt*u[1]/(u[1]+gg)

Srr=xxr/(xxr+1)
Srn=xxn/(xxn+1)

LEBr =(KLr*(xxr)^aa)-(Khunt*(0.5/Srr)*(u[1]/(u[1]+gg))*(u[2]/u[1]))
LEBn =(KLn*(xxn)^aa)

du[1]= u[1]*Freducr*csr*nobr*Srr*(LEBr-AoMr)/IBIr/LEBr - u[1]/LEBr  
du[2]= u[2]*Freducn*csn*nobn*Srn*(LEBn-AoMn)/IBIn/LEBn - u[2]/LEBn 
end

function LEBVenatorian_g(du,u,p,t)
nobr,csr,Freducr,IBIr,AoMr,nobn,csn,Freducn,IBIn,AoMn,aa,gg,hh,Ar,KLr,Khunt,KRtox,KLn = p
#u[1]=RRo
#u[2]=NNo
du[1]= 0.1*0.03*(u[1]+u[2])
du[2]= 0.1*0.003*(u[1]+u[2])
# 0.1*0.03 and 0.1*0.003 estimated to acheive 20% Wiener noise.
end


 function LEBVenatorian_SDDE_f(du,u,h,pdde,t)
nobr,csr,Freducr,IBIr,AoMr,nobn,csn,Freducn,IBIn,AoMn,aa,gg,hh,Ar,KLr,Khunt,KRtox,KLn,lagrepro,lagdeath = pdde
#u[1]=RRo
#u[2]=NNo

histrepro = h(pdde, t-lagrepro)
histdeath = h(pdde, t-lagdeath)

xxr=Ar-hh*u[1] 
xxn=KRtox*Khunt*u[1]/(u[1]+gg)

Srr=xxr/(xxr+1)
Srn=xxn/(xxn+1)

LEBr =(KLr*(xxr)^aa)-(Khunt*(0.5/Srr)*(u[1]/(u[1]+gg))*(u[2]/u[1]))
LEBn =(KLn*(xxn)^aa)

du[1]= u[1]*Freducr*csr*nobr*Srr*(LEBr-AoMr)/IBIr/LEBr - u[1]/LEBr  
du[2]= histrepro[2]*Freducn*csn*nobn*Srn*(LEBn-AoMn)/IBIn/LEBn - histdeath[2]/LEBn 
end


function LEBVenatorian_SDDE_g(du,u,h,pdde,t)
nobr,csr,Freducr,IBIr,AoMr,nobn,csn,Freducn,IBIn,AoMn,aa,gg,hh,Ar,KLr,Khunt,KRtox,KLn,lagrepro,lagdeath = pdde
#u[1]=RRo
#u[2]=NNo
du[1]= 0.1*0.03*(u[1]+u[2])
du[2]= 0.1*0.003*(u[1]+u[2])
# 0.1*0.03 and 0.1*0.003 estimated to acheive 20% Wiener noise.
end


lagrepro=25.0
#average age of reproduction for venatorian homo sapiens
lagdeath=30.0
#average age of death for venatorian homo sapiens


pdde = [nobr,csr,Freducr,IBIr,AoMr,nobn,csn,Freducn,IBIn,AoMn,aa,gg,hh,Ar,KLr,Khunt,KRtox,KLn,lagrepro,lagdeath]

p = [nobr,csr,Freducr,IBIr,AoMr,nobn,csn,Freducn,IBIn,AoMn,aa,gg,hh,Ar,KLr,Khunt,KRtox,KLn]

u0 = [0.9, 0.9]; #Initial populations
tspan = (0.0,dur);
ODEprob = ODEProblem(LEBVenatorian_f,u0,tspan,p);
sol = solve(ODEprob);

SDEprob = SDEProblem(LEBVenatorian_f,LEBVenatorian_g,u0,tspan,p);
SDEsol = solve(SDEprob);

h(t, pdde) = [0.9, 0.9]
DDEprob = DDEProblem(LEBVenatorian_SDDE_f, u0, h, tspan,pdde,constant_lags = [lagrepro,lagdeath])
DDEsol = solve(DDEprob)

SDDEprob = SDDEProblem(LEBVenatorian_SDDE_f,LEBVenatorian_SDDE_g, u0, h, tspan,pdde;constant_lags = [lagrepro,lagdeath])
SDDEsol = solve(SDDEprob,RKMil())


#################################
# write .csv of ODE, SDE, DDE, SDDE

CSV.write("venatorian01_ODE.csv",sol)
CSV.write("venatorian01_SDE.csv",SDEsol)
CSV.write("venatorian01_DDE.csv",DDEsol)
CSV.write("venatorian01_SDDE.csv",SDDEsol)

#####################################
# plot ODE, SDE, DDE, SDDE
#quiver, isoclines and trajectory 
#color

plt1 = plot(layout=(1,1),size=(600,400),title="LEBVenatorian")
quiver!(xR,yN,quiver=(0.5*uu,0.5*vv), color=:lightgray)
contour!(0:0.01:a, 0:0.01:b, dRRot, levels=[0], color=:cyan, colorbar = false, w=3)
contour!(0:0.01:a, 0:0.01:b, dNNot, levels=[0], color=:orange,colorbar = false, w=3)

plot!(plt1[1],sol,vars=(1,2),label="ODE",w=1)
plot!(plt1[1],SDEsol,vars=(1,2),label="SDE",w=2)
plot!(plt1[1],DDEsol,vars=(1,2),label="DDE", w=3)
plot!(plt1[1],SDDEsol,vars=(1,2),label="SDDE", w=1,
xlims = (0,a),ylims = (0,(b-0.5)),xaxis="[R/Ro]",yaxis="[N/No]",xticks = 0:0.25:a,yticks = 0:0.25:(b-0.5),grid= false, show = true)

#read(stdin, Char)
savefig("SDDE_LEBVenatorian01_2")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_LEBVenatorian01_2big")

#####################################
# plot ODE, SDE, DDE, SDDE
#quiver, isoclines and trajectory 
#BW

plt1 = plot(layout=(1,1),size=(600,400),title="LEBVenatorian")
quiver!(xR,yN,quiver=(0.5*uu,0.5*vv), color=:black)
contour!(0:0.01:a, 0:0.01:b, dRRot, levels=[0], color=:black, colorbar = false, w=4)
contour!(0:0.01:a, 0:0.01:b, dNNot, levels=[0], color=:black, colorbar = false, w=3)

plot!(plt1[1],sol,vars=(1,2),label="ODE",w=1, color=:black,linestyle=:dot) 
plot!(plt1[1],SDEsol,vars=(1,2),label="SDE",w=2, color=:black,linestyle=:dashdot)
plot!(plt1[1],DDEsol,vars=(1,2),label="DDE", w=2, color=:black,linestyle=:dash)
plot!(plt1[1],SDDEsol,vars=(1,2),label="SDDE", w=1, color=:black,linestyle=:solid,
xlims = (0,a),ylims = (0,(b-0.5)),xaxis="[R/Ro]",yaxis="[N/No]",xticks = 0:0.25:a,yticks = 0:0.25:(b-0.5),grid= false, show = true)

#read(stdin, Char)
savefig("SDDE_LEBVenatorian01_bw_2")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_LEBVenatorian01_bw_2big")

####################################
# plot ODE, SDE, DDE, SDDE
#Plot trajectories vs time
#color

plt2 = plot(layout=(1,1),title="LEBVenatorian",size=(600,400))
plot!(plt2[1],sol,vars=(0,1),label="ODE R/Ro")
plot!(plt2[1],sol,vars=(0,2), label="ODE N/No")

plot!(plt2[1],SDEsol,vars=(0,1),label="SDE R/Ro",w=2)
plot!(plt2[1],SDEsol,vars=(0,2), label="SDE N/No",w=2)

plot!(plt2[1],DDEsol,vars=(0,1),label="DDE R/Ro",w=2)
plot!(plt2[1],DDEsol,vars=(0,2), label="DDE N/No",w=2)

plot!(plt2[1],SDDEsol,vars=(0,1),label="SDDE R/Ro",w=2)
plot!(plt2[1],SDDEsol,vars=(0,2), label="SDDE N/No",w=2,
xaxis="Time[Years]",yaxis="[R/Ro and N/No]",xticks = 0:(floor(dur/20)):dur,yticks = 0:0.25:b,grid= false,ylims = (0,b),show = true)

#read(stdin, Char)
savefig("SDDE_LEBVenatorian01_3")

plot!(plt2[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_LEBVenatorian01_3big")

####################################
# plot ODE, SDE, DDE, SDDE
#Plot trajectories vs time
#BW

plt2 = plot(layout=(1,1),title="LEBVenatorian",size=(600,400))
plot!(plt2[1],sol,vars=(0,1),label="ODE R/Ro",w=1,color=:black,linestyle=:dot)
plot!(plt2[1],sol,vars=(0,2), label="ODE N/No",w=3,color=:black,linestyle=:dot)

plot!(plt2[1],SDEsol,vars=(0,1),label="SDE R/Ro",w=1,color=:black,linestyle=:dashdot)
plot!(plt2[1],SDEsol,vars=(0,2), label="SDE N/No",w=3,color=:black,linestyle=:dashdot)

plot!(plt2[1],DDEsol,vars=(0,1),label="DDE R/Ro",w=1,color=:black,linestyle=:dash)
plot!(plt2[1],DDEsol,vars=(0,2), label="DDE N/No",w=3,color=:black,linestyle=:dash)

plot!(plt2[1],SDDEsol,vars=(0,1),label="SDDE R/Ro",w=1,color=:black,linestyle=:solid)
plot!(plt2[1],SDDEsol,vars=(0,2), label="SDDE N/No",w=3,color=:black,linestyle=:solid,
xaxis="Time[Years]",yaxis="[R/Ro and N/No]",xticks = 0:(floor(dur/20)):dur,yticks = 0:0.25:b,grid= false,ylims = (0,b), show = true)

#read(stdin, Char)
savefig("SDDE_LEBVenatorian01_bw_3")

plot!(plt2[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_LEBVenatorian01_bw_3big")

#end

