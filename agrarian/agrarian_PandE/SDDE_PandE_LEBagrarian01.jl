#=
/*
01 Version.  Nominal Conditions 
*/
/* Stochastic DDE LEB Agrarian People (P) and Elite (E) version  P+E=N  (N= total Population)*/
/* build 1.0 */
/* all data normalized with supposed equilibrium at P/No=0.9 and E/No=0.1 (or very close) */

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
bb= 0.24 by choice

Strong Input/Output Targets:
at NNo=1 with PNo= 0.9 and ENo=0.1 then dPNot=0 and dENot=0

All Other parameters set with common sense, trying to match a Produce/People with Defenders/Exploiters/Elite 
early agrarian society.  

cshs  : 0.5; copulation system of homo sapiens
Freducr : 1.0 ; fertility reduction .  Considered 0 reduction at Herbivore conditions
nobhs: 1.0 ; number of offsping per birth for humans. set to 1.0.
Freducp : 1.0 ; fertility reduction for P (people) .  Considered 0 reduction at paleo conditions
Freduce : 0.6 ; fertility reduction for E (Elite) .  Considered 0 reduction at paleo conditions
AoMhs : 18.0 ; Age of menarche for humans
IBIhs  :  3.0 ;  Inter birth interval for humans, includes correction for 50% miscarriges and death 
(Includes pregnancy(9 month) +LAM (14 months)+ refecundization probabylity time (6 months, 4 moths for modern humans)). 
Srp : xxr (RRo,NNo)/(xxr (RRo,NNo)+1); Perinatal and juvenile survival rate for P (People)
Sre : xxn (RRo,NNo)/(xxn (RRo,NNo)+1); Perinatal and juvenile survival rate for E (Elite)

Hha : 1.0  ; habitat container size
Eor : 5.0 ; Elite opulence ratio , chosen by choice
Klhs : 30 ; constant linking x to LEB for humans.
Kev : 220 ; constant for elite violence.  The higher the number the more violence there is in the elite.
KtoX : 3.0; constant linking E and P numbers to the x output.

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

a=1.2 	# axis max */
b=0.14	# axis max */
dur=1000  # durat of t */

nobhs =1.0 
cshs  =0.5 
Freducp = 1.0 
Freduce = 0.6
IBIhs  = 3.0 
AoMhs =18.0  

aa=0.25
bb=0.24
HHa=1.0
EOr=5.0
KLhs=30.0
Kev=220.0
KtoX=3.0

#--------------------------------------------------

# isoclines and quiver plot preparation

NNo(PNo,ENo)=PNo+ENo

XXraw(PNo,ENo)= abs(KtoX*((ENo)^bb)*(PNo))
xxnraw(PNo,ENo)= XXraw(PNo,ENo)/(PNo+ENo)
XXo(PNo,ENo)= ((xxnraw(PNo,ENo)-0.2)^(1+bb))*XXraw(PNo,ENo)*(HHa/(xxnraw(PNo,ENo)*(PNo+ENo)))*(1/(((xxnraw(PNo,ENo)/1.8)^20)+1))

xxp(PNo,ENo)= XXo(PNo,ENo)/(PNo+EOr*ENo)
xxe(PNo,ENo)= EOr*xxp(PNo,ENo)

Srp(PNo,ENo)=xxp(PNo,ENo)/(xxp(PNo,ENo)+1)
Sre(PNo,ENo)=xxe(PNo,ENo)/(xxe(PNo,ENo)+1)

LEBp(PNo,ENo) =(KLhs*(xxp(PNo,ENo))^aa)

LEBe(PNo,ENo) =(KLhs*(xxe(PNo,ENo))^aa)-(Kev*(ENo/PNo)*(0.5/Sre(PNo,ENo)))

dPNot(PNo,ENo)= PNo*Freducp*cshs*nobhs*Srp(PNo,ENo)*(LEBp(PNo,ENo)-AoMhs)/IBIhs/LEBp(PNo,ENo) - PNo/LEBp(PNo,ENo)  
dENot(PNo,ENo)= ENo*Freduce*cshs*nobhs*Sre(PNo,ENo)*(LEBe(PNo,ENo)-AoMhs)/IBIhs/LEBe(PNo,ENo) - ENo/LEBe(PNo,ENo)   


#####################################

# plot Quiver and Isoclines

# quiver chart
meshgrid(xP, yE) = (repeat(xP, outer=length(yE)), repeat(yE, inner=length(xP)))
xP, yE = meshgrid(0.5:0.05:a, 0.05:0.005:b)
uu = @.dPNot(xP,yE)
vv = @.dENot(xP,yE)

# isoclines color ###############
plt1 = plot(layout=(1,1),title="P and E LEBagrarian",size=(600,400))
quiver!(xP,yE,quiver=(0.5*uu,0.5*vv), color=:lightgray)
contour!(0.5:0.001:a, 0.05:0.001:b, dPNot, levels=[0], color=:cyan, colorbar = false, w=3)
contour!(0.5:0.001:a, 0.05:0.001:b, dENot, levels=[0], color=:orange, colorbar = false, w=3)

plot!(plt1[1],
xlims = (0.6,a),ylims = (0.05,(b-0.01)),xaxis="[P/No]",yaxis="[E/No]",xticks = 0:0.1:a,yticks = 0:0.01:b,grid= false, show = true)

#read(stdin, Char)
savefig("SDDE_PandE_LEBagrarian01_1")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_PandE_LEBagrarian01_1big")

# isoclines BW ###############
plt1 = plot(layout=(1,1),title="P and E LEBagrarian",size=(600,400))
quiver!(xP,yE,quiver=(0.5*uu,0.5*vv), color=:black)
contour!(0.5:0.001:a, 0.05:0.001:b, dPNot, levels=[0], color=:black, colorbar = false, w=4)
contour!(0.5:0.001:a, 0.05:0.001:b, dENot, levels=[0], color=:black, colorbar = false, w=3)

plot!(plt1[1],
xlims = (0.6,a),ylims = (0.05,(b-0.01)),xaxis="[P/No]",yaxis="[E/No]",xticks = 0:0.1:a,yticks = 0:0.01:b,grid= false, show = true)

#read(stdin, Char)
savefig("SDDE_PandE_LEBagrarian01_bw_1")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_PandE_LEBagrarian01_bw_1big")

#-----------------------------------------------

 function PE_LEBagrarian_f(du,u,p,t)
nobhs,cshs,Freducp,Freduce,IBIhs,AoMhs,aa,bb,HHa,EOr,KLhs,Kev,KtoX = p

#u[1]=PNo
#u[2]=ENo

NNo=u[1]+u[2]

XXraw = abs(KtoX*((u[2])^bb)*(u[1]))
xxnraw = XXraw/(u[1]+u[2])
XXo = ((xxnraw-0.2)^(1+bb))*XXraw*(HHa/(xxnraw*(u[1]+u[2])))*(1/(((xxnraw/1.8)^20)+1))

xxp= XXo/(u[1]+EOr*u[2])
xxe= EOr*xxp

Srp=xxp/(xxp+1)
Sre=xxe/(xxe+1)

LEBp =(KLhs*(xxp)^aa)
LEBe =(KLhs*(xxe)^aa)-(Kev*(u[2]/u[1])*(0.5/Sre))

du[1]= u[1]*Freducp*cshs*nobhs*Srp*(LEBp-AoMhs)/IBIhs/LEBp - u[1]/LEBp  
du[2]= u[2]*Freduce*cshs*nobhs*Sre*(LEBe-AoMhs)/IBIhs/LEBe - u[2]/LEBe
end

function PE_LEBagrarian_g(du,u,p,t)
nobhs,cshs,Freducp,Freduce,IBIhs,AoMhs,aa,bb,HHa,EOr,KLhs,Kev,KtoX = p

du[1]= 0.1*0.006*(u[1]+u[2])
du[2]= 0.1*0.024*(u[1]+u[2])
# 0.1*0.03 and 0.1*0.003 estimated to acheive 20% Wiener noise.
end


 function PE_LEBagrarian_SDDE_f(du,u,h,pdde,t)
 
nobhs,cshs,Freducp,Freduce,IBIhs,AoMhs,aa,bb,HHa,EOr,KLhs,Kev,KtoX,lagrepro,lagdeath = pdde
 
#u[1]=PNo
#u[2]=ENo

histrepro = h(pdde, t-lagrepro)
histdeath = h(pdde, t-lagdeath)

NNo=u[1]+u[2]

XXraw = abs(KtoX*((histrepro[2])^bb)*(histrepro[1]))
xxnraw = XXraw/(u[1]+u[2])
XXo = ((xxnraw-0.2)^(1+bb))*XXraw*(HHa/(xxnraw*(u[1]+u[2])))*(1/(((xxnraw/1.8)^20)+1))

xxp= XXo/(u[1]+EOr*u[2])
xxe= EOr*xxp

Srp=xxp/(xxp+1)
Sre=xxe/(xxe+1)

LEBp =(KLhs*(xxp)^aa)
LEBe =(KLhs*(xxe)^aa)-(Kev*(u[2]/u[1])*(0.5/Sre))

du[1]= histrepro[1]*Freducp*cshs*nobhs*Srp*(LEBp-AoMhs)/IBIhs/LEBp - histdeath[1]/LEBp  
du[2]= histrepro[2]*Freduce*cshs*nobhs*Sre*(LEBe-AoMhs)/IBIhs/LEBe - histdeath[2]/LEBe
end


function PE_LEBagrarian_SDDE_g(du,u,h,pdde,t)
nobhs,cshs,Freducp,Freduce,IBIhs,AoMhs,aa,bb,HHa,EOr,KLhs,Kev,KtoX,lagrepro,lagdeath = pdde

du[1]= 0.1*0.006*(u[1]+u[2])
du[2]= 0.1*0.024*(u[1]+u[2])
# 0.1*0.03 and 0.1*0.003 estimated to acheive 20% Wiener noise.
end


lagrepro=20.0
#average age of reproduction for agrarian homo sapiens
lagdeath=40.0
#average age of death for agrarian homo sapiens

pdde = [nobhs,cshs,Freducp,Freduce,IBIhs,AoMhs,aa,bb,HHa,EOr,KLhs,Kev,KtoX,lagrepro,lagdeath]

p = [nobhs,cshs,Freducp,Freduce,IBIhs,AoMhs,aa,bb,HHa,EOr,KLhs,Kev,KtoX]

u0 = [0.90, 0.10]; #Initial populations
tspan = (0.0,dur);
ODEprob = ODEProblem(PE_LEBagrarian_f,u0,tspan,p);
sol = solve(ODEprob);


SDEprob = SDEProblem(PE_LEBagrarian_f,PE_LEBagrarian_g,u0,tspan,p);
SDEsol = solve(SDEprob);

h(t, pdde) = [0.90,0.10] #Hystory populations
DDEprob = DDEProblem(PE_LEBagrarian_SDDE_f, u0, h, tspan,pdde,constant_lags = [lagrepro,lagdeath])
DDEsol = solve(DDEprob)

SDDEprob = SDDEProblem(PE_LEBagrarian_SDDE_f,PE_LEBagrarian_SDDE_g, u0, h, tspan,pdde;constant_lags = [lagrepro,lagdeath])
SDDEsol = solve(SDDEprob,RKMil())


#################################
# write .csv of ODE, SDE, DDE, SDDE

CSV.write("PandE_agrarian01_ODE.csv",sol)
CSV.write("PandE_agrarian01_SDE.csv",SDEsol)
CSV.write("PandE_agrarian01_DDE.csv",DDEsol)
CSV.write("PandE_agrarian01_SDDE.csv",SDDEsol)


#####################################
# plot ODE, SDE, DDE, SDDE
#quiver, isoclines and trajectory 
#color

plt1 = plot(layout=(1,1),title="P and E LEBagrarian",size=(600,400))
quiver!(xP,yE,quiver=(0.5*uu,0.5*vv), color=:lightgray)
contour!(0.5:0.001:a, 0.05:0.001:b, dPNot, levels=[0], color=:cyan, colorbar = false, w=3)
contour!(0.5:0.001:a, 0.05:0.001:b, dENot, levels=[0], color=:orange, colorbar = false, w=3)

plot!(plt1[1],sol,vars=(1,2),label="ODE")
plot!(plt1[1],SDEsol,vars=(1,2),label="SDE",w=2)
plot!(plt1[1],DDEsol,vars=(1,2),label="DDE",w=3)
plot!(plt1[1],SDDEsol,vars=(1,2),label="SDDE",w=1,
xlims = (0.6,a),ylims = (0.05,(b-0.01)),xaxis="[P/No]",yaxis="[E/No]",xticks = 0:0.1:a,yticks = 0:0.01:b,grid= false, show = true)

#read(stdin, Char)
savefig("SDDE_PandE_LEBagrarian01_2")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_PandE_LEBagrarian01_2big")

#####################################
# plot ODE, SDE, DDE, SDDE
#quiver, isoclines and trajectory 
#BW

plt1 = plot(layout=(1,1),title="P and E LEBagrarian",size=(600,400))
quiver!(xP,yE,quiver=(0.5*uu,0.5*vv), color=:black)
contour!(0.5:0.001:a, 0.05:0.001:b, dPNot, levels=[0], color=:black, colorbar = false, w=4)
contour!(0.5:0.001:a, 0.05:0.001:b, dENot, levels=[0], color=:black, colorbar = false, w=3)

plot!(plt1[1],sol,vars=(1,2),label="ODE",w=1, color=:black,linestyle=:dot)
plot!(plt1[1],SDEsol,vars=(1,2),label="SDE",w=2, color=:black,linestyle=:dashdot)
plot!(plt1[1],DDEsol,vars=(1,2),label="DDE", w=2, color=:black,linestyle=:dash)
plot!(plt1[1],SDDEsol,vars=(1,2),label="SDDE", w=1, color=:black,linestyle=:solid,
xlims = (0.6,a),ylims = (0.05,(b-0.01)),xaxis="[P/No]",yaxis="[E/No]",xticks = 0:0.1:a,yticks = 0:0.01:b,grid= false, show = true)

#read(stdin, Char)
savefig("SDDE_PandE_LEBagrarian01_bw_2")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_PandE_LEBagrarian01_bw_2big")

##############################################
# plot ODE, SDE, DDE, SDDE
#Plot trajectories vs time
#color

plt2 = plot(layout=(1,1),title="P and E LEBagrarian",size=(600,400))
plot!(plt2[1],sol,vars=(0,1),label="ODE P/No")
plot!(plt2[1],sol,vars=(0,2), label="ODE E/No")

plot!(plt2[1],SDEsol,vars=(0,1),label="SDE P/No",w=2)
plot!(plt2[1],SDEsol,vars=(0,2), label="SDE E/No",w=2)

plot!(plt2[1],DDEsol,vars=(0,1),label="DDE P/No",w=3)
plot!(plt2[1],DDEsol,vars=(0,2), label="DDE E/No",w=3)

plot!(plt2[1],SDDEsol,vars=(0,1),label="SDDE P/No",w=1)
plot!(plt2[1],SDDEsol,vars=(0,2), label="SDDE E/No",w=1,
xaxis="Time[Years]",yaxis="[P/No and E/No]",xticks = 0:(floor(dur/20)):dur,yticks = 0:0.05:a,grid= false,ylims = (0,a), show = true)

#read(stdin, Char)
savefig("SDDE_PandE_LEBagrarian01_3")

plot!(plt2[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_PandE_LEBagrarian01_3big")

##############################################
# plot ODE, SDE, DDE, SDDE
#Plot trajectories vs time
#BW

plt2 = plot(layout=(1,1),title="P and E LEBagrarian",size=(600,400))
plot!(plt2[1],sol,vars=(0,1),label="ODE P/No",w=1,color=:black,linestyle=:dot)
plot!(plt2[1],sol,vars=(0,2), label="ODE E/No",w=3,color=:black,linestyle=:dot)

plot!(plt2[1],SDEsol,vars=(0,1),label="SDE P/No",w=1,color=:black,linestyle=:dashdot)
plot!(plt2[1],SDEsol,vars=(0,2), label="SDE E/No",w=3,color=:black,linestyle=:dashdot)

plot!(plt2[1],DDEsol,vars=(0,1),label="DDE P/No",w=1,color=:black,linestyle=:dash)
plot!(plt2[1],DDEsol,vars=(0,2), label="DDE E/No",w=3,color=:black,linestyle=:dash)

plot!(plt2[1],SDDEsol,vars=(0,1),label="SDDE P/No",w=1,color=:black,linestyle=:solid)
plot!(plt2[1],SDDEsol,vars=(0,2), label="SDDE E/No",w=3,color=:black,linestyle=:solid,
xaxis="Time[Years]",yaxis="[P/No and E/No]",xticks = 0:(floor(dur/20)):dur,yticks = 0:0.05:a,grid= false,ylims = (0,a), show = true)

#read(stdin, Char)
savefig("SDDE_PandE_LEBagrarian01_bw_3")

plot!(plt2[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("SDDE_PandE_LEBagrarian01_bw_3big")

#end

