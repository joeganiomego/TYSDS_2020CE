""" 
/* Simply Venatorian */
/* build 1.0 */
/* all data normalized with supposed equilibrium at R/Ro=1 and N/No=1 (or very close) */

/*
Conditions to calculate parameters:

at RRo=1 and NNo=1 then dRRot=0   (to have isocline cross at 1,1)
at RRo=1 (and NNo=1) then dNNot=0   (to have isocline cross at 1,1)
at RRo=2.2 and NNo=0 then dRRot=0  therefore Kc = 2.2 (to have carry capacity Kc = 2.2)

at RRo=inf(ite) and (RRo=inf/(gg+RRo=inf))*NNo*Khunt=1.2*((RRo=1/(gg+RRo=1))*NNo*Khunt)(1,1)  therefore  gg=0.2   (by feeling)
at RRo=1 and NNo=0 then dRRot about 0.08 therefore KFDr : 1.4667   (to have doubling population in 6 years.  Herbivore style)
all above generates Khunt= 0.96.

at RRo=2 and NNo=1 then dNNot about 0.04 therefore KFDn*Kxtoy = 5.5  (to have doubling population in 10 years.  HS style)
to copy LEBcombo Kfdn = 0.05  therefore Kxtoy = 110 .  Therefore Kdrn = 88.0 .

*/
"""

using Dierckx
using LsqFit
using DifferentialEquations
# using DiffEqParamEstim

using LeastSquaresOptim
using Optim

using Plots

# plotly()   # works AOK, can save, but opens a chain of browsers
# gr()  #simple and frequently do not work
# pgfplotsx()    # doesnt work
#  inspectdr()   #Works AOK , can't save!

# using PyPlot #doesent work at all.
# using Statistics
# using Loess
# using Polynomials
# using Calculus
# using ForwardDiff
# using Interpolations
# using DataInterpolations
# using SimplePCHIP

# import Pkg; Pkg.add("InspectDR")



# axis max
a= 2.5  
# axis max	 
b= 3.0 
# durat of t */
dur= 20.0   

Kc= 2.2  
gg= 0.2
KFDr = 1.4667
KFDn = 0.05
Khunt= 0.96 
Kndr= 88.0

"""/* OLD Kxtoy:110;   transfer hunted x into reproduced y, set to have RRo=1 */

/* transfer hunted x into reproduced y, set to have RRo=0.5 */
"""
Kxtoy=110   

#--------------------------------------------------

# isoclines and quiver plot preparation
dRRot(RRo,NNo)= RRo*KFDr*(1-RRo/Kc)-(RRo/(gg+RRo))*NNo*Khunt
dNNot(RRo,NNo)= NNo*KFDn*(((RRo/(gg+RRo))*Kxtoy*Khunt)-Kndr)

#####################################

# plot Quiver and Isoclines
# https://goropikari.github.io/PlotsGallery.jl/    Julia Plots Gallery
# https://gist.github.com/gizmaa/7214002   Various Julia plotting examples using PyPlot 


plt1 = plot(layout=(1,1),show = true)

# quiver chart
meshgrid(xR, yN) = (repeat(xR, outer=length(yN)), repeat(yN, inner=length(xR)))
xR, yN = meshgrid(0:0.15:a, 0:0.15:b)
u = @.dRRot(xR,yN)
v = @.dNNot(xR,yN)
quiver!(xR,yN,quiver=(0.03*u,0.03*v), color=:lightgray)

# isoclines
contour!(0:0.01:a, 0:0.01:b, dRRot, levels=[0], color=:cyan, colorbar = false, w=3)
contour!(0:0.01:a, 0:0.01:b, dNNot, levels=[0], color=:orange,colorbar = false, w=3)


#-----------------------------------------------


 function SimpVenatorian(du,u,p,t)
KFDr,Kc,gg,Khunt,KFDn,Kxtoy,Kndr = p
#u[1]=RRo
#u[2]=NNo
du[1]= u[1]*KFDr*(1-u[1]/Kc)-(u[1]/(gg+u[1]))*u[2]*Khunt
du[2]= u[2]*KFDn*(((u[1]/(gg+u[1]))*Kxtoy*Khunt)-Kndr)
#dRRot(RRo,NNo)= RRo*KFDr*(1-RRo/Kc)-(RRo/(gg+RRo))*NNo*Khunt
#dNNot(RRo,NNo)= NNo*KFDn*(((RRo/(gg+RRo))*Kxtoy*Khunt)-Kndr)
end
p = [KFDr,Kc,gg,Khunt,KFDn,Kxtoy,Kndr]
u0 = [0.8, 0.8]; #Initial populations
tspan = (0.0,dur);
prob = ODEProblem(SimpVenatorian,u0,tspan,p);
sol = solve(prob);

#####################################

# plot ODE
# https://goropikari.github.io/PlotsGallery.jl/    Julia Plots Gallery
# https://gist.github.com/gizmaa/7214002   Various Julia plotting examples using PyPlot 

plot!(plt1[1],sol,vars=(1,2),
xlims = (0,a),ylims = (0,b),xaxis="[R/Ro]",yaxis="[N/No]",xticks = 0:0.25:a,yticks = 0:0.25:b,grid= false,
title="SimpVenatorian",size=(1000,600),label=false, show = true)


read(stdin, Char)

plt2 = plot(layout=(1,1))
plot!(plt2[1],sol,vars=(0,1),label="R/Ro")
plot!(plt2[1],sol,vars=(0,2), label="N/No",
xaxis="Time[Years]",yaxis="[R/Ro and N/No]",xticks = 0:1.0:dur,yticks = 0:0.25:b,grid= false,ylims = (0,b),
title="SimpVenatorian",size=(1000,600), show = true)

read(stdin, Char)


