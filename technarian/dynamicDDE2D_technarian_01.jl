""" 
Code run on Julia (1.5.3, julialang.org).
Technarian Age 2D version (mov, x/mov)

NOTE: DDE and SDDE generated unpredictable instability.  It is not clear if that was due to numerical error or 
real technarian world extreme sensitivity to delays, but the ambition to model Technarian age also with DDE and SDDE
has been cancelled.

NOTE2:  An attempt to model technarian with 3 parameters (mov, x, x0/mov0) has been done but it landed solutions prone to instability.
All that activity has been cancelled.

All parameters set to get frequency, amplitudes and shapes somewhat similar to the employement (mov) , GDPPC (xmov*mov) trends seen 
in the world economy in the period CE 1960 , 2020.

DDE2D_technarian_01.jl = version with no habitat constraints (HHa = 1.7)
DDE2D_technarian_02.jl = version with habitat constraints (HHa = 0.7)
dynamicDDE2D_technarian_01.jl = dynamic version with habitat constraint going from 3.7 to 0.7 in about 50 years.

"""

#using Dierckx
#using LsqFit
using DifferentialEquations
#using StochasticDelayDiffEq
#using DataFrames
using CSV
using Plots

a=1.5 	# axis max */
b=1.6 	# axis max */
dur=150  # durat of t */

# Nominal Version Parameters

HHa1=3.7    # Habitat enveloppe.  Basically max amount of x (bioenergy) that can be cointained by the habitat without running out of space.
HHa2=0.7    # Habitat enveloppe at which population growth stops (Habitat Carring Capacity)
optmov = 1.0    # Optimal Movement.  Species movement that has maximal efficiency.  Moving more then optmov will generate fatigue and efficiency xmov (x/xmov) will decrease.  Moving less then optmov will create some spare movement that will be used to increase efficiency.
x0mov0 = 1.0    # Efficiency that will generate zero profit for the organizers.  If xmov is above x0mov0 there is profit and therefore the interest or possibility to increase movement.


#--------------------------------------------------
# mov,xmov isoclines and quiver plot preparation
#--------------(mov,xmov)

Haoptmov1(mov,xmov) = (HHa1/mov)*(1/(((mov/optmov)^20)+1))
dxmovt1(mov,xmov) = 0.20*(Haoptmov1(mov,xmov)-xmov) +0.02*(0.25*(xmov-x0mov0))*10.0^(-0.25*(xmov-x0mov0))
dmovt1(mov,xmov)= 0.35*(xmov-x0mov0)  +0.07*dxmovt1(mov,xmov)*10.0^(-dxmovt1(mov,xmov))

Haoptmov2(mov,xmov) = (HHa2/mov)*(1/(((mov/optmov)^20)+1))
dxmovt2(mov,xmov) = 0.20*(Haoptmov2(mov,xmov)-xmov) +0.02*(0.25*(xmov-x0mov0))*10.0^(-0.25*(xmov-x0mov0))
dmovt2(mov,xmov)= 0.35*(xmov-x0mov0)  +0.07*dxmovt2(mov,xmov)*10.0^(-dxmovt2(mov,xmov))


#####################################

# plot Quiver and Isoclines

# quiver chart
meshgrid(xP, yS) = (repeat(xP, outer=length(yS)), repeat(yS, inner=length(xP)))
xP, yS = meshgrid(0.3:0.1:a, 0.1:0.1:b)
uu = @.dmovt1(xP,yS)
vv = @.dxmovt1(xP,yS)


# isoclines color ###############
plt1 = plot(layout=(1,1),size=(600,400),title="dynamicDDE2D Technarian")
quiver!(xP,yS,quiver=(0.03*uu,0.03*vv), color=:lightgray)
contour!(0.3:0.001:a, 0.1:0.001:b, dmovt1, levels=[0], color=:cyan, colorbar = false, w=3)
contour!(0.3:0.001:a, 0.1:0.001:b, dxmovt1, levels=[0], color=:orange,colorbar = false, w=3)
contour!(0.3:0.001:a, 0.1:0.001:b, dmovt2, levels=[0], color=:cyan, colorbar = false, w=1)
contour!(0.3:0.001:a, 0.1:0.001:b, dxmovt2, levels=[0], color=:orange,colorbar = false, w=1)

plot!(plt1[1],
xlims = (0.3,a),ylims = (0.2,(b-0.01)),xaxis="[mov]",yaxis="[xmov]",xticks = 0:0.1:a,yticks = 0:0.1:(b-0.01),grid= false, show = true)


#read(stdin, Char)
savefig("dynamicDDE2D_technarian_01_1")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("dynamicDDE2D_technarian_01_1big")

# isoclines BW ###############
plt1 = plot(layout=(1,1),size=(600,400),title="dynamicDDE2D Technarian")
quiver!(xP,yS,quiver=(0.03*uu,0.03*vv), color=:black)
contour!(0.3:0.001:a, 0.1:0.001:b, dmovt, levels=[0], color=:black, colorbar = false, w=4)
contour!(0.3:0.001:a, 0.1:0.001:b, dxmovt, levels=[0], color=:black, colorbar = false, w=3)
contour!(0.3:0.001:a, 0.1:0.001:b, dmovt2, levels=[0], color=:black, colorbar = false, w=3,linestyle=:dot)
contour!(0.3:0.001:a, 0.1:0.001:b, dxmovt2, levels=[0], color=:black, colorbar = false, w=2,linestyle=:dot)

plot!(plt1[1],
xlims = (0.3,a),ylims = (0.2,(b-0.01)),xaxis="[mov]",yaxis="[xmov]",xticks = 0:0.1:a,yticks = 0:0.1:(b-0.01),grid= false, show = true)


#read(stdin, Char)
savefig("dynamicDDE2D_technarian_01_bw_1")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("dynamicDDE2D_technarian_01_bw_1big")

#-----------------------------------------------


 function dDE2D_Technarian_f(du,u,p,t)
HHa2,optmov,x0mov0 = p
#u[1]=mov
#u[2]=xmov
#u[3]=HHa

Haoptmov = (u[3]/u[1])*(1/(((u[1]/optmov)^20)+1))
du[2] = 0.20*(Haoptmov-u[2]) +0.02*(0.25*(u[2]-x0mov0))*10.0^(-0.25*(u[2]-x0mov0))
du[1]= 0.35*(u[2]-x0mov0)  +0.07*du[2]*10.0^(-du[2])
du[3]= 0.15*(HHa2-u[2]*u[1])

end

function dDE2D_Technarian_g(du,u,p,t)
HHa2,optmov,x0mov0 = p
#u[1]=mov
#u[2]=xmov
#u[3]=HHa

du[1]= 0.05*0.20*(u[1]+u[2]+u[3])
du[2]= 0.05*0.06*(u[1]+u[2]+u[3])
du[3]= 0.05*0.06*(u[1]+u[2]+u[3])

end

#=
NOTE: DDE and SDDE discontinued because of strange instability and results

 function DE2D_Technarian_SDDE_f(du,u,p,pdde,t)
HHa,optmov,x0mov0,laginventory,lagmov = pdde
#u[1]=mov
#u[2]=xmov
histinventory = h(pdde, t-laginventory)
histmov = h(pdde, t-lagmov)
Haoptmov = (HHa/histmov[1])*(1/(((histmov[1]/optmov)^20)+1))
du[2] = 0.25*(Haoptmov-u[2])
du[1]= 0.25*(histmov[2]-x0mov0)  +0.2*(u[1]*u[2]-histinventory[1]*histinventory[2])*10.0^(-1.0*(u[1]*u[2]-histinventory[1]*histinventory[2]))    
end
function DE2D_Technarian_SDDE_g(du,u,p,pdde,t)
HHa,optmov,x0mov0,laginventory,lagmov = pdde
#u[1]=mov
#u[2]=xmov
du[1]= 0.01*0.15*(u[1]+u[2])
du[2]= 0.01*0.06*(u[1]+u[2])
end

laginventory= 0.25
#average inventory reconciliation time
lagmov=0.05
#average movement efficiency lag
pdde = [HHa,optmov,x0mov0,laginventory,lagmov]
=#
p = [HHa2,optmov,x0mov0]

u0 = [1.0, 1.0, HHa1]; #Initial populations
tspan = (0.0,dur);
ODEprob = ODEProblem(dDE2D_Technarian_f,u0,tspan,p);
sol = solve(ODEprob);

SDEprob = SDEProblem(dDE2D_Technarian_f,dDE2D_Technarian_g,u0,tspan,p);
SDEsol = solve(SDEprob);

#=
NOTE: DDE and SDDE discontinued because of strange instability and results
h(t, pdde) = [0.90,0.90] #Hystory populations
DDEprob = DDEProblem(DE2D_Technarian_SDDE_f, u0, h, tspan,pdde,constant_lags = [laginventory,lagmov])
DDEsol = solve(DDEprob)

SDDEprob = SDDEProblem(DE2D_Technarian_SDDE_f,DE2D_Technarian_SDDE_g, u0, h, tspan,pdde;constant_lags = [laginventory,lagmov])
SDDEsol = solve(SDDEprob,RKMil())
=#


#################################
# write .csv of ODE, SDE, DDE, SDDE

CSV.write("dynamic_technarian01_ODE.csv",sol)
CSV.write("dynamic_technarian01_SDE.csv",SDEsol)
#CSV.write("dynamic_technarian01_DDE.csv",DDEsol)
#CSV.write("dynamic_technarian01_SDDE.csv",SDDEsol)


#####################################

# plot ODE, SDE, DDE, SDDE
#quiver, isoclines and trajectory
#color

plt1 = plot(layout=(1,1),size=(600,400),title="dynamicDDE2D Technarian")
quiver!(xP,yS,quiver=(0.03*uu,0.03*vv), color=:lightgray)
contour!(0.3:0.001:a, 0.1:0.001:b, dmovt1, levels=[0], color=:cyan, colorbar = false, w=3)
contour!(0.3:0.001:a, 0.1:0.001:b, dxmovt1, levels=[0], color=:orange,colorbar = false, w=3)
contour!(0.3:0.001:a, 0.1:0.001:b, dmovt2, levels=[0], color=:cyan, colorbar = false, w=1)
contour!(0.3:0.001:a, 0.1:0.001:b, dxmovt2, levels=[0], color=:orange,colorbar = false, w=1)

plot!(plt1[1],sol,vars=(1,2),label="ODE", w=3)
plot!(plt1[1],SDEsol,vars=(1,2),label="SDE", w=1,
xlims = (0.3,a),ylims = (0.2,(b-0.01)),xaxis="[mov]",yaxis="[xmov]",xticks = 0:0.1:a,yticks = 0:0.1:(b-0.01),grid= false, show = true)

#plot!(plt1[1],SDEsol,vars=(1,2),label="SDE", w=2)
#plot!(plt1[1],DDEsol,vars=(1,2),label="DDE", w=3)
#plot!(plt1[1],SDDEsol,vars=(1,2),label="SDDE",w=4,
#xlims = (0.3,a),ylims = (0.2,(b-0.01)),xaxis="[mov]",yaxis="[xmov]",xticks = 0:0.1:a,yticks = 0:0.1:(b-0.01),grid= false,
#title="mov_xxn_SimpTechnarian",size=(1000,600), show = true)

#read(stdin, Char)
savefig("dynamicDDE2D_technarian_01_2")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("dynamicDDE2D_technarian_01_2big")


#####################################
# plot ODE, SDE, DDE, SDDE
#quiver, isoclines and trajectory
#BW

plt1 = plot(layout=(1,1),size=(600,400),title="dynamicDDE2D Technarian")
quiver!(xP,yS,quiver=(0.03*uu,0.03*vv), color=:black)
contour!(0.3:0.001:a, 0.1:0.001:b, dmovt, levels=[0], color=:black, colorbar = false, w=4)
contour!(0.3:0.001:a, 0.1:0.001:b, dxmovt, levels=[0], color=:black, colorbar = false, w=3)
contour!(0.3:0.001:a, 0.1:0.001:b, dmovt2, levels=[0], color=:black, colorbar = false, w=3,linestyle=:dot)
contour!(0.3:0.001:a, 0.1:0.001:b, dxmovt2, levels=[0], color=:black, colorbar = false, w=2,linestyle=:dot)

plot!(plt1[1],sol,vars=(1,2),label="ODE", w=2, color=:black,linestyle=:solid)
plot!(plt1[1],SDEsol,vars=(1,2),label="SDE", w=1, color=:black,linestyle=:dash,
xlims = (0.3,a),ylims = (0.2,(b-0.01)),xaxis="[mov]",yaxis="[xmov]",xticks = 0:0.1:a,yticks = 0:0.1:(b-0.01),grid= false, show = true)

#read(stdin, Char)
savefig("dynamicDDE2D_technarian_01_bw_2")

plot!(plt1[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("dynamicDDE2D_technarian_01_bw_2big")

###########################################################à
# plot ODE, SDE, DDE, SDDE
#Plot trajectories vs time
#color

plt2 = plot(layout=(1,1),size=(600,400),title="dynamicDDE2D Technarian")
plot!(plt2[1],sol,vars=(0,1),label="ODE mov")
plot!(plt2[1],sol,vars=(0,2), label="ODE xmov")

plot!(plt2[1],SDEsol,vars=(0,1),label="SDE mov")
plot!(plt2[1],SDEsol,vars=(0,2), label="SDE xmov",
xaxis="Time[Years]", show = true)

#read(stdin, Char)
savefig("dynamicDDE2D_technarian_01_3")

plot!(plt2[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("dynamicDDE2D_technarian_01_3big")

###########################################################à
# plot ODE, SDE, DDE, SDDE
#Plot trajectories vs time
#BW

plt2 = plot(layout=(1,1),size=(600,400),title="dynamicDDE2D Technarian")
plot!(plt2[1],sol,vars=(0,1),label="ODE mov",w=1,color=:black,linestyle=:dash)
plot!(plt2[1],sol,vars=(0,2), label="ODE xmov",w=1,color=:black,linestyle=:dashdot)

plot!(plt2[1],SDEsol,vars=(0,1),label="SDE mov", w=1,color=:black,linestyle=:solid)
plot!(plt2[1],SDEsol,vars=(0,2), label="SDE xmov", w=1,color=:black,linestyle=:dot,
xaxis="Time[Years]", show = true)

#read(stdin, Char)
savefig("dynamicDDE2D_technarian_01_bw_3")

plot!(plt2[1],size=(1000,600),show = true)

#read(stdin, Char)
savefig("dynamicDDE2D_technarian_01_bw_3big")

#end

