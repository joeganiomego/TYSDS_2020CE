

![growth chart](https://github.com/joeganiomego/TYSDS_2020CE/blob/main/In_Caput_Evolution/TheSwitchTECHLEVELGrowth_01_1big.png)



MPADDE_inCaputEvolution_01.jl  :
Main Allometric model of homo sapiens species
MPADDE   Manual Parameter Adjustement DDE style



Best Retained Models : model5  and model6



tau = 3
lags = [tau]

function PopuGDPPCnorm5(dx,x,h,p,t)
> popu derivative with polynomial+exponential technical complexity + Popu density correction
dx[1] = paramp1[1]*x[1]^(paramp1[2])*(x[2]^(paramp1[3]))*2^(x[2]*(paramp1[4]))+(-0.055)*(1/(1+h(p, t-tau)[1]))
> GDPPC derivative with polynomial+exponential technical complexity + GDPPC level and sentiment correction
dx[2] = paramg1[1]*x[1]^(paramg1[2])*(x[2]^(paramg1[3]))*2^(x[2]*(paramg1[4]))+(-0.024)*(x[2]/(5.37+x[2]))*
(
((0.06)*h(p, t-tau)[2])/
(abs(x[2]-h(p, t-tau)[2])+((0.06)*h(p, t-tau)[2]))
)
end


function PopuGDPPCnorm6(dx,x,h,p,t)
> popu derivative with polynomial+stirling technical complexity + Popu density correction
dx[1] = paramp2[1]*x[1]^(paramp2[2])*(x[2]^(paramp2[3]))*(x[2]^(x[2]*(paramp2[4])))+(-0.0552)*(1/(1+h(p, t-tau)[1]))
> GDPPC derivative with polynomial+stirling technical complexity + GDPPC level and sentiment correction
dx[2] = paramg2[1]*x[1]^(paramg2[2])*(x[2]^(paramg2[3]))*(x[2]^(x[2]*(paramg2[4])))+(-0.0235)*(x[2]/(5.159+x[2]))*
(
((0.0602)*h(p, t-tau)[2])/
(abs(x[2]-h(p, t-tau)[2])+((0.0602)*h(p, t-tau)[2]))
)
end
