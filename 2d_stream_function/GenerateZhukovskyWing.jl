using Plots

c = 1.0
ξ₀ = -0.1*c
η₀ = 0.1*c

θ = [i*2π for i in 0:0.125/2.0:1]
r = sqrt(η₀^2 + (c-ξ₀)^2)

ξ = r*cos.(θ).+ξ₀
η = r*sin.(θ).+η₀

x = ξ.+(c^2*ξ)./(ξ.^2+η.^2)
y = η.-(c^2*η)./(ξ.^2+η.^2)

plot(ξ,η,label="\$\\xi\\,\\mathrm{space}\$",aspect_ratio=1.0,legendfontsize=16)
plot!(x,y,label="\$z\\,\\mathrm{space}\$",aspect_ratio=1.0,legendfontsize=16)
savefig("2d_stream_function/ZhukovskyWing.pdf")