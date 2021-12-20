using Plots, CSV, DataFrames
import Gmsh: gmsh

c = 1.0
ζ₀ = -0.1 * c + 0.1 * c * im
n = 36

θ = [i * 2π / n for i = 0:n]
ζ = ζ₀ .+ abs(c - ζ₀) * exp.(im * θ)

z = ζ + 1.0 ./ ζ

ℜ(z::Complex) = z.re
ℑ(z::Complex) = z.im

#write aerofoil shape to  CSV file
df = DataFrame(x=ℜ.(z), y=ℑ.(z))
CSV.write("2d_velocity_potential/ZhukovskyWing.csv",df)

#Plot
plot(ℜ.(ζ), ℑ.(ζ), label = "\$\\zeta\$")
plot!(ℜ.(z), ℑ.(z), label = "\$z\$")
plot!(;aspect_ratio = 1.0, xtickfontsize=16, ytickfontsize=16, legendfontsize = 16, grid = false)
savefig("2d_velocity_potential/ZhukovskyWing.pdf")
