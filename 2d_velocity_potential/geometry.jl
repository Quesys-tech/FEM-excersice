using Plots, CSV, DataFrames
import Gmsh: gmsh

c = 1.0
ζ₀ = -0.1 * c + 0.1 * c * im

θ = [i * 2π for i = 0:0.1/2.0:1]
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

#save gmsh geometry
gmsh.initialize()
gmsh.model.add("aerofoil")
lc = 1e-2

for i in 1:size(θ,1)
    gmsh.model.geo.addPoint(ℜ(z[i]), ℑ(z[i]), 0, lc, i)
end
gmsh.model.geo.addSpline([i for i=1:size(z,1)],1)

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)
gmsh.write("aerofoil.msh")
gmsh.finalize()
