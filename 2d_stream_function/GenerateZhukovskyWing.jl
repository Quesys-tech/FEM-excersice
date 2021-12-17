using Plots, CSV, DataFrames

c = 1.0
ζ₀ = -0.1 * c + 0.1 * c * im

θ = [i * 2π for i = 0:0.125/2.0:1]
ζ = ζ₀ .+ abs(c - ζ₀) * exp.(im * θ)

z = ζ + 1.0 ./ ζ

ℜ(z::Complex) = z.re
ℑ(z::Complex) = z.im

#write aerofoil shape to  CSV file
df = DataFrame(x=ℜ.(z), y=ℑ.(z))
CSV.write("2d_stream_function/ZhukovskyWing.csv",df)

#Plot
plot(ℜ.(ζ), ℑ.(ζ), label = "\$\\zeta\$")
plot!(ℜ.(z), ℑ.(z), label = "\$z\$")
plot!(;aspect_ratio = 1.0, xtickfontsize=16, ytickfontsize=16, legendfontsize = 16, grid = false)
savefig("2d_stream_function/ZhukovskyWing.pdf")