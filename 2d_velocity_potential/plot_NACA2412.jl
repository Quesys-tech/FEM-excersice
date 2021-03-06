using Plots

#NACA 2412 のパラメータ
m = 0.02
p = 0.4
t = 0.12

yₜ(x) = 5t * (0.2969 * sqrt(x) - 0.1260x - 0.3516x^2 + 0.2843x^3 - 0.1015x^4)
y_c(x) = ifelse(x <= p, m / (p^2) * (2p * x - x^2), m / ((1 - p)^2) * ((1 - 2p) + 2p * x - x^2))
dy_c_dx(x) = ifelse(x <= p, m / (p^2) * (p - x), 2m / ((1 - p)^2) * (p - x))
θ(x) = atan(dy_c_dx(x))

xᵤ = []; yᵤ = []; xₗ = []; yₗ = []
# multi resoluton range
# https://discourse.julialang.org/t/combine-multiple-ranges-with-different-resolution/36822/7
for x in sort(unique(vcat(range(0, 0.1, length = 20), range(0.1, 1, length = 20))))
    push!(xᵤ, x - yₜ(x) * sin(θ(x)))
    push!(yᵤ, y_c(x) + yₜ(x) * cos(θ(x)))
    push!(xₗ, x + yₜ(x) * sin(θ(x)))
    push!(yₗ, y_c(x) - yₜ(x) * cos(θ(x)))
end
x_foil = vcat(xᵤ, reverse(xₗ))
y_foil = vcat(yᵤ, reverse(yₗ))

x_mean = [x for x in sort(unique(vcat(range(0, 0.1, length = 20), range(0.1, 1, length = 20))))]
y_mean = y_c.(x_mean)

plot(x_mean, y_mean, label = "\$y_c\$")
plot!(x_foil, y_foil, label = "NACA2412")
plot!(; aspect_ratio = 1.0, legendfontsize = 16,
    axis = nothing, showaxis = false, grid = false)
savefig("2d_velocity_potential/NACA2412.pdf")