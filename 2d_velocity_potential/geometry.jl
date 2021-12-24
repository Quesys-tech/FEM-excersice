import Gmsh: gmsh

m = 0.02
p = 0.4
t = 0.12

yₜ(x) = 5t * (0.2969 * sqrt(x) - 0.1260x - 0.3516x^2 + 0.2843x^3 - 0.1015x^4)

y_c(x) = ifelse(x <= p, m / (p^2) * (2p * x - x^2), m / ((1 - p)^2) * ((1 - 2p) + 2p * x - x^2))
dy_c_dx(x) = ifelse(x <= p, m / (p^2) * (p - x), 2m / ((1 - p)^2) * (p - x))
θ(x) = atan(dy_c_dx(x))

xᵤ = []
yᵤ = []
xₗ = []
yₗ = []
n = 10
# multi resoluton range
# https://discourse.julialang.org/t/combine-multiple-ranges-with-different-resolution/36822/7
for x in sort(unique(vcat(range(0, 0.1, length = n), range(0.1, 1, length = n))))

    push!(xᵤ, x - yₜ(x) * sin(θ(x)))
    push!(yᵤ, y_c(x) + yₜ(x) * cos(θ(x)))
    push!(xₗ, x + yₜ(x) * sin(θ(x)))
    push!(yₗ, y_c(x) - yₜ(x) * cos(θ(x)))
end

lc = 1e-1 #mesh quality

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("NACA2412")

id_point = 1

for i = 1:size(xₗ, 1)
    gmsh.model.geo.addPoint(xₗ[i], yₗ[i], 0, lc, id_point)
    global id_point += 1
end
gmsh.model.geo.addSpline([i for i = 1:(id_point - 1)], 1)

id_point_u = id_point
for i = 2:size(xᵤ, 1) -1
    gmsh.model.geo.addPoint(xᵤ[i], yᵤ[i], 0, lc, id_point)
    global id_point += 1
end
gmsh.model.geo.addSpline([i for i = id_point_u:(id_point - 1)], 2)
gmsh.model.geo.addCurveLoop([1, -2], 1)

gmsh.model.geo.addPoint(-0.5, -0.5, 0, lc, id_point)
id_point += 1
gmsh.model.geo.addPoint(-0.5, 0.5, 0, lc, id_point)
id_point += 1
gmsh.model.geo.addPoint(1.5, 0.5, 0, lc, id_point)
id_point += 1
gmsh.model.geo.addPoint(1.5, -0.5, 0, lc, id_point)
id_point += 1

gmsh.model.geo.addLine(id_point - 4, id_point - 3, 3)
gmsh.model.geo.addLine(id_point - 3, id_point - 2, 4)
gmsh.model.geo.addLine(id_point - 2, id_point - 1, 5)
gmsh.model.geo.addLine(id_point - 4, id_point - 1, 6)

gmsh.model.geo.addCurveLoop([3, 4, 5, -6], 2)
gmsh.model.geo.addPlaneSurface([2], 1)
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)
gmsh.write("2d_velocity_potential/NACA2412.msh")
gmsh.finalize()
