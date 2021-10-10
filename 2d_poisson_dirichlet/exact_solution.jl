using CSV, DataFrames, GLMakie

nodes = CSV.File("2d_poisson_dirichlet/nodes.csv") |> DataFrame
connectivity = CSV.File("2d_poisson_dirichlet/connectivity.csv") |> DataFrame

N = size(nodes, 1) 

u(x,y) = sin(pi * x) * sin(pi * y)

vertices = zeros(Float64, N, 2)
color = zeros(Float64, N)

for i = 1:N
    vertices[i,1] = nodes[i, "x"]
    vertices[i,2] = nodes[i, "y"]
    color[i] = u(vertices[i,1], vertices[i,2])
end

faces = zeros(Int64, size(connectivity, 1), 3)

for i = 1:size(connectivity, 1)
    faces[i,1] = connectivity[i, "node1"]
    faces[i,2] = connectivity[i, "node2"] 
    faces[i,3] = connectivity[i, "node3"]
end

scene = mesh(vertices, faces, color=color, shading=false)

