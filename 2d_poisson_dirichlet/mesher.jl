using CSV, DataFrames, Plots

N_elements = 5 # 1辺あたりの同じ形状の要素数
N_nodes = N_elements + 1 # 1辺あたりの節点数
nodes = DataFrame(x=Float64[], y=Float64[])

boundary = DataFrame(node=Int64[], val=Float64[])

for iy = 1:N_nodes
    for ix = 1:N_nodes
        x = (ix - 1) / N_elements
        y = (iy - 1) / N_elements
        push!(nodes, (x, y))

        if x == 1 || x == 0 || y == 0 || y == 1
            push!(boundary, (size(nodes,1), 0.0))
        end
    end
end

connectivity = DataFrame(node1=Int64[], node2=Int64[], node3=Int64[])
for iy = 1:N_elements
    for ix = 1:N_elements
        i = (iy - 1) * N_nodes + ix
        push!(connectivity, (i, i + 1, i + N_nodes))
        push!(connectivity, (i + 1, i + N_nodes + 1, i + N_nodes))
    end
end

CSV.write("2d_poisson_dirichlet/connectivity.csv", connectivity)
CSV.write("2d_poisson_dirichlet/nodes.csv", nodes)
CSV.write("2d_poisson_dirichlet/bc.csv", boundary)