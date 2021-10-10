using LinearAlgebra, Match
using CSV, DataFrames
using GLMakie


function area(xᵉ::Vector{Float64}, yᵉ::Vector{Float64})
    u = [xᵉ[2] - xᵉ[1] xᵉ[3] - xᵉ[1]]
    v = [yᵉ[2] - yᵉ[1] yᵉ[3] - yᵉ[1]]

    return abs(u[1] * v[2] - u[2] * v[1]) / 2.0
end

function local_matrix(xᵉ::Vector{Float64}, yᵉ::Vector{Float64})
    Aᵉ = zeros(Float64, 3, 3)
    bᵉ = [yᵉ[2] - yᵉ[3] yᵉ[3] - yᵉ[1] yᵉ[1] - yᵉ[2]] / (2 * area(xᵉ, yᵉ))
    cᵉ = [xᵉ[3] - xᵉ[2] xᵉ[1] - xᵉ[3] xᵉ[2] - xᵉ[1]] / (2 * area(xᵉ, yᵉ))

    for ij in CartesianIndices(Aᵉ)
        i, j = Tuple(ij)
        Aᵉ[i,j] = (bᵉ[i]bᵉ[j] + cᵉ[i]cᵉ[j]) * area(xᵉ, yᵉ)
    end

    return Aᵉ
end

function local_vector(xᵉ::Vector{Float64}, yᵉ::Vector{Float64}, f::Function)
    bᵉ = zeros(Float64, 3)
    ξ = [1 / 3 1 / 5 3 / 5 1 / 5]
    η = [1 / 3 1 / 5 1 / 5 3 / 5]
    W = [-27 / 96 25 / 96 25 / 96 25 / 96]
    
    Nᵉ(i, ξₚ, ηₚ) = @match i begin
        1 => 1 - ξₚ - ηₚ
        2 => ξₚ
        3 => ηₚ
    end

    x(ξₚ, ηₚ) = [Nᵉ(i, ξₚ, ηₚ) for i = 1:3] ⋅ xᵉ
    y(ξₚ, ηₚ) = [Nᵉ(i, ξₚ, ηₚ) for i = 1:3] ⋅ yᵉ

    for p = 1:length(W)
        for i = 1:3
            bᵉ[i] += 2 * f(x(ξ[p], η[p]), y(ξ[p], η[p])) * Nᵉ(i, ξ[p], η[p]) * area(xᵉ, yᵉ) * W[p]
        end
    end

    return bᵉ 
end

f(x,y) = 2 * π^2 * sin(π * x) * sin(π * y)
# f(x,y) = 1

nodes = CSV.read("2d_poisson_dirichlet/nodes.csv", DataFrame)
connectivity = CSV.read("2d_poisson_dirichlet/connectivity.csv", DataFrame)

A = zeros(Float64, size(nodes, 1), size(nodes, 1))
b = zeros(Float64, size(nodes, 1))

for k = 1:size(connectivity, 1)
    xᵉ = zeros(Float64, 3)
    yᵉ = zeros(Float64, 3)
    
    local2global = [connectivity[k,"node" * string(j)] for j = 1:3]

    for j = 1:3
        xᵉ[j] = nodes[local2global[j],"x"] 
        yᵉ[j] = nodes[local2global[j],"y"]
    end

    Aᵉ = local_matrix(xᵉ, yᵉ)
    bᵉ = local_vector(xᵉ, yᵉ, f)
    
    for ij in CartesianIndices(Aᵉ)
        i, j = Tuple(ij)
        A[local2global[i],local2global[j]] += Aᵉ[i,j]
    end
    for i = 1:3
        b[local2global[i]] += bᵉ[i]
    end
end

for row in CSV.File("2d_poisson_dirichlet/bc.csv")
    i = row.node
    b[i] = Float64(row.val)
    
    for j = 1:size(A, 2)
        if j == i
            A[i,j] = 1    
        else
            A[i,j] = 0
        end
    end

    for j = 1:size(A, 1)
        if j != i
            b[j] -= A[j,i] * b[i]
            A[j,i] = 0
        end
    end
end

u = A \ b

result = DataFrame(x=Float64[], y=Float64[], u=Float64[])
for i = 1:size(nodes, 1)
    push!(result, (nodes[i,"x"], nodes[i,"y"], u[i]))
end
println("Finished!")
CSV.write("2d_poisson_dirichlet/result.csv",result)