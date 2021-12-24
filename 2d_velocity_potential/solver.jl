using LinearAlgebra
using StaticArrays


function area(xᵉ::Vector{Float64}, yᵉ::Vector{Float64})
    u = [xᵉ[2] - xᵉ[1] xᵉ[3] - xᵉ[1]]
    v = [yᵉ[2] - yᵉ[1] yᵉ[3] - yᵉ[1]]

    return abs(u[1] * v[2] - u[2] * v[1]) / 2.0
end

function local_matrix(xᵉ::Vector{Float64}, yᵉ::Vector{Float64})
    Aᵉ = zeros(SVector{3,3})
    bᵉ = zeros(SVector{3})
    cᵉ = zeros(SVector{3})

    perms = [1,2,3]
    for i = 1:3
        j = perms[(i+1) % 3]
        k = perms[(i+2) % 3]
        
        bᵉ[i] = yᵉ[j] - yᵉ[k]
        cᵉ[i] = xᵉ[j] - xᵉ[k]
    end
    bᵉ /= 2.0*area(xᵉ,yᵉ)
    cᵉ /= 2.0*area(xᵉ,yᵉ)

    for ij in CartesianIndices(Aᵉ)
        i, j = Tuple(ij)
        Aᵉ[i,j] = (bᵉ[i]bᵉ[j] + cᵉ[i]cᵉ[j]) * area(xᵉ, yᵉ)
    end

    return Aᵉ
end