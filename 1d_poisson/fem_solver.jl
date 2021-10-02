using CSV, DataFrames, LinearAlgebra,Plots, LaTeXStrings

const q = 2

function local_matrix(x₁, x₂)
    lᵉ = x₂ - x₁
    kᵉ = Array{Float64}(undef, 2, 2)
    kᵉ[1, 1] = kᵉ[2, 2] = 1 / lᵉ
    kᵉ[2, 1] = kᵉ[1, 2] = -1 / lᵉ
    return kᵉ
end

function local_vector(x₁, x₂)
    lᵉ = x₂ - x₁
    fᵉ = ones(Float64, 2)
    fᵉ *= q * lᵉ / 2.0
    return fᵉ
end


nodes = CSV.read("1d_poisson/nodes.csv", DataFrame)

node_coordinates = nodes.x
num_nodes = length(node_coordinates)

K = zeros(Float64, num_nodes, num_nodes)
F = zeros(Float64, num_nodes)

connectivity = CSV.read("1d_poisson/connectivity.csv", DataFrame)
relations = [connectivity.node1 connectivity.node2]

#Neumann B.C.
node_neumann = 1
neuman_value = -2 #外向き方向微分に注意

# 全体マトリックスとベクトルを構築
for i = 1:size(relations, 1)
    x₁ = node_coordinates[relations[i, 1]]
    x₂ = node_coordinates[relations[i, 2]]

    kᵉ = local_matrix(x₁, x₂)
    fᵉ = local_vector(x₁, x₂)

    if relations[i, 2] == node_neumann
        fᵉ += [0, neuman_value]
    elseif relations[i, 1] == node_neumann
        fᵉ += [neuman_value, 0]
    end

    for j = 1:2
        for k = 1:2
            K[relations[i, j], relations[i, k]] += kᵉ[j, k]
        end
        F[relations[i, j]] += fᵉ[j]
    end
end

# Dirichlet B.C.
i = num_nodes
F[i] = 1

for j = 1:size(K, 1) # 境界条件となる列の対角成分を1にしてそれ以外を0にする
    if j == i
        K[i, j] = 1
    else
        K[i, j] = 0
    end
end

#行列の対称性を維持するために対応する列の係数を消去
for j = 1:size(K, 1)
    if j != i
        F[j] -= K[j, i] * F[i]
        K[j, i] = 0
    end
end

# 連立方程式を解く
U = K \ F

exact_function(x) = -x^2 + 2x #解析解

#プロット
plot(node_coordinates, U, marker = :circle, label = "FEM")
plot!(exact_function,0,1, label = "exact",xlabel=L"x",ylabel=L"u",legend=:bottomright)
savefig(string("1d_poisson/", size(relations, 1), "-elements.png"))