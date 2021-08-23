using CSV, DataFrames, LinearAlgebra
using PyPlot

const q = 1

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

# 全体マトリックスとベクトルを構築
for i = 1:size(relations, 1)
    x₁ = node_coordinates[relations[i, 1]]
    x₂ = node_coordinates[relations[i, 2]]

    kᵉ = local_matrix(x₁, x₂)
    fᵉ = local_vector(x₁, x₂)
    for j = 1:2
        for k = 1:2
            K[relations[i, j], relations[i, k]] += kᵉ[j, k]
        end
        F[relations[i, j]] += fᵉ[j]
    end
end

# Dirichlet 境界条件を適用
i = 1
F[i] = 0

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
        K[j, i] = 0
    end
end

# 連立方程式を解く
U = K \ F

#解析解を計算
exact_function(x) = -0.5 * x^2 + x

x_exact = range(0, 1, step = 0.01)
u_exact = [exact_function(x) for x in x_exact]

#プロット
plot(node_coordinates, U, marker = "s", linestyle = "-", label = "FEM")
plot(x_exact, u_exact, label = "exact")
xlabel(raw"$x$")
ylabel(raw"$u$")
legend()
savefig(string("1d_poisson/", size(relations, 1), "-elements.png"))