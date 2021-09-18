using CSV, DataFrames

M = 5
if length(ARGS) > 1
    M = parse(UInt, ARGS[1]) #要素数
end
N = M + 1

connectivity = DataFrame(element = 1:M, node1 = 1:M, node2 = 2:N)
nodes = DataFrame(node = 1:N, x = [(i - 1) / M for i = 1:N])


CSV.write("1d_poisson/connectivity.csv", connectivity)
CSV.write("1d_poisson/nodes.csv", nodes)