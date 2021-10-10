import meshio
import csv
points = []
u = []
cells =[]

with open("2d_poisson_dirichlet/result.csv") as csvfile:
    reader = csv.reader(csvfile)
    header = next(reader)
    for row in reader:
        points.append([float(row[0]), float(row[1]), 0])
        u.append(float(row[2]))

with open("2d_poisson_dirichlet/connectivity.csv") as csvfile:
    reader = csv.reader(csvfile)
    header = next(reader)
    for row in reader:
        cells.append([int(row[i]) -1 for i in range(3)])


mesh = meshio.Mesh(
    points,
    cells = [("triangle",cells)],
    point_data={"u":u}
)

mesh.write("2d_poisson_dirichlet/result.vtu")