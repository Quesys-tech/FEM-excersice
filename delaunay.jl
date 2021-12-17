using GLMakie 
using GeometryBasics 
using Clustering
using Triangulate

# Construct a triangle 
outerpnts = [Point(0., 0), Point(1., 0.), Point(0.5, 1.)]
tri = Triangle(outerpnts...)

# Generate points that are inside of the triangle
innerpnts = Point[]
while length(innerpnts) ≤ 2000
    p = Point(rand(2)...)
    p ∈ tri && push!(innerpnts, p)
end
innerpntsmat = [getindex.(innerpnts, 1) getindex.(innerpnts, 2)]
innernodes = collect(kmeans(innerpntsmat', 100).centers')
outernodes = [getindex.(outerpnts, 1) getindex.(outerpnts, 2)]

# Triangulate the points.
triin = TriangulateIO() 
triin.pointlist = [outernodes; innernodes]'
msh, _ = triangulate("vcDQ", triin)

# Convert TriangleMesh.TriMesh to GeometryBasics.Mesh
pts = [Point(val[1], val[2]) for val in eachcol(msh.pointlist)]
fcs = [TriangleFace(val[1], val[2], val[3]) for val in eachcol(msh.trianglelist)]
msh = GeometryBasics.Mesh(pts, fcs)

println(fcs)

# # Plot mesh
# scn = Makie.mesh(msh, color = 1 : length(msh.position))
# Makie.wireframe!(msh) 