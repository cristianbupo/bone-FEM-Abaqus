import gmsh
gmsh.initialize()
gmsh.model.add("fixed_local_refine")

lc = 1

# Points
p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
p2 = gmsh.model.geo.addPoint(1, 0, 0, lc)
p3 = gmsh.model.geo.addPoint(1, 1, 0, lc)
p4 = gmsh.model.geo.addPoint(0, 1, 0, lc)

# Lines
l1 = gmsh.model.geo.addLine(p1, p2)  # bottom
l2 = gmsh.model.geo.addLine(p2, p3)  # right
l3 = gmsh.model.geo.addLine(p3, p4)  # top
l4 = gmsh.model.geo.addLine(p4, p1)  # left

# Surface
cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
s = gmsh.model.geo.addPlaneSurface([cl])

# Transfinite curves (matching number of points on opposite sides)
Nx = 40  # horizontal edges
Ny = 30  # vertical edges

# Apply progression to only one of each pair (refine bottom & left sides)
gmsh.model.geo.mesh.setTransfiniteCurve(l1, Nx, "Bump", 0.1)  # bottom: refined start
gmsh.model.geo.mesh.setTransfiniteCurve(l3, Nx, "Bump", 0.1)                 # top: uniform
gmsh.model.geo.mesh.setTransfiniteCurve(l2, Ny)                      # right: uniform
gmsh.model.geo.mesh.setTransfiniteCurve(l4, Ny)  # left: refined start

gmsh.model.geo.mesh.setTransfiniteSurface(s)
gmsh.model.geo.mesh.setRecombine(2, s)

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)
gmsh.write("fixed_local_refine.msh")
gmsh.fltk.run()
gmsh.finalize()
