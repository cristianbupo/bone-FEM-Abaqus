import gmsh
import numpy as np

gmsh.initialize()
gmsh.model.add("Example")

# Simple geometry
gmsh.model.geo.addPoint(0, 0, 0, 1, 1)
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(0)

nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()

# View creation
view = gmsh.view.add("ExampleView")
gmsh.view.addModelData(view, 0, "", "NodeData", nodeTags, [[42.0] for _ in nodeTags])

# Show
gmsh.fltk.run()

# Save IMMEDIATELY after seeing
gmsh.write("example_output.pos")

gmsh.finalize()
