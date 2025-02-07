import vtk
import numpy as np

# Parameters for the semicircles
outer_radius = 1.0
inner_radius = 0.5
num_points = 50
num_radial_divisions = 10

# Create a points object
points = vtk.vtkPoints()

# Generate points for the grid
for i in range(num_points + 1):
    angle = np.pi * i / num_points
    for j in range(num_radial_divisions + 1):
        radius = inner_radius + j * (outer_radius - inner_radius) / num_radial_divisions
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        points.InsertNextPoint(x, y, 0.0)

# Create cells for the grid
cells = vtk.vtkCellArray()

for i in range(num_points):
    for j in range(num_radial_divisions):
        quad = vtk.vtkQuad()
        quad.GetPointIds().SetId(0, i * (num_radial_divisions + 1) + j)
        quad.GetPointIds().SetId(1, (i + 1) * (num_radial_divisions + 1) + j)
        quad.GetPointIds().SetId(2, (i + 1) * (num_radial_divisions + 1) + j + 1)
        quad.GetPointIds().SetId(3, i * (num_radial_divisions + 1) + j + 1)
        cells.InsertNextCell(quad)

# Create a vtkUnstructuredGrid
unstructured_grid = vtk.vtkUnstructuredGrid()
unstructured_grid.SetPoints(points)
unstructured_grid.SetCells(vtk.VTK_QUAD, cells)

# Write the unstructured grid to a .vtu file in ASCII format
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(r"tests\semicircle_grid.vtu")
writer.SetInputData(unstructured_grid)
writer.SetDataModeToAscii()  # Set the writer to ASCII mode
writer.Write()
