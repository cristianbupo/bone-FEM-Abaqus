import vtk
import numpy as np

def createAndWritePolydata(points, normals, magnitudes):

    # Create vtkPolyData for normals visualization
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.GetPointData().SetNormals(normals)
    polydata.GetPointData().AddArray(magnitudes)

    # Write vtkPolyData to a file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName("tests/semicircle_normals.vtp")
    writer.SetInputData(polydata)
    writer.SetDataModeToAscii()  # Set the writer to ASCII mode
    writer.Write()

# Parameters for the semicircles
outer_radius = 1.0
inner_radius = 0.5
num_points = 50
num_radial_divisions = 10

# Calculate midpoints and normals for the outer surface
midpoints = vtk.vtkPoints()
global_normals = vtk.vtkDoubleArray()
load_magnitudes = vtk.vtkDoubleArray()
global_normals.SetNumberOfComponents(3)  # Normal vectors are 3D
global_normals.SetName("Load Orientation")
load_magnitudes.SetName("Load Magnitude")

for i in range(num_points):
    # Calculate midpoint between two consecutive points on the outer surface
    angle1 = np.pi * i / num_points
    angle2 = np.pi * (i + 1) / num_points
    midpoint_x = outer_radius * np.cos((angle1 + angle2) / 2)
    midpoint_y = outer_radius * np.sin((angle1 + angle2) / 2)
    midpoints.InsertNextPoint(midpoint_x, midpoint_y, 0.0)
    load_magnitudes.InsertNextValue(1 - (2 * i/(num_points-1) - 1) ** 2)

    # Normal vector (radially outward)
    normal = [np.cos((angle1 + angle2) / 2),
              np.sin((angle1 + angle2) / 2), 0.0]
    global_normals.InsertNextTuple(normal)

createAndWritePolydata(midpoints, global_normals, load_magnitudes)
