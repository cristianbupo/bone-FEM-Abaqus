#### INPUTS:
input0 = self.GetInput()
output = self.GetOutput()

# Get point coordinates
coords = input0.GetPoints().GetData()

# Get existing field arrays (example: 'Velocity')
velocities = input0.GetPointData().GetArray('Velocity')

# Prepare new array
import vtk
new_array = vtk.vtkDoubleArray()
new_array.SetName('Velocity_Magnitude_Squared')
new_array.SetNumberOfComponents(1)
new_array.SetNumberOfTuples(coords.GetNumberOfTuples())

# Calculate: u^2 + v^2
for i in range(coords.GetNumberOfTuples()):
    v = velocities.GetTuple3(i)
    mag_sq = v[0]**2 + v[1]**2  # Just u^2 + v^2 (no z component)
    new_array.SetValue(i, mag_sq)

# Attach new field
output.GetPointData().AddArray(new_array)
