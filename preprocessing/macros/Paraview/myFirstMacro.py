
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
defaultvtu = XMLUnstructuredGridReader(registrationName='default.vtu', FileName=['D:\\bone-FEM\\results\\analisis\\default.vtu'])

# Properties modified on defaultvtu
defaultvtu.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
defaultvtuDisplay = Show(defaultvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
defaultvtuDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

#changing interaction mode based on data extents
renderView1.CameraPosition = [6.64500121594358e-11, -1.5417170317728, 13.12049588712224]
renderView1.CameraFocalPoint = [6.64500121594358e-11, -1.5417170317728, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(defaultvtuDisplay, ('CELLS', 'S_Inv'))

# rescale color and/or opacity maps used to include current data range
defaultvtuDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
defaultvtuDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'S_Inv'
s_InvLUT = GetColorTransferFunction('S_Inv')

# get opacity transfer function/opacity map for 'S_Inv'
s_InvPWF = GetOpacityTransferFunction('S_Inv')

# get 2D transfer function for 'S_Inv'
s_InvTF2D = GetTransferFunction2D('S_Inv')

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=defaultvtu)

# set active source
SetActiveSource(cellDatatoPointData1)

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Surface'

# hide data in view
Hide(defaultvtu, renderView1)

# set scalar coloring
ColorBy(cellDatatoPointData1Display, ('POINTS', 'U', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
cellDatatoPointData1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
cellDatatoPointData1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'U'
uLUT = GetColorTransferFunction('U')

# get opacity transfer function/opacity map for 'U'
uPWF = GetOpacityTransferFunction('U')

# get 2D transfer function for 'U'
uTF2D = GetTransferFunction2D('U')

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
uLUT.ApplyPreset('Rainbow Desaturated', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
uLUT.ApplyPreset('Rainbow Desaturated', True)

# set scalar coloring
ColorBy(cellDatatoPointData1Display, ('POINTS', 'S_Inv'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(uLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
cellDatatoPointData1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
cellDatatoPointData1Display.SetScalarBarVisibility(renderView1, True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
s_InvLUT.ApplyPreset('Rainbow Desaturated', True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
s_InvLUT.ApplyPreset('Rainbow Desaturated', True)

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'UnstructuredGridRepresentation')

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

# hide data in view
Hide(defaultvtu, renderView1)

# show color bar/color legend
cellDatatoPointData1Display.SetScalarBarVisibility(renderView1, True)

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(556, 471)


#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [6.64500121594358e-11, -1.5417170317728, 8.17292696945905]
renderView1.CameraFocalPoint = [6.64500121594358e-11, -1.5417170317728, 0.0]
renderView1.CameraParallelScale = 2.11530915392803
