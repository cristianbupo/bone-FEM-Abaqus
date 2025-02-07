import os
import sys
import FreeCAD as App

current_dir = os.path.abspath(os.path.dirname(__file__))
main_dir = os.path.abspath(os.path.join(current_dir, '../..'))
sys.path.append(main_dir)

import geomdl2gmsh as g2g

p = App.ParamGet("User parameter:BaseApp/Preferences/Macro")
macros_dir = p.GetString("MacroPath")
os.chdir(macros_dir)
os.chdir('../..')

fileName='CADs/longBone.FCStd'
doc = App.openDocument(fileName)
sketch = doc.getObject('Sketch')

curvesMesh, curvesArea, curvesLength = g2g.processSketchNurbs(sketch)

g2g.container2gmsh(curvesMesh, 0, 1, 1, runFltk=True, skipWrite=True, numberElements=10)

App.closeDocument(fileName)
