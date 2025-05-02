import os
import sys
import FreeCAD as App

p = App.ParamGet("User parameter:BaseApp/Preferences/Macro")
macros_dir = p.GetString("MacroPath")
main_dir = os.path.abspath(os.path.join(macros_dir, '../..'))
sys.path.append(main_dir)

import geomdl2gmsh as g2g

sketch = App.ActiveDocument.getObject('Sketch')

curvesMesh, _, _ = g2g.processSketchNurbs(sketch)

g2g.container2gmsh(curvesMesh, 0, 1, 1, runFltk=True, skipWrite=True, numberElements=10)
