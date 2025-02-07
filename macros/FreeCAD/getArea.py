import sys
import os
import FreeCAD as App
import numpy as np
from scipy import integrate
from geomdl import operations

current_dir = os.path.abspath(os.path.dirname(__file__))
main_dir = os.path.abspath(os.path.join(current_dir, '../..'))
sys.path.append(main_dir)

import geomdl2gmsh as gn

def findIntegral(curve, nPoints=1000, invert=False):
    startTag = 0.0
    endTag = 1.0
    startPoint = curve.ctrlpts[0]
    endPoint = curve.ctrlpts[-1]
    if invert:
        endTag = 0.0
        startTag = 1.0
        startPoint = curve.ctrlpts[-1]
        endPoint = curve.ctrlpts[0]

    print("Start point:", startPoint)
    print("End point:", endPoint)

    # Evaluate derivatives at all points
    u_values = np.linspace(startTag, endTag, nPoints)
    derivatives = np.array([curve.derivatives(u, order=1) for u in u_values])

    # Extract x and y coordinates of the derivatives
    x_coords = derivatives[:, 0, 0]
    y_coords = derivatives[:, 0, 1]

    # Extract x and y coordinates of the derivatives
    x_der = derivatives[:, 1, 0]
    y_der = derivatives[:, 1, 1]

    integer = x_coords * y_der - y_coords * x_der

    integral = integrate.simps(integer, u_values)/2
    return integral

def findLenght(curve, startTag = 0, endTag = 1.0, nPoints = 1000):

    # Evaluate derivatives at all points
    u_values = np.linspace(startTag, endTag, nPoints)
    derivatives = np.array([curve.derivatives(u, order=1) for u in u_values])

    # Extract x and y coordinates of the derivatives
    x_der = derivatives[:, 1, 0]
    y_der = derivatives[:, 1, 1]

    lenght = integrate.simps(np.hypot(x_der + y_der), u_values)

    return lenght

if __name__ == "__main__":

    p=App.ParamGet("User parameter:BaseApp/Preferences/Macro")
    macros_dir=p.GetString("MacroPath")
    os.chdir(macros_dir)

    sketch = App.ActiveDocument.getObject('Sketch')
    _, curvesArea, curvesLenght = gn.processSketchNurbs(sketch)

    integral1 = findIntegral(curvesArea[1])
    integral2 = findIntegral(curvesArea[3])
    print("Integral 1:", integral1)
    print("Integral 2:", integral2)
    print("Area:", - integral1 + integral2)
    print("Length:", findLenght(curvesLenght[1]))
    print("Coarser length:", operations.length_curve(curvesLenght[1]))
