from math import pi
import numpy as np
from scipy import integrate

def setConstraintValue(sketch, constraint_name, rounded_val, min_val, max_val):
    if isinstance(min_val, str):
        min_val = sketch.getDatum(min_val).Value

    if isinstance(max_val, str):
        max_val = sketch.getDatum(max_val).Value

    changeConstraintValue(sketch, constraint_name, rounded_val)


def changeConstraintValue(sketch, var, desired_val, max_outer_attempts=10, max_inner_attempts=10):
    angle_factor = 1

    unitType = sketch.getDatum(var).Unit.Type
    if unitType == "Angle":
        angle_factor = pi / 180

    outer_attempts = 0

    while outer_attempts < max_outer_attempts:
        try:
            sketch.setDatum(var, desired_val * angle_factor) # Set in radians
            break  # If successful, break the loop
        except Exception:
            current_val = sketch.getDatum(var).Value # Get in degrees
            k = 2
            inner_attempts = 0

            while inner_attempts < max_inner_attempts:
                try:
                    try_val = ((k - 1) * current_val + desired_val) / k
                    sketch.setDatum(var, try_val * angle_factor) # Set in radians
                    break  # If successful, break the inner loop
                except Exception:
                    k *= 2
                    inner_attempts += 1

            outer_attempts += 1


def limit_val(rounded_val, min_val, max_val):
    min_space = 0
    if rounded_val < min_val + min_space:
        rounded_val = min_val + min_space
    elif rounded_val > max_val - min_space:
        rounded_val = max_val - min_space

    return rounded_val


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

    integral = integrate.simpson(y=integer, x=u_values)/2
    return integral

def findLenght(curve, startTag = 0, endTag = 1.0, nPoints = 1000):

    # Evaluate derivatives at all points
    u_values = np.linspace(startTag, endTag, nPoints)
    derivatives = np.array([curve.derivatives(u, order=1) for u in u_values])

    # Extract x and y coordinates of the derivatives
    x_der = derivatives[:, 1, 0]
    y_der = derivatives[:, 1, 1]

    lenght = integrate.simpson(y=np.hypot(x_der, y_der), x=u_values)

    return lenght

def findArea(curvesArea):
    return findIntegral(curvesArea[3]) - findIntegral(curvesArea[1])