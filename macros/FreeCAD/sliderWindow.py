
import os
from numpy import pi
import FreeCAD as App
import json

from PySide6.QtWidgets import *
from PySide6.QtCore import *
from PySide6.QtWidgets import QVBoxLayout, QWidget

n = 5
min_space = 0  # Minimum space to add when reaching min or max
default_decimals = 1


class AnimationSliderWindow(QWidget):
    num_sliders = 0  # Class attribute to keep track of the number of sliders

    def __init__(self, sketch, parent=None):
        super().__init__(parent)
        self.sketch = sketch
        self.setWindowFlags(Qt.WindowStaysOnTopHint)
        self.setWindowTitle("Animation Slider")

        layout = QVBoxLayout()

        # Create a slider for each item in the dictionary
        for var_name, var_info in variables.items():
            min_val, max_val = var_info["limits"]
            self.createSlider(layout, var_name, min_val, max_val)

        self.setLayout(layout)

        # Calculate the window height based on the number of sliders
        self.calculateWindowHeight()
        self.show()

    def createSlider(self, layout, constraint_name, min_val, max_val):
        label_text = constraint_name
        label = QLabel(label_text)
        layout.addWidget(label)

        factor = 10**n

        initial_slider_val = int(self.sketch.getDatum(constraint_name) * factor)

        slider = QSlider(Qt.Horizontal)
        slider.setObjectName(f"slider_{constraint_name}")
        slider.setMinimum(min_val * factor)
        slider.setMaximum(max_val * factor)

        # Create increment and decrement buttons
        increment_button = QPushButton("+")
        decrement_button = QPushButton("-")

        # Set a fixed size for the buttons
        increment_button.setFixedSize(30, 30)
        decrement_button.setFixedSize(30, 30)

        # Connect the buttons' clicked signals to slots that increment and decrement
        # the slider value
        increment_button.clicked.connect(lambda: slider.setValue(slider.value()
                                                                 + factor/10))
        decrement_button.clicked.connect(lambda: slider.setValue(slider.value()
                                                                 - factor/10))

        slider.setSingleStep(0.5)
        slider.setValue(initial_slider_val)
        slider.setFixedWidth(200)  # Set the fixed width to 400 pixels

        value_input_box = QLineEdit()
        value_input_box.setFixedWidth(60)
        value_input_box.setText(f"{initial_slider_val / factor:.{n}f}")

        decimals_input_box = QSpinBox()
        decimals_input_box.setMinimum(-5)
        decimals_input_box.setMaximum(5)
        decimals_input_box.setValue(default_decimals)

        rounded_val_box = QLabel()

        min_max_step_n_label = QLabel(f"Min: {min_val:.{n}f}  Max: {max_val:.{n}f}")

        slider_layout = QHBoxLayout()
        slider_layout.addWidget(slider)
        slider_layout.addWidget(decrement_button)
        slider_layout.addWidget(increment_button)
        slider_layout.addWidget(value_input_box)
        slider_layout.addWidget(decimals_input_box)
        slider_layout.addWidget(rounded_val_box)
        slider_layout.addWidget(min_max_step_n_label)

        layout.addLayout(slider_layout)

        slider.valueChanged.connect(lambda value, cn=constraint_name,
                                    vib=value_input_box, dbox=decimals_input_box,
                                    rbox=rounded_val_box:
                                    self.sliderMoved(value, cn, vib, dbox, rbox))

        min_slider_val = slider.minimum() / factor
        max_slider_val = slider.maximum() / factor
        value_input_box.editingFinished.connect(lambda vib=value_input_box,
                                                slider=slider, cn=constraint_name:
                                                self.inputBoxEdited(vib, slider, cn,
                                                                    min_slider_val,
                                                                    max_slider_val))

    def sliderMoved(self, value, constraint_name, value_input_box, decimals_input_box,
                    rounded_val_box):
        factor = 10**n
        real_val = float(value) / factor
        decimals = decimals_input_box.value()
        rounded_val = round(real_val, decimals)

        # Get the minimum and maximum values of the slider
        slider = value_input_box.parentWidget().findChild(QSlider,
                                                          f"slider_{constraint_name}")
        min_slider_val = slider.minimum() / factor
        max_slider_val = slider.maximum() / factor

        # Check if the rounded value is within 1e-4 of the minimum or maximum
        rounded_val = limit_val(rounded_val, min_slider_val, max_slider_val)

        value_input_box.setText(f"{real_val}")
        rounded_val_box.setText(f"Rounded: {rounded_val}")

        setConstraintValue(self.sketch, constraint_name, rounded_val,
                           min_slider_val, max_slider_val)

    def inputBoxEdited(self, value_input_box, slider, constraint_name, min_slider_val,
                       max_slider_val):
        try:
            new_val = float(value_input_box.text())
            factor = 10**n
            new_slider_val = int(new_val * factor)
            slider.setValue(new_slider_val)

            # Initialize rounded_val before using it in limit_val
            rounded_val = new_val

            # Check if the rounded value is within 1e-4 of the minimum or maximum
            rounded_val = limit_val(rounded_val, min_slider_val, max_slider_val)

            setConstraintValue(self.sketch, constraint_name, rounded_val,
                               min_slider_val, max_slider_val)
        except ValueError:
            # Handle invalid input (non-numeric)
            pass

    def calculateWindowHeight(self):
        # Calculate the height of the window based on the number of sliders
        window_height = 100 + AnimationSliderWindow.num_sliders * 50
        # Adjust the base height (100) and height per slider (50) as needed
        self.resize(600, window_height)
        # Adjust the width and set the calculated height


def setConstraintValue(sketch, constraint_name, rounded_val, min_val,
                       max_val):
    if isinstance(min_val, str):
        min_val = sketch.getDatum(min_val).Value

    if isinstance(max_val, str):
        max_val = sketch.getDatum(max_val).Value

    unitType = sketch.getDatum(constraint_name).Unit.Type

    if unitType == "Angle":
        rounded_val = rounded_val*pi/180

    changeConstraintValue(sketch, constraint_name, rounded_val)


def changeConstraintValue(sketch, var, desired_val, max_outer_attempts=5,
                          max_inner_attempts=10, angleFactor=1.0):
    outer_attempts = 0
    inner_attempts = 0
    while outer_attempts < max_outer_attempts and inner_attempts < max_inner_attempts:
        try:
            sketch.setDatum(var, desired_val)
            break  # If successful, break the loop
        except Exception as e:
            current_val = sketch.getDatum(var).Value*angleFactor
            print(f"An error occurred: {e}")
            k = 2
            inner_attempts = 0
            while inner_attempts < max_inner_attempts:
                try:
                    try_val = ((k-1)*current_val + desired_val) / k
                    print(f"Trying change {var} value with {try_val}...")
                    sketch.setDatum(var, try_val)
                    break
                except Exception as ex:  # Catching a more general exception
                    print(f"An error occurred: {ex}")  # Logging the error
                    k = 2*k
                    inner_attempts += 1
            outer_attempts += 1

def limit_val(rounded_val, min_val, max_val):
    if rounded_val < min_val + min_space:
        rounded_val = min_val + min_space
    elif rounded_val > max_val - min_space:
        rounded_val = max_val - min_space

    return rounded_val


if __name__ == "__main__":
    global_sketch = App.ActiveDocument.getObject('Sketch')
    p = App.ParamGet("User parameter:BaseApp/Preferences/Macro")
    macros_dir = p.GetString("MacroPath")
    os.chdir(macros_dir)

    # Load dictionary from the JSON file
    with open('../../combined_vars.json', 'r') as f:
        variables = json.load(f)['geom_vars']

    # Remove entries where "ignore" is true
    variables = {k: v for k, v in variables.items() if not v.get('ignore', False)}

    AnimationSliderWindow(global_sketch)
