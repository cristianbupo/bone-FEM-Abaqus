import FreeCAD as App
import argparse
import json
import macros.FreeCAD.sliderWindow as sw
import geomdl2gmsh as g2g
import FreeCADGui as Gui
import os
import shutil
import gmsh


def load_json_files():
    with open('pickle/geom_vars.json', 'r') as file1, \
            open('pickle/mesh_vars.json', 'r') as file2, \
            open('pickle/load_vars.json', 'r') as file3:
        geom_vars = json.load(file1)
        mesh_vars = json.load(file2)
        load_vars = json.load(file3)
    return geom_vars, mesh_vars, load_vars


def filter_geom_vars(geom_vars):
    return {k: v for k, v in geom_vars.items() if not v.get('ignore', False)}


def add_geom_args(parser, geom_vars):
    for var_name, var_info in geom_vars.items():
        min_val, max_val = var_info["limits"]
        def_val = (min_val + max_val) / 2
        short_name = var_info["short_name"]
        parser.add_argument(f"-{short_name}", f"--{var_name}", type=lambda x,
                            vn=var_name: check_limits(geom_vars, x, vn),
                            default=def_val,
                            help=f"{var_name.replace('_', ' ').capitalize()} (limit: "
                            f"{min_val} to {max_val}), default: {def_val}")


def add_mesh_args(parser, mesh_vars):
    for var_name, var_info in mesh_vars.items():
        if var_name == 'number_elements':
            min_val = var_info["min_value"]
            def_val = var_info["default"]
            short_name = var_info["short_name"]
            parser.add_argument(f"-{short_name}", f"--{var_name}", type=lambda x,
                                vn=var_name: check_min_limit(mesh_vars, x, vn),
                                default=def_val,
                                help=f"{var_name.replace('_', ' ').capitalize()} "
                                f"(limit: {min_val}), default: {def_val}")


def add_load_args(parser, load_vars):
    for var_name, var_info in load_vars.items():
        min_val = var_info["min_value"]
        def_val = var_info["default"]
        short_name = var_info["short_name"]
        parser.add_argument(f"-{short_name}", f"--{var_name}", type=lambda x,
                            vn=var_name: check_min_limit(load_vars, x, vn),
                            default=def_val,
                            help=f"{var_name.replace('_', ' ').capitalize()} "
                            f"(limit: {min_val}), default: {def_val}")


def add_misc_args(parser):
    parser.add_argument('--showLoad', action='store_true', default=False,
                        help="Show the graph of the load distribution if specified.")
    parser.add_argument('--renderRaw', action='store_true', default=False,
                        help="Render the raw curves if specified.")
    parser.add_argument('--renderMesh', action='store_true', default=False,
                        help="Render the mesh if specified.")
    parser.add_argument('--renderArea', action='store_true', default=False,
                        help="Render the area if specified.")
    parser.add_argument('--renderLength', action='store_true', default=False,
                        help="Render the length if specified.")
    parser.add_argument('--runFltk', action='store_true', default=False,
                        help="Run gmsh with the generated script if specified.")
    parser.add_argument('--skipWrite', action='store_true', default=False,
                        help="Skip writing the generated script if specified.")
    parser.add_argument('--writeVTK', action='store_true', default=False,
                        help="Write the generated VTK file if specified.")
    parser.add_argument('--runMeshing', action='store_true', default=False,
                        help="Mesh the model")
    parser.add_argument('--pickleIn', action='store_true', default=False,
                        help="Pickle out data")
    parser.add_argument('--pickleOut', action='store_true', default=False,
                        help="Pickle out data")
    parser.add_argument('--deleteBackUp', action='store_true', default=False,
                        help="Delete the backup files")


def setup_sketch(args, geom_vars):
    initial_values = {var: getattr(args, var) for var in geom_vars}
    Gui.setupWithoutGUI()
    fileName = 'CADs/longBone.FCStd'
    doc = App.openDocument(fileName)
    sketch = doc.getObject('Sketch')

    for var, var_info in geom_vars.items():
        desired_value = getattr(args, var)
        min_value, max_value = var_info["limits"]
        sw.setConstraintValue(sketch, var, desired_value, min_value, max_value)

    for var in geom_vars:
        current_value = sketch.getDatum(var)
        if current_value != initial_values[var]:
            print(f"Value for {var} was not set correctly. Expected {initial_values[var]}, "
                  f"but got {current_value}.")
        else:
            print(f"Value for {var} was set correctly: {current_value}")

    return sketch


def process_gmsh(args, sketch):
    gmsh.initialize()
    curvesMesh, _, _ = g2g.processSketchNurbs(
        sketch,
        renderRaw=args.renderRaw,
        renderMesh=args.renderMesh,
        renderArea=args.renderArea,
        renderLength=args.renderLength
    )

    if args.deleteBackUp:
        backup_dir = os.path.join("currentBackUp")
        shutil.rmtree(backup_dir, ignore_errors=True)
        os.makedirs(backup_dir, exist_ok=True)

        currentLoadHistory_dir = os.path.join("currentLoadHistory")
        shutil.rmtree(currentLoadHistory_dir, ignore_errors=True)
        os.makedirs(currentLoadHistory_dir, exist_ok=True)

    g2g.container2gmsh(
        args.load_center,
        args.load_amplitude,
        args.load_radius,
        args.kOI,
        curvesMesh,
        inputPath='current',
        numberElements=args.number_elements,
        runFltk=args.runFltk,
        skipWrite=args.skipWrite,
        writeVTK=args.writeVTK
    )

    if args.pickleOut:
        backup_files("current", "currentBackUp")

    if args.pickleIn:
        backup_files("currentBackUp", "current")

    gmsh.finalize()


def backup_files(source_dir, backup_dir):
    os.makedirs(source_dir, exist_ok=True)
    os.makedirs(backup_dir, exist_ok=True)

    for file_name in os.listdir(source_dir):
        if os.path.isfile(os.path.join(source_dir, file_name)):
            shutil.copy(os.path.join(source_dir, file_name), backup_dir)

    print(f"Files have been copied to {backup_dir}")


def parserFunction():
    geom_vars, mesh_vars, load_vars = load_json_files()
    geom_vars = filter_geom_vars(geom_vars)

    parser = argparse.ArgumentParser(description="This script parses and checks variable "
                                                 "limits. Use the short names (e.g., -bw) "
                                                 "or long names (e.g., --bone_width) to "
                                                 "specify values within their respective "
                                                 "limits.")
    add_geom_args(parser, geom_vars)
    add_mesh_args(parser, mesh_vars)
    add_load_args(parser, load_vars)
    add_misc_args(parser)

    args = parser.parse_args()
    sketch = setup_sketch(args, geom_vars)
    process_gmsh(args, sketch)


def parse_numeric(value):
    try:
        # First, try parsing the value as an integer.
        int_value = int(value)
        return int_value
    except ValueError:
        # If it fails, try parsing it as a float.
        try:
            float_value = float(value)
            return float_value
        except ValueError:
            # If both conversions fail, the value is not numeric.
            return None


def check_limits(dictionary, value, variable_name):
    min_limit, max_limit = dictionary[variable_name]["limits"]
    value = parse_numeric(value)
    if value < min_limit or value > max_limit:
        raise argparse.ArgumentTypeError(f"{variable_name} must be between {min_limit} "
                                         f"and {max_limit}")
    return value


def check_min_limit(dictionary, value, variable_name):
    min_limit = dictionary[variable_name]["min_value"]
    value = parse_numeric(value)
    if value < min_limit:
        raise argparse.ArgumentTypeError(f"{variable_name} must be greater than or "
                                         f"equal to {min_limit}")
    return value


if __name__ == "__main__":
    parserFunction()