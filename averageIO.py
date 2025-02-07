import os
import vtk
import glob
import numpy as np


def read_vtu_file(file_path):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_path)
    reader.Update()
    return reader.GetOutput()


def extract_data_array(vtk_data, array_name):
    # data_array = vtk_data.GetPointData().GetArray(array_name)
    data_array = vtk_data.GetCellData().GetArray(array_name)
    if data_array is None:
        raise ValueError(f"Data array '{array_name}' not found in the file.")
    num_tuples = data_array.GetNumberOfTuples()
    values = [data_array.GetValue(i) for i in range(num_tuples)]
    return np.array(values)


def write_to_file(data, output_file):
    with open(output_file, 'w') as f:
        for value in data:
            f.write(f"{value}\n")


def write_vtu_file(data, output_file, reference_vtk_data):

    # Create a new vtkUnstructuredGrid
    output_data = vtk.vtkUnstructuredGrid()

    # Copy points and cells from the reference data
    output_data.SetPoints(reference_vtk_data.GetPoints())
    output_data.SetCells(reference_vtk_data.GetCellTypesArray(), 
                         reference_vtk_data.GetCellLocationsArray(), 
                         reference_vtk_data.GetCells())

    # Preserve the Region array
    region_array = reference_vtk_data.GetCellData().GetArray("Region")
    if region_array:
        output_data.GetCellData().AddArray(region_array)

    # Create a new vtkDoubleArray to store the input data
    io_data_array = vtk.vtkDoubleArray()
    io_data_array.SetName("IO")
    io_data_array.SetNumberOfComponents(1)
    io_data_array.SetNumberOfTuples(len(data))
    for i, element in enumerate(data):
        io_data_array.SetValue(i, element)

    # Add the IO data array to the CellData of the output grid
    output_data.GetCellData().AddArray(io_data_array)

    # Write the output grid to a .vtu file in ASCII format
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(output_data)
    writer.Write()


def calculate_average(folder_path):
    # Path to the directory containing the .vtu files
    # folder_path = 'currentLoadHistory'

    # Path to the directory containing the .vtu files
    analysis_files = glob.glob(os.path.join(folder_path, 'analisis*.vtu'))
    all_data_arrays = []

    print(os.path.join(folder_path, 'analisis*.vtu'))
    print('Analysis files found:')
    for file_path in analysis_files:
        print(file_path)

    for file_path in analysis_files:
        vtk_data = read_vtu_file(file_path)
        data_array = extract_data_array(vtk_data, "IO")
        all_data_arrays.append(data_array)

    if not all_data_arrays:
        print("No data arrays found. Exiting.")
        return

    all_data_arrays = np.array(all_data_arrays)
    average_data = np.mean(all_data_arrays, axis=0)

    # Output file to write the average data
    output_file = 'averageOI.vtu'
    output_path = os.path.join(folder_path, output_file)
    write_vtu_file(average_data, output_path, vtk_data)
    print(f"Average data written to {output_path} and {output_path}")


if __name__ == "__main__":
    calculate_average('currentLoadHistory')
