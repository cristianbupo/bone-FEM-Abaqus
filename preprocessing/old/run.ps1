param (
    [switch]$debug
)

# Run Python script
Write-Output "Running Python script..."
try {
    & "C:\Users\crist\venvs\fcenv\Scripts\python.exe" "C:\Users\crist\git\bone-FEM-Abaqus\longBone.py"
    if ($LASTEXITCODE -ne 0) {
        throw "Python script failed with exit code $LASTEXITCODE"
    }
} catch {
    Write-Error "Failed to run Python script: $_"
    exit 1
}

# Change directory to Abaqus current folder
Set-Location "C:\Users\crist\git\bone-FEM-Abaqus\current"

# Determine which batch file to run

if ($debug) {
    Write-Output "Running Abaqus with debug.bat..."
    $batFile = ".\debug.bat"
} else {
    Write-Output "Running Abaqus with run.bat..."
    $batFile = ".\run.bat"
}

# Run the selected batch file
try {
    Start-Process -FilePath $batFile -Wait -NoNewWindow
    if ($LASTEXITCODE -ne 0) {
        throw "Batch file ($batFile) failed with exit code $LASTEXITCODE"
    }
} catch {
    Write-Error "Failed to run $batFile"
    exit 1
}

Set-Location ..
exit

