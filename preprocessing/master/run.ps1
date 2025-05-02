#Get the current file directory
$scriptPath = Split-Path -Parent $MyInvocation.MyCommand.Definition

# Define the path to the runCommands.txt file, on the scriptPath directory
$filename = "$scriptPath\runCommands.txt"

# Check if the file exists
if (-Not (Test-Path $filename)) {
    Write-Host "Error opening file: $filename"
    exit
}

# Read the file line by line
Get-Content $filename | ForEach-Object {
    # Split the line into folder and jobname
    $line = $_.Trim()
    $parts = $line -split ","
    $folder = $parts[0]
    $jobname = $parts[1]

    # Construct the command and lock file path
    $command = "$folder\$jobname.bat"
    $lockFilePath = "$folder\$jobname.lck"

    # Debug output
    Write-Host "Folder: $folder"
    Write-Host "Job name: $jobname"
    Write-Host "Command: $command"
    Write-Host "Lock file: $lockFilePath"

    # Change the directory to the specified folder
    Set-Location -Path $folder

    # Execute the command
    & $command
    $status = $LASTEXITCODE

    if ($status -ne 0) {
        Write-Host "Error executing command: $command, status code: $status"
    } else {
        Write-Host "Command executed successfully: $jobname.lck"
    }

    # Wait for the lock file to appear
    Write-Host "Waiting for lock file to appear: $jobname.lck"
    while (-Not (Test-Path $lockFilePath)) {
        Start-Sleep -Seconds 1
    }

    # Wait for the lock file to disappear
    Write-Host "Waiting for lock file to disappear: $jobname.lck"
    while (Test-Path $lockFilePath) {
        Start-Sleep -Seconds 1
    }
}