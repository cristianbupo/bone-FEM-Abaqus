param (
    [string]$fileName = "current\analisis.log"
)

# Extract the actual file name from the full path
# $actualFileName = [System.IO.Path]::GetFileName($fileName)

# Define a method to set the console title
Add-Type -TypeDefinition @"
using System;
using System.Runtime.InteropServices;

public class ConsoleTitle {
    [DllImport("kernel32.dll", SetLastError = true)]
    public static extern bool SetConsoleTitle(string lpConsoleTitle);
}
"@

# Set the console title to include the file name
[ConsoleTitle]::SetConsoleTitle("$fileName")

# Wait for the specified file to be deleted and monitor the log file in the same window
Write-Output "Monitoring $fileName..."

while ($true) {
    if (-not (Test-Path $fileName)) {
        Clear-Host
        Write-Output "$fileName deleted. Waiting for it to appear again..."
        while (-not (Test-Path $fileName)) {
            Start-Sleep -Seconds 1
        }
        Write-Output "$fileName appeared. Monitoring log file..."
    }
    Get-Content $fileName -Wait
    Start-Sleep -Seconds 1
}