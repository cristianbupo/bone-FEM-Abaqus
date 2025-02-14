# Wait for analisis.log to be deleted and monitor the log file in the same window
Write-Output "Monitoring Abaqus log file..."

while ($true) {
    if (-not (Test-Path "current\analisis.log")) {
        Clear-Host
        Write-Output "analisis.log deleted. Waiting for it to appear again..."
        while (-not (Test-Path "current\analisis.log")) {
            Start-Sleep -Seconds 1
        }
        Write-Output "analisis.log appeared. Monitoring log file..."
    }
    Get-Content 'current\analisis.log' -Wait
    Start-Sleep -Seconds 1
}