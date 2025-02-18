Clear-Host
Set-Location ".\current"
$batFile = ".\run.bat"
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
