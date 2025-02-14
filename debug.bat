@echo off
c:\Users\crist\venvs\fcenv\Scripts\python.exe c:\Users\crist\git\bone-FEM-Abaqus\longBone.py
cd C:\Users\crist\git\bone-FEM-Abaqus\current
call debug.bat
powershell -Command "Get-Content .\analisis.log -Wait"
cd ..