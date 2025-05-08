@echo off
set SOURCE="../"
set DESTINATION="gdtr:abaqusFortranBackUp"
rclone sync %SOURCE% %DESTINATION% --progress