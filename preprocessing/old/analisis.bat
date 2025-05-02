echo Running Abaqus

call runAbaqus.bat

echo Waiting for Abaqus to finish
:loop
IF EXIST current\analisis.lck (
    :innerloop
    IF NOT EXIST current\analisis.lck (
        echo Lock file removed
        GOTO end
    )
    TIMEOUT /T 1
    GOTO innerloop
)
TIMEOUT /T 1
GOTO loop
:end

echo Abaqus finished