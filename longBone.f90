module chdir_mod

  implicit none

  interface
    integer function c_chdir(path) bind(C,name="chdir")
      use iso_c_binding
      character(kind=c_char) :: path(*)
    end function
  end interface

contains

  subroutine chdir(path, err)
    use iso_c_binding
    character(*) :: path
    integer, optional, intent(out) :: err
    integer :: loc_err

    loc_err =  c_chdir(path//c_null_char)

    if (present(err)) err = loc_err
  end subroutine
end module chdir_mod


program run_bat

    implicit none

    integer :: status

    ! Define commands

    character(len=512) :: pythonVenvPath
    character(len=512) :: pythonScriptPath
    character(len=512) :: pythonCommand

    ! Define the full path to run.bat
    character(len=512) :: inputPath
    character(len=512) :: command

    ! Get the current working directory
    CHARACTER(len=255) :: cwd


    pythonVenvPath = "C:/Users/crist/venvs/fcenv/Scripts/python.exe"
    pythonScriptPath = "c:/Users/crist/git/bone-FEM-Abaqus/helloWorld.py"
    pythonCommand = pythonVenvPath // " " // pythonScriptPath

    inputPath = "C:\Users\crist\OneDrive\Documents\results\singleAnalysis\"
    command = "run.bat"


    ! Run python script
    ! call execute_command_line(pythonCommand, exitstat=status)

    ! Change the current working directory
    call chdir(inputPath)

    ! Get the current working directory
    CALL getcwd(cwd)

    ! Print the current working directory
    WRITE(*,*) TRIM(cwd)

    ! Execute the command
    call execute_command_line(command, exitstat=status)

    ! Check the exit status
    if (status /= 0) then
        print *, "Error executing run.bat, status code: ", status
    else
        print *, "Batch file executed successfully!"
    end if

end program run_bat
