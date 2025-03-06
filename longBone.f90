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
  end subroutine chdir

  subroutine run_commands_from_file(filename)
   character(len=*), intent(in) :: filename
   character(len=512) :: command, folder, lockFilePath, jobname
   integer :: unit, ios, status
   logical :: lock_exists

   open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
   if (ios /= 0) then
      print *, "Error opening file: ", filename
      return
   end if

   do
      read(unit, *, iostat=ios) folder, jobname
      if (ios /= 0) exit

      command = trim(folder) // "\\" // trim(jobname) // ".bat"
      lockFilePath = trim(folder) // "\\" // trim(jobname) // ".lck"

      ! Debug output
      print *, "Folder: ", trim(folder)
      print *, "Job name: ", trim(jobname)
      print *, "Command: ", trim(command)
      print *, "Lock file: ", trim(lockFilePath)

      ! Execute the command
      call chdir(trim(folder), err=status)
      if (status /= 0) then
         print *, "Error changing directory to: ", trim(folder), " status code: ", status
         cycle
      end if

      call execute_command_line(trim(command), wait=.true., exitstat=status)
      if (status /= 0) then
         print *, "Error executing command: ", trim(command), " status code: ", status
      else
         print *, "Command executed successfully: ", trim(jobname) // ".lck"
      end if

      ! Wait for the lock file to appear
      print *, "Waiting for lock file to appear: ", trim(jobname) // ".lck"
      do
         inquire(file=lockFilePath, exist=lock_exists)
         if (lock_exists) exit
         call sleep(1)  ! Wait for 1 second before checking again
      end do

      ! Wait for the lock file to disappear
      print *, "Waiting for lock file to disappear: ", trim(jobname) // ".lck"
      do
         inquire(file=lockFilePath, exist=lock_exists)
         if (.not. lock_exists) exit
         call sleep(1)  ! Wait for 1 second before checking again
      end do

   end do

   close(unit)
  end subroutine run_commands_from_file

  subroutine run_python_file()
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
   pythonScriptPath = "C:/Users/crist/git/bone-FEM-Abaqus/longBone.py"
   pythonCommand = trim(pythonVenvPath) // " " // trim(pythonScriptPath)

   inputPath = "C:/Users/crist/OneDrive/Documents/results/singleAnalysis/"
   command = "run.bat"

   write(*,*) trim(pythonCommand)
   
   ! Run python script
   call execute_command_line(pythonCommand, wait=.true., exitstat=status)

   ! Check the exit status
   if (status /= 0) then
      write(*,*) "Error executing run.bat, status code: ", status
   else
      write(*,*) "Batch file executed successfully!"
   end if

  end subroutine run_python_file

end module chdir_mod


program run_bat

   use chdir_mod

   ! call run_python_file()
   ! Run commands from runCommands.bat
   call run_commands_from_file('runCommands.bat')

end program run_bat