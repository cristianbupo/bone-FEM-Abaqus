      module debug_utils
         implicit none

         ! Global variable to control printing
         logical :: debugMode = .false.

      contains

      subroutine printStepInfo(lop, kstep, kinc, message)
         implicit none
         integer, intent(in) :: lop, kstep, kinc
         character(len=13), intent(in) :: message

         ! Check if debugMode is enabled
         if (debugMode) then
            print *, message, " LOP=", lop, "KSTEP=", kstep, "KINC=", kinc
            print *, ''
         end if
      end subroutine printStepInfo

      end module debug_utils