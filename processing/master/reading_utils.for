      module reading_utils
      implicit none
   ! Variables a ser leidas del archivo de entrada
      integer :: NUMNODE, NELEMS, dim, nnod
      integer :: nLoads
      integer :: nResNod, nResElem
      integer :: maxNElementLoads
      integer :: numProps, numMats
      integer :: axi, tipo_def
      real*8 :: stdWeight
      real*8 :: kOI
      integer :: ndofdiff
      integer, allocatable :: listNElementLoads(:,:)

   ! Otras variables globales
      real*8 :: CMIThreshold,CMIavg,CMIstd,CMImax

   ! Arreglos cuyos tamaños dependen de las variables leidas del archivo de entrada
      real*8, allocatable :: elementLoads(:,:) 
      integer, allocatable :: elementFaces(:,:)
      real*8, allocatable :: nodes(:,:)
      integer, allocatable :: conectividades(:,:)
      integer, allocatable :: adyacencias(:,:)
      integer, allocatable :: grupoFisico(:,:)
      integer, allocatable :: grupoFisicoN(:,:)
      logical, allocatable :: nodoBorde(:,:)
      logical, allocatable :: elementoBorde(:,:)
      real*8, allocatable :: propiedades(:,:)
      real*8, allocatable :: resNod(:,:), resElem(:,:)
      real*8, allocatable :: cumulativeResNod(:,:), cumulativeResElem(:,:)
      real*8, allocatable :: CMICriteria(:)

      contains

      subroutine allocate_arrays()
         implicit none

      ! Allocate arrays based on the retrieved integer values
         allocate(elementLoads(maxNElementLoads, 1))
         allocate(elementFaces(maxNElementLoads, 2))
         allocate(nodes(NUMNODE, dim))
         allocate(conectividades(NELEMS, nnod+1))
         allocate(adyacencias(NELEMS, nnod+1))
         allocate(grupoFisico(NELEMS, 2))
         allocate(grupoFisicoN(NUMNODE, 2))
         allocate(nodoBorde(NUMNODE,1))
         allocate(elementoBorde(NELEMS,1))
         allocate(propiedades(numMats, numProps))
         allocate(resNod(NUMNODE, nResNod))
         allocate(resElem(NELEMS, nResElem))
         allocate(cumulativeResNod(NUMNODE, nResNod))
         allocate(cumulativeResElem(NELEMS, nResElem))
         allocate(CMICriteria(NELEMS))

      end subroutine allocate_arrays

      subroutine initialize_arrays()
         implicit none
      ! Initialize arrays to default values
         print *, ''
         print *, 'Initializing arrays...'

         nodes = 0.d0
         conectividades = 0
         adyacencias = 0
         grupoFisico = 0
         grupoFisicoN = 0
         nodoBorde = .false.
         elementoBorde = .false.
         propiedades = 0.d0
      end subroutine initialize_arrays

      subroutine initialize_results()
         implicit none
         print *, '' 
         print *, 'Initializing results...'
         
         resNod = 0.d0
         resElem = 0.d0
         cumulativeResNod = 0.d0
         cumulativeResElem = 0.d0
         CMICriteria = 0.d0
         elementLoads = 0.d0
         elementFaces = 0
      end subroutine initialize_results

      subroutine get_global_variables(jobdir,file)
         implicit none
      ! Input file containing the variables
         character(len=*) :: jobdir
         character(len=*) :: file
         character(len=256) :: filename
         character(len=256) :: cargasLine
         integer :: i


         filename=trim(jobdir)//"\"//trim(file)

      ! Call the subroutine for each variable and print the results
         call get_variable_from_file(filename, "NUMNODE", .false., NUMNODE)
         print *, "NUMNODE =", NUMNODE

         call get_variable_from_file(filename, "NELEMS", .false., NELEMS)
         print *, "NELEMS =", NELEMS

         call get_variable_from_file(filename, "dim", .false., dim)
         print *, "dim =", dim

         call get_variable_from_file(filename, "nnod", .false., nnod)
         print *, "nnod =", nnod

         call get_variable_from_file(filename, "nResNod", .false., nResNod)
         print *, "nResNod =", nResNod

         call get_variable_from_file(filename, "nResElem", .false., nResElem)
         print *, "nResElem =", nResElem

         call get_variable_from_file(filename, "numProps", .false., numProps)
         print *, "numProps =", numProps

         call get_variable_from_file(filename, "numMats", .false., numMats)
         print *, "numMats =", numMats

         call get_variable_from_file(filename, "axi", .false., axi)
         print *, "axi =", axi

         call get_variable_from_file(filename, "tipo_def", .false., tipo_def)
         print *, "tipo_def =", tipo_def

         call get_variable_from_file(filename, "stdWeight", .true., real_value=stdWeight)
         print *, "stdWeight =", stdWeight

         call get_variable_from_file(filename, "kOI", .true., real_value=kOI)
         print *, "kOI =", kOI

         call get_variable_from_file(filename, "ndofdiff", .false., ndofdiff)
         print *, "ndofdiff =", ndofdiff

         if (dim == 2) then
            filename=trim(jobdir)//"\"//trim('carga.inp')
            call readNLoads(filename,cargasLine,nLoads)
            print *, "nLoads =", nLoads
            print *, "Numero de cargas =", cargasLine
         end if

         allocate(listNElementLoads(nLoads, 1))

         ! Call the readVLoads subroutine to fill the vector
         call readVLoads(cargasLine, nLoads, listNElementLoads)

         ! Print the parsed integers
         maxNElementLoads = 0
         do i = 1, nLoads
             print *, listNElementLoads(i, 1)
             if (listNElementLoads(i, 1) > maxNElementLoads) then
                 maxNElementLoads = listNElementLoads(i, 1)
             end if
         end do
         print *, "maxNElementLoads =", maxNElementLoads
      end
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------Searstr----------- 
!     Funcion que devuelve un valor logico que permite saber si la lectura es 
!     adecuadad.
! 
!     Busca Str en la unidad Lu y devuelve .true.// si lo encuentra. Si no
!     lo encuentra devuelve un .false.
!
!-----------------------------------------------------------------------------------------
      function Searstr(Lu, Str)
!
      integer            :: Lu, j, L1, ios
      character (len=*)  :: Str
      character (len=80) :: Rec
      logical            :: Salida, Searstr
!
      Searstr = .false.
      Salida  = .true.
      L1 = len(Str)
      j       = 1
!
      do while (Salida .and. (j <= 2))
         read (Lu, '(a80)', iostat=ios) Rec
         if (ios /= 0) then
            print *, 'Error reading file in Searstr'
            exit
         end if
         if (Rec(1:L1) == (Str)) then
            Salida  = .false.
            Searstr = .true.
         else if (Rec(1:4) .eq. '*END') then
            rewind (Lu)! Rebobina el archivo
            j = j + 1
         end if
      enddo
!
      if (.not. Searstr) then
         print *, '###..Error en la funcion Searstr(Unidad=)', Lu, Str
         stop
      end if
!
      return
      end function Searstr
!------------------------------------------------------------------------------------------
      function count_char_in_string(str, char) result(count)
         implicit none
         character(len=*), intent(in) :: str
         character, intent(in) :: char
         integer :: count, i
       
         count = 0
         do i = 1, len_trim(str)
           if (str(i:i) == char) count = count + 1
         end do
       end function
!-------------------------------------------------------------------------------------------
      subroutine get_variable_from_file(filename, var_name, is_real, int_value, real_value)
         implicit none
         character(len=*), intent(in) :: filename! Name of the input file
         character(len=*), intent(in) :: var_name! Name of the variable to retrieve
         logical, intent(in) :: is_real          ! Flag to indicate if the variable is real*8
         integer, intent(out), optional :: int_value ! Value of the variable (integer)
         real*8, intent(out), optional :: real_value ! Value of the variable (real*8)
         character(len=256) :: line              ! Line read from the file
         character(len=256) :: key               ! Key extracted from the line
         character(len=256) :: value_str         ! Value as a string
         integer :: ios, pos, start_pos          ! I/O status, position of '=', and start position
         logical :: found                        ! Flag to indicate if variable is found
         integer :: unit_number                  ! Unit number for the file

         found = .false.
         unit_number = 11  ! Arbitrary unit number for file operations

      ! Open the file for reading
         open(unit=unit_number, file=filename, status='old', action='read', iostat=ios)
         if (ios /= 0) then
            print *, "Error: Unable to open file: ", filename
            stop
         end if

      ! Read the file line by line
         do
            read(unit_number, '(A)', iostat=ios) line
            if (ios /= 0) then
               exit  ! Exit loop at end of file or error
            end if

            start_pos = 1
            do
            ! Find the position of '=' in the line
               pos = index(line(start_pos:), '=')
               if (pos == 0) then
                  exit  ! No more key=value pairs on this line
               end if

            ! Extract the key (variable name) from the line
               key = adjustl(trim(line(start_pos:start_pos + pos - 2)))

            ! Find the end of the value (delimited by ',' or end of line)
               start_pos = start_pos + pos
               pos = index(line(start_pos:), ',')
               if (pos == 0) then
                  value_str = adjustl(trim(line(start_pos:)))
                  start_pos = len(line) + 1
               else
                  value_str = adjustl(trim(line(start_pos:start_pos + pos - 2)))
                  start_pos = start_pos + pos
               end if

            ! Check if the key matches the requested variable name
               if (trim(key) == trim(var_name)) then
                  if (is_real) then
                  ! Replace ',' with '.' for decimal separator if needed
                     value_str = adjustl(trim(value_str))
                     call replace_comma_with_dot(value_str)
                     read(value_str, *, iostat=ios) real_value
                  else
                     read(value_str, *, iostat=ios) int_value
                  end if

                  if (ios == 0) then
                     if (ios == 0) then
                     found = .true.
                     exit
                     end if
                  end if
               end if
            end do

            if (found) then
               exit
            end if
         end do

      ! Close the file
         close(unit_number)

      ! Check if the variable was found
         if (.not. found) then
            print *, "Error: Variable not found: ", var_name
            stop
         end if
      end subroutine get_variable_from_file


      subroutine replace_comma_with_dot(value_str)
         implicit none
         character(len=*), intent(inout) :: value_str
         integer :: i

         do i = 1, len(value_str)
            if (value_str(i:i) == ',') then
                  value_str(i:i) = '.'
            end if
         end do
      end subroutine replace_comma_with_dot


      subroutine read_file_integer(jobdir,file,searchLine,fillVariable,m,n,gap,noTitle)

         character (len=*) :: file
         character (len=*) :: searchLine
         integer :: fillVariable(m,n)
         character (len=256) :: jobdir
         character (len=256) :: filename
         integer :: i,j,m,n
         integer :: dummyVariable
         integer, optional :: gap
         logical, optional :: noTitle

      ! allocate(dummyVariable(gap))
         
         
         filename=trim(jobdir)//"\"//trim(file)
!
         open(UNIT=14,file=filename, status='old')

         if (present(noTitle)) then
            call perform_read_int(14, gap, dummyVariable, fillVariable, m, n)
         else if (Searstr(14, trim(searchLine))) then
            call perform_read_int(14, gap, dummyVariable, fillVariable, m, n)
         else
            stop '###..Error en lectura'
         end if
         
         close(14)
      end subroutine read_file_integer
!
!
      subroutine read_file_real(jobdir,file,searchLine,fillVariable,m,n,gap,noTitle)

         character (len=*) :: file
         character (len=*) :: searchLine
         real*8 :: fillVariable(m,n)
         character (len=256) :: jobdir
         character (len=256) :: filename
         integer :: i,j,m,n
         real*8 :: dummyVariable
         integer, optional :: gap
         logical, optional :: noTitle
         
         filename=trim(jobdir)//"\"//trim(file)
!
         open(UNIT=14, file=filename, status='old')

         if (present(noTitle)) then
            call perform_read_real(14, gap, dummyVariable, fillVariable, m, n)
         else if (Searstr(14, trim(searchLine))) then
            call perform_read_real(14, gap, dummyVariable, fillVariable, m, n)
         else
            stop '###..Error en lectura'
         end if

         close(14)
      end subroutine read_file_real
!
      subroutine perform_read_int(unit, gap, dummyVariable, fillVariable, m, n)
         integer, intent(in) :: unit, m, n
         integer, intent(out) :: dummyVariable
         integer, intent(out) :: fillVariable(m, n)
         integer, optional, intent(in) :: gap
         integer :: i, j

         if (present(gap)) then
            READ(unit, *)(dummyVariable, (fillVariable(i, j), j = 1, n), i = 1, m)
         else
            READ(unit, *)((fillVariable(i, j), j = 1, n), i = 1, m)
         end if
      end subroutine perform_read_int
!
      subroutine perform_read_real(unit, gap, dummyVariable, fillVariable, m, n)
         integer, intent(in) :: unit, m, n
         real*8, intent(out) :: dummyVariable
         real*8, intent(out) :: fillVariable(m, n)
         integer, optional, intent(in) :: gap
         integer :: i, j

         if (present(gap)) then
            READ(unit, *)(dummyVariable, (fillVariable(i, j), j = 1, n), i = 1, m)
         else
            READ(unit, *)((fillVariable(i, j), j = 1, n), i = 1, m)
         end if
      end subroutine perform_read_real

      subroutine read_loads(jobdir,file,KINC,fillInt,fillReal,m,n1,n2)
         character (len=*) :: file
         character (len=256) :: stepString, loadString
         integer :: fillInt(m,n1)
         real*8 :: fillReal(m,n2)
         character (len=256) :: jobdir
         character (len=256) :: filename
         integer :: initialIndex, finalIndex
         integer :: i,j,m,n1,n2,KINC
         character(256) :: formatLine
         character(256) :: wholeLine
         
         write(loadString, *) KINC
         stepString = trim('*Step, name=LoadStep' // trim(adjustl(loadString)))

         formatLine = '('
         do i=1,n1
            formatLine = trim(formatLine)//'I5,'
         enddo
         do i=1,n2
            if (i==n2) then
               formatLine = trim(formatLine)//'F10.5)'
            else
               formatLine = trim(formatLine)//'F10.5,'
            end if
         end do
         
         filename=trim(jobdir)//"\"//trim(file)

      ! Find the initial index (first non-space character)
         initialIndex = 1
         do while (initialIndex <= len(stepString) .and. stepString(initialIndex:initialIndex) == " ")
            initialIndex = initialIndex + 1
         end do

      ! Find the final index (last non-space character)
         finalIndex = len(stepString)
         do while (finalIndex >= 1 .and. stepString(finalIndex:finalIndex) == " ")
            finalIndex = finalIndex - 1
         end do

         open(UNIT=15, file=filename, status='old')
         if (Searstr(15,stepString(initialIndex:finalIndex))) then
         ! Number of element loads per load, starts 2 lines after step line
            read(15, '(A)') wholeLine
            do i=1,listNElementLoads(KINC,1)
               read(15, formatLine) fillInt(i,1:n1), fillReal(i,1:n2)
            enddo
         else
            stop '###..Error en lectura'
         end if
         close(15)


      end subroutine read_loads

      subroutine readNLoads(filename, loadsLine, var)
         integer, intent(inout) :: var
         character(len=256), intent(in) :: filename
         integer :: ios, i
         character(len=256), intent(inout) :: loadsLine
     
         ! Open the file
         open(unit=15, file=filename, status='old', form='formatted', iostat=ios)
         if (ios /= 0) then
             print *, "Error: Unable to open file. IOSTAT =", ios
             stop
         end if
     
         ! Read the line containing the loads
         read(15, '(A)', iostat=ios) loadsLine
         if (ios /= 0) then
             print *, "Error: Unable to read from file. IOSTAT =", ios
             stop
         end if

         ! Close the file
         close(15)
     
         ! Count the number of integers (comma-separated values)
         var = count_char_in_string(loadsLine, ",") + 1

      end subroutine readNLoads

      
      subroutine readVLoads(loadsLine, var, vectorNLoads)
         integer, intent(in) :: var
         character(len=256), intent(in) :: loadsLine
         integer, intent(inout) :: vectorNLoads(var, 1)
         integer :: pos, startPos, endPos
         integer :: i
         character(len=256) :: tempStr
     
         ! Initialize the starting position
         startPos = 1
     
         ! Parse the integers from the line
         do i = 1, var
             pos = index(loadsLine(startPos:), ",")  ! Find the next comma
             if (pos == 0) then
                 endPos = len_trim(loadsLine)  ! Last value
             else
                 endPos = startPos + pos - 2
             end if
     
             ! Extract the substring
             tempStr = trim(loadsLine(startPos:endPos))
     
             ! Remove invalid characters (e.g., '**')
             tempStr = adjustl(tempStr)
             if (index(tempStr, "**") /= 0) then
                 tempStr = tempStr(index(tempStr, "**") + 2:)  ! Skip past '**'
             end if
     
             ! Validate the substring
             if (len_trim(tempStr) == 0 .or. verify(tempStr, '0123456789 ') /= 0) then
                 print *, "Error: Invalid input in loadsLine: '", tempStr, "'"
                 stop
             end if
     
             ! Convert the substring to an integer
             read(tempStr, *) vectorNLoads(i, 1)
     
             ! Update the start position for the next value
             startPos = endPos + 2
         end do
      end subroutine readVLoads
      
      end module reading_utils
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------UEXTERNALDB-------------------
!
!     Rutina para leer condiciones fuente y condiciones de archivos externos
!
!     Abre y cierra los archivos necesarios para el calculo
!
!-------------------------------------------------------------------------------------------
!
      subroutine UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
!
      use debug_utils
      use reading_utils
      include 'ABA_PARAM.INC'
!
      character(256)        JOBDIR
      character(256)        JOBNAME

      character*276         folderName
      character*276         line, key
      integer               ios, pos
      integer               i,j,k,elem

      call printStepInfo(LOP, KSTEP, KINC, 'UEXTERNALDB: ')

      if (LOP.eq.0) then ! Variables llamadas al comienzo del análisis
         call GETOUTDIR(JOBDIR,LENJOBDIR)
         call GETJOBNAME(JOBNAME,LENJOBNAME)
         call get_global_variables(jobdir,'parametros.txt')
         call allocate_arrays()
         call initialize_arrays()
         call initialize_results()
         call read_file_real(jobdir,'propiedades.csv','*UEL PROPERTY',propiedades,numMats,numProps)
         call read_file_integer(jobdir,'gruposFisicos.txt','Element Tag, Physical Group Tag',grupoFisico,NELEMS,2)
         call read_file_integer(jobdir,'gruposFisicosN.txt','Node Tag, Physical Group Tag',grupoFisicoN,NUMNODE,2)
         call read_file_real(jobdir,'nodos.inp','*NODE,NSET=N2',nodes,NUMNODE,dim,gap=1)
         call read_file_integer(jobdir,'conectividades.inp','*ELEMENT,TYPE=U1,ELSET=UEL',conectividades,NELEMS,nnod+1)
      ! call read_file_integer(jobdir,'adyacenciaElementos.inp','elem, face1, face2, face3, face4',adyacencias,NELEMS,nnod+1)
!
      else if (LOP.eq.1) then ! Variables llamadas al comienzo de cada paso
         ! Reinicialización de variables
         if (dim == 2) then
            call GETOUTDIR(JOBDIR,LENJOBDIR)
            call read_loads(jobdir,'carga.inp',KINC,elementFaces,elementLoads,maxNElementLoads,2,1)
         end if
      else if (LOP.eq.6) then
!
         ! Calculo de la condición de actualización de propiedades
         CMICriteria(:) = cumulativeResElem(:, 12)

         call calculateCMIthreshold()
         
         ! Actualización de propiedades
         call updateProps()

         ! Inicialización de resultados
         call initialize_results()
      else if (LOP.eq.2) then

      end if
      if (KINC.eq.1) then ! Variables llamadas al comienzo de cada paso.'
         call detectBorders(resElem(:,17),resNod(:,6))
      end if
      return
      end