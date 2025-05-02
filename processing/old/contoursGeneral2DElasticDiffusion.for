CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Elaborado por Diego Alexander Garzón Alvarado
C                 Oscar Rodrigo Lopez
C                 Profesores Asociados Ingenieria Mecánica y Mecatrónica
C                 Universidad Nacional de Colombia
C                 Sede Bogotá
C 
C   copy right - Universidad Nacional de Colombia
C
C   Software general para la solución de ecuaciones de reacción-convección-difusión
C   en 2D.
C
C   La aproximación se realiza mediante el método de Petrov-Galerkin en ABAQUS
C
C   SALIDA
C   La solución de este problema arroja un archivo *.dat que contiene los resultados
C
C   ENTRADA:
C   Requiere un archivo de conectividades.inp ==> archivo con las conectividades de los elementos
C                          nodos.inp          ==> archivo con los puntos nodales
C                          analisis.inp       ==> script de ABAQUS
C                          conec.txt          ==> variables globales
C--------------------------------------------------------------------------------------------------
      module common_utils
         implicit none

         ! Declare variables to hold the retrieved values
         integer :: NUMNODE, NELEMS, dim, nnod
         integer :: a2e, a3e, a4e, be0
         integer :: nLoads, nResNod, nResElem
         integer :: maxNElementLoads
         integer :: numProps, numMats
         integer :: axi, tipo_def
         integer :: nContornos
         integer :: filasContorno1
         integer :: filasContorno2
         integer :: velocidad
         ! real*8 :: OIthreshold
         real*8 :: stdWeight

         ! Declare array variables as allocatable
         real*8, allocatable :: elementLoads(:) 
         integer, allocatable :: elementFaces(:,:)
         integer, allocatable :: listNElementLoads(:)
         real*8, allocatable :: nodes(:,:)
         integer, allocatable :: conectividades(:,:)
         integer, allocatable :: grupoFisico(:,:)
         logical, allocatable :: cambioGrupoFisico(:)
         integer, allocatable :: sigGrupoFisico(:)
         integer, allocatable :: advanceElements(:,:)
         real*8, allocatable :: propiedades(:,:)
         integer, allocatable :: contorno1(:,:), contorno2(:,:)
         real*8, allocatable :: resNod(:,:), resElem(:,:)
         real*8, allocatable :: cumulativeResNod(:,:), cumulativeResElem(:,:)


      contains

      subroutine get_variable_from_file(filename, var_name, is_real, int_value, real_value)
         implicit none
         character(len=*), intent(in) :: filename   ! Name of the input file
         character(len=*), intent(in) :: var_name   ! Name of the variable to retrieve
         logical, intent(in) :: is_real             ! Flag to indicate if the variable is real*8
         integer, intent(out), optional :: int_value ! Value of the variable (integer)
         real*8, intent(out), optional :: real_value ! Value of the variable (real*8)
         character(len=256) :: line                 ! Line read from the file
         character(len=256) :: key                  ! Key extracted from the line
         character(len=256) :: value_str            ! Value as a string
         integer :: ios, pos, start_pos             ! I/O status, position of '=', and start position
         logical :: found                           ! Flag to indicate if variable is found
         integer :: unit_number                     ! Unit number for the file

         found = .false.
         unit_number = 10  ! Arbitrary unit number for file operations

         ! Open the file for reading
         open(unit=unit_number, file=filename, status='old', action='read', iostat=ios)
         if (ios /= 0) then
            print *, "Error: Unable to open file: ", filename
            stop
         end if

         ! Read the file line by line
         do
            read(unit_number, '(A)', iostat=ios) line
            if (ios /= 0) exit  ! Exit loop at end of file or error

            start_pos = 1
            do
                  ! Find the position of '=' in the line
                  pos = index(line(start_pos:), '=')
                  if (pos == 0) exit  ! No more key=value pairs on this line

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
                        found = .true.
                        exit
                     end if
                  end if
            end do

            if (found) exit
         end do

         ! Close the file
         close(unit_number)

         ! Check if the variable was found
         if (.not. found) then
            print *, "Error: Variable not found: ", var_name
            stop
         end if
      end subroutine get_variable_from_file
            
         subroutine test_get_variable_from_file(filename)
            implicit none

            ! Input file containing the variables
            character(len=256) :: filename

            ! Call the subroutine for each variable and print the results
            call get_variable_from_file(filename, "NUMNODE", .false., NUMNODE)
            print *, "NUMNODE =", NUMNODE

            call get_variable_from_file(filename, "NELEMS", .false., NELEMS)
            print *, "NELEMS =", NELEMS

            call get_variable_from_file(filename, "dim", .false., dim)
            print *, "dim =", dim

            call get_variable_from_file(filename, "nnod", .false., nnod)
            print *, "nnod =", nnod

            call get_variable_from_file(filename, "a2e", .false., a2e)
            print *, "a2e =", a2e

            call get_variable_from_file(filename, "a3e", .false., a3e)
            print *, "a3e =", a3e

            call get_variable_from_file(filename, "a4e", .false., a4e)
            print *, "a4e =", a4e

            call get_variable_from_file(filename, "be0", .false., be0)
            print *, "be0 =", be0

            call get_variable_from_file(filename, "nLoads", .false., nLoads)
            print *, "nLoads =", nLoads

            call get_variable_from_file(filename, "nResNod", .false., nResNod)
            print *, "nResNod =", nResNod

            call get_variable_from_file(filename, "nResElem", .false., nResElem)
            print *, "nResElem =", nResElem

            call get_variable_from_file(filename, "maxNElementLoads", .false., maxNElementLoads)
            print *, "maxNElementLoads =", maxNElementLoads

            call get_variable_from_file(filename, "numProps", .false., numProps)
            print *, "numProps =", numProps

            call get_variable_from_file(filename, "numMats", .false., numMats)
            print *, "numMats =", numMats

            call get_variable_from_file(filename, "axi", .false., axi)
            print *, "axi =", axi

            call get_variable_from_file(filename, "tipo_def", .false., tipo_def)
            print *, "tipo_def =", tipo_def

            call get_variable_from_file(filename, "nContornos", .false., nContornos)
            print *, "nContornos =", nContornos

            call get_variable_from_file(filename, "filasContorno1", .false., filasContorno1)
            print *, "filasContorno1 =", filasContorno1

            call get_variable_from_file(filename, "filasContorno2", .false., filasContorno2)
            print *, "filasContorno2 =", filasContorno2

            call get_variable_from_file(filename, "velocidad", .false., velocidad)
            print *, "velocidad =", velocidad

            ! call get_variable_from_file(filename, "OIthreshold", .true., real_value=OIthreshold)
            ! print *, "OIthreshold =", OIthreshold

            call get_variable_from_file(filename, "stdWeight", .true., real_value=stdWeight)
            print *, "stdWeight =", stdWeight
         end

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


         subroutine initialize_arrays()
            implicit none

            ! Allocate arrays based on the retrieved integer values
            allocate(elementLoads(maxNElementLoads))
            allocate(elementFaces(maxNElementLoads, 2))
            allocate(listNElementLoads(nLoads))
            allocate(nodes(NUMNODE, dim))
            allocate(conectividades(NELEMS, nnod+1))
            allocate(grupoFisico(NELEMS, 2))
            allocate(cambioGrupoFisico(NELEMS))
            allocate(sigGrupoFisico(NELEMS))
            allocate(advanceElements(a4e, a3e + 2*a2e))
            allocate(propiedades(numMats, numProps))
            allocate(contorno1(filasContorno1, 6))
            allocate(contorno2(filasContorno2, 6))
            allocate(resNod(NUMNODE, nResNod))
            allocate(resElem(NELEMS, nResElem))
            allocate(cumulativeResNod(NUMNODE, nResNod))
            allocate(cumulativeResElem(NELEMS, nResElem))

            ! Initialize arrays to default values
            ! elementLoads = 0.0
            ! elementFaces = 0
            ! listNElementLoads=(/27 ,26 ,26 ,26 ,27/) !TODO: Change to actual values
            ! nodes = 0.0
            ! conectividades = 0
            ! grupoFisico = 0
            ! cambioGrupoFisico = .false.
            ! sigGrupoFisico = 0
            ! advanceElements = 0
            ! propiedades = 0.0
            ! contorno1 = 0
            ! contorno2 = 0
            ! resNod = 0.0
            ! resElem = 0.0
            ! cumulativeResNod = 0.0
            ! cumulativeResElem = 0.0
            ! OIthreshold = 0.0
   
         end subroutine initialize_arrays

      end module common_utils
      
      
      subroutine UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS, NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3 NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4 PERIOD)
C     
      use common_utils
      include 'ABA_PARAM.INC'
	   ! include 'conec.for'
C
      dimension RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
C
	   real*8    x(dim,nnod)!,de(8)
      integer   k1,k2
C
      integer debugUnit
      character*276         filename
      character(256)        JOBDIR
      character(256)        JOBNAME
C
      logical change
C
      debugUnit = 15  ! or any unused Fortran unit number
C
C     Se llaman las propiedades del modelo (parametros)
C      call INPUT(props,de,NPROPS)
C
C     Inicializacion

      DO 6 K1=1,NDOFEL                      
		   RHS(K1,NRHS)=0.0
      DO 4 K2=1,NDOFEL
		   AMATRX(K2,K1)=0.0
    4 CONTINUE                                      
    6 CONTINUE   
C
      do k1=1,nnod
         do j=1,dim
            x(j,k1)=coords(j,k1)
         enddo
      enddo

      change = (KSTEP .ne. 1) .and. (KINC .eq. 1) .and. not(cambioGrupoFisico(jelem))

      if (change) then
         call calculateCMIthreshold()
         call updateProps()
      endif
C
C     Funcion que regresa RHS y AMATRX
C    
      call ENSAMBLE(de,x,du,u,v,nst,ndofel,KSTEP,KINC,
     1 nrhs,dtime,svars,nsvars,jelem,time,RHS,AMATRX)
C
C----- DEBUG OUTPUT SECTION ------------------------------------------
      call GETOUTDIR(JOBDIR,LENJOBDIR)
      call GETJOBNAME(JOBNAME,LENJOBNAME)
      filename=' '
      filename(1:lenjobdir)=jobdir(1:lenjobdir)
      filename(lenjobdir+1:lenjobdir+1)='\'
      filename(lenjobdir+2:lenjobdir+lenjobname+1)=jobname(1:lenjobname)
      filename(lenjobdir+lenjobname+2:lenjobdir+lenjobname+6)='.txt'

      if (.false.) THEN
C     1) Open a debug file
      OPEN(UNIT=debugUnit, FILE=filename, ACTION='WRITE', STATUS='UNKNOWN')
C     Write debug information including Step and Increment
      WRITE(debugUnit,*) 'Step: ', KSTEP, 'Increment: ', KINC
      WRITE(debugUnit,*) 'UEL subroutine executed at this increment.'

C     2) Write a header
      WRITE(debugUnit,*) '====================================================='
      WRITE(debugUnit,*) ' ELEMENT MATRIX (AMATRX) AND RHS DEBUG OUTPUT'
      WRITE(debugUnit,*) ' Element number = ', jelem
      WRITE(debugUnit,*) ' ndofel = ', ndofel
      WRITE(debugUnit,*) '====================================================='

C     3) Write the stiffness matrix
      WRITE(debugUnit,*) 'AMATRX:'
      do i = 1, NDOFEL 
      WRITE(debugUnit,'(100(1x, E14.6))') (AMATRX((i),(j)), j=1, NDOFEL)
      enddo

C     4) Write the RHS vector
      WRITE(debugUnit,*) 'RHS ='

      do i=1, MLVARX
      WRITE(debugUnit,'(I6, 2x, E14.6)') (i), RHS((i),1)
      enddo

C	  5) Write the U vector
      WRITE(debugUnit,*)
      WRITE(debugUnit,*) 'U ='
      do i=1, NDOFEL
      WRITE(debugUnit,'(I6, 2x, E14.6)') (i), U((i))
      enddo

C     6) Close the file
      close(debugUnit)
      END IF
C----- END DEBUG SECTION ---------------------------------------------

      return
      end
C---------------------------------------------------------------------------------
C------------------------------------------------------------------------- BGAUSS2 -------
C
C     Funcion bgauss(sg,wg)
C
C----------------------------------------------------------------------------------------
      SUBROUTINE bgauss(sg,wg)
      real*8 sg(2),wg(2)

      sg(1) = -0.577350269189626
      sg(2) =  0.577350269189626

      wg(1) =  1.0
      wg(2) =  1.0

      return
      end
C-----------------------------------------------------------------------------------------
C-----------------------------------------------------------------------Searstr----------- 
C     Funcion que devuelve un valor logico que permite saber si la lectura es 
C     adecuadad.
C 
C     Busca Str en la unidad Lu y devuelve .true.// si lo encuentra. Si no
C     lo encuentra devuelve un .false.
C
C-----------------------------------------------------------------------------------------
	function Searstr(Lu, Str)
C
	integer            :: Lu,j,L1
	character (len=*)  :: Str
	character (len=80) :: Rec
	logical            :: Salida,Searstr
C
	Searstr = .false.
	Salida  = .true.
	Nstrlen = len(Str)
	L1      = Nstrlen
	j       = 1
C
	do while (Salida .and. (j <= 2))
	 read (Lu, '(a80)',err=10) Rec
	 if (Rec(1:L1) == (Str)) then
        Salida  = .false.
        Searstr = .true.
     	 else if (Rec(1:4) .eq.('*END')) then
         rewind (Lu)   ! Rebobina el archivo
         j=j+1
       end if
	enddo
C
	return
   10 	Write (*,20) Lu,Str
	stop
C
   20 format ('###..Error en la funcion Searstr(Unidad=)',I3,2X,15A)
C
	end function Searstr
C---------------------------------------------------------------------------------------------
C----------------------------------------------------------------------Matrices_de_calculo2D--
C    
C     Devuelve las matrices de calculo requeridas en los analisis escalares
C
C---------------------------------------------------------------------------------------------
      subroutine Matrices_de_calculo2D(shp,N,B,Ne,Be)
C     
      use common_utils
      ! include 'conec.for'
C
      integer  col,n1
      real*8  shp(3,nnod),N(1,nnod),Ne(2,2*nnod), B(2,nnod),Be(4,2*nnod)
C
C     Puesta a cero
      N = 0.d0
      Ne= 0.d0
      B = 0.d0
      Be= 0.d0
C
C     Incializacion
      do n1 = 1,nnod
         col=2*(n1-1)+1
C
C        Funciones de forma para descripcion de variable escalar
         N(1,n1)    = shp(1,n1)
C
C        Funciones de forma para descripcion de variable vectorial	  
         Ne(1,col)  = shp(1,n1)
         Ne(2,col+1)= shp(1,n1)
C
C        Funciones de forma derivadas en el espacio para variable escalar	
         B(1,n1)    = shp(2,n1)
         B(2,n1)    = shp(3,n1)
C
C        Funciones de forma derivadas para el caso elastico
         Be(1,col)  = shp(2,n1)
         Be(2,col+1)= shp(3,n1)
         Be(3,col)  = shp(3,n1)
         Be(3,col+1)= shp(2,n1)
         Be(4,col)  = shp(1,n1)
	   enddo  
C
      return
      end
C----------------------------------------------------------------------------------------
C------------------------------------------------------------------------ F_FORMA2D------
C	  Funcion f_forma2D(chi,eta,x,shp,xjac,d2shp)
C
C----------------------------------------------------------------------------------------
      subroutine f_forma2D(chi,eta,x,shp,xjac)
C	
      use common_utils
      ! include 'conec.for'
C
      integer i
      real*8 chi,eta,xjac,shp(3,nnod),x(2,nnod),d2shp(3,nnod)
      real*8 dxchi,dxeta,dychi,dyeta,dchix,dchiy,detax,detay
      real*8 dNchi(nnod),dNeta(nnod)
      real*8 d2Nchi_eta(nnod)
      real*8 d2xchi_eta,d2ychi_eta,dxjacchi,dxjaceta,dxjacx,dxjacy
      real*8 d2eta_xx,d2chi_xx,d2chi_yy,d2eta_yy,d2eta_xy,d2chi_xy
C
      dxchi=0.d0
      dxeta=0.d0
      dychi=0.d0
      dyeta=0.d0
      dchix=0.d0
      dchiy=0.d0
      detax=0.d0
      detay=0.d0
C
      shp=0.d0
C
c     Primera fila de la matriz shp - Funciones de forma en los 8 nodos
      shp(1,1) = 0.25*(1.0 - chi)*(1.0 - eta)
      shp(1,2) = 0.25*(1.0 + chi)*(1.0 - eta)
      shp(1,3) = 0.25*(1.0 + chi)*(1.0 + eta)
      shp(1,4) = 0.25*(1.0 - chi)*(1.0 + eta)
C
C     Primeras derivadas de las funciones de forma con respecto a chi
      dNchi(1) = -0.25*(1.0 - eta)
      dNchi(2) =  0.25*(1.0 - eta)
      dNchi(3) =  0.25*(1.0 + eta)
      dNchi(4) = -0.25*(1.0 + eta)
C
C     Primeras derivadas de las funciones de forma con respecto a eta
      dNeta(1) = -0.25*(1.0 - chi)
      dNeta(2) = -0.25*(1.0 + chi)
      dNeta(3) =  0.25*(1.0 + chi)
      dNeta(4) =  0.25*(1.0 - chi)
C

C     Segundas derivadas de las funciones de forma con respecto a chi y eta
      d2Nchi_eta(1) =  0.25
      d2Nchi_eta(2) = -0.25
      d2Nchi_eta(3) =  0.25
      d2Nchi_eta(4) = -0.25
C
c     Calculo de la matriz jacobiana
      dxchi=0.25*((x(1,2)-x(1,1))*(1.0-eta)+(x(1,3)-x(1,4))*(1.0+eta))
      dxeta=0.25*((x(1,4)-x(1,1))*(1.0-chi)+(x(1,3)-x(1,2))*(1.0+chi))
      dychi=0.25*((x(2,2)-x(2,1))*(1.0-eta)+(x(2,3)-x(2,4))*(1.0+eta))
      dyeta=0.25*((x(2,4)-x(2,1))*(1.0-chi)+(x(2,3)-x(2,2))*(1.0+chi))
C
C     Calculo de las segundas derivadas de x y y
      d2xchi_eta=0.25*(-(x(1,2)-x(1,1))+(x(1,3)-x(1,4)))
      d2ychi_eta=0.25*(-(x(2,2)-x(2,1))+(x(2,3)-x(2,4)))
C
c     Calculo del determinante de la matriz jacobiana -Jacobiano -
C
      xjac = dxchi*dyeta - dxeta*dychi
C
c     Calculo de la matriz inversa de la matriz jacobiana
      dchix =   dyeta/xjac
      dchiy = - dxeta/xjac
      detax = - dychi/xjac
      detay =   dxchi/xjac
C
C     Calculo de las derivadas del jacobiano con respecto a chi
      dxjacchi = dxchi*d2ychi_eta - d2xchi_eta*dychi
      dxjaceta = d2xchi_eta*dyeta - dxeta*d2ychi_eta
C
C     Calculo de las derivadas del jacobiano con respecto a x
      dxjacx = dxjacchi*dchix + dxjaceta*detax
C     Calculo de las derivadas del jacobiano respecto a y
      dxjacy = dxjacchi*dchiy + dxjaceta*detay
C
C     Calculo de la segunda derivada de eta con respecto a x
      d2eta_xx=((-d2ychi_eta*detax*xjac)+(dxjacx*dychi))/(xjac**2)
C
C     Calculo de la segunda derivada de chi con respecto a x
      d2chi_xx=((d2ychi_eta*dchix*xjac)-(dxjacx*dyeta))/(xjac**2)
C 
C     Calculo de la segunda derivada de chi con respecto a y
      d2chi_yy=(-(d2xchi_eta*dchiy*xjac)+(dxjacy*dxeta))/(xjac**2)
C
C     Calculo de la segunda derivada de eta con respecto a y
      d2eta_yy=((d2xchi_eta*detay*xjac)-(dxjacy*dxchi))/(xjac**2)
C
C     Calculo de la segunda derivada de eta con respecto a x y luego y
      d2eta_xy=((-d2ychi_eta*detay*xjac)+(dxjacy*dychi))/(xjac**2)
C
C     Calculo de la segunda derivada de chi con respecto a x y luego y
      d2chi_xy=((d2ychi_eta*dchiy*xjac)-(dxjacy*dyeta))/(xjac**2)
C
C     Calculo de las segundas derivadas, la primera fila indica las segundas derivadas de
C     las funciones de forma con respecto a x, la segunda con respecto a y, y la tercera, los terminos cruzados
      do i =1,nnod     
	   d2shp(1,i)=2.d0*d2Nchi_eta(i)*detax*dchix+dNchi(i)*d2chi_xx+
     &             dNeta(i)*d2eta_xx
	   d2shp(2,i)=2.d0*d2Nchi_eta(i)*detay*dchiy+dNchi(i)*d2chi_yy+
     &             dNeta(i)*d2eta_yy
      d2shp(3,i)=d2Nchi_eta(i)*(detay*dchix+detax*dchiy)+
     &             dNchi(i)*d2chi_xy+dNeta(i)*d2eta_xy
	   enddo
C
C     Calculo de las derivadas de las funciones de forma con respecto a x
C
      shp(2,1) = - 0.25*((1.0-eta)*dchix + (1.0-chi)*detax)
      shp(2,2) =   0.25*((1.0-eta)*dchix - (1.0+chi)*detax)
      shp(2,3) =   0.25*((1.0+eta)*dchix + (1.0+chi)*detax)
      shp(2,4) = - 0.25*((1.0+eta)*dchix - (1.0-chi)*detax)

C     Calculo de las derivadas de las funciones de forma con respecto a y
C
      shp(3,1) = - 0.25*((1.0-eta)*dchiy + (1.0-chi)*detay)
      shp(3,2) =   0.25*((1.0-eta)*dchiy - (1.0+chi)*detay)
      shp(3,3) =   0.25*((1.0+eta)*dchiy + (1.0+chi)*detay)
      shp(3,4) = - 0.25*((1.0+eta)*dchiy - (1.0-chi)*detay)
C
      return
      end
C---------------------------------------------------------------------------------------------
C-------------------------------------------------------------updateProps-----------
C
C      Funcion updateProps
C	
C      Funcion que actualiza los grupos físicos.
C
C-------------------------------------------------------------------------------*/
C-------------------------------------------------------------------------------*/
C
      subroutine updateProps()
C
      use common_utils
      ! include 'conec.for'
C
C     Actualización de propiedades
C     
      logical :: ossByFront, ossBySOC, isEarlyOss
C
      ossByFront = .false.
      do i = 1, a2e + a3e + a2e
         do j = 1, velocidad
            if (advanceElements(velocidad*(KSTEP-2) + j, i) == jelem) then
               ossByFront = .true.
               exit
            end if
         end do
      end do

      ossBySOC = cumulativeResElem(jelem,15) >= CMIthreshold

      isEarlyOss = grupoFisico(jelem,2) == 3

C     Centro secundario de osificacion
      if (isEarlyOss) then
         grupoFisico(jelem,2) = sigGrupoFisico(jelem)
         cambioGrupoFisico(jelem) = .true.
      elseif (grupoFisico(jelem,2) == 2) then
         if (ossBySOC) then
            grupoFisico(jelem,2) = 3
            sigGrupoFisico(jelem) = 4
            cambioGrupoFisico(jelem) = .true.
         elseif (ossByFront) then
            grupoFisico(jelem,2) = 3
            sigGrupoFisico(jelem) = 1
            cambioGrupoFisico(jelem) = .true.
         endif
      endif

      return
      end
C---------------------------------------------------------------------------------------------
C-------------------------------------------------------------ENSAMBLE-----------
C
C      Funcion ENSAMBLE
C	
C      Funcion que regresa la matriz de rigidez tangente AMATRX y
C      el residuo RHS.
C
C-------------------------------------------------------------------------------*/
C-------------------------------------------------------------------------------*/
C
      subroutine ENSAMBLE(de,x,du,u,v,nst,ndofel,KSTEP,KINC,
     1 nrhs,dtime,svars,nsvars,jelem,time,p,m_k)
C
      
      use common_utils
      ! include 'conec.for'
C
C     Entradas
      integer nst,ndofel,nrhs,nsvars,jelem
      real*8  x(dim,nnod),du(ndofel,*),u(ndofel),v(ndofel)!,de(8)
	   real*8  dtime,time(2),svars(nsvars)
C
C     Salidas
      real*8  p(ndofel,nrhs),m_k(ndofel,ndofel)
C
C     Arreglos y variables de la subrutina
C
C     Generales
      integer k1,k2,k3,cc,ff,fil,col,tipo,colg,filg,colp,filp
      integer ndofnod,ndofdiff
	   real*8  paux(ndofel),Kelast(dim*nnod,dim*nnod) ! 2 gdl correspondientes a elasticidad
      real*8  Kdiff(2*nnod,2*nnod) ! 2 gdl correspondientes a difusion
      real*8  Def(dim*nnod)
      real*8  Dif(2)
C
C     Variables de la carga distribuida
      integer n1, n2, i, j, k
      real*8 mag
      real*8  x1(dim),x2(dim)
      real*8  f, fx, fy

      logical :: found
C     real*8  len, angulo, 
C     Inicializacion del tiempo
      if (dtime.eq.0.0) then
        dtime=1.e-15
      end if
C
C     Inicializacion de matrices y variables
      m_k   =   0.d0
      p     =   0.d0
      ndofnod = ndofel/nnod 
      ndofdiff = 2
C
C     Ensamblar la matriz de rigidez elastica en la matriz tangente global
C     Llamado de la matriz de rigidez para la expansion
      call matriz_rigidez_el(u,ndofel,de,x,jelem,time(2),Kelast)
      call matriz_rigidez_dif(u,ndofdiff,ndofel,de,x,jelem,time(2),Kdiff)
C
      do k1=1,nnod
         do k2=1,nnod
            fil = ndofnod*(k1-1)+1
            col = ndofnod*(k2-1)+1
            ff = dim*(k1-1)+1
            cc = dim*(k2-1)+1
            do i=1,dim
               do j=1,dim
                  colg = col+(j-1)
                  filg = fil+(i-1)
                  colp = cc +(j-1)
                  filp = ff +(i-1) 
                  m_k(filg,colg) = m_k(filg,colg) + Kelast(filp,colp)
               enddo
            enddo
            do i=1,ndofdiff
               do j=1,ndofdiff
                  colg = col+(j-1)+dim
                  filg = fil+(i-1)+dim
                  colp = cc +(j-1)
                  filp = ff +(i-1) 
                  m_k(filg,colg) = m_k(filg,colg) + Kdiff(filp,colp)
               enddo
            enddo
         enddo
      enddo
C
C     RHS:
	   paux=matmul(m_k,u)
C
C     Ensamble del vector residuo
      do k1=1,ndofel
	      p(k1,1)= - paux(k1)
      enddo
C
C     Insercion de carga distribuida
      found = .false.

      do i = 1, listNElementLoads(KINC)
         if (jelem == elementFaces(i, 1)) then
            found = .true.
            exit
         end if
      enddo

C     TODO: Use a similar approach to KDLOAD, NDLOAD and ADLMAG to load elements
C     DO KDLOAD = 1, NDLOAD
C         n1 = JDLTYP(KDLOAD,1) ! Cara de la carga distribuida
C         mag = ADLMAG(KDLOAD,1) ! Magnitud de la carga distribuida
C     ENDDO

      if (found) then ! Works just for 1 load per element
         n1 = elementFaces(i, 2) ! Cara de la carga distribuida
         mag = elementLoads(i) ! Magnitud de la carga distribuida
         n2 = mod(n1, nnod)+1

         x1 = x(:,n1)
         x2 = x(:,n2)

C         len = sqrt((x2(1)-x1(1))**2 + (x2(2)-x1(2))**2)
C         angulo = atan2(x2(2)-x1(2), x2(1)-x1(1))
         f = mag/2 ! *len
         fx = -f*(x2(2)-x1(2)) !/len
         fy = f*(x2(1)-x1(1)) !/len

         p(ndofnod*(n1-1)+1, 1) = p(ndofnod*(n1-1)+1, 1) + fx
         p(ndofnod*(n1-1)+2, 1) = p(ndofnod*(n1-1)+2, 1) + fy
         p(ndofnod*(n2-1)+1, 1) = p(ndofnod*(n2-1)+1, 1) + fx
         p(ndofnod*(n2-1)+2, 1) = p(ndofnod*(n2-1)+2, 1) + fy
      end if
C     Insercion del vector de deformacion
      Def = 0.0
C
      do k1=1,nnod
         fil = ndofnod*(k1-1)+1
         ff = dim*(k1-1)+1
         do i=1,dim
            filg = fil+(i-1)
            filp = ff +(i-1) 
C	       p(fil,1)=p(fil,1)+Def(ff)
            p(filg,1)=p(filg,1)+Def(filp)
         enddo
      enddo
C
      return
      end
C------------------------------------------------------------------------------------------
C----------------------------------------------------------------matriz_rigidez_el---------
C
C	 Funcion matriz_rigidez_el: matriz de rigidez elastica
c
C------------------------------------------------------------------------------------------
	subroutine matriz_rigidez_el(u,ndofel,de,x,jelem,t,Cmat)
C
      use common_utils
	   ! include 'conec.for'
C
C     Variables de entrada
      integer  ndofel,jelem
      real*8   u(ndofel),x(dim,nnod),t!,de(8)
C
C     variables de salida
      real*8   Cmat(dim*nnod,dim*nnod) !dim seria en este caso el número de grados de libertad por nodo
C
C     Variables de la subrutina
C     Generales
      integer  i,j,k
      real*8   wg(2),sg(2),matriz(dim*nnod,dim*nnod)
      real*8   chi,eta,dx,xjac
      real*8   Dmat(3*(dim-1)+axi,3*(dim-1)+axi)
C
C     En 2D
      real*8   shp2D(3,nnod),Nmat2D(1,nnod)
      real*8   Bmat2D(2,nnod),Nmatel2D(2,2*nnod),Bmatel2D(4,dim*nnod)
      real*8   Bmatplano2D(3,2*nnod),r,r0,Bmataxi2D(4,2*nnod)
C
C     Inicializacion de variables
      xjac	     =   0.d0
      chi	     =   0.d0
      eta        =   0.d0
      zita       =   0.d0
      dx	        =   0.d0
      Nmat2D     =   0.d0
      Bmat2D     =   0.d0
      Nmatel2D   =   0.d0
      Bmatel2D   =   0.d0
      Bmatplano2D=   0.d0
      Bmataxi2D  =   0.d0
      Cmat 	     =   0.d0
      Dmat       =   0.d0
      matriz     =   0.d0
      r0         =   0.d0
C
C	Calculo de los puntos de Gauss
	   call bgauss(sg,wg)
C
C	Bucle para cada punto de Gauss
      if(dim.eq.2)then
	      do i = 1,2
  	         chi = sg(i)
	      do j = 1,2
	         eta = sg(j)
C
C           Inicializacion de la matriz de constantes elasticas
            call matriz_constantes_el(chi,eta,zita,u,ndofel,de,x,
     &           jelem,t,Dmat)
C
C	      Se calculan las funciones de forma
	         call f_forma2D(chi,eta,x,shp2D,xjac)
C
C           Se obtienen las matrices de calculo
            call Matrices_de_calculo2D(shp2D,Nmat2D,Bmat2D,Nmatel2D,
     &           Bmatel2D)
C
C           Se calcula el diferencial de la integral
            dx = wg(i)*wg(j)*xjac
C
C           Tipo de analisis
            if(axi.eq.0)then ! Caso plano: esfuerzo o deformacion plana
	            Bmatplano2D = Bmatel2D(1:3,:)
C              Se calculan las matriz elemental correspondiente a reaccion
	            matriz = matmul(transpose(Bmatplano2D),
     &                   matmul(Dmat(1:3,1:3),Bmatplano2D))
	         r = 1
            elseif(axi.eq.1)then
C              Radio de giro
               r = dot_product(Nmat2D(1,:),x(1,:))+r0	        
C
C              Definicion de la matriz operador de deformaciones
               Bmataxi2D(1:3,:) = Bmatel2D(1:3,:)
               Bmataxi2D(4,:) = Bmatel2D(4,:)/r
C
C              Se calculan las matriz elemental correspondiente a reaccion
	          matriz = matmul(transpose(Bmataxi2D),
     &                   matmul(Dmat(1:4,1:4),Bmataxi2D))
            end if
C
C           Se llevan a cabo las operaciones de sumatoria
	      Cmat= Cmat + matriz * r * dx
	    enddo
	  enddo
	end if
C
	return 
	end
C------------------------------------------------------------------------------------------
C----------------------------------------------------------------matriz_rigidez_dif---------
C
C	 Funcion matriz_rigidez_dif: matriz de rigidez difusiva
C
C   CMat : Matriz de rigidez difusiva
C------------------------------------------------------------------------------------------
	subroutine matriz_rigidez_dif(u,ndofdiff,ndofel,de,x,jelem,t,Cmat)
C
      use common_utils
	   ! include 'conec.for'
C
C     Variables de entrada
      integer  ndofdiff,ndofel,jelem
      real*8   x(dim,nnod),t!,de(8)
C
C     variables de salida
      real*8   Cmat(ndofdiff*nnod,ndofdiff*nnod)
C     Variables de la subrutina
C     Generales
      integer  i,j,k,l,m,n
      real*8   wg(2),sg(2),matriz(ndofdiff*nnod,ndofdiff*nnod)
      real*8   matrizgdl(nnod,nnod)
      real*8   chi,eta,dx,xjac
      real*8   Dmat(3*(dim-1)+axi,3*(dim-1)+axi)
C
C     En 2D
      real*8   shp2D(3,nnod),Nmat2D(1,nnod)
      real*8   Bmat2D(2,nnod),Nmatel2D(2,2*nnod),Bmatel2D(4,2*nnod)
      real*8   Bmatplano2D(3,2*nnod),r,r0,Bmataxi2D(4,2*nnod)
C
C     Inicializacion de variables
      xjac	     =   0.d0
      chi	     =   0.d0
      eta        =   0.d0
      zita       =   0.d0
      dx	        =   0.d0
      Nmat2D     =   0.d0
      Bmat2D     =   0.d0
      Nmatel2D   =   0.d0
      Bmatel2D   =   0.d0
      Bmatplano2D=   0.d0
      Bmataxi2D  =   0.d0
      Cmat 	     =   0.d0
      Dmat       =   0.d0
      matriz     =   0.d0
      r0         =   0.d0
C
C	Calculo de los puntos de Gauss
	   call bgauss(sg,wg)
C
C	Bucle para cada punto de Gauss
      if(dim.eq.2)then
	      do i = 1,2
  	         chi = sg(i)
	      do j = 1,2
	         eta = sg(j)
C
C	         Se calculan las funciones de forma
	         call f_forma2D(chi,eta,x,shp2D,xjac)
C
C           Se obtienen las matrices de calculo
            call Matrices_de_calculo2D(shp2D,Nmat2D,Bmat2D,Nmatel2D,
     &           Bmatel2D)
C
C           Se calcula el diferencial de la integral
            dx = wg(i)*wg(j)*xjac
C
C           Se calcula la matriz de difusion para un solo gdl
            matrizgdl = matmul(transpose(Bmat2D),Bmat2D)
C
C           Se calcula el aporte a la matriz de rigidez elemental
            do k = 1,nnod
               do l = 1,nnod
                  m = 2*k
                  n = 2*l
                  matriz(m,n) = matrizgdl(k,l)
                  matriz(m-1,n-1) = matrizgdl(k,l)
               enddo
            enddo
C
C           
C           Se calculan las matriz elemental correspondiente a reaccion
C	         matriz = matmul(transpose(Bmatplano2D),
C     &                   matmul(Dmat(1:3,1:3),Bmatplano2D))
C           Se llevan a cabo las operaciones de sumatoria
	         Cmat = Cmat + matriz * dx
	      enddo
	      enddo
	end if
C
	return
	end
C------------------------------------------------------------------------------------------
C---------------------------------------------------------------MATRIZ_CONSTANTES_E________
C
C	 Funcion matriz_constantes_e de constantes elásticas
c
C------------------------------------------------------------------------------------------
	subroutine matriz_constantes_el(chi,eta,zita,u,ndofel,de,x,
     1 jelem,t,Dmat)
C
      use common_utils
	   ! include 'conec.for'
C
C     Variables de entrada
      integer  ndofel,jelem
      real*8   chi,eta,zita
	   real*8   u(ndofel),x(dim,nnod),t!,de(8)
C
C     variables de salida
      real*8   Dmat(3*(dim-1)+axi,3*(dim-1)+axi)
C
C     Variables de la subrutina
C     Generales
      integer  i,j,k
      real*8   E,nu,d11,d22,d33,d12
      real*8   Dmatb(3*(dim-1)+axi,3*(dim-1)+axi)

      character*276         filename
      character(256)        JOBDIR
      character(256)        JOBNAME
      character(256) :: logFileName
      integer :: logUnit, ierr

      logUnit = 14  ! Define a unit number for the log file
C
C     Definicion del tensor de constantes elasticas
      Dmat   = 0.d0
C
      E      = propiedades(grupoFisico(jelem,2),1)
      nu     = propiedades(grupoFisico(jelem,2),2)
C
      if(dim.eq.2)then
C       Analisis bidimensional
         if(axi.eq.0)then ! Analisis plano
            if(tipo_def.eq.1)then ! Esfuerzo plano
C            Constantes
               d11 = E/(1.d0-nu**2)
               d22 = d11
               d12 = nu*d11
               d33 = 0.5d0*E/(1.d0+nu)   
C
C            Definicion de la matriz de constantes
               Dmat(1,1) = d11
               Dmat(2,2) = d22
               Dmat(1,2) = d12
               Dmat(2,1) = Dmat(1,2)
               Dmat(3,3) = d33


C 
	         elseif(tipo_def.eq.2)then ! Deformacion plana
C            Constantes
               d11 = E*(1.d0-nu)/((1.d0-2.d0*nu)*(1.d0+nu))
               d22 = d11
               d12 = d11*nu/(1.d0-nu)
               d33 = 0.5d0*E/(1.d0+nu)   
C
C            Definicion de la matriz de constantes
               Dmat(1,1) = d11
               Dmat(2,2) = d22
               Dmat(1,2) = d12
               Dmat(2,1) = Dmat(1,2)
               Dmat(3,3) = d33
C 
            end if
C
C       Axisimetrico
         elseif(axi.eq.1)then ! Analisis axisimetrico
C         Constantes
            d11 = E*(1.d0-nu)/((1.d0-2.d0*nu)*(1.d0+nu))
            d22 = d11
            d12 = d11*nu/(1.d0-nu)
            d33 = 0.5d0*E/(1.d0+nu)   
C         Definicion de la matriz de constantes
            Dmat(1,1) = d11
            Dmat(2,2) = d22
            Dmat(1,2) = d12
            Dmat(2,1) = Dmat(1,2)
            Dmat(3,3) = d33
            Dmat(1,4) = d12
            Dmat(2,4) = d12
            Dmat(4,1) = Dmat(1,4)
            Dmat(4,2) = Dmat(2,4)
            Dmat(4,4) = d11
	      end if
	   end if

C----- DEBUG OUTPUT SECTION ------------------------------------------
      call GETOUTDIR(JOBDIR,LENJOBDIR)
      call GETJOBNAME(JOBNAME,LENJOBNAME)
      filename=' '
      filename(1:lenjobdir)=jobdir(1:lenjobdir)
      filename(lenjobdir+1:lenjobdir+1)='\'
      filename(lenjobdir+2:lenjobdir+8+1)='Dmat_log'
      filename(lenjobdir+8+2:lenjobdir+8+6)='.txt'
C     
      IF (.FALSE.) THEN
      open(unit=logUnit, FILE=filename, ACTION='WRITE', STATUS='UNKNOWN')
      
      ! Write the Dmat matrix to the log file
      WRITE(logUnit,*) 'Dmat matrix:'
      do i = 1, size(Dmat, 1)
      WRITE(logUnit,'(100(1x, E14.6))') (Dmat(i, j), j = 1, size(Dmat, 2))
      enddo

      ! Close the log file
      close(logUnit)
      END IF
C----- END DEBUG SECTION ---------------------------------------------

	return 
	end
C------------------------------------------------------------------------------------------
C---------------------------------------------------------------------------oute-----------
C       
C       Funcion outsigma()
C
C
C       Funcion de calculo de las deformaciones requeridas para los calculos de las 
C       matrices restantes, es llamado al inicio de cada incremento
C
C       Diccionario
C -----------------------------------------------------------------------------------------
C
C      Funcion sin argumentos:
C
C
C       Funcion sin argumentos que permite escribir los datos requeridos para 
C       calculos alternos, por ejemplo, las deformaciones
C
      subroutine outsigma()
C  
      use common_utils
      include 'ABA_PARAM.INC'
      ! include    'conec.for'
C
      integer    i,j,k,jelem,j2,n,M
      real*8     ESFUERZOS(3),bmat(3,nnod*2),
     1           shp(3,nnod),xjac,x(2,nnod),
     2           DEFORMACIONES(3),DESP(nnod*2)
      real*8     dmat2(3*(dim-1),3*(dim-1))
      real*8     chi,eta,sg(2),sigma_oct,thao_oct
      real*8     Nm(1,nnod),Ne(2,2*nnod),Bm(2,nnod),Be(4,2*nnod)
      real*8     I1,I2,I3,r,def_crec(4),dx
      real*8     zita,t,sigma_zz,Et(10),nut(10),Def2D(3*(dim-1)+axi)
      real*8     Esf_Hid,Esf_VM,OI

      real*8   dv,EsMax,EsMin,theta
      real*8   E,nu
C
C     Inicializacion
C      
      DESP= 0.d0
      ESFUERZOS    = 0.d0
      DEFORMACIONES= 0.d0
      n            = 1
      def2D        = 0.d0
C     Definicion del tensor de constantes elasticas
      Dmat   = 0.d0
      Esf_Hid= 0.D0
      Esf_VM = 0.D0
C     
      do jelem=1,NELEMS
         Dmat2 = 0.d0
         do J=1,nnod
            x(1,J) = nodes(conectividades(jelem,J+1),1)
            x(2,J) = nodes(conectividades(jelem,J+1),2)
         enddo
C           Funcion de forma
         call f_forma2D(0.d0,0.d0,x,shp,xjac)
C           Se obtienen las matrices de calculo
         call Matrices_de_calculo2D(shp,Nm,Bm,Ne,Be)
C           Calculo de los esfuerzos
         call matriz_constantes_el(a,b,zita,U,ndofel,de,x,jelem,
     1     t,Dmat2)

         do J=1,nnod
            DESP(2*(J-1)+1) = resNod(conectividades(jelem,J+1),1) 
            DESP(2*(J-1)+2) = resNod(conectividades(jelem,J+1),2) 
         enddo
C
         DEFORMACIONES = MATMUL(Be(1:3, :),DESP) !-def2D
C
         ESFUERZOS = 0.d0
         ESFUERZOS = MATMUL(Dmat2,DEFORMACIONES) !-def2D)   
C           Centro y radio del esfuerzo (circulo de Mohr)
         radio  = sqrt(((ESFUERZOS(1) - ESFUERZOS(2))/2.d0)**2 + 
     &             ESFUERZOS(3)**2)
         centro = 0.5d0*(ESFUERZOS(1) + ESFUERZOS(2))
C
         EsMax = centro+radio
         EsMin = centro-radio
C
         if(tipo_def.eq.2)then
            E      = propiedades(grupoFisico(jelem,2),1)
            nu     = propiedades(grupoFisico(jelem,2),2)
            sigma_zz = nu*(ESFUERZOS(1) + ESFUERZOS(2))
C           Esfuerzos equivalentes de Von Mises
         else if(tipo_def.eq.1)then
            sigma_zz = 0.d0
         end if
C
         Esf_VM = sqrt(((EsMax-EsMin)**2+(EsMin-sigma_zz)**2
     1      +(sigma_zz-EsMax)**2)/2.d0)
C             
         Esf_Hid=(EsMax+EsMin+sigma_zz)/3.d0
C
         thao_oct=Esf_VM*sqrt(2.d0)/3.d0
C     
         OI =  thao_oct + kOI * Esf_Hid

C
         resElem(jelem, 1) = DEFORMACIONES(1) !E11
         resElem(jelem, 2) = DEFORMACIONES(2) !E22
         resElem(jelem, 3) = 0.0 !E33 (Not implemented)
         resElem(jelem, 4) = DEFORMACIONES(3) !E12
         resElem(jelem, 5) = ESFUERZOS(1) !S11
         resElem(jelem, 6) = ESFUERZOS(2) !S22
         resElem(jelem, 7) = sigma_zz !S33
         resElem(jelem, 8) = ESFUERZOS(3) !S12
         resElem(jelem, 9) = Esf_VM !S_Mises
         resElem(jelem, 10) = Esf_Hid !S_Hyd
         resElem(jelem, 11) = thao_oct !S_Oct
         resElem(jelem, 12) = OI !OI
         resElem(jelem, 13) = 0.0 !BG (Not implemented)
         resElem(jelem, 14) = OI !CMI
C     
      enddo
      return
      end
C
C------------------------------------------------------------------------------------------
C---------------------------------------------------------------------------oute-----------
C       
C       Funcion accumulateResults()
C
C
C       Funcion para acumular los resultados de los elementos
C       
C
C -----------------------------------------------------------------------------------------
C
C      Argumentos: KINC, paso
C
      subroutine accumulateResults(KINC)
C  
      use common_utils
      include 'ABA_PARAM.INC'
      ! include    'conec.for'
C
      integer    i,j,KINC

c     resElem(i,1) = E11
c     resElem(i,2) = E22
c     resElem(i,3) = E33 (Not implemented)
c     resElem(i,4) = E12
c     resElem(i,5) = S11
c     resElem(i,6) = S22
c     resElem(i,7) = S33
c     resElem(i,8) = S12
c     resElem(i,9) = S_Mises
c     resElem(i,10) = S_Hyd
c     resElem(i,11) = S_Oct
c     resElem(i,12) = OI
c     resElem(i,13) = BG (Not implemented)
c     resElem(i,14) = CMI

      do i=1,NUMNODE
         do j=1,nResNod
            if (KINC == 1) then
               cumulativeResNod(i,j) = resNod(i,j)
            else
               cumulativeResNod(i,j) = cumulativeResNod(i,j) + resNod(i,j)
            endif
            if (KINC==nLoads) then
               cumulativeResNod(i,j) = cumulativeResNod(i,j) / DBLE(nLoads)
            endif
         enddo
      enddo

      do i=1,NELEMS
         do j=1,nResElem
            if (KINC == 1) then
               cumulativeResElem(i,j) = resElem(i,j)
            else
               cumulativeResElem(i,j) = cumulativeResElem(i,j) + resElem(i,j)
            endif
            if (KINC==nLoads) then
               cumulativeResElem(i,j) = cumulativeResElem(i,j) / DBLE(nLoads)
            endif
         enddo

         if (KINC == 1) then
            resElem(i,13) = resElem(i,10) !min Hyd
            resElem(i,14) = resElem(i,11) !max Oct
         else
            if (resElem(i,10) < resElem(i,13)) then
               resElem(i,13) = resElem(i,10)
            end if

            if (resElem(i,11) > resElem(i,14)) then
               resElem(i,14) = resElem(i,11)
            end if
         end if

         resElem(i,13) = resElem(i,13)
         resElem(i,14) = resElem(i,14)
         resElem(i,15) = resElem(i,14) + kOI * resElem(i,13)

         if (KINC==nLoads) then
            cumulativeResElem(i,13) = resElem(i,13)
            cumulativeResElem(i,14) = resElem(i,14)
            cumulativeResElem(i,15) = resElem(i,15)
         end if
      enddo
      return
      end
C
C------------------------------------------------------------------------------------------
C---------------------------------------------------------------------------oute-----------
C       
C       Funcion calculateCMIthreshold()
C
C
C       Funcion para acumular los resultados de los elementos
C       
C
C -----------------------------------------------------------------------------------------
C
C      Argumentos: 
C
      subroutine calculateCMIthreshold()
C
      use common_utils
      include 'ABA_PARAM.INC'
      ! include    'conec.for'
C
      integer    i,j,nCart

      CMIavg = 0.0
      CMIstd = 0.0
      nCart = 0

      do i=1,NELEMS
         if (grupoFisico(i,1) == 2) then
            CMIavg = CMIavg + cumulativeResElem(i, 14)
            nCart = nCart + 1
         endif
      enddo

      CMIavg = CMIavg / DBLE(nCart)

      do i=1,NELEMS
         if (grupoFisico(i,1) == 2) then
            CMIstd = CMIstd + (cumulativeResElem(i, 14) - CMIavg)**2
         endif
      enddo

      CMIstd = sqrt(CMIstd / DBLE(nCart)) ! Standard deviation of CMI

      CMIthreshold = CMIavg + stdWeight * CMIstd

      return
      end
C------------------------------------------------------------------------------------------
C---------------------------------------------------------VECTOR DE DEFORMACION------------
C
C	 Funcion vector de deformacion impuesta
c
C------------------------------------------------------------------------------------------
	subroutine deformacion(chi,eta,zita,u,ndofel,de,x,jelem,
	1t,def2D)
C
      use common_utils
      ! include 'conec.for'
C
C     Variables de entrada
      integer  ndofel,jelem
      real*8   chi,eta,zita
      real*8   u(ndofel),x(dim,nnod),t!,de(8)
C
C     Variables de salida
      real*8   def2D(4)
C
C     Variables de la subrutina
C     Generales
      integer  i,j,k
C	real*8   u1(nnod),u2(nnod),u3(nnod),u4(nnod),u5(nnod)
      real*8   u1(nnod),u2(nnod)
C
C     Definicion del tensor de constantes elasticas
      Def2D   = 0.d0
C
C     Para el caso 2D
C      if(dim.eq.2)then
C         Def2D(2)=0.3d0
C	    end if
C
	return 
	end
C------------------------------------------------------------------------------
C----------------------------------------------------------------URDFIL--------
C
C     Rutina utilizada para leer los datos de salida y escribirlos
C     en el archivo para TECPLOT o MATLAB
C
C     URDFIL sirve para leer los datos a la salida del archivo
C
C------------------------------------------------------------------------------
      SUBROUTINE URDFIL(LSTOP,LOVRWRT,KSTEP,KINC,DTIME,TIME)
C
      use common_utils
      use iso_c_binding
      INCLUDE 'ABA_PARAM.INC'
      ! include 'conec.for'
C
      DIMENSION ARRAY(513),JRRAY(NPRECD,513),TIME(2)
      EQUIVALENCE (ARRAY(1),JRRAY(1,1))
C
      LOGICAL :: firstKey = .TRUE., lastKey = .FALSE.
C      INTEGER :: locID = 0
      INTEGER :: flagOutput = 0
      INTEGER :: prevFlagOutput = 0
      INTEGER :: i, j, k
      real*8 S11, S22, S12
      real*8 S1, S2, S_avg, R
C
      character*276         filename
      character*276         folderName
      character(256)        JOBDIR
      character(256)        JOBNAME
      character(21)         loadString
      character(21)         stepString
C
C FIND CURRENT INCREMENT.
C
      j = 0
      k = 0
      CALL POSFIL(KSTEP,KINC,ARRAY,JRCD)
C
      DO K1=1,999999
         CALL DBFILE(0,ARRAY,JRCD)
         IF (JRCD .NE. 0) GO TO 110
         KEY=JRRAY(1,2)
C
C        RECORD 201
C
         IF (KEY.EQ.201) THEN
            k = k+1
            resNod(k, 1) = ARRAY(4)
            resNod(k, 2) = ARRAY(5)
            resNod(k, 3) = ARRAY(6)
            resNod(k, 4) = ARRAY(7)
         END IF
C
      END DO
C
 110  CONTINUE

C     Cálculo de los esfuerzos y las deformaciones
      call outsigma()
C
C     Acumulación de los resultados, TODO: Revisar que lo hace en el ultimo paso
      call accumulateResults(KINC)
C
      call GETOUTDIR(JOBDIR,LENJOBDIR)
      call GETJOBNAME(JOBNAME,LENJOBNAME)
      write(loadString, '(I3.3)') KINC
      write(stepString, '(I3.3)') KSTEP
C     
      folderName = trim(jobdir) // '\resultados\' // trim(jobname) // '\' //
     & 'step' // trim(stepString) // '\'

      call execute_command_line('if not exist ' // trim(folderName) //
     & ' mkdir ' // trim(folderName))

      filename = trim(folderName) // trim(jobname) //
     &'_step' // trim(stepString) // '_load' // trim(loadString) // '.vtu'
C
C     Escritura de los resultados en el archivo VTK
      call writeVTKFile(filename, resNod, resElem)
C
C     Escritura de los resultados acumulados
      
      if (KINC==nLoads) then
         folderName = trim(jobdir) // '\resultados\' // trim(jobname)
         filename = trim(folderName) // '\' // trim(jobname) // '_step' // trim(stepString)//'.vtu'
C
         call writeVTKFile(filename, cumulativeResNod, cumulativeResElem)
      end if
C
      return
      end
C-------------------------------------------------------------------------------------------
C-------------------------------------------------------------writeVTKFile-------------------
C
C     Rutina para escribir datos en VTK
C
C     
C
C-------------------------------------------------------------------------------------------
      subroutine writeVTKFile(filename, mNod, mElem)
C     
      use common_utils
      ! include 'conec.for'
C
      character*276         filename
      real*8 mNod(NUMNODE, nResNod), mElem(NELEMS, nResElem) ! Matrices de nodos y elementos, mismo tamaño que matrices de resultados
C
      open(UNIT=16,file=filename,action='write',status='UNKNOWN')
      write(16,'(a73)') '<VTKFile type="UnstructuredGrid" version="1,0" byte_order="LittleEndian">'
      write(16,'(a18)') '<UnstructuredGrid>'
      write(16,'(a23,i0,a17,i0,a2)') '<Piece NumberOfPoints="', NUMNODE, 
     & '" NumberOfCells="', NELEMS, '">'
      write(16,'(a8)') '<Points>'
      write(16,'(a64)') '<DataArray type="Float64" NumberOfComponents="3" format="ascii">'
C
      do i=1,NUMNODE
         write(16,'(3(E20.13,1X))') (nodes(i,j),j=1,2),0.0d0
      enddo 
C
      write(16,'(a12)') '</DataArray>'
      write(16,'(a9)') '</Points>'
C     
      write(16,'(a7)') '<Cells>'
      write(16,'(a59)')'<DataArray type="Int64" Name="connectivity" format="ascii">'
C
      do i=1,NELEMS
         write(16,'(I0,1X,I0,1X,I0,1X,I0)') (conectividades(i,j)-1,j=2,nnod+1)
      end do
C
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a54)')'<DataArray type="Int64" Name="offsets" format="ascii">'
      DO i = 1, NELEMS
         write(16,'(I0)') i*4
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a52)')'<DataArray type="Int64" Name="types" format="ascii">'
      DO i = 1, NELEMS
         write(16,'(I0)') 9
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a8)') '</Cells>'
C
      write(16,'(a30)') '<PointData Vectors="''U'', ''C''">'
      write(16,'(a58,a55)') '<DataArray type="Float64" Name="U" NumberOfComponents="2" ',
     & 'ComponentName0="U1" ComponentName1="U2" format="ascii">'
      DO i=1,NUMNODE
         write(16,'(2(E20.13,1X))') mNod(i, 1), mNod(i, 2)
      END DO
      write(16,'(a12)') '</DataArray>'

      write(16,'(a58,a55)') '<DataArray type="Float64" Name="C" NumberOfComponents="2" ',
     & 'ComponentName0="C1" ComponentName1="C2" format="ascii">'
      DO i=1,NUMNODE
         write(16,'(2(E20.13,1X))') mNod(i, 3), mNod(i, 4)
      END DO
      write(16,'(a12)') '</DataArray>'
      write(16,'(a12)') '</PointData>'
C
      write(16,'(a64,a40)') '<DataArray type="Float64" Name="Load" NumberOfComponents="1" ',
     & 'ComponentName0="Load" format="ascii">'
      DO i=1,NELEMS
         write(16,'(1(E20.13,1X))') mElem(i, 9)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16, '(a46)') '<CellData Tensors="''E_Centroid'',''S_Centroid''">'
      write(16,'(a67,a99)') '<DataArray type="Float64" Name="E_Centroid" NumberOfComponents="4" ',
     & 'ComponentName0="E11" ComponentName1="E22" ComponentName2="E33" ComponentName3="E12" format="ascii">'
      DO i=1,NELEMS
         write(16,'(4(E20.13,1X))') mElem(i, 1), mElem(i, 2),
     & mElem(i, 3), mElem(i, 4)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a67,a99)') '<DataArray type="Float64" Name="S_Centroid" NumberOfComponents="4" ',
     & 'ComponentName0="S11" ComponentName1="S22" ComponentName2="S33" ComponentName3="S12" format="ascii">'
      DO i=1,NELEMS
         write(16,'(4(E20.13,1X))') mElem(i, 5), mElem(i, 6),
     &         mElem(i, 7), mElem(i, 8)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a64,a40)') '<DataArray type="Float64" Name="S_Mises" NumberOfComponents="1" ',
     & 'ComponentName0="S_Mises" format="ascii">'
      DO i=1,NELEMS
         write(16,'(1(E20.13,1X))') mElem(i, 9)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a62,a38)') '<DataArray type="Float64" Name="S_Hyd" NumberOfComponents="1" ',
     & 'ComponentName0="S_Hyd" format="ascii">'
      DO i=1,NELEMS
         write(16,'(1(E20.13,1X))') mElem(i, 10)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a62,a38)') '<DataArray type="Float64" Name="S_Oct" NumberOfComponents="1" ',
     & 'ComponentName0="S_Oct" format="ascii">'
      DO i=1,NELEMS
         write(16,'(1(E20.13,1X))') mElem(i, 11)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a59,a35)') '<DataArray type="Float64" Name="OI" NumberOfComponents="1" ',
     & 'ComponentName0="OI" format="ascii">'
      DO i=1,NELEMS
         write(16,'(1(E20.13,1X))') mElem(i, 12)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a63,a35)') '<DataArray type="Float64" Name="minHyd" NumberOfComponents="1" ',
     & 'ComponentName0="CK" format="ascii">'
      DO i=1,NELEMS
         write(16,'(1(E20.13,1X))') mElem(i, 13)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a63,a35)') '<DataArray type="Float64" Name="maxOct" NumberOfComponents="1" ',
     & 'ComponentName0="CK" format="ascii">'
      DO i=1,NELEMS
         write(16,'(1(E20.13,1X))') mElem(i, 14)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a62,a35)') '<DataArray type="Float64" Name="OI_2" NumberOfComponents="1" ',
     & 'ComponentName0="CK" format="ascii">'
      DO i=1,NELEMS
         write(16,'(1(E20.13,1X))') mElem(i, 15)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a53)') '<DataArray type="Int32" Name="Region" format="ascii">'
      DO i=1,NELEMS
         write(16,*) grupoFisico(i, 2)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a11)') '</CellData>'
      write(16,'(a8)') '</Piece>'
      write(16,'(a19)') '</UnstructuredGrid>'
      write(16,'(a10)') '</VTKFile>'
      close(16)
C
      RETURN
      END
C-------------------------------------------------------------------------------------------
C-------------------------------------------------------------UEXTERNALDB-------------------
C
C     Rutina para leer condiciones fuente y condiciones de archivos externos
C
C     Abre y cierra los archivos necesarios para el calculo
C
C-------------------------------------------------------------------------------------------
C
      subroutine UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      use common_utils
      include 'ABA_PARAM.INC'
      ! include 'conec.for'
C
      logical            :: Searstr
      character(256)       :: faceString ! U letter for compatbility with *DLOAD
      character(256)       :: loadString ! up to 9 loads
      character(25)      :: stepString ! For "*Step, name=LoadStep1" detection
      character(256)     :: wholeLine
      character(256)        JOBDIR
      character(256)        JOBNAME
      character*276         filename
      character*276         folderName
      character*276         line, key
      integer               ios, pos
      integer               i,j,k,elem

C      integer NUMNODE, NELEMS, dim, nnod
C      integer a2e, a3e, a4e, be0
C      integer nLoads, nResNod, nResElem
C      integer maxNElementLoads
C      integer numProps, numMats
C      integer axi, tipo_def
C      integer filasContorno1
C      integer filasContorno2
C      integer velocidad


      print *, '######################################################'
      print *, ''
      print *, '###..Ejecutando UEXTERNALDB'
      print *, '###..LOP = ', LOP
      print *, '###..KSTEP = ', KSTEP
      print *, '###..KINC = ', KINC
      print *, ''
      print *, '######################################################'
C     Variables llamadas al comienzo del análisis
      if (LOP.eq.0) then
C        Extracción de la información de los archivos
         call GETOUTDIR(JOBDIR,LENJOBDIR)
         call GETJOBNAME(JOBNAME,LENJOBNAME)
C        Llenado de la matriz de avance

         elem = 2*a3e*be0 + 1
         
         do j = a2e+1 , a2e+a3e
            do i = 1 , a4e
               advanceElements(i,j) = elem
               elem = elem + 1
            enddo
         enddo
         
         do j = 1 , a2e
            do i = 1 , a4e
               advanceElements(i,j) = elem
               elem = elem + 1
            enddo
         enddo
            
         do j = a2e + a3e + 1, a2e + a3e + a2e
            do i = 1 , a4e
               advanceElements(i,j) = elem
               elem = elem + 1
            enddo
         enddo
C-------------------------------------
C        Llamada al archivo de parametros (reemplazo de conec.for)
         filename=trim(jobdir)//'\parametros.txt'
         call test_get_variable_from_file(filename)
         call initialize_arrays()
C-------------------------------------
C        Llamada al archivo de propiedades
         filename=' '
         filename(1:lenjobdir)=jobdir(1:lenjobdir)
         filename(lenjobdir+1:lenjobdir+17)='/propiedades.txt'
C
         open(UNIT=14,file=filename(1:lenjobdir+17), status='old')
            if (Searstr(14,'*UEL PROPERTY')) then
               READ(14,*)((propiedades(i,j),j=1,numProps),i=1,numMats)
            else
               stop '###..Error en lectura'
            end if
         close(14)
C-------------------------------------
C        Llamada al archivo de grupos físicos
         filename=' '
         filename(1:lenjobdir)=jobdir(1:lenjobdir)
         filename(lenjobdir+1:lenjobdir+19)='/gruposFisicos.txt'
C
         open(UNIT=14,file=filename(1:lenjobdir+19), status='old')
            if (Searstr(14,'Element Tag, Physical Group Tag')) then
               READ(14,*)((grupoFisico(i,j),j=1,2),i=1,NELEMS)
            else
               stop '###..Error en lectura'
            end if
         close(14)
C
C-------------------------------------
C        Llamada al archivo conectividades.inp
         filename=' '
         filename(1:lenjobdir)=jobdir(1:lenjobdir)
         filename(lenjobdir+1:lenjobdir+19)='/conectividades.inp'
C
         open(UNIT=15,file=filename(1:lenjobdir+19), status='old')
            if (Searstr(15,'*ELEMENT,TYPE=U1,ELSET=UEL')) then
               READ(15,*)((conectividades(i,j),j=1,nnod + 1),i=1,NELEMS)
            else
               stop '###..Error en lectura'
            end if
         close(15)
C
C-------------------------------------
C        Llamada al archivo nodos.inp

         filename=' '
	      filename(1:lenjobdir)=jobdir(1:lenjobdir)
         filename(lenjobdir+1:lenjobdir+10)='/nodos.inp'
C
         open(UNIT=16,file=filename(1:lenjobdir+10), status='old')
            if (Searstr(16,'*NODE,NSET=N2')) then
               READ(16,*) (k,(nodes(i,j),j=1,dim),i=1,NUMNODE)
            else
               stop '###..Error en lectura'
            end if
         close(16)
C-------------------------------------
C       Llamada al archivo contorno.inp

         filename=' '
         filename(1:lenjobdir)=jobdir(1:lenjobdir)
         filename(lenjobdir+1:lenjobdir+14)='/contorno.inp'
C
         open(UNIT=15,file=filename(1:lenjobdir+14), status='old')
            if (Searstr(15,'*NSET,NSET=contorno1')) then
               READ(15,*)((contorno1(i,j),j=1,6),i=1,filascontorno1)
            else
               stop '###..Error en lectura'
            end if
         close(15)
C
         filename=' '
         filename(1:lenjobdir)=jobdir(1:lenjobdir)
         filename(lenjobdir+1:lenjobdir+14)='/contorno.inp'
C
         open(UNIT=15,file=filename(1:lenjobdir+14), status='old')
            if (Searstr(15,'*NSET,NSET=contorno2')) then
               READ(15,*)((contorno2(i,j),j=1,6),i=1,filascontorno2)
            else
               stop '###..Error en lectura'
            end if
         close(15)
C        listNElementLoads
         filename=' '
         filename(1:lenjobdir)=jobdir(1:lenjobdir)
         filename(lenjobdir+1:lenjobdir+18)='/nFilasCargas.txt'
C
         open(UNIT=15,file=filename(1:lenjobdir+18), status='old')
            READ(15,*)(listNElementLoads(i),i=1,nLoads)
         close(15)

         do i=1,nLoads
            print*, listNElementLoads(i)
         enddo
C     Archivos llamados al comienzo de cada paso
      else if (LOP.eq.1) then
C     Reinicialización de variables
         cambioGrupoFisico = .false.
C-------------------------------------
C        Extracción de la información de los archivos
         call GETOUTDIR(JOBDIR,LENJOBDIR)
C        Llamada al archivo de cargas
         write(loadString, *) KINC
         stepString = trim('*Step, name=LoadStep' // trim(adjustl(loadString)))
C-------------------------------------
C        Llamada al archivo carga.inp
         filename=' '
         filename(1:lenjobdir)=jobdir(1:lenjobdir)
         filename(lenjobdir+1:lenjobdir+11)='/carga.inp'
C
         open(UNIT=15,file=filename(1:lenjobdir+11), status='old')
         if (Searstr(15,stepString)) then
            ! Number of element loads per load, starts 2 lines after step line
            read(15, '(A)') wholeLine
            do i=1,listNElementLoads(KINC)
               read(15, '(I, I, F)') elementFaces(i, 1),
     1            elementFaces(i, 2), elementLoads(i)
            enddo
         else
            stop '###..Error en lectura'
         end if
        close(15)
      end if
      return
      end