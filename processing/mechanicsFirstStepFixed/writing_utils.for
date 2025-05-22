      module writing_utils
      implicit none

      contains
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------writeVTKFile-------------------
!
!     Rutina para escribir datos en VTK
!
!     
!
!-------------------------------------------------------------------------------------------
      subroutine writeVTKFile(filename, mNod, mElem)
!  
      use reading_utils
      character*276         filename
      integer i,j
      real*8 mNod(NUMNODE, nResNod), mElem(NELEMS, nResElem) ! Matrices de nodos y elementos, mismo tamaño que matrices de resultados

      open(UNIT=16,file=filename,action='write',status='UNKNOWN')
      write(16,'(a73)') '<VTKFile type="UnstructuredGrid" version="1,0" byte_order="LittleEndian">'
      write(16,'(a18)') '<UnstructuredGrid>'
      write(16,'(a23,i0,a17,i0,a2)') '<Piece NumberOfPoints="', NUMNODE,'" NumberOfCells="', NELEMS, '">'

      write(16,'(a8)') '<Points>'
      call writeNodes()
      write(16,'(a9)') '</Points>'
!     
      write(16,'(a7)') '<Cells>'
      call writeDataArrayIntNoComp('connectivity',conectividades(:,2:nnod+1)-1, NELEMS, nnod)
      call writeDataArrayIntNoComp('offsets',[(i * nnod, i = 1, NELEMS)], NELEMS, 1)
      call writeDataArrayIntNoComp('types',[(3 * (nnod - 1), i = 1, NELEMS)], NELEMS, 1)
      write(16,'(a8)') '</Cells>'
!     
!     PointData
      write(16,'(a30)') '<PointData Vectors="''U'', ''C''">'
      call writeDataArrayReal('U', mNod(:, 1:dim), NUMNODE, dim)
      call writeDataArrayReal('C', mNod(:, 1+dim:ndofdiff+dim), NUMNODE, ndofdiff)
      call writeDataArrayInt('Region', grupoFisicoN(:, 2), NUMNODE, 1)
      call writeDataArrayReal('BorderNode', borderVectorNod, NUMNODE, 1)
      write(16,'(a12)') '</PointData>'

!     CellData
      write(16, '(a46)') '<CellData Tensors="''E_Centroid'',''S_Centroid''">'
      if (dim == 2) then
         call writeDataArrayReal('E_Centroid', mElem(:, 1:4), NELEMS, 4)
         call writeDataArrayReal('S_Centroid', mElem(:, 5:8), NELEMS, 4)
         call writeDataArrayReal('S_Mises', mElem(:, 9), NELEMS, 1)
         call writeDataArrayReal('S_Hyd', mElem(:, 10), NELEMS, 1)
         call writeDataArrayReal('S_Oct', mElem(:, 11), NELEMS, 1)
         call writeDataArrayReal('MG', mElem(:, 12), NELEMS, 1)
         call writeDataArrayReal('C', mElem(:, 13:17), NELEMS, ndofdiff)
         ! call writeDataArrayReal('BG', mElem(:, 15), NELEMS, 1)
         ! call writeDataArrayReal('CMI', mElem(:, 16), NELEMS, 1)
         call writeDataArrayReal('BorderElement', borderVectorElem, NELEMS, 1)
      end if
      call writeDataArrayInt('Region', grupoFisico(:, 2), NELEMS, 1)
      call writeProps()
      write(16,'(a11)') '</CellData>'

      write(16,'(a8)') '</Piece>'
      write(16,'(a19)') '</UnstructuredGrid>'
      write(16,'(a10)') '</VTKFile>'
      close(16)
!
      RETURN
      END subroutine writeVTKFile
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------writeVTKFilePre-------------------
!
!     Rutina para escribir datos en VTK
!
!     
!
!-------------------------------------------------------------------------------------------
      subroutine writeVTKFilePre(filename, mNod, mElem)
!  
      use reading_utils
      character*276         filename
      integer i,j
      real*8 mNod(NUMNODE, nResNod), mElem(NELEMS, nResElem) ! Matrices de nodos y elementos, mismo tamaño que matrices de resultados

      open(UNIT=16,file=filename,action='write',status='UNKNOWN')
      write(16,'(a73)') '<VTKFile type="UnstructuredGrid" version="1,0" byte_order="LittleEndian">'
      write(16,'(a18)') '<UnstructuredGrid>'
      write(16,'(a23,i0,a17,i0,a2)') '<Piece NumberOfPoints="', NUMNODE,'" NumberOfCells="', NELEMS, '">'

      write(16,'(a8)') '<Points>'
      call writeNodes()
      write(16,'(a9)') '</Points>'
!     
      write(16,'(a7)') '<Cells>'
      call writeDataArrayIntNoComp('connectivity',conectividades(:,2:nnod+1)-1, NELEMS, nnod)
      call writeDataArrayIntNoComp('offsets',[(i * nnod, i = 1, NELEMS)], NELEMS, 1)
      call writeDataArrayIntNoComp('types',[(3 * (nnod - 1), i = 1, NELEMS)], NELEMS, 1)
      write(16,'(a8)') '</Cells>'
!     
!     PointData
      write(16,'(a30)') '<PointData Vectors="''U'', ''C''">'
      call writeDataArrayInt('Region', grupoFisicoN(:, 2), NUMNODE, 1)
      call writeDataArrayReal('BorderNode', borderVectorNod, NUMNODE, 1)
      write(16,'(a12)') '</PointData>'

!     CellData
      write(16, '(a46)') '<CellData Tensors="''E_Centroid'',''S_Centroid''">'
      if (dim == 2) then
         call writeDataArrayReal('BorderElement', borderVectorElem, NELEMS, 1)
      end if
      call writeDataArrayInt('Region', grupoFisico(:, 2), NELEMS, 1)
      call writeProps()
      write(16,'(a11)') '</CellData>'

      write(16,'(a8)') '</Piece>'
      write(16,'(a19)') '</UnstructuredGrid>'
      write(16,'(a10)') '</VTKFile>'
      close(16)
!
      RETURN
      END subroutine writeVTKFilePre
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------writeDataArrayReal-------------------
!
!     Rutina para escribir datos en VTK
!
!-------------------------------------------------------------------------------------------
!     TODO: Read numComp from size
      subroutine writeDataArrayReal(name, dataArray, numElements, numComp)
      
      character(len=*) name
      integer numElements, numComp
      real*8 dataArray(numElements, numComp)
      integer i,j
      character(len=256) componentString

      componentString = getComponentString(name, numComp)

      write(16,'(a)')  '<DataArray type="Float64" Name="'//trim(name)//'"'//trim(componentString)//' format="ascii">'
   
      DO i=1,numElements
         DO j=1,numComp
            write(16,'(ES21.13E3)') dataArray(i,j)
         END DO
      END DO
      write(16,'(a)') ' </DataArray>'

      end subroutine writeDataArrayReal
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------writeDataArrayReal-------------------
!
!     Rutina para escribir datos en VTK
!
!-------------------------------------------------------------------------------------------
!     TODO: Read numComp from size
      subroutine writeDataArrayInt(name, dataArray, numElements, numComp)
      
      character(len=*) name
      integer numElements, numComp
      integer dataArray(numElements, numComp)
      integer i,j
      character(len=256) componentString
      
      componentString = getComponentString(name, numComp)

      write(16,'(a)')  '<DataArray type="Int64" Name="'//trim(name)//'"'// trim(componentString)//' format="ascii">'
      DO i=1,numElements
         DO j=1,numComp
            write(16,'(I10,1X)', advance="no") dataArray(i,j)
         END DO
         write(16, *)
      END DO
      write(16,'(a)') '</DataArray>'

      end subroutine writeDataArrayInt
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------writeNodes-------------------
!
!     Rutina para escribir datos en VTK
!
!-------------------------------------------------------------------------------------------
      subroutine writeNodes()
      
         use reading_utils
         integer i,j
   
         write(16,'(a)')  '<DataArray type="Float64" NumberOfComponents="3" format="ascii">'
         
         DO i=1,NUMNODE
            write(16,'(3E20.13,1X)') (nodes(i,j), j=1,dim), (0.d0, j=dim+1,3)
         END DO
         write(16, *)
         write(16,'(a)') ' </DataArray>'
   
         end subroutine writeNodes
!-------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------writeProps---------
!
      subroutine writeProps()
   
      use reading_utils
      integer i,j
      real*8 E, nu, D, f1, f2, g1, g2, h1, h2

      write(16,'(a)')  '<DataArray type="Float64" Name="props" NumberOfComponents="9"'
      write(16,'(a)') ' ComponentName0="E" ComponentName1="nu" ComponentName2="D" ComponentName3="f1"'
      write(16,'(a)') ' ComponentName4="f2" ComponentName5="g1" ComponentName6="g2" ComponentName7="h1"'
      write(16,'(a)') ' ComponentName8="h2" format="ascii">'
   
      DO i=1,NELEMS
         E     = propiedades(grupoFisico(i,2)+1,1)
         nu    = propiedades(grupoFisico(i,2)+1,2)
         D     = propiedades(grupoFisico(i,2)+1,3)
         f1    = propiedades(grupoFisico(i,2)+1,4)
         f2    = propiedades(grupoFisico(i,2)+1,5)
         g1    = propiedades(grupoFisico(i,2)+1,6)
         g2    = propiedades(grupoFisico(i,2)+1,7)
         h1    = propiedades(grupoFisico(i,2)+1,8)
         h2    = propiedades(grupoFisico(i,2)+1,9)

         write(16,'(E20.13,1X)', advance="no") E, nu, D, f1, f2, g1, g2, h1, h2
      ENDDO

      write(16, *)
      write(16,'(a)') ' </DataArray>'

      end subroutine writeProps
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------writeDataArrayReal-------------------
!
!     Rutina para escribir datos en VTK
!
!-------------------------------------------------------------------------------------------
!     TODO: Read numComp from size
      subroutine writeDataArrayIntNoComp(name, dataArray, numElements, numComp)

      character(len=*) name
      integer numElements, numComp
      integer dataArray(numElements, numComp)
      integer i,j
      character(len=256) componentString

      write(16,'(a)')  '<DataArray type="Int64" Name="'//trim(name)//'" format="ascii">'
      DO i=1,numElements
         DO j=1,numComp
            write(16,'(I10,1X)', advance="no") dataArray(i,j)
         END DO
         write(16, *)
      END DO
      write(16,'(a)') '</DataArray>'

      end subroutine writeDataArrayIntNoComp
!--------------------------------------------------------------------------------------------
!-------------------------------------------------------------itoa-------------------
!
!
!     Rutina para convertir un entero a cadena de caracteres
!
!
!-------------------------------------------------------------------------------------------
      function itoa(num) result(str)
         integer, intent(in) :: num
         character(len=32) :: str  ! Adjust length as needed
         write(str, '(I0)') num    ! Convert integer to string
      end function itoa
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------getComponentString-------------------
!
!
!     Rutina para crear el string de componentes
!
!
!--------------------------------------------------------------------------------------------
      function getComponentString(name, numComp) result(componentString)
         character(len=*) name
         integer numComp
         integer i,j
         character(len=256) componentString, tempString, itoaString

         componentString = " NumberOfComponents=""" // trim(adjustl(itoa(numComp))) // """"

         DO i = 1, numComp
            write(tempString, '(a,i0,a)') "ComponentName", i - 1, "="

            if (numComp == 1) then
               itoaString = trim(name)
            elseif (numComp == 4) then
               if (i == 4) then
                  itoaString = trim(name) // "12"
               else
                  itoaString = trim(name) // trim(adjustl(itoa(i))) // trim(adjustl(itoa(i)))
               end if
            else
               itoaString = trim(name) // trim(adjustl(itoa(i)))
            end if

            itoaString =  """" // trim(itoaString) // """"
            componentString = trim(componentString) // " " // trim(tempString) // itoaString
         END DO
      end function getComponentString

      end module writing_utils
!------------------------------------------------------------------------------
!----------------------------------------------------------------URDFIL--------
!
!     Rutina utilizada para leer los datos de salida y escribirlos
!     en el archivo para TECPLOT o MATLAB
!
!     URDFIL sirve para leer los datos a la salida del archivo
!
!------------------------------------------------------------------------------
      SUBROUTINE URDFIL(LSTOP,LOVRWRT,KSTEP,KINC,DTIME,TIME)
!
      use debug_utils
      use reading_utils
      use writing_utils
      use iso_c_binding
      INCLUDE 'ABA_PARAM.INC'
!
!
      DIMENSION ARRAY(513),JRRAY(NPRECD,513),TIME(2)
      EQUIVALENCE (ARRAY(1),JRRAY(1,1))
!
      INTEGER :: i, j, k, k1
!
      character*276         filename
      character*276         folderName
      character(256)        JOBDIR
      character(256)        JOBNAME
      character(21)         incString
      character(21)         stepString
!
! FIND CURRENT INCREMENT.
!
      call printStepInfo(-1, KSTEP, KINC, 'URDFIL     : ')
      k = 0
      CALL POSFIL(KSTEP,KINC,ARRAY,JRCD)
!
      DO K1=1,999999
         CALL DBFILE(0,ARRAY,JRCD)
         IF (JRCD .NE. 0) GO TO 110
         KEY=JRRAY(1,2)
!
!        RECORD 201
!
         IF (KEY.EQ.201) THEN
            k = k+1
            do i=1,dim+ndofdiff
               resNod(k, i) = ARRAY(3+i)
            end do
         END IF
!
      END DO
!
 110  CONTINUE
      call GETOUTDIR(JOBDIR,LENJOBDIR)
      call GETJOBNAME(JOBNAME,LENJOBNAME)
!     Cálculo de los esfuerzos y las deformaciones
      call outsigma()
      if (KSTEP.eq.1) then
         call accumulateResults(KINC)
      end if
      write(incString, '(I3.3)') KINC
      write(stepString, '(I3.3)') KSTEP
      folderName = trim(jobdir) // '\resultados\' // trim(jobname) // '\' // 'step' // trim(stepString) // '\'
      call execute_command_line('if not exist ' // trim(folderName) // ' mkdir ' // trim(folderName))
      filename = trim(folderName) // trim(jobname) // '_step' // trim(stepString) // '_inc' // trim(incString) // '.vtu'
!
!     Escritura de los resultados en el archivo VTK
      call writeVTKFile(filename, resNod, resElem)
!
!     Ultimo incremento del paso
!    
      if (KSTEP == 1 .and. KINC == nLoads) then
         folderName = trim(jobdir) // '\resultados\' // trim(jobname)
         filename = trim(folderName) // '\' // trim(jobname) // '_step' // trim(stepString)//'.vtu'
!        
!        Escritura de resultados
         call writeVTKFile(filename, cumulativeResNod, cumulativeResElem)
      end if
      return
      end