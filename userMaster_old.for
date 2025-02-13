      SUBROUTINE DLOAD(F,KSTEP,KINC,TIME,NOEL,NPT,LAYER,KSPT,
     1 COORDS,JLTYP,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TIME(2), COORDS (3)
      CHARACTER*80 SNAME
C
      RETURN
      END

      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
      include 'fortVars.for'
C
      DIMENSION TIME(2)
C
      logical            :: Searstr
      character(256)        JOBDIR
      character*276         filename
      integer               i,j
C
      if (LOP.eq.0) then
C       Llamada al archivo de grupos físicos
            call GETOUTDIR(JOBDIR,LENJOBDIR)
            filename=' '
            filename(1:lenjobdir)=jobdir(1:lenjobdir)
            filename(lenjobdir+1:lenjobdir+19)='/gruposFisicos.txt'
            open(UNIT=14,file=filename(1:lenjobdir+19), status='old')
C
            if (Searstr (14,'Element Tag, Physical Group Tag')) then
            READ(14,*)((grupoFisico(i,j),j=1,2),i=1,NELEMS)
            else
            stop '###..Error en lectura'
            end if
            close(14)
C
C       Llamada al archivo de conectividades 
            call GETOUTDIR(JOBDIR,LENJOBDIR)
            filename=' '
            filename(1:lenjobdir)=jobdir(1:lenjobdir)
            filename(lenjobdir+1:lenjobdir+19)='/conectividades.inp'
            open(UNIT=15,file=filename(1:lenjobdir+19), status='old')
C
            if (Searstr (15,'*ELEMENT, TYPE=CPE4, ELSET=UEL')) then
            READ(15,*)((conectividades(i,j),j=1,nnod + 1),i=1,NELEMS)
            else
            stop '###..Error en lectura'
            end if
            close(15)
C
C       Llamada al archivo de nodos
            filename=' '
            filename(1:lenjobdir)=jobdir(1:lenjobdir)
            filename(lenjobdir+1:lenjobdir+10)='/nodos.inp'
            open(UNIT=16,file=filename(1:lenjobdir+10), status='old')
C
            if (Searstr (16,'*NODE, NSET=N2')) then
                  READ(16,*) (k,(nodes(i,j),j=1,dim),i=1,NUMNODE)
            else
            stop '###..Error en lectura'
            end if
            close(16)
C
      ENDIF
      RETURN
      END

      SUBROUTINE URDFIL(LSTOP,LOVRWRT,KSTEP,KINC,DTIME,TIME)
C
      INCLUDE 'ABA_PARAM.INC'
      include 'fortVars.for'
C
      DIMENSION ARRAY(513),JRRAY(NPRECD,513),TIME(2)
      EQUIVALENCE (ARRAY(1),JRRAY(1,1))
      LOGICAL :: firstKey = .TRUE., lastKey = .FALSE.
      INTEGER :: locID = 0
      INTEGER :: flagOutput = 0
      INTEGER :: prevFlagOutput = 0
      REAL :: resNod(NUMNODE, 2)
      REAL :: resElem(NELEMS, 12)
      INTEGER :: i, j, k
      REAL :: S11, S22, S12
      REAL :: S1, S2, S_avg, R
C
      character*276         filename
      character(256)        JOBDIR
      character(256)        JOBNAME
C
      call GETOUTDIR(JOBDIR,LENJOBDIR)
      call GETJOBNAME(JOBNAME,LENJOBNAME)
      filename=' '
      filename(1:lenjobdir)=jobdir(1:lenjobdir)
      filename(lenjobdir+1:lenjobdir+1)='\'
      filename(lenjobdir+2:lenjobdir+lenjobname+1)=jobname(1:lenjobname)
      filename(lenjobdir+lenjobname+2:lenjobdir+lenjobname+6)='.vtu'
C
      open(UNIT=16,file=filename,action='write',status='replace')
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
C RECORD 1 CONTAINS VALUES FOR Element header record
C
            IF (KEY.EQ.1) THEN
C                  WRITE(*,*) 'Element number = ',JRRAY(1,3)
C                  WRITE(*,*) 'Integration point = ',JRRAY(1,4)
                  j = j + 1
                  locID = JRRAY(1,6)
            END IF

            IF (KEY .EQ. 1911) THEN

            END IF

            IF (KEY .EQ. 2001) THEN

            END IF
C
C RECORD 101 CONTAINS VALUES FOR U
C           
            IF (KEY.EQ.101 .AND. locID .EQ. 1) THEN
                  k = k+1
                  resNod(k, 1) = ARRAY(4)
                  resNod(k, 2) = ARRAY(5)

C                  WRITE(16,'(2(E20.13,1X))') ARRAY(4), ARRAY(5)
            END IF
C
C RECORD 21 CONTAINS VALUES FOR E
C
            IF (KEY.EQ.21 .AND. locID .EQ. 1) THEN
                  resElem(j, 1) = ARRAY(3)
                  resElem(j, 2) = ARRAY(4)
                  resElem(j, 3) = ARRAY(5)
                  resElem(j, 4) = ARRAY(6)
C                  WRITE(16,'(4(E20.13,1X))') ARRAY(3), ARRAY(4), 0.0d0, ARRAY(5)
            END IF
C
C RECORD 11 CONTAINS VALUES FOR S
C
            IF (KEY.EQ.11 .AND. locID .EQ. 1) THEN
                  resElem(j, 5) = ARRAY(3)
                  resElem(j, 6) = ARRAY(4)
                  resElem(j, 7) = ARRAY(5)
                  resElem(j, 8) = ARRAY(6)
C                  WRITE(*,'(2F8.4)') ARRAY(3), ARRAY(4)
            END IF
C
C RECORD 12 CONTAINS VALUES FOR SINV
C
            IF (KEY.EQ.12 .AND. locID .EQ. 1) THEN
                  resElem(j, 9) = ARRAY(3)
C                  WRITE(*,'(2F8.4)'), ARRAY(3)
            END IF
C
      END DO
C
 110  CONTINUE
C     
C     Calculo de indice osteogénico
C
      do i=1,NELEMS

      S11 = resElem(i, 5)
      S22 = resElem(i, 6)
      S33 = resElem(i, 7)
      S12 = resElem(i, 8)

      resElem(i, 10) = - (S11 + S22 + S33) / 3.d0
      resElem(i, 11) = sqrt((S11-S22)**2+(S22-S33)**2+(S33-S11)**2 + 6*S12**2)/3.d0
      resElem(i, 12) = resElem(i, 11) + k_OI * resElem(i, 10)

      end do

      write(16,'(a25)') '<PointData Vectors="''U''">'
      write(16,'(a58,a55)') '<DataArray type="Float32" Name="U" NumberOfComponents="2" ',
     & 'ComponentName0="U1" ComponentName1="U2" format="ascii">'
      DO i=1,NUMNODE
            write(16,'(2(E20.13,1X))') resNod(i, 1), resNod(i, 2)
      END DO
      write(16,'(a12)') '</DataArray>'
      write(16,'(a12)') '</PointData>'
C
      write(16,'(a64,a40)') '<DataArray type="Float32" Name="Load" NumberOfComponents="1" ',
     & 'ComponentName0="Load" format="ascii">'
      DO i=1,NELEMS
            write(16,'(1(E20.13,1X))') resElem(i, 9)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16, '(a46)') '<CellData Tensors="''E_Centroid'',''S_Centroid''">'
      write(16,'(a67,a99)') '<DataArray type="Float32" Name="E_Centroid" NumberOfComponents="4" ',
     & 'ComponentName0="E11" ComponentName1="E22" ComponentName2="E33" ComponentName3="E12" format="ascii">'
      DO i=1,NELEMS
            write(16,'(4(E20.13,1X))') resElem(i, 1), resElem(i, 2), resElem(i, 3), resElem(i, 4)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a67,a99)') '<DataArray type="Float32" Name="S_Centroid" NumberOfComponents="4" ',
     & 'ComponentName0="S11" ComponentName1="S22" ComponentName2="S33" ComponentName3="S12" format="ascii">'
      DO i=1,NELEMS
            write(16,'(4(E20.13,1X))') resElem(i, 5), resElem(i, 6), resElem(i, 7), resElem(i, 8)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a64,a40)') '<DataArray type="Float32" Name="S_Mises" NumberOfComponents="1" ',
     & 'ComponentName0="S_Mises" format="ascii">'
      DO i=1,NELEMS
            write(16,'(1(E20.13,1X))') resElem(i, 9)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a62,a38)') '<DataArray type="Float32" Name="S_Hyd" NumberOfComponents="1" ',
     & 'ComponentName0="S_Hyd" format="ascii">'
      DO i=1,NELEMS
            write(16,'(1(E20.13,1X))') resElem(i, 10)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a62,a38)') '<DataArray type="Float32" Name="S_Oct" NumberOfComponents="1" ',
     & 'ComponentName0="S_Oct" format="ascii">'
      DO i=1,NELEMS
            write(16,'(1(E20.13,1X))') resElem(i, 11)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a59,a35)') '<DataArray type="Float32" Name="IO" NumberOfComponents="1" ',
     & 'ComponentName0="IO" format="ascii">'
      DO i=1,NELEMS
            write(16,'(1(E20.13,1X))') resElem(i, 12)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a59,a35)') '<DataArray type="Float32" Name="CK" NumberOfComponents="1" ',
     & 'ComponentName0="CK" format="ascii">'
      DO i=1,NELEMS
            write(16,'(1(E20.13,1X))') resElem(i, 11) - sqrt(2.d0) * resElem(i, 9) / 3.d0
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
C
      close(16)
C
      RETURN
      END

C-----------------------------------------------------------------------Searstr----------- 
C     Funcion que devuelve un valor logico que permite saber si la lectura es 
C     adecuadad.
C 
C     Busca Str en la unidad Lu y devuelve .true.// si lo encuentra. Si no
C     lo encuentra devuelve un .false.
C
C-----------------------------------------------------------------------------------------
      function Searstr (Lu, Str)
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
   10       Write (*,20) Lu,Str
      stop
C
   20 format ('###..Error en la funcion Searstr (Unidad=)',I3,2X,15A)
C
      end function Searstr
C---------------------------------------------------------------------------------------------