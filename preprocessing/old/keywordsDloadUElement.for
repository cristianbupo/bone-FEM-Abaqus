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

      subroutine UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS, NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3 NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4 PERIOD)
C
      include 'ABA_PARAM.INC'
	include 'conec.for'
C
      dimension RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
C
	real*8    x(dim,nnod),de(8)
      integer   k1,k2
C
      integer debugUnit
      character*276         filename
      character(256)        JOBDIR
      character(256)        JOBNAME
C
      debugUnit = 15  ! or any unused Fortran unit number
C
C     Se llaman las propiedades del modelo (parametros)
	call INPUT(props,de,NPROPS)
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
C
C     Funcion que regresa RHS y AMATRX
C    
      call ENSAMBLE(de,x,du,u,v,nst,ndofel,MDLOAD,NDLOAD,JDLTYP,ADLMAG,
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
C-------------------------------------------------------------INPUT---------------
C															
C     Funcion de entrada INPUT
C      
C     Obtiene los parametros del modelo
C     A definir
C									
C----------------------------------------------------------------------------------
C
      subroutine INPUT(propr,de,NPROPS)
C
	include 'conec.for'
C
      real*8 propr(*),de(8)
C
C     Almacenamos las propiedades enteras y reales en el vector de 
C
	do i=1,8
	  de(i) = propr(i)
        myprops(i) = propr(i)
	enddo
C
      return
      end
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
   10 	Write (*,20) Lu,Str
	stop
C
   20 format ('###..Error en la funcion Searstr (Unidad=)',I3,2X,15A)
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
      include 'conec.for'
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
C       Funciones de forma para descripcion de variable escalar
	  N(1,n1)    = shp(1,n1)
C
C       Funciones de forma para descripcion de variable vectorial	  
	  Ne(1,col)  = shp(1,n1)
	  Ne(2,col+1)= shp(1,n1)
C
C       Funciones de forma derivadas en el espacio para variable escalar	
	  B(1,n1)    = shp(2,n1)
	  B(2,n1)    = shp(3,n1)
C
C       Funciones de forma derivadas para el caso elastico
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
      include 'conec.for'
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
	1             dNeta(i)*d2eta_xx
	  d2shp(2,i)=2.d0*d2Nchi_eta(i)*detay*dchiy+dNchi(i)*d2chi_yy+
	1             dNeta(i)*d2eta_yy
        d2shp(3,i)=d2Nchi_eta(i)*(detay*dchix+detax*dchiy)+
	1             dNchi(i)*d2chi_xy+dNeta(i)*d2eta_xy
	enddo
C
c     Calculo de las derivadas de las funciones de forma con respecto a x
C
       shp(2,1) = - 0.25*((1.0-eta)*dchix + (1.0-chi)*detax)
       shp(2,2) =   0.25*((1.0-eta)*dchix - (1.0+chi)*detax)
       shp(2,3) =   0.25*((1.0+eta)*dchix + (1.0+chi)*detax)
       shp(2,4) = - 0.25*((1.0+eta)*dchix - (1.0-chi)*detax)

c     Calculo de las derivadas de las funciones de forma con respecto a y

       shp(3,1) = - 0.25*((1.0-eta)*dchiy + (1.0-chi)*detay)
       shp(3,2) =   0.25*((1.0-eta)*dchiy - (1.0+chi)*detay)
       shp(3,3) =   0.25*((1.0+eta)*dchiy + (1.0+chi)*detay)
       shp(3,4) = - 0.25*((1.0+eta)*dchiy - (1.0-chi)*detay)
C
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
      subroutine ENSAMBLE(de,x,du,u,v,nst,ndofel,MDLOAD,NDLOAD,JDLTYP,ADLMAG,
     1 nrhs,dtime,svars,nsvars,jelem,time,p,m_k)
C
      include 'conec.for'
C
C     Entradas
      integer nst,ndofel,nrhs,nsvars,jelem
      integer KDLOAD,MDLOAD,NDLOAD,JDLTYP(MDLOAD,*)
      real*8  ADLMAG(MDLOAD,*)
      real*8  x(dim,nnod),du(ndofel,*),u(ndofel),v(ndofel),de(8)
	real*8  dtime,time(2),svars(nsvars)
C
C     Salidas
      real*8  p(ndofel,nrhs),m_k(ndofel,ndofel)
C
C     Arreglos y variables de la subrutina
C
C     Generales
      integer k1,k2,k3,cc,ff,fil,col,tipo,colg,filg,colp,filp
      real*8  du1(nnod),du2(nnod)
	real*8  u1(nnod),u2(nnod)
C     real*8  du1(nnod),du2(nnod),du3(nnod),du4(nnod),du5(nnod)
C	real*8  u1(nnod),u2(nnod),u3(nnod),u4(nnod),u5(nnod)
	real*8  paux(ndofel),Kelast(dim*nnod,dim*nnod)
	real*8  masa_el(dim*nnod,dim*nnod),vector_carga(dim*nnod)
	real*8  Def(dim*nnod)
C
C     Variables de la carga distribuida
      integer n1, n2
      real*8 mag
      real*8  x1(dim),x2(dim)
      real*8  f, fx, fy
C      real*8  len, angulo, 
C     Inicializacion del tiempo
      if (dtime.eq.0.0) then
        dtime=1.e-15
      endif
C
C     Inicializacion de los vectores solucion
      do i=1,nnod
	  du1(i) = du(2*(i-1)+1,1)
	  du2(i) = du(2*(i-1)+2,1)
C	  du3(i) = du(5*(i-1)+3,1)
C	  du4(i) = du(5*(i-1)+4,1)
C	  du5(i) = du(5*(i-1)+5,1)
	  u1(i)  = u(2*(i-1)+1)
	  u2(i)  = u(2*(i-1)+2)
C	  u3(i)  = u(5*(i-1)+3)
C	  u4(i)  = u(5*(i-1)+4)
C	  u5(i)  = u(5*(i-1)+5)
	enddo
C
C     Inicializacion de matrices y variables
      m_k   =   0.d0
      p     =   0.d0
C
C     Ensamblar la matriz de rigidez elastica en la matriz tangente global
C     Llamado de la matriz de rigidez para la expansion
      call matriz_rigidez_el(u,ndofel,de,x,jelem,time(2),Kelast)
      do k1=1,nnod
        do k2=1,nnod
           fil=2*(k1-1)+1
           col=2*(k2-1)+1
           ff =dim*(k1-1)+1
           cc =dim*(k2-1)+1
	     do i=1,dim
	       do j=1,dim
	         colg = col+(j-1)
		   filg = fil+(i-1)
	         colp = cc +(j-1)
		   filp = ff +(i-1) 
               m_k(filg,colg) = m_k(filg,colg)+Kelast(filp,colp)
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

      DO KDLOAD = 1, NDLOAD
            n1 = JDLTYP(KDLOAD,1) ! Cara de la carga distribuida
            mag = ADLMAG(KDLOAD,1) ! Magnitud de la carga distribuida

            n2 = mod(n1, nnod)+1

            x1 = x(:,n1)
            x2 = x(:,n2)

C            len = sqrt((x2(1)-x1(1))**2 + (x2(2)-x1(2))**2)
C            angulo = atan2(x2(2)-x1(2), x2(1)-x1(1))
            f = mag/2 ! *len
            fx = -f*(x2(2)-x1(2)) !/len
            fy = f*(x2(1)-x1(1)) !/len

            p(2*(n1-1)+1, 1) = p(2*(n1-1)+1, 1) + fx
            p(2*(n1-1)+2, 1) = p(2*(n1-1)+2, 1) + fy
            p(2*(n2-1)+1, 1) = p(2*(n2-1)+1, 1) + fx
            p(2*(n2-1)+2, 1) = p(2*(n2-1)+2, 1) + fy

      enddo
C
C     Insercion del vector de deformacion
	call vector_deformacion(u,ndofel,de,x,jelem,time(2),Def)
      do k1=1,nnod
	  fil=2*(k1-1)+1 !CUIDADO
	  ff =dim*(k1-1)+1
	  do i=1,dim
	    filg = fil+(i-1)
	    filp = ff +(i-1) 
C	    p(fil,1)=p(fil,1)+Def(ff)
          p(filg,1)=p(filg,1)+Def(filp)
	  enddo
	enddo
C
      return
      end
C------------------------------------------------------------------------------------------
C---------------------------------------------------------------MATRIZ_DIFUSION_SCH________
C
C	 Funcion matriz_rigidez_el: matriz de rigidez elastica
c
C------------------------------------------------------------------------------------------
	subroutine matriz_rigidez_el(u,ndofel,de,x,jelem,t,Cmat)
C
	include 'conec.for'
C
C     Variables de entrada
      integer  ndofel,jelem
	real*8   u(ndofel),x(dim,nnod),t,de(8)
C
C     variables de salida
      real*8   Cmat(dim*nnod,dim*nnod)
C
C     Variables de la subrutina
C     Generales
      integer  i,j,k
	real*8   wg(2),sg(2),matriz(dim*nnod,dim*nnod)
C	real*8   u1(nnod),u2(nnod),u3(nnod),u4(nnod),u5(nnod)
	real*8   u1(nnod),u2(nnod)
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
	dx	     =   0.d0
	Nmat2D     =   0.d0
	Bmat2D     =   0.d0
	Nmatel2D   =   0.d0
	Bmatel2D   =   0.d0
	Bmatplano2D=   0.d0
	Bmataxi2D  =   0.d0
	Cmat 	     =   0.d0
	Dmat       =   0.d0
	matriz     =   0.d0
	r0         =   de(5)
C
C     Extraccion de los valores de u y v
      do i=1,nnod
	  u1(i)  = u(2*(i-1)+1)
	  u2(i)  = u(2*(i-1)+2)
C	  u3(i)  = u(5*(i-1)+3)
C	  u4(i)  = u(5*(i-1)+4)
C	  u5(i)  = u(5*(i-1)+5)
	enddo
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
	1           jelem,t,Dmat)
C
C	      Se calculan las funciones de forma
	      call f_forma2D(chi,eta,x,shp2D,xjac)
C
C           Se obtienen las matrices de calculo
            call Matrices_de_calculo2D(shp2D,Nmat2D,Bmat2D,Nmatel2D,
	1           Bmatel2D)
C
C           Se calcula el diferencial de la integral
            dx = wg(i)*wg(j)*xjac
C
C           Tipo de analisis
            if(axi.eq.0)then ! Caso plano: esfuerzo o deformacion plana
	         Bmatplano2D = Bmatel2D(1:3,:)
C              Se calculan las matriz elemental correspondiente a reaccion
	          matriz = matmul(transpose(Bmatplano2D),
	1                   matmul(Dmat(1:3,1:3),Bmatplano2D))
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
	1                   matmul(Dmat(1:4,1:4),Bmataxi2D))
            endif
C
C           Se llevan a cabo las operaciones de sumatoria
	      Cmat= Cmat + matriz * r * dx
	    enddo
	  enddo
	endif
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
	include 'conec.for'
C
C     Variables de entrada
      integer  ndofel,jelem
      real*8   chi,eta,zita
	real*8   u(ndofel),x(dim,nnod),t,de(8)
C
C     variables de salida
      real*8   Dmat(3*(dim-1)+axi,3*(dim-1)+axi)
C
C     Variables de la subrutina
C     Generales
      integer  i,j,k
C	real*8   u1(nnod),u2(nnod),u3(nnod),u4(nnod),u5(nnod)
	real*8   u1(nnod),u2(nnod)
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
      E      = myprops(3)
	nu     = myprops(4) 
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
          endif
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
	  endif
	endif

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
C---------------------------------------------------------------VECTOR DE DEFORMACION______
C
C	 Funcion VECTOR DE DEFORMACION INICIAL IMPUESTA
c
C------------------------------------------------------------------------------------------
	subroutine vector_deformacion(u,ndofel,de,x,jelem,t,Cmat)
C
	include 'conec.for'
C
C     Variables de entrada
      integer  ndofel,jelem
	real*8   u(ndofel),x(dim,nnod),t,de(8)
C
C     variables de salida
      real*8   Cmat(dim*nnod)
C
C     Variables de la subrutina
C     Generales
      integer  i,j,k
	real*8   wg(2),sg(2),vector(dim*nnod)
C	real*8   u1(nnod),u2(nnod),u3(nnod),u4(nnod),u5(nnod)
	real*8   u1(nnod),u2(nnod)
	real*8   chi,eta,dx,xjac
	real*8   Dmat(3*(dim-1)+axi,3*(dim-1)+axi)
C
C     En 2D
      real*8   shp2D(3,nnod),Nmat2D(1,nnod)
	real*8   Bmat2D(2,nnod),Nmatel2D(2,2*nnod),Bmatel2D(4,2*nnod)
	real*8   Bmatplano2D(3,2*nnod),r,r0,Bmataxi2D(4,2*nnod)
	real*8   Def2D(4)
C
C     Inicializacion de variables
	xjac	     =   0.d0
	chi	     =   0.d0
	eta        =   0.d0
	zita       =   0.d0
	dxchi      =   0.d0
	Nmat2D     =   0.d0
	Bmat2D     =   0.d0
	Nmatel2D   =   0.d0
	Bmatel2D   =   0.d0
	Bmatplano2D=   0.d0
	Bmataxi2D  =   0.d0
	Cmat 	     =   0.d0
	Dmat       =   0.d0
	vector     =   0.d0
	r0         =   de(5)
C
C     Extraccion de los valores de u y v
      do i=1,nnod
	  u1(i)  = u(2*(i-1)+1)
	  u2(i)  = u(2*(i-1)+2)
C	  u3(i)  = u(5*(i-1)+3)
C	  u4(i)  = u(5*(i-1)+4)
C	  u5(i)  = u(5*(i-1)+5)
	enddo
C
C	Calculo de los puntos de Gauss
	call bgauss(sg,wg)
C
      if(dim.eq.2)then
	  do i = 1,2
  	    chi = sg(i)
	    do j = 1,2
	      eta = sg(j)
C
C           Inicializacion de la matriz de constantes elasticas
            call matriz_constantes_el(chi,eta,zita,u,ndofel,de,x,
	1           jelem,t,Dmat)
C
C	      Se calculan las funciones de forma
	      call f_forma2D(chi,eta,x,shp2D,xjac)
C
C           Se obtienen las matrices de calculo
            call Matrices_de_calculo2D(shp2D,Nmat2D,Bmat2D,Nmatel2D,
	1           Bmatel2D)
C
C           Llamado al vector deformacion
            call deformacion(chi,eta,zita,u,ndofel,de,x,jelem,
	1      t,def2D)
C
C           Se calcula el diferencial de la integral
            dx = wg(i)*wg(j)*xjac
C
C           Tipo de analisis
            if(axi.eq.0)then ! Caso plano: esfuerzo o deformacion plana
	         Bmatplano2D = Bmatel2D(1:3,:)
C              Se calculan las matriz elemental correspondiente a reaccion
	          vector = matmul(transpose(Bmatplano2D),
	1                   matmul(Dmat(1:3,1:3),Def2D(1:3)))
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
	          vector = matmul(transpose(Bmataxi2D),
	1                   matmul(Dmat(1:4,1:4),Def2D))
            endif
C
C           Se llevan a cabo las operaciones de sumatoria
	      Cmat= Cmat + vector * r * dx
	    enddo
	  enddo
	endif
C
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
      include 'ABA_PARAM.INC'
      include    'conec.for'
C
      integer    i,J,k,jelem,j2,n,M
      real*8     ESFUERZOS(3),bmat(3,nnod*2),
     1           shp(3,nnod),xjac,x(2,nnod),
     2           DEFORMACIONES(3),DESP(nnod*2)
      real*8     dmat2(3*(dim-1),3*(dim-1))
      real*8     chi,eta,sg(2),puntos(2),a,b,aux_OI(2),
     1           puntos2(2),a8(4),b8(4),sigma_oct,thao_oct
      real*8     Nm(1,nnod),Ne(2,2*nnod),Bm(2,nnod),Be(4,2*nnod)
      real*8     I1,I2,I3,r,def_crec(4),dx
      real*8     zita,t,sigma_zz,Et(10),nut(10),Def2D(3),deforma
      real*8     Esf_Hid,Esf_VM
      real*8     OI, kOI

      real*8   def3D(6)
      real*8   dv,EsMax,EsMin,theta,thetadf
      real*8   E,nu
C
      data       puntos /-0.557350, 0.557350/
      data       puntos2/ 0.557350, -0.557350/
      data       a8/ 0.0, 1.0, 0.0, -1.0/
      data       b8/-1.0, 0.0, 1.0,  0.0/
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

      print *, ''
      print *, 'Elemento: ', jelem

        Dmat2 = 0.d0
        do J=1,nnod
          x(1,J)       = nodes(conectividades(jelem,J+1),1)
          x(2,J)       = nodes(conectividades(jelem,J+1),2)
          print *, 'Coordenadas: ', x(1,J), x(2,J)
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
           print *, 'Desplazamientos: ', DESP(2*(J-1)+1), DESP(2*(J-1)+2)
      enddo
C
      DEFORMACIONES = MATMUL(Be(1:3, :),DESP) !-def2D
C
      ESFUERZOS = 0.d0
      ESFUERZOS = MATMUL(Dmat2,DEFORMACIONES) !-def2D)   
C           Centro y radio del esfuerzo (circulo de Mohr)
      radio  = sqrt(0.25d0*(ESFUERZOS(1) - ESFUERZOS(2))**2)
      centro = 0.5d0*(ESFUERZOS(1) + ESFUERZOS(2))
C
      EsMax = centro+radio
      EsMin = centro-radio
C
      if(tipo_def.eq.2)then
            E      = myprops(3)
            nu     = myprops(4) 
            sigma_zz = (nu*E/((1.d0+nu)*(1.d0-2.d0*nu)))
     1      *(DEFORMACIONES(1)+DEFORMACIONES(2))
C           Esfuerzos equivalentes de Von Mises
      else if(tipo_def.eq.1)then
            sigma_zz = 0.d0
      endif
C
      Esf_VM = sqrt(((EsMax-EsMin)**2+(EsMin-sigma_zz)**2
     1      +(sigma_zz-EsMax)**2)/2.d0)
             
      Esf_Hid=-(1.d0/3.d0)*(EsMax+EsMin+sigma_zz)
C
      thao_oct=Esf_VM*sqrt(2.d0)/3.d0
C     
      kOI = 0.5d0
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


      print *, 'Desplazamientos: ', DESP
      print *, 'Deformaciones: ', DEFORMACIONES
      print *, 'Esfuerzos: ', ESFUERZOS
C
      enddo
C
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
	include 'conec.for'
C
C     Variables de entrada
      integer  ndofel,jelem
      real*8   chi,eta,zita
	real*8   u(ndofel),x(dim,nnod),t,de(8)
C
C     variables de salida
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
C        Def2D(2)=0.3d0
C	endif
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
      INCLUDE 'ABA_PARAM.INC'
      INCLUDE 'conec.for'
C
      DIMENSION ARRAY(513),JRRAY(NPRECD,513),TIME(2)
      EQUIVALENCE (ARRAY(1),JRRAY(1,1))
C
      LOGICAL :: firstKey = .TRUE., lastKey = .FALSE.
      INTEGER :: locID = 0
      INTEGER :: flagOutput = 0
      INTEGER :: prevFlagOutput = 0
      INTEGER :: i, j, k
      real*8 S11, S22, S12
      real*8 S1, S2, S_avg, R
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
                  j = j + 1
                  locID = JRRAY(1,6)
            END IF
C
C RECORD 101 CONTAINS VALUES FOR U
C           
            IF (KEY.EQ.101 .AND. locID .EQ. 1) THEN
                  k = k+1
                  resNod(k, 1) = ARRAY(4)
                  resNod(k, 2) = ARRAY(5)
            END IF
C
C RECORD 21 CONTAINS VALUES FOR E
C
            IF (KEY.EQ.21 .AND. locID .EQ. 1) THEN
                  resElem(j, 1) = ARRAY(3)
                  resElem(j, 2) = ARRAY(4)
                  resElem(j, 3) = ARRAY(5)
                  resElem(j, 4) = ARRAY(6)
            END IF
C
C RECORD 11 CONTAINS VALUES FOR S
C
            IF (KEY.EQ.11 .AND. locID .EQ. 1) THEN
                  resElem(j, 5) = ARRAY(3)
                  resElem(j, 6) = ARRAY(4)
                  resElem(j, 7) = ARRAY(5)
                  resElem(j, 8) = ARRAY(6)
            END IF
C
C RECORD 12 CONTAINS VALUES FOR SINV
C
            IF (KEY.EQ.12 .AND. locID .EQ. 1) THEN
                  resElem(j, 9) = ARRAY(3)
            END IF
C
C RECORD 201 CONTAINS ...
C
            IF (KEY.EQ.201) THEN
                  k = k+1
                  print *, 'I''m saving displacement info!!!!!!!!!!!!!!'
                  resNod(k, 1) = ARRAY(4)
                  resNod(k, 2) = ARRAY(5)
                  print *, '' 
                  print *, 'Displacement: ', ARRAY(4), ARRAY(5)
                  print *, 'Displacement: ', resNod(k, 1), resNod(k, 2)
            END IF
C
      END DO
C
 110  CONTINUE

C     Cálculo de los esfuerzos y las deformaciones
      call outsigma()
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
      write(16,'(a25)') '<PointData Vectors="''U''">'
      write(16,'(a58,a55)') '<DataArray type="Float64" Name="U" NumberOfComponents="2" ',
     & 'ComponentName0="U1" ComponentName1="U2" format="ascii">'
      DO i=1,NUMNODE
            write(16,'(2(E20.13,1X))') resNod(i, 1), resNod(i, 2)
      END DO
      write(16,'(a12)') '</DataArray>'
      write(16,'(a12)') '</PointData>'
C
      write(16,'(a64,a40)') '<DataArray type="Float64" Name="Load" NumberOfComponents="1" ',
     & 'ComponentName0="Load" format="ascii">'
      DO i=1,NELEMS
            write(16,'(1(E20.13,1X))') resElem(i, 9)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16, '(a46)') '<CellData Tensors="''E_Centroid'',''S_Centroid''">'
      write(16,'(a67,a99)') '<DataArray type="Float64" Name="E_Centroid" NumberOfComponents="4" ',
     & 'ComponentName0="E11" ComponentName1="E22" ComponentName2="E33" ComponentName3="E12" format="ascii">'
      DO i=1,NELEMS
            write(16,'(4(E20.13,1X))') resElem(i, 1), resElem(i, 2), resElem(i, 3), resElem(i, 4)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a67,a99)') '<DataArray type="Float64" Name="S_Centroid" NumberOfComponents="4" ',
     & 'ComponentName0="S11" ComponentName1="S22" ComponentName2="S33" ComponentName3="S12" format="ascii">'
      DO i=1,NELEMS
            write(16,'(4(E20.13,1X))') resElem(i, 5), resElem(i, 6), resElem(i, 7), resElem(i, 8)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a64,a40)') '<DataArray type="Float64" Name="S_Mises" NumberOfComponents="1" ',
     & 'ComponentName0="S_Mises" format="ascii">'
      DO i=1,NELEMS
            write(16,'(1(E20.13,1X))') resElem(i, 9)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a62,a38)') '<DataArray type="Float64" Name="S_Hyd" NumberOfComponents="1" ',
     & 'ComponentName0="S_Hyd" format="ascii">'
      DO i=1,NELEMS
            write(16,'(1(E20.13,1X))') resElem(i, 10)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a62,a38)') '<DataArray type="Float64" Name="S_Oct" NumberOfComponents="1" ',
     & 'ComponentName0="S_Oct" format="ascii">'
      DO i=1,NELEMS
            write(16,'(1(E20.13,1X))') resElem(i, 11)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a59,a35)') '<DataArray type="Float64" Name="OI" NumberOfComponents="1" ',
     & 'ComponentName0="OI" format="ascii">'
      DO i=1,NELEMS
            write(16,'(1(E20.13,1X))') resElem(i, 12)
      END DO
      write(16,'(a12)') '</DataArray>'
C
      write(16,'(a59,a35)') '<DataArray type="Float64" Name="CK" NumberOfComponents="1" ',
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
      close(16)

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
      include 'ABA_PARAM.INC'
      include 'conec.for'
C
      logical            :: Searstr
      character(256)        JOBDIR
	character*276         filename
	integer               i,j,k
C
      if (LOP.eq.0) then
C       Llamada al archivo de grupos físicos
            call GETOUTDIR(JOBDIR,LENJOBDIR)
            filename=' '
            filename(1:lenjobdir)=jobdir(1:lenjobdir)
            filename(lenjobdir+1:lenjobdir+19)='/gruposFisicos.txt'
C
            open(UNIT=14,file=filename(1:lenjobdir+19), status='old')
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
C
        filename(lenjobdir+1:lenjobdir+19)='/conectividades.inp'
C
        open(UNIT=15,file=filename(1:lenjobdir+19), status='old')
          if (Searstr (15,'*ELEMENT,TYPE=U1')) then
            READ(15,*)((conectividades(i,j),j=1,nnod + 1),i=1,NELEMS)
          else
            stop '###..Error en lectura'
          end if
        close(15)
C
C       Llamada al archivo de nodos
        call GETOUTDIR(JOBDIR,LENJOBDIR)
        filename=' '
	  filename(1:lenjobdir)=jobdir(1:lenjobdir)
        filename(lenjobdir+1:lenjobdir+10)='/nodos.inp'
C
	  open(UNIT=16,file=filename(1:lenjobdir+10), status='old')
	    if (Searstr (16,'*NODE,NSET=N2')) then
  	      READ(16,*) (k,(nodes(i,j),j=1,dim),i=1,NUMNODE)
	    else
	      stop '###..Error en lectura'
	    end if
	  close(16)
C       Se lee el archivo de entrada inp del analisis
	  call GETOUTDIR(JOBDIR,LENJOBDIR)	
        filename=' '
        filename(1:lenjobdir)=jobdir(1:lenjobdir)
        filename(lenjobdir+1:lenjobdir+14)='/contorno.inp'
C
        open(UNIT=15,file=filename(1:lenjobdir+14), status='old')
        if (Searstr (15,'*NSET,NSET=contorno1')) then
         READ(15,*)((contorno1(i,j),j=1,6),i=1,filascontorno1)
        else
        stop '###..Error en lectura'
        end if
        close(15)
C
C       Se lee el archivo de entrada inp del analisis
	  call GETOUTDIR(JOBDIR,LENJOBDIR)	
        filename=' '
        filename(1:lenjobdir)=jobdir(1:lenjobdir)
        filename(lenjobdir+1:lenjobdir+14)='/contorno.inp'
C
        open(UNIT=15,file=filename(1:lenjobdir+14), status='old')
        if (Searstr (15,'*NSET,NSET=contorno2')) then
         READ(15,*)((contorno2(i,j),j=1,6),i=1,filascontorno2)
        else
        stop '###..Error en lectura'
        end if
        close(15)
      endif
	return
	end