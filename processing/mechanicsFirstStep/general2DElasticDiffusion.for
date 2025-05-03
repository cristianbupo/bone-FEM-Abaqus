!-------------------------------------------------------------------------------------------
!   Elaborado por Diego Alexander Garzón Alvarado
!                 Oscar Rodrigo Lopez
!                 Profesores Asociados Ingenieria Mecánica y Mecatrónica
!                 Universidad Nacional de Colombia
!                 Sede Bogotá
! 
!   copy right - Universidad Nacional de Colombia
!
!   Software general para la solución de ecuaciones de reacción-convección-difusión
!   en 2D.
!
!   La aproximación se realiza mediante el método de Petrov-Galerkin en ABAQUS
!
!   SALIDA
!   La solución de este problema arroja un archivo *.dat que contiene los resultados
!
!   ENTRADA:
!   Requiere un archivo de conectividades.inp ==> archivo con las conectividades de los elementos
!                          nodos.inp          ==> archivo con los puntos nodales
!                          analisis.inp       ==> script de ABAQUS
!                          conec.txt          ==> variables globales
!-------------------------------------------------------------------------------------------
      include 'debug_utils.for'
      include 'reading_utils.for'
      include 'writing_utils.for'

      
      subroutine UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS, NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3 NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4 PERIOD)
!     
      use reading_utils
      include 'ABA_PARAM.INC'
!
      dimension RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
!
      real*8    x(dim,nnod)!,de(8)
      integer   k1,k2
!
      integer debugUnit
      character*276         filename
      character(256)        JOBDIR
      character(256)        JOBNAME
!
      logical change
!
      debugUnit = 15  ! or any unused Fortran unit number
!
!     Se llaman las propiedades del modelo (parametros)
!      call INPUT(props,de,NPROPS)
!
!     Inicializacion

      DO 6 K1=1,NDOFEL                      
         RHS(K1,NRHS)=0.0
      DO 4 K2=1,NDOFEL
         AMATRX(K2,K1)=0.0
    4 CONTINUE                                      
    6 CONTINUE   
!
      do k1=1,nnod
         do j=1,dim
            x(j,k1)=coords(j,k1)
         enddo
      enddo
!
!     Funcion que regresa RHS y AMATRX
!    
      call ENSAMBLE(de,x,du,u,v,nst,ndofel,KSTEP,KINC,nrhs,dtime,svars,nsvars,jelem,time,RHS,AMATRX)
!     
      return
      end
!---------------------------------------------------------------------------------
!------------------------------------------------------------------------- BGAUSS2 -------
!
!     Funcion bgauss(sg,wg)
!
!----------------------------------------------------------------------------------------
      SUBROUTINE bgauss(sg,wg)
      real*8, intent(out) :: sg(2),wg(2)

      sg(1) = -0.577350269189626
      sg(2) =  0.577350269189626

      wg(1) =  1.0
      wg(2) =  1.0

      return
      end
!---------------------------------------------------------------------------------------------
!----------------------------------------------------------------------Matrices_de_calculo2D--
!    
!     Devuelve las matrices de calculo requeridas en los analisis escalares
!
!---------------------------------------------------------------------------------------------
      subroutine Matrices_de_calculo2D(shp,N,B,Ne,Be)
!     
      use reading_utils
!
      integer  col,n1
      real*8, intent(in) ::  shp(3,nnod)
      real*8, intent(out) :: N(1,nnod),Ne(2,2*nnod), B(2,nnod),Be(4,2*nnod)
!
!     Puesta a cero
      N = 0.d0
      Ne= 0.d0
      B = 0.d0
      Be= 0.d0
!
!     Incializacion
      do n1 = 1,nnod
         col=2*(n1-1)+1
!
!        Funciones de forma para descripcion de variable escalar
         N(1,n1)    = shp(1,n1)
!
!        Funciones de forma para descripcion de variable vectorial     
         Ne(1,col)  = shp(1,n1)
         Ne(2,col+1)= shp(1,n1)
!
!        Funciones de forma derivadas en el espacio para variable escalar   
         B(1,n1)    = shp(2,n1)
         B(2,n1)    = shp(3,n1)
!
!        Funciones de forma derivadas para el caso elastico
         Be(1,col)  = shp(2,n1)
         Be(2,col+1)= shp(3,n1)
         Be(3,col)  = shp(3,n1)
         Be(3,col+1)= shp(2,n1)
         Be(4,col)  = shp(1,n1)
      enddo  
      !
      return
      end
!---------------------------------------------------------------------------------------------
!----------------------------------------------------------------------Matrices_de_calculo1D--
!
!     Devuelve las matrices de calculo requeridas en los analisis escalares en 1D
!
!---------------------------------------------------------------------------------------------
      subroutine Matrices_de_calculo1D(shp, N, Ne, B, Be)
!
      use reading_utils
!
      integer n1, col
      real*8, intent(in) :: shp(2, nnod)
      real*8, intent(out) :: N(1, nnod), Ne(1, nnod), B(1, nnod), Be(2, nnod)
!
!     Puesta a cero
      N  = 0.d0
      Ne = 0.d0
      B  = 0.d0
      Be = 0.d0
!
!     Inicializacion
      do n1 = 1, nnod
         col = n1
!
!        Funciones de forma para descripcion de variable escalar
         N(1, n1) = shp(1, n1)
!
!        Funciones de forma para descripcion de variable vectorial (simple en 1D)
         Ne(1, col) = shp(1, n1)
!
!        Funciones de forma derivadas en el espacio para variable escalar
         B(1, n1) = shp(2, n1)
!
!        Funciones de forma derivadas para el caso elastico en 1D
         Be(1, col) = shp(2, n1)  ! Derivada respecto a x
         Be(2, col) = shp(1, n1)  ! Valor de la funcion de forma
      enddo  
!
      return
      end
!----------------------------------------------------------------------------------------
!------------------------------------------------------------------------ F_FORMA2D------
!     Funcion f_forma2D(chi,eta,x,shp,xjac,d2shp)
!
!----------------------------------------------------------------------------------------
      subroutine f_forma2D(chi,eta,x,shp,xjac)
!   
      use reading_utils
!
!
      integer i
      real*8 chi,eta,xjac,shp(3,nnod),x(2,nnod),d2shp(3,nnod)
      real*8 dxchi,dxeta,dychi,dyeta,dchix,dchiy,detax,detay
      real*8 dNchi(nnod),dNeta(nnod)
      real*8 d2Nchi_eta(nnod)
      real*8 d2xchi_eta,d2ychi_eta,dxjacchi,dxjaceta,dxjacx,dxjacy
      real*8 d2eta_xx,d2chi_xx,d2chi_yy,d2eta_yy,d2eta_xy,d2chi_xy
!
      dxchi=0.d0
      dxeta=0.d0
      dychi=0.d0
      dyeta=0.d0
      dchix=0.d0
      dchiy=0.d0
      detax=0.d0
      detay=0.d0
!
      shp=0.d0
!
!     Primera fila de la matriz shp - Funciones de forma en los 8 nodos
      shp(1,1) = 0.25*(1.0 - chi)*(1.0 - eta)
      shp(1,2) = 0.25*(1.0 + chi)*(1.0 - eta)
      shp(1,3) = 0.25*(1.0 + chi)*(1.0 + eta)
      shp(1,4) = 0.25*(1.0 - chi)*(1.0 + eta)
!
!     Primeras derivadas de las funciones de forma con respecto a chi
      dNchi(1) = -0.25*(1.0 - eta)
      dNchi(2) =  0.25*(1.0 - eta)
      dNchi(3) =  0.25*(1.0 + eta)
      dNchi(4) = -0.25*(1.0 + eta)
!
!     Primeras derivadas de las funciones de forma con respecto a eta
      dNeta(1) = -0.25*(1.0 - chi)
      dNeta(2) = -0.25*(1.0 + chi)
      dNeta(3) =  0.25*(1.0 + chi)
      dNeta(4) =  0.25*(1.0 - chi)
!

!     Segundas derivadas de las funciones de forma con respecto a chi y eta
      d2Nchi_eta(1) =  0.25
      d2Nchi_eta(2) = -0.25
      d2Nchi_eta(3) =  0.25
      d2Nchi_eta(4) = -0.25
!
!     Calculo de la matriz jacobiana
      dxchi=0.25*((x(1,2)-x(1,1))*(1.0-eta)+(x(1,3)-x(1,4))*(1.0+eta))
      dxeta=0.25*((x(1,4)-x(1,1))*(1.0-chi)+(x(1,3)-x(1,2))*(1.0+chi))
      dychi=0.25*((x(2,2)-x(2,1))*(1.0-eta)+(x(2,3)-x(2,4))*(1.0+eta))
      dyeta=0.25*((x(2,4)-x(2,1))*(1.0-chi)+(x(2,3)-x(2,2))*(1.0+chi))
!
!     Calculo de las segundas derivadas de x y y
      d2xchi_eta=0.25*(-(x(1,2)-x(1,1))+(x(1,3)-x(1,4)))
      d2ychi_eta=0.25*(-(x(2,2)-x(2,1))+(x(2,3)-x(2,4)))
!
!     Calculo del determinante de la matriz jacobiana -Jacobiano -
!
      xjac = dxchi*dyeta - dxeta*dychi
!
!     Calculo de la matriz inversa de la matriz jacobiana
      dchix =   dyeta/xjac
      dchiy = - dxeta/xjac
      detax = - dychi/xjac
      detay =   dxchi/xjac
!
!     Calculo de las derivadas del jacobiano con respecto a chi
      dxjacchi = dxchi*d2ychi_eta - d2xchi_eta*dychi
      dxjaceta = d2xchi_eta*dyeta - dxeta*d2ychi_eta
!
!     Calculo de las derivadas del jacobiano con respecto a x
      dxjacx = dxjacchi*dchix + dxjaceta*detax
!     Calculo de las derivadas del jacobiano respecto a y
      dxjacy = dxjacchi*dchiy + dxjaceta*detay
!
!     Calculo de la segunda derivada de eta con respecto a x
      d2eta_xx=((-d2ychi_eta*detax*xjac)+(dxjacx*dychi))/(xjac**2)
!
!     Calculo de la segunda derivada de chi con respecto a x
      d2chi_xx=((d2ychi_eta*dchix*xjac)-(dxjacx*dyeta))/(xjac**2)
! 
!     Calculo de la segunda derivada de chi con respecto a y
      d2chi_yy=(-(d2xchi_eta*dchiy*xjac)+(dxjacy*dxeta))/(xjac**2)
!
!     Calculo de la segunda derivada de eta con respecto a y
      d2eta_yy=((d2xchi_eta*detay*xjac)-(dxjacy*dxchi))/(xjac**2)
!
!     Calculo de la segunda derivada de eta con respecto a x y luego y
      d2eta_xy=((-d2ychi_eta*detay*xjac)+(dxjacy*dychi))/(xjac**2)
!
!     Calculo de la segunda derivada de chi con respecto a x y luego y
      d2chi_xy=((d2ychi_eta*dchiy*xjac)-(dxjacy*dyeta))/(xjac**2)
!
!     Calculo de las segundas derivadas, la primera fila indica las segundas derivadas de
!     las funciones de forma con respecto a x, la segunda con respecto a y, y la tercera, los terminos cruzados
      do i =1,nnod     
      d2shp(1,i)=2.d0*d2Nchi_eta(i)*detax*dchix+dNchi(i)*d2chi_xx+dNeta(i)*d2eta_xx
      d2shp(2,i)=2.d0*d2Nchi_eta(i)*detay*dchiy+dNchi(i)*d2chi_yy+dNeta(i)*d2eta_yy
      d2shp(3,i)=d2Nchi_eta(i)*(detay*dchix+detax*dchiy)+dNchi(i)*d2chi_xy+dNeta(i)*d2eta_xy
      enddo
!
!     Calculo de las derivadas de las funciones de forma con respecto a x
!
      shp(2,1) = - 0.25*((1.0-eta)*dchix + (1.0-chi)*detax)
      shp(2,2) =   0.25*((1.0-eta)*dchix - (1.0+chi)*detax)
      shp(2,3) =   0.25*((1.0+eta)*dchix + (1.0+chi)*detax)
      shp(2,4) = - 0.25*((1.0+eta)*dchix - (1.0-chi)*detax)

!     Calculo de las derivadas de las funciones de forma con respecto a y
!
      shp(3,1) = - 0.25*((1.0-eta)*dchiy + (1.0-chi)*detay)
      shp(3,2) =   0.25*((1.0-eta)*dchiy - (1.0+chi)*detay)
      shp(3,3) =   0.25*((1.0+eta)*dchiy + (1.0+chi)*detay)
      shp(3,4) = - 0.25*((1.0+eta)*dchiy - (1.0-chi)*detay)
!
      return
      end
!----------------------------------------------------------------------------------------
!------------------------------------------------------------------------ F_FORMA1D------
!     Funcion f_forma1D(chi,x,shp,xjac)
!----------------------------------------------------------------------------------------
      subroutine f_forma1D(chi,x,shp,xjac)
         !
               use reading_utils
         !
               integer i
               real*8, intent(in) :: chi, x(1, nnod)
               real*8, intent(out) :: xjac, shp(2, nnod)
               real*8 dNchi(nnod)
               real*8 dxchi, dchix, detax
         !
               dxchi = 0.d0
         !
               shp = 0.d0
         !
         !     Funciones de forma lineales en el elemento 1D (dos nodos)
               shp(1,1) = 0.5*(1.0 - chi)
               shp(1,2) = 0.5*(1.0 + chi)
         !
         !     Derivadas de las funciones de forma respecto a chi
               dNchi(1) = -0.5
               dNchi(2) =  0.5
         !
         !     Cálculo del Jacobiano (dx/dchi)
               do i = 1, nnod
                  dxchi = dxchi + dNchi(i)*x(1,i)
               enddo
         !
         !     Determinante del Jacobiano
               xjac = dxchi
         !
         !     Inversa del Jacobiano
               dchix = 1.d0 / xjac
         !
         !     Derivadas de las funciones de forma respecto a x
               do i = 1, nnod
                  shp(2,i) = dNchi(i) * dchix
               enddo
         !
               return
               end
!------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------oute-----------
!       
!       Funcion calculateCMIthreshold()
!
!
!       Funcion para acumular los resultados de los elementos
!       
!
! -----------------------------------------------------------------------------------------
!
!      Argumentos: 
!
      subroutine calculateCMIthreshold()
!
      use reading_utils
      include 'ABA_PARAM.INC'
!
      integer  i,j,k,jelem,nCart
      real*8   x,minx,CMIcurrent

      CMIavg = 0.0
      CMIstd = 0.0
      CMImax = -huge(1.0d0)
      nCart = 0

      do jelem=1,NELEMS
         if (grupoFisico(jelem,2) .gt. 0) then
            CMIcurrent = CMICriteria(jelem)
            CMIavg = CMIavg + CMIcurrent
            nCart = nCart + 1

            minx = huge(1.0d0)
            do j = 1, nnod
               k = conectividades(jelem, j+1)
               x = abs(nodes(k,1))
               if (x < minx) then
                  minx = x
               endif
               
            enddo

            if (minx < 1.0d-5) then
               if (CMIcurrent>CMImax) then
                  CMImax = CMIcurrent
               endif
            endif
         endif
      enddo

      CMIavg = CMIavg / DBLE(nCart)

      do jelem=1,NELEMS
         if (grupoFisico(i,2) .gt. 0) then
            CMIstd = CMIstd + (CMICriteria(jelem) - CMIavg)**2
         endif
      enddo

      CMIstd = sqrt(CMIstd / DBLE(nCart)) ! Standard deviation of CMI

      CMIThreshold = CMIavg + stdWeight * CMIstd

      print *, ''
      print *, 'CMIavg:', CMIavg
      print *, 'CMIstd:', CMIstd
      print *, 'CMImax:', CMImax
      print *, 'CMIthreshold:', CMIThreshold

      return
      end
!---------------------------------------------------------------------------------------------
!-------------------------------------------------------------updateProps-----------
!
!      Funcion updateProps
!   
!      Funcion que actualiza los grupos físicos.
!
!-------------------------------------------------------------------------------*/
!-------------------------------------------------------------------------------*/
      subroutine updateProps()
!
      use reading_utils
      logical :: update(nelems)
      integer :: i, grupoFisicoActual
      integer :: newGroup(nelems)
!
!     Actualización de propiedades
!     
      update = .false.
      newGroup = 0

      print *, 'Updating properties...'

      ! Decidir si se va a cambiar la propiedad
      do i=1,NELEMS
         grupoFisicoActual = grupoFisico(i,2)

         if (grupoFisicoActual.gt.1) then
            if (CMICriteria(i) >= CMIThreshold) then
               update(i) = .true.
               newGroup(i) = 1
            end if
         else if (grupoFisicoActual==1) then
            update(i) = .true.
            newGroup(i) = 0
         endif
      enddo

      ! Cambiar la propiedad

      do i=1,NELEMS
         if (update(i)) then
            grupoFisico(i,2) = newGroup(i)
         end if
      enddo

      return
      end subroutine updateProps
!---------------------------------------------------------------------------------------------
!-------------------------------------------------------------ENSAMBLE-----------
!
!      Funcion ENSAMBLE
!   
!      Funcion que regresa la matriz de rigidez tangente AMATRX y
!      el residuo RHS.
!
!-------------------------------------------------------------------------------*/
!-------------------------------------------------------------------------------*/
      subroutine ENSAMBLE(de,x,du,u,v,nst,ndofel,KSTEP,KINC,
     1 nrhs,dtime,svars,nsvars,jelem,time,p,m_k)
!
      
      use reading_utils
!
!
!     Entradas
      integer nst,ndofel,nrhs,nsvars,jelem,jnod
      real*8  x(dim,nnod),du(ndofel,*),u(ndofel),v(ndofel)!,de(8)
      real*8  dtime,time(2),svars(nsvars)
!
!     Salidas
      real*8  p(ndofel,nrhs),m_k(ndofel,ndofel)
!
!     Arreglos y variables de la subrutina
!
!     Generales
      integer k1,k2,k3,cc,ff,fil,col,tipo,colg,filg,colp,filp
      integer ndofnod
      real*8  paux(ndofel),Kelast(dim*nnod,dim*nnod)
      real*8  Kdiff(nnod,nnod)
      real*8  Masa(nnod,nnod)
      real*8  Mec(dim*nnod)
      real*8  Reac(1,nnod)
      real*8 , allocatable :: myU(:,:), myDU(:,:)
!
!     Variables de la carga distribuida
      integer n1, n2, i, j, k
      real*8 mag
      real*8  x1(dim),x2(dim)
      real*8  force, fx, fy
!
      integer  grupo
      real*8  E, nu, D(2), f(2), g(2), h(2)

      logical :: found
!     Inicializacion de matrices y variables
      m_k   =   0.d0
      p     =   0.d0
      ndofnod = ndofel/nnod 
      paux = 0.d0
      Kelast = 0.d0
      Kdiff = 0.d0
      Masa = 0.d0
      Mec = 0.d0
      Reac = 0.d0
!     
      allocate(myU(ndofnod,nnod), myDU(ndofnod,nnod))

!     Inicializacion del tiempo
      if (dtime.eq.0.0) then
        dtime=1.e-15
      endif
!
!    Inicializacion de los residuos
      do i=1,nnod
         do j=1,ndofnod
            myU(j,i) = u(ndofnod*(i-1)+j)
            myDU(j,i) = du(ndofnod*(i-1)+j,1)
         enddo
	   enddo
!
!     Ensamblar la matriz de rigidez elastica en la matriz tangente global
!     Llamado de la matriz de rigidez para la expansion
      
      call matriz_rigidez_dif(u,ndofel,de,x,jelem,time(2),Kdiff)
      call matriz_masa(u,ndofel,de,x,jelem,time(2),Masa)
      call vector_reaccion(u,ndofel,de,x,jelem,time(2),Reac)
      call matriz_rigidez_el(u,ndofel,de,x,jelem,time(2),Kelast)
      call carga_distribuida(i, x, KINC, jelem, Mec)

      grupo = grupoFisico(jelem,2)+1
      E  = propiedades(grupo,1)
      nu = propiedades(grupo,2)
      D(1)  = propiedades(grupo,3)
      D(2)  = propiedades(grupo,4)
      f(1) = propiedades(grupo,5)
      f(2) = propiedades(grupo,6)
      g(1) = propiedades(grupo,7)
      g(2) = propiedades(grupo,8)
      h(1) = propiedades(grupo,9)
      h(2) = propiedades(grupo,10)
!
!      print '(F6.1, 1X, F6.1, 1X, F6.1, 1X, F6.1,'//
!     & '1X, F6.1, 1X, F6.1, 1X, F6.1, 1X, F6.1, 1X, F6.1)'
!     & , E, nu, D, f(1), f(2), g(1), g(2), h(1), h(2)
!
!      print *, (Reac(1,j), j=1,nnod)
      do k1=1,nnod
         ff = dim*(k1-1)+1
         do k2=1,nnod
            cc = dim*(k2-1)+1
            fil = ndofnod*(k1-1)+1
            col = ndofnod*(k2-1)+1
            
            ! Ensamble de la matriz de rigidez elastica
            do i=1,dim
               filg = fil+(i-1)
               filp = ff +(i-1)
               ! Ensamble del vector de carga
               p(filg,1) = p(filg,1) + Mec(filp)
               do j=1,dim
                  colg = col+(j-1)
                  colp = cc +(j-1)
                  ! Ensamble de la matriz de rigidez elastica
                  m_k(filg,colg) = m_k(filg,colg) + Kelast(filp,colp)
               enddo
            enddo

            ! Actualizacion de los indices 
            fil = fil + dim
            col = col + dim

            do i=1,2
               filg = fil + (i-1)
               colg = col + (i-1)
               ! Ensamble del vector de reaccion
               p(filg,1) = p(filg,1) - h(i)*Reac(1,k1) ! revisar
               ! Ensamble de la matriz de rigidez difusiva
               m_k(filg,colg) = m_k(filg,colg) - D(i)*Kdiff(k1,k2)
               m_k(filg,colg) = m_k(filg,colg) + g(i)*Masa(k1,k2)
            enddo

            ! Ensamble de las matrices de masa
            m_k(fil,col+1) = m_k(fil,col+1) + f(1)*Masa(k1,k2) 
            m_k(fil+1,col) = m_k(fil+1,col) + f(2)*Masa(k1,k2)

            ! Actualizacion de los indices 
            filg = fil + 2
            colg = col + 2
            m_k(filg,colg) = m_k(filg,colg) - Kdiff(k1,k2)
         enddo
      enddo
!
!     RHS:
      paux=matmul(m_k,u)
!
!     Ensamble del vector residuo
      do k1=1,ndofel
         p(k1,1) = p(k1,1) - paux(k1)
      enddo
!      
      return
      end
!------------------------------------------------------------------------------------------
!-------------------------------------------------------------carga_distribuida-----------
!
!      Funcion carga_distribuida
!
!
!      Funcion que regresa la carga distribuida en el elemento
!      Se utiliza para el analisis elastico
!
!-----------------------------------------------------------------------------------------
      subroutine carga_distribuida(i, x, KINC, jelem, vec)
      
      use reading_utils

      logical :: found
      integer :: n1, n2, i, KINC, jelem
      real*8 :: x1(2), x2(2), mag, vec(dim*nnod), x(dim,nnod)
      real*8 :: force, fx, fy

      vec = 0.0d0
!     TODO: Use a similar approach to KDLOAD, NDLOAD and ADLMAG to load elements
!     DO KDLOAD = 1, NDLOAD
!         n1 = JDLTYP(KDLOAD,1) ! Cara de la carga distribuida
!         mag = ADLMAG(KDLOAD,1) ! Magnitud de la carga distribuida
!     ENDDO

!     Insercion de carga distribuida

      if (KINC .le. nLoads) then 
         found = .false.

         do i = 1, listNElementLoads(KINC,1)
            if (jelem == elementFaces(i, 1)) then
               found = .true.
               exit
            end if
         enddo


         if (found) then ! Works just for 1 load per element
            n1 = elementFaces(i, 2) ! Cara de la carga distribuida
            mag = elementLoads(i, 1) ! Magnitud de la carga distribuida
            n2 = mod(n1, nnod)+1

            x1 = x(:,n1)
            x2 = x(:,n2)

            force = mag/2 ! *len
            fx = -force*(x2(2)-x1(2)) !/len
            fy = force*(x2(1)-x1(1)) !/len

            vec(2*n1-1) = fx
            vec(2*n1) = fy
            vec(2*n2-1) = fx
            vec(2*n2) = fy
         end if
      end if
      
      end subroutine carga_distribuida

!------------------------------------------------------------------------------------------
!----------------------------------------------------------------matriz_rigidez_el---------
!
!    Funcion matriz_rigidez_el: matriz de rigidez elastica
!
!------------------------------------------------------------------------------------------
      subroutine matriz_rigidez_el(u,ndofel,de,x,jelem,t,Cmat)
!
      use reading_utils
!
!     Variables de entrada
      integer  ndofel,jelem
      real*8   u(ndofel),x(dim,nnod),t!,de(8)
!
!     variables de salida
      real*8   Cmat(dim*nnod,dim*nnod) !dim seria en este caso el número de grados de libertad por nodo
!
!     Variables de la subrutina
!     Generales
      integer  i,j,k
      real*8   wg(2),sg(2),matriz(dim*nnod,dim*nnod)
      real*8   chi,eta,zita,dx,xjac
      real*8   Dmat(dim*(dim+1)/2+axi,dim*(dim+1)/2+axi)
!
!     En 2D
      real*8   shp2D(3,nnod),Nmat2D(1,nnod)
      real*8   Bmat2D(2,nnod),Nmatel2D(2,2*nnod),Bmatel2D(4,dim*nnod)
      real*8   Bmatplano2D(3,2*nnod),r,r0,Bmataxi2D(4,2*nnod)
!     En 1D
      real*8   shp1D(2,nnod), Nmat1D(1,nnod)
      real*8   Bmat1D(1,nnod), Nmatel1D(1,nnod), Bmatel1D(2,nnod)
!
!     Inicializacion de variables
      xjac       =   0.d0
      chi        =   0.d0
      eta        =   0.d0
      zita       =   0.d0
      dx         =   0.d0
      shp2D      =   0.d0
      Nmat2D     =   0.d0
      Bmat2D     =   0.d0
      Nmatel2D   =   0.d0
      Bmatel2D   =   0.d0
      Bmatplano2D=   0.d0
      Bmataxi2D  =   0.d0
      shp1D      =   0.d0
      Nmat1D     =   0.d0
      Bmat1D     =   0.d0
      Nmatel1D   =   0.d0
      Bmatel1D   =   0.d0
      Cmat       =   0.d0
      Dmat       =   0.d0
      matriz     =   0.d0
      r0         =   0.d0
!
!   Calculo de los puntos de Gauss
      call bgauss(sg,wg)
!
!   Bucle para cada punto de Gauss
      if(dim.eq.2)then
         do i = 1,2
            chi = sg(i)
         do j = 1,2
            eta = sg(j)
!
!           Inicializacion de la matriz de constantes elasticas
            call matriz_constantes_el(u,ndofel,de,x,jelem,t,Dmat)
!
!           Se calculan las funciones de forma
            call f_forma2D(chi,eta,x,shp2D,xjac)
!
!           Se obtienen las matrices de calculo
            call Matrices_de_calculo2D(shp2D,Nmat2D,Bmat2D,Nmatel2D,Bmatel2D)
!
!           Se calcula el diferencial de la integral
            dx = wg(i)*wg(j)*xjac
!
!           Tipo de analisis
            if(axi.eq.0)then ! Caso plano: esfuerzo o deformacion plana
               Bmatplano2D = Bmatel2D(1:3,:)
!              Se calculan las matriz elemental correspondiente a reaccion
               matriz = matmul(transpose(Bmatplano2D),matmul(Dmat(1:3,1:3),Bmatplano2D))
            r = 1
            elseif(axi.eq.1)then
!              Radio de giro
               r = dot_product(Nmat2D(1,:),x(1,:))+r0
!
!              Definicion de la matriz operador de deformaciones
               Bmataxi2D(1:3,:) = Bmatel2D(1:3,:)
               Bmataxi2D(4,:) = Bmatel2D(4,:)/r
!
!              Se calculan las matriz elemental correspondiente a reaccion
             matriz = matmul(transpose(Bmataxi2D),matmul(Dmat(1:4,1:4),Bmataxi2D))
            end if
!
!           Se llevan a cabo las operaciones de sumatoria
         Cmat= Cmat + matriz * r * dx
         enddo
         enddo
      end if
!
      if(dim.eq.1)then
         do i = 1,2
            chi = sg(i)
!
!           Inicializacion de la matriz de constantes elásticas
            call matriz_constantes_el(u,ndofel,de,x,jelem,t,Dmat)
!
!           Se calculan las funciones de forma
            call f_forma1D(chi,x,shp1D,xjac)
!
!           Se obtienen las matrices de cálculo
            call Matrices_de_calculo1D(shp1D,Nmat1D,Bmat1D,Nmatel1D,Bmatel1D)
!
!           Se calcula el diferencial de la integral
            dx = wg(i) * xjac
!
!           Se calculan la matriz elemental de rigidez
!           matriz = B^T * D * B
            do j = 1, nnod
               do k = 1, nnod
                  matriz(j,k) = Dmat(1,1) * Bmatel1D(1,j) * Bmatel1D(1,k)
               enddo
            enddo
!
!           Se llevan a cabo las operaciones de sumatoria
            Cmat = Cmat + matriz * dx
         enddo
      end if
      return 
      end
!------------------------------------------------------------------------------------------
!----------------------------------------------------------------matriz_rigidez_dif---------
!
!    Funcion matriz_rigidez_dif: matriz de rigidez difusiva
!
!   CMat : Matriz de rigidez difusiva
!------------------------------------------------------------------------------------------
      subroutine matriz_rigidez_dif(u,ndofel,de,x,jelem,t,Cmat)
!
      use reading_utils
!
!     Variables de entrada
      integer  ndofel,jelem
      real*8   x(dim,nnod),t!,de(8)
!
!     variables de salida
      real*8   Cmat(nnod,nnod)
!     Variables de la subrutina
!     Generales
      integer  i,j,k,l,m,n
      real*8   wg(2),sg(2),matriz(nnod,nnod)
      real*8   chi,eta,dx,xjac
      real*8   Dmat(dim*(dim+1)/2+axi,dim*(dim+1)/2+axi)
!
!     En 2D
      real*8   shp2D(3,nnod),Nmat2D(1,nnod)
      real*8   Bmat2D(2,nnod),Nmatel2D(2,2*nnod),Bmatel2D(4,2*nnod)
      real*8   Bmatplano2D(3,2*nnod),r,r0,Bmataxi2D(4,2*nnod)
!     En 1D
      real*8   shp1D(2,nnod), Nmat1D(1,nnod)
      real*8   Bmat1D(1,nnod), Nmatel1D(1,nnod), Bmatel1D(2,nnod)
!
!     Inicializacion de variables
      xjac       =   0.d0
      chi        =   0.d0
      eta        =   0.d0
      zita       =   0.d0
      dx         =   0.d0
      shp2D      =   0.d0
      Nmat2D     =   0.d0
      Bmat2D     =   0.d0
      Nmatel2D   =   0.d0
      Bmatel2D   =   0.d0
      Bmatplano2D=   0.d0
      Bmataxi2D  =   0.d0
      shp1D      =   0.d0
      Nmat1D     =   0.d0
      Bmat1D     =   0.d0
      Nmatel1D   =   0.d0
      Bmatel1D   =   0.d0
      Cmat       =   0.d0
      Dmat       =   0.d0
      matriz     =   0.d0
      r0         =   0.d0
!
!   Calculo de los puntos de Gauss
      call bgauss(sg,wg)
!
!   Bucle para cada punto de Gauss
      if(dim.eq.2)then
         do i = 1,2
            chi = sg(i)
            do j = 1,2
               eta = sg(j)
!
!            Se calculan las funciones de forma
               call f_forma2D(chi,eta,x,shp2D,xjac)
!
!           Se obtienen las matrices de calculo
               call Matrices_de_calculo2D(shp2D,Nmat2D,Bmat2D,Nmatel2D,Bmatel2D)
!
!           Se calcula el diferencial de la integral
               dx = wg(i)*wg(j)*xjac
!
!           Se calcula la matriz de difusion para un solo gdl
               matriz = matmul(transpose(Bmat2D),Bmat2D)
!
!           Se llevan a cabo las operaciones de sumatoria
               Cmat = Cmat + matriz * dx
            enddo
         enddo
!   Bucle para cada punto de Gauss en 1D
      else if (dim.eq.1) then
         do i = 1, 2
            chi = sg(i)
!
!           Se calculan las funciones de forma en 1D
            call f_forma1D(chi, x, shp1D, xjac)
!
!           Se obtienen las matrices de cálculo en 1D
            call Matrices_de_calculo1D(shp1D, Nmat1D, Nmatel1D, Bmat1D, Bmatel1D)
!
!           Se calcula el diferencial de la integral
            dx = wg(i) * xjac
!
!           Se calcula la matriz de difusion para un solo gdl
            matriz = matmul(transpose(Bmat1D), Bmat1D)
!
!           Se llevan a cabo las operaciones de sumatoria
            Cmat = Cmat + matriz * dx
         enddo
      end if
!
      return
      end
!------------------------------------------------------------------------------------------
!----------------------------------------------------------------matriz_masa---------
!
!    Funcion matriz_masa: matriz de masa de la absorción cruzada entre reactivos
!
!-------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
      subroutine matriz_masa(u,ndofel,de,x,jelem,t,Cmat)
!
      use reading_utils
!
!     Variables de entrada
      integer  ndofel,jelem
      real*8   x(dim,nnod),t!,de(8)
!
!     variables de salida
      real*8   Cmat(nnod,nnod)
!     Variables de la subrutina
!     Generales
      integer  i,j,k,l,m,n
      real*8   wg(2),sg(2),matriz(nnod,nnod)
      real*8   chi,eta,dx,xjac
      real*8   Dmat(dim*(dim+1)/2+axi,dim*(dim+1)/2+axi)
!
!     En 2D
      real*8   shp2D(3,nnod),Nmat2D(1,nnod)
      real*8   Bmat2D(2,nnod),Nmatel2D(2,2*nnod),Bmatel2D(4,2*nnod)
      real*8   Bmatplano2D(3,2*nnod),r,r0,Bmataxi2D(4,2*nnod)
!     En 1D
      real*8   shp1D(2,nnod), Nmat1D(1,nnod)
      real*8   Bmat1D(1,nnod), Nmatel1D(1,nnod), Bmatel1D(2,nnod)
!
!     Inicializacion de variables
      xjac       =   0.d0
      chi        =   0.d0
      eta        =   0.d0
      zita       =   0.d0
      dx         =   0.d0
      shp2D      =   0.d0
      Nmat2D     =   0.d0
      Bmat2D     =   0.d0
      Nmatel2D   =   0.d0
      Bmatel2D   =   0.d0
      Bmatplano2D=   0.d0
      Bmataxi2D  =   0.d0
      shp1D      =   0.d0
      Nmat1D     =   0.d0
      Bmat1D     =   0.d0
      Nmatel1D   =   0.d0
      Bmatel1D   =   0.d0
      Cmat       =   0.d0
      Dmat       =   0.d0
      matriz     =   0.d0
      r0         =   0.d0      
!
!   Calculo de los puntos de Gauss
      call bgauss(sg,wg)
!
!   Bucle para cada punto de Gauss
      if(dim.eq.2)then
         do i = 1,2
            chi = sg(i)
            do j = 1,2
               eta = sg(j)
!
!            Se calculan las funciones de forma
               call f_forma2D(chi,eta,x,shp2D,xjac)
!
!           Se obtienen las matrices de calculo
               call Matrices_de_calculo2D(shp2D,Nmat2D,Bmat2D,Nmatel2D,Bmatel2D)
!
!           Se calcula el diferencial de la integral
               dx = wg(i)*wg(j)*xjac
!
!           Se calcula la matriz de masa para un solo gdl
               matriz = matmul(transpose(Nmat2D),Nmat2D)
!
!           Se llevan a cabo las operaciones de sumatoria
               Cmat = Cmat + matriz * dx
            enddo
         enddo
!   Bucle para cada punto de Gauss en 1D
      else if (dim.eq.1) then
         do i = 1, 2
            chi = sg(i)
!
!           Se calculan las funciones de forma en 1D
            call f_forma1D(chi, x, shp1D, xjac)
!
!           Se obtienen las matrices de cálculo en 1D
            call Matrices_de_calculo1D(shp1D, Nmat1D, Bmat1D, Nmatel1D, Bmatel1D)
!
!           Se calcula el diferencial de la integral
            dx = wg(i) * xjac
!
!           Se calcula la matriz de masa para un solo gdl
            matriz = matmul(transpose(Nmat1D), Nmat1D)
!
!           Se llevan a cabo las operaciones de sumatoria
            Cmat = Cmat + matriz * dx
         enddo
      end if
!
      return
      end

      subroutine vector_reaccion(u,ndofel,de,x,jelem,t,Cvec)
!
      use reading_utils
!
!     Variables de entrada
      integer  ndofel,jelem
      real*8   x(dim,nnod),t!,de(8)
!
!     variables de salida
      real*8   Cvec(1,nnod)
!     Variables de la subrutina
!     Generales
      integer  i,j,k,l,m,n
      real*8   wg(2),sg(2),vector(1,nnod)
      real*8   chi,eta,dx,xjac
      real*8   Dmat(dim*(dim+1)/2+axi,dim*(dim+1)/2+axi)
!
!     En 2D
      real*8   shp2D(3,nnod),Nmat2D(1,nnod)
      real*8   Bmat2D(2,nnod),Nmatel2D(2,2*nnod),Bmatel2D(4,2*nnod)
      real*8   Bmatplano2D(3,2*nnod),r,r0,Bmataxi2D(4,2*nnod)
!     En 1D
      real*8   shp1D(2,nnod), Nmat1D(1,nnod)
      real*8   Bmat1D(1,nnod), Nmatel1D(1,nnod), Bmatel1D(2,nnod)
!
!     Inicializacion de variables
      xjac       =   0.d0
      chi        =   0.d0
      eta        =   0.d0
      zita       =   0.d0
      dx         =   0.d0
      shp2D      =   0.d0
      Nmat2D     =   0.d0
      Bmat2D     =   0.d0
      Nmatel2D   =   0.d0
      Bmatel2D   =   0.d0
      Bmatplano2D=   0.d0
      Bmataxi2D  =   0.d0
      shp1D      =   0.d0
      Nmat1D     =   0.d0
      Bmat1D     =   0.d0
      Nmatel1D   =   0.d0
      Bmatel1D   =   0.d0
      Cvec       =   0.d0
      Dmat       =   0.d0
      vector     =   0.d0
      r0         =   0.d0      
!
!   Calculo de los puntos de Gauss
      call bgauss(sg,wg)
!
!   Bucle para cada punto de Gauss
      if(dim.eq.2)then
         do i = 1,2
            chi = sg(i)
            do j = 1,2
               eta = sg(j)
!
!            Se calculan las funciones de forma
               call f_forma2D(chi,eta,x,shp2D,xjac)
!
!           Se obtienen las matrices de calculo
               call Matrices_de_calculo2D(shp2D,Nmat2D,Bmat2D,Nmatel2D,Bmatel2D)
!
!           Se calcula el diferencial de la integral
               dx = wg(i)*wg(j)*xjac
!
!           Se calcula el vector reaccion
               vector = Nmat2D
!
!           Se llevan a cabo las operaciones de sumatoria
               Cvec = Cvec + vector * dx
            enddo
         enddo
!   Bucle para cada punto de Gauss en 1D
      else if (dim.eq.1) then
         do i = 1, 2
            chi = sg(i)
!
!           Se calculan las funciones de forma en 1D
            call f_forma1D(chi, x, shp1D, xjac)
!
!           Se obtienen las matrices de cálculo en 1D
            call Matrices_de_calculo1D(shp1D, Nmat1D, Bmat1D, Nmatel1D, Bmatel1D)
!
!           Se calcula el diferencial de la integral
            dx = wg(i) * xjac
!
!           Se calcula el vector de reaccion
            vector = Nmat1D
!
!           Se llevan a cabo las operaciones de sumatoria
            Cvec = Cvec + vector * dx
         enddo
      end if
!
      return
      end
!-----------------------------------------------------------------------------------------
!---------------------------------------------------------------MATRIZ_CONSTANTES_EL------
!
!    Funcion matriz_constantes_e de constantes elásticas
!
!------------------------------------------------------------------------------------------
      subroutine matriz_constantes_el(u,ndofel,de,x,jelem,t,Dmat)
!
      use reading_utils
!
!     Variables de entrada
      integer  ndofel,jelem
      real*8   chi,eta,zita
      real*8   u(ndofel),x(dim,nnod),t,de(8)
!
!     variables de salida
      real*8   Dmat(dim*(dim+1)/2+axi,dim*(dim+1)/2+axi)
!
!     Variables de la subrutina
!     Generales
      integer  i,j,k
      integer  grupo
      real*8   E,nu,d11,d22,d33,d12
      real*8   Dmatb(dim*(dim+1)/2+axi,dim*(dim+1)/2+axi)

      character*276         filename
      character(256)        JOBDIR
      character(256)        JOBNAME
      character(256) :: logFileName
      integer :: logUnit, ierr

      logUnit = 14  ! Define a unit number for the log file
!
!     Definicion del tensor de constantes elasticas
      Dmat   = 0.d0
!
      grupo  = grupoFisico(jelem,2)+1
      E      = propiedades(grupo,1)
      nu     = propiedades(grupo,2)
!
      ! print *, 'jelem', jelem, 'E:', E, 'nu:', nu
      if(dim.eq.2)then
!       Analisis bidimensional
         if(axi.eq.0)then ! Analisis plano
            if(tipo_def.eq.1)then ! Esfuerzo plano
!            Constantes
               d11 = E/(1.d0-nu**2)
               d22 = d11
               d12 = nu*d11
               d33 = 0.5d0*E/(1.d0+nu)   
!
!            Definicion de la matriz de constantes
               Dmat(1,1) = d11
               Dmat(2,2) = d22
               Dmat(1,2) = d12
               Dmat(2,1) = Dmat(1,2)
               Dmat(3,3) = d33
! 
            elseif(tipo_def.eq.2)then ! Deformacion plana
!            Constantes
               d11 = E*(1.d0-nu)/((1.d0-2.d0*nu)*(1.d0+nu))
               d22 = d11
               d12 = d11*nu/(1.d0-nu)
               d33 = 0.5d0*E/(1.d0+nu)   
!
!            Definicion de la matriz de constantes
               Dmat(1,1) = d11
               Dmat(2,2) = d22
               Dmat(1,2) = d12
               Dmat(2,1) = Dmat(1,2)
               Dmat(3,3) = d33
! 
            end if
!
!       Axisimetrico
         elseif(axi.eq.1)then ! Analisis axisimetrico
!         Constantes
            d11 = E*(1.d0-nu)/((1.d0-2.d0*nu)*(1.d0+nu))
            d22 = d11
            d12 = d11*nu/(1.d0-nu)
            d33 = 0.5d0*E/(1.d0+nu)   
!         Definicion de la matriz de constantes
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
!
      elseif(dim.eq.1)then
!       Analisis unidimensional
            Dmat(1,1) = E
      end if
!
      return 
      end
!------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------oute-----------
!       
!       Funcion outsigma()
!
!
!       Funcion de calculo de las deformaciones requeridas para los calculos de las 
!       matrices restantes, es llamado al inicio de cada incremento
!
!       Diccionario
! -----------------------------------------------------------------------------------------
!
!      Funcion sin argumentos:
!
!
!       Funcion sin argumentos que permite escribir los datos requeridos para 
!       calculos alternos, por ejemplo, las deformaciones
!
! -----------------------------------------------------------------------------------------
      subroutine outsigma()

         use reading_utils
         include 'ABA_PARAM.INC'
      
         integer :: jelem, j, k, grupo
         real*8  :: xjac, radio, centro
         real*8  :: EsMax, EsMin, sigma_zz, Esf_VM, Esf_Hid, thao_oct, OI
         real*8  :: E, nu, BG, CMI
         real*8  :: x(dim, nnod)                   ! nodal coordinates
         real*8  :: shp(3, nnod)                  ! shape functions
         real*8  :: Nm(1, nnod), Bm(2, nnod)
         real*8  :: Ne(2, 2*nnod), Be(4, 2*nnod)
         real*8  :: DESP2D(2*nnod)                ! displacements
         real*8  :: Dmat2(3, 3)
         real*8  :: DEFORM(3), ESF(3)
         real*8  :: CEN(ndofdiff)
   
         ! integer :: jelem, j, k
         ! real*8  :: x(1, nnod)              ! Nodal coordinates (only x)
         real*8  :: shp1D(2, nnod)
         real*8  :: N1D(1, nnod), B1D(1, nnod)
         real*8  :: Ne1D(1, nnod), Be1D(1, nnod)
         real*8  :: DESP1D(nnod)            ! 1D displacements
         real*8  :: eps1D, sig1D
         real*8  :: Dmat(1,1)
         real*8  :: de(8)
         real*8  :: t
      
         
         if (dim == 2) then
            resElem = 0.0d0
         
            do jelem = 1, NELEMS
               ! Reset per element
               x     = 0.0d0
               Dmat2 = 0.0d0
               DESP2D = 0.0d0
               CEN   = 0.0d0
         
               do j = 1, nnod
                  k = conectividades(jelem, j+1)
                  do i = 1,dim
                     x(i,j) = nodes(k,i)
                  end do

                  do i = 1,ndofdiff
                     CEN(i) = CEN(i) + Nm(1,j)*resNod(k,ndofdiff+i-1)
                  end do
               end do
         
               call f_forma2D(0.0d0, 0.0d0, x, shp, xjac)
               call Matrices_de_calculo2D(shp, Nm, Bm, Ne, Be)
               call matriz_constantes_el(U, ndofel, de, x, jelem, t, Dmat2)
         
               do j = 1, nnod
                  k = conectividades(jelem, j+1)
                  
                  do i=1,dim
                     DESP2D(2*j-1+i-1) = resNod(k,i)
                  end do
               end do
         
               DEFORM = matmul(Be(1:3,:), DESP2D)
               ESF    = matmul(Dmat2, DEFORM)
         
               ! Mohr-circle
               radio  = sqrt( ((ESF(1)-ESF(2))/2.d0)**2 + ESF(3)**2 )
               centro = 0.5d0*(ESF(1)+ESF(2))
               EsMax  = centro + radio
               EsMin  = centro - radio
         
               if (tipo_def == 2) then
                  grupo = grupoFisico(jelem,2)+1
                  E = propiedades(grupo,1)
                  nu = propiedades(grupo,2)
                  sigma_zz = nu*(ESF(1) + ESF(2))
               else
                  sigma_zz = 0.0d0
               end if
         
               Esf_VM   = sqrt( ((EsMax-EsMin)**2 + (EsMin-sigma_zz)**2 + (sigma_zz-EsMax)**2) / 2.d0 )
               Esf_Hid  = (EsMax + EsMin + sigma_zz) / 3.d0
               thao_oct = Esf_VM * sqrt(2.d0)/3.d0
               OI  = thao_oct + kOI * Esf_Hid
               BG  = CEN(1) - 0.5*CEN(2)
               CMI = 0.5*BG + OI
         
               resElem(jelem, 1) = DEFORM(1) !E11
               resElem(jelem, 2) = DEFORM(2) !E22
               resElem(jelem, 3) = 0.0 !E33 (Not implemented)
               resElem(jelem, 4) = DEFORM(3) !E12
               resElem(jelem, 5) = ESF(1) !S11
               resElem(jelem, 6) = ESF(2) !S22
               resElem(jelem, 7) = sigma_zz !S33
               resElem(jelem, 8) = ESF(3) !S12
               resElem(jelem, 9) = Esf_VM !S_Mises
               resElem(jelem, 10) = Esf_Hid !S_Hyd
               resElem(jelem, 11) = thao_oct !S_Oct
               resElem(jelem, 12) = OI !OI
               resElem(jelem, 13) = CEN(1) !C1
               resElem(jelem, 14) = CEN(2) !C2
               resElem(jelem, 15) = CEN(3) !C3
               ! resElem(jelem, 15) = BG
               ! resElem(jelem, 16) = CMI
            end do
         else if (dim == 1) then
            do jelem = 1, NELEMS
               x      = 0.0d0
               DESP1D = 0.0d0
               matriz = 0.0d0
         
               ! Get nodal coordinates for the element
               do j = 1, nnod
                  k = conectividades(jelem, j+1)
                  x(1, j) = nodes(k, 1)
               end do
      
               ! Initialize material properties
               call matriz_constantes_el(U, ndofel, de, x, jelem, t, Dmat)
      
               ! Compute shape functions and Jacobian
               call f_forma1D(0.d0, x, shp1D, xjac)
      
               ! Compute calculation matrices
               call Matrices_de_calculo1D(shp1D, N1D, Ne1D, B1D, Be1D)
         
               ! Compute strains and stresses
               eps1D = dot_product(Be1D(1, :), DESP1D(1:nnod))
               sig1D = Dmat(1, 1) * eps1D
         
               ! Store results
               ! resElem(jelem, 1) = eps1D     ! Strain
               ! resElem(jelem, 2) = sig1D     ! Stress
            end do
         end if
   
         end subroutine outsigma
!
!------------------------------------------------------------------------------------------
!-------------------------------------------------------------accumulateResults------------
!       
!       Funcion accumulateResults()
!
!
!       Funcion para acumular los resultados de los elementos
!       
!
! -----------------------------------------------------------------------------------------
!
!      Argumentos: KINC, paso
!
      subroutine accumulateResults(KINC)
!  
      use reading_utils
      include 'ABA_PARAM.INC'
!
      integer    i,j,KINC
!
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

      enddo
      return
      end
!------------------------------------------------------------------------------------------
!-------------------------------------------------------------detectBorders--------
!
!
!       Funcion detectBorders()
!
!
!----------------------------------------------------------------------------
      subroutine detectBorders(borderVectorElem,borderVectorNod)
!
      use reading_utils
!
      logical :: belongsToLimit, condition
      integer :: limitGroup, limitElement, index1, index2
      real*8 :: borderVectorElem(NELEMS,1), borderVectorNod(NUMNODE,1)
      integer :: i,j

      nodoBorde = .false.
      elementoBorde = .false.
      
      do i=1,NELEMS
         belongsToLimit = .false.
         
         if (grupoFisico(i,2) == 2) then
            do j=1,nnod
               limitElement = adyacencias(i,j+1)
               limitGroup = grupoFisico(limitElement,2)
               condition = (limitGroup /= 2) .and. (limitElement/= 0) .and. (limitGroup /= 1)
               
               if (condition) then
                  index1 = conectividades(i,j+1)
                  index2 = conectividades(i,mod(j,nnod)+2)
                  nodoBorde(index1,1) = nodoBorde(index1,1) .or. condition
                  nodoBorde(index2,1) = nodoBorde(index2,1) .or. condition
               end if

               belongsToLimit = belongsToLimit .or. condition

            enddo
         end if
         
         elementoBorde(i,1) = belongsToLimit

         if (belongsToLimit) then
            borderVectorElem(i, 1) = 1.0d0
         else
            borderVectorElem(i, 1) = 0.0d0
         end if

      end do

      do i = 1, NUMNODE
         if (nodoBorde(i,1)) then
            borderVectorNod(i, 1) = 1.0d0
         else
            borderVectorNod(i, 1) = 0.0d0
         end if
      end do
      
      end subroutine detectBorders
!------------------------------------------------------------------------------------------
!-------------------------------------------------------------MPC---------------------
! Degree of freedom constraint subroutine
!--------------------------------------------------------------------------------------------
      SUBROUTINE MPC(UE,A,JDOF,MDOF,N,JTYPE,X,U,UINIT,MAXDOF,
     1 LMPC,KSTEP,KINC,TIME,NT,NF,TEMP,FIELD,LTRAN,TRAN)
!
      use reading_utils
      INCLUDE 'ABA_PARAM.INC'
!
      DIMENSION A(N),JDOF(N),X(6,N),U(MAXDOF,N),UINIT(MAXDOF,N),
     1 TIME(2),TEMP(NT,N),FIELD(NF,NT,N),LTRAN(N),TRAN(3,3,N)

      real :: factor

!      user coding to define UE, A, JDOF, and, optionally, LMPC
      IF (.false.) THEN
         JDOF(2) = 20 ! Make unit variable the independent degree of freedom
         A(1) = 1.0d0 ! Coefficient for unit variable
         JDOF(1) = 13 ! Make C1 the dependent degree of freedom
         factor = 2.0d0
         A(2) = -factor ! Coefficient for C2
         UE = factor*U(20,1)
   
      ELSEIF (.false.) THEN
         JDOF(2) = 13 ! Make C1 the independent degree of freedom
         A(1) = 1.0d0 ! Coefficient for C1
         JDOF(1) = 14 ! Make C2 the dependent degree of freedom
         factor = 1.0d0
         A(2) = -factor ! Coefficient for C2
         UE = factor*U(13,1)
      ELSE
         LMPC = 0
      ENDIF

      RETURN
      END