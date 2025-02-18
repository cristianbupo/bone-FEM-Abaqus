      integer, parameter :: NUMNODE=5802, NELEMS=5680, dim=2, nnod=4
      integer, parameter :: listNElementLoads(6)=(/25 ,25 ,25 ,25 ,25 ,25/)
      integer, parameter :: maxNElementLoads=25
      integer, parameter :: order2(2) = (/ 2, 1 /)
      real*8, parameter :: propiedades(2,2) = reshape((/500.0, 0.2, 6.0, 0.47/),
     1 (/2, 2/),  order=order2)
      integer, parameter :: axi=0, tipo_def=2 
      integer, parameter :: filasContorno1=14
      integer, parameter :: filasContorno2=14

C     E, nu
C     500.0, 0.2
C     6.0, 0.47

C     Contains the magnitude of the loads
      real*8 elementLoads(maxNElementLoads) 
C     Contains the tag of the elements and the face where the load is applied
      integer elementFaces(maxNElementLoads, 2)

      real*8 nodes(NUMNODE, dim)
      integer conectividades(NELEMS, nnod+1)
      integer grupoFisico(NELEMS, 2)
      integer contorno1(filasContorno1, 6), contorno2(filasContorno2, 6)
      real*8 resNod(NUMNODE, 2)
      real*8 resElem(NELEMS, 12)

      common resNod, resElem
      common nodes
      common conectividades
      common grupoFisico
      common contorno1, contorno2
      common elementLoads, elementFaces