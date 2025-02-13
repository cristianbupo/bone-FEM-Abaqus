      integer, parameter :: NUMNODE=26, NELEMS=18, dim=2, nnod=4
      real*8 nodes(NUMNODE, dim)
      parameter(axi=0,tipo_def=2)
      real*8, parameter :: k_OI=0.5
      integer conectividades(NELEMS, nnod+1)
      integer grupoFisico(NELEMS, 2)
      integer, parameter :: filasContorno1=1, filasContorno2=1
      integer, parameter :: filasContorno3=1, filasContorno4=1
      integer contorno1(filasContorno1, 6), contorno2(filasContorno2, 6)
      real*8 resNod(NUMNODE, 2)
      real*8 resElem(NELEMS, 12)
      real*8 myprops(8)
      common resNod, resElem
      common myprops
      common nodes
      common conectividades
      common grupoFisico
      common contorno1, contorno2
      common contorno3, contorno4
