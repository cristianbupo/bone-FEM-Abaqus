      integer, parameter :: NUMNODE=5802, NELEMS=5680, dim=2, nnod=4
      real*8 nodes(NUMNODE, dim)
      parameter(axi=0,tipo_def=2)
      real*8, parameter :: kOI=0.5
      integer conectividades(NELEMS, nnod+1)
      integer grupoFisico(NELEMS, 2)
      parameter(filasContorno1=14)
      parameter(filasContorno2=14)
      integer contorno1(filasContorno1,6),contorno2(filasContorno2,6)
      real*8 resNod(NUMNODE, 2)
      real*8 resElem(NELEMS, 12)
      real*8 myprops(8)
      common resNod, resElem
      common myprops
      common nodes
      common conectividades
      common grupoFisico
      common contorno1, contorno2
