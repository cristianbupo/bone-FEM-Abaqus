      integer, parameter :: NUMNODE=566, NELEMS=528, dim=2, nnod=4
      real*8 nodes(NUMNODE, dim)
      real*8, parameter :: k_OI=0.5
      integer conectividades(NELEMS, nnod+1)
      integer grupoFisico(NELEMS, 2)
      common nodes
      common conectividades
      common grupoFisico
