      parameter(NUMNODE=591, NELEMS=550, dim=2, nnod=4)
      real*8 nodes(NUMNODE, dim)
      integer conectividades(NELEMS, nnod+1)
      integer grupoFisico(NELEMS, 2)
      common nodes
      common conectividades
      common grupoFisico
