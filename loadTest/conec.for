      parameter(NUMNODE=4, NELEMS=1, dim=2, nnod=4)
      real*8 nodes(NUMNODE, dim)
      integer conectividades(NELEMS, nnod+1)
      common nodes
      common conectividades
