c program to test the function 'correlation' from correlation.f
c MM 1999-2-11

      program testcorrelation

      include 'parameters.f'
      integer*2 iz0s(mxxy), nin(mxxy)
      real*8  zsbnd(mxxy, maxbnd)
      common/mats/ zsbnd, iz0s, nin
      
      real  c
      integer sgn
      
      nin(1) = 2
      zsbnd(1, 1) = 0.5
      zsbnd(1, 2) = 0.75
      iz0s(1) = 1
      
      nin(2) = 2
      zsbnd(2, 1) = 0.5
      zsbnd(2, 2) = 0.60000
      iz0s(2) = 1
      

      sgn = iz0s(1) * iz0s(2)
      
      c = correlation(1,0, nin(1), 0.,  
     ,                1,1, nin(1), 0.5,   
     ,                sgn, 0.5)

 100  format ('correlation: ', f10.6)
      write (6, 100) c
      
      
      end ! program testcorrelation
      
      include 'correlation.f'
