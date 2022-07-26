c
c 1999-2-9: Began to add the correlation-calculation in subroutine core.
c
c
c
c
c 1999-2-8: Dramaticly decreased the number of 'continue's.
c           in other words, replaced a lot of goto-loops
c           with for example a do while(.true.) construction
c           Goto's are only used to break out of loops, except in the procedure mcw, which is too complex.           
c
c Version: 1999-2-4
c  
c   MM:
c   Hardly changed anything until now. Mainly added commence and lay out.
c   
c
c
c for wolff cluster 3d s=1/2 ising model anisotropic limit
c
c select neighbors via continuous intervals to be skipped
c nrannr,irn in common block 
c periodic boundary conditions                                          
c program anlim.f (based on w1n.f)
c
c core: wn in common; vz=zsbnd ...
c a little bit cryptic (MM)      
c
      program anlim 
      
      implicit real*8 (a-h, o-z)                                         
      parameter (nprlll = 1)                                              
c      dimension itask(nprlll)                             

      do while(.true.)      
c        read from standard input
         print *,' n1,n2,a3,ncycle,ntoss,numint,nint,nr,itmax'             
         read *, n1, n2, a3, ncycle, ntoss, numint, nint, nr, itmax                 
         
         itmax=itmax * 3600
         
         if (n1.eq.0) stop      ! that's how we break out 
         
c        read more from standard input 
         print *,' tpar,fac'                                               
         read *, tpar, fac        ! temperature    
         
c        open the output file
         open(8,file='anli.dat',access='append')                          
         
c        main loop, which is not a loop right now, because nprlll = 1 (MM)
         do iprlll = 1 , nprlll                 
            call bain(n1, n2, a3, ncycle, ntoss, numint, nint, nr,itmax,
     ,           tpar, fac)    
         enddo ! iprlll
         
         close(8) !  close the output file
      enddo    ! while(.true.)  
      end   ! program anlim                                             


c----------------------------------------------------------                 
c  main program
c  (why it's bain, not main, I don't know.(MM)
c----------------------------------------------------------
      subroutine bain(n1, n2, a3, ncycle, ntoss, numint, nint, nr,itmax,
     ,                tpar, fac)
c     program for simulation of anisotropic limit of ising model 
      implicit real*8(a-h, o-z)                                          
      save
      character*80 ident                                                
c      dimension idate(14)                                               
      parameter (npmax = 2000)                       
                    
      parameter (mxnx = 40, mxny = 40, maxbnd = 150)                     
      parameter (mxxy = mxnx * mxny)         
      parameter (mxrn = mxxy * maxbnd * 2)
      common /lookuc/ t2p, maxr, maxi
      data ident /'anli v0'/                                            

      print *, ident                                                     
      xt0   = ctijd(0.) ! 0. is dummy, doesn't have meaning 
      tmax  = itmax - 10
      ntemp = 0
                                                           
      if(n1.gt.mxnx.or.n2.gt.mxny) then                   
         print *,' main: at least one lattice dimension too big'        
         stop                                                           
      endif                                                             

      mcstps = numint * nint                                                
      npart  = numint / 2                                                    
      if (npart.gt.npmax) npart = npmax                                   
      if (npart.le.10) npart = 1                                          

      npp   = mcstps / npart                                                  
      nip   = npp / nint                                                      
      npp   = nint * nip                                                      
      mcstps= npart * npp                                                  
      ncount= mcstps / nint                                                
c     104   continue                                                         
      ntemp = ntemp + 1                                                     

c     initialise: 
      call carlo(tpar,fac,n1,n2,a3,ntemp,nr)       

      vcu = n1 * n2 * a3                                                      

c     write to standard output what we're doing (MM)
      write(6,900) n1,n2,a3                                             
 900  format('program for mc simulation  on a',  ! added quotes and took away 32h(MM) 
     *     i4, 1h*, i4,1h*, f12.6, 8h lattice)                                      
      print 901,tpar,fac                                              
 901  format(' tpar= ', f11.6,
     ,     ' factor=', f11.6)
      print 902, ncycle, mcstps, ntoss                                     
 902  format(1x, i3, ' consecutive runs of ', i10, ' Wolff steps', 1x,
     *     'after', i10, ' times nint steps for equilibration'/)
      print 924, nr                                                      
 924  format(' nr=', i8/)                                                
      print 903, nint                                                    
 903  format(' data are taken at intervals of ', i3, ' mcs/site'/)        


      do itoss = 1, ntoss
         call monte(nint)
      enddo ! itoss

      do nsigh = 1, ncycle                                               
         write(8,120) ident,n1,numint,nint,tpar,nr
 120    format(a8,'; n1,numint,nint,knn,s:     ',i3,i10,i3,f12.6,12x,i4)

         call setpqr                                                       

c     in xt? are stored several times
         xt3 = ctijd(0.)                                                   
         xt1 = xt3 - xt0                                                       
         print 554,xt1                                                     
 554     format(20x,'setup time (sec.)=',f10.6)                            
         tmc = 0                                                             
         tan = 0                                                             
         print 999                                                         
 999     format(/70('*')/)                                             
         nc=0                                                              
         do ip = 1, npart                                                  
            inp = ip                                                            
            do iip = 1, nip                                                    
               nc = nc + 1
                                                           
               call monte(nint)                          

               xt2 = ctijd(0.)
               tmc = tmc + xt2 - xt3 

               call addpqr(n1, n2, ncup)                        

               xt3 = ctijd(0.)                                            
               tan = tan + xt3 - xt2 
               if ( (tmax - xt3 + xt0).lt.0) goto 1111

            enddo ! iip
            call part(inp, nc, ncup)    
                                        
         enddo ! ip                                                         

         goto 1112                                                         

 1111    ncount = nc                                                         
         
         call part(inp, nc, ncup)                                            
         
         mcstps = ncount * nint     
         write(6,*) 'dat gaat niet goed'

         write(6,1120)mcstps, xt3 - xt0                                           
 1120    format(' ****** emergency stop,mcs=',i9, 3h = ,f11.3,' sec.') ! added quotes, took away 37h (MM)

 1112    continue                                                          

         write(6,1130) xt3 - xt0                                           
 1130    format(' mc stop at t=',f12.2,' sec.')    

c     call lpsb                                                         
         tmc = (tmc * 1.d + 6 / mcstps)/ vcu
         tan = (tan * 1.d + 6 / ncount)/ vcu

         call pqr(ncup, ncount, npart, n1, n2, a3)                             

         if ( (tmax - xt3 + xt0).lt.0) stop                                         

      enddo ! do nsigh = 1 , ncycle 

      print 999             ! format(/70('*')/)

      print 553,maxr, mxrn, maxi, maxbnd
 553  format('max nr random numbers used:',i10,' allowed:',i10/
     /     'max nr  cluster boundaries:',i10,' allowed:',i10)

      xt4 = ctijd(0.)                                                   
      anal= xt4 - xt3                                                      
      t4  = xt4 - xt0                                                        

      print 555, tmc, tan, anal, t4                                         
 555  format('0mc simulation  time (microseconds/mc step/site)=',f20.6/ 
     *     /' summation over lattice (microseconds/time/site)=',f20.6//    
     *     ' analysis (seconds) =',f10.2// ' total time (sec.)=',f10.2)    

      return                                                            
      end    ! bain


c------------------------------------------------------ 
c carlo
c Initialisation of Monte-Carlo algorithm.(MM)
c
      subroutine carlo(tpar, fac, m1, m2, aa, ntmp, nr)
c     version for periodic ising model                                  
      implicit real*8(a-h, o-z)                                          
      save
      parameter (eps = 2.d-15)       

c     maximum sizes(MM): 
      parameter (mxnx = 40, mxny = 40, maxbnd = 150)                     
      parameter (mxxy = mxnx * mxny)                                   
      parameter (mxrn = mxxy * maxbnd * 2)
      
      parameter (ilut=6)                                    
      common /lookup/ rmul, radd, rmul2, radd2, rmul3, radd3
      common /lookuc/ t2p, maxr, maxi
      common/nran/ nrannr, irn(mxrn)
      integer*2 iz0s(mxxy), nin(mxxy)
      common/mats/ zsbnd(mxxy, maxbnd), iz0s, nin
      integer*2 ngc(4, mxxy)
      common /ngcb/ ngc
      parameter (tp16 = 2**16, tm32 = 1.0d0/(tp16 * tp16))
      parameter (tp32m1=tp16 * tp16 - 1, maxint = 2147483647)
      common/dats/a3, wn, n1, n2, n12

c     check lattice size:
      if ((m1.gt.mxnx) .or. (m2.gt.mxny)) then         
         write(6,9) m1, m2                                                  
 9       format('carlo disagrees with specified lattice size',2i3)       
         stop 
      endif

      a3  = aa ! size in z-direction
      t2p = tpar / 2.0d0  ! K/2
      maxr= 0 ! # of used random numbers
      maxi= 0 ! maximum number of interfaces, of course we begin at 0.

c     size of system(MM)
      n1 = m1                                                             
      n2 = m2                  
      n12= n1 * n2              

      wn = 1.0d0 / (a3 * n12) ! 1/size of system

c     This loop fills ngc.
c     this is a table in which the neigbors of every site are stored.
      do index = 1 , n12                                                

         ! x and y coordinate of the site
         iy   = (index - 1) / n1  + 1
         ix   = index - (iy - 1)* n1  

         ! x+1 and y+1
         ixp  = ix + 1 - (ix / n1)* n1 
         iyp  = iy + 1 - (iy / n2)* n2

         ! x-1 and y-1
         ixm  = ix - 1 + ( (n1 - ix + 1)/n1 )* n1 
         iym  = iy - 1 + ( (n2 - iy + 1)/n2 )* n2

         ! 4 neighbors:
         ngc( 1,index) = (iy  - 1)* n1 + ixm
         ngc( 3,index) = (iy  - 1)* n1 + ixp
         ngc( 2,index) = (iym - 1)* n1 + ix 
         ngc( 4,index) = (iyp - 1)* n1 + ix 
         
      enddo ! index

c     these things are needed by the random-number generator:
      radd  = 0.5d0 * n12 + 1.0d0               ! N/2 + 1
      rmul  = tm32  * (1.0d0 - 1.0d-15 ) * n12  ! 2^-32 N
      radd2 = 0.5d0 * a3                        ! l_z / 2
      rmul2 = tm32  * (1.0d0 - 1.0d-15 ) * a3   ! 2^-32 l_z
      radd3 = 0.5d0                             ! 1/2
      rmul3 = tm32  * (1.0d0 - 1.0d-15 )        ! 2^-32


      write(6,*)' wolff mc of s=1/2 anisotropic ising model'   
      
      if(ntmp.le.1) then         
c     `ntmp' probably means `tempory integer'
c     it seems to be always .gt.1, exept for the first time.(MM)

         iseed = iabs(nr) + 1    
c        so 'nr' is what is called the seed from outside (MM)
         
c        initialize random number generator:
         call ransi(iseed)
         
         if(nr.lt.0) then       ! nr = `seed'
            do ij = 1, n12
               iz0s(ij) = 1
               nin(ij)  = 0
            enddo ! ij
         else                   ! nr >= 0
            nrannr = n12
            
            call ransr()    
            
            do ij = 1, n12                                                   
               rnd = irn(ij) * rmul3  + radd3
               ipm = rnd + rnd
               ipm =-1 + ipm + ipm                                                  
               iz0s(ij) = ipm                                              
               nin(ij)  = 0
            enddo ! ij                                                        
         endif  ! nr .lt.0      

      endif ! ntmp .le. 1

      nrannr = n12 * maxbnd * 2

      return                                                            
      end ! carlo                          

c----------------------------------------------------- 
c monte
c This subroutine only calls mcw (MM).
c     
      subroutine monte(nsteps)                
c     perform nsteps mc steps per site                                      
c     version for arbitrary model                                       

c     everything is implicitly a double except i,j,k..n.
      implicit real*8 (a-h, o-z)    
      
c     again the maximal sizes (MM):
      parameter (mxnx = 40, mxny = 40, maxbnd = 150)                     
      parameter (mxxy = mxnx * mxny)   
      
c     do 2 nn=1,nsteps                                 
      call mcw(nsteps)                           
c     2    continue                                                          
      return                                                            
      end  ! monte                                                             

c----------------------------------------------------------
c mcw  = monte
c This is the heart of the program, here happen
c the real Monte-Carlo simulations (MM)

      subroutine mcw(nsteps)                                     
c     anisotropic wolff cluster
      implicit real*8 (a-h, o-z)
      
c     yet another time the maximal sizes:(MM)
      parameter (mxnx = 40, mxny = 40, maxbnd = 150)                     
      parameter (mxxy = mxnx * mxny)                                   
      parameter (mxrn = mxxy * maxbnd * 2)

      common/nran/ nrannr, irn(mxrn)
      integer*2 iz0s(mxxy), nin(mxxy)
      common/mats/ zsbnd(mxxy, maxbnd), iz0s, nin
      common/dats/ a3, wn, n1, n2, n12
      common /lookuc/ t2p, maxr, maxi
      integer*2 ngc(4, mxxy)
      common /ngcb/ ngc
      common /lookup/ rmul, radd, rmul2, radd2, rmul3, radd3
      parameter (tp16 = 2**16, tm32 = 1.0d0 / (tp16 * tp16))

c     stacks are a kind of arrays, seemingly (MM)
      integer*2 ijstack(mxxy * maxbnd)
      real*8     zstack(mxxy * maxbnd)
      real*8     wstack(mxxy * maxbnd)

c     logical prt
      save

c     the main loop of this subroutine
      do isteps = 1, nsteps
c     prt=(isteps.eq.17)
         
         call ransr() ! draw a set of random numbers

c     choose origin
         nrannr = 1 ! we use one randomnumber 
         ik     = irn(nrannr) * rmul + radd
         nrannr = 2 ! and the second for the z-direction
         zk     = irn(nrannr) * rmul2 + radd2


c     determine sign of cluster (before flip)
         icsp = iz0s(ik) ! spin on z = 0
         jin  = 0        ! count of interfaces 

         do while (.true.)
            if (jin.ge.nin(ik)) goto 181 !break, we are at the end of the sample
            jp1 = jin + 1 ! why is there a need for jp1?
            
            if (zsbnd(ik, jp1).gt.zk) goto 181 !break, this interface is too far
            jin  = jp1
            
            icsp = -icsp     ! we found an interface, flip the spin
         enddo  ! while(.true.)
 181     continue ! break out

c     now we know the sign of the cluster, it's icsp

         nstack = 0 ! we start with an empty stack
         cw1    = 0.0d0 ! a cluster width, before we do anything, it's 0

         nrannr = 3 ! and random number 3 is not yet used.
         brd =  rmul3 * irn(nrannr) + 0.5d0 ! 0.5 .. 1.5 ??
         brd = -dlog(brd)  ! jump in z-direction

         nrannr = 4
         bod = rmul3 * irn(nrannr) + 0.5d0 
         bod = -dlog(bod) * t2p  ! t2p is the paramater tpar / 2 
         cw4 = 4.0d0 * cw1 ! 0, a cluster width times 4. In the beginning it's simply zero.
         abd = 0.0d0       ! 0  like cw4 this is needed for 
c     zk,cw4,ik,abd,icsp,brd,bod defined


c ----------do while loop----------------------------------------------
 104     continue ! This is the start of a while-like loop

c This loop is still rather messy, but this may be because the algororithm simply
c is pretty complex.        
         
c     begin determination horiz. bounds of cluster
c     Two cases are handled seperately:
c       - there are already interfaces on this z-line
c       - there are aren't
         
            if (nin(ik).gt.0) goto 194 ! there are already interfaces, go handling that

            if (brd.ge.a3) then ! there was no interface, and we jumped to whole distance
               cl = 0.0d0
               cw = a3
               brd= brd - a3  ! (a new random number)
               iz0s(ik) = -iz0s(ik) ! so we can simple flip the sign on z=0
               goto 193 ! go to add it to the stack.
            endif

            ! so the jump in our, by the way empty z-direction, was smaller then a3
            ! what happens now?
            dlc = brd 
            dr  = a3 - brd

            nrannr = nrannr + 1 ! we need a new one
            brd = rmul3 * irn(nrannr) + 0.5d0
            brd =-dlog(brd)  ! again a jump in the z-direction

            ! so we have two of them now: dlc and brd (right and left?)

c     write(6,889) irn(nrannr),brd 
            if (brd.ge.dr) then
               cl = 0.0d0      ! cluster-left
               cw = a3         ! cluster-width
               brd= brd - dr
               iz0s(ik) = -iz0s(ik)
               goto 193 ! add to stack
            endif

            nin(ik) = 2 ! We now have two interfaces. Is it possible to have one?
                        ! No, we have periodic boundary conditions, also in the 
                        ! z-direction 

            cl = zk - dlc 
                          ! zk is the z-coordinate of our starting point
                          ! dlc is also a jump

            ! cl is a possible new location of an interface

            if (cl.lt.0.0d0) then ! periodic boundary conditions
               cl = cl + a3  ! c left
               iz0s(ik) = -iz0s(ik)
               zsbnd(ik,1) = zk + brd
               zsbnd(ik,2) = cl
            else
               cr = zk + brd  ! c right
               if (cr.ge.a3) then  ! periodic boundary conditions
                  iz0s(ik) = -iz0s(ik)
                  zsbnd(ik, 1) = cr-a3
                  zsbnd(ik, 2) = cl
               else
                  zsbnd(ik, 1) = cl
                  zsbnd(ik, 2) = cr
               endif               
            endif            
            ! We have now for every case decided where the 2 new interfaces are.

            cw = dlc + brd  ! the size of our new region (isn't it?)
                            ! `w' is probably of `width'
                            ! `c' could stand for `cluster'

c     calculate new break length
            ! so in `brd' there's always a new jump ready to be used.
            nrannr = nrannr + 1
            brd = rmul3 * irn(nrannr) + 0.5d0
            brd = -dlog(brd)

c     write(6,889) irn(nrannr),brd 
            goto 193 ! dealt with case of no interfaces

            
 194        continue ! there where already interfaces
c     what are we going to do in this case:
c     The case is a little different now because
c     - interfaces can also be removed now
c     - 
            
c     determine size Ising cluster

c     dl: distance to left interface
c     dr: distance to right interface
c     pl: position of left interface
c     pr: position of right interface
            
c     LEFT
            if (jin.eq.0) then    ! the point was taken in the first region
               pl = zsbnd(ik, nin(ik)) ! so pl is the location of the most right boundary, which is the one left of the point we've chosen
               dl = zk - pl + a3    ! consider periodic boundary conditions                              
            else
               pl = zsbnd(ik, jin)  ! it now is simply the boundary left of the chosed point
               dl = zk - pl
            endif

c     RIGHT
            if (jin.eq.nin(ik)) then  ! in the last region
               pr = zsbnd(ik, 1)      
               dr = pr + a3 - zk        ! consider periodic boundary conditions
            else
               pr = zsbnd(ik, jin + 1)
               dr = pr - zk
            endif

c     first find left random cluster reach
            if (brd.ge.dl) then  ! jumped over the interface
               cl = pl           
               brd = brd - dl    ! we can reuse the random number like this
c     eliminate left interface
               if (jin.ne.0)  then
                  do ii = jin + 1, nin(ik)
                     zsbnd(ik, ii-1) = zsbnd(ik, ii) ! so eliminating is done by shifting.
                  enddo ! ii
                  jin = jin - 1 ! one interface less.
               else
                  iz0s(ik) = -iz0s(ik) ! of course, but why is the most right interface not removed?
                                       ! there is no need see \/
               endif
               nin(ik) = nin(ik) - 1   ! so the number of interfaces int decreased, possibly throwing out the most right one.
            else ! so not jumped over the interface, there should be added one:
               dl = brd    ! new distance to the left interface.
c     add left interface
               cl = zk - brd   ! its's new position.

               nin(ik) = nin(ik) + 1  ! increase the number of interfaces. (why this is done first? in contradiction to with the elimination)

               if (cl.lt.0.0d0) then  ! consider periodic boundary conditions.
                  cl = cl + a3
                  zsbnd(ik, nin(ik)) = cl  ! add one interface on the right.
                  iz0s(ik) = -iz0s(ik)
               else                  ! didn't jump over boundary
                  jin = jin + 1
                  do ii = nin(ik), jin + 1, -1
                     zsbnd(ik, ii) = zsbnd(ik, ii - 1)  ! insert one interface, they all have to shift.
                  enddo ! ii
                  zsbnd(ik, jin) = cl       ! the value of the new interface
               endif
c              because we jumped over the interface, we have used our random number, make a new one:
               nrannr = nrannr + 1
               brd = rmul3 * irn(nrannr) + 0.5d0
               brd = -dlog(brd) ! `jump' distance.
c     write(6,889) irn(nrannr),brd 
            endif ! jumped over left interface/ not

c--------------------------------
c     analogous for right hand interval
            if (brd.gt.dr) then  ! jumped over right interface
               cr = pr
               brd = brd - dr
               cw = dl + dr
c              eliminate right interface
               if (jin.ne.nin(ik))  then ! consider periodic boundary conditions
                  do ii = jin + 1, nin(ik) - 1
                     zsbnd(ik, ii) = zsbnd(ik, ii + 1)
                  enddo ! ii
               else
                  do ii = 1, nin(ik) - 1
                     zsbnd(ik, ii) = zsbnd(ik, ii + 1)
                  enddo ! ii
                  iz0s(ik) = -iz0s(ik)
               endif
               nin(ik) = nin(ik) - 1  ! one interface left
c     write(6,352) ik,nin(ik),(zsbnd(ik,iii),iii = 1,nin(ik))
c     352  format('el r int: k,nin,zsbnd',2i3,20f6.2)
            else    ! not jumped over right interface
c            add right interface
               dr = brd
               cw = dl + dr 
               cr = zk + brd
               nin(ik) = nin(ik) + 1
               if (cr.gt.a3) then
                  cr = cr - a3
                  do ii = nin(ik), 2, -1
                     zsbnd(ik, ii) = zsbnd(ik, ii - 1)
                  enddo ! ii
                  zsbnd(ik, 1) = cr
                  iz0s(ik) = -iz0s(ik)
               else
                  do ii = nin(ik), jin + 2, -1
                     zsbnd(ik, ii) = zsbnd(ik, ii - 1)
                  enddo ! ii
                  zsbnd(ik, jin + 1) = cr
               endif
c              used random number, make a new one:
               nrannr = nrannr + 1
               brd = rmul3 * irn(nrannr) + 0.5d0
               brd = -dlog(brd)
c     write(6,889) irn(nrannr),brd  
c88   9  format('breakd.:'i12,2f20.14)
            endif ! eliminating or adding right interface.
            
c     end determination horiz. bounds of flipped part
c     write in stack
 193       continue
           nstack = nstack + 1
c     so, we have 3 stacks, all with same size, because we want to store 3 different quantities.
c     namely:
           ijstack(nstack) = ik ! which interface is left of our point 
           zstack(nstack) = cl  ! where it is
           wstack(nstack) = cw  ! and the width of the spin-region in which our point was
c     write(6,888) nstack,ik,cl,cw
c888  format('stackwr:'2i6,2f20.14)

c     Now we are going to deal with the `jumps' in the x-y plane

 192        continue ! return from read from stack

c  THIS is the factual loop of the program.
            do while(.true.)
c    abd is the part of the area we already dealt with.
c    it goes to 4 times the width of the area (cw4) because
c    we have to find perpendicular bounds in 4 direction.
c    bod is the 'jump' size               
               abd = abd + bod  ! abd? bod? cw4?
               
               if (abd.ge.cw4) then
                  bod = abd - cw4
                  goto 161      ! goto read from stack
               endif
               
c  idr: in wich of the four direction we are  making our perp. bounds.
c  cw1: width of the cluster/area

               idr = abd / cw1  ! idr? cw1? how cryptic!

               zk = cl1 + abd - idr * cw1 ! the z of the perp. bound.
               
               if (zk.gt.a3) zk = zk - a3 ! periodic boundary conditions.

               ik = ngc(idr + 1, ij)      ! neighbor 'site' from the table
c     find sign of nb
               isnb = iz0s(ik)
               jin = 0                 
               do while (.true.)
                  if (jin.ge.nin(ik)) goto 189 ! break
                  jp1 = jin + 1
                  if (zsbnd(ik, jp1).gt.zk) goto 189 ! break
                  jin = jp1
                  isnb = -isnb
               enddo            ! while (.true.)            
 189           continue         ! break out                        
c     end find sign

c     first determine interval to next perp bond
c     = draw new random number               
               nrannr = nrannr + 1
               bod = rmul3 * irn(nrannr) + 0.5d0
               bod = -dlog(bod) * t2p
            
c     if it is the right sign then we can go to redetermine the size and to add it to the stack
               if (isnb.eq.icsp) goto 104  
                                ! isnb = sign of neigboring cluster
                                ! icsp = sign of cluster
c     104: just before determining the horizontal bound

            enddo  ! while(.true.)
         
c----------------- Read from stack  ----------------------------------------    
 161     continue 
         if (nstack.ne.0) then
            ij = ijstack(nstack)
            cl1 = zstack(nstack)
            cw1 = wstack(nstack)
            cw4 = 4.0d0 * cw1   ! 
            nstack = nstack - 1
c     write(6,981) isteps,nstack,ij,cl1,cw1
c     981  format ('str isteps,nstack ij cl cw',3i4,2f8.3)
            abd = 0.0d0
            goto 192   ! return
         endif
! else continue to next Wolff iteration  because stack was empty 
c----------------------------------------------------------------------------      

c check on number of used random numbers
         if (nrannr.gt.maxr) then
            maxr = nrannr
            if (nrannr.gt.mxrn) then
               print *,' **** nrannr > mxrn:',nrannr,mxrn
               stop
            endif
         endif
         
      enddo   ! isteps = 1, nsteps , end of main loop of this subroutine 

      return          

      end ! mcw
c---------------------------------------------------------------
c correlation

      include 'correlation.f' ! see in that file

c-----------------------------------------------------------------------------------------
c  zk intoduced everywhere; OK for both par and perp bonds?
c parallel parms: append 1
c vermijd interf parms in par and perp sectie
c interfaces (grain boundaries) everywhere updated?
c check precise bounds of intervals generated by bod and brd

c---------------------------------------------------------------
c core
c
c
      subroutine core(sn,ncor)                                
c     version for cylinder ising model                                  
      implicit real*8 (a-h,o-z)                                         

c     The maximal sizes...(MM)
      parameter (mxnx = 40, mxny = 40, maxbnd = 150)                     
      parameter (mxxy = mxnx * mxny)                                   
      parameter (mxrn = mxxy * maxbnd * 2)
      
      common/nran/ nrannr, irn(mxrn)
      integer*2 iz0s(mxxy), nin(mxxy)
      common/mats/ zsbnd(mxxy, maxbnd), iz0s, nin
      common/dats/a3, wn, n1, n2, n12
      common /lookuc/ t2p, maxr, maxi
      integer*2 ngc(4, mxxy)
      
c     the neighbors: (MM)
      common /ngcb/ ngc

      parameter (mxcr = 9)                                      

      dimension  sn(mxcr)
      save

c     Strange things happen here. 
c     Perhaps this is indeed a good spot to calculate the correlation function.
      ncor = 9                                        
c---------------------------------------------------------------------------
c correlation function (MM)

c calculating the correlation function is only possible with even system sizes:

      if(2 * n1 / 2 .ne. n1 .or. 2 * n2 / 2 .ne. n2) then
         ! correlation function is not possible
      else
         ! it is possible, we're going to do it:

         correlation2 = 0  ! we start at zero.  (2 of 'half system size')

         ! first in the z-direction
         halfa3 = a3 / 2
         do ij = 1, n12
            
                                ! find the sign and interface number half way:
            nij = nin(ij)       ! number of interfaces here
            i_sign_z0 = iz0s(ij) ! sign on z = 0
            
            int_halfway = 0     ! number of the interface left of halfway
            i_sign_halfway = i_sign_z0 ! sign on halfway
                                ! we're going to search these two:
            do while (zsbnd(ij,int_halfway) .le. halfa3 .and. 
     *           int_halfway .lt. nij)
               int_halfway    =   int_halfway + 1 
               i_sign_halfway = - i_sign_halfway
            enddo
                       
c            correlation2 = correlation2 + 
c     +           correlation(ij, 0,           int_halfway,       0d0, 
c     ,                       ij, int_halfway, nij - int_halfway, halfa3,
c     ,                       i_sign_z0 * i_sign_halfway, halfa3)
         enddo                  ! ij

         ! in y-direction
         do ij = 1, n12 / 2 
            correlation2 = correlation2 + something
         enddo
         
         ! in x-direction
         do ij = 1, n12            
            correlation2 = correlation2 + somthing
         enddo
            
        
      endif      
c     if (ncor.eq.9) goto 10
c     sample correlation function over half the system size 
c 10   continue

c-----------------------------------------------------
c     sample magnetization and nr of interfaces
      nit = 0  ! total number of interfaces (so not just at one (x,y), but of the whole sample
      amt = 0  ! probably `magnetisation', beginning with an `a' to make it a real... 

c     n12 is size of system (only x.y of course)
      do ij = 1, n12
         nij = nin(ij)          ! number of interfaces at this (x,y)
         if (nij.gt.maxi) maxi = nij ! so maxi = maximum number of interfaces found until now
c                                      I suppose it's just a check
         nit = nit + nij
         izs = iz0s(ij)        ! spin at z = 0
         amt = amt + a3 * izs  ! multiplied by a3?? Why?

         if (nij.gt.0) vz = zsbnd(ij, nij) 

         do int = 1, nij ! over all interfaces
            z   =  zsbnd(ij, int)
            amt =  amt + (z - vz) * izs
            izs = -izs
            vz  =  z
         enddo    ! int, over all interfaces

      enddo ! ij,  over all (x,y)

c     I suppose that `amt' now is the magnetisation
      
c     another loop of the whole system
      sps = 0.0d0 ! sps?

      do ij1 = 1, n12
         ni1 = nin(ij1)
         do inb = 1, 2 ! 2 neigbors?

            in1 = 0
            in2 = 0

            ij2 = ngc(inb, ij1) ! a neighbor of ij1
            ni2 = nin(ij2)      ! number of interfaces at this neighbor

            isp2 = iz0s(ij1) * iz0s(ij2) ! s(ij,z=0)s(i+-1,z=0)
            sps = sps + a3 * isp2
            isp2 = isp2 + isp2

            if (0.eq.ni1) goto 21 ! may be the use of goto is acceptable here.
            if (0.eq.ni2) goto 22

            do while(.true.)
               if (zsbnd(ij1, in1 + 1).gt.zsbnd(ij2, in2 + 1)) then
                  in2 = in2 + 1
                  sps = sps + isp2 * zsbnd(ij2, in2)
                  isp2 = -isp2
                  if (in2.eq.ni2) goto 22
               else
                  in1 = in1 + 1
                  sps = sps + isp2 * zsbnd(ij1, in1)
                  isp2 = -isp2
                  if (in1.eq.ni1) goto 21
               endif
            enddo ! while(.true.)

 21         continue
            do jn2 = in2 + 1, ni2
               sps = sps + isp2 * zsbnd(ij2,  jn2)
               isp2 = -isp2
            enddo
            goto 25

 22         continue
            do jn1 = in1 + 1, ni1
               sps = sps + isp2 * zsbnd(ij1, jn1)
               isp2 = -isp2
            enddo ! jn1

 25         continue

         enddo ! inb
      enddo ! ij1
c     I've changed some <nr> continue to enddo

c     the things this program calculates: 
c     wn: 1/ size of system. (1/(lx ly az))
      bsn = nit * wn  ! 
      sps = sps * wn  !
      st  = amt * wn  ! magnetisation per `site': m
      s2  = st  * st  ! m^2                                                  

      sn(1) = s2            ! m^2
      sn(2) = s2 * s2       ! m^4
      sn(3) = bsn           ! ?
      sn(4) = bsn * bsn     ! ?^2
      sn(5) = sn(3) * sn(4) ! ?^3
      sn(6) = sn(4) * sn(4) ! ?^4
      sn(7) = sn(1) * sn(3)
      sn(8) = sn(2) * sn(3)
      sn(9) = sps

      return                                                            
      end ! core


c --------------------------------------------------------------
c setpqr
c this subrouting probably sets p q (to zero?) (MM)

      subroutine setpqr                                                 
      implicit real*8 (a - h, o - z)                                         

c     npmax?(MM)
      parameter (npmax = 2000)               

c     the maximal sizes
      parameter (mxnx = 40, mxny = 40, maxbnd = 150)                     
      parameter (mxcr = 9)    
                                  
      double precision p(mxcr),q(mxcr,mxcr),pp(mxcr,npmax)              
      dimension        ipr(npmax)                                              
      common/pqrh/     p, q, pp, ipr                                            

      do k = 1, mxcr                                                     
         p(k) = 0.0                                                          
         do l = 1, mxcr                                                     
            q(k,l) = 0.0                                                        
         enddo ! l
      enddo ! k

      do  k = 1, npmax                                                    
         ipr(k) = 0                                                          
      enddo ! k
      
      return                                                            
      end ! setpqr  

c----------------------------------------------------------
c addpqr
c p,q are things which are greatened every iterations, by means 
c of this funcion.
      subroutine addpqr(n1, n2, ntop) 

      implicit real*8 (a - h, o - z)                                          
      parameter (npmax = 2000)                                           
c     again. 

      parameter (mxnx = 40, mxny = 40, maxbnd = 150)                     
      parameter (mxcr = 9)                                            
      dimension        ipr(npmax)                                              
      double precision sn (mxcr)                                         
      double precision p(mxcr), q(mxcr, mxcr), pp(mxcr, npmax)              
      common/pqrh/     p, q, pp, ipr                                            

      call core(sn, ntop)
      
      do j = 1, ntop                                                     
         p(j) = p(j) + sn(j)                                                   
         do k = 1, ntop
            q(j, k) = q(j, k) + sn(j) * sn(k)                                         
         enddo ! k
      enddo ! j
     
      return                                                            
      end ! addpqr

c--------------------------------------------------------------
c part
c      
      subroutine part(ip,nc,ntop)                                       
 
      implicit real*8(a-h,o-z)                                          

      parameter (npmax = 2000)                                           
      parameter (mxnx = 40, mxny = 40, maxbnd = 150)                     
      parameter (mxcr = 9)                                      

      dimension ipr(npmax)                                              
      double precision p(mxcr), q(mxcr,mxcr), pp(mxcr,npmax)              
      common/pqrh/ p, q, pp, ipr                                            
c in part:
      save

      do i = 1, ntop                                                   
         pp(i, ip) = p(i)                                                     
         p(i) = 0.0d0                                                     
      enddo ! i

      if (ip.eq.1) nv = 0      
      ipr(ip) = nc - nv
      nv = nc
c end
      return                                                            
      end ! part

c -----------------------------------------------------
c pqr
c
c
      subroutine pqr(ncup,ncount,npart,n1,n2,a3)                       
 
      implicit real*8(a-h,o-z)                                          
      logical prt                                                       

      parameter (npmax = 2000)                                            
      parameter (npmxh = npmax / 2)                                          
      parameter (mxnx = 40,mxny = 40,maxbnd = 150)                     
      parameter (mxcr = 9)                                      

      double precision p(mxcr), q(mxcr, mxcr), pp(mxcr, npmax)              
      double precision b(mxcr, mxcr), rp(mxcr, npmxh)                      

      dimension ipr(npmax), jpr(npmxh), ndiv(12)                          
      common/pqrh/ p, q, pp, ipr                                            
      common/cbh/  b                                                      
      data ndiv /2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000,5000,10000/         

      prt = .false.                                                       

      do j = 1, ncup                                                     
         p(j) = 0.0d0
         do i = 1, npart                                                  
            p(j) = p(j) + pp(j, i)                                                
         enddo ! i
         p(j) = p(j) / ncount                                                  
      enddo ! j
      
      do j = 1, ncup                                                     
         do k = 1, ncup                                                     
            q(j, k) = q(j, k) / ncount                                              
            b(j, k) = q(j, k) - p(j) * p(k)                                           
         enddo ! k
      enddo ! j
         
      do i = 1, npart 
         if (ipr(i).ne.0) then            
            do j = 1, ncup                                                    
               pp(j, i) = pp(j, i) / ipr(i)
            enddo ! j
         endif
      enddo ! i
                                                          
      call avdv(npart, ncup, ipr, pp, p, prt, n1, n2, a3)                      

      nrdiv = 1                                                           

      do while(.true.)         
         npdiv = npart / ndiv(nrdiv)                                           
         prt = (nrdiv.eq.1)                                        
         if (npdiv.le.3) return                               
c     this is the only way to exit this subroutine
c     the following is done as long as necessary (MM)
         
         do i = 1, npdiv                                                   
            jpr(i) = 0                                                          
            do j = 1, ncup                                                    
               rp(j,i) = 0.0d0                                                     
            enddo ! j
         enddo ! i
         
         npd = 0                                                             
         ind = 0                                                             
         do i = 1, npdiv                                                   
            ntog = ndiv(nrdiv)                                                  
            do k = 1, ntog                                                    
               ind = ind + 1                                                         
               jpr(i) = jpr(i) + ipr(ind)                                            
               do j = 1, ncup                                                    
                  rp(j, i) = rp(j, i) + pp(j, ind) * ipr(ind)                                
               enddo ! j
            enddo ! k
            if (jpr(i).ne.0) then                       
               npd = npd + 1                                                         
               do j = 1, ncup                                                    
                  rp(j, i) = rp(j, i) / jpr(i)                                            
               enddo ! j
            endif
         enddo   ! i                                                       
         
         call avdv(npdiv, ncup, jpr, rp, p, prt, n1, n2, a3)                      
         
         nrdiv = nrdiv + 1                                                              
      enddo ! while(.true.)
                                                           
      end  ! pqr                                                             

c--------------------------------------------------------
c avdv
c `average' dv?
c
c
      subroutine avdv(npart, ncup, ipr, pp, p, prt, n1, n2, a3)
      implicit real*8(a-h, o-z)                     
      logical prt                                 
      
      parameter (mxnx = 40, mxny = 40, maxbnd = 150)
      parameter (mxcr = 9)                        
      parameter (npmax = 2000)                   
      parameter (npmxh = npmax / 2)               
c     real*8 pp(mxcr,npmax),p(mxcr),sdev(mxcr)

      real*8 pp(mxcr, 1), p(mxcr), sdev(mxcr)  
      real*8 corr(mxcr), aver(20)
c     double precision wn,dev,vdev,sk,sc    
      common /lookuc/ t2p, maxr, maxi
      dimension ipr(npmxh)                 

      n1h = n1 / 2                            
      n2h = n2 / 2                           
      nctr = ncup                         
      vcu = n1 * n2 * a3

      if (nctr.gt.6) nctr = 6            

      if (npart.ge.100) then          
         write (6,70) npart          
 70      format('0npart = ',i8,' : partial results suppressed'//)         
      else                                 
         write (6,60) npart          
 60      format('0npart = ',i8,' partial results for averages'//)            
         write (6,61)(i,i = 1,nctr)   
 61      format(' sum terms',9(i8,4x))             
         write (6,62)                             
 62      format(//' ')                           
         
         do k = 1,npart                        
            write (6,64) k,ipr(k),(pp(j,k),j = 1,nctr)
 64         format(1h ,i3,i8,10f12.6)                 
         enddo ! k
         if (npart.le.3) return     

      endif

      do j = 1, ncup                         

         np = 0                                  
         sk = 0.0d0                             
         sc = 0.0d0                           
         dev = 0.0d0                         

         do k = 1, npart                              
            if (ipr(k).eq.0) goto 12  !break out of this do-loop (MM)
            vdev = dev                       
            dev = pp(j,k)-p(j)              
            sk = sk + dev*dev                
            sc = sc + dev*vdev              
            np = np + 1                    
         enddo ! k
 12      continue ! break out

         if (np.lt.3) return               
         sdev(j) = dsqrt(sk/(np*(np-1)))             

         if (sk.ne.0.0d0) then
            corr(j) = np * sc /(sk * (np - 1))  
         else
            corr(j) = 0.0d0
         endif

      enddo ! j

      write(6, 82)                            
 82   format(/'  nc     average    st. dev.    corr. c.')

      do j = 1, ncup                        
         write(6,80) (j), p(j), sdev(j), corr(j)                     
         if (prt) write(8,80) (j), p(j), sdev(j), corr(j)  
 80      format(1h ,i3, 3f14.8, i3,3f12.6, 3x, 3f12.6)                         
      enddo ! j

      b3 = vcu
      b4 = b3 * vcu
      aver(1) = (p(1)**2) / p(2) 
      aver(2) = vcu * (p(4) - p(3)**2) ! c
      aver(3) = b3 * (p(5) - 3.0d0 * p(4) * p(3) + 2.0d0 * p(3)**3) ! c'
      aver(4) = b4 * (p(6) - 4.0d0 * p(5) * p(3) + 
     +                12.0d0 * p(4) * p(3)**2 - 3.0d0 * p(4)**2 -
     -                6.0d0 * p(3)**4 )       ! c''
      aver(5) = p(7) - p(3) * p(1)  ! em2
      aver(6) = p(8) - p(3) * p(2)  ! em4
      aver(7) = 2.0d0 * p(7) / p(1) - p(8) / p(2) - p(3) ! qp
      
      do i = 1, 7
         j = ncup + i
         np = 0                   
         sk = 0.0d0              
         sc = 0.0d0             
         dev = 0.0d0           
         do k = 1, npart                  
            vdev = dev                   
    
            goto (30, 31, 32, 33, 34, 35, 36) i  ! switch(i)
c i=1
 30         if (pp(1,k).eq.0.0d0) then
               dev = 1.0d0/3.0d0-aver(1) 
            else
               dev = (pp(1, k)**2)/pp(2, k) - aver(1) 
            endif     
            goto 50
c i=2
 31         dev = vcu * (pp(4, k) - pp(3, k)**2) - aver(2)
            goto 50
c i=3
 32         dev = b3 * ( pp(5, k) - 3.0d0 * pp(4, k) * pp(3, k)  + 
     +                   2.0d0 * pp(3, k)**3 ) -
     -            aver(3)
            goto 50
c i=4
 33         dev = b4 * ( pp(6, k) - 4.0d0 * pp(5, k) * pp(3, k)  + 
     +                   1.2d1 * pp(4, k) * pp(3, k)**2 -
     -                   3.0d0 * pp(4, k)**2 - 6.0d0 * pp(3, k)**4 )-
     -            aver(4)
            goto 50
c i=5
 34         dev = (pp(7, k)-pp(3, k)*pp(1, k))-aver(5)
            goto 50
c i=6
 35         dev = (pp(8, k)-pp(3,k)*pp(2, k))-aver(6)
            goto 50
c i=7
 36         if (pp(1,k).eq.0.0d0) then
               dev = -aver(7) 
            else
               dev = 2*pp(7, k)/pp(1,k)-pp(8,k)/pp(2, k)-pp(3,k)-aver(7)  
            endif

 50         continue
            sk = sk + dev *  dev                 
            sc = sc + dev * vdev    !vdev? `vorige dev'?
            np = np + 1                     
         enddo ! do k = 1, npart

         sdevj1 = dsqrt(sk / (np * (np - 1)))              
         if (sk.ne.0.0d0) then
            corrj1 = np * sc / ( sk * (np - 1))                 
         else
            corrj1 = 0.0d0
         endif
         write(6,80) (j), aver(i), sdevj1, corrj1       
         if (prt) write(8,80) (j), aver(i), sdevj1, corrj1  
      enddo ! do i = 1, 7

c 20   continue  ! where is this one for ? 

c     write(6,81) aver(1),aver(2),aver(3),aver(4)
c81   format(' q = ',f10.6,' c = ',f10.5,' em2 = ',f12.7,' em4 = ',f12.7) 
      return                                                            
      end  ! avdv 

c---------------------------------------------------------------
c ransi
c random numbers 
      subroutine ransi(iseed)                                 

c sequential version                                           
c     shift register random generator with very long period     

      implicit real*8 (a-h, o-z)            
      save

      parameter (mult = 32781)               
      parameter (mod2 = 2796203, mul2 = 125)
      parameter (two  = 2d0,  tm32 = two**(-32))   
      parameter (len1 = 9689, ifd1 = 471)         
      parameter (len2 = 127,  ifd2 = 30)         
      integer*2 inxt1(len1)
      integer*2 inxt2(len2)
      common/ransrb/ ir1(len1), ir2(len2), inxt1, inxt2
      parameter (mxnx = 40, mxny = 40, maxbnd = 150)                     
      parameter (mxxy = mxnx * mxny)         
      parameter (mxrn = mxxy * maxbnd * 2)

      common/nran/ n, irn(mxrn) ! n random numbers are stored in array irn. (MM)

c     look out with iseed!!
c     may be I should build some control here.

      k = 3**18 + 2 * iseed     
      k1 =  1313131 * iseed
      k1 =  k1 - (k1 / mod2) * mod2

      do i = 1, len1      
         k =   k  * mult               
         k1 =  k1 * mul2
         k1 =  k1 - (k1 / mod2) * mod2
         ir1(i) = k + k1 * 8193        
      enddo ! i
             
      do i = 1, len2      
         k = k  * mult               
         k1 =  k1 * mul2
         k1 =  k1 - (k1 / mod2) * mod2
         ir2(i) = k + k1 * 4099        
      enddo ! i
             
      do i = 1, len1           
         inxt1(i) = i + 1
      enddo ! i
     
      inxt1(len1) = 1
      ipnt1       = 1
      ipnf1       = ifd1 + 1

      do i = 1, len2           
         inxt2(i) = i + 1                
      enddo ! i

      inxt2(len2) = 1                
      ipnt2       = 1                     
      ipnf2       = ifd2 + 1                 

      return          ! so ransi itself does not give numbers
      
c -----------------entry ransr--------------------------------
      entry ransr()
      
c     calculate n random numbers  (n is in common block)
      do i = 1, n                     

         l     = ieor( ir1(ipnt1) , ir1(ipnf1) )   
         k     = ieor( ir2(ipnt2) , ir2(ipnf2) )   
         irn(i)= ieor( k          , l          )

         ir1(ipnt1) = l                    
         ipnt1 = inxt1(ipnt1)             
         ipnf1 = inxt1(ipnf1)            

         ir2(ipnt2) = k                    
         ipnt2 = inxt2(ipnt2)             
         ipnf2 = inxt2(ipnf2)            

      enddo ! i

      return                       
      end ! ransi/ransr

c-----------------------------------------------------------------
c ctijd
c
      function ctijd(a)
c     function ctijd returns time in seconds; a is a dummy argument
      real*4 etime, ta(2)
      save
      t = etime(ta)
c     ctijd = ta(1)
      ctijd = t
      return
      end   ! ctijd

c -----------END OF FILE ani.f --------------------------- 
