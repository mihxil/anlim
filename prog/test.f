c for wolff cluster 3d s=1/2 ising model anisotropic limit
c select neighbors via continuous intervals to be skipped
c nrannr,irn in common block 
c periodic boundary conditions                                          
c program anlim.f (based on w1n.f)
c core: wn in common; vz=zsbnd ...
      implicit real*8 (a-h,o-z)                                         
      parameter (nprlll=1)                                              
      dimension itask(nprlll)                                           
 102   continue                                                         
      print *,' n1,n2,a3,ncycle,ntoss,numint,nint,nr,itmax'             
      read *,n1,n2,a3,ncycle,ntoss,numint,nint,nr,itmax                 
      itmax=itmax*3600
      if (n1.eq.0) stop                                                 
      print *,' tpar,fac'                                               
      read *,tpar,fac                                                   
      open(8,file='anli.dat',access='append')                          
      do 20 iprlll=1,nprlll                                             
      call bain(n1,n2,a3,ncycle,ntoss,numint,nint,nr,itmax,
     ,tpar,fac)    
20    continue                                                          
      close(8)                          
      goto 102                                                          
      end                                                               
      subroutine bain(n1,n2,a3,ncycle,ntoss,numint,nint,nr,itmax,
     ,tpar,fac)
c     program for simulation of anisotropic limit of ising model 
      implicit real*8(a-h,o-z)                                          
      save
      character*80 ident                                                
      dimension idate(14)                                               
      parameter (npmax=2000)                                           
      parameter (mxnx=10,mxny=10,maxbnd=100)                     
      parameter (mxxy=mxnx*mxny)         
      parameter (mxrn=mxxy*maxbnd*2)
      common /lookuc/ t2p,maxr,maxi
      data ident /'anli v0'/                                            
      print *,ident                                                     
      xt0=ctijd(0.)                                                   
      tmax=itmax-10
      ntemp=0                                                           
      if(n1.gt.mxnx.or.n2.gt.mxny) then                   
         print *,' main: at least one lattice dimension too big'        
         stop                                                           
      endif                                                             
      mcstps=numint*nint                                                
      npart=numint/2                                                    
      if (npart.gt.npmax) npart=npmax                                   
      if (npart.le.10) npart=1                                          
      npp=mcstps/npart                                                  
      nip=npp/nint                                                      
      npp=nint*nip                                                      
      mcstps=npart*npp                                                  
      ncount=mcstps/nint                                                
c104   continue                                                         
      ntemp=ntemp+1                                                     
      call carlo(tpar,fac,n1,n2,a3,ntemp,nr)       
      vcu=n1*n2*a3                                                      
      write(6,900) n1,n2,a3                                             
900   format(32h program for mc simulation  on a,                       
     *i4,1h*,i4,1h*,f12.6,8h lattice)                                      
      print 901,tpar,fac                                              
901   format(' tpar= ',f11.6,
     ,' factor=',f11.6)
      print 902,ncycle,mcstps,ntoss                                     
902   format(1x,i3,' consecutive runs of ',i10,' Wolff steps',1x,
     *'after',i10,' times nint steps for equilibration'/)
      print 924,nr                                                      
924   format(' nr=',i8/)                                                
      print 903,nint                                                    
903   format(' data are taken at intervals of ',i3,' mcs/site'/)        
      do 6 itoss=1,ntoss
      call monte(nint)
6     continue
      do 7 nsigh=1,ncycle                                               
      write(8,120) ident,n1,numint,nint,tpar,nr
120   format(a8,'; n1,numint,nint,knn,s:     ',i3,i10,i3,f12.6,12x,i4)
      call setpqr                                                       
      xt3=ctijd(0.)                                                   
      xt1=xt3-xt0                                                       
      print 554,xt1                                                     
554   format(20x,'setup time (sec.)=',f10.6)                            
      tmc=0                                                             
      tan=0                                                             
      print 999                                                         
999   format(/70('*')/)                                             
      nc=0                                                              
      do 11 ip=1,npart                                                  
      inp=ip                                                            
      do 1 iip=1,nip                                                    
      nc=nc+1                                                           
      call monte(nint)                          
      xt2=ctijd(0.)                                                   
      tmc=tmc+xt2-xt3                                                   
      call addpqr(n1,n2,ncup)                                       
      xt3=ctijd(0.)                                                   
      tan=tan+xt3-xt2                                                   
      if ((tmax-xt3+xt0).lt.0) goto 1111                                    
1     continue                                                          
      call part(inp,nc,ncup)                                            
11    continue                                                          
      goto 1112                                                         
1111  ncount=nc                                                         
      call part(inp,nc,ncup)                                            
      mcstps=ncount*nint                                                
      write(6,1120)mcstps,xt3-xt0                                           
1120  format(27h ****** emergency stop,mcs=,i9, 3h\= ,f11.3,' sec.')    
1112  continue                                                          
      write(6,1130) xt3-xt0                                           
1130  format(' mc stop at t=',f12.2,' sec.')    
c     call lpsb                                                         
      tmc=(tmc*1.d+6/mcstps)/vcu                                       
      tan=(tan*1.d+6/ncount)/vcu                                        
      call pqr(ncup,ncount,npart,n1,n2,a3)                             
      if ((tmax-xt3+xt0).lt.0) stop                                         
7     continue                                                          
      print 999                                                         
      print 553,maxr,mxrn,maxi,maxbnd
553   format('max nr random numbers used:',i10,' allowed:',i10/
     /'max nr  cluster boundaries:',i10,' allowed:',i10)
      xt4=ctijd(0.)                                                   
      anal=xt4-xt3                                                      
      t4=xt4-xt0                                                        
      print 555,tmc,tan,anal,t4                                         
555   format('0mc simulation  time (microseconds/mc step/site)=',f20.6/ 
     *  /' summation over lattice (microseconds/time/site)=',f20.6//    
     *  ' analysis (seconds) =',f10.2// ' total time (sec.)=',f10.2)    
      return                                                            
      end                                                               
      subroutine carlo(tpar,fac,m1,m2,aa,ntmp,nr)
c     version for periodic ising model                                  
      implicit real*8(a-h,o-z)                                          
      save
      parameter (eps=2.d-15)                                            
      parameter (mxnx=10,mxny=10,maxbnd=100)                     
      parameter (mxxy=mxnx*mxny)                                   
      parameter (mxrn=mxxy*maxbnd*2)
      parameter (ilut=6)                                    
      common /lookup/ rmul,radd,rmul2,radd2,rmul3,radd3
      common /lookuc/ t2p,maxr,maxi
      common/nran/ nrannr,irn(mxrn)
      integer*2 iz0s(mxxy),nin(mxxy)
      common/mats/ zsbnd(mxxy,maxbnd),iz0s,nin
      integer*2 ngc(4,mxxy)
      common /ngcb/ ngc
      parameter (tp16=2**16,tm32=1.0d0/(tp16*tp16))
      parameter (tp32m1=tp16*tp16-1,maxint=2147483647)
      common/dats/a3,wn,n1,n2,n12
4     if (m1-mxnx) 8,8,7                                                
7     write(6,9) m1,m2                                                  
9     format(44h1carlo disagrees with specified lattice size,2i3)       
      stop                                                              
8     continue                                                          
      if (m2-mxny) 6,6,7                                                
6     continue                                                          
      a3=aa
      t2p=tpar/2.0d0
      maxr=0 
      maxi=0
      facc=fac
      n1=m1                                                             
      n2=m2                                                             
      n12=n1*n2                                                         
      wn=1.0d0/(a3*n12)
      istep=2                                                           
      itel=0                                                            
      do 210 index=1,n12                                                
c     write(6,202) index,itel,istep                          
c02   format(' index,itel,istep',30i3)                       
      itel=itel+1                                                       
      iy=(index-1)/n1+1                                      
      ix=index-(iy-1)*n1                                     
      ixp=ix+1-(ix/n1)*n1                                               
      iyp=iy+1-(iy/n2)*n2                                               
      ixm=ix-1+((n1-ix+1)/n1)*n1                                        
      iym=iy-1+((n2-iy+1)/n2)*n2                                        
      ngc( 1,index)=(iy -1)*n1+ixm
      ngc( 3,index)=(iy -1)*n1+ixp
      ngc( 2,index)=(iym-1)*n1+ix 
      ngc( 4,index)=(iyp-1)*n1+ix 
c     write(6,950)index,ngc(1,index),ngc(2,index),ngc(3,index),
c    ,ngc(4,index)
c950  format('carlo index,ngc',i4,2x,4i4)
210   continue                                                          
      radd=0.5d0*n12+1.0d0
      rmul=tm32*(1.0d0-1.0d-15)*n12
      radd2=0.5d0*a3
      rmul2=tm32*(1.0d0-1.0d-15)*a3
      radd3=0.5d0
      rmul3=tm32*(1.0d0-1.0d-15)
      write(6,*)' wolff mc of s=1/2 anisotropic ising model'   
      if(ntmp.gt.1) goto 20                                             
      iseed=iabs(nr)+1                                                  
      call ransi(iseed)                                  
      if(nr.ge.0) goto 320                                              
      do 112 ij=1,n12
      iz0s(ij)=1
112   nin(ij)=0
      goto 20                                                           
320   continue                                                          
      nrannr=n12
      call ransr()                         
      do 330 ij=1,n12                                                   
      rnd=irn(ij)*rmul3+radd3
      ipm=rnd+rnd
      ipm=-1+ipm+ipm                                                  
      iz0s(ij)=ipm                                              
      nin(ij)=0
330   continue                                                          
20    continue                                                          
      nrannr=n12*maxbnd*2
      return                                                            
      end                                                               
      subroutine monte(nsteps)                
c perform nsteps mc steps per site                                      
c     version for arbitrary model                                       
      implicit real*8 (a-h,o-z)                                         
      parameter (mxnx=10,mxny=10,maxbnd=100)                     
      parameter (mxxy=mxnx*mxny)                                   
c     do 2 nn=1,nsteps                                                  
      call mcw(nsteps)                           
c2    continue                                                          
      return                                                            
      end                                                               
      subroutine mcw(nsteps)                                     
c     anisotropic wolff cluster
      implicit real*8 (a-h,o-z)                                         
      parameter (mxnx=10,mxny=10,maxbnd=100)                     
      parameter (mxxy=mxnx*mxny)                                   
      parameter (mxrn=mxxy*maxbnd*2)
      common/nran/ nrannr,irn(mxrn)
      integer*2 iz0s(mxxy),nin(mxxy)
      common/mats/ zsbnd(mxxy,maxbnd),iz0s,nin
      common/dats/a3,wn,n1,n2,n12
      common /lookuc/ t2p,maxr,maxi
      integer*2 ngc(4,mxxy)
      common /ngcb/ ngc
      common /lookup/ rmul,radd,rmul2,radd2,rmul3,radd3
      parameter (tp16=2**16,tm32=1.0d0/(tp16*tp16))
      integer*2 ijstack(mxxy*maxbnd)
      real*8     zstack(mxxy*maxbnd)
      real*8     wstack(mxxy*maxbnd)
c     logical prt
      save
      do 150 isteps=1,nsteps
c     prt=(isteps.eq.17)
      call ransr()
c choose origin
      nrannr=1
      ik=irn(nrannr)*rmul+radd
      nrannr=2
      zk=irn(nrannr)*rmul2+radd2
c determine sign of cluster (before flip)
      icsp=iz0s(ik)
      jin=0         
 180  continue
      if (jin.ge.nin(ik)) goto 181
      jp1=jin+1
      if (zsbnd(ik,jp1).gt.zk) goto 181
      jin=jp1
      icsp=-icsp
      goto 180
 181  continue
      nstack=0
      cw1=0.0d0
      nrannr=3
      brd=rmul3*irn(nrannr)+0.5d0
      brd=-dlog(brd)
      nrannr=4
      bod=rmul3*irn(nrannr)+0.5d0
      bod=-dlog(bod)*t2p
      cw4=4.0d0*cw1
      abd=0.0d0
c  zk,cw4,ik,abd,icsp,brd,bod defined
 104  continue  
c begin determination horiz. bounds of cluster
      if (nin(ik).gt.0) goto 194
      if (brd.ge.a3) then
         cl=0.0d0
         cw=a3
         brd=brd-a3
         iz0s(ik)=-iz0s(ik)
         goto 193
      endif
      dlc=brd
      dr=a3-brd
      nrannr=nrannr+1
      brd=rmul3*irn(nrannr)+0.5d0
      brd=-dlog(brd)
      write(6,889) irn(nrannr),brd  
      if (brd.ge.dr) then
         cl=0.0d0
         cw=a3
         brd=brd-dr
         iz0s(ik)=-iz0s(ik)
         goto 193
      endif
      nin(ik)=2   
      cl=zk-dlc
      if (cl.lt.0.0d0) then
         cl=cl+a3
         iz0s(ik)=-iz0s(ik)
         zsbnd(ik,1)=zk+brd
         zsbnd(ik,2)=cl
      else
         cr=zk+brd
         if (cr.ge.a3) then
           iz0s(ik)=-iz0s(ik)
           zsbnd(ik,1)=cr-a3
           zsbnd(ik,2)=cl
         else
           zsbnd(ik,1)=cl
           zsbnd(ik,2)=cr
         endif
      endif
      cw=dlc+brd
c calculate new break length
      nrannr=nrannr+1
      brd=rmul3*irn(nrannr)+0.5d0
      brd=-dlog(brd)
      write(6,889) irn(nrannr),brd
      goto 193
  194 continue
c  determine size Ising cluster
      if (jin.eq.0) then 
         pl=zsbnd(ik,nin(ik))
         dl=zk-pl+a3   
      else
         pl=zsbnd(ik,jin)
         dl=zk-pl
      endif
      if (jin.eq.nin(ik)) then 
         pr=zsbnd(ik,1)
         dr=pr+a3-zk
      else
         pr=zsbnd(ik,jin+1)
         dr=pr-zk
      endif
c first find left random cluster reach
      if (brd.ge.dl) then
         cl=pl
         brd=brd-dl
c  eliminate left interface
         if (jin.ne.0)  then
           do 183 ii=jin+1,nin(ik)
           zsbnd(ik,ii-1)=zsbnd(ik,ii)
 183       continue
           jin=jin-1
         else
           iz0s(ik)=-iz0s(ik)
         endif
         nin(ik)=nin(ik)-1
      else
         dl=brd
c  add left interface
         cl=zk-brd
         nin(ik)=nin(ik)+1
         if (cl.lt.0.0d0) then
           cl=cl+a3
           zsbnd(ik,nin(ik))=cl
           iz0s(ik)=-iz0s(ik)
         else
           jin=jin+1
           do 184 ii=nin(ik),jin+1,-1
           zsbnd(ik,ii)=zsbnd(ik,ii-1)
  184      continue
           zsbnd(ik,jin)=cl
         endif
c     write(6,351) ik,nin(ik),(zsbnd(ik,iii),iii=1,nin(ik))
c351  format('ad l int: k,nin,zsbnd',2i3,20f6.2)
         nrannr=nrannr+1
         brd=rmul3*irn(nrannr)+0.5d0
         brd=-dlog(brd)
      write(6,889) irn(nrannr),brd
      endif
c  analogous for right hand interval
      if (brd.gt.dr) then
         cr=pr
         brd=brd-dr
         cw=dl+dr
c  eliminate right interface
         if (jin.ne.nin(ik))  then
           do 185 ii=jin+1,nin(ik)-1
           zsbnd(ik,ii)=zsbnd(ik,ii+1)
 185       continue
         else
           do 195 ii=1,nin(ik)-1
           zsbnd(ik,ii)=zsbnd(ik,ii+1)
 195       continue
           iz0s(ik)=-iz0s(ik)
         endif
         nin(ik)=nin(ik)-1
c     write(6,352) ik,nin(ik),(zsbnd(ik,iii),iii=1,nin(ik))
c352  format('el r int: k,nin,zsbnd',2i3,20f6.2)
      else
c  add right interface
         dr=brd
         cw=dl+dr 
         cr=zk+brd
         nin(ik)=nin(ik)+1
         if (cr.gt.a3) then
           cr=cr-a3
           do 187 ii=nin(ik),2,-1
           zsbnd(ik,ii)=zsbnd(ik,ii-1)
 187       continue
           zsbnd(ik,1)=cr
           iz0s(ik)=-iz0s(ik)
         else
           do 186 ii=nin(ik),jin+2,-1
           zsbnd(ik,ii)=zsbnd(ik,ii-1)
 186       continue
           zsbnd(ik,jin+1)=cr
         endif
         nrannr=nrannr+1
         brd=rmul3*irn(nrannr)+0.5d0
         brd=-dlog(brd)
      write(6,889) irn(nrannr),brd  
 889  format('breakd.:'i12,2f20.14)
      endif
c  end determination horiz. bounds of flipped part
c  write in stack
 193  nstack=nstack+1
      ijstack(nstack)=ik
      zstack(nstack)=cl
      wstack(nstack)=cw
      write(6,888) nstack,ik,cl,cw
 888  format('stackwr:'2i6,2f20.14)
 192  abd=abd+bod
      if (abd.ge.cw4) then
         bod=abd-cw4
         goto 161
      endif
      idr=abd/cw1
      zk=cl1+abd-idr*cw1
      if (zk.gt.a3) zk=zk-a3
      ik=ngc(idr+1,ij)
c find sign of nb
      isnb=iz0s(ik)
      jin=0         
 188  continue
      if (jin.ge.nin(ik)) goto 189
      jp1=jin+1
      if (zsbnd(ik,jp1).gt.zk) goto 189
      jin=jp1
      isnb=-isnb
      goto 188
 189  continue
c end find sign
c first determine interval to next perp bond
      nrannr=nrannr+1
      bod=rmul3*irn(nrannr)+0.5d0
      bod=-dlog(bod)*t2p
      if (isnb.eq.icsp) goto 104
      goto 192
 161  continue
      if (nstack.eq.0) goto 105
c read from stack
      ij=ijstack(nstack)
      cl1=zstack(nstack)
      cw1=wstack(nstack)
      cw4=4.0d0*cw1
      nstack=nstack-1
c     write(6,981) isteps,nstack,ij,cl1,cw1
c981  format ('str isteps,nstack ij cl cw',3i4,2f8.3)
      abd=0.0d0
      goto 192
 105  continue     
      if (nrannr.gt.maxr) then
         maxr=nrannr
         if (nrannr.gt.mxrn) then
           print *,' **** nrannr>mxrn:',nrannr,mxrn
           stop
         endif
      endif
 150  continue     
      return          
      end
c  zk intoduced everywhere; OK for both par and perp bonds?
c parallel parms: append 1
c vermijd interf parms in par and perp sectie
c interfaces (grain boundaries) everywhere updated?
c check precise bounds of intervals generated by bod and brd
      subroutine core(sn,ncor)                                
c     version for cylinder ising model                                  
      implicit real*8 (a-h,o-z)                                         
      parameter (mxnx=10,mxny=10,maxbnd=100)                     
      parameter (mxxy=mxnx*mxny)                                   
      parameter (mxrn=mxxy*maxbnd*2)
      common/nran/ nrannr,irn(mxrn)
      integer*2 iz0s(mxxy),nin(mxxy)
      common/mats/ zsbnd(mxxy,maxbnd),iz0s,nin
      common/dats/a3,wn,n1,n2,n12
      common /lookuc/ t2p,maxr,maxi
      integer*2 ngc(4,mxxy)
      common /ngcb/ ngc
      parameter (mxcr=9)                                      
      dimension  sn(mxcr)
      save
      ncor=9                                             
c     if (ncor.eq.9) goto 10
c     sample correlation function over half the system size 
 10   continue
c     sample magnetization and nr of interfaces
      nit=0
      amt=0
      do 18 ij=1,n12                                                    
      nij=nin(ij)
      if (nij.gt.maxi) maxi=nij
      nit=nit+nij
      izs=iz0s(ij)
      amt=amt+a3*izs
      if (nij.gt.0) vz=zsbnd(ij,nij)
      do 19 int=1,nij
        z=zsbnd(ij,int)
        amt=amt+(z-vz)*izs
        izs=-izs
        vz=z
 19   continue
 18   continue                                                          
      sps=0.0d0
      do 26 ij1=1,n12
      ni1=nin(ij1)
      do 26 inb=1,2
      in1=0
      in2=0
      ij2=ngc(inb,ij1)
      ni2=nin(ij2)
      isp2=iz0s(ij1)*iz0s(ij2)
      sps=sps+a3*isp2
      isp2=isp2+isp2
      if (0.eq.ni1) goto 21
      if (0.eq.ni2) goto 22
 20   continue
      if (zsbnd(ij1,in1+1).gt.zsbnd(ij2,in2+1)) then
         in2=in2+1
         sps=sps+isp2*zsbnd(ij2,in2)
         isp2=-isp2
         if (in2.eq.ni2) goto 22
      else
         in1=in1+1
         sps=sps+isp2*zsbnd(ij1,in1)
         isp2=-isp2
         if (in1.eq.ni1) goto 21
      endif
      goto 20
 21   continue
      do 23 jn2=in2+1,ni2
      sps=sps+isp2*zsbnd(ij2,jn2)
      isp2=-isp2
 23   continue
      goto 25
 22   continue
      do 24 jn1=in1+1,ni1
      sps=sps+isp2*zsbnd(ij1,jn1)
      isp2=-isp2
 24   continue
 25   continue
 26   continue
      bsn=nit*wn
      sps=sps*wn
      st=amt*wn
      s2=st*st                                                         
      sn(1)=s2
      sn(2)=s2*s2
      sn(3)=bsn
      sn(4)=bsn*bsn
      sn(5)=sn(3)*sn(4)
      sn(6)=sn(4)*sn(4)
      sn(7)=sn(1)*sn(3)
      sn(8)=sn(2)*sn(3)
      sn(9)=sps
      return                                                            
      end                                                               
      subroutine setpqr                                                 
      implicit real*8 (a-h,o-z)                                         
      parameter (npmax=2000)                                           
      parameter (mxnx=10,mxny=10,maxbnd=100)                     
      parameter (mxcr=9)                                      
      double precision p(mxcr),q(mxcr,mxcr),pp(mxcr,npmax)              
      dimension ipr(npmax)                                              
      common/pqrh/p,q,pp,ipr                                            
      do 1 k=1,mxcr                                                     
      p(k)=0.0                                                          
      do 1 l=1,mxcr                                                     
      q(k,l)=0.0                                                        
1     continue                                                          
      do 4 k=1,npmax                                                    
4     ipr(k)=0                                                          
      return                                                            
      end                                                               
      subroutine addpqr(n1,n2,ntop)                                 
      implicit real*8(a-h,o-z)                                          
      parameter (npmax=2000)                                           
      parameter (mxnx=10,mxny=10,maxbnd=100)                     
      parameter (mxcr=9)                                      
      dimension ipr(npmax)                                              
      double precision sn(mxcr)                                         
      double precision p(mxcr),q(mxcr,mxcr),pp(mxcr,npmax)              
      common/pqrh/p,q,pp,ipr                                            
      call core(sn,ntop)                                      
      do 1 j=1,ntop                                                     
      p(j)=p(j)+sn(j)                                                   
      do 1 k=1,ntop                                                     
      q(j,k)=q(j,k)+sn(j)*sn(k)                                         
1     continue                                                          
      return                                                            
      end                                                               
      subroutine part(ip,nc,ntop)                                       
      implicit real*8(a-h,o-z)                                          
      parameter (npmax=2000)                                           
      parameter (mxnx=10,mxny=10,maxbnd=100)                     
      parameter (mxcr=9)                                      
      dimension ipr(npmax)                                              
      double precision p(mxcr),q(mxcr,mxcr),pp(mxcr,npmax)              
      common/pqrh/p,q,pp,ipr                                            
c in part:
      save
      do 100 i=1,ntop                                                   
      pp(i,ip)=p(i)                                                     
 100  p(i)=0.0d0                                                     
      if (ip.eq.1) nv=0
      ipr(ip)=nc-nv                                                     
      nv=nc
c end
      return                                                            
      end                                                               
      subroutine pqr(ncup,ncount,npart,n1,n2,a3)                       
      implicit real*8(a-h,o-z)                                          
      logical prt                                                       
      parameter(npmax=2000)                                            
      parameter(npmxh=npmax/2)                                          
      parameter (mxnx=10,mxny=10,maxbnd=100)                     
      parameter (mxcr=9)                                      
      double precision p(mxcr),q(mxcr,mxcr),pp(mxcr,npmax)              
      double precision b(mxcr,mxcr),rp(mxcr,npmxh)                      
      dimension ipr(npmax),jpr(npmxh),ndiv(12)                          
      common/pqrh/p,q,pp,ipr                                            
      common/cbh/b                                                      
      data ndiv /2,5,10,20,50,100,200,500,1000,2000,5000,10000/         
      prt=.false.                                                       
      do 5 j=1,ncup                                                     
      p(j)=0.0d0
      do 11 i=1,npart                                                  
      p(j)=p(j)+pp(j,i)                                                
11    continue                                                          
      p(j)=p(j)/ncount                                                  
5     continue                                                          
      do 1 j=1,ncup                                                     
      do 1 k=1,ncup                                                     
      q(j,k)=q(j,k)/ncount                                              
      b(j,k)=q(j,k)-p(j)*p(k)                                           
1     continue                                                          
      do 20 i=1,npart                                                   
      if (ipr(i).eq.0) goto 20                                          
      do 21 j=1,ncup                                                    
21    pp(j,i)=pp(j,i)/ipr(i)                                            
20    continue                                                          
      call avdv(npart,ncup,ipr,pp,p,prt,n1,n2,a3)                      
      nrdiv=1                                                           
30    npdiv=npart/ndiv(nrdiv)                                           
      prt=(nrdiv.eq.1)                                         
      if (npdiv.le.3) return                                            
      do 31 i=1,npdiv                                                   
      jpr(i)=0                                                          
      do 31 j=1,ncup                                                    
      rp(j,i)=0.0d0                                                     
31    continue                                                          
      npd=0                                                             
      ind=0                                                             
      do 32 i=1,npdiv                                                   
      ntog=ndiv(nrdiv)                                                  
      do 33 k=1,ntog                                                    
      ind=ind+1                                                         
      jpr(i)=jpr(i)+ipr(ind)                                            
      do 34 j=1,ncup                                                    
34    rp(j,i)=rp(j,i)+pp(j,ind)*ipr(ind)                                
33    continue                                                          
      if (jpr(i).eq.0) goto 32                                          
      npd=npd+1                                                         
      do 35 j=1,ncup                                                    
35    rp(j,i)=rp(j,i)/jpr(i)                                            
32    continue                                                          
      call avdv(npdiv,ncup,jpr,rp,p,prt,n1,n2,a3)                      
      nrdiv=nrdiv+1                                                     
      goto 30                                                           
      end                                                               
      subroutine avdv(npart,ncup,ipr,pp,p,prt,n1,n2,a3)
      implicit real*8(a-h,o-z)                     
      logical prt                                 
      parameter (mxnx=10,mxny=10,maxbnd=100)
      parameter (mxcr=9)                        
      parameter (npmax=2000)                   
      parameter (npmxh=npmax/2)               
c     real*8 pp(mxcr,npmax),p(mxcr),sdev(mxcr)
      real*8 pp(mxcr,1),p(mxcr),sdev(mxcr)  
      real*8 corr(mxcr),aver(20)
c     double precision wn,dev,vdev,sk,sc    
      common /lookuc/ t2p,maxr,maxi
      dimension ipr(npmxh)                 
      n1h=n1/2                            
      n2h=n2/2                           
      nctr=ncup                         
      vcu=n1*n2*a3
      if (nctr.gt.6) nctr=6            
      if (npart.ge.100) then          
         write (6,70) npart          
 70      format('0npart=',i8,' : partial results suppressed'//)         
         goto 71                    
      endif                        
      write (6,60) npart          
 60   format('0npart=',i8,' partial results for averages'//)            
      write (6,61)(i,i=1,nctr)   
 61   format(' sum terms',9(i8,4x))             
      write (6,62)                             
 62   format(//' ')                           
      do 63 k=1,npart                        
      write (6,64) k,ipr(k),(pp(j,k),j=1,nctr)
 64   format(1h ,i3,i8,10f12.6)                 
 63   continue                                 
      if (npart.le.3) return                  
 71   do 10 j=1,ncup                         
      np=0                                  
      sk=0.0d0                             
      sc=0.0d0                           
      dev=0.0d0                         
      do 11 k=1,npart                  
      if (ipr(k).eq.0) goto 12        
      vdev=dev                       
      dev=pp(j,k)-p(j)              
      sk=sk+dev*dev                
      sc=sc+dev*vdev              
      np=np+1                    
 11   continue                  
 12   if (np.lt.3) return      
      sdev(j)=dsqrt(sk/(np*(np-1)))             
      if (sk.ne.0.0d0) then
      corr(j)=np*sc/(sk*(np-1))  
      else
      corr(j)=0.0d0
      endif
 10   continue                  
      write(6,82)                            
 82   format(/'  nc     average    st. dev.    corr. c.')
      do 15 j=1,ncup                        
      write(6,80) (j),p(j),sdev(j),corr(j)                     
      if (prt) write(8,80) (j),p(j),sdev(j),corr(j)  
 80   format(1h ,i3,3f14.8,i3,3f12.6,3x,3f12.6)                         
 15   continue                                 
      b3=vcu
      b4=b3*vcu
      aver(1)=(p(1)**2)/p(2)                             ! q
      aver(2)=vcu*(p(4)-p(3)**2)                        ! c
      aver(3)=b3*(p(5)-3.0d0*p(4)*p(3)+2.0d0*p(3)**3)    ! c'
      aver(4)=b4*(p(6)-4.0d0*p(5)*p(3)+
     +12.0d0*p(4)*p(3)**2-3.0d0*p(4)**2-6.0d0*p(3)**4)   ! c''
      aver(5)=p(7)-p(3)*p(1)                             ! em2
      aver(6)=p(8)-p(3)*p(2)                             ! em4
      aver(7)=2.0d0*p(7)/p(1)-p(8)/p(2)-p(3)             ! qp
      do 22 i=1,7
      j=ncup+i
      np=0                   
      sk=0.0d0              
      sc=0.0d0             
      dev=0.0d0           
      do 21 k=1,npart                  
      vdev=dev                       
      goto (30,31,32,33,34,35,36) i
 30   if (pp(1,k).eq.0.0d0) then
         dev=1.0d0/3.0d0-aver(1) 
      else
         dev=(pp(1,k)**2)/pp(2,k)-aver(1) 
      endif
      goto 50
 31   dev=vcu*(pp(4,k)-pp(3,k)**2)-aver(2)
      goto 50
 32   dev=b3*(pp(5,k)-3.0d0*pp(4,k)*pp(3,k)+2.0d0*pp(3,k)**3)-aver(3)
      goto 50
 33   dev=b4*(pp(6,k)-4.0d0*pp(5,k)*pp(3,k)+1.2d1*pp(4,k)*pp(3,k)**2-
     -3.0d0*pp(4,k)**2-6.0d0*pp(3,k)**4)-aver(4)
      goto 50
 34   dev=(pp(7,k)-pp(3,k)*pp(1,k))-aver(5)
      goto 50
 35   dev=(pp(8,k)-pp(3,k)*pp(2,k))-aver(6)
      goto 50
 36   if (pp(1,k).eq.0.0d0) then
         dev=-aver(7) 
      else
         dev=2*pp(7,k)/pp(1,k)-pp(8,k)/pp(2,k)-pp(3,k)-aver(7)  
      endif
 50   continue
      sk=sk+dev*dev                 
      sc=sc+dev*vdev               
      np=np+1                     
 21   continue                                   
      sdevj1=dsqrt(sk/(np*(np-1)))              
      if (sk.ne.0.0d0) then
      corrj1=np*sc/(sk*(np-1))                 
      else
      corrj1=0.0d0
      endif
      write(6,80) (j),aver(i),sdevj1,corrj1       
      if (prt) write(8,80) (j),aver(i),sdevj1,corrj1  
 22   continue
 20   continue
c     write(6,81) aver(1),aver(2),aver(3),aver(4)
c81   format(' q=',f10.6,' c=',f10.5,' em2=',f12.7,' em4=',f12.7) 
      return                                                            
      end                                                               
      subroutine ransi(iseed)                                 
c sequential version                                           
c     shift register random generator with very long period     
      implicit real*8(a-h,o-z)            
      save
      parameter (mult=32781)               
      parameter (mod2=2796203,mul2=125)
      parameter (two=2d0,tm32=two**(-32))   
      parameter (len1=9689,ifd1=471)         
      parameter (len2=127,ifd2=30)         
      integer*2 inxt1(len1)
      integer*2 inxt2(len2)
      common/ransrb/ ir1(len1),ir2(len2),inxt1,inxt2
      parameter (mxnx=10,mxny=10,maxbnd=100)                     
      parameter (mxxy=mxnx*mxny)         
      parameter (mxrn=mxxy*maxbnd*2)
      common/nran/ n,irn(mxrn)
      k=3**18+2*iseed     
      k1=1313131*iseed
      k1=k1-(k1/mod2)*mod2
      do 100 i=1,len1      
      k=k*mult               
      k1=k1*mul2
      k1=k1-(k1/mod2)*mod2
      ir1(i)=k+k1*8193        
 100  continue                 
      do 102 i=1,len2      
      k=k*mult               
      k1=k1*mul2
      k1=k1-(k1/mod2)*mod2
      ir2(i)=k+k1*4099        
 102  continue                 
      do 101 i=1,len1           
 101  inxt1(i)=i+1                
      inxt1(len1)=1                
      ipnt1=1                     
      ipnf1=ifd1+1                 
      do 103 i=1,len2           
 103  inxt2(i)=i+1                
      inxt2(len2)=1                
      ipnt2=1                     
      ipnf2=ifd2+1                 
      return                         
      entry ransr()
c     calculate n random numbers      
      do 200 i=1,n                     
      l=ieor(ir1(ipnt1),ir1(ipnf1))   
      k=ieor(ir2(ipnt2),ir2(ipnf2))   
      irn(i)=ieor(k,l)
      ir1(ipnt1)=l                    
      ipnt1=inxt1(ipnt1)             
      ipnf1=inxt1(ipnf1)            
      ir2(ipnt2)=k                    
      ipnt2=inxt2(ipnt2)             
      ipnf2=inxt2(ipnf2)            
 200  continue                      
      return                       
      end                         
      function ctijd(a)
c     function ctijd returns time in seconds; a is a dummy argument
      implicit real*8 (a-h,o-z)
      real*4 etime,ta(2)
      save
      t = etime(ta)
c     ctijd = ta(1)
      ctijd = t
      return
      end

