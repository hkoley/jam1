      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15
      character*80 fnin,fnout
      logical dump 
      data dump/.false./

c...Default Parameters
      mstc(1)=20001224  ! random seed. (Merry Christmas!)
      nev=100            ! total simulation event
      bmin= 0.0D0        ! minimum impact parameter
      bmax=-3.3D0        ! maximum impact parameter 350mb
      dt=100.0D0         ! collision time(fm/c)
      nstep=1            ! # of time step (=1: Cascade)
      frame='collider'   ! Calculation frame
      cwin='200gev'      ! incident energy
      proj='197Au'       ! projectile
      targ='197Au'       ! target
      fname(6)='0'
      fname(7)='0'
      fname(8)='0'

1000  continue
      fnin='jam.cfg'
      read(*,'(a)',err=999) fnin
c     if(fnin(1:1).eq.'-') goto 998
c     if(fnin(1:1).eq.'0') goto 998
      fname(1)=fnin
c
c     read(*,'(a)',err=999) fnout
c     if(fnout(1:1).ne.'0') then
c        dump=.true.
c        open(33,file=fnout,status='unknown')
c     endif

c....Initialize JAM.
      call jaminit(nev,bmin,bmax,dt,nstep,frame,proj,targ,cwin)
	nev=mstc(2)
	bmin=parc(3)
	bmax=parc(4)

c     read(*,'(a)',err=999) fnout
      fnout=fname(6)
      if(fnout(1:1).ne.'0') then
         dump=.true.
         open(33,file=fnout,status='unknown')
      endif

c...[ Event Loop Start.
      do iev=1,mstc(2)
        call jamevt(iev)
c......[ Phase Space Output
        if(dump) then
          write(33,801)iev,nv,nbary,nmeson,pard(2)
          do i=1,nv
            write(33,802) k(2,i),p(5,i),(p(j,i),j=1,3),(v(j,i),j=1,4)
          end do
        endif
 801  format('#',4(1x,i8),f6.2)
 802  format(1x,i8,8(1x,1pe11.4))
c......] Phase Space Output
      end do
c...] Event Loop Start.

c...Final output.
      call jamfin

c     goto 1000

 998  write(*,*) 'Good Bye!'
 999  continue
      end
c***********************************************************************

      subroutine jamanaus0(indd,nadd)

c...Analysize collision spectra.
      include 'jam1.inc'
      include 'jam2.inc'
      dimension indd(100)
      end
c***********************************************************************

      subroutine jamanaus(indd,nadd)

c...Analysize collision spectra.
      include 'jam1.inc'
      include 'jam2.inc'
      dimension indd(100)

c...Local Variable
      parameter(mab=5)
      dimension dns(mab)
      irecdec=0

      if(irecdec.ne.1) return

c...Useful vectors for analysis.
      ichanel=mste(1)
      icltyp=mste(2)

c...Decay of resonance or string fragmentation.
      if(icltyp.eq.-1) then
      kf1=kcp(2,1)   ! PDG code for Mother
      i1=mste(21)    ! line number of the ingoing particle 1
c
c     write(36,901)pard(1),mste(21),kcp(1,1),kcp(2,1),pcp(5,1),nadd
c    $  ,chaf(mste(22),(3-isign(1,kcp(2,1)))/2)
c     write(36,'(''r='',5(g10.3,1x),''v5='',g15.5)')
c    $ (rcp(j,1),j=1,5),vcp(5,1)
c901  format(/f8.3,'fm/c Decay i=',i6,' ks=',i4,' kf=',i7,
c    $  ' em=',g10.3,' nadd=',i4,1x,a8)
c
c     write(36,801) k(2,i1),nadd,p(5,i1),(p(j,i1),j=1,3),(r(j,i1),j=1,4)
c
	if(nadd.ge.2) then
c     kf1=kcp(1,1)   ! PDG code for colliding particle 1
c     kf2=kcp(1,2)   ! PDG code for colliding particle 2
c     kf3=kcp(2,1)   ! PDG code for outgoing  particle 1
c     kf4=kcp(2,2)   ! PDG code for outgoing  particle 2
c     i2=mste(23)    ! line number of the ingoing particle 2
c     i3=mste(25)    ! line number of outgoing particle 1
c     i4=mste(27)    ! line number of outgoing particle 2
c...
      x0=rcp(1,1)
      y0=rcp(2,1)
      z0=rcp(3,1)
      time=rcp(4,1)
      call jampdns(i1,mab,dns,x0,y0,z0,time)
      write(35,801) kf1,nadd,pcp(5,1),(pcp(j,1),j=1,3),(rcp(j,1),j=1,4)
     &              ,(dns(j),j=1,5)
801   format(i7,1x,i3,1x,f6.2,7(1x,f7.2),5(1x,1pe10.3))
c...
         endif
      endif

      end
c***********************************************************************

      subroutine jampdns(i0,mab,dns,x0,y0,z0,time)

c...Density and Temperature at the point of particle i0
      include 'jam1.inc'
      include 'jam2.inc'
       
      parameter(widg=2*1.0d0)                   ! Gaussian width
      dimension dns(mab)
c
      dimension gmetric(4,4)
      dimension cu(4),cub(4),tens(4,4)
      data gmetric/ -1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,1/

c-----------------------------------------------------------------------
c...dns(1)  Lorentz invariant scalar number density
c...dns(2)  Lorentz invariant baryon density
c...dns(3)  Lorentz invariant energy density
c...dns(4)  Lorentz invariant pressure
c...dns(5)  Lorentz invariant temperature
c-----------------------------------------------------------------------
      widcof=(1.d0/(3.14159d0*widg))**1.5d0
c     x0=r(1,i0)
c     y0=r(2,i0)
c     z0=r(3,i0)
c     time=r(4,i0)
c=======================================================================

c...Clear Array
        do im=1,4
          cu(im)=0.0d0
          cub(im)=0.0d0
        do in=1,4
          tens(im,in)=0.0d0
        end do
        end do

c---------------------------------------------------------------------
c...[ Loop over all particles
      do 100 i=1,nv
	if(i.eq.i0) goto 100
 
        k1=k(1,i)
        if(k1.gt.10) goto 100   ! dead particle
        if(p(5,i).le.1d-5) goto 100
c       if(r(5,i).gt.time) goto 100
        if(abs(k(7,i)).eq.1) goto 100   ! not yet interaction

        dt=time-r(4,i)
        if(dt.lt.0.0d0) goto 100
c
        x1=r(1,i)+dt*p(1,i)/p(4,i)-x0
        y1=r(2,i)+dt*p(2,i)/p(4,i)-y0
        z1=r(3,i)+dt*p(3,i)/p(4,i)-z0
c
        bar=k(9,i)/3.0d0
        if(r(5,i).gt.time) then
          iq=mod(abs(k(1,i))/10,10)
          if(iq.eq.3) iq=2
          bar=iq/3d0*isign(1,k(2,i))
          if(k(9,i).eq.0)bar=abs(bar)
        endif
c
        xtra=x1**2+y1**2+z1**2
     $          +((x1*p(1,i)+y1*p(2,i)+z1*p(3,i))/p(5,i))**2
        if(xtra/widg.gt.30.d0) goto 100
        gam=p(4,i)/p(5,i)
        den=widcof*gam*exp(-xtra/widg)
c
        do im=1,4
          cu(im)   =cu(im)    + p(im,i)/p(4,i)*den
          cub(im)  =cub(im)   + p(im,i)/p(4,i)*den*bar
        do in=1,4
          tnsmn=p(im,i)*p(in,i)/p(4,i)
          tens(im,in)  =tens(im,in)  + tnsmn*den
        end do
        end do
100   continue
c...] Loop over all particles

c...Output
c---------------------------------------------------------------------
      do i=1,mab
       dns(i)=0.0d0
      end do

c...[
c...Lorentz invariant Scalar Number Density
      dns(1)=0.0d0
      cc=cu(4)**2-(cu(1)**2+cu(2)**2+cu(3)**2)
      if(cc.le.0.0d0) goto 250
      dns(1)=sqrt(cc)

c...Lorentz invariant Baryon density
      dns(2)=0.0d0
      bnorm=cub(4)**2 -( cub(1)**2+cub(2)**2+cub(3)**2 )
      if(bnorm.gt.0.0d0) dns(2)=sqrt(bnorm)

c...Lorentz invariant pressure and energy density
      dns(3)=0.0d0
      dns(4)=0.0d0
      do i=1,4
      do j=1,4
       tmp=gmetric(i,j)*cu(i)*cu(j)/cc
       dns(3)=dns(3)+tens(i,j)*tmp                    ! energy density
       dns(4)=dns(4)-1.d0/3.d0*tens(i,j)*(gmetric(i,j)-tmp)     ! pressure
      end do
      end do

c...Lorentz invariant temperature
      dns(5)=0.0d0
      if(dns(1).gt.0.0d0) dns(5)=dns(4)/dns(1)
c...]

250   continue
      end
c***********************************************************************
