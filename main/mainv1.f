c...A main program for calculation of v1(y) and v2(y).

      include 'jam1.inc'
      include 'jam2.inc'
      external pjdat

      character frame*8,proj*8,targ*8,cwin*15
      character feos*80
      logical dump 
c     data bmin/0.0/,bmax/-1.08/,dt/100.0/,mevent/1/
c     data dump/.false./
      data dump/.true./


c      print *,'paru(1)',paru(1),mstc(2)

      mstc(52)=2
c     parc(5)=0d0

      if(dump)
     $  open(33,file='phase.dat',status='new')

c     mdcy(jamcomp(321),1)=1    ! K+
c     mdcy(jamcomp(-321),1)=1   ! K-

c....Initialize JAM.
c     mstc(6)=101         ! RQMD/S mode
c     mstc(108)=1         ! feel potential after scatt.
c     mstc(109)=2    ! 1:include potential in p^0  2:use effective mass
c     mstc(7)=0      ! no collision
c     mstc(104)=12 
c     mstc(79)=1   ! no s-channel string
c     mstc(58)=1

      fname(1)='../jam.cfg'  ! input file name.
c     fname(1)='jam.cfg'  ! input file name.
      mstc(1)=0   ! random seed.
      mevent=1 ! total simulation event

c     bmin=4.0D0          ! minimum impact parameter
c     bmax=-8.0D0        ! maximum impact parameter

c...Star
c     bmin=4.6D0          ! minimum impact parameter
c     bmax=-9.4D0        ! maximum impact parameter

      bmin=0.0
c     bmax=-3.4
      bmax=-4.0

      dt=100.0D0          ! collision time(fm/c)
      nstep=1

c     dt=0.1D0          ! collision time(fm/c)
c     nstep=400

      if(mstc(6).eq.101) then
        dt=0.1d0            ! time step size (fm/c)
        nstep=400           ! total number of timestep
c       mstc(106)=14        !=1 mh, =2 mm, =3 ms, =4 ah, =5 as  par. set
c       mstc(107)=0         !=0 mom+sky(mh,mm,ms), =1 sky(ah,as), =2 cut mom+sky
c       mstc(103)=0        ! RQMD/S t-development approx.
      endif
      cwin='2gev'
c     cwin='4gev'
c     cwin='6gev'
c     cwin='8gev'
c     cwin='8gev'
c     cwin='10.7gev         '  ! incident energy
c     cwin='20gev         '  ! incident energy
c     cwin='30gev         '  ! incident energy
      frame='nn      '        ! comp. frame
c     cwin='11.7gev         '  ! incident energy
c     cwin='158gev         '  ! incident energy
c     cwin='7.7gev         '  ! incident energy
c     cwin='11.5gev         '  ! incident energy
c     cwin='19.6gev         '  ! incident energy
c     frame='collider'        ! comp. frame
      proj='197Au   '         ! projectile
      targ='197Au   '         ! projectile
c     proj='208Pb   '         ! projectile
c     targ='208Pb   '         ! target
c     proj='32S'
c     targ='32S'
c     mstc(8)=3  ! job mode.
c     mstc(156)=1  ! analysis of collision distribution
c     mstc(155)=0  ! flow anal.
c     mstc(163)=0  ! time evolution of directed transverse flow
c     parc(7)= 1.0D0    ! Output time interval (fm/c)
c     mstc(51)=0  !   only BB collisions
c     mstc(55)=1  !   frozen resonance
c     mstc(63)=0  ! frozen delta
c     mstc(111)=1 !  econ


      call jaminit(mevent,bmin,bmax,dt,nstep,
     $                             frame,proj,targ,cwin)
      nevent=mstc(2)
c     if(dump)write(33)nevent,pard(17),pard(5),pard(6),mstc(4)
      if(dump)write(33,800)nevent,pard(17),pard(5),pard(6),mstc(4)

c     if(mstc(50).ge.21) then
c       call readeos(5,fname(9))
c       iopt=mod(mstc(50),10)
c       print *,'iopt=',iopt
c       call readeos(iopt,fname(9))
c     else if(mstc(50).ge.11) then
c        feos="hgeos.dat"
c       call readEOStable(fname(9))
c     endif





c...Initialize analysis.
      call anal1

c...Simulation start.
      do iev=1,nevent

c...Simulate one event.
        call jamevt(iev)

c...Dump phase space data.
        if(dump) then

         ncol=mstd(41)+mstd(42)
         ncolbb=mstd(44)
         nonpart=0
         do i=1,nv
          if(abs(k(7,i)).eq.1) nonpart=nonpart+1
         end do
         npart=mstd(11)-nonpart

c         write(33,801)iev,nv,nbary,nmeson,pard(2)
          write(33,811)iev,nv,nbary,nmeson,pard(2),npart,ncol,ncolbb
          do i=1,nv
c           write(33,802)(k(j,i),j=1,11),(r(j,i),j=1,5),(p(j,i),j=1,5)
c           write(33,802) k(2,i),p(5,i),(p(j,i),j=1,3),(v(j,i),j=1,4)
c           write(33,802) k(2,i),p(5,i),(p(j,i),j=1,4)
c           write(33,803) k(1,i),k(2,i),p(5,i),(p(j,i),j=1,4),
c    & (r(j,i),j=1,5)
            write(33,813) k(1,i),k(2,i),k(7,i),p(5,i),(p(j,i),j=1,4),
     & (r(j,i),j=1,5)

c            if(abs(k(2,i)).eq.321) then
c              call kaondecay(i)
c            endif

          end do
        endif

        if(mod(iev,1000).eq.0) then
          p1=0.0
          p2=0.0
          p3=0.0
          p4=0.0
          do i=1,nv
           p1=p1+p(1,i) 
           p2=p2+p(2,i) 
           p3=p3+p(3,i) 
           p4=p4+p(4,i) 
          end do
          write(6,*)'event=',iev,p1,p2,p3,p4
        endif

c...Data analysis.
        call anal2

c...List phase space data.
c       call jamlist(1)

      end do

 800  format('#',1x,i8,3(1x,f11.5),1x,i4)
 801  format('#',4(1x,i8),f6.2)
 811  format('#',4(1x,i8),f6.2,1x,3(i6,1x))
c802  format(1x,i8,8(1x,1pe11.4))
 802  format(1x,i8,8(1x,1pe16.8))
c803  format(1x,i8,1x,i12,10(1x,1pe16.8))
 813  format(1x,i8,1x,i12,1x,i5,10(1x,1pe16.8))

      if(dump) close(33)

c...Final output.
      call jamfin

c...Print analysis results.
      call anal3

      print *,'mstd110=',mstd(110),mstd(111),mstd(112)
      end

c***********************************************************************

      subroutine anal1

      include 'jam1.inc'
      include 'jam2.inc'
      dimension dn(6,100)
      save ymin, ymax
      save wy,wp
      save ylab,wevt
      save nymx,npmx


c....Rapidity distribution.
      ylab=pard(17)
c     ymin=-7.0D0
c     ymax=7.0D0
c     wy=0.25D0
c     nymx=(ymax-ymin)/wy
c     print *,ylab

      ymin=-ylab*1.2
      ymax=ylab*1.2
      nymx=50
      wy=(ymax-ymin)/nymx

      wevt=1.D0/(dble(mstc(2))*wy)

      pmin=0.0D0
      pmax=5.0D0
      wp=0.1D0
      npmx=(pmax-pmin)/wp


      call vbook1(11,'n(y) protons   ',nymx,ymin,ymax)
      call vbook1(12,'n(y) antiproton',nymx,ymin,ymax)
      call vbook1(13,'n(y) pion-     ',nymx,ymin,ymax)
      call vbook1(14,'n(y) pion+     ',nymx,ymin,ymax)
      call vbook1(15,'n(y) K-        ',nymx,ymin,ymax)
      call vbook1(16,'n(y) K+        ',nymx,ymin,ymax)

c....Directed and flow as a function of rapidity.
      call vbook1(21,'v1(y) protons   ',nymx,ymin,ymax)
      call vbook1(22,'v1(y) antiproton',nymx,ymin,ymax)
      call vbook1(23,'v1(y) pion-     ',nymx,ymin,ymax)
      call vbook1(24,'v1(y) pion+     ',nymx,ymin,ymax)
      call vbook1(25,'v1(y) K-        ',nymx,ymin,ymax)
      call vbook1(26,'v1(y) K+        ',nymx,ymin,ymax)

c...Ellipse flow as a function of y.
      call vbook1(31,'v2(y) protons   ',nymx,ymin,ymax)
      call vbook1(32,'v2(y) antiproton',nymx,ymin,ymax)
      call vbook1(33,'v2(y) pion-     ',nymx,ymin,ymax)
      call vbook1(34,'v2(y) pion+     ',nymx,ymin,ymax)
      call vbook1(35,'v2(y) K-        ',nymx,ymin,ymax)
      call vbook1(36,'v2(y) K+        ',nymx,ymin,ymax)

c....Directed and flow as a function of rapidity.
      call vbook1(41,'<v1(y)> protons   ',nymx,ymin,ymax)
      call vbook1(42,'<v1(y)> antiproton',nymx,ymin,ymax)
      call vbook1(43,'<v1(y)> pion-     ',nymx,ymin,ymax)
      call vbook1(44,'<v1(y)> pion+     ',nymx,ymin,ymax)
      call vbook1(45,'<v1(y)> K-        ',nymx,ymin,ymax)
      call vbook1(46,'<v1(y)> K+        ',nymx,ymin,ymax)

c...Ellipse flow as a function of y.
      call vbook1(51,'<v2(y)> protons   ',nymx,ymin,ymax)
      call vbook1(52,'<v2(y)> antiproton',nymx,ymin,ymax)
      call vbook1(53,'<v2(y)> pion-     ',nymx,ymin,ymax)
      call vbook1(54,'<v2(y)> pion+     ',nymx,ymin,ymax)
      call vbook1(55,'<v2(y)> K-        ',nymx,ymin,ymax)
      call vbook1(56,'<v2(y)> K+        ',nymx,ymin,ymax)


c....Directed and flow as a function of rapidity.
      call vbook1(61,'<<v1(y)>> protons   ',nymx,ymin,ymax)
      call vbook1(62,'<<v1(y)>> antiproton',nymx,ymin,ymax)
      call vbook1(63,'<<v1(y)>> pion-     ',nymx,ymin,ymax)
      call vbook1(64,'<<v1(y)>> pion+     ',nymx,ymin,ymax)
      call vbook1(65,'<<v1(y)>> K-        ',nymx,ymin,ymax)
      call vbook1(66,'<<v1(y)>> K+        ',nymx,ymin,ymax)

c...Ellipse flow as a function of y.
      call vbook1(71,'<<v2(y)>> protons   ',nymx,ymin,ymax)
      call vbook1(72,'<<v2(y)>> antiproton',nymx,ymin,ymax)
      call vbook1(73,'<<v2(y)>> pion-     ',nymx,ymin,ymax)
      call vbook1(74,'<<v2(y)>> pion+     ',nymx,ymin,ymax)
      call vbook1(75,'<<v2(y)>> K-        ',nymx,ymin,ymax)
      call vbook1(76,'<<v2(y)>> K+        ',nymx,ymin,ymax)

      return

c***********************************************************************

      entry anal2

c...Loop over all particles.
      do i=1,6
      do j=1,100
        dn(i,j)=0.0d0
      end do
      end do

      do 300 i=1,nv
        if(abs(k(7,i)).eq.1) goto 300
        kf=k(2,i)
        id=0
        if(kf.eq.2212) then
          id=1
        else if(kf.eq.-2212) then
          id=2
        else if(kf.eq.-211) then
          id=3
        else if(kf.eq.211) then
          id=4
        else if(kf.eq.-321) then
          id=5
        else if(kf.eq.321) then
          id=6
        endif
        if(id.gt.0) then
        y=0.5D0*log( max(p(4,i)+p(3,i),1.D-8)/max(p(4,i)-p(3,i),1.D-8) )
          iy=(y-ymin)/wy
          if(iy.ge.1.and.iy.le.nymx) dn(id,iy)=dn(id,iy)+1.0d0
        endif
 300  continue

      do 3000 i=1,nv

c...Exclude spectetor.
        if(abs(k(7,i)).eq.1) goto 3000

        kf=k(2,i)
        kc=jamcomp(kf)
        if(kc.le.0.or.kc.gt.mstu(6)) then
           write(6,*)'Invalid code i kf kc',i,kf,kc,nv,nbary,nmeson
           goto 3000
        endif

        y=0.5D0*log( max(p(4,i)+p(3,i),1.D-8)/max(p(4,i)-p(3,i),1.D-8) )

        ptsq=p(1,i)**2+p(2,i)**2
        pt=sqrt(ptsq)
        pp=sqrt(ptsq+p(3,i)**2)
        pt=max(pt,1.D-8)
        eta=0.5D0*log( max(pp+p(3,i),1.D-8)/max(pp-p(3,i),1.D-8) )
        et=p(4,i)*pt/max(pp,1.D-8)
        px=p(1,i)
        if(pt.gt.1D-8) then
          cos1=px/pt
          cos2 = (p(1,i)**2-p(2,i)**2)/ptsq
        else
          cos1=0.0D0
          cos2=0.0d0
        endif
c       phi=acos(cos1)
c       cos2=cos(2*phi)

        id=0
        if(kf.eq.2212) then
          id=11
        else if(kf.eq.-2212) then
          id=12
        else if(kf.eq.-211) then
          id=13
        else if(kf.eq.211) then
          id=14
        else if(kf.eq.-321) then
          id=15
        else if(kf.eq.321) then
          id=16
        endif

        if(id.gt.0) then
          call vfill1(id,y,1.0d0)
          call vfill1(id+10,y,cos1)
          call vfill1(id+20,y,cos2)
          iy=(y-ymin)/wy
          if(iy.ge.1.and.iy.le.nymx) then
            idd=id-10
            call vfill1(id+50,y,cos1/dn(idd,iy))
            call vfill1(id+60,y,cos2/dn(idd,iy))
          endif
        endif

3000  continue

      return

c***********************************************************************

      entry anal3

c...Output histograms.

      fac=1.D0/dble(mstc(2))
      mnorm=0
      mform=1

c....Directed and ellipse flow.
      do i=1,6
        call vopera(20+i,'/',10+i,40+i,1.0D0,1.0D0)
        call vopera(30+i,'/',10+i,50+i,1.0D0,1.0D0)
      end do

      do i=1,6
       call vscale(10+i,wevt)
       call vprint(10+i,mnorm,mform)
       call vprint(40+i,mnorm,mform)
       call vprint(50+i,mnorm,mform)

       call vscale(60+i,fac)
       call vscale(70+i,fac)
       call vprint(60+i,mnorm,mform)
       call vprint(70+i,mnorm,mform)
      end do

      end

c***********************************************************************
      subroutine kaondecay(ip)
      include 'jam1.inc'
      include 'jam2.inc'
      common/jyjets/njet,npad,kjet(1000,5),pjet(1000,5),vjet(1000,5)
      save  /jyjets/

      mdcy(jamcomp(321),1)=1    ! K+
      mdcy(jamcomp(-321),1)=1   ! K-
      njet=1
      kjet(1,1)=1
      kjet(1,2)=k(2,ip)
      kjet(1,3)=0
      kjet(1,4)=0
      kjet(1,5)=0
      do j=1,5
          pjet(1,j)=p(j,ip)
      end do

      call pjdecy(1,icon)
      call pjlist(1)

      mdcy(jamcomp(321),1)=0    ! K+
      mdcy(jamcomp(-321),1)=0   ! K-

      end

