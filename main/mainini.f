c...A main program to use the initial condition of hadronic cascade
c...from prepared file.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15
      data ioptbinary/0/

c=========Set input values and switches ========================
c....Initialize JAM
      mstc(1)=48127      ! random seed.
      mstc(8)=0          ! job mode.
      mstc(16)=0         ! display on/off.
      parc(6)=5.0        ! scale of display
      mstc(54)=0         !avoid first coll inside the same nucleus off
      mstc(42)=0         ! weak decay on

c....Switch on some analysis.
c     mstc(156)=1        ! analysis of collision distribution
c     mstc(162)=1        ! Output collision history
c     mstc(165)=1        ! 
c     parc(7)= 1.0D0     ! Output time interval (fm/c)

c....Initial setting for JAM.
      mevent=1           ! total simulation event
      frame='user'       ! comp. frame in this case, user defined 
      bmin=0.0d0         ! minimum impact parameter (dummy)
      bmax=0.0d0         ! maximum impact parameter (dummy)
      dt=100.d0          ! collision time(fm/c)
      nstep=1            ! time step (i.e. no time step)
      cwin='10gev'       ! initial c.m. energy per nucl, in this case,
                         ! most energitic two-body collisions expected
c     proj='209Pb'       ! projectile
c     targ='209Pb'       ! target
c     proj='179Au'       ! projectile
c     targ='179Au'       ! target
      proj='32S'         ! projectile
      targ='32S'         ! target

c================ end input section ==================================

      fname(1)='jam1.cfg'

        if(ioptbinary.eq.0) then
          open(33,file='phase.dat',status='new')
        else
          open(33,file='phase.dat',status='new',
     &           access='stream',action='write')
c    &           form='unformatted',access='stream')
       endif
      
c...Initialize jam.
      call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)
      nevent=mstc(2)

        if(ioptbinary.eq.0) then
          write(33,800)nevent,pard(17),pard(5),pard(6),mstc(4)
 800  format('#',1x,i8,3(1x,f11.5),1x,i4)
        else
          write(33)nevent,pard(17),pard(5),pard(6),mstc(4)
        endif

c....Total neutron number.
      in=mstd(12)-mstd(13)

c....Total charge.
      iz=mstd(13)

c...Simulation start.
      do iev=1,nevent

c...Sampling particle momentum and coordinate.
        call init_mom

        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...Simulate one event.
        call jamevt(iev)    !...generate hadronic cascade
c       call jamfdec

c...Dump phase space data.
         ncol=mstd(41)+mstd(42)
         ncolbb=mstd(44)
         nonpart=0
         do i=1,nv
          if(abs(k(7,i)).eq.1) nonpart=nonpart+1
         end do
         npart=mstd(11)-nonpart
          if(ioptbinary.eq.0) then
            write(33,811)iev,nv,mstd(79),mstd(80),pard(2),
     &            npart,ncol,ncolbb
            do i=1,nv
              write(33,813) k(1,i),k(2,i),k(7,i),p(5,i),(p(j,i),j=1,4),
     &         (r(j,i),j=1,5)
            end do

          else

            write(33)iev,nv,mstd(79),mstd(80),pard(2),npart,ncol,ncolbb
            do i=1,nv
              write(33) k(1,i),k(2,i),k(7,i),p(5,i),(p(j,i),j=1,4),
     &         (r(j,i),j=1,5)
            end do

          endif

        flush(33)

      end do

c...Final output.
      call jamfin

      close(33)
      if(ioptbinary.eq.0) call system('gzip -f '//'phase.dat')

 811  format('#',4(1x,i8),f6.2,1x,3(i6,1x))
 813  format(1x,i8,1x,i12,1x,i5,10(1x,1pe16.8))

      end

c***********************************************************************

      subroutine init_mom

c...reading  particle momentum and position from a file.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      real*8 jamdtim
      character cfile*80

      cfile='../AngantyrInitCond.dat'
      iunit=10
c     leng=index(fname(8),' ')-1
c     cfile=chfile
c     if(leng.gt.1) cfile=fname(8)(1:leng)//chfile
c     open(unit=iunit,file=cfile,status='unknown')
      open(unit=iunit,file=cfile,status='old')

      read(iunit,*) nv1
      print *,'nv=',nv1
c     mstd(12)=0
c     mstd(13)=0
c     mstd(14)=0
      nbary=0
      nmeson=0

      do ip=1,nv1

c...Zero the vector.
        call jamzero(ip)

        read(iunit,*,end=900) kf,p(5,ip),p(1,ip),p(2,ip),p(3,ip),
     &   p(4,ip),
     &   r(1,ip),r(2,ip),r(3,ip),r(4,ip)

        nv=ip
        kc=jamcomp(kf)
        ibary=isign(kchg(kc,6),kf)

c       print *,kf,(p(i,ip),i=1,5)

c       mstd(12)=mstd(12)+ibary/3
c       mstd(13)=mstd(13)+jamchge(kf)/3
c       mstd(14)=mstd(14)+kchg(kc,7)*isign(1,kf)
        if(ibary.eq.0) then
          nmeson=nmeson+1
          xm=xm+1
        else
          nbary=nbary+1
          xb=xb+1
        endif

c....Particle mass.
        if(pmas(kc,2).le.1d-7.or.mdcy(kc,1).eq.0
     $              .or.mdcy(kc,2).eq.0.or.mdcy(kc,3).eq.0)then
          k(1,ip)=1
        else
          k(1,ip)=2
        endif
        k(2,ip)=kf
        k(3,ip)=0
        k(4,ip)=0
        k(5,ip)=-1
        k(6,ip)=0
        k(7,ip)=1
        k(8,ip)=1
        k(9,ip)=ibary
        k(10,ip)=0
        k(11,ip)=0
        r(5,ip)=r(4,ip)

c...Vertex
        v(1,ip)=r(1,ip)
        v(2,ip)=r(2,ip)
        v(3,ip)=r(3,ip)
        v(4,ip)=r(4,ip)

c.....Set resonance decay time.
        v(5,ip)=1.d+35
        if(k(1,ip).eq.2)
     $  v(5,ip)=r(4,ip)+jamdtim(1,kf,kc,k(1,ip),p(5,ip),p(4,ip))

      end do
900   continue
      end
