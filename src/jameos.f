c***********************************************************************
      function jamrev(einv,pressure,rhob,dot,gtime)
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      parameter (ecut=1.88)
      parameter (slope=1.0/3.5)
      parameter(rho0=0.168)

c     iopt=4
      iopt=1
      jamrev=0

      if(iopt.eq.1) then

c...Testing by time slice.
c*********************************************************************

c     if(gtime.ge.1.0d0.and.gtime.le.2.0d0) then
c       jamrev=0
c       if(dot.gt.0.0d0) jamrev=1
c       return
c     endif

c     jamrev=1
c     if(gtime.ge.2.0d0.and.dot.gt.0.0) return

c     if(gtime.le.1.0d0.and.dot.gt.0d0) return

c.....repulsive for time < 4
c     if(gtime.le.4.0d0.and.dot.lt.0d0) then

c.....repulsive for time < 3
c     if(gtime.le.3.0d0.and.dot.lt.0d0) then

c.....repulsive for time < 2
c     if(gtime.le.2.0d0.and.dot.lt.0d0) then

c.....repulsive for time < 7
c     if(gtime.le.7.0d0.and.dot.lt.0d0) then

c.....repulsive for time > 4
c     if(gtime.gt.4.0d0.and.dot.lt.0d0) then

c.....attractive for time < 3
c     if(gtime.le.3.0d0.and.dot.gt.0d0) then

c....attractive for time < 4 fm/c and repulsive for time > 4 fm/c
c     if(gtime.le.4.0d0) then
c          if(dot.gt.0d0) jamrev=1
c     else if(gtime.gt.4.0d0) then
c        if(dot.lt.0d0) jamrev=1

c....repulsive for time < 4 fm/c and attractive for time > 4 fm/c
      if(gtime.le.4.0d0) then
           if(dot.lt.0d0) jamrev=1
      else if(gtime.gt.4.0d0) then
         if(dot.gt.0d0) jamrev=1

c       jamrev=1
        return
      endif

      return

c*********************************************************************
      else if(iopt.eq.2) then

      if(einv.ge.3.0d0) then
        jamrev=0
        if(rn(0).le.0.3) jamrev=1
      else if(einv.ge.1.5d0) then
        jamrev=0 
        if(rn(0).ge.0.5) jamrev=1
      else if(einv.ge.0.5) then
        jamrev=1
      else
        jamrev=0
      endif
      return

c*********************************************************************
      else if(iopt.eq.3) then

      jamrev=0
c     if(einv.ge.1.5.and.rn(0).le.0.5) then
      if(rhob/rho0.ge.4.0.and.rn(0).le.0.5) then
        jamrev=1
      endif
      return

c*********************************************************************
      else

        jamrev=0
        pf=geteos(einv,rhob)
        prob=(pressure-pf)/pf

c       prob=(pressure*1.2-pf)/pf
c       prob=1.3*prob
        prob=parc(50)*prob

        if(prob.gt.0.0) then
          if(rn(0).le.prob.and.dot.gt.0.0) jamrev=1
        endif

c       if(prob.gt.0.0) then
c         if(rn(0).le.prob.and.dot.gt.0.0d0) jamrev=1
c       else 
c         if(rn(0).le.abs(prob).and.dot.lt.0.0d0) jamrev=1
c       endif

        return

      endif


c     if(einv.gt.1.5) then
c       if(pressure.ge.slope*(einv-1.64)) then
c       if(pressure.ge.0.286*(einv-1.5)+0.0675) then
c         jamrev=1
c       else
c         jamrev=0
c       endif
c       return

c       prob=min(0.8,0.5*(evinv-ecut))
c       if(rn(0).le.prob) return
c       jamrev=1
c       return

c     endif

c     jamrev=1

c    if(einv.gt.2.0.and.rn(0).le.0.3333333) return
c    if(einv.gt.2.0.and.rn(0).ge.0.5) return

      end

c***********************************************************************
      function jamclck(ibar1,ibar2,mstc87)

      jamclck=0
      if(mstc87.eq.0) then
        jamclck=1
        return
      else if(mstc87.eq.1) then  ! BB only
        if(ibar1*ibar2.eq.0) return
      else if(mstc87.eq.2) then  ! MB only
        if(ibar1*ibar2.ne.0) return
        if(ibar1.eq.0.and.ibar2.eq.0) return
      else if(mstc87.eq.3) then  ! MB and MM only
        if(ibar1*ibar2.ne.0) return
      endif
      jamclck=1

      end

c***********************************************************************
      subroutine jameosscat(i1,i2,pxr,pyr,pyz,pt,pr0)

c...Purpose: EoS modified collision by choosing attractive or repulsive orbit.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      if(mstc(59).ge.30.and.mste(42).ne.0) return
      if(jamclck(mste(32),mste(33),mstc(87)).eq.0) return

      irev=0
      if(mstc(59).eq.30) then
        if(mste(42).ne.0) return
        irev=jamrev(pare(31),pare(34),pare(33),pare(20),pard(1))
        if(irev.eq.0) return
      else if(mod(mstc(59),10).eq.1) then ! Repulsive
         if(pare(20).lt.0.0d0) irev=1
      else if(mod(mstc(59),10).eq.2) then ! Attractive 
         if(pare(20).gt.0.0d0) irev=1
      else if(mod(mstc(59),10).eq.3) then ! Attractive/Repulsive
         if(pare(20).gt.0.0d0) irev=1
         if(rn(0).ge.parc(83)) irev=1-irev ! attractive parc(83)=70%
      endif

      if(irev.eq.1) then
        pxr=-pxr
        pyr=-pyr
      else if(irev.eq.2) then
c...Find azimuth angle which gives smallest of (p'_i-p_j) dot (r_i-r_j).
        pare(18)=paru(1)+atan2(pare(26),pare(25))
        pxr=pt*cos(pare(18))
        pyr=pt*sin(pare(18))
      endif

      pare(20)=pxr*pare(25)+pyr*pare(26)+(pzr-pr0)*pare(27)
      pare(28)=pxr
      pare(29)=pyr

      end

c***********************************************************************
      subroutine jameosscats(i1,i2,pxr,pyr,pyz,pt,pr0)

c...Purpose: EoS modified collision by choosing attractive or repulsive orbit
c...for s-channel scattering.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      if(mstc(59).ge.30) then
c       i1=mste(21)
        call jamtensor(icon,i1,i2,cc,rho,rhob,einv,pfree,pxx,pyy,pzz)
        if(icon.ne.0) return
        mste(42)=icon
       endif

       pare(31)=einv
       pare(32)=rho
       pare(33)=rhob
       pare(34)=pfree
       pare(35)=(pxx+pyy)/2
       pare(36)=pzz

      irev=0
      if(mstc(59).eq.30) then
        irev=jamrev(pare(31),pare(34),pare(33),pare(20),pard(1))
        if(irev.eq.0) return
      else if(mod(mstc(59),10).eq.1) then ! Repulsive
         if(pare(20).lt.0.0d0) irev=1
      else if(mod(mstc(59),10).eq.2) then ! Attractive 
         if(pare(20).gt.0.0d0) irev=1
      endif

      if(irev.eq.1) then
        pxr=-pxr
        pyr=-pyr
      else if(irev.eq.2) then
c...Find azimuth angle which gives smallest of (p'_i-p_j) dot (r_i-r_j).
        pt=sqrt(pxr**2+pyr**2)
        pare(18)=paru(1)+atan2(pare(26),pare(25))
        pxr=pt*cos(pare(18))
        pyr=pt*sin(pare(18))
      endif

      pare(20)=pxr*pare(25)+pyr*pare(26)+(pzr-pr0)*pare(27)
      pare(28)=pxr
      pare(29)=pyr

      end

c***********************************************************************
      subroutine jameosscat2(p1,p2,pkc,iq1,iq2,pr0)

c...Purpose: EoS modified collision by choosing attractive or repulsive orbit.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension p1(15),p2(15),qm1(2),qm2(2)

      if(mstc(59).ge.30.and.mste(42).ne.0) return
c     k7c=k(7,mste(21))
c     k7b=k(7,mste(23))
c     if(mstc(130).eq.1.and.abs(k7a*k7b).eq.1) return
c     if(jamclck(k(9,mste(21)),k(9,mste(23)),mstc(87)).eq.0) return
      if(jamclck(mste(32),mste(33),mstc(87)).eq.0) return

      if(p1(1)**2+p1(2)**2.lt.1d-10) return

      irev=0
      if(mstc(59).eq.30) then
        if(mste(42).ne.0) return
c       i1=mste(21)
c       i2=mste(23)
c       call jamtensor(icon,i1,i2,cc,rho,rhob,einv,pfree,pxx,pyy,pzz)
c       if(icon.ne.0) return
        irev=jamrev(pare(31),pare(34),pare(33),pare(20),pard(1))
        if(irev.eq.0) return
      else if(mod(mstc(59),10).eq.1.and.pare(20).lt.0d0) then
          irev=1
      else if(mod(mstc(59),10).eq.2.and.pare(20).gt.0d0) then
          irev=1
      endif

      if(irev.eq.1) then 
        p1(1)=-p1(1)
        p1(2)=-p1(2)
        p1(6)=-p1(6)
        p1(7)=-p1(7)
        p1(8)=-p1(8)
        p1(9)=-p1(9)

        p2(1)=-p2(1)
        p2(2)=-p2(2)
        p2(6)=-p2(6)
        p2(7)=-p2(7)
        p2(8)=-p2(8)
        p2(9)=-p2(9)

      else if(irev.eq.2) then

c       t1=paru(1)+pjangl(pxcm,pycm)
        t1=paru(1)+atan2(pycm,pxcm)
        pare(18)=t1
        pkc11=pkc*cos(t1)
        pkc12=pkc*sin(t1)
        if(iq1.eq.1) then
          p1(6)=0d0
          p1(7)=0d0
          p1(8)=pkc11
          p1(9)=pkc12
        else
          p1(6)=pkc11
          p1(7)=pkc12
          p1(8)=0d0
          p1(9)=0d0
        endif
        if(iq2.eq.1) then
          p2(6)=0d0
          p2(7)=0d0
          p2(8)=-pkc11
          p2(9)=-pkc12
        else
          p2(6)=-pkc11
          p2(7)=-pkc12
          p2(8)=0d0
          p2(9)=0d0
        endif
        p1(1)=p1(6)+p1(8)
        p1(2)=p1(7)+p1(9)
        p2(1)=p2(6)+p2(8)
        p2(2)=p2(7)+p2(9)

      endif

      pare(20)=p1(1)*pare(25)+p1(2)*pare(26)+(p1(3)-pr0)*pare(27)
      pare(28)=p1(1)
      pare(29)=p1(2)

      end

c***********************************************************************
      subroutine jameosscat3(i1,i2,pg)

c...Purpose: EoS modified collision by choosing attractive or repulsive orbit.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common/jyjets/n,npad,kjet(1000,5),pjet(1000,5),vjet(1000,5)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jyjets/,/pjpars/,/pjint1/
      dimension pg(6,3)

      if(mstc(59).ge.30.and.mste(42).ne.0) return
      if(jamclck(mste(32),mste(33),mstc(87)).eq.0) return

      irev=0
      if(mstc(59).eq.30) then
        if(mste(42).ne.0) return
        irev=jamrev(pare(31),pare(34),pare(33),pare(20),pard(1))
        if(irev.eq.0) return
      else if(mod(mstc(59),10).eq.1) then ! Repulsive
        if(pare(20).lt.0.0d0) irev=1
      else if(mod(mstc(59),10).eq.2) then ! Attractive 
        if(pare(20).gt.0.0d0) irev=1
      endif

      if(irev.eq.0) return

      if(irev.eq.1) then

        do i=mint(84)+1,n
          pjet(i,1)=-pjet(i,1)
          pjet(i,2)=-pjet(i,2)
        end do

      else if(irev.eq.2) then

        ii=1
        do i=mint(84)+1,n
          pjet(i,1)=pg(ii,1)
          pjet(i,2)=pg(ii,2)
          pjet(i,3)=pg(ii,3)
          ii=ii+1
        end do
C...Rotate outgoing partons/particles using cos(theta).
        vint(24)=paru(1)+atan2(pare(26),pare(25))
        pare(18)=vint(24)
        if(vint(23).lt.0.9d0) then
          call pjrobo(n1,n,acos(vint(23)),vint(24),0d0,0d0,0d0)
        else
          call pjrobo(n1,n,asin(vint(59)),vint(24),0d0,0d0,0d0)
        endif

      endif

c...Total momentum from projectile side.
      px1=0d0
      py1=0d0
      pz1=0d0
      do i=mint(84)+1,n
        if(kjet(i,3).eq.3) then
          px1=px1+pjet(i,1)
          py1=py1+pjet(i,2)
          pz1=pz1+pjet(i,3)
        endif
      end do
      pare(28)=px1
      pare(29)=py1
      pare(30)=pz1-sqrt(px1**2+py1**2+pz1**2)
      pare(20)=px1*pare(25) + py1*pare(26) + pare(20)*pare(27)


c     if(n-n1.eq.1) then
c       dot=p(n-1,1)*dx+p(n-1,2)*dy
c       pare(28)=p(n-1,1)
c       pare(29)=p(n-1,2)
c       pare(30)=p(n-1,3)-pz
c       pare(20)=dot
c       irev=0
c       if(mstc(59).ge.1) then
c         if(mod(mstc(59),10).eq.1) then ! Repulsive
c          if(dot.lt.0.0d0) irev=1
c         else if(mod(mstc(59),10).eq.2) then ! Attractive 
c          if(dot.gt.0.0d0) irev=1
c         endif
c       endif
c       if(irev.eq.1) then
c         p(n-1,1)=-p(n-1,1)
c         p(n-1,2)=-p(n-1,2)
c         p(n,1)=-p(n,1)
c         p(n,2)=-p(n,2)
c         pare(20)=-dot
c         pare(28)=p(n,1)
c         pare(29)=p(n,2)
c       endif
c     endif

      end

c***********************************************************************

      subroutine jameoscols(p1,p2,pr0,iq1,iq2,emin1,emin2)

c...Purpose: EoS modified collision.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      dimension p1(15),p2(15)
c     dimension cur(4),tens(4,4)
c p1(1)-p1(5)=(p_x,p_y,p_z,E,m) four momentum in the two-body c.m.
c frame and the invariant mass of projectile hadron.

      if(mste(42).ne.0) return
      if(jamclck(mste(32),mste(33),mstc(87)).eq.0) return
      if(p1(1)**2+p1(2)**2.lt.1d-10) return
      mstd(110)=mstd(110)+1

      einv=pare(31)
      rho=pare(32)
      rhob=pare(33)
      pid=pare(34)
      if(mstc(48).eq.2) pid=pare(35)

      dot=pare(20)
      if(mstc(49).eq.0) then
        if(abs(dot).lt.1d-10) return
        if(jamrev(einv,pid,rhob,pare(20),pard(1)).eq.1) then
          p1(1)=-p1(1)
          p1(2)=-p1(2)
          p2(1)=-p2(1)
          p2(2)=-p2(2)
          dot = - dot
          do i=6,9
            p1(i)=-p1(i)
            p2(i)=-p2(i)
          end do
        endif
        pare(20)=dot
        pare(28)=p1(1)
        pare(29)=p1(2)
        pare(30)=p1(3)-pr0
        return
      endif

c...Compute pressure difference: Delta P = P - P_free
      pressure=geteos(einv,rhob)
      dt1=pare(49)
      dt2=pare(50)
      if(mstc(86).eq.1) then
        dx=pare(25)
        dy=pare(26)
        dz=pare(27)
        dpr=3*(pressure-pid)*(dt1+dt2)/rho 
c       dpr=dpr/2
      else
        dx=pare(43)/dt1 - pare(46)/dt2
        dy=pare(44)/dt1 - pare(47)/dt2
        dz=pare(45)/dt1 - pare(48)/dt2
        dpr=3*(pressure-pid)/rho
      endif
      if(dpr.lt.0d0) mstd(112)=mstd(112)+1


      if(mstc(49).eq.1) then
        call jameosang(0d0,dpr,dx,dy,dz,p1(1),p1(2),p1(3),t1,theta)
      else if(mstc(49).eq.2) then
        call jameost1(dpr,dx,dy,dz,p1(1),p1(2),p1(3),t1,theta,1)
      else
        if(dot*dpr.lt.0d0) then
          p1(1)=-p1(1)
          p1(2)=-p1(2)
        endif
        pr2=p1(1)**2+p1(2)**2+p1(3)**2
        em1=sqrt(p1(4)**2-pr2)
        em2=sqrt(p2(4)**2-pr2)
        call jameospr(dpr,p1(1),p1(2),p1(3),em1,em2,emin1,emin2,1)
        p1(5)=em1
        p2(5)=em2
        pr2=p1(1)**2+p1(2)**2+p1(3)**2
        p1(4)=sqrt(em1**2+pr2)
        p2(4)=sqrt(em2**2+pr2)
      endif

      p2(1)=-p1(1)
      p2(2)=-p1(2)
      p2(3)=-p1(3)

      if(iq1.eq.1) then
        p1(6)=0d0
        p1(7)=0d0
        p1(8)=p1(1)
        p1(9)=p1(2)
      else
        p1(6)=p1(1)
        p1(7)=p1(2)
        p1(8)=0d0
        p1(9)=0d0
      endif
      if(iq2.eq.1) then
        p2(6)=0d0
        p2(7)=0d0
        p2(8)=p2(1)
        p2(9)=p2(2)
      else
        p2(6)=p2(1)
        p2(7)=p2(2)
        p2(8)=0d0
        p2(9)=0d0
      endif

      end

c***********************************************************************

      subroutine jameoscol(i1,i2,pr0,pxr,pyr,pzr,em1,em2)

c...Purpose: EoS modified scattering
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
c---------------------------------------------------------------------

      if(mste(42).ne.0) return
      if(jamclck(mste(32),mste(33),mstc(87)).eq.0) return

      mstd(110)=mstd(110)+1

      einv=pare(31)
      rho=pare(32)
      rhob=pare(33)

c...Lorentz invariant pressure.
      pid=pare(34)

c...Use transverse part of pressure.
      if(mstc(48).eq.2) pid=pare(35)

      if(mstc(49).eq.0) then
        if(jamrev(einv,pid,rhob,pare(20),pard(1)).eq.1) then
          pxr=-pxr
          pyr=-pyr
          pare(20)=pxr*pare(25)+pyr*pare(26)+(pzr-pr0)*pare(27)
          pare(28)=pxr
          pare(29)=pyr
          pare(30)=pzr-pr0
        endif
        return
      endif

c....Proper time intervals between collisions.
      dt1=pare(49)
      dt2=pare(50)

c...Compute pressure difference: Delta P = P - P_free
      pressure=geteos(einv,rhob)

      if(mstc(86).eq.1) then
        dx=pare(25)
        dy=pare(26)
        dz=pare(27)
        dpr=3*(pressure-pid)*(dt1+dt2)/rho
c       dpr=dpr/2
      else
        dx=pare(43)/dt1 - pare(46)/dt2
        dy=pare(44)/dt1 - pare(47)/dt2
        dz=pare(45)/dt1 - pare(48)/dt2
        dpr=3*(pressure-pid)/rho
      endif

      if(dpr.lt.0d0) mstd(112)=mstd(112)+1

c...Find an angle to match the EoS.
      if(mstc(49).eq.1) then

        call jameosang(0d0,dpr,dx,dy,dz,pxr,pyr,pzr,t1,theta)

      else if(mstc(49).eq.2) then
       
        call jameost1(dpr,dx,dy,dz,pxr,pyr,pzr,t1,theta,2)

c...Modify effective masses of particle instead of changing angle.
      else

        dot=pare(20)
        if(abs(dot).lt.1d-10) return
        if(abs(pressure-pid).lt.1d-4) return

c...Set the orbit consistent with the sign of pressure difference.
        if(dpr*dot.lt.0.0d0) then
          pxr=-pxr
          pyr=-pyr
          pare(20)=pxr*pare(25)+pyr*pare(26)+(pzr-pr0)*pare(27)
          pare(20)=dot
          pare(28)=pxr
          pare(29)=pyr
          pare(30)=pzr-pr0
        endif

        emin1=0.134
        emin2=0.134
        call jameospr(dpr,pxr,pyr,pzr,em1,em2,emin1,emin2,2)

      endif

      end

c***********************************************************************

      subroutine jameoscolsch(i1,i2,pr0,pxr,pyr,pzr,em1,em2)

c...Purpose: EoS modified scattering for s-channel scatterings.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
c---------------------------------------------------------------------

c     if(pare(25)**2+pare(26)**2+pare(27)**2.lt.1d-7) return
c     i1=mste(21)
      mste(42)=0
      call jamtensor(icon,i1,i2,cc,rho,rhob,einv,pfree,pxx,pyy,pzz)
       mste(42)=icon
      if(mste(42).ne.0) return
       pare(31)=einv
       pare(32)=rho
       pare(33)=rhob
       pare(34)=pfree
       pare(35)=(pxx+pyy)/2
       pare(36)=pzz

      mstd(110)=mstd(110)+1

      einv=pare(31)
      rho=pare(32)
      rhob=pare(33)

c...Lorentz invariant pressure.
      pid=pare(34)

c...Use transverse part of pressure.
      if(mstc(48).eq.2) pid=pare(35)

      if(mstc(49).eq.0) then
        if(jamrev(einv,pid,rhob,pare(20),pard(1)).eq.1) then
          pxr=-pxr
          pyr=-pyr
          pare(20)=pxr*pare(25)+pyr*pare(26)+(pzr-pr0)*pare(27)
          pare(28)=pxr
          pare(29)=pyr
          pare(30)=pzr-pr0
        endif
        return
      endif

c....Proper time intervals between collisions.
c     dt1=pare(49)
c     dt2=pare(50)
      dt1=pare(49)
      dt2=0.0d0

c...Compute pressure difference: Delta P = P - P_free
      pressure=geteos(einv,rhob)

c     if(mstc(86).eq.1) then
        dx=pare(25)
        dy=pare(26)
        dz=pare(27)
        dpr=3*(pressure-pid)*(dt1+dt2)/rho
c     endif

      if(dpr.lt.0d0) mstd(112)=mstd(112)+1

c...Find an angle to match the EoS.
      if(mstc(49).eq.1) then

        call jameosang(0d0,dpr,dx,dy,dz,pxr,pyr,pzr,t1,theta)

      else if(mstc(49).eq.2) then
       
        call jameost1(dpr,dx,dy,dz,pxr,pyr,pzr,t1,theta,4)

c...Modify effective masses of particle instead of changing angle.
      else

        dot=pare(20)
        if(abs(dot).lt.1d-10) return
        if(abs(pressure-pid).lt.1d-4) return

c...Set the orbit consistent with the sign of pressure difference.
        if(dpr*dot.lt.0.0d0) then
          pxr=-pxr
          pyr=-pyr
          pare(20)=pxr*pare(25)+pyr*pare(26)+(pzr-pr0)*pare(27)
          pare(20)=dot
          pare(28)=pxr
          pare(29)=pyr
          pare(30)=pzr-pr0
        endif

        emin1=0.134
        emin2=0.134
        call jameospr(dpr,pxr,pyr,pzr,em1,em2,emin1,emin2,4)

      endif

      end


c***********************************************************************
      subroutine jameosdec(i1,i2,mtype)

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'

c     if(pare(25)**2+pare(26)*2+pare(27)**2.lt.1d-18) return

        px1=p(1,i1)
        py1=p(2,i1)
        pz1=p(3,i1)
        pe1=p(4,i1)
        em1=sqrt(max(0d0,pe1**2-px1**2-py1**2-pz1**2))
        px2=p(1,i2)
        py2=p(2,i2)
        pz2=p(3,i2)
        pe2=p(4,i2)
        em2=sqrt(max(0d0,pe2**2-px2**2-py2**2-pz2**2))

        ee=pe1+pe2
        be1=(px1+px2)/ee
        be2=(py1+py2)/ee
        be3=(pz1+pz2)/ee
        gg=ee/sqrt(ee**2 - (px1+px2)**2-(py1+py2)**2-(pz1+pz2)**2)
        call jamrobo(0d0,0d0,-be1,-be2,-be3,gg,px1,py1,pz1,pe1)
        pr2=px1**2+py1**2+pz1**2
        pr=sqrt(pr2)
        pt=sqrt(px1**2+py1**2)


c       iopt=0
c       if(iopt.eq.0) then
        if(mstc(91).eq.1) then
        rx1=r(1,i1)
        ry1=r(2,i1)
        rz1=r(3,i1)
        rt1=r(4,i1)
        rx2=r(1,i2)
        ry2=r(2,i2)
        rz2=r(3,i2)
        rt2=r(4,i2)
        call jamrobo(0d0,0d0,-be1,-be2,-be3,gg,rx1,ry1,rz1,rt1)
        call jamrobo(0d0,0d0,-be1,-be2,-be3,gg,rx2,ry2,rz2,rt2)
c       the=pjangl(px1,py1)
c       phi=pjangl(pz1,pt)
c       call jamrobo(0d0,-the,0d0,0d0,0d0,1d0,rx1,ry1,rz1,rt1)
c       call jamrobo(-phi,0d0,0d0,0d0,0d0,1d0,rx1,ry1,rz1,rt1)
c       call jamrobo(0d0,-the,0d0,0d0,0d0,1d0,rx2,ry2,rz2,rt2)
c       call jamrobo(-phi,0d0,0d0,0d0,0d0,1d0,rx2,ry2,rz2,rt2)

        pare(25)=rx1-rx2
        pare(26)=ry1-ry2
        pare(27)=rz1-rz2
        endif

        if(pare(25)**2+pare(26)**2+pare(27)**2.lt.1d-25) then
           print *,'jameosdec dx=',pare(25),pare(26),pare(27)
           print *,k(2,i1),(r(j,i1),j=1,5)
           print *,k(2,i2),(r(j,i2),j=1,5)
           print *,'parc42',parc(42),'mtype=',mtype
           return
        endif
c       pare(49)=taucl ! this should be set before call this!
        pare(18)=pjangl(px1,py1)
        pare(19)=pz1/pr
        pare(20)=px1*pare(25) + py1*pare(26) + pz1*pare(27)

        if(mstc(59).ge.101) then
          call jameoscolsch(i1,i2,pr,px1,py1,pz1,em1,em2)
          pr2=px1**2+py1**2+pz1**2
          pe1=sqrt(em1*em1 + pr2)
          pe2=sqrt(em2*em2 + pr2)
        else
          call jameosscats(i1,i2,px1,py1,pz1,pt,pr)
        endif

        px2=-px1
        py2=-py1
        pz2=-pz1
        pe2=sqrt(em2*em2 + px2**2+py2**2+pz2**2)
        call jamrobo(0d0,0d0,be1,be2,be3,gg,px1,py1,pz1,pe1)
        p(1,i1)=px1
        p(2,i1)=py1
        p(3,i1)=pz1
        p(4,i1)=sqrt(em1*em1 + px1**2+py1**2+pz1**2)
        call jamrobo(0d0,0d0,be1,be2,be3,gg,px2,py2,pz2,pe2)
        p(1,i2)=px2
        p(2,i2)=py2
        p(3,i2)=pz2
        p(4,i2)=sqrt(em2*em2 + px2**2+py2**2+pz2**2)

c       print *,'e1 e2',pe1,pe2
c       print *,'p',px1+px2,py1+py2,pz1+pz2
c       print *,'px1=',px1+px2,py1+py2,pz1+pz2,p(4,i1)+p(4,i2)

      end

c***********************************************************************
      subroutine jameoscold(i1,i2,n1,pxr,pyr,pzr,pg)

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common/jyjets/n,npad,kjet(1000,5),pjet(1000,5),vjet(1000,5)
      common/pjpars/mstp(200),parp(200),msti(200),pari(200)
      common/pjint1/mint(400),vint(400)
      save /jyjets/,/pjpars/,/pjint1/
      dimension pg(6,3)

c     call jamtensor(icon,i1,i2,cc,rho,rhob,einv,pfree,pxx,pyy,pzz)
c     if(icon.ne.0) return
c     pare(31)=einv
c     pare(32)=rho
c     pare(33)=pfree

      if(mste(42).ne.0) return
      if(jamclck(mste(32),mste(33),mstc(87)).eq.0) return

      einv=pare(31)
      rho=pare(32)
      rhob=pare(33)
      pid=pare(34)
      if(mstc(48).eq.2) pid=pare(35)

      dx=pare(25)
      dy=pare(26)
      dz=pare(27)
      pr=sqrt(pxr**2 + pyr**2 + pzr**2)
      dot=pare(20)

c...Set the orbit consistent with the sign of pressure difference.
      if(mstc(49).eq.0) then
        if(jamrev(einv,pid,rhob,pare(20),pard(1)).eq.1) then
            pxr=-pxr
            pyr=-pyr
            dot=pxr*dx+pyr*dy+(pzr-pr)*dz
            do i=n1,n
              pjet(i,1)=-pjet(i,1)
              pjet(i,2)=-pjet(i,2)
            end do
        endif
        pare(20)=dot
        pare(28)=pxr
        pare(29)=pyr
        pare(30)=pzr-pr
        return
      endif

c...Compute pressure difference: Delta P = P - P_free
      dt1=pare(49)
      dt2=pare(50)
      pressure=geteos(einv,rhob)
      if(mstc(86).eq.1) then
        dx=pare(25)
        dy=pare(26)
        dz=pare(27)
        dpr=3*(pressure-pid)*(dt1+dt2)/rho
c       dpr=dpr/2
      else
        dx=pare(43)/dt1 - pare(46)/dt2
        dy=pare(44)/dt1 - pare(47)/dt2
        dz=pare(45)/dt1 - pare(48)/dt2
        dpr=3*(pressure-pid)/rho
      endif
      if(dpr.lt.0d0) mstd(112)=mstd(112)+1

c...Find an angle to match the EoS.
      if(mstc(49).le.2) then

      if(mstc(49).eq.1) then
c       pr=sqrt(pxr**2 + pyr**2 + pzr**2)
c       phi=pjangl(pxr,pyr)
c       the=pjangl(pzr,sqrt(pxr**2+pyr**2))
        phi=vint(24)
        the=acos(vint(23))
c       print *,'the=',the,the1
c       print *,'phi',phi,phi1
        call jameosang(0d0,dpr,dx,dy,dz,pxr,pyr,pzr,t1,theta)
      else if(mstc(49).eq.2) then
        call jameost1(dpr,dx,dy,dz,pxr,pyr,pzr,t1,theta,3)
      endif

      px1=0d0
      py1=0d0
      pz1=0d0
      e1=0d0
      l=1
      px2=0d0
      py2=0d0
      pz2=0d0
      e2=0d0
      ptotx=0d0
      ptoty=0d0
      ptotz=0d0
      do i=n1,n
        pjet(i,1)=pg(l,1)
        pjet(i,2)=pg(l,2)
        pjet(i,3)=pg(l,3)
        ptotx=ptotx+pjet(i,1)
        ptoty=ptoty+pjet(i,2)
        ptotz=ptotz+pjet(i,3)
        l=l+1
        call jamrobo(theta,t1,0d0,0d0,0d0,1d0,
     &         pjet(i,1),pjet(i,2),pjet(i,3),0d0)

        if(mstc(49).eq.1)
     &    call jamrobo(the,phi,0d0,0d0,0d0,1d0,
     &         pjet(i,1),pjet(i,2),pjet(i,3),0d0)

        if(kjet(i,3).eq.3) then
          px1=px1+pjet(i,1)
          py1=py1+pjet(i,2)
          pz1=pz1+pjet(i,3)
          e1=e1+pjet(i,4)
        else
          px2=px2+pjet(i,1)
          py2=py2+pjet(i,2)
          pz2=pz2+pjet(i,3)
          e2=e2+pjet(i,4)
        endif
      end do

      dpx=px1+px2-ptotx
      dpy=py1+py2-ptoty
      dpz=pz1+pz2-ptotz

      if(abs(pxr-px1).gt.1e-3.or.abs(pyr-py1).gt.1e-3
     & .or.abs(pzr-pz1).gt.1e-3) then
          print *,'eoscold pr?',t1,theta
          print *,pxr,pyr,pzr
          print *,px1,py1,pz1
          stop
      endif

      if(abs(dpx).gt.1e-3.or.abs(dpy).gt.1e-3.or.abs(dpz).gt.1e-3) then
        mstu(11)=6
        call pjlist(1)
        print *,'jameoscold momentum does not conserve',dpx,dpy,dpz
        print *,'theta=',theta,t1
        do i=1,n-n1+1
        print *,'pg=',pg(i,1),pg(i,2),pg(i,3)
        end do
      endif

      pare(20)=px1*pare(25)+py1*pare(26)+(pz1-pr)*pare(27)
      pare(28)=px1
      pare(29)=py1
      pare(30)=pz1-pr

      if(mstc(49).eq.1) then
C...Rotate outgoing partons/particles using cos(theta).
        vint(24)=paru(2)*pjr(0)
      if(vint(23).lt.0.9d0) then
        call pjrobo(n1,n,acos(vint(23)),vint(24),0d0,0d0,0d0)
      else
        call pjrobo(n1,n,asin(vint(59)),vint(24),0d0,0d0,0d0)
      endif

      endif


c...Modify effective masses of particle instead of changing angle.
      else

        if(dpr*dot.lt.0.0d0) then
            pxr=-pxr
            pyr=-pyr
            dot=pxr*pare(25)+pyr*pare(26)+(pzr-pr)*pare(27)
            do i=n1,n
              pjet(i,1)=-pjet(i,1)
              pjet(i,2)=-pjet(i,2)
            end do
        endif
        pare(20)=dot
        pare(28)=pxr
        pare(29)=pyr
        pare(30)=pzr-pr
        if(abs(pressure-pid).lt.1d-4) return

      px1=0d0
      py1=0d0
      pz1=0d0
      e1=0d0
      px2=0d0
      py2=0d0
      pz2=0d0
      e2=0d0
      do i=n1,n
        if(kjet(i,3).eq.3) then
          px1=px1+pjet(i,1)
          py1=py1+pjet(i,2)
          pz1=pz1+pjet(i,3)
          e1 =e1 +pjet(i,4)
        else
          px2=px2+pjet(i,1)
          py2=py2+pjet(i,2)
          pz2=pz2+pjet(i,3)
          e2 =e2 +pjet(i,4)
        endif
      end do
      em1=sqrt(e1**2-px1**2-py1**2-pz1**2)
      em2=sqrt(e2**2-px2**2-py2**2-pz2**2)

        emin1=0.134
        emin2=0.134
        call jameospr(dpr,px1,py1,pz1,em1,em2,emin1,emin2,3)
        if(n-n1.eq.1) then
          pjet(n1,1)=px1
          pjet(n1,2)=py1
          pjet(n1,3)=pz1
          pjet(n,1)=-px1
          pjet(n,2)=-py1
          pjet(n,3)=-pz1
          pr2=px1**2+py1**2+pz1**2
          pjet(n1,4)=sqrt(em1**2+pr2)
          pjet(n,4)=sqrt(em2**2+pr2)
        endif

      endif

c     mstu(11)=6
c     call pjlist(1)

      end

c***********************************************************************
      subroutine jamang(icol,pr,c1)

      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      if(icol.eq.2) then

      srt=pare(2)
      pr0=pare(7)
      ibar1=mste(32)
      ibar2=mste(33)
      kc3=mste(26)
      kc4=mste(28)
      kc1=mste(22)
      kc2=mste(22)
      kf1=kchg(kc1,4)
      kf2=kchg(kc2,4)
      em1=pare(41)
      em2=pare(42)

      if(mste(1).eq.1) then
          sig=pare(4)
          call jamangel(pr,sig,kf1,kf2,ibar1,ibar2,t1,c1) 
      else
        snew=sqrt(em1**2+pr*pr)+sqrt(em2**2+pr*pr)
        call jamangin(snew,pr,pr0,em1,em2,ibar1,ibar2,
     $                kf1,kf2,kc1,kc2,kc3,kc4,t1,c1)
       endif

      else if(icol.eq.1) then
         print *,'under construction icol=',icol
         stop
        call jamsftpt(0,pare(40),pkc)

      else if(icol.eq.3) then
         print *,'under construction icol=',icol
         stop

      endif

      end

c***********************************************************************
      subroutine jameost1(dpr,dx,dy,dz,px,py,pz,t1,the,icol)

c...Find phi angle to modify the orbit. Scattering angle will not be
c...changed.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      pr=sqrt(px**2 + py**2 + pz**2)
      pt=sqrt(px**2+py**2)

c....Scattering angle.
      the=pjangl(pz,sqrt(px**2+py**2))
      cos1=pare(19)
      if(cos1.ge.1d0) then
        print *,'jameost1(1) cos1=',cos1
        return
      endif
      if(pr.le.0d0) then
        print *,'jameost1 pr=',pr
        return
      endif

      sin1=sqrt(1.d0-cos1**2)

      if(sin1.le.0d0) then
        print *,'jameost1 sin1=',sin1,cos1
        return
      endif


      if(abs(cos1-pz/pr).gt.1d-3) then
        print *,'jameost1 ? cos1=',cos1,pz/pr
      endif
      if(abs(cos1).eq.1d0) then
        print *,'scatt. ang. jameos1 cos1=',cos1,pz/pr,
     &         mste(1),mste(2),icol
        return
      endif

c     print *,'icol=',icol
c     print *,'pt=',pt
c     print *,'pr*sin',pr*sin1
c     print *,'t1=',t1,phi,atan2(px,py),atan2(py,px)
c     the=atan2(pt,pz)
c     cos2=pz/pr
c     print *,'cos(t1)=',cos1,cos(the),cos2
c     print *,'sin(t1)=',sin1,sin(the),sqrt(1d0-cos2**2)
c     print *,'theta=',c1,cos1

      mste41=0
      ntry=0
100   continue
      ntry=ntry+1
      c=(dpr/pr + dz - dz*cos1)/sin1
      d=c/sqrt(dx**2 + dy**2)
      d0=d

c...Try different scattering angle.
      if(abs(d).gt.1.0d0) then

        if(icol.eq.2.and.ntry.le.30) then
          call jamang(icol,pr,cos1)
          sin1=sqrt(1.d0-cos1**2)
          if(sin1.gt.0.0d0) then
          pz=pr*cos1
          pt=pr*sin1
          goto 100
c         else
c           print *,'jameost1 cos1=',cos1
          endif
        endif

        d=1.0*sign(1.0d0,d)
        mstd(111)=mstd(111)+1
        mste41=1
      endif
c     if(ntry.ge.100) then
c       print *,'scattering angle ntry',ntry,cos1
c     endif

c.....Find phi angle.
c       if(abs(dy).gt.1d-5) then
c         alpha=atan2(dx,dy)
c         t1=asin(d)-alpha
c      else 
          alpha=atan2(dy,dx)
          t1=acos(d)+alpha
c       endif

      px=pt*cos(t1)
      py=pt*sin(t1)

c     px=pr*sin1*cos(t1)
c     py=pr*sin1*sin(t1)
c     pz=pr*cos1

c     dot0=pare(20)

      pare(28)=px
      pare(29)=py
      pare(30)=pz-pr
      pare(20)=pare(28)*dx + pare(29)*dy + pare(30)*dz
      pare(18)=t1

c     print *,'icol time',icol,pard(1)
c     print *,'sin d',sin(t1+alpha),dz
c     t2=acos(d)+atan2(dy,dx)
c     px2=pt*cos(t2)
c     py2=pt*sin(t2)
c     dot2=px2*dx + py2*dy + pare(30)*dz
c     print *,'t1=',t1,t2
c     print *,'px=',px,px2
c     print *,'py=',py,py2
c     print *,'d=',d,'d0=',d0,'dpr=',dpr
c     print *,'dot',dot0,pare(20)

      if(mste41.eq.0.and.abs(dpr-pare(20)).gt.1d-5) then
         print *,'icol=',icol,dpr,pare(20)
         stop
      endif

      end

c***********************************************************************
      subroutine jameosang(dot,dpr,dx,dy,dz,px,py,pz,t1,theta)

c...Find theta to modify the orbit by elastic collision.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      pr=sqrt(px**2 + py**2 + pz**2)
      phi=pjangl(px,py)
      the=pjangl(pz,sqrt(px**2+py**2))

      b=dz
      c=dpr/pr+dz-dot

c...Solve a*sin(theta)+b*cos(thea)=c
c...s*cos(theta+alpha)=c, s=sqrt(a**2+b**2), tan(alpha)=a/b
      do i=1,50
        t1=2.0d0*paru(1)*rn(0)
        if(i.eq.50) then
          if(abs(dx).ge.1e-7) then
            t1=atan(dy/dx)
          else
            t1=paru(1)/2 - atan(dx/dy)
          endif
        endif
        a=dx*cos(t1)+dy*sin(t1)
        s=c/sqrt(a**2+b**2)
        if(abs(s).le.1d0) goto 300
      end do
      s=1.0*sign(1.0d0,s)
      mstd(111)=mstd(111)+1
 300  continue

c.....Rotate momentum vector to new position.
      if(abs(a).gt.1d-5) then
        alpha=atan2(b,a)
        theta=asin(s)-alpha
      else
        alpha=atan2(a,b)
        theta=acos(s)+alpha
      endif

c     cos1=abs(cos(theta))
c     a=cos(the)*sin(theta)*cos(t1)+cos1*sin(the)
c     b=sin(theta)*sin(t1)
c     px=pr*(cos(phi)*a - sin(phi)*b)
c     py=pr*(sin(phi)*a + cos(phi)*b)
c     pz=pr*(cos(the)*cos1 - sin(the)*sin(theta)*cos(t1))

      pt=pr*sin(theta)
      px=pt*cos(t1)
      py=pt*sin(t1)
      pz=pr*cos(theta)

      pare(28)=px
      pare(29)=py
      pare(30)=pz-pr
      pare(20)=pare(28)*dx + pare(29)*dy + pare(30)*dz

      call jamrobo(the,phi,0d0,0d0,0d0,1d0,px,py,pz,0d0)


c     print *,'dxt time',icol,pard(1),the,phi
c     print *,pare(28),pare(29),pare(30)
c     print *,'px',px,py,pz
c     print *,'dot',pare(20),dpr,d

      end

c***********************************************************************
      subroutine jameospr(dpr,pxr,pyr,pzr,em1,em2,pm1,pm2,msel)
c...Change the momentum and mass modification factor to match the EoS.

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
   
      srt=pare(2)
      s=srt**2
      pr=sqrt(pxr**2 + pyr**2 + pzr**2)
      dr=sqrt(pare(25)**2+pare(26)**2+pare(27)**2)
      dot=pxr*pare(25) + pyr*pare(26) + (pzr-pr)*pare(27)

      if(abs(dot).le.1d-10) then
        print *,'dot=0?',dot,'srt=',srt,'msel=',msel
        print *,'p=',pxr,pyr,pzr,pr
        print *,'r=',pare(25),pare(26),pare(27)
        return
      endif

      cos1=dot/(pr*dr)
      prnew=dpr/cos1

      if(prnew.le.0d0) then
         print *,'prnew<0?',prnew,msel
         print *,'dpr=',dpr,'dot=',dot
         print *,'pr=',pr,'dr=',dr
         print *,'p=',pxr,pyr,pzr
         print *,'pzr-pr',pzr-pr,pare(27),mstc(52)
         return
      endif

      prmin=sqrt(max(0d0,(s-pm1**2-pm2**2)**2-(2d0*pm1*pm2)**2))
     &/(2d0*srt)

      if(prnew.ge.prmin) then
c        print *,'prnew too big?',prnew,pr,dpr
         prnew = prmin-0.003
      endif


      if(abs(em1-em2).lt.1d-10) then
          s2=(s-4*prnew**2)/(2*(em1**2+em2**2))
      else
c         eem=(em1+em2)**2*(em1-em2)**2
          eem=(em1**2 - em2**2)**2
          d=s*(em1**2 + em2**2)**2 - eem*(s-4*prnew*prnew)
          if(d.lt.0d0) then
            print *,'no solution? ',d
            return
          endif
          s2a=(s*(em1**2+em2**2)+srt*sqrt(d))/eem
          s2=(s*(em1**2+em2**2)-srt*sqrt(d))/eem

      endif

      if(s2.le.0.0d0) then
        print *,'s2<0?',s2
        print *,'s2a=',s2a
        e1=(s2a*em1**2+prnew**2)
        e2=(s2a*em2**2+prnew**2)
        print *,'srt=',srt,e1+e2
        print *,'pm1 pm2',pm1,pm2
        print *,'em1 em2',em1,em2
        return
      endif

      pxr=prnew*pxr/pr
      pyr=prnew*pyr/pr
      pzr=prnew*pzr/pr

      em1=sqrt(s2)*em1
      em2=sqrt(s2)*em2

      pare(28)=pxr
      pare(29)=pyr
      pare(30)=pzr-pr
      pare(20)=pxr*pare(25) + pyr*pare(26) + (pzr-pr)*pare(27)

c     pr2=pxr**2+pyr**2+pzr**2
c     e1=sqrt(em1**2+pr2)
c     e2=sqrt(em2**2+pr2)
c     if(abs(srt-e1-e2).gt.1d-4) then
c       print *,'jameospr energy does not conserve',srt,e1+e2
c       stop
c     endif

      end

c***********************************************************************
      subroutine jamtensor(icon,i1,i2,cc,rho,rhob,einv,
     & pfree,pxx,pyy,pzz)

c...Purpose: compute energy-momentum tensor and hydrodynamics velocity.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      common/jamdecmom/po(5),xo(4)
      save /jamdecmom/
c     parameter(widg=2*1.0d0)                   ! Gaussian width
c     parameter(widg=2*1.2d0**2)                   ! Gaussian width
c     parameter(widcof=(1.d0/(3.14159d0*widg))**1.5d0)
      dimension pv(4)
      dimension g(4,4),cur(4),curb(4),tens(4,4),u(4)
      dimension gg(4,4),dl(4,4),z(4)
      data g/ -1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,1/
      data gg/ 1,1,1,-1,  1,1,1,-1, 1,1,1,-1, -1,-1,-1,1/
      data dl/ 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1/

      widcof=pard(90)
      widg=2*parc(84)**2

      icon=1
      einv=0.0d0
      rho=0d0   ! 2018/1/13
      rhob=0d0  ! 2018/1/13
      pfree=0.0d0
      pxx=0.0d0
      pyy=0.0d0
      pzz=0.0d0
      pperp=0d0
      plong=0d0
      isel=0

      if(i2.eq.0) then
c       x2=r(1,i1)
c       y2=r(2,i1)
c       z2=r(3,i1)
c       time = r(4,i1)
c       bex=p(1,i1)/p(4,i1)
c       bey=p(2,i1)/p(4,i1)
c       bez=p(3,i1)/p(4,i1)
c       emdsq=p(4,i1)**2-p(1,i1)**2-p(2,i1)**2-p(3,i1)**2
c       gam=p(4,i1)/sqrt(emdsq)
        x2=xo(1)
        y2=xo(2)
        z2=xo(3)
        time = xo(4)
        bex=po(1)/po(4)
        bey=po(2)/po(4)
        bez=po(3)/po(4)
        gam=po(4)/po(5)

        if(mstc(88).ge.1) then
          isel=2
          if(abs(k(9,i1)).eq.3) isel=1
        endif
      else

        if(mstc(91).ge.1) then
        if(abs(k(7,i1)).eq.1) return   ! not yet interact
        if(abs(k(7,i2)).eq.1) return   ! not yet interact
        endif

        x2=(r(1,i1)+r(1,i2))/2
        y2=(r(2,i1)+r(2,i2))/2
        z2=(r(3,i1)+r(3,i2))/2
        time = (r(4,i1)+r(4,i2))/2
        bex=pare(12)
        bey=pare(13)
        bez=pare(14)
        gam=pare(15)
        if(mstc(88).ge.1) then
          if(mste(2).eq.1) isel=1 ! BB collision
          if(mste(2).eq.2) isel=1 ! MB collision
          if(mste(2).eq.3) isel=2 ! MM collision
          if(mste(2).eq.4) isel=1 ! antiBB collision
        endif
      endif

c     gamz=1d0/sqrt(1d0-bez**2)
      ycm=0.5*log((1d0+bez)/(1d0-bez))

      do i=1,4
        cur(i)=0.0d0
        curb(i)=0.0d0
        do j=1,4
          tens(i,j)=0.0d0
      end do
      end do

      ncount=0
c....Loop over all particles
      do i=1,nv
 
        if(i.eq.i1) goto 100
        if(i.eq.i2) goto 100
        k1=k(1,i)
        if(k1.gt.10) goto 100   ! dead particle
c       if(p(5,i).le.1d-5) goto 100


c.....Exclude mesons.
        if(isel.eq.1.and.k(9,i).eq.0) then
          goto 100
c.....Exclude baryons.
        else if(isel.eq.2.and.abs(k(9,i)).eq.3) then
          goto 100
        endif

c.....pre-formed hadrons.
        if(mstc(89).eq.1.and.r(5,i).gt.time) goto 100

        if(abs(k(7,i)).eq.1) goto 100   ! not yet interact

        dt=time-r(4,i)
        if(dt.lt.0.0d0) goto 100


        if(k(2,i).eq.92) then
          do j=1,4
            pv(j)=vq(j,i)+vq(j+5,i)
          end do
        else
          pv(1)=p(1,i)
          pv(2)=p(2,i)
          pv(3)=p(3,i)
          pv(4)=p(4,i)
        endif

        y1=0.5*log((pv(4)+pv(3))/(pv(4)-pv(3)))
        if(abs(y1-ycm).gt.parc(30)) goto 100

        vx=pv(1)/pv(4)
        vy=pv(2)/pv(4)
        vz=pv(3)/pv(4)
        x1=r(1,i) + dt*vx - x2
        y1=r(2,i) + dt*vy - y2
        z1=r(3,i) + dt*vz - z2



c       bar=kfprop(k(2,i),2)/3.0d0
        bar=k(9,i)/3.0d0
        iq=1
        facq=1.0d0
        if(r(5,i).gt.time) then
          iq=mod(abs(k(1,i))/10,10)
          if(iq.eq.3) iq=2
          bar=iq/3d0*isign(1,k(2,i))
          if(k(9,i).eq.0)bar=abs(bar)  ! meson can have one const.q.
          if(abs(k(9,i)).eq.3) qnum1=3.d0
          if(k(9,i).eq.0) qnum1=2.d0
          facq=iq/qnum1
        endif

        emd=pv(4)**2-pv(1)**2-pv(2)**2-pv(3)**2
        if(emd.le.1d-5) goto 100
        emd=sqrt(emd)

        xtra=x1**2+y1**2+z1**2
     $          +((x1*pv(1)+y1*pv(2)+z1*pv(3))/emd)**2
        if(xtra/widg.gt.30.d0) goto 100

        gam=pv(4)/emd
        den=widcof*gam*exp(-xtra/widg)
        bden=den*bar
        den=den*facq

c...Compute current and energy-momentum tensor.
        do im=1,4
          cur(im)  = cur(im)  + pv(im)/pv(4)*den
          curb(im) = curb(im) + pv(im)/pv(4)*bden
        do ik=1,4
          tnsmn = pv(im)*pv(ik)/pv(4)
          tens(im,ik) = tens(im,ik)  + tnsmn*den
        end do
        end do
        ncount=ncount+1

100   end do   ! End loop over particles.

      if(ncount.eq.0) return
      if(cur(4).le.1d-7) return

      cc=cur(4)**2-(cur(1)**2+cur(2)**2+cur(3)**2)
      if(cc.le.1d-7) return
      cc=sqrt(cc)
      u(1)=cur(1)/cc
      u(2)=cur(2)/cc
      u(3)=cur(3)/cc
      u(4)=cur(4)/cc

      if(mstc(47).eq.2) call landauframe(u,tens)

c...Lorentz invariant Scalar Number Density
      rho=cur(4)*u(4)-cur(1)*u(1)-cur(2)*u(2)-cur(3)*u(3)

c...Lorentz invariant Baryon density
      rhob=curb(4)*u(4)-curb(1)*u(1)-curb(2)*u(2)-curb(3)*u(3)

c...Lorentz invariant pressure and energy density
      gam=u(4)
      vz=u(3)/u(4)
      gamz=1.0d0/sqrt(1.0d0-vz**2)
      z(4)=gamz*vz
      z(3)=gamz
      z(2)=0d0
      z(1)=0d0
      do i=1,4
      do j=1,4
       tmp=gg(i,j)*u(i)*u(j)
       einv=einv + tens(i,j)*tmp                    ! energy density
       pfree=pfree - 1.d0/3.d0*tens(i,j)*(g(i,j)-tmp) ! pressure

       pl=gg(i,j)*z(i)*z(j)
       plong=plong+tens(i,j)*pl
       pperp=pperp-0.5d0*tens(i,j)*(g(i,j)-tmp+pl)

c...Txx,Tyy,Tzz at the local rest frame.
       ci=-1.0d0
       cj=-1.0d0
       if(i.le.3) ci=u(i)/(1+gam)
       if(j.le.3) cj=u(j)/(1+gam)
       pxx=pxx+tens(i,j)*(dl(1,i)+u(1)*ci)*(dl(1,j)+u(1)*cj)
       pyy=pyy+tens(i,j)*(dl(2,i)+u(2)*ci)*(dl(2,j)+u(2)*cj)
       pzz=pzz+tens(i,j)*(dl(3,i)+u(3)*ci)*(dl(3,j)+u(3)*cj)

      end do
      end do

      icon=0

c     print *,'p=',pfree,sqrt(pxx**2+pyy**2),pzz
c     print *,'p=',pfree,pperp,plong

c...Lorentz transform to the local rest frame to check.
c     be1=-u(1)/u(4)
c     be2=-u(2)/u(4)
c     be3=-u(3)/u(4)
c     call jamlemt(be1,be2,be3,gam,tens)
c     print *,'gam=',gam
c     print *,'p=',pfree,(pxx+pyy+pzz)/3,
c    &(tens(1,1)+tens(2,2)+tens(3,3))/3
c     print *,'pl=',pl,-einv+pl
c     print *,'pzz=',pzz,-einv-pzz
c     print *,'pt=',pperp
c     print *,'pzz=',(pxx+pyy)/2,-einv-(pxx+pyy)/2
c     print *,'einv=',einv,tens(4,4)
c     print *,'T=',tens(4,1),tens(4,2),tens(4,3)
c     print *,'T=',tens(1,4),tens(2,4),tens(3,4)
c     do i=1,4
c     do j=1,4
c     print *,'T(',i,j,')=',tens(i,j)
c     end do
c     end do

      end

c***********************************************************************
      subroutine landauframe(u,tens)
c...PRC82,044904(2010)
      implicit double precision(a-h, o-z)
      parameter(eps=1d-6)
      dimension tens(4,4),u(4)
      dimension v(3),v1(3),v2(3)
      dimension gg(4,4)
      data gg/ 1,1,1,-1,  1,1,1,-1, 1,1,1,-1, -1,-1,-1,1/

c...Initial condition taken from Eckart's definition.

      u1s=u(1)
      u2s=u(2)
      u3s=u(3)
      u4s=u(4)
      v(1)=u(1)/u(4)
      v(2)=u(2)/u(4)
      v(3)=u(3)/u(4)
      vx0=v(1)
      vy0=v(2)
      vz0=v(3)

      do it=1,1000

c     e=0d0
c     do i=1,4
c     do j=1,4
c      e = e + tens(i,j)*gg(i,j)*u(i)*u(j)
c     end do
c     end do


c     u1(1)=tens(1,4)*u(4)-tens(1,1)*u(1)-tens(1,2)*u(2)-tens(1,3)*u(3)
c     u1(2)=tens(2,4)*u(4)-tens(2,1)*u(1)-tens(2,2)*u(2)-tens(2,3)*u(3)
c     u1(3)=tens(3,4)*u(4)-tens(3,1)*u(1)-tens(3,2)*u(2)-tens(3,3)*u(3)
c     u1(4)=tens(4,4)*u(4)-tens(4,1)*u(1)-tens(4,2)*u(2)-tens(4,3)*u(3)

c      u1(1)=u1(1)/e
c      u1(2)=u1(2)/e
c      u1(3)=u1(3)/e
c      u1(4)=u1(4)/e

c      if(abs(u1(1)-u(1)).le.1d-4.and.
c    &    abs(u1(2)-u(2)).le.1d-4.and.
c    &    abs(u1(3)-u(3)).le.1d-4)then
c    &    abs(u1(4)-u(4)).le.1d-4) then
c         go to 1000
c      endif
c     print *,'it=',it
c     print *,'ux=',u1(1),u(1),u1(1)-u(1)
c     print *,'uy=',u1(2),u(2),u1(2)-u(2)
c     print *,'uz=',u1(3),u(3),u1(3)-u(3)
c     print *,' e=',e
c     print *,'ue=',u1(4),u(4),u1(4)-u(4)

c      do k=1,3
c        u2(k)=u(k)
c        u(k)=u1(k)
c      end do

      e=tens(4,4)-tens(4,1)*v(1)-tens(4,2)*v(2)-tens(4,3)*v(3)
      v1(1)=tens(1,4)-tens(1,1)*v(1)-tens(1,2)*v(2)-tens(1,3)*v(3)
      v1(2)=tens(2,4)-tens(2,1)*v(1)-tens(2,2)*v(2)-tens(2,3)*v(3)
      v1(3)=tens(3,4)-tens(3,1)*v(1)-tens(3,2)*v(2)-tens(3,3)*v(3)
      if(abs(e).lt.1d-8) goto 2000
c     v1(1)=v1(1)/e
c     v1(2)=v1(2)/e
c     v1(3)=v1(3)/e

       if(abs(v1(1)-v(1)).le.eps.and.
     &    abs(v1(2)-v(2)).le.eps.and.
     &    abs(v1(3)-v(3)).le.eps)then
        goto 1000
       endif

       do k=1,3
         v2(k)=v(k)
         v(k)=v1(k)
       end do

      end do

      print *,'does not converge u1 u=',v2(1),v(1),v2(1)-v(1)
      print *,'u1y uy=',v2(2),v(2),v2(2)-v(2)
      print *,'u1z uz=',v2(3),v(3),v2(3)-v(3)
      print *,'u1s',u1s,u2s,u3s,u4s
c     call jamlist(5)
      goto 2000

1000  continue
      g=1.0d0/sqrt(1d0-(v(1)**2+v(2)**2+v(3)**2))
      u(1)=v(1)*g
      u(2)=v(2)*g
      u(3)=v(3)*g
      u(4)=g
      return

2000  continue
      u(1)=u1s
      u(2)=u2s
      u(3)=u3s
      u(4)=u4s

c      print *,'it=',it
c      print *,vx0,vy0,vz0
c      print *,v2(1),v2(2),v2(3)
c      print *,v1(1),v1(2),v1(3)

      end

c***********************************************************************
c...EOS
      function geteos(e,b)
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      if(mstc(50).eq.1) then  ! EOS-Q
        geteos=eosq(e)
      else if(mstc(50).eq.2) then  ! EOS-Q modified
        geteos=eosqmod(e)

      else if(mstc(50).eq.3) then  ! lattice EoS
        geteos=s95pv1(e)

      else if(mstc(50).eq.4) then  ! lattice EoS
        geteos=s95pv11(e)

      else if(mstc(50).eq.5) then  ! super soft EoS
        geteos=eoss(e)

      else if(mstc(50).eq.6) then  ! hadron gas
        geteos=0.15*e

      else if(mstc(50).eq.7) then
        a1=0.15
        e0=2.0
        a2=0.17
        a2=0.1
        geteos=a1*e
        if(e.ge.e0) geteos = a2*e + (a1-a2)*e0
      else if(mstc(50).eq.11) then
        geteos=getP(b,e)
      else if(mstc(50).ge.21.and.mstc(50).le.25) then
        fac = 146.51751415742*1d-3
        b0 = 0.15891d0 
        geteos=press(e/fac,b/b0)*fac
      else  
        print *,'geteos::invalid option mstc(50)=',mstc(50)
        stop
      endif

      end

c***********************************************************************
c...EOS super soft
      function eoss(e)
      implicit double precision(a-h, o-z)

      hc=0.1972
      eh=0.45   ! GeV/fm^3
      B4=0.23   ! GeV
      B=B4**4/hc**3
      eq=0.15*3*eh + 4*B
      pc=0.15*eh

      if(e .le. eh) then
         eoss = 0.15*e
      else
         eoss = pc
      endif

      end

c***********************************************************************
c...EOS Q
      function eosq(e)
      implicit double precision(a-h, o-z)

      hc=0.1972
      eh=0.45   ! GeV/fm^3
      B4=0.23   ! GeV
      B=B4**4/hc**3
      eq=0.15*3*eh + 4*B
      pc=0.15*eh

      if(e .le. eh) then
         eosq = 0.15*e
      else if( e .le.eq) then
         eosq = pc
      else
         eosq = (e - 4*B)/3.0
      endif

      end

c***********************************************************************
c...EOS Qmod
      function eosqmod(e)
      implicit double precision(a-h, o-z)

      a=1.0/3.5

      hc=0.1972
      eh=0.45   ! GeV/fm^3
      B4=0.23   ! GeV
      B=B4**4/hc**3
      eq=0.15*eh/a + (1+a)/a*B
      pc=0.15*eh

      if(e .le. eh) then
         eosqmod = 0.15*e
      else if( e .le.eq) then
         eosqmod = pc
      else
         eosqmod = a*e - (a+1)*B
      endif

      end

c***********************************************************************
      
      function s95pv1(eden)
      implicit double precision(a-h, o-z)
      parameter(nn=356)
      dimension e(nn),p(nn)

      data (e(i),i=1,nn)/
     & 0.000389186,0.000432619,0.000479651,0.000530493,0.000585368,
     & 0.000644507,0.000708154,0.000776564,0.000850005,0.000928757,
     & 0.00101312,0.00110339,0.00119991,0.00130302,0.00141308,
     & 0.00153046,0.00165557,0.00178883,0.00193068,0.00208159,
     & 0.00224206,0.0024126,0.00259377,0.00278614,0.00299032,
     & 0.00320695,0.00343673,0.00368036,0.0039386,0.00421225,
     & 0.00450216,0.00480923,0.00513438,0.00547863,0.00584301,
     & 0.00622864,0.00663669,0.00706839,0.00752506,0.00800808,
     & 0.00851889,0.00905904,
     & 0.00963016,0.010234,0.0108722,0.0115469,0.01226,0.0130136,
     & 0.0138099,0.0146514,0.0155405,0.0164799,0.0174723,0.0185206,
     & 0.0196278,0.0207973,0.0220324,0.0233367,0.0247139,0.0261681,
     & 0.0277033,0.0293239,0.0310346,0.0328401,0.0347455,0.0367562,
     & 0.0388776,0.0411157,0.0434766,0.0459667,0.0485926,0.0513616,
     & 0.0542808,0.0573581,0.0606015,0.0640194,0.0676206,0.0714144,
     & 0.0754103,0.0796183,0.084049,0.0887132,0.0936223,0.098788,
     & 0.104223,0.10994,0.115952,0.122273,0.128918,0.135901,
     & 0.143238,0.150946,0.159042,0.167542,0.176464,0.185828,
     & 0.195654,0.20596,0.216768,0.2281,0.239978,0.252425,
     & 0.265465,0.279123,0.293423,0.308394,0.324061,0.340452,
     & 0.357598,0.375527,0.394271,0.413861,0.43433,0.455712,
     & 0.478041,0.501353,0.525685,0.551074,0.57756,0.605181,
     & 0.63398,0.663999,0.69528,0.727868,0.761808,0.797148,
     & 0.833935,0.872219,0.912049,0.953478,0.996558,1.04134,
     & 1.08789,1.13626,1.1865,1.23858,1.29232,1.34751,
     & 1.404,1.46165,1.52037,1.58006,1.64065,1.70208,
     & 1.7643,1.82728,1.89099,1.95541,2.02052,2.08632,
     & 2.1528,2.21996,2.2878,2.35633,2.42554,2.49546,
     & 2.56608,2.63742,2.7095,2.78232,2.8559,2.93025,
     & 3.00539,3.08133,3.15809,3.23568,3.31412,3.39342,
     & 3.47361,3.55469,3.63669,3.71961,3.80347,3.8883,
     & 3.9741,4.06089,4.14869,4.23752,4.32738,4.41829,
     & 4.51028,4.60334,4.69751,4.79279,4.8892,4.98675,
     & 5.08547,5.18536,5.28643,5.38871,5.49221,5.59694,
     & 5.70292,5.81015,5.91867,6.02847,6.13958,6.252,
     & 6.36576,6.48087,6.59733,6.71517,6.8344,6.95503,
     & 7.07708,7.20056,7.32549,7.45187,7.57973,7.70907,
     & 7.83991,7.97227,8.10616,8.24158,8.37857,8.51712,
     & 8.65726,8.79899,8.94234,9.08731,9.23392,9.38218,
     & 9.53211,9.68372,9.83703,9.99204,10.1488,10.3073,
     & 10.4675,10.6295,10.7932,10.9588,11.1262,11.2954,
     & 11.4664,11.6393,11.814,11.9906,12.1691,12.3495,
     & 12.5319,12.7161,12.9023,13.0905,13.2807,13.4728,
     & 13.667,13.8632,14.0614,14.2617,14.464,14.6684,
     & 14.8749,15.0835,15.2942,15.5071,15.7221,15.9393,
     & 16.1587,16.3803,16.6041,16.8301,17.0584,17.2889,
     & 17.5217,17.7568,17.9941,18.2338,18.4759,18.7202,
     & 18.967,19.2161,19.4676,19.7216,19.9779,20.2367,
     & 20.498,20.7617,21.0279,21.2967,21.5679,21.8417,
     & 22.118,22.3969,22.6784,22.9625,23.2492,23.5385,
     & 23.8305,24.1251,24.4224,24.7224,25.0251,25.3306,
     & 25.6388,25.9497,26.2635,26.58,26.8994,27.2215,
     & 27.5465,27.8744,28.2052,28.5388,28.8754,29.2148,
     & 29.5573,29.9027,30.251,30.6024,30.9568,31.3142,
     & 31.6746,32.0381,32.4047,32.7744,33.1472,33.5232,
     & 33.9023,34.2845,34.67,35.0586,35.4505,35.8456,
     & 36.244,36.6456,37.0505,37.4588,37.8703,38.2852,
     & 38.7035,39.1252,39.5502,39.9787,40.4106,40.846,
     & 41.2848,41.7271,42.173,42.6223,43.0752,43.5317,
     & 43.9917,44.4554,44.9227,45.3936,45.8681,46.3464,
     & 46.8283,47.314,47.8034,48.2965,48.7934,49.2941,
     & 49.7986,50.3069/

      data (p(i),i=1,nn)/
     & 8.0974d-05,9.08028d-05,0.000101518,0.000113174,0.000125827,
     & 0.000139534,0.000154358,0.00017036,0.000187605,0.000206162,
     & 0.000226101,0.000247494,0.000270418,0.000294951,0.000321175,
     & 0.000349175,0.00037904,0.000410861,0.000444735,0.00048076,
     & 0.000519041,0.000559686,0.000602806,0.000648519,0.000696948,
     & 0.000748218,0.000802464,0.000859823,0.00092044,0.000984465,
     & 0.00105206,0.00112338,0.0011986,0.0012779,0.00136147,0.0014495,
     & 0.0015422,0.00163978,0.00174247,0.00185049,0.00196409,0.00208353,
     & 0.00220907,0.002341,0.0024796,0.00262518,0.00277807,0.00293859,
     & 0.0031071,0.00328396,0.00346956,0.0036643,0.0038686,0.0040829,
     & 0.00430766,0.00454336,0.00479051,0.00504963,0.00532128,
     & 0.00560604,0.0059045,0.0062173,0.00654511,0.0068886,
     & 0.00724852,0.0076256,
     & 0.00802064,0.00843447,0.00886795,0.00932197,0.00979749,0.0102955,
     & 0.010817,0.011363,0.0119348,0.0125334,0.01316,0.013816,
     & 0.0145027,0.0152214,0.0159736,0.0167607,0.0175844,0.0184462,
     & 0.0193479,0.0202911,0.0212778,0.0223099,0.0233892,0.0245179,
     & 0.0256981,0.026932,0.028222,0.0295704,0.0309796,0.0324524,
     & 0.0339913,0.0355992,0.0372788,0.0390332,0.0408655,0.0427789,
     & 0.0447767,0.0468622,0.0490391,0.051311,0.0536817,0.0561552,
     & 0.0587354,0.0614267,0.0642332,0.0671595,0.0702102,0.07339,
     & 0.0767039,0.0801569,0.0837542,0.0875013,0.0914037,0.0954671,
     & 0.0996974,0.104101,0.108683,0.113451,0.118412,0.123571,
     & 0.128936,0.134515,0.140314,0.146342,0.152606,0.159113,
     & 0.165873,0.172894,0.180184,0.187752,0.195606,0.203753,
     & 0.212199,0.220948,0.230005,0.239374,0.249058,0.25906,
     & 0.269382,0.280026,0.290996,0.302292,0.313917,0.325873,
     & 0.338161,0.350783,0.363741,0.377036,0.390671,0.404647,
     & 0.418966,0.433629,0.448638,0.463996,0.479704,0.495764,
     & 0.512179,0.52895,0.546079,0.563569,0.581421,0.599639,
     & 0.618225,0.63718,0.656507,0.67621,0.696289,0.716749,
     & 0.737591,0.758819,0.780434,0.80244,0.82484,0.847636,
     & 0.870831,0.894429,0.918432,0.942843,0.967666,0.992903,
     & 1.01856,1.04463,1.07113,1.09806,1.12542,1.15321,
     & 1.18144,1.2101,1.23922,1.26878,1.29879,1.32925,
     & 1.36017,1.39156,1.42341,1.45573,1.48852,1.52178,
     & 1.55553,1.58976,1.62447,1.65968,1.69538,1.73158,
     & 1.76828,1.80549,1.84321,1.88144,1.92019,1.95946,
     & 1.99926,2.03959,2.08045,2.12185,2.16379,2.20628,
     & 2.24932,2.29291,2.33706,2.38177,2.42705,2.4729,
     & 2.51933,2.56634,2.61392,2.6621,2.71087,2.76024,
     & 2.81021,2.86078,2.91196,2.96376,3.01618,3.06921,
     & 3.12288,3.17718,3.23211,3.28769,3.34391,3.40079,
     & 3.45831,3.5165,3.57535,3.63487,3.69506,3.75593,
     & 3.81748,3.87972,3.94265,4.00628,4.07061,4.13565,
     & 4.20139,4.26786,4.33504,4.40295,4.47159,4.54096,
     & 4.61108,4.68194,4.75355,4.82591,4.89904,4.97293,
     & 5.04759,5.12302,5.19923,5.27623,5.35402,5.4326,
     & 5.51199,5.59218,5.67317,5.75499,5.83762,5.92108,
     & 6.00538,6.0905,6.17647,6.26329,6.35096,6.43948,
     & 6.52887,6.61913,6.71025,6.80226,6.89515,6.98893,
     & 7.0836,7.17918,7.27566,7.37305,7.47136,7.57058,
     & 7.67074,7.77183,7.87385,7.97682,8.08074,8.18562,
     & 8.29145,8.39825,8.50602,8.61477,8.72451,8.83523,
     & 8.94694,9.05966,9.17338,9.28811,9.40386,9.52063,
     & 9.63843,9.75727,9.87714,9.99806,10.12,10.2431,
     & 10.3672,10.4923,10.6186,10.7459,10.8743,11.0038,
     & 11.1344,11.2661,11.3989,11.5329,11.6679,11.8041,
     & 11.9414,12.0799,12.2195,12.3603,12.5022,12.6453,
     & 12.7895,12.935,13.0816,13.2294,13.3784,13.5286,
     & 13.6801,13.8327,13.9866,14.1417,14.298,14.4556,
     & 14.6144,14.7744/


      do i=1,nn
        if(e(i).ge.eden) then
          ii=i
          goto 100
        endif
      end do
      ii=nn
100   continue
      
      if(ii.ge.2) then
        s95pv1=(p(ii)-p(ii-1))*(eden-e(ii-1))/(e(ii)-e(ii-1))+p(ii-1)
      else
        s95pv1=p(ii)
      endif

      end

c***********************************************************************
      
      function s95pv11(eden)
      implicit double precision(a-h, o-z)
      parameter(nn=664)
      dimension e(nn),p(nn)

      data (e(i),i=1,nn)/
     & 0.000389675,0.000433249,0.000480455,0.000531511
     & ,0.000586647,0.000646101,
     & 0.000710127,0.00077899,0.000852969,0.000932356,
     & 0.00101746,0.00110861,
     & 0.00120614,0.00131042,0.00142183,0.00154076,
     & 0.00166765,0.00180293,
     & 0.00194707,0.00210057,0.00226395,0.00243777,0.0026226,0.00281906,
     & 0.00302779,0.00324948,0.00348484,0.00373464,0.00399968,0.0042808,
     & 0.0045789,0.00489491,0.00522983,0.0055847,0.00596064,0.0063588,
     &0.00678042,0.00722679,0.00769927,0.00819932,0.00872843,0.00928823,
     & 0.00988039,0.0105067,0.011169,0.0118693,0.0126097,0.0133923,
     & 0.0142196,0.0150937,0.0160175,0.0169935,0.0180245,0.0191136,
     & 0.0202639,0.0214786,0.0227612,0.0241154,0.025545,0.0270539,
     & 0.0286464,0.0303269,0.0321,0.0339707,0.0359439,0.0380251,
     & 0.0402198,0.0425339,0.0449735,0.0475451,0.0502554,0.0531115,
     & 0.0561207,0.0592907,0.0626295,0.0661455,0.0698476,0.0737449,
     & 0.0778468,0.0821635,0.0867052,0.0914829,0.0965077,0.101791,
     & 0.107346,0.113185,0.119321,0.125768,0.13254,0.139652,
     & 0.147119,0.154958,0.163185,0.171817,0.180873,0.190371,
     & 0.200329,0.210769,0.22171,0.233174,0.245183,0.25776,
     & 0.270928,0.284712,0.299137,0.314229,0.330016,0.346523,
     & 0.363782,0.38182,0.400669,0.420359,0.440923,0.462395,
     & 0.484808,0.508199,0.532602,0.558057,0.5846,0.612273,
     & 0.641114,0.671167,0.702473,0.735077,0.769024,0.80436,
     & 0.841134,0.879393,0.919188,0.960569,1.00359,1.0483,
     & 1.09477,1.14303,1.19316,1.24521,1.29912,1.35463,
     & 1.41152,1.46963,1.52881,1.58896,1.64997,1.71179,
     & 1.77437,1.83765,1.90162,1.96625,2.03153,2.09746,
     & 2.16402,2.23123,2.29908,2.36759,2.43676,2.5066,
     & 2.57713,2.64837,2.72032,2.793,2.86643,2.94062,
     & 3.01559,3.09136,3.16793,3.24534,3.3236,3.40271,
     & 3.48271,3.5636,3.64541,3.72814,3.81183,3.89647,
     & 3.98209,4.0687,4.15633,4.24498,4.33467,4.42542,
     & 4.51724,4.61015,4.70417,4.7993,4.89556,4.99297,
     & 5.09155,5.1913,5.29225,5.3944,5.49777,5.60238,
     & 5.70824,5.81537,5.92377,6.03347,6.14447,6.2568,
     & 6.37046,6.48547,6.60185,6.7196,6.83874,6.95929,
     & 7.08126,7.20467,7.32952,7.45583,7.58362,7.71289,
     & 7.84367,7.97597,8.10979,8.24516,8.38209,8.52059,
     & 8.66068,8.80236,8.94566,9.09058,9.23715,9.38537,
     & 9.53526,9.68683,9.8401,9.99508,10.1518,10.3102,
     & 10.4704,10.6324,10.7961,10.9616,11.129,11.2982,
     & 11.4692,11.642,11.8167,11.9933,12.1718,12.3522,
     & 12.5345,12.7188,12.905,13.0931,13.2833,13.4754,
     & 13.6696,13.8657,14.0639,14.2642,14.4665,14.6709,
     & 14.8774,15.086,15.2967,15.5096,15.7246,15.9418,
     & 16.1612,16.3828,16.6065,16.8326,17.0608,17.2913,
     & 17.5241,17.7592,17.9966,18.2363,18.4783,18.7227,
     & 18.9694,19.2186,19.4701,19.724,19.9804,20.2392,
     & 20.5004,20.7642,21.0304,21.2991,21.5704,21.8442,
     & 22.1205,22.3994,22.6809,22.965,23.2517,23.541,
     & 23.833,24.1276,24.4249,24.725,25.0277,25.3331,
     & 25.6414,25.9523,26.2661,26.5826,26.902,27.2241,
     & 27.5491,27.877,28.2078,28.5414,28.878,29.2175,
     & 29.5599,29.9053,30.2537,30.6051,30.9595,31.3169,
     & 31.6773,32.0409,32.4075,32.7772,33.15,33.5259,
     & 33.905,34.2873,34.6728,35.0614,35.4533,35.8484,
     & 36.2468,36.6484,37.0534,37.4616,37.8732,38.2881,
     & 38.7064,39.1281,39.5531,39.9816,40.4135,40.8489,
     & 41.2878,41.7301,42.1759,42.6253,43.0782,43.5347,
     & 43.9948,44.4584,44.9257,45.3966,45.8712,46.3495,
     & 46.8314,47.3171,47.8065,48.2996,48.7965,49.2972,
     & 49.8017,50.3101,50.8223,51.3383,51.8582,52.3821,
     & 52.9099,53.4416,53.9772,54.5169,55.0605,55.6082,
     & 56.1599,56.7157,57.2756,57.8395,58.4076,58.9798,
     & 59.5562,60.1368,60.7215,61.3105,61.9037,62.5012,
     & 63.103,63.709,64.3194,64.9341,65.5532,66.1767,
     & 66.8046,67.4369,68.0736,68.7148,69.3605,70.0107,
     & 70.6654,71.3247,71.9885,72.6569,73.3299,74.0076,
     & 74.6899,75.3769,76.0686,76.765,77.4661,78.172,
     & 78.8827,79.5982,80.3185,81.0436,81.7736,82.5085,
     & 83.2483,83.9931,84.7427,85.4974,86.2571,87.0217,
     & 87.7914,88.5662,89.3461,90.131,90.9211,91.7164,
     & 92.5168,93.3224,94.1332,94.9492,95.7706,96.5972,
     & 97.4291,98.2663,99.1089,99.9569,100.81,101.669,
     & 102.533,103.403,104.278,105.159,106.045,106.937,
     & 107.834,108.737,109.646,110.56,111.48,112.405,
     & 113.337,114.274,115.216,116.165,117.119,118.08,
     & 119.046,120.018,120.996,121.98,122.969,123.965,
     & 124.967,125.975,126.989,128.008,129.034,130.067,
     & 131.105,132.149,133.2,134.257,135.32,136.389,
     & 137.465,138.547,139.635,140.73,141.831,142.939,
     & 144.052,145.173,146.3,147.433,148.573,149.72,
     & 150.873,152.032,153.199,154.372,155.552,156.738,
     & 157.931,159.131,160.338,161.551,162.772,163.999,
     & 165.233,166.474,167.722,168.977,170.239,171.508,
     & 172.784,174.067,175.357,176.655,177.959,179.271,
     & 180.59,181.916,183.249,184.59,185.938,187.293,
     & 188.656,190.026,191.403,192.788,194.181,195.58,
     & 196.988,198.403,199.825,201.255,202.693,204.138,
     & 205.591,207.052,208.521,209.997,211.481,212.973,
     & 214.473,215.98,217.496,219.019,220.551,222.09,
     & 223.638,225.193,226.756,228.328,229.908,231.496,
     & 233.092,234.696,236.308,237.929,239.558,241.196,
     & 242.841,244.495,246.158,247.829,249.508,251.196,
     & 252.892,254.597,256.311,258.033,259.763,261.503,
     & 263.251,265.008,266.773,268.547,270.33,272.122,
     & 273.923,275.733,277.551,279.379,281.215,283.061,
     & 284.915,286.779,288.652,290.534,292.425,294.325,
     & 296.234,298.153,300.081,302.018,303.964,305.92,
     & 307.885,309.86,311.844,313.838,315.841,317.854,
     & 319.876,321.908,323.949,326,328.061,330.132,
     & 332.212,334.302,336.402,338.512,340.631,342.761,
     & 344.9,347.05,349.209,351.379,353.558,355.748,
     & 357.947,360.157,362.377,364.608,366.848,369.099,
     & 371.36,373.631,375.913,378.205,380.508,382.821,
     & 385.145,387.479,389.824,392.179,394.545,396.921,
     & 399.309,401.707,404.115,406.535,408.965,411.406,
     & 413.858,416.321,418.795,421.28,423.776,426.283,
     & 428.801,431.33,433.87,436.421,438.984,441.558,
     & 444.143,446.739,449.347,451.966,454.597,457.238,
     & 459.892,462.557,465.233,467.921,470.621,473.332,
     & 476.055,478.789,481.536,484.294,487.064,489.845,
     & 492.639,495.444,498.262,501.091/

      data (p(i),i=1,nn)/
     & 8.10124e-05,9.08532e-05,0.000101583,0.000113258,0.000125934,
     & 0.00013967,0.000154528,0.000170572,0.000187869,0.000206487,
     & 0.000226498,0.000247977,
     & 0.000271003,0.000295655,0.000322018,0.00035018,
     & 0.000380232,0.000412269,
     & 0.000446391,0.0004827,0.000521305,0.000562317,
     & 0.000605854,0.000652038,
     & 0.000700996,0.000752861,0.000807772,0.000865874,
     & 0.000927318,0.000992261,
     &0.00106087,0.00113331,0.00120978,0.00129044,0.00137551,0.00146518,
     &0.00155968,0.00165921,0.00176403,0.00187437,0.00199049,0.00211266,
     & 0.00224116,0.00237627,0.00251832,0.0026676,0.00282447,0.00298927,
     & 0.00316236,0.00334413,0.00353498,0.00373532,0.0039456,0.00416626,
     & 0.0043978,0.0046407,0.00489549,0.00516271,0.00544293,0.00573676,
     & 0.00604481,0.00636773,0.00670621,0.00706094,0.00743268,0.0078222,
     & 0.0082303,0.00865783,0.00910567,0.00957474,0.010066,0.0105804,
     & 0.0111191,0.0116831,0.0122735,0.0128916,0.0135386,0.0142156,
     & 0.0149242,0.0156657,0.0164414,0.017253,0.0181021,0.0189901,
     & 0.0199189,0.0208902,0.0219058,0.0229677,0.0240778,0.0252382,
     & 0.0264511,0.0277185,0.029043,0.0304268,0.0318724,0.0333825,
     & 0.0349596,0.0366065,0.0383262,0.0401216,0.0419956,0.0439517,
     & 0.0459929,0.0481228,0.0503449,0.0526628,0.0550803,0.0576013,
     & 0.0602297,0.0629699,0.065826,0.0688026,0.0719041,0.0751354,
     & 0.0785012,0.0820067,0.085657,0.0894574,0.0934135,0.097531,
     & 0.101816,0.106273,0.110911,0.115733,0.120749,0.125963,
     & 0.131383,0.137017,0.142871,0.148953,0.15527,0.161831,
     & 0.168645,0.175718,0.18306,0.19068,0.198586,0.206787,
     & 0.215288,0.224095,0.233212,0.242643,0.252392,0.262461,
     & 0.272851,0.283567,0.294608,0.305978,0.317677,0.329708,
     & 0.342071,0.354769,0.367803,0.381174,0.394885,0.408936,
     & 0.42333,0.438067,0.453151,0.468582,0.484363,0.500495,
     & 0.516981,0.533822,0.551021,0.568579,0.5865,0.604784,
     & 0.623436,0.642456,0.661848,0.681614,0.701756,0.722278,
     & 0.743181,0.764469,0.786144,0.808208,0.830666,0.853519,
     & 0.876771,0.900425,0.924483,0.948949,0.973825,0.999115,
     & 1.02482,1.05095,1.0775,1.10448,1.13189,1.15973,
     & 1.188,1.21672,1.24588,1.27549,1.30555,1.33606,
     & 1.36703,1.39846,1.43036,1.46272,1.49556,1.52887,
     & 1.56266,1.59693,1.63169,1.66694,1.70269,1.73893,
     & 1.77567,1.81293,1.85069,1.88896,1.92775,1.96707,
     & 2.00691,2.04728,2.08818,2.12962,2.1716,2.21413,
     & 2.25721,2.30084,2.34503,2.38978,2.4351,2.48099,
     & 2.52746,2.57451,2.62214,2.67035,2.71916,2.76857,
     & 2.81857,2.86919,2.92041,2.97225,3.0247,3.07778,
     & 3.13148,3.18582,3.24079,3.29641,3.35267,3.40958,
     & 3.46715,3.52537,3.58426,3.64382,3.70405,3.76496,
     & 3.82655,3.88882,3.95179,4.01546,4.07983,4.1449,
     & 4.21069,4.27719,4.34441,4.41236,4.48104,4.55045,
     & 4.6206,4.6915,4.76315,4.83555,4.90871,4.98264,
     & 5.05734,5.13281,5.20906,5.2861,5.36392,5.44255,
     & 5.52197,5.60219,5.68323,5.76508,5.84776,5.93126,
     & 6.01559,6.10075,6.18676,6.27361,6.36132,6.44988,
     & 6.53931,6.6296,6.72077,6.81281,6.90574,6.99956,
     & 7.09427,7.18989,7.2864,7.38383,7.48218,7.58145,
     & 7.68164,7.78277,7.88483,7.98784,8.0918,8.19671,
     & 8.30258,8.40942,8.51723,8.62602,8.73579,8.84655,
     & 8.95831,9.07106,9.18482,9.29959,9.41538,9.53219,
     & 9.65003,9.76891,9.88882,10.0098,10.1318,10.2549,
     & 10.379,10.5042,10.6305,10.7578,10.8863,11.0158,
     & 11.1465,11.2782,11.4111,11.545,11.6801,11.8164,
     & 11.9537,12.0922,12.2319,12.3727,12.5147,12.6578,
     & 12.8021,12.9476,13.0942,13.2421,13.3911,13.5414,
     & 13.6929,13.8455,13.9994,14.1546,14.3109,14.4685,
     & 14.6274,14.7875,14.9489,15.1115,15.2754,15.4406,
     & 15.6071,15.7748,15.9439,16.1142,16.2859,16.4589,
     & 16.6332,16.8089,16.9859,17.1642,17.3439,17.5249,
     & 17.7073,17.8911,18.0763,18.2628,18.4507,18.64,
     & 18.8308,19.0229,19.2165,19.4115,19.6079,19.8057,
     & 20.005,20.2058,20.408,20.6117,20.8168,21.0234,
     & 21.2315,21.4411,21.6522,21.8649,22.079,22.2946,
     & 22.5118,22.7305,22.9508,23.1726,23.3959,23.6208,
     & 23.8473,24.0754,24.3051,24.5363,24.7691,25.0036,
     & 25.2397,25.4773,25.7167,25.9576,26.2002,26.4444,
     & 26.6903,26.9379,27.1871,27.438,27.6906,27.9449,
     & 28.2009,28.4586,28.718,28.9792,29.2421,29.5067,
     & 29.773,30.0411,30.311,30.5826,30.8561,31.1312,
     & 31.4082,31.687,31.9676,32.25,32.5343,32.8203,
     & 33.1082,33.3979,33.6895,33.983,34.2783,34.5755,
     & 34.8746,35.1756,35.4784,35.7832,36.0899,36.3985,
     & 36.709,37.0215,37.3359,37.6523,37.9707,38.291,
     & 38.6133,38.9375,39.2638,39.5921,39.9224,40.2547,
     & 40.589,40.9253,41.2637,41.6042,41.9467,42.2913,
     & 42.6379,42.9867,43.3375,43.6904,44.0454,44.4026,
     & 44.7618,45.1232,45.4868,45.8525,46.2203,46.5903,
     & 46.9625,47.3369,47.7134,48.0922,48.4732,48.8563,
     & 49.2417,49.6294,50.0192,50.4114,50.8057,51.2024,
     & 51.6013,52.0025,52.406,52.8118,53.2199,53.6304,
     & 54.0431,54.4582,54.8756,55.2954,55.7176,56.1421,
     & 56.569,56.9982,57.4299,57.864,58.3005,58.7394,
     & 59.1808,59.6245,60.0708,60.5195,60.9706,61.4243,
     & 61.8804,62.339,62.8001,63.2637,63.7299,64.1986,
     & 64.6698,65.1435,65.6198,66.0987,66.5802,67.0642,
     & 67.5508,68.0401,68.5319,69.0264,69.5235,70.0232,
     & 70.5256,71.0306,71.5383,72.0487,72.5618,73.0775,
     & 73.596,74.1171,74.641,75.1676,75.697,76.2291,
     & 76.764,77.3016,77.842,78.3852,78.9312,79.48,
     & 80.0316,80.5861,81.1434,81.7035,82.2665,82.8323,
     & 83.401,83.9726,84.5471,85.1245,85.7048,86.288,
     & 86.8741,87.4632,88.0553,88.6503,89.2482,89.8492,
     & 90.4531,91.0601,91.67,92.283,92.8989,93.518,
     & 94.14,94.7652,95.3934,96.0246,96.659,97.2965,
     & 97.937,98.5807,99.2275,99.8775,100.531,101.187,
     & 101.846,102.509,103.175,103.844,104.516,105.191,
     & 105.87,106.552,107.237,107.925,108.617,109.312,
     & 110.01,110.712,111.416,112.124,112.836,113.551,
     & 114.269,114.99,115.715,116.444,117.175,117.91,
     & 118.649,119.391,120.136,120.885,121.637,122.393,
     & 123.152,123.915,124.681,125.451,126.224,127,
     & 127.781,128.564,129.352,130.143,130.937,131.735,
     & 132.537,133.342,134.151,134.964,135.78,136.6,
     & 137.424,138.251,139.082,139.917,140.755,141.597,
     & 142.443,143.292,144.146,145.003,145.863,146.728,
     & 147.597,148.469,149.345,150.225,151.108,151.996,
     & 152.887,153.783,154.682,155.585,156.492,157.403,
     & 158.318,159.237,160.16,161.086/

      do i=1,nn
        if(e(i).ge.eden) then
          ii=i
          goto 100
        endif
      end do
      ii=nn
100   continue
      
      if(ii.ge.2) then
        s95pv11=(p(ii)-p(ii-1))*(eden-e(ii-1))/(e(ii)-e(ii-1))+p(ii-1)
      else
        s95pv11=p(ii)
      endif

      end

c**************************************************************************
      double precision function getP(b,e)
      implicit double precision(a-h, o-z)
      include 'jameos.inc'

      ie=int((e-emin)/de)
      if(ie.ge.maxe) then
c       e=maxe-1
c       Bag = 0.39693 ! GeV/fm**3 = 235 MeV^{1/4}
        getP = (e - 4*Bagconst)/3.0
        return
      endif

      if(ie.lt.0) ie=0

      ib=int(b/dn)
      if(ib.lt.0) then
        ib=0
      else if(ib.ge.maxn) then
        ib=maxn-1
      endif

      if(ie.eq.maxe-1) then
        if(ib.eq.maxn-1) then
          getP=eostable(ib,ie)
        else
          val1=eostable(ib,ie)
          val2=eostable(ib+1,ie)
          getP=rlinter(ib*dn,(ib+1)*dn,val1,val2,b)
        endif
        return
      endif

      ! first interpolation in e
      y1=eostable(ib,ie)
      y2=eostable(ib,ie+1)
      x1=emin+ie*de
      x2=emin+(ie+1)*de
      val1=rlinter(x1,x2,y1,y2,e)
      if(ib.eq.maxn-1) then
        getP=val1
        return
      endif
      y1=eostable(ib+1,ie)
      y2=eostable(ib+1,ie+1)

      val2=rlinter(x1,x2,y1,y2,e)
      getP=rlinter(ib*dn,(ib+1)*dn,val1,val2,b)

      end

      function rlinter(x1,x2,y1,y2,x)
      implicit double precision(a-h, o-z)
      rlinter = (y2-y1)*(x-x1)/(x2-x1)+y1;
      end

c**************************************************************************
      subroutine geteos11(e,b)
      implicit double precision(a-h, o-z)
      include 'jameos.inc'

      ie=int((e-emin)/de)
      if(ie.ge.maxe) then
c       e=maxe-1
c       Bag = 0.39693 ! GeV/fm**3
        getP = (e - 4*Bagconst)/3.0
        return
      endif

      if(ie.lt.0) ie=0

      ib=int(b/dn)
      if(ib.lt.0) then
        ib=0
      else if(ib.ge.maxn) then
        ib=maxn-1
      endif

      if(ie.eq.maxe-1) then
        if(ib.eq.maxn-1) then
          pnow=eostable(ib,ie)
          tnow=eostabt(ib,ie)
          bmnow=eostabmu(ib,ie)
          smnow=eostabsmu(ib,ie)
          fnow=eostablam(ib,ie)
          snow=eostabs(ib,ie)
        else
          val1=eostable(ib,ie)
          val2=eostable(ib+1,ie)
          pnow=rlinter(ib*dn,(ib+1)*dn,val1,val2,b)

          val1=eostabt(ib,ie)
          val2=eostabt(ib+1,ie)
          tnow=rlinter(ib*dn,(ib+1)*dn,val1,val2,b)

          val1=eostabmu(ib,ie)
          val2=eostabmu(ib+1,ie)
          bmunow=rlinter(ib*dn,(ib+1)*dn,val1,val2,b)

c         val1=eostabsmu(ib,ie)
c         val2=eostabsmu(ib+1,ie)
c         smunow=rlinter(ib*dn,(ib+1)*dn,val1,val2,b)

        endif
        return
      endif

      ! first interpolation in e
      y1=eostable(ib,ie)
      y2=eostable(ib,ie+1)
      x1=emin+ie*de
      x2=emin+(ie+1)*de
      val1=rlinter(x1,x2,y1,y2,e)

      y1m=eostabmu(ib,ie)
      y2m=eostabmu(ib,ie+1)
      val1m=rlinter(x1,x2,y1m,y2m,e)

      y1t=eostabt(ib,ie)
      y2t=eostabt(ib,ie+1)
      val1t=rlinter(x1,x2,y1t,y2t,e)

      if(ib.eq.maxn-1) then
        pnow=val1
        bmnow=val1m
        tnow=val1t
        return
      endif

      y1=eostable(ib+1,ie)
      y2=eostable(ib+1,ie+1)
      val2=rlinter(x1,x2,y1,y2,e)
      pnow=rlinter(ib*dn,(ib+1)*dn,val1,val2,b)

      y1=eostabmu(ib+1,ie)
      y2=eostabmu(ib+1,ie+1)
      val2m=rlinter(x1,x2,y1,y2,e)
      bmnow=rlinter(ib*dn,(ib+1)*dn,val1m,val2m,b)

      y1=eostabt(ib+1,ie)
      y2=eostabt(ib+1,ie+1)
      val2t=rlinter(x1,x2,y1,y2,e)
      tnow=rlinter(ib*dn,(ib+1)*dn,val1t,val2t,b)

      end

c***************************************************************************
      subroutine readEOStable(eosfname)
      implicit none
      character  eosfname*(*),tmp*26
      character tmpl(3)*80
      real*8 bden,e,p,t,bmu,smu,xqgp,s,pdens,bag
      include 'jam2.inc'
      include 'jameos.inc'
      integer ie,im,mstc106
      integer maxe1,maxb1

c     Bagconst=parc(145)


      open(81,file=eosfname,status='old')
      read(81,'(a)')tmpl(1)
      read(81,'(a)')tmpl(2)
      read(81,'(a)')tmpl(3)

      pard(142)=0.0d0
c     read(81,*)pard(142)   ! potential parameter K
      read(81,*)mstc106,bag     ! potential parameter
      if(mstc106.ne.mstc(106)) then
        print *,'readEOStable inconsistent potential parameter'
        print *,' mstc(106)=',mstc(106),'mstc106=',mstc106
        stop
      endif

      Bagconst=bag

      read(tmpl(2)(14:17),*) maxe1
      read(tmpl(3)(14:17),*) maxb1

      de = emax/(maxe-1)
      dn = bmax/(maxn-1)


      do im=0,maxn-1
          do ie=1,maxe
            if(mstc(145).eq.2) then
              read(81,*) bden,e,p,t,bmu,smu,xqgp,s,pdens
            else
              read(81,*) bden,e,p,t,bmu,smu,xqgp,s
            endif
c           read(81,*) e,p,bden,t,bmu
            eostable(im,maxe-ie)=p
            eostabt(im,maxe-ie)=t
            eostabmu(im,maxe-ie)=bmu
            eostabsmu(im,maxe-ie)=smu
            eostablam(im,maxe-ie)=xqgp
            eostabs(im,maxe-ie)=s
            if(mstc(145).eq.2) eostabp(im,maxe-ie)=pdens
        end do
        read(81,'(a)')tmp
      end do
      close(81)

      end

c***********************************************************************
      double precision function getpden(e,b)
c...Input   e (1/fm^4), b (1/fm^3)
      implicit none
      include 'fluid.inc'
      include 'jameos.inc'
      real*8 e,b,t,e0,b0,ev,bv,getEoStab
      parameter(e0=0.14651751415742, b0=0.15891d0)

      if(eos_mode.eq.11.or.eos_mode.eq.12) then ! EoS-Q
        ev=e*hbc
        getpden=getEoStab(b,ev,eostabp)
      else if(eos_mode.ge.1.and.eos_mode.le.5) then
        print *,'getpden not implemented for eos_mode=',eos_mode
        print *,'you should not set mstc(145)=2'
        stop
c       ev=e*hbc/e0
c       bv=b/b0
c       getpden=entro(ev,bv)*b0
      else  
        print *,'getpden::invalid option eos_mode=',eos_mode
        stop
      endif

      end

c***********************************************************************
      double precision function getentropy(e,b)
c...Input   e (1/fm^4), b (1/fm^3)
      implicit none
      include 'fluid.inc'
      include 'jameos.inc'
      real*8 e,b,t,e0,b0,ev,bv,getEoStab
      real*8 mu,entro,getP,pre
      real*8 temp,chem,press
      parameter(e0=0.14651751415742, b0=0.15891d0)

      if(eos_mode.eq.11.or.eos_mode.eq.12) then ! EoS-Q
        ev=e*hbc
c       getentropy=getEoStab(b,ev,eostabs)

        mu=getEoStab(b,ev,eostabmu)
        t=getEoStab(b,ev,eostabt)
        pre=getP(b,ev)
        if(t.gt.0d0) then
          getentropy=(e + pre - mu*b)/t
        else
          getentropy=0d0
        endif

      else if(eos_mode.ge.1.and.eos_mode.le.5) then
        ev=e*hbc/e0
        bv=b/b0
        getentropy=entro(ev,bv)*b0
c       t=temp(ev,bv)*0.001
c       mu=3*chem(ev,bv)*0.001
c       pre=press(ev,bv)*e0
c       if(t.gt.0d0) then
c         getentropy=(e + pre - mu*b)/t
c       else
c         getentropy=0d0
c       endif

      else  
        print *,'getentropy::invalid option eos_mode=',eos_mode
        stop
      endif

      end

c***********************************************************************
      subroutine geteos2(e,b,t,mu,smu,s)
c...Input   e (1/fm^4), b (1/fm^3)
c...Output  t (1/fm),  mu (1/fm), s (1/fm^3)
      implicit none
      include 'fluid.inc'
      include 'jameos.inc'
      real*8 e,b,t,mu,smu,s,e0,b0,ev,bv,getEoStab
      real*8 temp,chem,entro,schem
      parameter(e0=0.14651751415742, b0=0.15891d0)

      if(eos_mode.eq.11.or.eos_mode.eq.12) then ! EoS-Q
        ev=e*hbc
        mu=getEoStab(b,ev,eostabmu)/hbc
        t=getEoStab(b,ev,eostabt)/hbc
        s=getEoStab(b,ev,eostabs)
        smu=getEoStab(b,ev,eostabsmu)/hbc
      else if(eos_mode.ge.1.and.eos_mode.le.5) then
        ev=e*hbc/e0
        bv=b/b0
        t=temp(ev,bv)*0.001/hbc
        mu=3*chem(ev,bv)*0.001/hbc
        smu=schem(ev,bv)*0.001/hbc
        smu=mu/3.0-smu
        s=entro(ev,bv)*b0
      else  
        print *,'geteos2::invalid option eos_mode=',eos_mode
        stop
      endif

      end

c***********************************************************************
c     subroutine geteos(e,b,pnow,tnow,bmnow,csnow)
      subroutine geteos3(e,b,pnow,csnow)
c...Compute EOS from the input of energy density and baryon density.
c....e: enegy density in 1/fm^4,  b: baryon density in 1/fm^3
      implicit none
      include 'fluid.inc'
      include 'jameos.inc'
      real*8 e,b,e0,b0,ev,bv,getP,getEoStab,pnow,csnow
      real*8 press,temp,chem,schem,lambda,entro,css
      parameter(e0=0.14651751415742, b0=0.15891d0)

c       smnow=0.0d0
c       fnow=0.0d0
c       snow=0.0d0
      if(eos_mode.eq.11.or.eos_mode.eq.12) then ! EoS-Q
        ev=e*hbc

         if(ev.ge.20d0) then
c          print *,'geteos3 energy density too high?',ev
         endif

c       pnow=getEoStab(b,ev,eostable)/hbc
        pnow=getP(b,ev)/hbc

c       bmnow=getEoStab(b,ev,eostabmu)/hbc
c       tnow=getEoStab(b,ev,eostabt)/hbc
        csnow=1.0/sqrt(3.0)
      else if(eos_mode.ge.1.and.eos_mode.le.5) then
        ev=e*hbc/e0
        bv=b/b0
        pnow=press(ev,bv)*e0/hbc
c       tnow=temp(ev,bv)*0.001/hbc
c       bmnow=3*chem(ev,bv)*0.001/hbc
c       smnow=bmnow/3-schem(ev,bv)*0.01/hbc
c       fnow=lambda(ev,bv)
c       snow=entro(ev,bv)*b0
c       csnow=css(ev,bv)
        csnow=1.0/sqrt(3.0)
c       print *,'e=',e,' n= ',b,'P=',pnow,' t=',tnow
      else  
        print *,'geteos3::invalid option eos_mode=',eos_mode
        stop
      endif

      end

c**************************************************************************
      double precision function getEoStab(b,e,eostab)
      implicit none
      include 'jameos.inc'
      integer ie,ib
      real*8 b,e,Bag,val1,val2,y1,y2,x1,x2
      real*8 eostab(0:maxn,0:maxe),rlinter

      getEoStab=0.0d0
      ie=int((e-emin)/de)
      if(ie.ge.maxe) then
        ie=maxe-1
      endif
      if(ie.lt.0) ie=0
      ib=int(b/dn)
      if(ib.lt.0) then
        ib=0
      else if(ib.ge.maxn) then
        ib=maxn-1
      endif
      if(ib.eq.0.and.ie.eq.0) return

      if(ie.eq.maxe-1) then
        if(ib.eq.maxn-1) then
          getEoStab=eostab(ib,ie)
        else
          val1=eostab(ib,ie)
          val2=eostab(ib+1,ie)
          getEoStab=max(0d0,rlinter(ib*dn,(ib+1)*dn,val1,val2,b))
        endif
        return
      endif

      ! first interpolation in e
      y1=eostab(ib,ie)
      y2=eostab(ib,ie+1)
      x1=emin+ie*de
      x2=emin+(ie+1)*de
      val1=rlinter(x1,x2,y1,y2,e)
      if(ib.eq.maxn-1) then
        getEoStab=max(0d0,val1)
        return
      endif
      y1=eostab(ib+1,ie)
      y2=eostab(ib+1,ie+1)

      val2=rlinter(x1,x2,y1,y2,e)
      getEoStab=max(0d0,rlinter(ib*dn,(ib+1)*dn,val1,val2,b))

      end

c**************************************************************************
      real*8 function eoslookup(e,b,eostab)
      implicit none
      include 'jameos.inc'
      integer i,j
      real*8 e,b,x2,y2,x1,y1
      real*8 eostab(0:maxn,0:maxe)

      i=min(max(0,int(b/dn)),maxn-1)
      j=min(max(0,int((e-emin)/de)),maxe-1)

      x2=(b-(i*dn))/dn
      y2=(e-(emin+j*de))/de
      x1=1d0-x2
      y1=1d0-y2

      eoslookup=x1*y1*eostab(i,j)
     &       +x2*y1*eostab(i+1,j)
     &       +x1*y2*eostab(i,j+1)
     &       +x2*y2*eostab(i+1,j+1)


      end

c**************************************************************************

      subroutine compute_epd(tf,mub,mus,eden,pre,pden,bden)
      implicit none
      include 'jam2.inc'
      real*8 eden,pre,pden,ed,pd,pe,bden
      integer kcmin,kcmax,kc,kf,ispin,ns,ibary,ch,nc,nbot,id
      integer istat
      real*8 pm,width,mu,mus,mub,tf,fac
      character*16 chap
      parameter(kcmin=132, kcmax=366)

      fac=1.0d0/(2*paru(1)**2*paru(3)**3)

      eden=0d0
      pre=0d0
      pden=0d0
      bden=0d0

      do kc=kcmin,kcmax
        kf=kchg(kc,4)
        if(kf.eq.130.or.kf.eq.310) cycle
        if(kc.ge.195.and.kc.le.251) cycle
        ispin=max(1,mod(kf,10))   ! spin 2J+1
        ns=kchg(kc,7)             ! strangeness
        ibary=kchg(kc,6)/3        ! baryon number times 3.
        pm=pmas(kc,1)             ! particle mass
        mu=ns*mus + ibary*mub
        istat=1
        if(ibary.eq.0) istat=-1

        call edensity(pm,tf,mu,istat,ed,pe,pd)
        eden=eden+ed*fac*ispin
        pre=pre+pe*fac*ispin
        pden=pden+pd*fac*ispin
        bden=bden+ibary*pd*fac*ispin

        if(kchg(kc,3).eq.1) then
          call edensity(pm,tf,-mu,istat,ed,pe,pd)
          eden=eden+ed*fac*ispin
          pre=pre+pe*fac*ispin
          pden=pden+pd*fac*ispin
          bden=bden-ibary*pd*fac*ispin
        endif

      end do

      end

c**************************************************************************
      subroutine edensity(pm,tf,mu,istat,eden,pre,pdensity)
      implicit none
      integer istat,i,j,nb,opt_integ
      double precision pm,tf,mu,stat,pden,dis,p,e,dw
      real*8 eden,pre,pdensity

      integer NGP
      parameter(NGP=38)
      real*8 xg(NGP),wg(NGP)
      common/gaussgrid/xg,wg
      real*8 xg1(NGP),wg1(NGP)
      common/gaussgrid1/xg1,wg1

c...Numerical integration.
      pdensity=0.0d0
      eden=0.0d0
      pre=0.0d0
      do i=1, NGP
        p=xg1(i)
        dw=wg1(i)
        e=sqrt(p*p+pm**2)
        dis=exp(-(e-mu)/tf)
        dis = dis/(1.0+istat*dis)
        pdensity = pdensity + p*p*dw*dis
        eden = eden + e*p*p*dw*dis
        pre = pre + p*p*p*p/(3*e)*dw*dis
      end do

      end

c**************************************************************************
      double precision function edminQGP(rhob)
c...Compute minimum energy density (at T=0) at the baryon density rhob with
c...massless QGP bag model.
      implicit none
      include 'jameos.inc'
      real*8 rhob,nq,pisq,hc,x
      parameter(nq=2.5, pisq=9.86960440108936,hc=0.1973269788)
c     x=max(0d0,rhob)
      x=abs(rhob)
      edminQGP=nq/(108*pisq)*(81*pisq/nq*x)**(4.0/3.0)*hc + Bagconst

      end

