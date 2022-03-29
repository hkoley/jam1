c**********************************************************************
c  This code is based on 1D PPM scheme for (3+1)Dimension
c  U(1) = g^2*(e+p)*v_x
c  U(2) = g^2*(e+p)*v_y
c  U(3) = g^2*(e+p)*v_z 
c  U(0) = g^2*(e+p)-p
c  U(4) = g*n_B
c...Original '98 6 4  BY TETSU
c...Rewritten by Y. Nara (Oct. 2017) but algorithm is exactly
c...the same as originally written by T.Hirano.
c T. Hirano and Y, Nara, Prog. Theor. Exp. Phys. 2012, 01A203
c DOI: 10.1093/ptep/pts007
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     subroutine ppm1d(tau,mxx,dl,dt,uu,vx,sv,pr,g,opt_timelike,iz)
      subroutine ppm1d(tau,mxx,dl,dt,uu,vx,sv,pr,fh,opt_timelike,iz)

c...tau: current time for tau-eta coordinate.
c...mxx: dimension of uu
c...dl:  = dx spacial size in fm
c...dt:  time step size
c...vx:  velocity: vx or vy or yz
c...vs:  sound velocity
c...pr:  pressure
c...g:   flux to update uu, i.e.  uu += g
c...iz:  =1: include source term for tau-eta coordinate.

      implicit none
      integer mxx,i,j,opt_timelike,iz
      real*8 dl,dt,tau
      real*8 uu(0:4,0:mxx),vx(0:mxx),sv(0:mxx),pr(0:mxx),g(0:4,0:mxx)
      real*8 dmu(0:4,0:mxx),ddu(0:4,0:mxx),s(0:4,0:mxx)
      real*8 ur(0:4,0:mxx),ul(0:4,0:mxx),uh(0:4)
      real*8 vr(0:mxx),svr(0:mxx),vl(0:mxx),svl(0:mxx)
      real*8 urh(0:4,0:mxx),ulh(0:4,0:mxx)
      real*8 fr(0:4,0:mxx),fl(0:4,0:mxx)
      real*8 difa,difb,du,difr,difl,mim,det,u2
      real*8 ene,nba,pre,sva,tplc,mulc
      real*8 vvh,svh,rv1,rv2,rv3,rv4,avrl,difrl
      real*8 fh(0:4,0:mxx),fe
      real*8 bplus,bminus,fg
      real*8 lam,va,adu,gam,vel
      parameter(mim = 0.00001d0)

      lam = dt/dl

c  HALF STEP
c  R: 0-->MAX-1
c  L: 1-->MAX      
c...Calculation of delta_m U#
c...Compute dmu and ddu
      do i=1,mxx-1
      do j=0,4
        difa = uu(j,i+1) - uu(j,i)
        difb = uu(j,i)    -uu(j,i-1)
        du = difa/2d0+ difb/2d0
        if(difa*difb .gt. 0d0) then
          adu=abs(du)
          dmu(j,i) = (min(adu,2d0*abs(difa),2d0*abs(difb)))*du/adu
        else 
          dmu(j,i) = 0d0
        endif
        ddu(j,i) = difa/6d0/dl/dl - difb/6d0/dl/dl
      end do
      end do

c...Average value at right and left state. Compute ur and ul
       ur(:,0) = 0.5*(uu(:,1)+uu(:,0))
       ur(:,1) = 0.5*(uu(:,2)+uu(:,1))-dmu(:,2)/6d0+dmu(:,1)/6d0
       ur(:,mxx-1) = 0.5*(uu(:,mxx)+uu(:,mxx-1))
       ur(:,mxx) = uu(:,mxx)

       ul(:,0) = uu(:,0)
       ul(:,1) = 0.5*(uu(:,1)+uu(:,0))
       ul(:,mxx-1) = 0.5*(uu(:,mxx-2)+uu(:,mxx-1))
     &  -dmu(:,mxx-1)/6d0+dmu(:,mxx-2)/6d0
       ul(:,mxx) = 0.5*(uu(:,mxx-1)+uu(:,mxx))

      do i = 2,mxx-2
        ur(:,i) = 0.5*(uu(:,i+1)+uu(:,i))
     &           - dmu(:,i+1)/6d0 + dmu(:,i)/6d0
        ul(:,i) = 0.5*(uu(:,i-1)+uu(:,i))
     &           - dmu(:,i)/6d0   + dmu(:,i-1)/6d0
      end do

c...Monotonicity condition
      do i = 0,mxx
      do j=0,4
        difr = ur(j,i) - uu(j,i)
        difl = uu(j,i) - ul(j,i) 
        difrl = difr+difl
        avrl = (ur(j,i)+ul(j,i))/2d0
c....Eq.(2.33)
        if(difr*difl .lt. 0d0) then
         ur(j,i) = uu(j,i)
         ul(j,i) = uu(j,i)
        endif
c...Eq.(2.34)
        if(difrl*(uu(j,i)-avrl).gt. difrl*difrl/6d0) then
          ul(j,i) = 3d0*uu(j,i)-2d0*ur(j,i)
         endif
c...Eq.(2.35)
        if(-difrl*difrl/6d0 .gt.difrl*(uu(j,i)-avrl)) then
         ur(j,i) = 3d0*uu(j,i)-2d0*ul(j,i)
        endif
      end do
      end do

c********************************************************************
c...Density MUST BE non-negative 
      do i = 0,mxx
        if((ur(0,i) .lt. 0d0).or. (ul(0,i) .lt. 0d0))then
         ur(0,i) = uu(0,i)
         ul(0,i) = uu(0,i)
        endif

c....Is this OK?  What about anti-baryon?
        if((ur(4,i) .lt.0d0) .or. (ul(4,i) .lt. 0d0))then
         ur(4,i) = uu(4,i)
         ul(4,i) = uu(4,i)
        endif

      end do

c...If a cell is in the vacuum,a value at cell faces is taken to be zero.
      do i = 0,mxx
        if(uu(0,i) .eq. 0d0) then
          ur(:,i) = 0d0
          ul(:,i) = 0d0
        endif
       end do

c....Determine v_r, v_l
c     call vrvl(vr,svr,vl,svl,ur,ul,tau,mxx)


c...Compute UL,UR and FL, FR  original subroutine avlr(vv,sr,...)
      do i = 0,mxx-1

        fr(:,i) = 0d0
        vr(i) = 0d0
        svr(i) = 0d0
        rv2=(vx(i)+sv(i))/(1d0+vx(i)*sv(i))
        bplus  = max(0d0,rv2)
c...Eq.(2.43)
        urh(:,i) = ur(:,i)-bplus*lam/2d0
     &   *(ur(:,i)-ul(:,i)-(1d0-2d0/3d0*bplus*lam)
     &   *6d0*(uu(:,i)-(ur(:,i)+ul(:,i))/2d0))

c     u2=sqrt(urh(1,i)**2+urh(2,i)**2+urh(3,i)**2)
c     if(urh(0,i).lt.u2) then
cc      urh(0,i)=u2+1d-1r
c       urh(1,i)=urh(0,i)*(1d0-1d-15)*urh(1,i)/u2
c       urh(2,i)=urh(0,i)*(1d0-1d-15)*urh(2,i)/u2
c       urh(3,i)=urh(0,i)*(1d0-1d-15)*urh(3,i)/u2
c     endif

      uh(:) = urh(:,i)/tau
      if (uh(0) .lt. 0d0) uh(0) = 0d0
      det = uh(0)*uh(0)-uh(1)*uh(1)-uh(2)*uh(2)-uh(3)*uh(3)
      if (det .gt. 0d0) then
         call thermal(uh,ene,nba,pre,sva,tplc,mulc,vel)
         uh=uh*tau
         fg = uh(0)+tau*pre
         if(fg .ne. 0d0) then
           fr(1,i) = uh(1)*uh(1)/fg + tau*pre
           fr(2,i) = uh(1)*uh(2)/fg
           fr(3,i) = uh(1)*uh(3)/fg
           fr(0,i) = uh(1)
           fr(4,i) = uh(1)*uh(4)/fg
           vr(i)   = uh(1)/fg
           svr(i)  = sva
         endif
      endif

        fl(:,i+1) = 0d0
        vl(i+1) = 0d0
        svl(i+1) = 0d0
        rv4=-(vx(i+1)-sv(i+1))/(1d0-vx(i+1)*sv(i+1))
        bminus = max(0d0,rv4)
c...Eq.(2.42)
        ulh(:,i+1)= ul(:,i+1)+bminus*lam/2d0
     & *(ur(:,i+1)-ul(:,i+1)+(1d0-2d0/3d0*bminus*lam)
     & *6d0*(uu(:,i+1)-(ur(:,i+1)+ul(:,i+1))/2d0))

c     u2=sqrt(ulh(1,i)**2+ulh(2,i)**2+ulh(3,i)**2)
c     if(ulh(0,i).lt.u2) then
cc      ulh(0,i)=u2+1d-12
c       ulh(1,i)=ulh(0,i)*(1d0-1d-15)*ulh(1,i)/u2
c       ulh(2,i)=ulh(0,i)*(1d0-1d-15)*ulh(2,i)/u2
c       ulh(3,i)=ulh(0,i)*(1d0-1d-15)*ulh(3,i)/u2
c     endif

      uh(:) = ulh(:,i+1)/tau
      if (uh(0) .lt. 0d0) uh(0) = 0d0
      det = uh(0)*uh(0)-uh(1)*uh(1)-uh(2)*uh(2)-uh(3)*uh(3)
      if(det.gt.0d0) then
         call thermal(uh,ene,nba,pre,sva,tplc,mulc,vel)
         uh=uh*tau
         fg = uh(0)+tau*pre
         if(fg .ne. 0d0) then
         fl(1,i+1) = uh(1)*uh(1)/fg + tau*pre
         fl(2,i+1) = uh(1)*uh(2)/fg
         fl(3,i+1) = uh(1)*uh(3)/fg
         fl(0,i+1) = uh(1)
         fl(4,i+1) = uh(1)*uh(4)/fg
         vl(i+1)   = uh(1)/fg
         svl(i+1)  = sva
         endif
      endif

      end do

c...Compute flux.
      do i = 0,mxx-1

        fh(:,i) = 0d0
        vvh = (vr(i)+vl(i+1))/2d0
        svh = (svr(i)+svl(i+1))/2d0
        rv1 = (vvh+svh)/(1d0+vvh*svh)
        rv2 = (vl(i+1)+svl(i+1))/(1d0+vl(i+1)*svl(i+1))
        rv3 = (vvh-svh)/(1d0-vvh*svh)
        rv4=(vr(i)-svr(i))/(1d0-vr(i)*svr(i))
        bplus = max(0d0,rv1,rv2)
        bminus = min(0d0,rv3,rv4)

c....Eq.(2.49)
        if(bplus-bminus.ne.0d0) then
          fh(:,i) = lam*( bplus*fr(:,i)-bminus*fl(:,i+1)
     &                +bplus*bminus*(ulh(:,i+1)-urh(:,i)) )
     &               /(bplus-bminus)
        endif

      end do

c...Source term in the case of tau-eta coordinate.
c     s=0d0
c     if(iz.ne.0) then
c       do i=0,mxx
c         fe=uu(0,i)+tau*pr(i)
c         if(fe.ne.0d0) then
c         s(0,i) = -uu(iz,i)*uu(iz,i)/fe/tau*dt - pr(i)*dt
c         s(0,i) = -uu(0,i)*vx(i)**2/tau*dt - pr(i)*(1d0+vx(i)**2)*dt
c         s(iz,i) = -uu(iz,i)/tau*dt
c         uu(0,i) = uu(0,i) -uu(iz,i)*uu(iz,i)/fe/tau*dt - pr(i)*dt
c         uu(iz,i) = uu(iz,i) -uu(iz,i)/tau*dt
c         endif
c       end do
c     end if

c...K. Murase method.
c     if(opt_timelike.eq.1) call limiter_timelike(mxx,fh,uu,vx,lam)

c     g(:,0)   = fh(:,0)-fh(:,1)
c     g(:,mxx) = fh(:,mxx-2)-fh(:,mxx-1)

      g(:,0)   = -fh(:,0)
      g(:,mxx) = fh(:,mxx-1)
      fh(:,mxx)=0d0

      do i=1,mxx-1
        g(:,i) = fh(:,i-1)-fh(:,i)
      end do


c...Advance the next time step. Eq.(2.50).
c     uu(:,0)   = uu(:,0)  -lam*(fh(:,1)-fh(:,0))
c     uu(:,mxx) = uu(:,mxx)-lam*(fh(:,mxx-1)-fh(:,mxx-2))
c     do i = 1,mxx-1
c       uu(:,i) = uu(:,i)-lam*(fh(:,i)-fh(:,i-1))
c     end do

CCC new terms in tau-eta coordinate 
CCC (devided by 3 due to operater splitting.)
c     if(opt_taueta.eq.1.and.is1 .eq. 2) then
c       uu(3,i) = uu(3,i)-uu(3,i)/tau*dt
c       if(fe(i) .ne. 0d0)then
c         uu(0,i) = uu(0,i)
c    &              -uu(3,i)*uu(3,i)/fe(i)/tau*dt-pressu(i)*dt
c       endif
c     endif


      end

c*******************************************************************
      subroutine limiter_timelike(mxx,f,uu,vx,lam)

c...Flux limiter to preserve timelikeness of energy-momentum tensor
c...by K.Murase method.
      implicit none
      integer mxx,i
      double precision f(0:4,0:mxx),uu(0:4,0:mxx),vx(0:mxx),lam,v1,v2
      double precision ur(0:4),ul(0:4),ur1(0:4),ul1(0:4),ur0,ul0
      double precision ur4,ul4,alpha,eps,eps1,eps3,scalef,fac
c     parameter(alpha=0.5d0,eps=0.02d0,eps3=1d-10)
      parameter(alpha=0.5d0,eps=0.00d0,eps3=1d-10)
      parameter(eps1=0d0)
      real*8 alp,alp1,alp2

c...Check.
      real*8 facsave(0:mxx)

      do i=0,mxx-1

c....Murase original.
c       ur(:)=(alpha-eps)*uu(:,i)
c       ul(:)=(1.0-alpha-eps)*uu(:,i+1)

c...Murase new.
        ur(:)=min(1.0d0,max(0.0d0,0.5+vx(i)*lam))*uu(:,i)
        ul(:)=min(1.0d0,max(0.0d0,0.5-vx(i+1)*lam))*uu(:,i+1)

c...Murase new II use the velocity of the total energy density.
c       v1=0.0; v2=0.0;
c       if(u(0,i).ne.0d0) v1=uu(1,i)/uu(0,i)
c       ur(:)=min(1.0d0,max(0.0d0,0.5+v1*lam))*uu(:,i)
c       if(u(0,i+1).ne.0d0) v2=uu(1,i+1)/uu(0,i+1)
c       ul(:)=min(1.0d0,max(0.0d0,0.5-v2*lam))*uu(:,i+1)

c       ur(:)=alp1*uu(:,i)
c       ul(:)=alp2*uu(:,i+1)


        ur0=ur(0)-f(0,i)
        ul0=ul(0)+f(0,i)
        if(ur0.lt.eps1) then
          f(:,i) = ur(:)
        else if(ul0.lt.eps1) then
          f(:,i) = -ul(:)
        endif

        ur1(:)=ur(:)-f(:,i)
        ul1(:)=ul(:)+f(:,i)
        ur4=ur1(0)**2-ur1(1)**2-ur1(2)**2-ur1(3)**2
        ul4=ul1(0)**2-ul1(1)**2-ul1(2)**2-ul1(3)**2
        if(ur4.ge.0d0.and.ul4.ge.0d0) cycle

        fac=max(0d0,min(scalef(f(:,i),-ur),scalef(f(:,i),ul))-eps3)
c       f(:,i) = fac*f(:,i)
        f(0,i) = fac*f(0,i)
        f(1,i) = fac*f(1,i)
        f(2,i) = fac*f(2,i)
        f(3,i) = fac*f(3,i)

c....Check.
        ur1(:)=ur(:)-f(:,i)
        ul1(:)=ul(:)+f(:,i)
        ur4=ur1(0)**2-ur1(1)**2-ur1(2)**2-ur1(3)**2
        ul4=ul1(0)**2-ul1(1)**2-ul1(2)**2-ul1(3)**2
        if(ur4.lt.0d0.or.ul4.lt.0d0) then
         print *,'ur<0 or ul<0?',ur4,ul4
        endif

c....Check.
        facsave(i)=fac
        if(i.gt.0) then
        ur1(:) = uu(:,i) + f(:,i-1)-f(:,i)
        ur4=ur1(0)**2-ur1(1)**2-ur1(2)**2-ur1(3)**2
        if(ur4.lt.0d0) then
         print *,'u<0 ?',i,ur4,ur1(0),ur1(1),ur1(2),ur1(3)
         print *,'fac=',facsave(i-1),fac
c        print *,'alp1 alp2=',alp1,alp2
        endif
        endif

      end do

      return

c...check
      do i=1,mxx-1
        ur1(:) = uu(:,i) + f(:,i-1)-f(:,i)
        ur4=ur1(0)**2-ur1(1)**2-ur1(2)**2-ur1(3)**2
        if(ur1(0).eq.0d0) cycle
c       if(ur1(0).le.1d-15.or.ur4.lt.0d0) then
        if(ur4.lt.0d0) then
         print *,'timelike: u<0 ?',i,ur4,ur1(0),ur1(1),ur1(2),ur1(3)
         print *,'u=',uu(0,i),uu(1,i),uu(2,i),uu(3,i),uu(4,i)
         print *,'f=',f(0,i-1),f(1,i-1),f(2,i-1),f(3,i-1),f(4,i-1)
         print *,'f2=',f(0,i),f(1,i),f(2,i),f(3,i),f(4,i)
        endif
      end do

      end

c*******************************************************************
      double precision function scalef(f,u)
      implicit none
      double precision f(0:4),u(0:4)
      double precision a,b,c,d4,s1,s2

      a=f(0)*f(0)-f(1)*f(1)-f(2)*f(2)-f(3)*f(3)
      b=f(0)*u(0)-f(1)*u(1)-f(2)*u(2)-f(3)*u(3)
      c=u(0)*u(0)-u(1)*u(1)-u(2)*u(2)-u(3)*u(3)

      if(a+2*b+c .ge.0d0) then
        scalef=1.0d0
      else if(abs(a).lt.1d-10) then
        if(abs(b).lt.1d-10) then
          scalef=0.0d0
        else
          scalef = c /(2*b)
        endif
      else
        b = b/a
        c = c/a
        d4=sqrt(max(0d0,b*b-c))
        s1=-b-d4
        s2=-b+d4
        if(s2 .le.1.0d0) then
          scalef=s2
        else
          scalef=s1
        endif
      endif

      end

c*******************************************************************
c...This routine is not used.
c....Determine v_r, v_l
      subroutine vrvl(vr,svr,vl,svl,ur,ul,tau,mxx)
      implicit none
      integer mxx,ix,io,j,il
      real*8 vr(0:mxx),svr(0:mxx),vl(0:mxx),svl(0:mxx)
      real*8 ur(0:4,0:mxx),ul(0:4,0:mxx)
      real*8 uh(0:4),tau,ene,nba,pre,sva,tplc,mulc,fg,det
      real*8 mim,vel
      parameter(mim = 0.00001d0)

      do ix = 0,mxx
        vr(ix) = 0d0
        svr(ix) = 0d0
        vl(ix) = 0d0
        svl(ix) = 0d0
      do il=1,2

      io=0
      if(il.eq.1) then
        if (ix .eq. mxx) io=-1
        uh(:) = ur(:,ix+io)
      else
        if (ix .eq. 0) io=1
        uh(:) = ul(:,ix+io)
      endif

      if (sqrt(uh(0)*uh(0)) .lt. mim) then
         uh(:)=0d0
      endif
      if (sqrt(uh(1)*uh(1)) .lt. mim) uh(1) = 0d0
      if (uh(0) .lt. 0d0) then
         uh(:)=0d0
      endif
      det = uh(0)*uh(0)-uh(1)*uh(1)-uh(2)*uh(2)-uh(3)*uh(3)
      if (det .gt. 0d0) then
        uh = uh/tau
        call thermal(uh,ene,nba,pre,sva,tplc,mulc,vel)
        fg = uh(0)+tau*pre
        if(fg. ne.0d0) then
          if(il.eq.1) then
            vr(ix) = uh(1)/fg
            svr(ix) = sva
          else
            vl(ix) = uh(1)/fg
            svl(ix) = sva
          endif
        endif
      endif

      end do
      end do

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine ppm1d2(tau,mxx,dl,dt,uu,vv,sv,ix,iy,iz,g)

      implicit none
      integer mxx,i,j,k,ix,iy,iz
      real*8 dl,dt,tau
      real*8 uu(0:4,0:mxx), vv(0:mxx),sv(0:mxx),g(0:4,0:mxx)
      real*8 dmu(0:4,0:mxx),ddu(0:4,0:mxx)
      real*8 ur(0:4,0:mxx),ul(0:4,0:mxx),uh(0:4)
      real*8 vr(0:mxx),svr(0:mxx),vl(0:mxx),svl(0:mxx)
      real*8 urh(0:4,0:mxx),ulh(0:4,0:mxx)
      real*8 fr(0:4,0:mxx),fl(0:4,0:mxx)
      real*8 difa,difb,du,difr,difl,mim,det
      real*8 ene,nba,pre,sva,tplc,mulc
      real*8 vvh,svh,rv1,rv2,rv3,rv4,avrl,difrl
      real*8 fh(0:4,0:mxx)
      real*8 bplus,bminus,fg
      real*8 lam,va,adu,gam,vel
      parameter(mim = 0.00001d0)

      lam = dt/dl

c  HALF STEP
c  R: 0-->MAX-1
c  L: 1-->MAX      
c...Calculation of delta_m U#
c...Compute dmu and ddu
      do i=1,mxx-1
      do j=0,4
        difa = uu(j,i+1) - uu(j,i)
        difb = uu(j,i)    -uu(j,i-1)
        du = 0.5d0*(difa + difb)
        if(difa*difb .gt. 0d0) then
          adu=abs(du)
          dmu(j,i) = (min(adu,2d0*abs(difa),2d0*abs(difb)))*du/adu
        else 
          dmu(j,i) = 0d0
        endif
        ddu(j,i) = difa/6d0/dl/dl - difb/6d0/dl/dl
      end do
      end do

c...Average value at right and left state. Compute ur and ul
       ur(:,0) = 0.5*(uu(:,1)+uu(:,0))
       ul(:,0) = uu(:,0)

       ur(:,1) = 0.5*(uu(:,2)+uu(:,1))-dmu(:,2)/6d0+dmu(:,1)/6d0
       ul(:,1) = 0.5*(uu(:,1)+uu(:,0))

       ur(:,mxx-1) = 0.5*(uu(:,mxx)+uu(:,mxx-1))
       ul(:,mxx-1) = 0.5*(uu(:,mxx-2)+uu(:,mxx-1))
     &  -dmu(:,mxx-1)/6d0+dmu(:,mxx-2)/6d0

       ur(:,mxx) = uu(:,mxx)
       ul(:,mxx) = 0.5*(uu(:,mxx-1)+uu(:,mxx))

      do i = 2,mxx-2
        ur(:,i) = 0.5*(uu(:,i+1)+uu(:,i))
     &           - dmu(:,i+1)/6d0 + dmu(:,i)/6d0
        ul(:,i) = 0.5*(uu(:,i-1)+uu(:,i))
     &           - dmu(:,i)/6d0   + dmu(:,i-1)/6d0
      end do

c...Monotonicity condition
      do i = 0,mxx
      do j=0,4
        difr = ur(j,i) - uu(j,i)
        difl = uu(j,i) - ul(j,i) 
        difrl = difr+difl
        avrl = (ur(j,i)+ul(j,i))/2d0
c....Eq.(2.33)
        if(difr*difl .lt. 0d0) then
         ur(j,i) = uu(j,i)
         ul(j,i) = uu(j,i)
        endif
c...Eq.(2.34)
        if(difrl*(uu(j,i)-avrl).gt. difrl*difrl/6d0) then
          ul(j,i) = 3d0*uu(j,i)-2d0*ur(j,i)
         endif
c...Eq.(2.35)
        if(-difrl*difrl/6d0 .gt.difrl*(uu(j,i)-avrl)) then
         ur(j,i) = 3d0*uu(j,i)-2d0*ul(j,i)
        endif
      end do
      end do

c*******************************************************************::
c...Density MUST BE non-negative 
      do i = 0,mxx
        if((ur(0,i) .lt. 0d0).or. (ul(0,i) .lt. 0d0))then
         ur(0,i) = uu(0,i)
         ul(0,i) = uu(0,i)
        endif

c....Is this OK?  What about anti-baryon?
        if((ur(4,i) .lt.0d0) .or. (ul(4,i) .lt. 0d0))then
         ur(4,i) = uu(4,i)
         ul(4,i) = uu(4,i)
        endif

      end do

c...If a cell is in the vacuum,a value at cell faces is taken to be zero.
      do i = 0,mxx
        if(uu(0,i) .eq. 0d0) then
          ur(:,i) = 0d0
          ul(:,i) = 0d0
        endif
       end do

c....Determine v_r, v_l
c     call vrvl(vr,svr,vl,svl,ur,ul,tau,mxx)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      subroutine avlr(vv,sv,vr,svr,vl,svl
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      k=ix + 2*iy + 3*iz

      do i = 0,mxx-1

        fr(:,i) = 0d0
        fl(:,i+1) = 0d0
        vr(i) = 0d0
        svr(i) = 0d0
        vl(i+1) = 0d0
        svl(i+1) = 0d0

        rv2=(vv(i)+sv(i))/(1d0+vv(i)*sv(i))
        rv4=-(vv(i+1)-sv(i+1))/(1d0-vv(i+1)*sv(i+1))
        bplus  = max(0d0,rv2)
        bminus = max(0d0,rv4)

c...Eq.(2.43)
        urh(:,i) = ur(:,i)-bplus*lam/2d0
     &   *(ur(:,i)-ul(:,i)-(1d0-2d0/3d0*bplus*lam)
     &   *6d0*(uu(:,i)-(ur(:,i)+ul(:,i))/2d0))

c...Eq.(2.42)
        ulh(:,i+1)= ul(:,i+1)+bminus*lam/2d0
     & *(ur(:,i+1)-ul(:,i+1)+(1d0-2d0/3d0*bminus*lam)
     & *6d0*(uu(:,i+1)-(ur(:,i+1)+ul(:,i+1))/2d0))

      uh(:) = urh(:,i)/tau
      if (uh(0) .lt. 0d0) uh(0) = 0d0
      det = uh(0)*uh(0)-uh(1)*uh(1)-uh(2)*uh(2)-uh(3)*uh(3)
      if (det .gt. 0d0) then
         call thermal(uh,ene,nba,pre,sva,tplc,mulc,vel)
         fg = uh(0)+tau*pre
         if(fg .ne. 0d0) then
           fr(1,i) = uh(k)*uh(1)/fg + tau*pre*ix
           fr(2,i) = uh(k)*uh(2)/fg + tau*pre*iy
           fr(3,i) = uh(k)*uh(3)/fg + tau*pre*iz
           fr(0,i) = uh(k)
           fr(4,i) = uh(k)*uh(4)/fg
           vr(i) = uh(k)/fg
           svr(i) = sva
         endif
      endif

      uh(:) = ulh(:,i+1)/tau
      if (uh(0) .lt. 0d0) uh(0) = 0d0
      det = uh(0)*uh(0)-uh(1)*uh(1)-uh(2)*uh(2)-uh(3)*uh(3)
      if(det.gt.0d0) then
         call thermal(uh,ene,nba,pre,sva,tplc,mulc,vel)
         fg = uh(0)+tau*pre
         if(fg .ne. 0d0) then
         fl(1,i+1) = uh(k)*uh(1)/fg + tau*pre*ix
         fl(2,i+1) = uh(k)*uh(2)/fg + tau*pre*iy
         fl(3,i+1) = uh(k)*uh(3)/fg + tau*pre*iz
         fl(0,i+1) = uh(k)
         fl(4,i+1) = uh(k)*uh(4)/fg
         vl(i+1) = uh(k)/fg
         svl(i+1) = sva
         endif
      endif

      end do

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c       subroutine flux(vr,vl,svr,svl
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do i = 0,mxx-1

        fh(:,i) = 0d0
        vvh = (vr(i)+vl(i+1))/2d0
        svh = (svr(i)+svl(i+1))/2d0
        rv1 = (vvh+svh)/(1d0+vvh*svh)
        rv2 = (vl(i+1)+svl(i+1))/(1d0+vl(i+1)*svl(i+1))
        rv3 = (vvh-svh)/(1d0-vvh*svh)
        rv4=(vr(i)-svr(i))/(1d0-vr(i)*svr(i))
        bplus = max(0d0,rv1,rv2)
        bminus = min(0d0,rv3,rv4)

c....Eq.(2.49)
        if(bplus-bminus.ne.0d0) then
          fh(:,i) = ( bplus*fr(:,i)-bminus*fl(:,i+1)
     &                +bplus*bminus*(ulh(:,i+1)-urh(:,i)) )
     &               /(bplus-bminus)
        endif

      end do

c...Nextstep
      g(:,0)   = lam*(fh(:,0)-fh(:,1))
      g(:,mxx) = lam*(fh(:,mxx-2)-fh(:,mxx-1))
      do i=1,mxx-1
        g(:,i) = lam*(fh(:,i-1)-fh(:,i))
      end do

c....Eq.(2.50)
c     uu(:,0)   = uu(:,0)  -lam*(fh(:,1)-fh(:,0))
c     uu(:,mxx) = uu(:,mxx)-lam*(fh(:,mxx-1)-fh(:,mxx-2))
c     do i = 1,mxx-1
c       uu(:,i) = uu(:,i)-lam*(fh(:,i)-fh(:,i-1))
c     end do

CCC new terms in tau-eta coordinate 
CCC (devided by 3 due to operater splitting.)
c     if(opt_taueta.eq.1.and.is1 .eq. 2) then
c       uu(3,i) = uu(3,i)-uu(3,i)/tau*dt
c       if(fe(i) .ne. 0d0)then
c         uu(0,i) = uu(0,i)
c    &              -uu(3,i)*uu(3,i)/fe(i)/tau*dt-pressu(i)*dt
c       endif
c     endif

c     call rescaleu(mxx,uu)

      end

c********************************************************************

      subroutine ppm1dp(tau,mxx,dl,dt,uu,vx,sv,g)

c....PPM with periodic boundary condition.
      implicit none
      integer mxx,i,j,ip,im
      real*8 dl,dt,tau
      real*8 uu(0:4,0:mxx), vx(0:mxx),sv(0:mxx),g(0:4,0:mxx)
      real*8 dmu(0:4,0:mxx),ddu(0:4,0:mxx)
      real*8 ur(0:4,0:mxx),ul(0:4,0:mxx),uh(0:4)
      real*8 vr(0:mxx),svr(0:mxx),vl(0:mxx),svl(0:mxx)
      real*8 urh(0:4,0:mxx),ulh(0:4,0:mxx)
      real*8 fr(0:4,0:mxx),fl(0:4,0:mxx)
      real*8 difa,difb,du,difr,difl,mim,det
      real*8 ene,nba,pre,sva,tplc,mulc
      real*8 vvh,svh,rv1,rv2,rv3,rv4,avrl,difrl
      real*8 fh(0:4,0:mxx)
      real*8 bplus,bminus,fg
      real*8 lam,va,adu,gam,vel
      parameter(mim = 0.00001d0)

      lam = dt/dl

c...Calculation of delta_m U#
c...Compute dmu and ddu
c     do i=1,mxx-1
      do i=0,mxx
        ip=i+1
        if(ip.eq.mxx+1) ip=0
        im=i-1
        if(im.eq.-1) im=mxx
      do j=0,4
        difa = uu(j,ip) - uu(j,i)
        difb = uu(j,i)  - uu(j,im)
        du = difa/2d0+ difb/2d0
        if(difa*difb .gt. 0d0) then
          adu=abs(du)
          dmu(j,i) = (min(adu,2d0*abs(difa),2d0*abs(difb)))*du/adu
        else 
          dmu(j,i) = 0d0
        endif
        ddu(j,i) = difa/6d0/dl/dl - difb/6d0/dl/dl
      end do
      end do

c...Average value at right and left state. Compute ur and ul
c     ur(:,0) = 0.5*(uu(:,1)+uu(:,0))
c     ur(:,1) = 0.5*(uu(:,2)+uu(:,1))-dmu(:,2)/6d0+dmu(:,1)/6d0
c     ur(:,mxx-1) = 0.5*(uu(:,mxx)+uu(:,mxx-1))
      ur(:,mxx) = 0.5*(uu(:,0)+uu(:,mxx))
     &           - dmu(:,0)/6d0 + dmu(:,mxx)/6d0
      do i = 0,mxx-1
        ur(:,i) = 0.5*(uu(:,i+1)+uu(:,i))
     &           - dmu(:,i+1)/6d0 + dmu(:,i)/6d0
      end do

c      ul(:,0) = uu(:,0)
c      ul(:,1) = 0.5*(uu(:,1)+uu(:,0))
c      ul(:,mxx-1) = 0.5*(uu(:,mxx-2)+uu(:,mxx-1))
c    &  -dmu(:,mxx-1)/6d0+dmu(:,mxx-2)/6d0
c      ul(:,mxx) = 0.5*(uu(:,mxx-1)+uu(:,mxx))

        ul(:,0) = 0.5*(uu(:,mxx)+uu(:,0))
     &           - dmu(:,0)/6d0   + dmu(:,mxx)/6d0
      do i = 1,mxx
        ul(:,i) = 0.5*(uu(:,i-1)+uu(:,i))
     &           - dmu(:,i)/6d0   + dmu(:,i-1)/6d0
      end do

c...Monotonicity condition
      do i = 0,mxx
      do j=0,4
        difr = ur(j,i) - uu(j,i)
        difl = uu(j,i) - ul(j,i) 
        difrl = difr+difl
        avrl = (ur(j,i)+ul(j,i))/2d0
c....Eq.(2.33)
        if(difr*difl .lt. 0d0) then
         ur(j,i) = uu(j,i)
         ul(j,i) = uu(j,i)
        endif
c...Eq.(2.34)
        if(difrl*(uu(j,i)-avrl).gt. difrl*difrl/6d0) then
          ul(j,i) = 3d0*uu(j,i)-2d0*ur(j,i)
         endif
c...Eq.(2.35)
        if(-difrl*difrl/6d0 .gt.difrl*(uu(j,i)-avrl)) then
         ur(j,i) = 3d0*uu(j,i)-2d0*ul(j,i)
        endif
      end do
      end do

c...Density MUST BE non-negative 
      do i = 0,mxx
        if((ur(0,i) .lt. 0d0).or. (ul(0,i) .lt. 0d0))then
         ur(0,i) = uu(0,i)
         ul(0,i) = uu(0,i)
        endif

c....Is this OK?  What about anti-baryon?
        if((ur(4,i) .lt.0d0) .or. (ul(4,i) .lt. 0d0))then
         ur(4,i) = uu(4,i)
         ul(4,i) = uu(4,i)
        endif

      end do

c...If a cell is in the vacuum,a value at cell faces is taken to be zero.
      do i = 0,mxx
        if(uu(0,i) .eq. 0d0) then
          ur(:,i) = 0d0
          ul(:,i) = 0d0
        endif
       end do

c....Determine v_r, v_l
c     call vrvl(vr,svr,vl,svl,ur,ul,tau,mxx)


c...Compute UL,UR and FL, FR  original subroutine avlr(vv,sr,...)
c     do i = 0,mxx-1
      do i = 0,mxx

        fr(:,i) = 0d0
        vr(i) = 0d0
        svr(i) = 0d0
        rv2=(vx(i)+sv(i))/(1d0+vx(i)*sv(i))
        bplus  = max(0d0,rv2)
c...Eq.(2.43)
        urh(:,i) = ur(:,i)-bplus*lam/2d0
     &   *(ur(:,i)-ul(:,i)-(1d0-2d0/3d0*bplus*lam)
     &   *6d0*(uu(:,i)-(ur(:,i)+ul(:,i))/2d0))

      uh(:) = urh(:,i)/tau
      if (uh(0) .lt. 0d0) uh(0) = 0d0
      det = uh(0)*uh(0)-uh(1)*uh(1)-uh(2)*uh(2)-uh(3)*uh(3)
      if (det .gt. 0d0) then
         call thermal(uh,ene,nba,pre,sva,tplc,mulc,vel)
         fg = uh(0)+tau*pre
         if(fg .ne. 0d0) then
           fr(1,i) = uh(1)*uh(1)/fg + tau*pre
           fr(2,i) = uh(1)*uh(2)/fg
           fr(3,i) = uh(1)*uh(3)/fg
           fr(0,i) = uh(1)
           fr(4,i) = uh(1)*uh(4)/fg
           vr(i)   = uh(1)/fg
           svr(i)  = sva
         endif
      endif

       ip=i+1
       if(ip.eq.mxx+1) ip=0
        fl(:,ip) = 0d0
        vl(ip) = 0d0
        svl(ip) = 0d0
        rv4=-(vx(ip)-sv(ip))/(1d0-vx(ip)*sv(ip))
        bminus = max(0d0,rv4)
c...Eq.(2.42)
        ulh(:,ip)= ul(:,ip)+bminus*lam/2d0
     & *(ur(:,ip)-ul(:,ip)+(1d0-2d0/3d0*bminus*lam)
     & *6d0*(uu(:,ip)-(ur(:,ip)+ul(:,ip))/2d0))

      uh(:) = ulh(:,ip)/tau
      if (uh(0) .lt. 0d0) uh(0) = 0d0
      det = uh(0)*uh(0)-uh(1)*uh(1)-uh(2)*uh(2)-uh(3)*uh(3)
      if(det.gt.0d0) then
         call thermal(uh,ene,nba,pre,sva,tplc,mulc,vel)
         fg = uh(0)+tau*pre
         if(fg .ne. 0d0) then
         fl(1,ip) = uh(1)*uh(1)/fg + tau*pre
         fl(2,ip) = uh(1)*uh(2)/fg
         fl(3,ip) = uh(1)*uh(3)/fg
         fl(0,ip) = uh(1)
         fl(4,ip) = uh(1)*uh(4)/fg
         vl(ip)   = uh(1)/fg
         svl(ip)  = sva
         endif
      endif

      end do

c...Compute flux.
c     do i = 0,mxx-1
      do i = 0,mxx

        ip=i+1
        if(ip.eq.mxx+1) ip=0

        fh(:,i) = 0d0
        vvh = (vr(i)+vl(ip))/2d0
        svh = (svr(i)+svl(ip))/2d0
        rv1 = (vvh+svh)/(1d0+vvh*svh)
        rv2 = (vl(ip)+svl(ip))/(1d0+vl(ip)*svl(ip))
        rv3 = (vvh-svh)/(1d0-vvh*svh)
        rv4=(vr(i)-svr(i))/(1d0-vr(i)*svr(i))
        bplus = max(0d0,rv1,rv2)
        bminus = min(0d0,rv3,rv4)

c....Eq.(2.49)
        if(bplus-bminus.ne.0d0) then
          fh(:,i) = ( bplus*fr(:,i)-bminus*fl(:,ip)
     &                +bplus*bminus*(ulh(:,ip)-urh(:,i)) )
     &               /(bplus-bminus)
        endif

      end do

      g(:,0)   = lam*(fh(:,mxx)-fh(:,0))
      do i=1,mxx
        g(:,i) = lam*(fh(:,i-1)-fh(:,i))
      end do

      end

