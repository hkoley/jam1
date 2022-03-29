c***********************************************************************
c 2019/5/23                                                            *
c        PART  : RQMD/Scalar potential new implementation              *
c                                                                      *
c   List of Subprograms in rough order of relevance with main purpose  *
c      (s = subroutine, f = function, b = block data, e = entry)       *
c  f fenga      to pre-factor from the baryon current                  *
c  s rqmddev1c  to calculate derivatives                               *
c  s jamvpotdot to compute time-derivatives of the baryon current.     *
c                                                                      *
************************************************************************

      subroutine jamrqmd3

c...Purpose: to calculate force in RQMD/S
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,NrqmdCount,icheckMF
      real*8 diff,fengi,fengj,emi,emj
      real*8 dt,emfsq,diffv,gfacs,gfac1,gfac2,gfac
      real*8 fi,fj,em1,em2,pp1,pp2

c...Check which particle feels potentials.
      NrqmdCount=icheckMF(0,pard(1))

c....set p0 to be the free kinetic energy for the calculation of
c....potential.
      if(mstc(112).ge.1) call setenergy(0)
      call jamrqmm3
      call jamepart3(diff,diffv)
      if(mstc(112).eq.0) call setenergy(1)

c...Compute scalar potential self-consistently.
c     if(mstc(109).ge.3) then
c       call jamsvpot(diff,diffv)
c     endif

c...Loop over all particles.
      forcer(:,:) = 0.0d0
      force(:,:) = 0.0d0
      rhog(:)=0.0
      do i=1,nv
        if(MF_on(i)==0) cycle
        if(rho(i).lt.1d-8) cycle
        rhog(i)=rho(i)**(pard(103)-1.0d0)
      end do
 
c...Loop over all particles i.
      do i=1,nv  ! Sum for ith particle (sigma_i) 

         if(MF_on(i).eq.0) cycle
         gfac1=gfacs(i)
         pp1=p(1,i)**2+p(2,i)**2+p(3,i)**2
         emi=sqrt(p(4,i)**2-pp1)
c...Factor in the EOM for scalar part.
c        em1=p(5,i)+pots(i)
c        fengi=em1/sqrt(em1**2+pp1)
         fengi=emi/p(4,i)


c....Factor for the scalar density.
         fi=1.0
         if(mstc(114).ge.1) then
           fi=emi/p(4,i)
           if(mstc(112).eq.1) fi=p(5,i)/(p(5,i)**2+pp1)
         endif

c...Loop over particles j.
         do j=i+1,nv     ! Sum for j (sigma_j(not=i))
           if(MF_on(j)==0) cycle
           pp2=p(1,j)**2+p(2,j)**2+p(3,j)**2
           emj=sqrt(p(4,j)**2-pp2)

c          em2=p(5,j)+pots(j)
c          fengj=em2/sqrt(em2**2+pp2)
           fengj=emj/p(4,j)

           fj=1.0
           if(mstc(114).ge.1) then
             fj=emj/p(4,j)
             if(mstc(112).eq.1) fj=p(5,j)/(p(5,j)**2+pp2)
           endif

           gfac2=gfacs(j)
           gfac=gfac1*gfac2

          if(mstc(113).eq.0) then
            call rqmddev0d(i,j,fengi,fengj,fi,fj,gfac)
          else if(mstc(113).le.1) then
            call rqmddev1d(i,j,fengi,fengj,fi,fj,emi,emj,gfac)
          else
            call rqmddev2d(i,j,fengi,fengj,emi,emj,gfac)
          endif

        end do
      end do

c....set p0 using effective mass.
c     if(mstc(112).ge.1) call setenergy(1)

        dt=parc(2)
        do i=1,nv
          if(MF_on(i)==0) cycle
c           emfsq=p(4,i)**2 - p(1,i)**2 - p(2,i)**2 - p(3,i)**2
            emfsq=(p(5,i)+pots(i))**2

            p(1,i)=p(1,i)+dt*force(1,i)
            p(2,i)=p(2,i)+dt*force(2,i)
            p(3,i)=p(3,i)+dt*force(3,i)

            p(4,i)=sqrt(emfsq+p(1,i)**2+p(2,i)**2+p(3,i)**2)
c displacement only with interaction
            r(1,i)=r(1,i)+dt*forcer(1,i)
            r(2,i)=r(2,i)+dt*forcer(2,i)
            r(3,i)=r(3,i)+dt*forcer(3,i)
        enddo

c       call jamrqmm3
c       call jamepart3(diff,diffv)
c       call setenergy(1)

c....set p0 using effective mass.
c     if(mstc(109).ge.2) call setenergy(1)

      end

**********************************************************************

      subroutine rqmddev0d(i,j,feng1,feng2,fi,fj,gfacts)

c...Purpose: to compute derivatives of the non-relativistic distance for
c...scalar Skyrme.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,n
      real*8 wg,feng1,feng2,gfacts
      real*8 fsky1,fsky2,fsky
      real*8 fsgam2i,fsgam2j,fmgam2i,fmgam2j
      real*8 fi,fj,p2,fac1,fac2,fengij,fmome,fmomd,facmom

c...Width parameter.
c     wg = 1.0/(4.d0*pard(104))
      wg = pard(109)  != 1.0/(4.d0*pard(104))

c...Scalar potential.
      fsky1=feng1*(t1+t3f*rhog(i))*rhom(i,j)*gfacts
      fsky2=feng2*(t1+t3f*rhog(j))*rhom(j,i)*gfacts

c...Total.
      fsky= -wg*(fsky1*fj + fsky2*fi)

      do n=1,3
        force(n,i)  = force(n,i)  - 2*fsky*(r(n,i)-r(n,j)) ! dP_i/dt
        force(n,j)  = force(n,j)  - 2*fsky*(r(n,j)-r(n,i)) ! dP_j/dt
      end do

c...Derivative of m_j/p^0_j in the scalar part.
      if(mstc(111).eq.1) then
        fsgam2i=-mstc(116)*fsky2*fi/p(4,i)
        fsgam2j=-mstc(116)*fsky1*fj/p(4,j)
        do n=1,3
          forcer(n,i) = forcer(n,i) + fsgam2i*p(n,i)/p(4,i)
          forcer(n,j) = forcer(n,j) + fsgam2j*p(n,j)/p(4,j)
        end do
      endif


c...Done if momentum dependent potential is off.
      if (mstd(101).eq.0) return
 
c....Momentum dependent potential part start.
      p2 = (p(1,i)-p(1,j))**2 + (p(2,i)-p(2,j))**2 + (p(3,i)-p(3,j))**2
      fac1 = 1.d0 + p2/pmu1
      fac2 = 1.d0 + p2/pmu2

c...p Derivative term.
      fengij=feng1*fj*rhom(i,j) + feng2*fi*rhom(j,i)
      fmome=-fengij*(vex1/pmu1/fac1**2 + vex2/pmu2/fac2**2)

c...r Derivative term.
      facmom=vex1/fac1 + vex2/fac2
      fmomd=-wg*fengij*facmom

      do n=1,3
        force(n,i)  = force(n,i)  - 2*fmomd*(r(n,i)-r(n,j))
        force(n,j)  = force(n,j)  - 2*fmomd*(r(n,j)-r(n,i))
        forcer(n,i) = forcer(n,i) + 2*fmome*(p(n,i)-p(n,j))
        forcer(n,j) = forcer(n,j) + 2*fmome*(p(n,j)-p(n,i))
      enddo

c...Derivative of  m_i/p^0_i part from the scalar density.
      if(mstc(111).eq.1.and.mstc(114).ne.0) then
        fmgam2i=-mstc(116)*facmom*feng2*fi/p(4,i)*rhom(j,i)
        fmgam2j=-mstc(116)*facmom*feng1*fj/p(4,j)*rhom(i,j)
        do n=1,3
          forcer(n,i) = forcer(n,i) + p(n,i)/p(4,i)*fmgam2i
          forcer(n,j) = forcer(n,j) + p(n,j)/p(4,j)*fmgam2j
        enddo
      endif

      end

**********************************************************************

      subroutine rqmddev1d(i,j,feng1,feng2,gi,gj,emi,emj,gfacts)

c...Purpose: to compute derivatives of the squared four-vector distance
c...q_{Tij}^2 and p_{Tij}^2 in the two body c.m. of particle i and j.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      real*8 feng1,feng2,fengi,fengj
      integer i,j,n
      real*8 bi(0:3),bj(0:3),vi(0:3),vj(0:3)
      real*8 dr2ri(3),dr2pi(3),dp2pi(3)
      real*8 dr2rj(3),dr2pj(3),dp2pj(3)
      real*8 pc1(5),pc2(5),pcm(4),bet(3),rk(3),pk(4)
      real*8 s,rbij,pma,fengij,fsky1,fsky2,fsky
      real*8 gi,gj,gfacts
      real*8 fsgam1,fsgam2i,fsgam2j,fmgam1,facsk
      real*8 fsgam3i,fsgam3j,fmgam2i,fmgam2j
      real*8 p2,fac1,fac2,facmom,fmomd,fmome
      real*8 bbi(3),bbj(3),wg,emi,emj,pe1,pe2

      pc1(:)=p(:,i)
      pc2(:)=p(:,j)
      pcm(4)=pc1(4)+pc2(4)
      bi(0)=1.0
      bj(0)=1.0
      vi(0)=1.0
      vj(0)=1.0
      pe1=p(4,i)
      pe2=p(4,j)
      if(mstc(114).eq.3) then
        pe1=emi
        pe2=emj
        vi(0)=p(4,i)/emi
        vj(0)=p(4,j)/emj
      endif

      do n=1,3
        pcm(n)=pc1(n)+pc2(n)
        bet(n)=pcm(n)/pcm(4)
        rk(n)=r(n,i)-r(n,j)
        pk(n)=pc1(n)-pc2(n)
        bi(n)=pc1(n)/pc1(4)
        bj(n)=pc2(n)/pc2(4)
        vi(n)=-p(n,i)/pe1
        vj(n)=-p(n,j)/pe2
        bbi(n)=bet(n)-mstc(116)*bi(n)
        bbj(n)=bet(n)-mstc(116)*bj(n)
      end do

      s=pcm(4)**2-pcm(1)**2-pcm(2)**2-pcm(3)**2
      rbij=(rk(1)*pcm(1)+rk(2)*pcm(2)+rk(3)*pcm(3))/s
      do n=1,3
        dr2ri(n)= rk(n)+rbij*pcm(n)         ! 1/2*dR~^2_ij/dR_i
        dr2rj(n)= -dr2ri(n)
        dr2pi(n)=(rk(n) + pcm(4)*rbij*bbi(n))*rbij  ! 1/2*dR~^2_ij/dP_i 
        dr2pj(n)=(rk(n) + pcm(4)*rbij*bbj(n))*rbij
      enddo

      fengi=1.0
      fengj=1.0

c...Width parameter.
      wg = 1.0/(4.d0*pard(104))

c...Scalar potential.
      fsky1= feng1*(t1+t3f*rhog(i))*rhom(i,j)*gfacts
      fsky2= feng2*(t1+t3f*rhog(j))*rhom(j,i)*gfacts
      fsky= -wg*(fsky1*gj + fsky2*gi)

c...Coulomb, Symmetry e, and Yukawa potentials are treated as scalar.
c     fsky=fsky + feng1*gj*rhc(i,j)+feng2*gi*rhc(j,i)  ! Add Coulomb potential
c     fsky=fsky + feng1*gj*rhs(i,j)+feng2*gi*rhs(j,i)  ! Add symmetry energy
c     fsky=fsky + feng1*gj*rhy(i,j)+feng2*gi*rhy(j,i)  ! Add Yukawa potential

      fsgam1=0.0  ! p_i(n) + p_j(n) term
      fsgam2i=0.0 ! bi(n) term
      fsgam2j=0.0 ! bi(n) term
      fsgam3i=0.0 ! j_i(n) term
      fsgam3j=0.0 ! j_i(n) term
c...for now option mstc(114)=2,3 is only for Skyrme + mom. dep. force.

c.....Derivatives of gamma_{ij}
      if(mstc(114).eq.1.or.mstc(114).eq.2) then
         facsk = gj*fsky1 + gi*fsky2
         fsgam1= facsk/s
         fsgam2i=mstc(116)*facsk*(1/pcm(4) - pcm(4)/s)
         fsgam2j=fsgam2i
      endif

      if(mstc(111).eq.1) then
      if(mstc(114).eq.2) then

c...Derivative of m_j/p^0_j in the scalar part.
         fsgam2i=fsgam2i-mstc(116)*fsky2*gi/p(4,i)
         fsgam2j=fsgam2j-mstc(116)*fsky1*gj/p(4,j)

c...Derivative of p^\mu_j/m_j
      else if(mstc(114).eq.3) then

         fsgam3i=fsgam3i -rhom(j,i)/emi/p(4,i)
         fsgam3j=fsgam3j -rhom(i,j)/emj/p(4,j)

      endif
      endif

      do n=1,3
        forcer(n,i) = forcer(n,i) + 2*fsky*dr2pi(n) ! dR_i/dt
     &           + pcm(n)*fsgam1 + bi(n)*fsgam2i
        forcer(n,j) = forcer(n,j) + 2*fsky*dr2pj(n) ! dR_j/dt
     &           + pcm(n)*fsgam1 + bj(n)*fsgam2j
        force(n,i) = force(n,i) - 2*fsky*dr2ri(n) ! dP_i/dt
        force(n,j) = force(n,j) - 2*fsky*dr2rj(n) ! dP_j/dt
      end do

c...Done if momentum dependent potential is off.
      if (mstd(101).eq.0) return
 
c....Momentum dependent potential part start.
      pma = (emi**2- emj**2)**2 / s
      pk(4)=pc1(4)-pc2(4)
      p2 = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + pma
      fac1 = 1.d0 + p2/pmu1
      fac2 = 1.d0 + p2/pmu2

      fengij=feng1*gj*rhom(i,j)+feng2*gi*rhom(j,i)
      fmome=-fengij*(vex1/pmu1/fac1**2+vex2/pmu2/fac2**2)

c...Derivative of rho_{ij} term
      facmom=vex1/fac1 + vex2/fac2
      fmomd=-wg*facmom*fengij

      do n=1,3
        dp2pi(n) = pk(n) - mstc(116)*pk(4)*bi(n) + pcm(4)/s*pma*bbi(n)
        dp2pj(n) =-pk(n) + mstc(116)*pk(4)*bj(n) + pcm(4)/s*pma*bbj(n)
        force(n,i)  = force(n,i)  - 2*fmomd*dr2ri(n)
        force(n,j)  = force(n,j)  - 2*fmomd*dr2rj(n)
        forcer(n,i) = forcer(n,i) + 2*fmomd*dr2pi(n) + 2*fmome*dp2pi(n)
        forcer(n,j) = forcer(n,j) + 2*fmomd*dr2pj(n) + 2*fmome*dp2pj(n)
      enddo

      fmgam1=0.0
      fmgam2i=0.0
      fmgam2j=0.0
c...Derivatives of gamma_{ij} part.
      if(mstc(114).ge.1) then
        fmgam1=fengij*facmom/s
        fmgam2i=mstc(116)*fengij*facmom*(1/pcm(4) - pcm(4)/s)
        fmgam2j=fmgam2i
      endif

c...Derivative of  m_j/p^0_j part.
      if(mstc(111).eq.1.and.mstc(114).ge.1) then
        fmgam2i=fmgam2i-mstc(116)*facmom*feng2*gi/pc1(4)*rhom(j,i)
        fmgam2j=fmgam2j-mstc(116)*facmom*feng1*gj/pc2(4)*rhom(i,j)
      endif

      do n=1,3
        forcer(n,i) = forcer(n,i) + pcm(n)*fmgam1
     &                            + bi(n)*fmgam2i
 
        forcer(n,j) = forcer(n,j) + pcm(n)*fmgam1
     &                            + bj(n)*fmgam2j
      enddo

      end

**********************************************************************

      subroutine rqmddev2d(i,j,feng1,feng2,emi,emj,gfacts)

c...Purpose: to compute derivatives of the squared four-vector distance
c...q_{Tij}^2 and p_{Tij}^2 in the rest frame of a particle i or j.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,n
      real*8 feng1,feng2,emi,emj
      real*8 dr2ri(3),dr2rj(3),dr2pi(3),dr2pj(3)
      real*8 dpij2pi(3),dpji2pi(3),dpij2pj(3),dpji2pj(3)
      real*8 fskyi,fskyj,fmomdi,fmomdj,fmomei,fmomej,pp2i,pp2j
      real*8 rk(3),pk(4),bi(3),bj(3),gfacts
      real*8 rbi,rbj,pij,ppk,wg,fac1,fac2
      real*8 vi(0:4),vj(0:4)
      real*8 facm1,facm2,bbi,bbj,bb

c...Width parameter.
      wg = 1.0/(4.d0*pard(104))
      fskyi=-wg*(gfacts*feng1*(t1+t3f*rhog(i)) )*rhom(i,j)
      fskyj=-wg*(gfacts*feng2*(t1+t3f*rhog(j)) )*rhom(j,i)

c     fskyi=fskyi + feng1*rhc(i,j)  ! add Coulomb potential
c     fskyj=fskyj + feng2*rhc(j,i)  ! add Coulomb potential
c     fskyi=fskyi + feng1*rhs(i,j)  ! add symmetry energy
c     fskyj=fskyj + feng2*rhs(j,i)  ! add symmetry energy
c     fskyi=fskyi + feng1*rhy(i,j)  ! add Yukawa potential
c     fskyj=fskyj + feng2*rhy(j,i)  ! add Yukawa potential

      vi(0)=p(4,i)/emi
      vj(0)=p(4,j)/emj
      do n=1,3
        vi(n)=p(n,i)/emi
        vj(n)=p(n,j)/emj
      end do

      rk(1)=r(1,i)-r(1,j)
      rk(2)=r(2,i)-r(2,j)
      rk(3)=r(3,i)-r(3,j)
      rbi=rk(1)*vi(1)+rk(2)*vi(2)+rk(3)*vi(3)
      rbj=rk(1)*vj(1)+rk(2)*vj(2)+rk(3)*vj(3)
      do n=1,3
        dr2ri(n) =  2*(rk(n)+rbj*vj(n))    ! dR~^2_ij/dR_i
        dr2rj(n) =  2*(rk(n)+rbi*vi(n))    ! dR~^2_ji/dR_i
        dr2pi(n) =  2*rk(n)*rbi/emi
        dr2pj(n) =  2*rk(n)*rbj/emj
      end do

      do n=1,3
        force(n,i)  = force(n,i)  - fskyi*dr2ri(n) - fskyj*dr2rj(n)
        force(n,j)  = force(n,j)  + fskyi*dr2ri(n) + fskyj*dr2rj(n)
        forcer(n,i) = forcer(n,i) + fskyj*dr2pi(n)
        forcer(n,j) = forcer(n,j) + fskyi*dr2pj(n)
      end do

 
      ! Momentum dependent potential off.
      if(mstd(101).eq.0) return

      pk(1)=p(1,i)-p(1,j)
      pk(2)=p(2,i)-p(2,j)
      pk(3)=p(3,i)-p(3,j)
      pk(4)=p(4,i)-p(4,j)
      pij=p(4,i)*p(4,j)-p(1,i)*p(1,j)-p(2,i)*p(2,j)-p(3,i)*p(3,j)
      ppk = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2
      pp2i= ppk + (pij-emi**2)**2/emi**2
      pp2j= ppk + (pij-emj**2)**2/emj**2
      facm1 = rhom(i,j)*feng1*gfacts
      facm2 = rhom(j,i)*feng2*gfacts

c     pp2j = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + p2j
      fac1 = 1.d0 + pp2j/pmu1
      fac2 = 1.d0 + pp2j/pmu2
      fmomei=-(vex1/pmu1/fac1**2+vex2/pmu2/fac2**2)*facm1
      fmomdi=-(vex1/fac1+vex2/fac2)*wg*facm1

c     pp2i = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + p2i
      fac1 = 1.d0 + pp2i/pmu1
      fac2 = 1.d0 + pp2i/pmu2
      fmomej=-(vex1/pmu1/fac1**2+vex2/pmu2/fac2**2)*facm2
      fmomdj=-(vex1/fac1+vex2/fac2)*wg*facm2

c     do n=1,3
c       bi(n)=p(n,i)/p(4,i)
c       bj(n)=p(n,j)/p(4,j)
c       bbi=mstc(116)*bi(n)-bj(n)
c       dp2pi(n) = pk(n) - mstc(116)*pk(4)*bi(n) + p(4,j)*p2j*bbi
c       dp2pj(n) =-pk(n) - mstc(116)*pk(4)*bi(n) + p(4,j)*p2i*bbi
c       forcer(n,i) = forcer(n,i) + 2*fmomei*dp2pi(n)
c       forcer(n,j) = forcer(n,j) + 2*fmomej*dp2pj(n)
c       force(n,i) = force(n,i) - 2*fmomdi*dr2ri(n) - 2*fmomdj*dr2rj(n)
c       force(n,j) = force(n,j) + 2*fmomdi*dr2ri(n) + 2*fmomdj*dr2rj(n)
c     enddo

      bbi=2*(pij-emi**2)/emi**2
      bbj=2*(pij-emj**2)/emj**2
      do n=1,3
        bi(n)=p(n,i)/p(4,i)
        bj(n)=p(n,j)/p(4,j)
        bb=mstc(116)*p(4,j)*bi(n)-p(n,j)
        dpij2pi(n) = 2*pk(n) - 2*mstc(116)*pk(4)*bi(n) + bb*bbj ! dpij^2/dp_i
        dpji2pi(n) = 2*pk(n) - 2*mstc(116)*pk(4)*bi(n) + bb*bbi ! dpji^2/pd_i

        bb=mstc(116)*p(4,i)*bj(n)-p(n,i)
        dpij2pj(n) = -2*pk(n) + 2*mstc(116)*pk(4)*bj(n) + bb*bbj ! dpij^2/dp_j
        dpji2pj(n) = -2*pk(n) + 2*mstc(116)*pk(4)*bj(n) + bb*bbi ! dpji^2/pd_j
      end do

      do n=1,3
        forcer(n,i) = forcer(n,i) + fmomei*dpij2pi(n) 
     &          + fmomdj*dr2pi(n) + fmomej*dpji2pi(n)
        forcer(n,j) = forcer(n,j) + fmomej*dpji2pj(n)
c    &          - fmomdi*dr2pj(n) + fmomei*dpij2pj(n)
     &          + fmomdi*dr2pj(n) + fmomei*dpij2pj(n)
        force(n,i) = force(n,i) - fmomdi*dr2ri(n) - fmomdj*dr2rj(n)
        force(n,j) = force(n,j) + fmomdi*dr2ri(n) + fmomdj*dr2rj(n)
      enddo

      end

**********************************************************************

      subroutine jamrqmm3

c...Purpose: to prepare matrix in calculating force
c...Option mstc(113)=2 is not implemented for Coulomb, Yukawa Symmetry force 

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j
      real*8 fac,wg,emisq,emjsq,emi,emj,bi,bj
      real*8 px,py,pz,pe,rx,ry,rz,s,r2,ps,p2,rs
      real*8 qfacr
      real*8 gfac,pmom2ij,pmom2ji
      real*8 gfacs,gfacs1,gfacts,vfac,vfac1,vfacs
      real*8 den,gam,r2i,r2j,p2i,p2j,pij

c     fac  = 1.0/(4.d0*paru(1)*pard(104))**1.5d0 !    [(4*pi*L)^3/2]
c     wg   = 1.0/4.d0/pard(104)

      fac = pard(110)
      wg   = pard(109)

      rho=0.d0   !...  <rho_i> = sigma_j(not= i) rho_ij
      vmoms=0.0
      vmom=0.0

      rhoc=0.d0  ! Coulomb
      rhc=0.d0   ! Coulomb
      rhos=0.d0  ! Symmetry energy
      rhs=0.d0   ! Symmetry energy
      rhoy=0.d0  ! Yukawa potential
      rhy=0.d0   ! Yukawa potential

      do 100 i=1,nv             ! sum for i-th particle (sigma_i) 

        if(MF_on(i).eq.0) goto 100
        bi=k(9,i)/3
        emisq=p(4,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2
        emi=1.0
        if(mstc(114).ge.1) then
          emi=sqrt(emisq)/p(4,i)
          if(mstc(112).eq.1) then
            emi=p(5,i)/sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
          endif
        endif

        gfacs1=gfacs(i)
        vfac1=vfacs(i)

        do 110 j = i+1 , nv      ! Sum for j (sigma_j(>i))

          if(MF_on(j).eq.0) goto 110
          bj=k(9,j)/3
          emjsq = p(4,j)**2 - p(1,j)**2 - p(2,j)**2 - p(3,j)**2
          emj=1.0
          if(mstc(114).ge.1) then
            emj=sqrt(emjsq)/p(4,j)
            if(mstc(112).eq.1) then
              emj=p(5,j)/sqrt(p(5,j)**2+p(1,j)**2+p(2,j)**2+p(3,j)**2)
            endif
          endif

          rx = r(1,i) - r(1,j)
          ry = r(2,i) - r(2,j)
          rz = r(3,i) - r(3,j)
          ps=(p(1,i)-p(1,j))**2+(p(2,i)-p(2,j))**2+(p(3,i)-p(3,j))**2
          rs= rx**2 + ry**2 + rz**2
          gfac=qfacr(i)*qfacr(j)
          gfacts=gfacs1*gfacs(j)
          vfac=vfac1*vfacs(j)

c...Non-relativistic.
          if(mstc(113).eq.0) then

          rhom(i,j) = fac *exp(-rs*wg) * gfac
          rhom(j,i) = rhom(i,j)
          pmom2ij = vex1/(1.d0+ps/pmu1)+vex2/(1.d0+ps/pmu2)
          pmom2ji = pmom2ij
          p2i=ps
          p2j=ps

c...Two body distance is defined by the two particle c.m. frame.
          else if(mstc(113).eq.1) then

          px = p(1,i) + p(1,j)
          py = p(2,i) + p(2,j)
          pz = p(3,i) + p(3,j)
          pe = p(4,i) + p(4,j)
          s = pe**2 - px**2 - py**2 - pz**2
          r2= rs + (rx*px + ry*py + rz*pz)**2/s
          den  = fac * exp(-r2*wg) * gfac
          gam=1.0
          if(mstc(114).ge.1) gam = pe/sqrt(s)
          rhom(i,j) = den*gam
          rhom(j,i) = den*gam

          p2 = ps -(p(4,i)-p(4,j))**2 + (emisq- emjsq)**2/s
          pmom2ij = vex1/(1.d0+p2/pmu1)+vex2/(1.d0+p2/pmu2)
          pmom2ji = pmom2ij
          p2i=p2
          p2j=p2

c...distance is measured from the rest-frame of particle j.
          else

          emi=1.0
          emj=1.0
          r2i=rs + (rx*p(1,i)+ry*p(2,i)+rz*p(3,i))**2/emisq
          r2j=rs + (rx*p(1,j)+ry*p(2,j)+rz*p(3,j))**2/emjsq
          rhom(i,j)= fac * exp(-r2j*wg) * gfac
          rhom(j,i)= fac * exp(-r2i*wg) * gfac
          pij=p(4,i)*p(4,j)-p(1,i)*p(1,j)-p(2,i)*p(2,j)-p(3,i)*p(3,j)
          p2i=ps - (p(4,i)-p(4,j))**2 + (pij-emisq)**2/emisq
          p2j=ps - (p(4,i)-p(4,j))**2 + (pij-emjsq)**2/emjsq
          pmom2ij = vex1/(1.d0+p2j/pmu1)+vex2/(1.d0+p2j/pmu2)
          pmom2ji = vex1/(1.d0+p2i/pmu1)+vex2/(1.d0+p2i/pmu2)

          endif

          if(mstc(118).ge.1) then
            rhom(i,j)=rhom(i,j)/(1.0+p2j/parc(105)**2)
            rhom(j,i)=rhom(j,i)/(1.0+p2i/parc(105)**2)
          endif

          rho(i) = rho(i) + rhom(i,j)*emj*gfacts
          rho(j) = rho(j) + rhom(j,i)*emi*gfacts

          if (mstd(101).eq.0) cycle
          vmoms(i) = vmoms(i) + pmom2ij*rhom(i,j)*emj*gfacts
          vmoms(j) = vmoms(j) + pmom2ji*rhom(j,i)*emi*gfacts

 110     continue               ! Loop end of sigma_j(>i)
 100  continue                  ! Loop end of sigma_i

      end

c***********************************************************************

      subroutine jamepart3(difm,difv)  ! energy

c...Purpose: to calculate single particle potential energy in RQMD/S
c H = sigma_i=1^N  1/(2*E_i) *[E_i^2 -vec{p}_i^2 - m_i^2 - 2m_i*V_i ]
c here V_i = t1 * <rho_i> + t3 *<rho_i>^gamma : Skyrme Potential
c...last revised: 2015/2/20
c...last revised: 2019/3/6

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
c...Local Variable
      integer i
      real*8 difm,difv,emf,pp2,em0

      difm=0.d0   ! Difference of effective mass
      difv=0.d0   ! Difference of vector potential

      do 100 i=1,nv

       potv(:,i)=0.0
       if(MF_on(i)==0) then
         pots(i)=0.0d0
         goto 100
       endif

c...Skyrme potential
       pots(i) = t1*rho(i) + t3*rho(i)**pard(103)
     $       +   rhoc(i) + rhos(i) + rhoy(i) + vmoms(i)

       emf=p(5,i) + pots(i)
       pp2 = p(1,i)**2 + p(2,i)**2 + p(3,i)**2
       em0=sqrt(max(0d0,p(4,i)**2 - pp2))
       difm=difm+abs(em0-emf)
c      p(4,i) = sqrt(emf**2 + pp2)
c      print *,'pot=',pots(i),vmoms(i),t1*rho(i),t3*rho(i)**pard(103)

100   continue

      mste(44)=1
      end

************************************************************************
      subroutine rqmdpotparam3
      implicit none
      real*8 hc,rho0,conv,alpha,beta,gam,C1,C2,mu1,mu2
      include 'jam2.inc'

      hc=paru(3)
      rho0=parc(21)
      conv=hc*hc*hc

      mstd(101)=1  ! Momentum dependent potential is used.

c...New scalar potential implementation.
      if (mstc(106).eq.1) then ! Hard K=380 MeV
        alpha= -0.12378466181689009
        beta=  0.06869780680663275
        gam=  2.1251621959276754
        C1=0.0
        C2=0.0
        mu1=1.0
        mu2=1.0
        mstd(101)=0  ! momentum dependent potential is not used.

      elseif (mstc(106).eq.2) then ! Soft K=210 MeV
        alpha= -0.26209321011086983
        beta=  0.20643197370756433
        gam=  1.2652712441818348
        C1=0.0
        C2=0.0
        mu1=1.0
        mu2=1.0
        mstd(101)=0  ! momentum dependent potential is not used.

      elseif (mstc(106).eq.3) then ! K=380 MeV
        alpha = 0.045166262659641
        beta = 0.03824971687635796
        gam = 2.509835588271163
        C1 = -0.1787443736909779
        mu1 = 3.0754425376340975*hc
        C2 = 0.0
        mu2 = 1.0

        !pard(101)= 0.03488602030826164 
        !pard(102)= 0.040113979691738355
        !pard(103)= 2.312948066341799
        !pard(105)=  2.992730755827161*hc ! pmu1 dummy
        !pard(106)=1.0  ! pmu2 dummy
        !pard(107)=  -0.16701187250281035 ! vex1
        !pard(108)= 0.0d0   ! vex2

      elseif (mstc(106).eq.4) then ! K=210 MeV
        alpha = -0.0896591573281629
        beta = 0.17249044531651028
        gam = 1.202859798739601
        C1 = -0.1787443736909779
        mu1 = 3.0754425376340975*hc
        C2 = 0.0
        mu2 = 1.0

cc--------------------------------------------------------------------------
cc optical potential is defined by sqrt{(m_N+S)^2+p^2} - sqrt{m_N^2 + p^2}
c // m*/m=0.88
      elseif (mstc(106).eq.11) then ! NH2 Skyrme + mom.dep Nara2021 K=380
c   // U( Pinf= 1.7 )= 0.05  Einf= 1.0036086114353737 Pinf= 1.7
c   // U( 0.685 )= 0.0 Elab= 0.22349429615474214
c   // meff=  0.8823864683629444  M*/M= 0.9407105206427979 rho0= 0.168
        alpha = -0.14746618076876553
        beta = 0.28231864726575806
        gam = 1.3086087855271107
        C1 = -0.8500057544566262
        C2= 1.0508877826863514
        mu1 = 2.02*hc
        mu2= 1.0*hc

      elseif (mstc(106).eq.12) then ! MS2 Skyrme + mom.dep Nara2021 K=210
c   // U( Pinf= 1.7 )= 0.05  Einf= 1.0036086114353737 Pinf= 1.7
c   // U( 0.685 )= 0.0 Elab= 0.22349429615474214
c   // meff=  0.8823864683629444  M*/M= 0.9407105206427979 rho0= 0.168
       alpha = -1.7403675255660511
       beta = 1.8746769579585103
       gam = 1.0351617222744887
       C1 = -0.8500057544566262
       C2= 1.0508877826863514
       mu1 = 2.02*hc
       mu2= 1.0*hc

      else
        write(6,*)'RQMDs wrong mode number mstc(106)=',mstc(106)
        stop
      endif

      pard(101)=alpha
      pard(102)=beta
      pard(103)=gam
      pard(105)= mu1
      pard(106)= mu2
      pard(107)= C1
      pard(108)= C2

      end

