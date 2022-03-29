c***********************************************************************
c 2020/2/8                                                             *
c        PART  : RQMD/S scalar implementation Skyrme potential (old)   *
c               Isse et.al Phys. Rev. C (2005)                         *
c                                                                      *
c   List of Subprograms in rough order of relevance with main purpose  *
c      (s = subroutine, f = function, b = block data, e = entry)       *
c                                                                      *
************************************************************************

      subroutine jamrqmd2

c...Purpose: to calculate force in RQMD/S
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,NrqmdCount,icheckMF
      real*8 diff,fengi,fengj,pc1(5),pc2(5)
      real*8 dt,emfsq,diffv,gfacs,gfac1,gfac2,gfac

c...Check which particle feels potentials.
      NrqmdCount=icheckMF(0,pard(1))

      call jamrqmm2
      call jamepart2(diff,diffv)
      call setenergy(1)

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
         pc1(:)=p(:,i)
         if(mstc(112).eq.1) then
           pc1(4)=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
         else
           pc1(5)=sqrt(p(5,i)**2 + 2*p(5,i)*pots(i))
         endif

c...Factor in the EOM for scalar part.
c        fengi=p(5,i)/p(4,i)
         fengi=pc1(5)/pc1(4)

c...Loop over particles j.
         do j=i+1,nv     ! Sum for j (sigma_j(not=i))
           if(MF_on(j)==0) cycle

           pc2(:)=p(:,j)
           if(mstc(112).eq.1) then
             pc2(4)=sqrt(p(5,j)**2+p(1,j)**2+p(2,j)**2+p(3,j)**2)
           else
             pc2(5)=sqrt(p(5,j)**2 + 2*p(5,j)*pots(j))
           endif

c          fengj=p(5,j)/p(4,j)
           fengj=pc2(5)/pc2(4)

           gfac2=gfacs(j)
           gfac=gfac1*gfac2

          if(mstc(113).eq.0) then
            call rqmddev0d2(i,j,fengi,fengj,gfac)
          else if(mstc(113).le.1) then
            call rqmddev1d2(i,j,fengi,fengj,pc1,pc2,gfac)
          else
            call rqmddev2d2(i,j,fengi,fengj,pc1,pc2,gfac)
          endif

        end do
      end do

      dt=parc(2)
      do i=1,nv
        if(MF_on(i)==0) cycle
        emfsq=p(5,i)**2 + 2*p(5,i)*pots(i)

        p(1,i)=p(1,i)+dt*force(1,i)
        p(2,i)=p(2,i)+dt*force(2,i)
        p(3,i)=p(3,i)+dt*force(3,i)

        p(4,i)=sqrt(emfsq+p(1,i)**2+p(2,i)**2+p(3,i)**2)
c displacement only with interaction
        r(1,i)=r(1,i)+dt*forcer(1,i)
        r(2,i)=r(2,i)+dt*forcer(2,i)
        r(3,i)=r(3,i)+dt*forcer(3,i)
      enddo

c       call jamrqmm2
c       call jamepart2(diff,diffv)
c       call setenergy(1)

      end

**********************************************************************

      subroutine rqmddev0d2(i,j,feng1,feng2,gfacts)

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
      real*8 p2,fac1,fac2,fengij,fmome,fmomd,facmom

c...Width parameter.
      wg = pard(109)

c...Scalar potential.
      fsky1=feng1*(t1+t3f*rhog(i))*rhom(i,j)*gfacts
      fsky2=feng2*(t1+t3f*rhog(j))*rhom(j,i)*gfacts

c...Total.
      fsky= -wg*(fsky1 + fsky2)

c...Coulomb, Symmetry e, and Yukawa potentials are treated as scalar.
c     fsky=fsky + feng1*rhc(i,j)+feng2*rhc(j,i)  ! Add Coulomb potential
c     fsky=fsky + feng1*rhs(i,j)+feng2*rhs(j,i)  ! Add symmetry energy
c     fsky=fsky + feng1*rhy(i,j)+feng2*rhy(j,i)  ! Add Yukawa potential

      do n=1,3
        force(n,i)  = force(n,i)  - 2*fsky*(r(n,i)-r(n,j)) ! dP_i/dt
        force(n,j)  = force(n,j)  - 2*fsky*(r(n,j)-r(n,i)) ! dP_j/dt
      end do

c...Done if momentum dependent potential is off.
      if (mstd(101).eq.0) return
 
c....Momentum dependent potential part start.
      p2 = (p(1,i)-p(1,j))**2 + (p(2,i)-p(2,j))**2 + (p(3,i)-p(3,j))**2
      fac1 = 1.d0 + p2/pmu1
      fac2 = 1.d0 + p2/pmu2

c...p Derivative term.
      fengij=feng1*rhom(i,j) + feng2*rhom(j,i)
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

      end

**********************************************************************

      subroutine rqmddev1d2(i,j,feng1,feng2,pc1,pc2,gfacts)

c...Purpose: to compute derivatives of the squared four-vector distance
c...q_{Tij}^2 and p_{Tij}^2 in the two body c.m. of particle i and j.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      real*8 feng1,feng2
      integer i,j,n
      real*8 bi(0:3),bj(0:3)
      real*8 dr2ri(3),dr2pi(3),dp2pi(3)
      real*8 dr2rj(3),dr2pj(3),dp2pj(3)
      real*8 pc1(5),pc2(5),pcm(4),bet(3),rk(3),pk(4)
      real*8 s,rbij,pma,fengij,fsky1,fsky2,fsky
      real*8 gfacts
      real*8 p2,fac1,fac2,facmom,fmomd,fmome
      real*8 bbi(3),bbj(3),wg

      pcm(4)=pc1(4)+pc2(4)
      bi(0)=1.0
      bj(0)=1.0
      do n=1,3
        pcm(n)=pc1(n)+pc2(n)
        bet(n)=pcm(n)/pcm(4)
        rk(n)=r(n,i)-r(n,j)
        pk(n)=pc1(n)-pc2(n)
        bi(n)=pc1(n)/pc1(4)
        bj(n)=pc2(n)/pc2(4)
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

c...Width parameter.
      wg = pard(109)

c...Scalar potential.
      fsky1= feng1*(t1+t3f*rhog(i))*rhom(i,j)*gfacts
      fsky2= feng2*(t1+t3f*rhog(j))*rhom(j,i)*gfacts
      fsky= -wg*(fsky1 + fsky2)

c...Coulomb, Symmetry e, and Yukawa potentials are treated as scalar.
c     fsky=fsky + feng1*rhc(i,j)+feng2*rhc(j,i)  ! Add Coulomb potential
c     fsky=fsky + feng1*rhs(i,j)+feng2*rhs(j,i)  ! Add symmetry energy
c     fsky=fsky + feng1*rhy(i,j)+feng2*rhy(j,i)  ! Add Yukawa potential

      do n=1,3
        forcer(n,i) = forcer(n,i) + 2*fsky*dr2pi(n) ! dR_i/dt
        forcer(n,j) = forcer(n,j) + 2*fsky*dr2pj(n) ! dR_j/dt
        force(n,i) = force(n,i) - 2*fsky*dr2ri(n) ! dP_i/dt
        force(n,j) = force(n,j) - 2*fsky*dr2rj(n) ! dP_j/dt
      end do

c...Done if momentum dependent potential is off.
      if (mstd(101).eq.0) return
 
c....Momentum dependent potential part start.
      pma = (pc1(5)**2- pc2(5)**2)**2 / s
      pk(4)=pc1(4)-pc2(4)
      p2 = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + pma
      fac1 = 1.d0 + p2/pmu1
      fac2 = 1.d0 + p2/pmu2

      fengij=feng1*rhom(i,j)+feng2*rhom(j,i)
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

      end

**********************************************************************

      subroutine rqmddev2d2(i,j,feng1,feng2,pc1,pc2,gfacts)

c...Purpose: to compute derivatives of the squared four-vector distance
c...q_{Tij}^2 and p_{Tij}^2 in the rest frame of a particle i or j.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,n
      real*8 feng1,feng2,pc1(5),pc2(5)
      real*8 dr2ri(3),dr2rj(3),dr2pi(3),dr2pj(3)
      real*8 dpij2pi(3),dpji2pi(3),dpij2pj(3),dpji2pj(3)
      real*8 fskyi,fskyj,fmomdi,fmomdj,fmomei,fmomej,pp2i,pp2j
      real*8 rk(3),pk(5),bi(3),bj(3),gfacts
      real*8 rbi,rbj,pij,p2i,p2j,wg,fac1,fac2
      real*8 vi(5),vj(5)
      real*8 facm1,facm2,bb

c...Width parameter.
      wg = pard(109)
      fskyi=-wg*(gfacts*feng1*(t1+t3f*rhog(i)) )*rhom(i,j)
      fskyj=-wg*(gfacts*feng2*(t1+t3f*rhog(j)) )*rhom(j,i)

c     fskyi=fskyi + feng1*rhc(i,j)  ! add Coulomb potential
c     fskyj=fskyj + feng2*rhc(j,i)  ! add Coulomb potential
c     fskyi=fskyi + feng1*rhs(i,j)  ! add symmetry energy
c     fskyj=fskyj + feng2*rhs(j,i)  ! add symmetry energy
c     fskyi=fskyi + feng1*rhy(i,j)  ! add Yukawa potential
c     fskyj=fskyj + feng2*rhy(j,i)  ! add Yukawa potential

      vi(:)=pc1(:)/pc1(5)
      vj(:)=pc2(:)/pc2(5)

      rk(1)=r(1,i)-r(1,j)
      rk(2)=r(2,i)-r(2,j)
      rk(3)=r(3,i)-r(3,j)
      rbi=rk(1)*vi(1)+rk(2)*vi(2)+rk(3)*vi(3)
      rbj=rk(1)*vj(1)+rk(2)*vj(2)+rk(3)*vj(3)
      do n=1,3
        dr2ri(n) =  2*(rk(n)+rbj*vj(n))    ! dR~^2_ij/dR_i
        dr2rj(n) =  2*(rk(n)+rbi*vi(n))    ! dR~^2_ji/dR_i
        dr2pi(n) =  2*rk(n)*rbi/pc1(5)
        dr2pj(n) =  2*rk(n)*rbj/pc2(5)
      end do

      do n=1,3
        force(n,i)  = force(n,i)  - fskyi*dr2ri(n) - fskyj*dr2rj(n)
        force(n,j)  = force(n,j)  + fskyi*dr2ri(n) + fskyj*dr2rj(n)
        forcer(n,i) = forcer(n,i) + fskyj*dr2pi(n)
        forcer(n,j) = forcer(n,j) + fskyi*dr2pj(n)
      end do

 
      ! Momentum dependent potential off.
      if(mstd(101).eq.0) return

      pk(:)=pc1(:)-pc2(:)
      pij=pc1(4)*pc2(4)-pc1(1)*pc2(1)-pc1(2)*pc2(2)-pc1(3)*pc2(3)
      p2i=(pij-pc1(5)**2)**2/pc1(5)**2
      p2j=(pij-pc2(5)**2)**2/pc2(5)**2
      facm1 = rhom(i,j)*feng1*gfacts
      facm2 = rhom(j,i)*feng2*gfacts

      pp2j = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + p2j
      fac1 = 1.d0 + pp2j/pmu1
      fac2 = 1.d0 + pp2j/pmu2
      fmomei=-(vex1/pmu1/fac1**2+vex2/pmu2/fac2**2)*facm1
      fmomdi=-(vex1/fac1+vex2/fac2)*wg*facm1

      pp2i = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + p2i
      fac1 = 1.d0 + pp2i/pmu1
      fac2 = 1.d0 + pp2i/pmu2
      fmomej=-(vex1/pmu1/fac1**2+vex2/pmu2/fac2**2)*facm2
      fmomdj=-(vex1/fac1+vex2/fac2)*wg*facm2

      do n=1,3
        bi(n)=pc1(n)/pc1(4)
        bj(n)=pc2(n)/pc2(4)
        bb=mstc(116)*pc2(4)*bi(n)-pc2(n)
        dpij2pi(n) = 2*pk(n) - 2*mstc(116)*pk(4)*bi(n) + bb*p2j ! dpij^2/dp_i
        dpji2pi(n) = 2*pk(n) - 2*mstc(116)*pk(4)*bi(n) + bb*p2i ! dpji^2/pd_i

        bb=mstc(116)*pc1(4)*bj(n)-pc1(n)
        dpij2pj(n) = -2*pk(n) + 2*mstc(116)*pk(4)*bj(n) + bb*p2j ! dpij^2/dp_j
        dpji2pj(n) = -2*pk(n) + 2*mstc(116)*pk(4)*bj(n) + bb*p2i ! dpji^2/pd_j
      end do

      do n=1,3
        forcer(n,i) = forcer(n,i) + fmomei*dpij2pi(n) 
     &          + fmomdj*dr2pi(n) + fmomej*dpji2pi(n)
        forcer(n,j) = forcer(n,j) + fmomej*dpji2pj(n)
     &          + fmomdi*dr2pj(n) + fmomei*dpij2pj(n)
        force(n,i) = force(n,i) - fmomdi*dr2ri(n) - fmomdj*dr2rj(n)
        force(n,j) = force(n,j) + fmomdi*dr2ri(n) + fmomdj*dr2rj(n)
      enddo

      end

**********************************************************************

      subroutine jamrqmm2

c...Purpose: to prepare matrix in calculating force
c...Option mstc(113)=2 is not implemented for Coulomb, Yukawa Symmetry force 

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j
      real*8 fac,wg,emisq,emjsq,bi,bj
      real*8 px,py,pz,pe,rx,ry,rz,s,r2,ps,p2,rs
      real*8 qfacr,pe1,pe2
      real*8 gfac,pmom2ij,pmom2ji
      real*8 gfacs,gfacs1,gfacts,vfac,vfac1,vfacs
      real*8 den,r2i,r2j,p2i,p2j,pij

c     fac  = 1.0/(4.d0*paru(1)*pard(104))**1.5d0 !    [(4*pi*L)^3/2]
c     wg   = 1.0/4.d0/pard(104)

      fac = pard(110)
      wg  = pard(109)

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
        if(mstc(112).eq.1) then
          emisq=p(5,i)**2
          pe1=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
        else
          emisq=p(4,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2
          pe1=p(4,i)
        endif
        gfacs1=gfacs(i)
        vfac1=vfacs(i)

        do 110 j = i+1 , nv      ! Sum for j (sigma_j(>i))

          if(MF_on(j).eq.0) goto 110
          bj=k(9,j)/3
          if(mstc(112).eq.1) then
            emjsq = p(5,j)**2
            pe2=sqrt(p(5,j)**2+p(1,j)**2+p(2,j)**2+p(3,j)**2)
          else
            emjsq = p(4,j)**2 - p(1,j)**2 - p(2,j)**2 - p(3,j)**2
            pe2=p(4,j)
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
          pe = pe1    + pe2
          s = pe**2 - px**2 - py**2 - pz**2
          r2= rs + (rx*px + ry*py + rz*pz)**2/s
          den  = fac * exp(-r2*wg) * gfac
          rhom(i,j) = den
          rhom(j,i) = den

          p2 = ps -(pe1-pe2)**2 + (emisq- emjsq)**2/s
          pmom2ij = vex1/(1.d0+p2/pmu1)+vex2/(1.d0+p2/pmu2)
          pmom2ji = pmom2ij
          p2i=p2
          p2j=p2

c...distance is measured from the rest-frame of particle j.
          else

          r2i=rs + (rx*p(1,i)+ry*p(2,i)+rz*p(3,i))**2/emisq
          r2j=rs + (rx*p(1,j)+ry*p(2,j)+rz*p(3,j))**2/emjsq
          rhom(i,j)= fac * exp(-r2j*wg) * gfac
          rhom(j,i)= fac * exp(-r2i*wg) * gfac
          pij=pe1*pe2-p(1,i)*p(1,j)-p(2,i)*p(2,j)-p(3,i)*p(3,j)
          p2i=ps - (pe1-pe2)**2 + (pij-emisq)**2/emisq
          p2j=ps - (pe1-pe2)**2 + (pij-emjsq)**2/emjsq
          pmom2ij = vex1/(1.d0+p2j/pmu1)+vex2/(1.d0+p2j/pmu2)
          pmom2ji = vex1/(1.d0+p2i/pmu1)+vex2/(1.d0+p2i/pmu2)

          endif

          if(mstc(118).ge.1) then
            rhom(i,j)=rhom(i,j)/(1.0+p2j/parc(105)**2)
            rhom(j,i)=rhom(j,i)/(1.0+p2i/parc(105)**2)
          endif

          rho(i) = rho(i) + rhom(i,j)*gfacts
          rho(j) = rho(j) + rhom(j,i)*gfacts

          if (mstd(101).eq.0) cycle
          vmoms(i) = vmoms(i) + pmom2ij*rhom(i,j)*gfacts
          vmoms(j) = vmoms(j) + pmom2ji*rhom(j,i)*gfacts

 110     continue               ! Loop end of sigma_j(>i)
 100  continue                  ! Loop end of sigma_i

      end

c***********************************************************************

      subroutine jamepart2(difm,difv)  ! energy

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

c      emf=p(5,i) + pots(i)
       emf = sqrt(p(5,i)**2 + 2*p(5,i)*pots(i))
       pp2 = p(1,i)**2 + p(2,i)**2 + p(3,i)**2
       em0=sqrt(max(0d0,p(4,i)**2 - pp2))
       difm=difm+abs(em0-emf)
c      p(4,i) = sqrt(emf**2 + pp2)
c      print *,'pot=',pots(i),vmoms(i),t1*rho(i),t3*rho(i)**pard(103)

100   continue

      mste(44)=1
      end

************************************************************************
      subroutine rqmdpotparam2
      implicit none
      real*8 hc,rho0,conv,alpha,beta,gam,C1,C2,mu1,mu2
      include 'jam2.inc'

      hc=paru(3)
      rho0=parc(21)
      conv=hc*hc*hc

      mstd(101)=1  ! Momentum dependent potential is used.

c...Old scalar potential implementation.
      if (mstc(106).eq.1) then ! Hard K=380 MeV
        alpha = -0.124838736062805
        beta  =  0.0707373967100532
        gam   =  2.00261915156202
        C1=0.0
        C2=0.0
        mu1=1.0
        mu2=1.0
        mstd(101)=0  ! momentum dependent potential is not used.

      elseif (mstc(106).eq.2) then ! Soft K=210 MeV
        alpha = -0.310514700691548
        beta  =  0.256413361338797
        gam   =  1.20292928384297
        C1=0.0
        C2=0.0
        mu1=1.0
        mu2=1.0
        mstd(101)=0  ! momentum dependent potential is not used.

      elseif (mstc(106).eq.8) then ! MH2 Skyrme + mom.dep Nara2019 K=380 m*/m=0.848
        alpha = -0.0065641639085685915
        beta = 0.08219105302398698
        gam = 1.7228705539199192
        C1= -0.386579535843693
        C2= 0.3439279165932326
        mu1= 2.02*hc
        mu2= 1.0*hc
 
      elseif (mstc(106).eq.9) then ! MS2 Skyrme + mom.dep Nara2019 K=210
        alpha = -0.31501938078689967
        beta = 0.3874745511416311
        gam = 1.1127892793311491
        C1= -0.38657953584369287
        C2= 0.34392791659323213
        mu1= 2.02*hc
        mu2= 1.0*hc


      elseif (mstc(106).eq.10) then ! MH3 K=380 MeV
        alpha = 0.029643038556651874
        beta = 0.04896806260897878
        gam = 2.121401412903513
        C1= -0.20263202109735332
        mu1= 2.8*hc
        C2= 0.06002279810547833
        mu2= 1.0*hc

      elseif (mstc(106).eq.11) then ! MS3 K=210 MeV
        alpha = -0.6765491747331438
        beta = 0.7551602758987743
        gam = 1.047703735031756
        C1= -0.20263202109735337
        C2= 0.06002279810547853
        mu1= 2.8*hc
        mu2= 1.0 *hc

      elseif (mstc(106).eq.12) then ! MH4 K=380 MeV
       alpha = 0.03862935265120862
       beta = 0.04169348484309908
       gam = 2.280403522511386
       C1= -0.16976586945727515
       C2= 0.0
       mu1= 3.172915754474445*hc
       mu2= 1.0*hc

      elseif (mstc(106).eq.13) then ! MS4 K=210 MeV

        alpha = -0.20784369609866815
        beta = 0.2881665335936613
        gam = 1.1197071552352167
        C1= -0.16976586945850788
        mu1= 3.1729157544574456*hc
        C2= 0.0
        mu2= 1.0*hc

      else
        write(6,*)'RQMD/S wrong mode number mstc(106)=',mstc(106)
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

