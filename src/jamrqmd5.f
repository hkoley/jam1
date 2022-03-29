c***********************************************************************
c                                                                      *
c        PART  : RQMD/Scalar+Vector potential                          *
c                                                                      *
c   List of Subprograms in rough order of relevance with main purpose  *
c      (s = subroutine, f = function, b = block data, e = entry)       *
c  f fenga      to pre-factor from the baryon current                  *
c  s rqmddev1c  to calculate derivatives                               *
c  s jamvpotdot to compute time-derivatives of the baryon current.     *
c                                                                      *
************************************************************************

      subroutine jamrqmd5

c...Purpose: to calculate force in RQMD.RMF mode
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,NrqmdCount,icheckMF,bar
      real*8 diff,fengi,fengj
      real*8 dt,emfsq,diffv,gfacs,gfac1,gfac2,gfac,vfac1,vfac,vfacs
      real*8 A1(0:3),A2(0:3),Ai(0:3),Aj(0:3),pkin1(5),pkin2(5)
      real*8 ms2,gs,g2,g3,g4,sigma,sig1,sig2,g24,g34,gs4,gg4,rscalar
      real*8 vdot
      real*8 fi,fj,pc1(5),pc2(5)
      common/jamvdot/vdot(3,mxv)
      real*8 potv0(0:3,mxv)
      real*8 x,dsigma
c     dsigma(x)=gs/(ms2+g2*x+g3*x**2)
      dsigma(x)=gs*(1.0+gg4*x**2)/(ms2 + g2*x + g3*x**2
     &       + g24*x**3 + g34*x**4 - 2*(1.0+gg4*x**2)*gs4*x*rscalar)

      ms2=(paru(3)*pard(112))**2  ! sigma mass square (GeV^2)
      gs=pard(101)*paru(3)**3     ! Scalar part g_sigma
      g2=2*pard(102)*paru(3)      ! Scalar part g2 (GeV)
      g3=3*pard(103)              ! Scalar part g3
      g4=pard(116)
      g24=2.0/3.0*g2*g4/paru(3)
      g34=5.0/4.0*g3*g4/paru(3)**2
      gs4=g4*gs/paru(3)**2
      gg4=0.5*g4/paru(3)**2

c...Check which particle feels potentials.
      NrqmdCount=icheckMF(0,pard(1))

      potv0(:,:)=potv(:,:)
      vdot(:,:)=0.0d0
      call jamrqmm5
      call jamepart5(diff,diffv)
      call setenergy(1)

c...Compute scalar and vector potential self-consistently.
      if(mstc(109).ge.3) then
        call jamsvpot5(diff,diffv)
      endif

c...Loop over all particles.
      rhog(:)=0.0
      forcer(:,:) = 0.0d0
      force(:,:) = 0.0d0
      do i=1,nv
        if(MF_on(i)==0) cycle
        if(rhoi(i).lt.1d-8) cycle
        rhog(i)=rhoi(i)**(pard(113)-1.0d0)
      end do
 
c...Loop over all particles i.
      do i=1,nv  ! Sum for ith particle (sigma_i) 

         if(MF_on(i).eq.0) cycle
         pkin1(:)=p(:,i)
         pkin1(5)=p(5,i)+pots(i)
         pkin1(4)=sqrt(pkin1(5)**2+pkin1(1)**2+pkin1(2)**2+pkin1(3)**2)
         if(mstc(119).eq.0) call getpkin(i,pkin1)

c...Factor in the EOM for scalar and vector part.
         fengi=pkin1(5)/pkin1(4)
c        call fenga5(i,pkin1,Ai)

c....Factor for the scalar density.
         pc1=p(:,i)
         if(mstc(112).eq.1) call getmom5(i,pc1)
         pc1(5)=sqrt(pc1(4)**2-pc1(1)**2-pc1(2)**2-pc1(3)**2)
         fi=pc1(5)/pc1(4)

         gfac1=gfacs(i)
         vfac1=vfacs(i)
c...For non-linear sigma-field.
         sig1=1.0
         if(mstd(102).eq.1) then
           sigma=(pots(i)-vmoms(i))/(t1*gfac1)
           rscalar=rho(i)
           sig1=dsigma(sigma)
c          sig1=gs/(ms2+g2*sigma+g3*sigma**2)
c          sig1=gs*(1.0+gg4*sigma**2)/(ms2+g2*sigma+g3*sigma**2
c    &         +g24*sigma**3 + g34*sigma**4
c    &         -2*(1.0+gg4*sigma**2)*gs4*sigma*rho(i))
         endif

c...Loop over particles j.
         do j=i+1,nv     ! Sum for j (sigma_j(not=i))
           if(MF_on(j)==0) cycle
           pkin2(:)=p(:,j)
           pkin2(5)=p(5,j)+pots(j)
          pkin2(4)=sqrt(pkin2(5)**2+pkin2(1)**2+pkin2(2)**2+pkin2(3)**2)
           if(mstc(119).eq.0) call getpkin(j,pkin2)
           fengj=pkin2(5)/pkin2(4)
           call fenga5(i,pkin2,Ai)
           call fenga5(j,pkin1,Aj)

           pc2=p(:,j)
           if(mstc(112).eq.1) call getmom5(j,pc2)
           pc2(5)=sqrt(pc2(4)**2-pc2(1)**2-pc2(2)**2-pc2(3)**2)
           fj=pc2(5)/pc2(4)

           gfac2=gfacs(j)
           gfac=gfac1*gfac2
           vfac=vfac1*vfacs(j)
           sig2=1.0
           if(mstd(102).eq.1) then
             sigma=(pots(j)-vmoms(j))/(t1*gfac2)
c            sig2=gs/(ms2+g2*sigma+g3*sigma**2)
             rscalar=rho(j)
             sig2=dsigma(sigma)
           endif

           bar=k(9,i)*k(9,j)/9
           A1(:)=bar*Ai(:)
           A2(:)=bar*Aj(:)

          if(mstc(113).eq.0) then
            call rqmddev0c5(i,j,fengi,fengj,fi,fj,A1,A2,pkin1,pkin2,
     &      pc1,pc2,gfac,vfac,sig1,sig2,bar)
          else if(mstc(113).le.1) then
            call rqmddev1c5(i,j,fengi,fengj,fi,fj,A1,A2,
     &      pkin1,pkin2,pc1,pc2,gfac,vfac,sig1,sig2,bar)
          else
            call rqmddev2c5(i,j,fengi,fengj,A1,A2,
     &       pkin1,pkin2,pc1,pc2,gfac,vfac,sig1,sig2,bar)
          endif

        end do
      end do

c...Compute the time derivative of the vector potential for the kinetic
c...momentum update in case vector potential is fully included.
      if(mstc(115).ge.2) then
        if(mstc(119).eq.0.and.mstc(124).eq.1) then
          call jamvpotdot5
        else if(mstc(119).eq.1) then
          call jamvpotdot5
        endif
      endif

      dt=parc(2)
      do i=1,nv
        if(MF_on(i)==0) cycle

        emfsq=(p(5,i)+pots(i))**2
        if(mstc(119).eq.1.or.mstc(119).eq.2) then
          p(1,i)=p(1,i)+dt*(force(1,i) - vdot(1,i))
          p(2,i)=p(2,i)+dt*(force(2,i) - vdot(2,i))
          p(3,i)=p(3,i)+dt*(force(3,i) - vdot(3,i))
        else
          p(1,i)=p(1,i)+dt*force(1,i)
          p(2,i)=p(2,i)+dt*force(2,i)
          p(3,i)=p(3,i)+dt*force(3,i)
        endif

        p(4,i)=sqrt(emfsq+p(1,i)**2+p(2,i)**2+p(3,i)**2)
c displacement only with interaction
        r(1,i)=r(1,i)+dt*forcer(1,i)
        r(2,i)=r(2,i)+dt*forcer(2,i)
        r(3,i)=r(3,i)+dt*forcer(3,i)
      enddo

c     potv(:,:)=potv0(:,:)
      call jamrqmm5
      call jamepart5(diff,diffv)

      if(mstc(119).eq.3) then
        do i=1,nv
          if(MF_on(i)==0) cycle
          p(1,i)=p(1,i)-dt*vdot(1,i)
          p(2,i)=p(2,i)-dt*vdot(2,i)
          p(3,i)=p(3,i)-dt*vdot(3,i)
        end do
c       call jamrqmm5
c       call jamepart5(diff,diffv)
      endif

c...vector potential is updated by using  vdot.
      if(mstc(124).eq.1.and.mstc(115).ge.2) then
        do i=1,nv
          if(MF_on(i)==0) cycle
          potv(1,i)=potv0(1,i)+vdot(1,i)*dt
          potv(2,i)=potv0(2,i)+vdot(2,i)*dt
          potv(3,i)=potv0(3,i)+vdot(3,i)*dt
        end do
      endif

c...Recover total energy conservation.
      if(mstc(123).eq.1) then
        call recovere5
      elseif(mstc(123).eq.2) then
        call recovere5all
      endif

      end

**********************************************************************
      subroutine getpkin(i,pk)

c...Purpose: to compute kinetic momentum from canonical momentum.
      implicit none
      include 'jam1.inc'
      include 'jam3.inc'
      integer i
      real*8 pk(5)

      pk(1)=p(1,i)-potv(1,i)
      pk(2)=p(2,i)-potv(2,i)
      pk(3)=p(3,i)-potv(3,i)
      pk(5)=p(5,i)+pots(i)
      pk(4)=sqrt(pk(5)**2+pk(1)**2+pk(2)**2+pk(3)**2)

      end

**********************************************************************
      subroutine getmom5(i,pk)

c...Purpose: to compute canonical momentum from kinetic momentum.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      integer i
      real*8 pk(5)

      if(mstc(119).ge.1) then
        pk(1)=p(1,i)+potv(1,i)
        pk(2)=p(2,i)+potv(2,i)
        pk(3)=p(3,i)+potv(3,i)
      else
        pk(1)=p(1,i)
        pk(2)=p(2,i)
        pk(3)=p(3,i)
      endif
      pk(5)=p(5,i)
      pk(4)=sqrt(pk(5)**2+pk(1)**2+pk(2)**2+pk(3)**2)

      end
************************************************************************
      subroutine fenga5(i,pkin,Ai)

c....Pre-factor from the baryon current.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i
      real*8 bi(0:3),Ai(0:3),vj,vv,dv
      real*8 pkin(5)

      Ai=0.0
c     if(abs(rhoi(i)).lt.1d-7) return
      if(rhoi(i).lt.1d-8) return

c...Vector potential V_i^\mu is fully included.
      if(mstc(115).eq.2) then
        bi(0)=1.0
        bi(1)=pkin(1)/pkin(4)
        bi(2)=pkin(2)/pkin(4)
        bi(3)=pkin(3)/pkin(4)

c...Only time-component V_i^0 term is included.
      else if(mstc(115).eq.1) then
        bi(0)=1.0
        bi(1)=0.0
        bi(2)=0.0
        bi(3)=0.0

c...Only time-component V^0 term is included with the form of V^i(rho_B).
      else
        dv = (t2 + t3f*rhog(i) + t4f*rhoi(i)**(pard(115)-1))/rhoi(i)
        Ai(:)=dv*rhoj(:,i)
        return
      endif

      vj=rhoj(0,i)-bi(1)*rhoj(1,i)-bi(2)*rhoj(2,i)-bi(3)*rhoj(3,i)
      vv = t2 + t3*rhog(i) + t4*rhoi(i)**(pard(115)-1)  ! V/rho_B
      dv = (pard(113)-1)*t3*rhoi(i)**(pard(113)-3)
     &   + (pard(115)-1)*t4*rhoi(i)**(pard(115)-3)
      Ai(:)=vj*dv*rhoj(:,i) + vv*bi(:)

      end

**********************************************************************

      subroutine rqmddev0c5(i,j,feng1,feng2,fi,fj,Ai,Aj,pkin1,pkin2,
     &  pc1,pc2,gfacts,vfac,sig1,sig2,bar)

c...Purpose: to compute derivatives of the non-relativistic distance
c...Here scalar and vector potentials are implemented as
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,n,bar
      real*8 pc1(5),pc2(5),bi(0:3),bj(0:3),Ai(0:3),Aj(0:3)
      real*8 wg,feng1,feng2,fengi,fengj,gfacts,vfac
      real*8 fsky1,fsky2,fsky3,fsky4,fsky,sig1,sig2
      real*8 fsgam2i,fsgam2j,fsgam3i,fsgam3j,fmgam2i,fmgam2j
      real*8 fi,fj,pkin1(5),pkin2(5),vf1,vf2,fm1,fm2,fm3,fm4
      real*8 p2,fac1,fac2,fengij,fmome,fmomd,facmom,facmomv

      bi(0)=1.0
      bj(0)=1.0
      do n=1,3
        bi(n)=pc1(n)/pc1(4)
        bj(n)=pc2(n)/pc2(4)
      end do

c...Factors for vector potential
c...V_i^\mu = V_i/rho_B J_i^\mu
      fengi=Ai(0)*bi(0)
      fengj=Aj(0)*bj(0)
      do n=1,3
        fengi = fengi - Ai(n)*bi(n)
        fengj = fengj - Aj(n)*bj(n)
      end do

c...Width parameter.
      wg = pard(109)  ! wg = 1.0/(4.d0*pard(104))

c...Scalar potential.
      fsky1=feng1*fj*t1*sig1*rhom(i,j)*gfacts
      fsky2=feng2*fi*t1*sig2*rhom(j,i)*gfacts

c...Vector part.
      fsky3= fengi*rhom(i,j)*vfac
      fsky4= fengj*rhom(j,i)*vfac

c...Total.
      fsky= -wg*(fsky1 + fsky2 + fsky3 + fsky4)

      do n=1,3
        force(n,i)  = force(n,i)  - 2*fsky*(r(n,i)-r(n,j)) ! dP_i/dt
        force(n,j)  = force(n,j)  - 2*fsky*(r(n,j)-r(n,i)) ! dP_j/dt
      end do

      if(mstc(111).eq.1) then
c...Derivative of m_j/p^0_j in the scalar part.
      do n=1,3
        forcer(n,i) = forcer(n,i) - bi(n)*fsky2/pc1(4)*mstc(116)
        forcer(n,j) = forcer(n,j) - bj(n)*fsky1/pc2(4)*mstc(116)
      end do

c....Derivatives of p^\mu/_i/p^0_i term in the vector part.
        call fenga5(i,pkin1,Ai)
        call fenga5(j,pkin2,Aj)
        fsgam2i= mstc(116)*rhom(j,i)*
     &         ( Aj(1)*bi(1) + Aj(2)*bi(2) + Aj(3)*bi(3) )/pc1(4)
        fsgam3i=-rhom(j,i)/pc1(4)

        fsgam2j= mstc(116)*rhom(i,j)*
     &         ( Ai(1)*bj(1) + Ai(2)*bj(2) + Ai(3)*bj(3) )/pc2(4)
        fsgam3j=-rhom(i,j)/pc2(4)

      do n=1,3
        forcer(n,i) = forcer(n,i) + bi(n)*fsgam2i + Aj(n)*fsgam3i
        forcer(n,j) = forcer(n,j) + bj(n)*fsgam2j + Ai(n)*fsgam3j
      end do

      endif

c...Done if momentum dependent potential is off.
      if (mstd(101).eq.0) return

      vf1=1.0
      vf2=1.0
c...V^mu(p) is fully included.
      if(mstc(115).eq.2) then
        do n=1,3
          vf1 = vf1 - pkin1(n)/pkin1(4)*bj(n)
          vf2 = vf2 - pkin2(n)/pkin2(4)*bi(n)
        end do
      endif
      vf1 = vf1*bar
      vf2 = vf2*bar

c....Momentum dependent potential part start.
      p2 = (pc1(1)-pc2(1))**2 + (pc1(2)-pc2(2))**2 + (pc1(3)-pc2(3))**2
      fac1 = 1.d0 + p2/pmu1
      fac2 = 1.d0 + p2/pmu2

c...p Derivative term (scalar).
      fengij=feng1*fj*rhom(i,j) + feng2*fi*rhom(j,i)
      fmome=-fengij*vex1/pmu1/fac1**2

c...r Derivative term (scalar).
      facmom=vex1/fac1
      fmomd=-wg*fengij*facmom

c...p Derivative term (vector).
      fengij = vf1*rhom(i,j) + vf2*rhom(j,i)
      fmome=fmome - fengij*vex2/pmu2/fac2**2

c...r Derivative term (vector).
      facmomv=vex2/fac2
      fmomd=fmomd - wg*fengij*facmomv

      do n=1,3
        force(n,i)  = force(n,i)  - 2*fmomd*(r(n,i)-r(n,j))
        force(n,j)  = force(n,j)  - 2*fmomd*(r(n,j)-r(n,i))
        forcer(n,i) = forcer(n,i) + 2*fmome*(pc1(n)-pc2(n))
        forcer(n,j) = forcer(n,j) + 2*fmome*(pc2(n)-pc1(n))
      enddo

c...Derivative of  m_i/p^0_i part from the scalar density.
      if(mstc(111).eq.1) then
        fmgam2i=-mstc(116)*facmom*feng2*fi/pc1(4)*rhom(j,i)
        fmgam2j=-mstc(116)*facmom*feng1*fj/pc2(4)*rhom(i,j)
        do n=1,3
          forcer(n,i) = forcer(n,i) + pc1(n)/pc1(4)*fmgam2i
          forcer(n,j) = forcer(n,j) + pc2(n)/pc2(4)*fmgam2j
        enddo
      endif

c...Derivative of  p_i/p^0_i part.
      if(mstc(111).eq.1.and.mstc(115).eq.2) then
        fm1=-rhom(j,i)*facmomv/pc1(4)/pkin2(4)
        fm2=-rhom(i,j)*facmomv/pc2(4)/pkin1(4)
c       vf2 = (pkin2(1)*bi(1)+pkin2(2)*bi(2)+pkin2(3)*bi(3))/pkin2(4)
c       vf1 = (pkin1(1)*bj(1)+pkin1(2)*bj(2)+pkin1(3)*bj(3))/pkin1(4)
        vf1 = bar - vf1
        vf2 = bar - vf2
        fm3=mstc(116)*facmomv/pc1(4)*rhom(j,i)*vf1
        fm4=mstc(116)*facmomv/pc2(4)*rhom(i,j)*vf2
        do n=1,3
          forcer(n,i) = forcer(n,i) + fm1*pkin2(n) + bi(n)*fm3
          forcer(n,j) = forcer(n,j) + fm2*pkin1(n) + bj(n)*fm4
        enddo
      endif


      end

**********************************************************************

      subroutine rqmddev1c5(i,j,feng1,feng2,gi,gj,Ai,Aj,
     &      pkin1,pkin2,pc1,pc2,gfacts,vfac,sig1,sig2,bar)

c...Purpose: to compute derivatives of the squared four-vector distance
c...q_{Tij}^2 and p_{Tij}^2 in the two body c.m. of particle i and j.
c...Here attractive part of the Skyrme potential is treated as Lorentz scalar,
c...repulsive part is treated as the Lorentz vector.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      real*8 feng1,feng2,fengi,fengj,sig1,sig2
      integer i,j,n,bar
      real*8 bi(0:3),bj(0:3)
      real*8 dr2ri(3),dr2pi(3),dp2pi(3)
      real*8 dr2rj(3),dr2pj(3),dp2pj(3)
      real*8 pc1(5),pc2(5),pcm(4),bet(3),rk(3),pk(4)
      real*8 s,rbij,pma,fengij,fsky1,fsky2,fsky
      real*8 fsky3,fsky4,gi,gj,gfacts,vfac
      real*8 fsgam1,fsgam2i,fsgam2j,fmgam1,facsk
      real*8 fsgam3i,fsgam3j,fmgam2i,fmgam2j
      real*8 p2,fac1,fac2,facmom,fmomd,fmome,facmomv,fengv,facij
      real*8 bbi(3),bbj(3),wg
      real*8 Ai(0:3),Aj(0:3),pkin1(5),pkin2(5)
      real*8 vf1,vf2,fm1,fm2,fm3,fm4,facm1,facm2

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

c...Factors for vector potential
c...V_i^\mu = V_i/rho_B J_i^\mu
      fengi=Ai(0)*bj(0)
      fengj=Aj(0)*bi(0)
      do n=1,3
        fengi = fengi - Ai(n)*bi(n)
        fengj = fengj - Aj(n)*bj(n)
      end do

c...Width parameter.
c     wg = 1.0/(4.d0*pard(104))
      wg = pard(109)

c...Scalar potential.
      fsky1=feng1*t1*sig1*rhom(i,j)*gfacts
      fsky2=feng2*t1*sig2*rhom(j,i)*gfacts

c...Vector part.
      fsky3= fengi*rhom(i,j)*vfac
      fsky4= fengj*rhom(j,i)*vfac
      fsky= -wg*(fsky1*gj + fsky2*gi + fsky3 + fsky4)

      fsgam1=0.0  ! p_i(n) + p_j(n) term
      fsgam2i=0.0 ! bi(n) term
      fsgam2j=0.0 ! bi(n) term
      fsgam3i=0.0 ! j_i(n) term
      fsgam3j=0.0 ! j_i(n) term
c...for now option mstc(114)=2,3 is only for Skyrme + mom. dep. force.

c.....Derivatives of gamma_{ij}
      if(mstc(114).ne.3) then
         facsk = gj*fsky1 + gi*fsky2 + fsky3 + fsky4
         fsgam1= facsk/s
         fsgam2i=mstc(116)*facsk*(1/pcm(4) - pcm(4)/s)
         fsgam2j=fsgam2i
      endif

      if(mstc(111).eq.1) then
      if(mstc(114).eq.2) then

c...Derivative of m_j/p^0_j in the scalar part.
         fsgam2i=fsgam2i-mstc(116)*fsky2*gi/pc1(4)
         fsgam2j=fsgam2j-mstc(116)*fsky1*gj/pc2(4)

c....Derivatives of p^\mu/_i/p^0_i term in the vector part.
        call fenga5(i,pkin1,Ai)
        call fenga5(j,pkin2,Aj)
         fsgam2i=fsgam2i + mstc(116)*rhom(j,i)*
     &         ( Aj(1)*bi(1) + Aj(2)*bi(2) + Aj(3)*bi(3) )/pc1(4)
         fsgam3i=-rhom(j,i)/pc1(4)

         fsgam2j=fsgam2j + mstc(116)*rhom(i,j)*
     &         ( Ai(1)*bj(1) + Ai(2)*bj(2) + Ai(3)*bj(3) )/pc2(4)
         fsgam3j=-rhom(i,j)/pc2(4)

c...Derivative of p^\mu_j/m_j
      else if(mstc(114).eq.3) then

         fsgam2i=fsgam2i + mstc(116)*rhom(j,i)*Aj(0)/pc1(5)
         fsgam3i=fsgam3i -rhom(j,i)/pc1(5)/pc1(4)
         fsgam2j=fsgam2j + mstc(116)*rhom(j,i)*Ai(0)/pc2(5)
         fsgam3j=fsgam3j -rhom(i,j)/pc2(5)/pc2(4)

      endif
      endif

      do n=1,3
        forcer(n,i) = forcer(n,i) + 2*fsky*dr2pi(n) ! dR_i/dt
     &           + pcm(n)*fsgam1 + bi(n)*fsgam2i + Aj(n)*fsgam3i
        forcer(n,j) = forcer(n,j) + 2*fsky*dr2pj(n) ! dR_j/dt
     &           + pcm(n)*fsgam1 + bj(n)*fsgam2j + Ai(n)*fsgam3j
        force(n,i) = force(n,i) - 2*fsky*dr2ri(n) ! dP_i/dt
        force(n,j) = force(n,j) - 2*fsky*dr2rj(n) ! dP_j/dt
      end do

c...Done if momentum dependent potential is off.
      if (mstd(101).eq.0) return

c...Factors when vector potential is fully included.
      vf1=1.0
      vf2=1.0
      if(mstc(115).eq.2) then
        do n=1,3
          vf1 = vf1 - pkin1(n)/pkin1(4)*bj(n)
          vf2 = vf2 - pkin2(n)/pkin2(4)*bi(n)
        end do
      endif
      vf1=vf1*bar
      vf2=vf2*bar
 
c....Momentum dependent potential part start.
      pma =  (pc1(5)**2- pc2(5)**2)**2 / s
      pk(4)= (pc1(4)-pc2(4))
      p2 = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + pma
      fac1 = 1.d0 + p2/pmu1
      fac2 = 1.d0 + p2/pmu2

c...Derivative of p_{Tij} term (scalar).
      fengij=feng1*gj*rhom(i,j)+feng2*gi*rhom(j,i)
      fmome=-fengij*vex1/pmu1/fac1**2

c...Derivative of rho_{ij} term (scalar)
      facmom=vex1/fac1
      fmomd=-wg*facmom*fengij

c...Derivative of p_{Tij} term (vector).
      fengv=vf1*rhom(i,j) + vf2*rhom(j,i)
      fmome=fmome-fengv*vex2/pmu2/fac2**2

c...Derivative of rho_{ij} term (vector)
      facmomv=vex2/fac2
      fmomd=fmomd-wg*facmomv*fengv

      do n=1,3
        dp2pi(n) = pk(n) - mstc(116)*pk(4)*bi(n) + pcm(4)/s*pma*bbi(n)
        dp2pj(n) =-pk(n) + mstc(116)*pk(4)*bj(n) + pcm(4)/s*pma*bbj(n)
        force(n,i)  = force(n,i)  - 2*fmomd*dr2ri(n)
        force(n,j)  = force(n,j)  - 2*fmomd*dr2rj(n)
        forcer(n,i) = forcer(n,i) + 2*fmomd*dr2pi(n) + 2*fmome*dp2pi(n)
        forcer(n,j) = forcer(n,j) + 2*fmomd*dr2pj(n) + 2*fmome*dp2pj(n)
      enddo

c     fmgam1=0.0
c     fmgam2i=0.0
c     fmgam2j=0.0
      if(mstc(111).eq.1) then

c...Derivatives of gamma_{ij} part.
        facij=fengij*facmom+fengv*facmomv
        fmgam1=facij/s
        fmgam2i=mstc(116)*facij*(1/pcm(4) - pcm(4)/s)
        fmgam2j=fmgam2i

c...Derivative of  m_j/p^0_j part.
        fmgam2i=fmgam2i-mstc(116)*facmom*feng2*gi/pc1(4)*rhom(j,i)
        fmgam2j=fmgam2j-mstc(116)*facmom*feng1*gj/pc2(4)*rhom(i,j)

      do n=1,3
        forcer(n,i) = forcer(n,i) + pcm(n)*fmgam1
     &                            + bi(n)*fmgam2i
        forcer(n,j) = forcer(n,j) + pcm(n)*fmgam1
     &                            + bj(n)*fmgam2j
      enddo

c...Derivative of  p_i/p^0_i part.
      if(mstc(115).eq.2) then
        facm1=facmomv*rhom(j,i)
        facm2=facmomv*rhom(i,j)
        fm1=-facm1/pc1(4)/pkin2(4)
        fm2=-facm2/pc2(4)/pkin1(4)
        vf1 = bar - vf1
        vf2 = bar - vf2
        fm3=mstc(116)*facm1/pc1(4)*vf2
        fm4=mstc(116)*facm2/pc2(4)*vf1
        do n=1,3
          forcer(n,i) = forcer(n,i) + fm1*pkin2(n) + bi(n)*fm3
          forcer(n,j) = forcer(n,j) + fm2*pkin1(n) + bj(n)*fm4
        enddo

      endif
      endif

      end

**********************************************************************

      subroutine rqmddev2c5(i,j,feng1,feng2,Ai,Aj,pkin1,pkin2,
     & pc1,pc2,gfacts,vfac,sig1,sig2,bar)

c...Purpose: to compute derivatives of the squared four-vector distance
c...q_{Tij}^2 and p_{Tij}^2 in the rest frame of a particle i or j.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,n,bar
      real*8 feng1,feng2,fengi,fengj,emi,emj,sig1,sig2,s,pcm(5),pma
      real*8 pc1(5),pc2(5),dr2ri(3),dr2rj(3),dr2pi(3),dr2pj(3)
      real*8 fskyi,fskyj,fmomdi,fmomdj,fmomei,fmomej,pp2i,pp2j
      real*8 bet(3),rk(3),pk(4),bi(3),bj(3),gfacts,vfac
      real*8 rbi,rbj,pij,bbi(3),bbj(3),wg,fac1,fac2,facmom
      real*8 vi(0:4),vj(0:4),dvi(3),dvj(3),bb1,bb2,ppk
      real*8 Ai(0:3),Aj(0:3),pkin1(5),pkin2(5)
      real*8 vf1,vf2,facm1,facm2,fm1,fm2,fm3,fm4,bb
      real*8 dpij2pi(3),dpji2pi(3),dpij2pj(3),dpji2pj(3)

c...Factors for vector potential
      emi=pc1(5)
      emj=pc2(5)
      vi(0)=pc1(4)/emi
      vj(0)=pc2(4)/emj
      do n=1,3
        vi(n)=pc1(n)/emi
        vj(n)=pc2(n)/emj
        bi(n)=pc1(n)/pc1(4)
        bj(n)=pc2(n)/pc2(4)
      end do

c...V_i^\mu = V_i/rho_B J_i^\mu
        fengi=Ai(0)*vi(0)
        fengj=Aj(0)*vj(0)
        do n=1,3
          fengi = fengi - Ai(n)*vi(n)
          fengj = fengj - Aj(n)*vj(n)
        end do

c...Width parameter.
c     wg = 1.0/(4.d0*pard(104))
      wg = pard(109)
      fskyi=-wg*(gfacts*feng1*t1*sig1 + vfac*fengi)*rhom(i,j)
      fskyj=-wg*(gfacts*feng2*t1*sig2 + vfac*fengj)*rhom(j,i)

c...Derivative of p^\mu_j/m_j in the vector potential.
      dvi=0.0
      dvj=0.0
      if(mstc(111).eq.1) then
        call fenga5(i,pkin1,Ai)
        call fenga5(j,pkin2,Aj)
      do n=1,3
        dvi(n)= ( mstc(116)*Aj(0)*vi(n)/pc1(4) - Aj(n)/emi )*rhom(j,i)
        dvj(n)= ( mstc(116)*Ai(0)*vj(n)/pc2(4) - Ai(n)/emj )*rhom(i,j)
      end do
      endif

      rk(1)=r(1,i)-r(1,j)
      rk(2)=r(2,i)-r(2,j)
      rk(3)=r(3,i)-r(3,j)
      rbi=rk(1)*vi(1)+rk(2)*vi(2)+rk(3)*vi(3)
      rbj=rk(1)*vj(1)+rk(2)*vj(2)+rk(3)*vj(3)
      do n=1,3
        dr2ri(n) = 2*(rk(n)+rbj*vj(n))    ! dR~^2_ij/dR_i
        dr2rj(n) = 2*(rk(n)+rbi*vi(n))    ! dR~^2_ji/dR_i
        dr2pi(n) = 2*rk(n)*rbi/emi
        dr2pj(n) = 2*rk(n)*rbj/emj
      end do

      do n=1,3
        force(n,i)  = force(n,i)  - fskyi*dr2ri(n) - fskyj*dr2rj(n)
        force(n,j)  = force(n,j)  + fskyi*dr2ri(n) + fskyj*dr2rj(n)
        forcer(n,i) = forcer(n,i) + fskyj*dr2pi(n) + dvi(n)
        forcer(n,j) = forcer(n,j) + fskyi*dr2pj(n) + dvj(n)
      end do

      ! Momentum dependent potential off.
      if(mstd(101).eq.0) return

      vf1=vj(0)
      vf2=vi(0)
c...V^mu(p) is fully included.
      if(mstc(115).eq.2) then
        do n=1,3
          vf1 = vf1 - pkin1(n)/pkin1(4)*vj(n)
          vf2 = vf2 - pkin2(n)/pkin2(4)*vi(n)
        end do
      endif
      vf1=vf1*bar
      vf2=vf2*bar

      pk(1)=pc1(1)-pc2(1)
      pk(2)=pc1(2)-pc2(2)
      pk(3)=pc1(3)-pc2(3)
      pk(4)=pc1(4)-pc2(4)
      ppk=pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2
      if(mstc(120).eq.1) then
        pcm(:)=pc1(:)+pc2(:)
        s=pcm(4)**2-pcm(1)**2-pcm(2)**2-pcm(3)**2
        pp2i = ppk + (emi**2- emj**2)**2 / s
        pp2j = pp2i
      else
        pij=pc1(4)*pc2(4)-pc1(1)*pc2(1)-pc1(2)*pc2(2)-pc1(3)*pc2(3)
        pp2i = ppk + (pij-emi**2)**2/emi**2
        pp2j = ppk + (pij-emj**2)**2/emj**2
      endif

      fac1 = 1.d0 + pp2j/pmu1
      fac2 = 1.d0 + pp2j/pmu2
      ! Scalar
      facmom = -feng1*rhom(i,j)
      fmomei=(vex1/pmu1/fac1**2)*facmom
      fmomdi=(vex1/fac1)*wg*facmom
 
      ! Vector
      fmomei=fmomei-(vex2/pmu2/fac2**2)*vf1*rhom(i,j)
      fmomdi=fmomdi-(vex2/fac2)*wg*vf1*rhom(i,j)
      facm2 = vex2/fac2*rhom(i,j)/emj * bar

      fac1 = 1.d0 + pp2i/pmu1
      fac2 = 1.d0 + pp2i/pmu2

      ! Scalar
      facmom = -feng2*rhom(j,i)
      fmomej=(vex1/pmu1/fac1**2)*facmom
      fmomdj=(vex1/fac1)*wg*facmom

      ! Vector
      fmomej=fmomej-(vex2/pmu2/fac2**2)*vf2*rhom(j,i)
      fmomdj=fmomdj-(vex2/fac2)*wg*vf2*rhom(j,i)
      facm1 = vex2/fac2*rhom(j,i)/emi * bar

c...Distance in the two-body cm is used for the argument of MD pot.
      if(mstc(120).eq.1) then

      do n=1,3
        bet(n)=pcm(n)/pcm(4)
        bbi(n)=bet(n)-mstc(116)*bi(n)
        bbj(n)=bet(n)-mstc(116)*bj(n)
      end do
      s=pcm(4)**2-pcm(1)**2-pcm(2)**2-pcm(3)**2
      pma =  (emi**2- emj**2)**2 / s
      do n=1,3
        dpij2pi(n) = pk(n) - mstc(116)*pk(4)*bi(n) + pcm(4)/s*pma*bbi(n)
        dpij2pj(n) =-pk(n) + mstc(116)*pk(4)*bj(n) + pcm(4)/s*pma*bbj(n)
      enddo
      dpji2pi(:)=dpij2pi(:)
      dpji2pj(:)=dpij2pj(:)

      else

      bb1=2*(pij-emi**2)/emi**2
      bb2=2*(pij-emj**2)/emj**2
      do n=1,3
        bb=mstc(116)*pc2(4)*bi(n)-pc2(n)
        dpij2pi(n) = 2*pk(n) - 2*mstc(116)*pk(4)*bi(n) + bb*bb2 ! dpij^2/dp_i
        dpji2pi(n) = 2*pk(n) - 2*mstc(116)*pk(4)*bi(n) + bb*bb1 ! dpji^2/pd_i

        bb=mstc(116)*pc1(4)*bj(n)-pc1(n)
        dpij2pj(n) = -2*pk(n) + 2*mstc(116)*pk(4)*bj(n) + bb*bb2 ! dpij^2/dp_j
        dpji2pj(n) = -2*pk(n) + 2*mstc(116)*pk(4)*bj(n) + bb*bb1 ! dpji^2/pd_j
      end do

      endif

      do n=1,3
        forcer(n,i) = forcer(n,i) + fmomei*dpij2pi(n)
     &          + fmomdj*dr2pi(n) + fmomej*dpji2pi(n)
        forcer(n,j) = forcer(n,j) + fmomej*dpji2pj(n)
     &          + fmomdi*dr2pj(n) + fmomei*dpij2pj(n)
        force(n,i) = force(n,i) - fmomdi*dr2ri(n) - fmomdj*dr2rj(n)
        force(n,j) = force(n,j) + fmomdi*dr2ri(n) + fmomdj*dr2rj(n)
      enddo

c...Derivative of  p_i/m_i part.
      if(mstc(111).eq.1.and.mstc(115).eq.2) then
        fm1=-1.0/pkin2(4)*facm1
        fm2=-1.0/pkin1(4)*facm2
        fm3=mstc(116)*facm1
        fm4=mstc(116)*facm2
        do n=1,3
          forcer(n,i) = forcer(n,i) + fm1*pkin2(n) + bi(n)*fm3
          forcer(n,j) = forcer(n,j) + fm2*pkin1(n) + bj(n)*fm4
        enddo
      endif


      end

**********************************************************************
      subroutine jamvpotdot5

c...Compute time-derivatives of the vector potential.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      real* 8 Ai(0:3),Aj(0:3),Bi(0:3),Bj(0:3),vi(0:3),vj(0:3)
      integer i,j,n,nn
      real*8 pcm(5),rk(5),pk(5),bbi(5),bbj(5),emi,emj
      real*8 s,rbij,dr2ri(3),dr2rj(3),dr2pi(3),dr2pj(3)
      real*8 drji2ri(3),drji2rj(3),drji2pi(3),drji2pj(3)
      real*8 wg,dvi,dvj,vvi,vvj,vij,vji,doti,dotj,fai,faj,fvi,fvj
      real*8 xdoti,pdoti,xdotj,pdotj,p1(5),p2(5)
      real*8 vdott(3),fai2,faj2,fai3,faj3,fvi2,fvj2
      real*8 dgami(3),dgamj(3),rbi,rbj,b1(5),b2(5)
      real*8 fmomdi,fmomdj,dp2ijpi(3),dp2ijpj(3),dp2jipi(3),dp2jipj(3)
      real*8 bb1,bb2,bb,p2i,p2j,psq1,psq2,fac1,fac2,fmomei,fmomej
      real*8 pma,pij,vk1(5),vk2(5)

      real*8 vdot
      common/jamvdot/vdot(3,mxv)

c...Width parameter.
c     wg = 1.0/(4.d0*pard(104))
      wg = pard(109)
      vdot(:,:)=0.0d0

      do i=1,nv
        if(MF_on(i).eq.0) cycle

        vvi = t2 + t3*rhog(i)  ! V_i/rho_i
        dvi=0.0
        if(abs(rhoi(i)).gt.1d-7) ! del(V_i/rho_i)/rho_i
     &  dvi = (pard(113)-1)*t3*rhoi(i)**(pard(113)-3)
        Bi(:)=dvi*rhoj(:,i)

        p1=p(:,i)
        vk1(:)=p1(:)/p1(4)
        if(mstc(112).eq.1) call getmom5(i,p1)
        p1(5)=sqrt(p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2)
        emi=p1(4)
        if(mstc(114).eq.3.or.mstc(113).eq.2) emi=p1(5)

        vi(0)=p1(4)/emi
        vi(1)=p1(1)/emi
        vi(2)=p1(2)/emi
        vi(3)=p1(3)/emi
        b1(:)=p1(:)/p1(4)

      do j=i+1,nv
        if(MF_on(j)==0) cycle

        vvj = t2 + t3*rhog(j)
        dvj=0.0
        if(abs(rhoi(j)).gt.1d-7)
     &  dvj = (pard(113)-1)*t3*rhoi(j)**(pard(113)-3)
        Bj(:)=dvj*rhoj(:,j)

        p2=p(:,j)
        vk2(:)=p2(:)/p2(4)
        if(mstc(112).eq.1) call getmom5(j,p2)
        p2(5)=sqrt(p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2)
        emj=p2(4)
        if(mstc(114).eq.3.or.mstc(113).eq.2) emj=p2(5)

        vj(0)=p2(4)/emj
        vj(1)=p2(1)/emj
        vj(2)=p2(2)/emj
        vj(3)=p2(3)/emj
        b2(:)=p2(:)/p2(4)

        vji=vj(0)*rhoj(0,i)
        vij=vi(0)*rhoj(0,j)
        do n=1,3
          vji=vji-vj(n)*rhoj(n,i)
          vij=vij-vi(n)*rhoj(n,j)
        end do
        Ai(:)=vji*Bi(:)
        Aj(:)=vij*Bj(:)

c...Non-relativistic.
        pk(:)=p1(:)-p2(:)
        rk(:)=r(:,i)-r(:,j)
        dgami(:)=0.0
        dgamj(:)=0.0
        dr2ri(1)= 2*rk(1)
        dr2ri(2)= 2*rk(2)
        dr2ri(3)= 2*rk(3)
        dr2rj(:)   = -dr2ri(:)
        drji2ri(:) =  dr2ri(:)
        drji2rj(:) =  dr2rj(:)
        dr2pi(:)=0.0
        dr2pj(:)=0.0
        drji2pi(:)=0.0
        drji2pj(:)=0.0

c...Compute derivatives of interaction densities.
      if(mstc(113).eq.1) then

      pcm(:)=p1(:)+p2(:)
      bbi(:)=pcm(:)/pcm(4)-mstc(116)*b1(:)
      bbj(:)=pcm(:)/pcm(4)-mstc(116)*b2(:)
      s=pcm(4)**2-pcm(1)**2-pcm(2)**2-pcm(3)**2
      rbij=(rk(1)*pcm(1)+rk(2)*pcm(2)+rk(3)*pcm(3))/s
      do n=1,3
        dr2ri(n)= 2*(rk(n)+rbij*pcm(n))               ! dR~^2_ij/dR_i
        dr2rj(n)= -dr2ri(n)
        dr2pi(n)=2*(rk(n) + pcm(4)*rbij*bbi(n))*rbij  ! dR~^2_ij/dP_i 
        dr2pj(n)=2*(rk(n) + pcm(4)*rbij*bbj(n))*rbij
      enddo
        drji2ri(:)= dr2ri(:)
        drji2rj(:)= -drji2ri(:)
        drji2pi(:)=dr2pi(:)
        drji2pj(:)=dr2pj(:)

c...derivatives from gamma_{ij}.
      do n=1,3
        dgami(n)=mstc(116)*(1.0/pcm(4)-pcm(4)/s)*b1(n)+pcm(n)/s
        dgamj(n)=mstc(116)*(1.0/pcm(4)-pcm(4)/s)*b2(n)+pcm(n)/s
      end do

c...Rest frame of particle.
      else if(mstc(113).eq.2) then

        rbi=rk(1)*vi(1)+rk(2)*vi(2)+rk(3)*vi(3)
        rbj=rk(1)*vj(1)+rk(2)*vj(2)+rk(3)*vj(3)
        do n=1,3
          dr2ri(n) =  2*(rk(n)+rbj*vj(n) )   ! 1/2*dR~^2_ij/dR_i
          dr2rj(n) =  -dr2ri(n)              ! 1/2*dR~^2_ij/dR_j
          dr2pi(n) =  0.0                    ! 1/2*dR~^2_ij/dP_i
          dr2pj(n) =  2*rk(n)*rbj/p2(5)     ! 1/2*dR~^2_ij/dP_j

          drji2rj(n) = 2*(rk(n)+rbi*vi(n))
          drji2rj(n) =  -drji2ri(n)
          drji2pi(n) =  2*rk(n)*rbi/p1(5)
          drji2pj(n) =  0.0
        end do

      endif

      xdoti=0.0
      pdoti=0.0
      xdotj=0.0
      pdotj=0.0
      do n=1,3
        xdoti=xdoti + (vk1(n)+forcer(n,i))*dr2ri(n) 
     &              + (vk2(n)+forcer(n,j))*dr2rj(n)
        xdotj=xdotj + (vk2(n)+forcer(n,j))*drji2rj(n)
     &              + (vk1(n)+forcer(n,i))*drji2ri(n)
        pdoti=pdoti + force(n,i)*(dr2pi(n)-dgami(n)/wg)
     &              + force(n,j)*(dr2pj(n)-dgamj(n)/wg)
        pdotj=pdotj + force(n,j)*(drji2pj(n)-dgamj(n)/wg)
     &              + force(n,i)*(drji2pi(n)-dgami(n)/wg)
      end do
      doti=-wg*(xdoti+pdoti)*rhom(i,j)
      dotj=-wg*(xdotj+pdotj)*rhom(j,i)
      do n=1,3
        vdot(n,i) = vdot(n,i) + doti*(Ai(n)+vvi*vj(n))
        vdot(n,j) = vdot(n,j) + dotj*(Aj(n)+vvj*vi(n))
      end do

c------------------------------------------------------------------------
      fmomdi=0.0
      fmomdj=0.0
      ! Momentum dependent potential on.
      if(mstd(101).ne.0) then

      ! Non-relativistic.
      if(mstc(113).eq.0) then
        psq1 = pk(1)**2 + pk(2)**2 + pk(3)**2
        psq2 = psq1
        do n=1,3
        dp2ijpi(n)=  pk(n)
        dp2ijpj(n)= -pk(n)
        end do
        dp2jipi(:) = dp2ijpi(:)
        dp2jipj(:) = dp2ijpj(:)

      ! Two-body C.M.
      else if(mstc(113).eq.1.or.mstc(120).eq.1) then

        pcm(:)=p1(:)+p2(:)
        bbi(:)=pcm(:)/pcm(4)-mstc(116)*b1(:)
        bbj(:)=pcm(:)/pcm(4)-mstc(116)*b2(:)
        s=pcm(4)**2-pcm(1)**2-pcm(2)**2-pcm(3)**2
        pma = (p1(5)**2- p2(5)**2)**2 / s
        psq1 = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + pma
        psq2 = psq1
        do n=1,3
        dp2ijpi(n) = pk(n) - mstc(116)*pk(4)*b1(n) + pcm(4)/s*pma*bbi(n)
        dp2ijpj(n) =-pk(n) + mstc(116)*pk(4)*b2(n) + pcm(4)/s*pma*bbj(n)
        end do
        dp2jipi(:)=dp2ijpi(:)
        dp2jipj(:)=dp2ijpj(:)

      ! Rest-frame of particle.
      else
        pij=p1(4)*p2(4)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
        p2i=(pij-p1(5)**2)**2/p1(5)**2
        p2j=(pij-p2(5)**2)**2/p2(5)**2
        psq2 = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + p2j
        psq1 = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + p2i
        bb1=2*(pij-p1(5)**2)/p1(5)**2
        bb2=2*(pij-p2(5)**2)/p2(5)**2
        do n=1,3
          bb=mstc(116)*p2(4)*b1(n)-p2(n)
          dp2ijpi(n) = 2*pk(n) - 2*mstc(116)*pk(4)*b1(n) + bb*bb2 ! dpij^2/dp_i
          dp2jipi(n) = 2*pk(n) - 2*mstc(116)*pk(4)*b1(n) + bb*bb1 ! dpji^2/pd_i

          bb=mstc(116)*p1(4)*b2(n)-p1(n)
          dp2ijpj(n) = -2*pk(n) + 2*mstc(116)*pk(4)*b2(n) + bb*bb2 ! dpij^2/dp_j
          dp2jipj(n) = -2*pk(n) + 2*mstc(116)*pk(4)*b2(n) + bb*bb1 ! dpji^2/pd_j
        end do
      endif

c...Derivative of rho_{ij} and gamma_{ij} part.
      fac1 = 1.d0 + psq1/pmu2
      fac2=  1.d0 + psq2/pmu2
      fmomdi = vex2/fac2
      fmomdj = vex2/fac1
      do n=1,3
        vdot(n,i) = vdot(n,i) + doti*fmomdi*vj(n)
        vdot(n,j) = vdot(n,j) + dotj*fmomdj*vi(n)
      end do

c...Derivative of D(p_{ij}) part.
      pdoti=0.0
      pdotj=0.0
      do n=1,3
        pdoti=pdoti + force(n,i)*dp2ijpi(n) + force(n,j)*dp2ijpj(n)
        pdotj=pdotj + force(n,j)*dp2jipj(n) + force(n,i)*dp2jipi(n)
      end do
      fmomei = -pdoti*(vex2/pmu2/fac2**2)*rhom(i,j)
      fmomej = -pdotj*(vex2/pmu2/fac1**2)*rhom(j,i)
      do n=1,3
        vdot(n,i) = vdot(n,i) + fmomei*vj(n)
        vdot(n,j) = vdot(n,j) + fmomej*vi(n)
      end do

      endif  ! End momentum dependent potential part.
c------------------------------------------------------------------------

      if(mstc(111).eq.0) cycle

c...Compute derivatives of pre-factors v_j^\mu in front of vector potential.
      fai=0.0
      faj=0.0
      fai2=0.0
      faj2=0.0
      fai3=0.0
      faj3=0.0
      do n=1,3
        fai=fai - force(n,j)*Bi(n)/emj
        faj=faj - force(n,i)*Bj(n)/emi
        fai2=fai2 + force(n,j)*vj(n)
        faj2=faj2 + force(n,i)*vi(n)
        fai3=fai3 + vj(n)*Bi(n)/p2(4)
        faj3=faj3 + vi(n)*Bj(n)/p1(4)
      end do

      fai=fai*rhom(i,j)
      faj=faj*rhom(j,i)

      if(mstc(114).eq.3.or.mstc(113).eq.2) then
        fai3=Bi(0)/p2(4)
        faj3=Bj(0)/p1(4)
        fvi2=0.0
        fvj2=0.0
      else
        fvi2=-mstc(116)*fai2*rhom(i,j)*vvi/p2(4)
        fvj2=-mstc(116)*faj2*rhom(j,i)*vvj/p1(4)
      endif

      fai=fai+mstc(116)*fai2*fai3*rhom(i,j)
      faj=faj+mstc(116)*faj2*faj3*rhom(j,i)

      fvi=mstc(116)*vvi*rhom(i,j)/emj
      fvj=mstc(116)*vvj*rhom(j,i)/emi


      do n=1,3
        vdot(n,i) = vdot(n,i) + fai*rhoj(n,i)
     &                        + fvi*force(n,j)
     &                        + fvi2*vj(n)

        vdot(n,j) = vdot(n,j) + faj*rhoj(n,j)
     &                        + fvj*force(n,i)
     &                        + fvj2*vi(n)
      end do

c---------------------------------------------------------------------------
      ! Momentum dependent potential on.
      if(mstd(101).ne.0) then
        do n=1,3
          vdot(n,i) = vdot(n,i) + fmomdi*rhom(i,j)*force(n,j)/emj
          vdot(n,j) = vdot(n,j) + fmomdj*rhom(j,i)*force(n,i)/emi
        end do
        if(mstc(114).ne.3.and.mstc(113).ne.2) then
          do n=1,3
        vdot(n,i)=vdot(n,i)-mstc(116)*fai2*fmomdi*rhom(i,j)/p2(4)*vj(n)
        vdot(n,j)=vdot(n,j)-mstc(116)*faj2*fmomdj*rhom(j,i)/p1(4)*vi(n)
          end do
        endif
      endif
c---------------------------------------------------------------------------

      end do
      end do

      if(mstc(124).eq.1) return

      vdott=0.0
      nn=0
      do i=1,nv
       if(MF_on(i).eq.0) cycle
       vdott(:) = vdott(:) + vdot(:,i)
       nn=nn+1
      end do
      if(nn.ge.1) then
      do i=1,nv
       if(MF_on(i).eq.0) cycle
       vdot(:,i) = vdot(:,i) - vdott(:)/nn
      end do
      endif

      end

**********************************************************************
      subroutine kineticmomit5
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,it
      real*8 diff0,diff,vv,pp
      real*8 pcan(5,mxv)

c....Save canonical momenta.
      pcan=p
      diff0=0.0
      do i=1,nv
        if(MF_on(i).eq.0) cycle
        diff0=diff0+sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+p(4,i)**2)
      end do

      do it=1,100

      call jamrqmm5
      diff=0.0
      do i=1,nv
        if(MF_on(i).eq.0) cycle
        vv=k(9,i)/3*(t2+t3*rhog(i))
        pots(i) = t1*rho(i) 
        p(1,i)=pcan(1,i) - vv*rhoj(1,i)
        p(2,i)=pcan(2,i) - vv*rhoj(2,i)
        p(3,i)=pcan(3,i) - vv*rhoj(3,i)
c       em=p(5,i)+pots(i)
        pp=p(1,i)**2+p(2,i)**2+p(3,i)**2
        p(4,i)=sqrt(p(5,i)**2+pp)
        diff=diff+sqrt(pp+p(4,i)**2)
      end do
        print *,it,'diff=',diff0-diff
      if(abs(diff0-diff).lt.1d-5) return
      diff0=diff

      end do

      print *,'kinetic momentum does not converge',diff0,diff
      stop

      end

**********************************************************************

      subroutine jamrqmm5

c...Purpose: to prepare matrix in calculating force
c...Option mstc(113)=2 is not implemented for Coulomb, Yukawa Symmetry force 

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j
      real*8 facg,wg,emisq,emjsq,fi,fj,bi,bj,bij
      real*8 px,py,pz,pe,rx,ry,rz,s,r2,ps,psq,rs
      real*8 qfacr,p1(5),p2(5),v1(5),v2(5),emi,emj
      real*8 gfac,pmom2ijs,pmom2jis,pmom2ijv,pmom2jiv
      real*8 gfacs,gfacs1,gfacts,vfac,vfac1,vfacs
      real*8 den,gam,r2i,r2j,p2i,p2j,pij

c     fac  =(4.d0*paru(1)*pard(104))**1.5d0 !    [(4*pi*L)^3/2]
c     wg   = 1.0/4.d0/pard(104)
      wg = pard(109)
      facg = pard(110)

      rho=0.d0   !...  <rho_i> = sigma_j(not= i) rho_ij
      rhoi=0.d0
      vmoms=0.0
      vmom=0.0
      rhoj=0.d0
      rhom=0.d0

      rhoc=0.d0  ! Coulomb
      rhc=0.d0   ! Coulomb
      rhos=0.d0  ! Symmetry energy
      rhs=0.d0   ! Symmetry energy
      rhoy=0.d0  ! Yukawa potential
      rhy=0.d0   ! Yukawa potential

      do 100 i=1,nv             ! sum for i-th particle (sigma_i) 

        if(MF_on(i).eq.0) goto 100
        bi=k(9,i)/3
        p1(:)=p(:,i)
        if(mstc(112).eq.1)  call getmom5(i,p1)
        emisq=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
        emi=sqrt(emisq)
        if(mstc(113).ne.2) then
          fi=emi/p1(4)
          v1(:)=p1(:)/p1(4)
        else
          fi=1.0
          v1(:)=p1(:)/emi
        endif
        gfacs1=gfacs(i)
        vfac1=vfacs(i)

        do 110 j = i+1 , nv      ! Sum for j (sigma_j(>i))

          if(MF_on(j).eq.0) goto 110
          bj=k(9,j)/3
          bij=bi*bj
c         bij=1.0
          p2(:)=p(:,j)
          if(mstc(112).eq.1) call getmom5(j,p2)
          emjsq = p2(4)**2 - p2(1)**2 - p2(2)**2 - p2(3)**2
          emj=sqrt(emjsq)
          if(mstc(113).ne.2) then
            fj=emj/p2(4)
            v2(:)=p2(:)/p2(4)
          else
            fj=1.0
            v2(:)=p2(:)/emj
          endif

          rx = r(1,i) - r(1,j)
          ry = r(2,i) - r(2,j)
          rz = r(3,i) - r(3,j)
          ps=(p1(1)-p2(1))**2+(p1(2)-p2(2))**2+(p1(3)-p2(3))**2
          rs= rx**2 + ry**2 + rz**2
          gfac=qfacr(i)*qfacr(j)
          gfacts=gfacs1*gfacs(j)
          vfac=vfac1*vfacs(j)

c...Non-relativistic.
          if(mstc(113).eq.0) then

          rhom(i,j) = facg * exp(-rs*wg) * gfac
          rhom(j,i) = rhom(i,j)
          pmom2ijs = vex1/(1.d0+ps/pmu1)
          pmom2ijv = vex2/(1.d0+ps/pmu2)
          pmom2jis = pmom2ijs
          pmom2jiv = pmom2ijv
          p2i=ps
          p2j=ps

c...Two body distance is defined by the two particle c.m. frame.
          else if(mstc(113).eq.1) then

          px = p1(1) + p2(1)
          py = p1(2) + p2(2)
          pz = p1(3) + p2(3)
          pe = p1(4) + p2(4)
          s = pe**2 - px**2 - py**2 - pz**2
          r2= rs + (rx*px + ry*py + rz*pz)**2/s
          den  = facg * exp(-r2*wg) * gfac
          gam = pe/sqrt(s)
          rhom(i,j) = den*gam
          rhom(j,i) = den*gam

          psq = ps - (p1(4)-p2(4))**2 + (emisq- emjsq)**2/s
          pmom2ijs = vex1/(1.d0+psq/pmu1)
          pmom2ijv = vex2/(1.d0+psq/pmu2)
          pmom2jis = pmom2ijs
          pmom2jiv = pmom2ijv
          p2i=psq
          p2j=psq

c...distance is measured from the rest-frame of particle j.
          else

          r2i=rs + (rx*p(1,i)+ry*p(2,i)+rz*p(3,i))**2/emisq
          r2j=rs + (rx*p(1,j)+ry*p(2,j)+rz*p(3,j))**2/emjsq
          rhom(i,j)= facg * exp(-r2j*wg) * gfac
          rhom(j,i)= facg * exp(-r2i*wg) * gfac

c...Distance in the two-body cm is used for the argument of MD pot.
          if(mstc(120).eq.1) then
            s=(p1(4)+p2(4))**2
     &          -(p1(1)+p2(1))**2-(p1(2)+p2(2))**2-(p1(3)+p2(3))**2
            p2i = ps - (p1(4)-p2(4))**2 + (emisq- emjsq)**2/s
            p2j=p2i
          else
            pij=p1(4)*p2(4)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
            p2i=ps - (p1(4)-p2(4))**2 + (pij-emisq)**2/emisq
            p2j=ps - (p1(4)-p2(4))**2 + (pij-emjsq)**2/emjsq
          endif
          pmom2ijs = vex1/(1.d0+p2j/pmu1)
          pmom2jis = vex1/(1.d0+p2i/pmu1)
          pmom2ijv = vex2/(1.d0+p2j/pmu2)
          pmom2jiv = vex2/(1.d0+p2i/pmu2)

          endif

          if(mstc(118).ge.1) then
            rhom(i,j)=rhom(i,j)/(1.0+p2j/parc(105)**2)
            rhom(j,i)=rhom(j,i)/(1.0+p2i/parc(105)**2)
          endif

          ! Scalar density.
          rho(i) = rho(i) + rhom(i,j)*fj*gfacts
          rho(j) = rho(j) + rhom(j,i)*fi*gfacts

          ! Baryon current.
          rhoj(0,i) = rhoj(0,i) + rhom(i,j)*v2(4)*vfac*bj
          rhoj(1,i) = rhoj(1,i) + rhom(i,j)*v2(1)*vfac*bj
          rhoj(2,i) = rhoj(2,i) + rhom(i,j)*v2(2)*vfac*bj
          rhoj(3,i) = rhoj(3,i) + rhom(i,j)*v2(3)*vfac*bj
          rhoj(0,j) = rhoj(0,j) + rhom(j,i)*v1(4)*vfac*bi
          rhoj(1,j) = rhoj(1,j) + rhom(j,i)*v1(1)*vfac*bi
          rhoj(2,j) = rhoj(2,j) + rhom(j,i)*v1(2)*vfac*bi
          rhoj(3,j) = rhoj(3,j) + rhom(j,i)*v1(3)*vfac*bi

          ! Momentum dependent potential off.
          if(mstd(101).eq.0) goto 110

          vmoms(i) = vmoms(i) + pmom2ijs*rhom(i,j)*fj*vfac
          vmoms(j) = vmoms(j) + pmom2jis*rhom(j,i)*fi*vfac

          vmom(0,i) = vmom(0,i) + pmom2ijv*rhom(i,j)*v2(4)*vfac*bij
          vmom(1,i) = vmom(1,i) + pmom2ijv*rhom(i,j)*v2(1)*vfac*bij
          vmom(2,i) = vmom(2,i) + pmom2ijv*rhom(i,j)*v2(2)*vfac*bij
          vmom(3,i) = vmom(3,i) + pmom2ijv*rhom(i,j)*v2(3)*vfac*bij
          vmom(0,j) = vmom(0,j) + pmom2jiv*rhom(j,i)*v1(4)*vfac*bij
          vmom(1,j) = vmom(1,j) + pmom2jiv*rhom(j,i)*v1(1)*vfac*bij
          vmom(2,j) = vmom(2,j) + pmom2jiv*rhom(j,i)*v1(2)*vfac*bij
          vmom(3,j) = vmom(3,j) + pmom2jiv*rhom(j,i)*v1(3)*vfac*bij

 110     continue               ! Loop end of sigma_j(>i)
 100  continue                  ! Loop end of sigma_i
 
c...Compute invariant baryon density.
      if(mstc(114).ne.0) then
        do i=1,nv
          if(MF_on(i).eq.0) cycle
          rhoi(i)=sqrt(max(0d0,rhoj(0,i)**2
     &        -rhoj(1,i)**2-rhoj(2,i)**2-rhoj(3,i)**2))
        enddo
      endif

      end

c***********************************************************************
      subroutine sigmafield5
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,it,maxit
      real*8 diff,hc,pot0
      real*8 gfacs,dens,gs,gss
      real*8 sigma1,sigma_bisec5

      hc=paru(3)
      gs=pard(101)        ! Scalar part g_sigma

      maxit=5
      if(mstc(112).eq.1.or.mstc(119).eq.0) maxit=1

      do it=1,maxit

        pot0=sum(pots)
        do i=1,nv  ! Loop over particles.
        if(MF_on(i)==0) cycle

        dens=gs*rho(i)/gfacs(i)  ! factor gr is already contained in rho(i) 
        gss=t1*gfacs(i)*hc

        sigma1=abs((p(5,i) + vmoms(i))/gss)
        pots(i)=gss*sigma_bisec5(dens,sigma1) + vmoms(i)
        p(4,i)=sqrt((p(5,i)+pots(i))**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)

        end do ! end particle loop

        diff=abs(pot0-sum(pots))
        if(diff.lt.1d-5) goto 100

c       do i=1,nv
c         if(MF_on(i)==0) cycle
c         em=p(5,i)+pots(i)
c         p(4,i)=sqrt(em**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
c       end do

        call jamrqmm5

      end do  ! end iteration loop
c     write(6,*)'sigma does not converge',diff
100   continue
c     write(6,*)'sigma ',diff,it


      end

c***********************************************************************
      real*8 function sigma_bisec5(dens,sigma1)
c...find sigma field in 1/fm
      implicit none
      include 'jam2.inc'
      integer it
      real*8 ms2,sigma0,sigma1,sigma,g2,g3,g4,dens,f
      real*8 func,x,f0,f1,sig0
c     func(x)= -dens + ms2*x + g2*x**2 + g3*x**3
      func(x)= -dens*(1.0+0.5*g4*x*x)**2
     &    + ms2*x + g2*x*x + g3*x*x*x
     &    + 1.0/6.0*g2*g4*x**4+0.25*g3*g4*x**5

      ms2=pard(112)**2    ! sigma mass (1/fm^2)
      g2=pard(102)     ! Scalar part g2 (1/fm)
      g3=pard(103)    ! Scalar part g3
      g4=pard(116)    ! Scalar part g4

      sigma0=0.0
      f0=func(sigma0)
      f1=func(sigma1)
      sig0=sigma1
c...Look for a new initial condition.
      if(f0*f1 > 0.0) then
        do while(sigma1 .gt. sigma0)
        sigma1 = sigma1 - 0.01
        if(sigma1 < sigma0) then
          print *,'sigma_bisec5::bisec no solution sigma1=',sigma1,sig0
          sigma_bisec5=sigma1-1e-5
          stop
          return
        endif
        f1=func(sigma1)
        if(f0*f1 .lt.0.0) exit
        end do
      endif

      if(f0*f1 > 0.0 .or. sigma1 < 0.0) then
          print *,'bisec no solution sigma=',sigma0,sigma,dens
          print *,'func1 func2',func(sigma0),func(sigma1)
          sigma_bisec5=0.0
          stop
          return
      endif

c...bisection method starts.
      do it=1,50
        sigma=0.5*(sigma0+sigma1)
        f=  func(sigma)
        if(f.gt.0.0) then
          sigma1=sigma
        else
          sigma0=sigma
        endif
        if(abs(f).lt.1d-5) goto 100
        end do
        print *,'does not converge make_sigma_table',it,sigma
        stop
100   continue

      sigma_bisec5=sigma

      end

c***********************************************************************

      subroutine jamepart5(difm,difv)  ! energy

c...Purpose: to calculate single particle potential energy in RQMD/S
c H = sigma_i=1^N  1/(2*E_i) *[E_i^2 -vec{p}_i^2 - m_i^2 - 2m_i*V_i ]
c here V_i = t1 * <rho_i> + t3 *<rho_i>^gamma : Skyrme Potential
c...last revised: 2015/2/20
c...last revised: 2019/3/6
c...last revised: 2019/12/16

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
c...Local Variable
      integer i,nn,nn2
      real*8 difm,difv,bar,emf,pp2,em0,vv,potv0(0:3),t3p,t4p
      real*8 vdott(3),vdot,vtot(0:3),vmtot(0:3)
      common/jamvdot/vdot(3,mxv)

c....Sigma-omega model. spots() is also computed.
      if(mstd(102).eq.1) call sigmafield5

      difm=0.d0   ! Difference of effective mass
      difv=0.d0   ! Difference of vector potential
      vdott=0.0
      vtot=0.0
      vmtot=0.0
      nn=0
      nn2=0

      do 100 i=1,nv

       if(MF_on(i)==0) then
         pots(i)=0.0d0
         potv(:,i)=0.0d0
         vdot(:,i)=0.0d0
         goto 100
       endif

       bar = k(9,i)/3
       potv0(:)=potv(:,i)
       potv(:,i)=0.0
       pp2 = p(1,i)**2 + p(2,i)**2 + p(3,i)**2
       emf=p(5,i)+pots(i)

       if(mstd(102).eq.0) then
         pots(i) = t1*rho(i) + vmoms(i)
         if(p(5,i)+pots(i).le.0d0) then
           print *,'effective mass=',p(5,i)+pots(i)
           pots(i)=-p(5,i)+0.001
         endif
         emf=p(5,i)+pots(i)
c        p(4,i)=sqrt(emf**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
       endif

       ! Option for vector potential
       if(mstc(115).ge.1.and.rhoi(i).gt.0d0) then
           t3p=0.0
           t4p=0.0
           if(t3.ne.0d0) t3p=t3*rhoi(i)**(pard(113)-1.0)
           if(t4.ne.0d0) t4p=t4*rhoi(i)**(pard(115)-1.0)
           vv = bar*(t2 + t3p + t4p)
           if(mstc(115).eq.2) then ! Fully included
             potv(:,i)=vv*rhoj(:,i) + vmom(:,i)
             vtot(:) = vtot(:) + potv(:,i)
             vmtot(:) = vmtot(:) + vmom(:,i)
             nn2 = nn2 + 1
           else
             potv(0,i)=vv*rhoj(0,i) + vmom(0,i)
           endif

c...Only  zero-th  component of the vector potential is included.
       else
         t3p=0.0d0
         t4p=0.0d0
         if(t3.ne.0d0) t3p=t3*rhoi(i)**pard(113)
         if(t4.ne.0d0) t4p=t4*rhoi(i)**pard(115)
         potv(0,i) = bar*(t2*rhoi(i) +  t3p + t4p) + vmom(0,i)
       endif

       em0=sqrt(max(0d0,p(4,i)**2 - pp2))
       difm=difm+abs(em0-emf)
       difv=difv+sqrt((potv0(0)-potv(0,i))**2+(potv0(1)-potv(1,i))**2
     &    +(potv0(2)-potv(2,i))**2+(potv0(3)-potv(3,i))**2)

       nn = nn + 1
       if(mstc(119).ge.2) then
       vdot(1,i) = potv(1,i) - potv0(1)
       vdot(2,i) = potv(2,i) - potv0(2)
       vdot(3,i) = potv(3,i) - potv0(3)
       vdott(:) = vdott(:) + vdot(:,i)
       endif

100   continue

      
c     if(nn2.ge.1) then
c     do i=1,nv
c      if(MF_on(i).eq.0) cycle
c       potv(1,i)=potv(1,i) - vtot(1)/nn2
c       potv(2,i)=potv(2,i) - vtot(2)/nn2
c       potv(3,i)=potv(3,i) - vtot(3)/nn2
c       vmom(1,i)=vmom(1,i) - vmtot(1)/nn2
c       vmom(2,i)=vmom(2,i) - vmtot(2)/nn2
c       vmom(3,i)=vmom(3,i) - vmtot(3)/nn2
c     end do
c     endif

      if(mstc(119).ge.2.and.nn.ge.1) then
      vdott(:) = vdott(:)/nn
      do i=1,nv
       if(MF_on(i).eq.0) cycle
       vdot(:,i) = ( vdot(:,i) - vdott(:) ) / parc(2)
      end do
      endif

      mste(44)=1
      end

c***********************************************************************

      subroutine recovere5all

c...Purpose: to recover the energy conservation by changing all particle
c...momenta.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      integer i,itry
      real*8 ptot(5),psys(4),etot0,etot,vx,vy,vz,gam,eps
      real*8 a,diff,diffv ,efinal,ee,pp2
      parameter(eps=1e-3)

      ptot(:)=0.0
      do i=1,nv
        ptot(:) = ptot(:) + p(:,i)
      end do

      psys(1)=pard(9)
      psys(2)=pard(10)
      psys(3)=pard(11)
      psys(4)=pard(12)
      etot0 = sqrt(psys(4)**2 - psys(1)**2 - psys(2)**2 - psys(3)**2)
      vx=ptot(1)/ptot(4)
      vy=ptot(2)/ptot(4)
      vz=ptot(3)/ptot(4)
      gam=ptot(4)/sqrt(ptot(4)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2)

      if(vx**2+vy**2+vz**2.gt.1e-5) then
c...Go to the C.M. frame of all interacting particles.
      do i=1,nv
       call jamrobo(0d0,0d0,-vx,-vy,-vz,gam,p(1,i),p(2,i),p(3,i),p(4,i))
      end do
      endif

      a = 1.0
      do itry=1,30
        call jamrqmm5
        call jamepart5(diff,diffv)
        etot=0.0
        efinal=0.0
        do i=1,nv
          ee=p(4,i)
          if(MF_on(i).ne.0) then
          if(mstc(119).eq.0) then
            pp2=(p(1,i)-potv(1,i))**2
     &        + (p(2,i)-potv(2,i))**2 + (p(3,i)-potv(3,i))**2
            ee=sqrt((p(5,i)+pots(i))**2 + pp2)
          else
          endif
          endif

          etot = etot + ee + potv(0,i)
          efinal = efinal + p(4,i)
        end do
        if(abs(etot-etot0).le.eps*etot0) exit
        a = etot0/etot*a
        do i=1,nv
          p(1,i)=p(1,i)*a
          p(2,i)=p(2,i)*a
          p(3,i)=p(3,i)*a
          p(4,i)=sqrt((p(5,i)+pots(i))**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
        end do
      end do
      if(itry.eq.30) print *,'a=',a,itry,abs(etot-etot0)/etot0
c     print *,'a=',a,itry,abs(etot-etot0)/etot0

      vx=psys(1)/efinal
      vy=psys(2)/efinal
      vz=psys(3)/efinal
      gam=1.0
      if(vx**2+vy**2+vz**2.gt.1e-5) then
c...Go back to the computational frame.
      do i=1,nv
       call jamrobo(0d0,0d0,vx,vy,vz,gam,p(1,i),p(2,i),p(3,i),p(4,i))
      end do
      endif

      call jamrqmm5
      call jamepart5(diff,diffv)

      end

c***********************************************************************

      subroutine recovere5

c...Purpose: to recover the energy conservation by changing interacting particle
c...momenta with potentials.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      integer i,itry,nn
      real*8 ptot(5),pfree(5),psys(4),etot0,etot,vx,vy,vz,gam,eps
      real*8 a,diff,diffv ,efinal,ee,pp2
      parameter(eps=1e-3)

      ptot(:)=0.0
      pfree(:)=0.0
      nn=0
      do i=1,nv
       if(MF_on(i).eq.0) then
         pfree = pfree + p(:,i)
         if(potv(0,i).ne.0d0) print *,potv(0,i)
       else
         ptot(:) = ptot(:) + p(:,i)
         nn=nn+1
       endif
      end do

      if(nn.eq.0) return

      psys(1)=pard(9) -pfree(1)
      psys(2)=pard(10)-pfree(2)
      psys(3)=pard(11)-pfree(3)
      psys(4)=pard(12)-pfree(4)
      etot0 = sqrt(psys(4)**2 - psys(1)**2 - psys(2)**2 - psys(3)**2)
      vx=ptot(1)/ptot(4)
      vy=ptot(2)/ptot(4)
      vz=ptot(3)/ptot(4)
      gam=ptot(4)/sqrt(ptot(4)**2-ptot(1)**2-ptot(2)**2-ptot(3)**2)

c...Go to the C.M. frame of all interacting particles.
      do i=1,nv
       if(MF_on(i).eq.0) cycle
       call jamrobo(0d0,0d0,-vx,-vy,-vz,gam,p(1,i),p(2,i),p(3,i),p(4,i))
      end do

      a = 1.0
      do itry=1,30
        call jamrqmm5
        call jamepart5(diff,diffv)
        etot=0.0
        efinal=0.0
        do i=1,nv
          if(MF_on(i).eq.0) cycle
          if(mstc(119).eq.0) then
            pp2=(p(1,i)-potv(1,i))**2
     &        + (p(2,i)-potv(2,i))**2 + (p(3,i)-potv(3,i))**2
            ee=sqrt((p(5,i)+pots(i))**2 + pp2)
          else
            ee=p(4,i)
          endif
          etot = etot + ee + potv(0,i)
          efinal = efinal + p(4,i)
        end do
        if(abs(etot-etot0).le.eps*etot0) exit
        a = etot0/etot*a
        do i=1,nv
          if(MF_on(i).eq.0) cycle
          p(1,i)=p(1,i)*a
          p(2,i)=p(2,i)*a
          p(3,i)=p(3,i)*a
          p(4,i)=sqrt((p(5,i)+pots(i))**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
        end do
      end do
      if(itry.eq.30) print *,'a=',a,itry,abs(etot-etot0)/etot0
c     print *,'a=',a,itry,abs(etot-etot0)/etot0

      vx=psys(1)/efinal
      vy=psys(2)/efinal
      vz=psys(3)/efinal
      gam=1.0
c...Go back to the computational frame.
      do i=1,nv
       if(MF_on(i).eq.0) cycle
       call jamrobo(0d0,0d0,vx,vy,vz,gam,p(1,i),p(2,i),p(3,i),p(4,i))
      end do
      call jamrqmm5
      call jamepart5(diff,diffv)

      end

**********************************************************************

      subroutine jamsvpot5(difm,difv)

c...Purpose: to compute scalar and vector potentials self-consistently.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      integer i,it,iopt
      real*8 difm,difv,pkin(3,mxv)
c     logical jamrqpb

      iopt=1
      if(mstc(112).eq.1) iopt=0

c     do i=1,nv
c        MF_on(i)=1
c        if(jamrqpb(i)) MF_on(i)=0
c     enddo

c     pots(:)=0.0d0
c     potv(:,:)=0.0d0

c...Initial condition.
c     call jamrqmm2
c     call jamepart2(difm,difv)

      pkin(1,:)=p(1,:)
      pkin(2,:)=p(2,:)
      pkin(3,:)=p(3,:)

      do it=1,50
        if(difm+difv.le.1d-5) goto 3000

        if(mstc(112).eq.1) then
          do i=1,nv
            if(MF_on(i)==0) cycle
             p(1,i)=pkin(1,i) + potv(1,i)
             p(2,i)=pkin(2,i) + potv(2,i)
             p(3,i)=pkin(3,i) + potv(3,i)
             p(4,i)=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
          end do
         endif

         call jamrqmm5
         call jamepart5(difm,difv)
         call setenergy(iopt)
      end do
      write(6,800)it,difm,difv,difm+difv
3000  continue

800   format('jamsvpot5: not conserved diff',i3,3(e15.8,1x))
      end

************************************************************************
      subroutine rqmdpotparam5
      implicit none
      real*8 hc,rho0,conv
      real*8 cs,cv,cd,ce,ms,mv,mn,gs,gv,g1,g2,g3,g4
      real*8 gsp,gvp,mu1,mu2
      include 'jam2.inc'
      include 'jam4.inc'

      hc=paru(3)
      rho0=parc(21)
      conv=hc*hc*hc

      pard(105)=1.0  ! pmu1 dummy
      pard(106)=1.0  ! pmu2 dummy
      pard(107)=0.0d0   ! vex1
      pard(108)=0.0d0   ! vex2
      mstd(101)=0  ! momentum dependent potential is not used.
      mstd(102)=0  ! nonlinear sigma-field is not used.
      g4 = 0.0
      pard(116)= g4  ! Scalar part g4

c...mean field with scalar + vector model by 
c...M.I.Gorenstein, D.H.Rischke, H.Stoecker, W.Greiner, J.Phys.G19 (1993) L69
      if (mstc(106).eq.101) then  ! K=553 MeV M*/M=0.546
        pard(101)= -355.1913106253926*conv  ! Scalar part alpha
        pard(102)=  0.0                     ! Scalar part beta
        pard(103)= 1.0/3.0                  ! Scalar part gamma
        pard(111)=  266.591142296637*conv   ! Vector part alpha
        pard(112)= 0.0                      ! Vector part beta
        pard(113)=  1.0/3.0                 ! Vector part gamma
      elseif (mstc(106).eq.102) then  ! K=380 MeV m*/m=0.603
        pard(101)= -306.72053219515277*conv  ! Scalar part alpha
        pard(102)=  0.0                 ! Scalar part beta
        pard(103)= 1.0/3.0              ! Scalar part gamma
        pard(111)=  239.85427189549642*conv   ! Vector part alpha
        pard(112)= -0.1231975684547163*hc   ! Vector part beta
        pard(113)=  1.0/3.0              ! Vector part gamma
      elseif (mstc(106).eq.103) then   ! K=300 MeV m*m=0.639
        pard(101)= -277.65372307188386*conv  ! Scalar part alpha
        pard(102)=  0.0                 ! Scalar part beta
        pard(103)= 1.0/3.0              ! Scalar part gamma
        pard(111)=  221.49524746939045*conv  ! Vector part alpha
        pard(112)= -0.182394720297761*hc   ! Vector part beta
        pard(113)=  1.0/3.0             ! Vector part gamma
      elseif (mstc(106).eq.104) then   ! K=210 MeV m*/m=0.691
        pard(101)= -235.38649645521997*conv  ! Scalar part alpha
        pard(102)=  0.0                 ! Scalar part beta
        pard(103)= 1.0/3.0              ! Scalar part gamma
        pard(111)= 192.23574594599043*conv   ! Vector part alpha
        pard(112)= -0.2524095316796683*hc  ! Vector part beta
        pard(113)=  1.0/3.0             ! Vector part gamma

      elseif (mstc(106).eq.105) then   ! K=170 MeV m*/m=0.722
        pard(101)= -211.08806268066482*conv  ! Scalar part alpha
        pard(102)=  0.0                 ! Scalar part beta
        pard(103)= 1.0/3.0              ! Scalar part gamma
        pard(111)= 174.2489492426502*conv    ! Vector part alpha
        pard(112)= -0.2854076048897029*hc  ! Vector part beta
        pard(113)=  1.0/3.0             ! Vector part gamma

      elseif (mstc(106).eq.106) then   ! K=150 MeV m*/m=0.74
        pard(101)= -196.92126626182198*conv  ! Scalar part alpha
        pard(102)=  0.0                 ! Scalar part beta
        pard(103)= 1.0/3.0              ! Scalar part gamma
        pard(111)= 163.4220206404303*conv    ! Vector part alpha
        pard(112)= -0.30254981062092773*hc  ! Vector part beta
        pard(113)= 1.0/3.0             ! Vector part gamma

      elseif (mstc(106).eq.107) then   ! K=80 MeV m*/m=0.8298
        pard(101)= -127.88090096098931*conv  ! Scalar part alpha
        pard(102)=  0.0                 ! Scalar part beta
        pard(103)= 1.0/3.0              ! Scalar part gamma
        pard(111)= 107.69357981235275*conv    ! Vector part alpha
        pard(112)= -0.36803632376707573*hc  ! Vector part beta
        pard(113)=  1.0/3.0             ! Vector part gamma

      elseif (mstc(106).eq.108) then   ! K=210 MeV m*/m=0.83

        pard(101)= -127.7884489283939*conv  ! Scalar part alpha
        pard(102)=  0.0                 ! Scalar part beta
        pard(103)= 1.0/3.0              ! Scalar part gamma
        pard(111)= 149.60102152417969*conv    ! Vector part alpha
        pard(112)= -2.359080953010444*hc  ! Vector part beta
        pard(113)=  1.0/3.0             ! Vector part gamma
        pard(114)=0.6150484714503821*hc**(3.0/5.0)  ! Vector part beta_2
        pard(115)=  1.0/5.0             ! Vector part gamma_2

      elseif (mstc(106).eq.109) then   ! K=380 MeV m*/m=0.83
        pard(101)= -127.78844892839392*conv  ! Scalar part alpha
        pard(102)=  0.0                 ! Scalar part beta
        pard(103)= 1.0/3.0              ! Scalar part gamma
        pard(111)= 204.47532213760095*conv    ! Vector part alpha
        pard(112)= -4.961285574309068*hc  ! Vector part beta
        pard(113)=  1.0/3.0             ! Vector part gamma
        pard(114)=1.4189173178206818*hc**(3.0/5.0)  ! Vector part beta_2
        pard(115)=  1.0/5.0             ! Vector part gamma_2

      elseif (mstc(106).eq.110) then   ! K=210 MeV m*/m=0.84
        pard(101)=  -120.17905846641514*conv
        pard(102)=  0.0                 ! Scalar part beta
        pard(103)= 1.0/3.0              ! Scalar part gamma
        pard(111)= 144.96973938530886*conv
        pard(112)= -2.448539846239129*hc
        pard(113)=  1.0/3.0             ! Vector part gamma
        pard(114)= 0.6409296393815799*hc**(3.0/5.0)
        pard(115)=  1.0/5.0             ! Vector part gamma_2

      elseif (mstc(106).eq.111) then   ! K=210 MeV m*/m=0.80
        Cv= 162.50463896111825
        Cs=  150.71138537875936
        Cd= -2.0515258045017397
        Ce= 0.5258414437591274
        pard(101)= -Cs*conv
        pard(102)=  0.0                 ! Scalar part beta
        pard(103)= 1.0/3.0              ! Scalar part gamma
        pard(111)= Cv*conv
        pard(112)= Cd*hc
        pard(113)=  1.0/3.0             ! Vector part gamma
        pard(114)= Ce*hc**(3.0/5.0)
        pard(115)=  1.0/5.0             ! Vector part gamma_2

      elseif (mstc(106).eq.112) then   ! K=380 MeV m*/m=0.80
        Cv= 217.37893957453923
        Cs=  150.71138537875936
        Cd= -4.653730425800347
        Ce= 1.3297102901294222
        pard(101)= -Cs*conv
        pard(102)=  0.0                 ! Scalar part beta
        pard(103)= 1.0/3.0              ! Scalar part gamma
        pard(111)= Cv*conv
        pard(112)= Cd*hc
        pard(113)=  1.0/3.0             ! Vector part gamma
        pard(114)= Ce*hc**(3.0/5.0)
        pard(115)=  1.0/5.0             ! Vector part gamma_2



c...Non-linear sigma-omega model.
c...A. Lang, B. Blattel, W. Cassing, V. Koch, U. Mosel, and K. Weber,
c...Z. Phys. A 340 (1991).
      else if (mstc(106).eq.1001) then   ! K=380 MeV NL1 m*/m=0.83
        pard(101)=  6.91               ! Scalar part g_sigma
        pard(102)= -40.6               ! Scalar part g2 (1/fm)
        pard(103)= 384.4               ! Scalar part g3
        pard(111)= 7.54                ! Vector part g_omega
        pard(112)= 2.79                ! sigma mass (1/fm)
        pard(113)= 3.97                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.
      else if (mstc(106).eq.1002) then   ! K=210 MeV NL2 m*/m=0.83
        pard(101)=  8.5                ! Scalar part g_sigma
        pard(102)=  50.37              ! Scalar part g2 (1/fm)
c       pard(102)=  50.57              ! Scalar part g2 (1/fm)
        pard(103)= -6.26               ! Scalar part g3
        pard(111)= 7.54                ! Vector part g_omega
        pard(112)= 2.79                ! sigma mass (1/fm)
        pard(113)= 3.97                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.
      else if (mstc(106).eq.1003) then   ! K=380 MeV NL3 m*/m=0.7
        pard(101)=  9.5                ! Scalar part g_sigma
        pard(102)= 1.589               ! Scalar part g2 (1/fm)
        pard(103)= 34.23               ! Scalar part g3
        pard(111)= 10.95               ! Vector part g_omega
        pard(112)= 2.79                ! sigma mass (1/fm)
        pard(113)= 3.97                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.
c...G.Q.Li PRC49(1994)1139
      else if (mstc(106).eq.1004) then   ! K=380 MeV m*/m=0.83 rho0=0.16
         ms=0.55/hc;
         mv=0.789/hc;
         mn=0.933/hc;
         gs=11.27*ms/mn;
         gv=8.498*mv/mn;
         g1=-0.0283*gs*gs*gs*mn;
         g2=0.1859*gs*gs*gs*gs;
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

c...G.Q.Li PRC49(1994)1139
      else if (mstc(106).eq.1005) then   ! K=200 MeV m*/m=0.83 rho0=0.16
         ms=0.55/hc;
         mv=0.789/hc;
         mn=0.933/hc;
         gs=13.95*ms/mn;
         gv=8.498*mv/mn;
         g1=0.0199*gs*gs*gs*mn;
         g2=-0.00296*gs*gs*gs*gs;
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.114) then   ! K=380 MeV m*/m=0.83 rho0=0.168
        gs=  6.448930391435675
        g1= -38.28483359158677
        g2= 342.08655094514916
        gv= 6.8548080859628415
         ms=2.79
         mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.


      else if (mstc(106).eq.115) then   ! K=210 MeV m*/m=0.83 rho0=0.168
        pard(101)= 7.907608564684195  ! Scalar part g_sigma
        pard(102)= 44.48219033563673        ! Scalar part g2 (1/fm)
        pard(103)= 22.691331489528103             ! Scalar part g3
        pard(111)= 6.8548080859628415               ! Vector part g_omega
        pard(112)= 2.79                ! sigma mass (1/fm)
        pard(113)= 3.97                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.116) then   ! K=380 MeV m*/m=0.63 rho0=0.168
        gs=  9.749707183693499
        g1= 5.166066279162845
        g2= -4.294987644673203
        gv= 11.40087186371087
        ms=2.79
        mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

c...NS2
      else if (mstc(106).eq.117) then   ! K=380 MeV m*/m=0.8 rho0=0.168
        pard(101)=  7.2110738469868885  ! Scalar part g_sigma
        pard(102)=  -17.888646553233265         ! Scalar part g2 (1/fm)
        pard(103)=  197.63846964767876          ! Scalar part g3
        pard(111)= 7.721187454057742              ! Vector part g_omega
        pard(112)= 2.79                ! sigma mass (1/fm)
        pard(113)= 3.97                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.


      else if (mstc(106).eq.118) then   ! K=210 MeV m*/m=0.85 rho0=0.168
        gs=  7.556630299429732
        g1= 41.130723737259125
        g2= 152.831761048523
        gv= 6.208251499369285
        ms=2.79
        mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.119) then   ! K=210 MeV m*/m=0.84 rho0=0.168
        gs=  7.744256707930581
        g1= 43.97572528567212
        g2= 73.94214333170736
        gv= 6.539723688876368
        ms=2.79
        mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.120) then   ! K=380 MeV m*/m=0.84 rho0=0.168
        gs=  6.132541126432369
        g1= -48.38145498219124
        g2= 406.53070136804973
        gv= 6.539723688876368
        ms=2.79
        mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.


      else if (mstc(106).eq.121) then   ! K=300 MeV m*/m=0.8 rho0=0.168
         ms=2.79
         mv=3.97
         gs=  7.716205383388513
         g1= 5.566120385681115
         g2= 110.40216192621635
         gv= 7.721187454057742

        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.122) then   ! K=240 MeV m*/m=0.8 rho0=0.168
         gs=  8.113799757099452
         g1= 27.53953936936618
         g2= 15.107648844227855
         gv= 7.721187454057742
         ms=2.79
         mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.


      else if (mstc(106).eq.123) then   ! K=300 MeV m*/m=0.722 rho0=0.168
        gs=  8.866001991402603
        g1= 11.808041964217319
        g2= 0.7801044590843295
        gv= 9.601382237013418

       ! K=240MeV does not work
       !gs=  9.096120377501608
       !g1= 21.203911530591597
       !g2= -37.92230322131382
       !gv= 9.601382237013418


         ms=2.79
         mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      ! NS3
      else if (mstc(106).eq.124) then   ! K=380 MeV m*/m=0.722 rho0=0.168
        gs=  8.562176910455312
        g1= 0.4429217593591168
        g2= 44.70405920222624
        gv= 9.601382237013418
         ms=2.79
         mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.125) then   ! K=300 MeV m*/m=0.76 rho0=0.168
          gs=  8.353120244786156
          g1= 10.960505005588097
          g2= 31.645061114519116
          gv= 8.739031343572359
        ms=2.79
        mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      ! NS1
      else if (mstc(106).eq.126) then   ! K=230 MeV m*/m=0.8 rho0=0.168
           ! K=230
           gs=  8.18183693589196
           g1= 31.622719664403174
           g2= -3.797739825658857
           gv= 7.721187454057742
         ms=2.79
         mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.127) then

           ! K=240 MeV m*/m=0.81 rho0=0.168
             gs=  7.964158855956515
             g1= 27.44775180881574
             g2= 38.11312917147808
             gv= 7.443987481918361


           ! K=210 MeV m*/m=0.80 rho0=0.168 does not work
         !gs=  8.319547680493411
         !g1= 40.18346413042249
         !g2= -44.529808144659235
         !gv= 7.721187454057742

         ms=2.79
         mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.
      else if (mstc(106).eq.128) then   ! K=380 MeV m*/m=0.7 rho0=0.168

         gs=  8.868721394180294
         g1= 2.193148626223312
         g2= 27.282385203907644
         gv= 10.064436308666291
         ms=2.79
         mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.129) then   ! K=300 MeV m*/m=0.75 rho0=0.168
        gs=  8.494244414361036
        g1= 11.39564805131095
        g2= 20.997831670872042
        gv= 8.974562635302759
         ms=2.79
         mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.130) then   ! K=340 MeV m*/m=0.678 rho0=0.168
       gs=  9.280947350664993
       g1= 7.370142764197598
       g2= -0.15036613894853254
       gv= 10.505168892531362
         ms=2.79
         mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.


      else if (mstc(106).eq.131) then   ! K=330 MeV m*/m=0.678 rho0=0.168
        gs=  9.311632864199382
        g1= 8.384449418073308
        g2= -3.968590286215991
        gv= 10.505168892531362

         ms=2.79
         mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.132) then   ! K=380 MeV m*/m=0.678 rho0=0.168
        gs=  9.158264292694614
        g1= 3.412164867373429
        g2= 14.494782498563072
        gv= 10.505168892531362
         ms=2.79
         mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.133) then   ! K=380 MeV m*/m=0.63 rho0=0.168
        gs=  9.749707183693499
        g1= 5.166066279162845
        g2= -4.294987644673203
        gv= 11.40087186371087

         ms=2.79
         mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.134) then   ! K=210 MeV m*/m=0.9 rho0=0.168
       gs=  5.529692800948896
       g1= -81.57448852407006
       g2= 1316.19097535084
       gv= 4.164372575012643
         ms=2.79
         mv=3.97
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)= ms               ! sigma mass (1/fm)
        pard(113)= mv                ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

c..2020/1/12 NS2'change Mnucl=0.939 and m_sigma, m_omega from 126
      else if (mstc(106).eq.150) then !NS1' K=210 MeV m*/m=0.8 rho0=0.168
        gs=  8.176015492139769
        g1= 31.48609618904411
        g2= -3.9926278758544633
        gv= 7.724939376757135
        pard(101)=  gs         ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv         ! Vector part g_omega
        pard(112)= 0.55/hc     ! sigma mass (1/fm)
        pard(113)= 0.783/hc    ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

c..2020/1/12 NS1'change Mnucl=0.939 and m_sigma, m_omega from 117
      else if (mstc(106).eq.151) then !NS2' K=380 MeV m*/m=0.8 rho0=0.168
        gs=  7.208208087191382
        g1= -17.740099004710263
        g2= 196.10588999399084
        gv= 7.724939376757135
        pard(101)=  gs         ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv         ! Vector part g_omega
        pard(112)= 0.55/hc     ! sigma mass (1/fm)
        pard(113)= 0.783/hc    ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

c..2020/1/12 NS3'change Mnucl=0.939 and m_sigma, m_omega from 124  NS3 in the paper 2020
      else if (mstc(106).eq.152) then !NS2' K=380 MeV m*/m=0.7 rho0=0.168
        gs=  8.863811172925521
        g1= 2.190629372074654
        g2= 27.070648189159165
        gv= 10.067747992111274
        pard(101)=  gs         ! Scalar part g_sigma
        pard(102)=  g1         ! Scalar part g2 (1/fm)
        pard(103)=  g2         ! Scalar part g3
        pard(111)=  gv         ! Vector part g_omega
        pard(112)= 0.55/hc     ! sigma mass (1/fm)
        pard(113)= 0.783/hc    ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.153) then !New NS1 K=380 MeV m*/m=0.83
        gs=  6.447696268747421
        g2= -38.001619264755874
        g3= 339.59550852599887
        gv= 6.858858568503943
        pard(101)=  gs         ! Scalar part g_sigma
        pard(102)=  g2         ! Scalar part g2 (1/fm)
        pard(103)=  g3         ! Scalar part g3
        pard(111)=  gv         ! Vector part g_omega
        pard(112)= 0.55/hc     ! sigma mass (1/fm)
        pard(113)= 0.783/hc    ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.154) then !New NS2 K=210 MeV m*/m=0.83
        gs=  7.901789322846869
        g2= 44.31272165001529
        g3= 21.98880833814041
        gv= 6.858858568503943
        pard(101)=  gs         ! Scalar part g_sigma
        pard(102)=  g2         ! Scalar part g2 (1/fm)
        pard(103)=  g3         ! Scalar part g3
        pard(111)=  gv         ! Vector part g_omega
        pard(112)= 0.55/hc     ! sigma mass (1/fm)
        pard(113)= 0.783/hc    ! omega mass (1/fm)
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.200) then ! K=380 MeV m*/m=0.65 rho0=0.168 SE
        gs= 8.446240226735835
        g2= 3.448041549527271
        g3= 12.667394877335811
        gv= 4.325553092079072
        gsp =  4.669818923094996
        gvp =  10.400937587189995
        mu1 = 0.6051656607160294 ! GeV
        mu2 = 1.5582862936704895 ! GeV
        ms = 0.55
        mv = 0.783
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g2         ! Scalar part g2 (1/fm)
        pard(103)=  g3         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)=  ms/hc            ! sigma mass (1/fm)
        pard(113)=  mv/hc               ! omega mass (1/fm)

        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv  ! vex2

        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.201) then ! K=380 MeV m*/m=0.65 Uopt->0 SE
        gs=  8.69894808411858
        g2= 3.3123747312338723
        g3= 10.151504907242417
        gv= 3.0515019650861466
        gsp =  4.147458199005445
        gvp =  10.805320816731728
        mu1 = 0.5377037526882253
        mu2 = 1.784448816788669
        ms = 0.55
        mv = 0.783
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g2         ! Scalar part g2 (1/fm)
        pard(103)=  g3         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)=  ms/hc            ! sigma mass (1/fm)
        pard(113)=  mv/hc               ! omega mass (1/fm)

        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2

        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.

      else if (mstc(106).eq.202) then ! K=380 MeV m*/m=0.65  omega = 0.0 SE
        gs=  4.239088246904197
        g2= -83.38119751325448
        g3= 1039.0072461869995
        gv= 0.0
        gsp =  8.299983097235295
        gvp =  11.415550399348907
        mu1 = 1.2
        mu2 = 1.3

        ms = 0.55
        mv = 0.783
        pard(101)=  gs ! Scalar part g_sigma
        pard(102)=  g2         ! Scalar part g2 (1/fm)
        pard(103)=  g3         ! Scalar part g3
        pard(111)=  gv             ! Vector part g_omega
        pard(112)=  ms/hc            ! sigma mass (1/fm)
        pard(113)=  mv/hc               ! omega mass (1/fm)

        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2

        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.

      ! MD1
      else if (mstc(106).eq.203) then ! K=380 MeV m*/m=0.65
      !} else if(eosType==14) {  // K=380 MeV m*/m=0.65
        gs=  9.02956219322301
        g2= 4.218109600022776
        g3= 6.667054503467467
        gv= 6.7402687603431675
        gsp =  3.1855887031669394
        gvp =  8.89548883191133
        mu1 = 0.6411295603159707
        mu2 = 1.841348887047491

        ms = 0.55
        mv = 0.783
        pard(101)=  gs       ! Scalar part g_sigma
        pard(102)=  g2       ! Scalar part g2 (1/fm)
        pard(103)=  g3       ! Scalar part g3
        pard(111)=  gv       ! Vector part g_omega
        pard(112)=  ms/hc    ! sigma mass (1/fm)
        pard(113)=  mv/hc    ! omega mass (1/fm)

        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2

        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.

c....small  MD2
c  } else if(eosType==17) {  // K=380 MeV m*/m=0.65 U(1.7)=0.06  U(6)=0.04
      else if (mstc(106).eq.204) then
c      gs=  9.283794468334284;
c      g2= 3.9415832307299548;
c      g3= 5.341587689277067;
c      gv= 4.892556211086802;
c      gsp =  2.306806519850028;
c      gvp =  9.993473536946235;
c      mu1 = 0.4268055932276202;
c      mu2 = 2.4686246424999148;

c   // K=380 MeV m*/m=0.65 U(1.7)=0.056  U(6)=0.025
      gs=  9.232564820463068
      g2= 4.012495169512666
      g3= 5.520251163159458
      gv= 3.8884399009614437
      gsp =  2.5019753284300768
      gvp =  10.431917516759935
      mu1 = 0.4896732480502826
      mu2 = 2.4886149663285533


        ms = 0.55
        mv = 0.783
        pard(101)=  gs       ! Scalar part g_sigma
        pard(102)=  g2       ! Scalar part g2 (1/fm)
        pard(103)=  g3       ! Scalar part g3
        pard(111)=  gv       ! Vector part g_omega
        pard(112)=  ms/hc    ! sigma mass (1/fm)
        pard(113)=  mv/hc    ! omega mass (1/fm)
        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2
        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.


c... MD3 smallest
c  } else if(eosType==18) { // Hama2 K=380 MeV m*/m=0.65 U(1.6)=0.052  U(3)=0.04
      else if (mstc(106).eq.205) then
        gs=  9.152032162363772
        g2= 4.1091984804037205
        g3= 5.886441070518481
        gv= 3.568789983273673
        gsp =  2.7878500218565083
        gvp =  10.555177917527832
        mu1 = 0.5668789032308825
        mu2 = 2.3990906874902547

        ms = 0.55
        mv = 0.783
        pard(101)=  gs       ! Scalar part g_sigma
        pard(102)=  g2       ! Scalar part g2 (1/fm)
        pard(103)=  g3       ! Scalar part g3
        pard(111)=  gv       ! Vector part g_omega
        pard(112)=  ms/hc    ! sigma mass (1/fm)
        pard(113)=  mv/hc    ! omega mass (1/fm)
        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2
        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.

c....  MD4 largest Upot(infty) = 200 MeV
c  } else if(eosType==22) { // K=380 MeV m*/m=0.65 U(1.7)=0.09 U(9)=0.18
      else if (mstc(106).eq.206) then
        gs=  8.943540882257514
        g2= 3.917672746640286
        g3= 8.728327369798889
        gv= 9.87112877312085
        gsp =  3.636579180635797
        gvp =  5.481831782950353
        mu1 = 0.5039127940027004
        mu2 = 0.6998101782971953
        ms = 0.55
        mv = 0.783
        pard(101)=  gs       ! Scalar part g_sigma
        pard(102)=  g2       ! Scalar part g2 (1/fm)
        pard(103)=  g3       ! Scalar part g3
        pard(111)=  gv       ! Vector part g_omega
        pard(112)=  ms/hc    ! sigma mass (1/fm)
        pard(113)=  mv/hc    ! omega mass (1/fm)
        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2
        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.

c....  MD4'  Upot(infty) = 175 MeV
c  } else if(eosType==21) { // K=380 MeV m*/m=0.65 U(1.7)=0.08 U(9)=0.15
      else if (mstc(106).eq.207) then
        gs=  9.304844831497922
        g2= 3.835117791722945
        g3= 5.5579843719039905
        gv= 9.199390296353101
        gsp =  2.277148182541522
        gvp =  6.277561735662023

        mu1 = 0.3652112126911749
        mu2 = 1.4301429116848454
        ms = 0.55
        mv = 0.783
        pard(101)=  gs       ! Scalar part g_sigma
        pard(102)=  g2       ! Scalar part g2 (1/fm)
        pard(103)=  g3       ! Scalar part g3
        pard(111)=  gv       ! Vector part g_omega
        pard(112)=  ms/hc    ! sigma mass (1/fm)
        pard(113)=  mv/hc    ! omega mass (1/fm)
        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2
        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.

c....  MD4'  Upot(infty) = 150 MeV
c  } else if(eosType==20) { // K=380 MeV m*/m=0.65 U(1.7)=0.072 U(9)=0.13
      else if (mstc(106).eq.208) then
        gs=  9.15499996067294
        g2= 4.0336963331563656
        g3= 6.143870343679022
        gv= 8.620286411780615
        gsp =  2.81531192821113
        gvp =  7.084464641127519
        mu1 = 0.5269385823575612
        mu2 = 1.4655608399069255
        ms = 0.55
        mv = 0.783
        pard(101)=  gs       ! Scalar part g_sigma
        pard(102)=  g2       ! Scalar part g2 (1/fm)
        pard(103)=  g3       ! Scalar part g3
        pard(111)=  gv       ! Vector part g_omega
        pard(112)=  ms/hc    ! sigma mass (1/fm)
        pard(113)=  mv/hc    ! omega mass (1/fm)
        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2
        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.

c....  MD5
c  } else if(eosType==24) { // K=210 MeV m*/m=0.8 U(1.7)=0.06 U(5)=0.075
      else if (mstc(106).eq.209) then
        gs=  6.594093254363855
        g2= 89.86702225927918
        g3= 108.51032370440886
        gv= 6.3580918751932325
        gsp =  5.277392175886514
        gvp =  4.4396533278183234
        mu1 = 0.6983616678464525
        mu2 = 2.207541181176542
        ms = 0.55
        mv = 0.783
        pard(101)=  gs       ! Scalar part g_sigma
        pard(102)=  g2       ! Scalar part g2 (1/fm)
        pard(103)=  g3       ! Scalar part g3
        pard(111)=  gv       ! Vector part g_omega
        pard(112)=  ms/hc    ! sigma mass (1/fm)
        pard(113)=  mv/hc    ! omega mass (1/fm)
        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2
        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.

c....  MD6   -> MD4 in the paper 2020
c } else if(eosType==25) { // K=210 MeV m*/m=0.83 U(1.7)=0.06 U(5)=0.07
      else if (mstc(106).eq.210) then

       gs=  4.0586271572892025
       g2= -160.31303134630116
       g3= 2684.378830396555
       gv= 5.63247647754081
       gsp =  5.543944245008469
       gvp =  3.9266304417034466
       mu1 = 0.7035464013689996
       mu2 = 4.252421078959087
        ms = 0.55
        mv = 0.783
        pard(101)=  gs       ! Scalar part g_sigma
        pard(102)=  g2       ! Scalar part g2 (1/fm)
        pard(103)=  g3       ! Scalar part g3
        pard(111)=  gv       ! Vector part g_omega
        pard(112)=  ms/hc    ! sigma mass (1/fm)
        pard(113)=  mv/hc    ! omega mass (1/fm)
        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2
        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.




      ! MD2 old
c} else if(eosType==15) {  // K=210 MeV m*/m=0.8 U(1.7)=0.058 U(3)=0.065
      else if (mstc(106).eq.211) then
        gs=  6.551389416446063
        g2= 90.3843215960433
        g3= 150.1299910657155
        gv= 5.757865467744954
        gsp =  5.324936778292298
        gvp =  5.198211638397418
        mu1 = 0.7053607242477306
        mu2 = 2.4912130056684383

        ms = 0.55
        mv = 0.783
        pard(101)=  gs       ! Scalar part g_sigma
        pard(102)=  g2       ! Scalar part g2 (1/fm)
        pard(103)=  g3       ! Scalar part g3
        pard(111)=  gv       ! Vector part g_omega
        pard(112)=  ms/hc    ! sigma mass (1/fm)
        pard(113)=  mv/hc    ! omega mass (1/fm)
        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2
        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.

c....  MD6 omega=0.0 extremely small !   MD3 in the paper 2020
c } else if(eosType==26) { // K=380 MeV m*/m=0.65 U(1.6)=0.05 U(10)=0.0
      else if (mstc(106).eq.212) then
        gs=  5.4390747082116775
        g2= -15.587197606469752
        g3= 391.86527371750225
        gv= 0.0
        gsp =  7.711348424066859
        gvp =  11.22012602741299
        mu1 = 1.7024475104128964
        mu2 = 1.898341169747188

        ms = 0.55
        mv = 0.783
        pard(101)=  gs       ! Scalar part g_sigma
        pard(102)=  g2       ! Scalar part g2 (1/fm)
        pard(103)=  g3       ! Scalar part g3
        pard(111)=  gv       ! Vector part g_omega
        pard(112)=  ms/hc    ! sigma mass (1/fm)
        pard(113)=  mv/hc    ! omega mass (1/fm)
        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2
        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.

c...eos27
      else if (mstc(106).eq.213) then
       gs=  9.048354922793935
       g2= 4.20958829512232
       g3= 6.5093840573668444
       gv= 5.900891294185006
       gsp =  3.1245160106443817
       gvp =  9.467708263559453
       mu1 = 0.6357372647745847
       mu2 = 2.005796795262418

        ms = 0.55
        mv = 0.783
        pard(101)=  gs       ! Scalar part g_sigma
        pard(102)=  g2       ! Scalar part g2 (1/fm)
        pard(103)=  g3       ! Scalar part g3
        pard(111)=  gv       ! Vector part g_omega
        pard(112)=  ms/hc    ! sigma mass (1/fm)
        pard(113)=  mv/hc    ! omega mass (1/fm)
        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2
        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.

c....K=380MeV m*/m=0.78 U( 1.7 )= 0.058,U( 0.65 )= 0.0,U( 3.0 )= 0.065 
      else if (mstc(106).eq.214) then

        gsp =  5.125888811529892
        gvp =  5.79518230378181
        mu1 = 0.7013309013578113
        mu2 = 2.267161800762232
        gs=  5.33703361590903
        g2= -54.695147262252604
        g3= 483.26127733316866
        gv= 5.935397496403753

        ms = 0.55
        mv = 0.783
        pard(101)=  gs       ! Scalar part g_sigma
        pard(102)=  g2       ! Scalar part g2 (1/fm)
        pard(103)=  g3       ! Scalar part g3
        pard(111)=  gv       ! Vector part g_omega
        pard(112)=  ms/hc    ! sigma mass (1/fm)
        pard(113)=  mv/hc    ! omega mass (1/fm)
        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2
        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.

c K=210MeV m*/m=0.65 U( Pinf= 1.7 )= 0.056 c U( 0.65 )= 0.0 U( 6.0 )= 0.025
      else if (mstc(106).eq.215) then

        gsp =  2.5019753284300768
        gvp =  10.431917516759935
        mu1 = 0.4896732480502826
        mu2 = 2.4886149663285533   ! in GeV
        gs=  9.715848996596957
        g2= 19.640169429599684
        g3= 0.0
        g4= 5.739407608095525      ! in fm^2
        gv= 3.8884399009614437

       gs=  9.534423675188515;
       g2= 21.44720328182629;
       g3= 0.0;
       g4= 6.379088354903728;
       gv= 6.7402687603431675;
       gsp =  3.1855887031669394;
       gvp =  8.89548883191133;
       mu1 = 0.6411295603159707;
       mu2 = 1.841348887047491;


        ms = 0.55
        mv = 0.783
        pard(101)=  gs       ! Scalar part g_sigma
        pard(102)=  g2       ! Scalar part g2 (1/fm)
        pard(103)=  g3       ! Scalar part g3
        pard(116)=  g4       ! Scalar part g4
        pard(111)=  gv       ! Vector part g_omega
        pard(112)=  ms/hc    ! sigma mass (1/fm)
        pard(113)=  mv/hc    ! omega mass (1/fm)
        pard(105)=mu1  ! pmu1 dummy
        pard(106)=mu2  ! pmu2 dummy
        pard(107)= -(gsp/ms)**2*conv   ! vex1
        pard(108)=  (gvp/mv)**2*conv   ! vex2
        mstd(101)=1  ! momentum dependent potential is used.
        mstd(102)=1  ! nonlinear sigma-field is used.

      else
        write(6,*)'RQMD.RMF wrong number mstc(106)=',mstc(106)
        stop
      endif

      if(mstc(109).eq.0) then
        print *,'inconsistent mstc(109) should be 2',mstc(109)
        stop
      endif

c...RQMD/NV non-linear vector potential version
      if(mstd(102).eq.0) then
        t1=0.5*pard(101)
        t2=0.5*pard(111)
        t3=pard(112)/(pard(113)+1.0)
        t3f=pard(113)*t3
        t4=pard(114)/(pard(115)+1.0)
        t4f=pard(115)*t4
        if(mstc(6).ge.100) write(mstc(37),*)'RQMD/NV mode'

c       t1=pard(101)
c       t2=pard(111)
c       t3=pard(112)
c       t3f=pard(113)*t3

c...RQMD/NS  non-linear sigma
      else if(mstd(102).eq.1) then
        t1=-0.5*pard(101)
        t2=0.5*(pard(111)/pard(113))**2*paru(3) ! GeVfm^3

        t3=0.0d0
        t3f=0.0d0
        t4=0.0d0
        t4f=0.0d0

c       g24=2.0/3.0*g2*g4/hc
c       g34=5.0/4.0*g3*g4/hc/hc
c       gs4=g4*gs/hc/hc
c       gg4=0.5*g4/hc/hc


        if(mstc(6).ge.100) write(mstc(37),*)'RQMD/NS mode'
      endif

      pmu1=pard(105)**2
      pmu2=pard(106)**2
      vex1=pard(107)/2
      vex2=pard(108)/2

      if(mstc(6).ge.100) then
      write(mstc(38),*)'t1=',t1,'t2=',t2,'t3f=',t3f
      write(mstc(38),*)'t3=',t3
      write(mstc(38),*)'pard112=',pard(111),pard(112),pard(113)
      write(mstc(38),*)'vex1=',vex1,'vex2=',vex2,
     & 'pmu1=',pmu1,'pmu2=',pmu2
      endif

c...Width parameter.
c     pard(109)= 1.0/(4.d0*pard(104))
c     pard(110)=1.0d0/(4.d0*paru(1)*pard(104))**1.5d0 ! 1/[(4*pi*L)^3/2]

      pard(109)= 1.0/(2.d0*pard(104))
      pard(110)=1.0d0/(2.d0*paru(1)*pard(104))**1.5d0 ! 1/[(4*pi*L)^3/2]

      if(mstc(6).ge.2.and.mstc(119).eq.0.and.mstc(115).eq.2) then
        mstd(91)=1
c       mstd(91)=2
      endif

c     if(mstd(102)==1) then  ! nonlinear sigma-field is used.
c      call make_sigma_table
c     endif

      end
