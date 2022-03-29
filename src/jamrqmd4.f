c***********************************************************************
c                                                                      *
c        PART  : RQMD  Vector potential implementation of Skyrme force *
c                                                                      *
c   List of Subprograms in rough order of relevance with main purpose  *
c      (s = subroutine, f = function, b = block data, e = entry)       *
c  f fenga      to pre-factor from the baryon current                  *
c  s rqmddev1c  to calculate derivatives                               *
c  s jamvpotdot4 to compute time-derivatives of the baryon current.    *
c                                                                      *
************************************************************************

      subroutine jamrqmd4

c...Purpose: to calculate force in RQMD/S
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,NrqmdCount,icheckMF
      real*8 diff,emi,emj,emf1,emf2
      real*8 dt,diffv,vfac1,vfac,vfacs
      real*8 Ai(0:3),Aj(0:3),pkin1(5),pkin2(5)
      real*8 vdot
c     real*8 vtot(0:3)
      common/jamvdot/vdot(3,mxv)

c...Check which particle feels potentials.
      NrqmdCount=icheckMF(0,pard(1))

c....set p0 to be the free kinetic energy for the calculation of
c....potential.
c     if(mstc(112).ge.1.and.mstc(109).ge.2) call setenergy(0)

      forcer(:,:) = 0.0d0
      force(:,:) = 0.0d0
      vdot(:,:)=0.0d0

      call jamrqmm4
      call jamepart4(diff,diffv)
c     if(mstc(112).eq.0) call setenergy(1)

c...Compute scalar and vector potential self-consistently.
c     if(mstc(109).ge.3) then
c       call jamsvpot4(diff,diffv)
c     endif

c...Loop over all particles.
      rhog(:)=0.0
      do i=1,nv
        if(MF_on(i)==0) cycle
        if(rhoi(i).lt.1d-8) cycle
        rhog(i)=rhoi(i)**(pard(103)-1.0d0)
      end do
 
c...Loop over all particles i.
      do i=1,nv  ! Sum for ith particle (sigma_i) 

         if(MF_on(i).eq.0) cycle
         emi=sqrt(p(4,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2)
         pkin1(:)=p(:,i)

c...Kinetic momentum in case p(,i) is canonical momentum.
         if(mstc(119).eq.0.and.mstc(115).eq.2) then
           pkin1(1)=p(1,i)-potv(1,i)
           pkin1(2)=p(2,i)-potv(2,i)
           pkin1(3)=p(3,i)-potv(3,i)
         endif
         emf1=p(5,i)
         pkin1(4)=sqrt(emf1**2+pkin1(1)**2+pkin1(2)**2+pkin1(3)**2)

c        call fengd(i,pkin1,Ai)
         vfac1=vfacs(i)

c...Loop over particles j.
         do j=i+1,nv     ! Sum for j (sigma_j(not=i))
           if(MF_on(j)==0) cycle
           emj=sqrt(p(4,j)**2-p(1,j)**2-p(2,j)**2-p(3,j)**2)
           pkin2(:)=p(:,j)
           if(mstc(119).eq.0.and.mstc(115).eq.2) then
             pkin2(1)=p(1,j)-potv(1,j)
             pkin2(2)=p(2,j)-potv(2,j)
             pkin2(3)=p(3,j)-potv(3,j)
           endif
           emf2=p(5,j)
           pkin2(4)=sqrt(emf2**2+pkin2(1)**2+pkin2(2)**2+pkin2(3)**2)
           call fengd(i,pkin2,Ai)
           call fengd(j,pkin1,Aj)

           vfac=vfac1*vfacs(j)

          if(mstc(113).eq.0) then
            call rqmddev0e(i,j,Ai,Aj,pkin1,pkin2,vfac)
          else if(mstc(113).le.1) then
            call rqmddev1e(i,j,Ai,Aj,emi,emj,pkin1,pkin2,vfac)
          else
            call rqmddev2e(i,j,Ai,Aj,emi,emj,pkin1,pkin2,vfac)
          endif

        end do
      end do


c...Compute the time derivative of the vector potential to compute the kinetic
c...momentum in case vector potential is fully included.
      if(mstc(119).eq.1.and.mstc(115).ge.2) call jamvpotdot4

c....set p0 using effective mass.
c     if(mstc(112).ge.1) call setenergy(1)

c      vtot=0.0

        dt=parc(2)
        do i=1,nv
          if(MF_on(i)==0) cycle
            if(mstc(119).eq.1.or.mstc(119).eq.2) then
              p(1,i)=p(1,i)+dt*(force(1,i) - vdot(1,i))
              p(2,i)=p(2,i)+dt*(force(2,i) - vdot(2,i))
              p(3,i)=p(3,i)+dt*(force(3,i) - vdot(3,i))
            else
              p(1,i)=p(1,i)+dt*force(1,i)
              p(2,i)=p(2,i)+dt*force(2,i)
              p(3,i)=p(3,i)+dt*force(3,i)
            endif
            p(4,i)=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
c displacement only with interaction
            r(1,i)=r(1,i)+dt*forcer(1,i)
            r(2,i)=r(2,i)+dt*forcer(2,i)
            r(3,i)=r(3,i)+dt*forcer(3,i)

c           vtot = vtot + potv(:,i)
        enddo

        call jamrqmm4
        call jamepart4(diff,diffv)

      if(mstc(119).eq.3) then
        do i=1,nv
          if(MF_on(i)==0) cycle
          p(1,i)=p(1,i)-dt*vdot(1,i)
          p(2,i)=p(2,i)-dt*vdot(2,i)
          p(3,i)=p(3,i)-dt*vdot(3,i)
        end do
      endif


c....set p0 using effective mass.
c     if(mstc(109).ge.2) call setenergy(1)

      end

c***********************************************************************

      subroutine jamepart4(difm,difv)  ! energy

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
      integer i,nn
      real*8 difm,difv,potv0(0:3),vsky
      real*8 vdott(3),vdot
      common/jamvdot/vdot(3,mxv)

      difm=0.d0   ! Difference of effective mass
      difv=0.d0   ! Difference of vector potential
      vdott=0.0
      nn=0

      do 100 i=1,nv

       if(MF_on(i)==0) then
         pots(i)=0.0d0
         potv(:,i)=0.0d0
         goto 100
       endif

       potv0(:)=potv(:,i)
       potv(1,i)=0.0
       potv(2,i)=0.0
       potv(3,i)=0.0

       pots(i)=0.0
       vsky=k(9,i)/3*(t1*rhoi(i) + t3*rhoi(i)**pard(103))
       potv(0,i) = vsky + rhoc(i) + rhos(i) + rhoy(i) + vmom(0,i)
       if(mstc(115).eq.2.and.rhoi(i).gt.1d-8) then
         potv(:,i)=vsky*rhoj(:,i)/rhoi(i) + vmom(:,i)
       endif

       difv=difv+sqrt((potv0(0)-potv(0,i))**2+(potv0(1)-potv(1,i))**2
     &    +(potv0(2)-potv(2,i))**2+(potv0(3)-potv(3,i))**2)

       if(mstc(119).ge.2) then
       vdot(1,i) = (potv(1,i) - potv0(1))/parc(2)
       vdot(2,i) = (potv(2,i) - potv0(2))/parc(2)
       vdot(3,i) = (potv(3,i) - potv0(3))/parc(2)
       vdott(:) = vdott(:) + vdot(:,i)
       nn = nn + 1
       endif

100   continue

      if(mstc(119).ge.2.and.nn.ge.1) then
      do i=1,nv
       if(MF_on(i).eq.0) cycle
       vdot(:,i) = vdot(:,i) - vdott(:)/nn
      end do
      endif

      mste(44)=1
      end

**********************************************************************

      subroutine jamrqmm4

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
      real*8 qfacr,pk1(5),pk2(5)
      real*8 pmom2ij,pmom2ji
      real*8 gfac,vfac,vfac1,vfacs,ei,ej,gi,gj
      real*8 den,gam,r2i,r2j,p2i,p2j,pij

      fac  =(4.d0*paru(1)*pard(104))**1.5d0 !    [(4*pi*L)^3/2]
      wg   = 1.0/4.d0/pard(104)

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
        emisq=p(4,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2
        vfac1=vfacs(i)
        emi=sqrt(emisq)
        ei = p(4,i)
        gi = 1.0

        do 110 j = i+1 , nv      ! Sum for j (sigma_j(>i))

          if(MF_on(j).eq.0) goto 110
          bj=k(9,j)/3
          emjsq = p(4,j)**2 - p(1,j)**2 - p(2,j)**2 - p(3,j)**2
          ej = p(4,j)
          gj = 1.0
          emj=sqrt(emjsq)

          rx = r(1,i) - r(1,j)
          ry = r(2,i) - r(2,j)
          rz = r(3,i) - r(3,j)
          ps=(p(1,i)-p(1,j))**2+(p(2,i)-p(2,j))**2+(p(3,i)-p(3,j))**2
          rs= rx**2 + ry**2 + rz**2
          gfac=qfacr(i)*qfacr(j)
          vfac=vfac1*vfacs(j)

c...Non-relativistic.
          if(mstc(113).eq.0) then

          rhom(i,j) = exp(-rs*wg) * gfac/fac
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
          den  = exp(-r2*wg) * gfac /fac
          gam=1.0
          if(mstc(114).ne.3) gam = pe/sqrt(s)
          rhom(i,j) = gam * den !* emj*bj
          rhom(j,i) = gam * den !* emi*bi

          p2 = ps - (p(4,i)-p(4,j))**2 + (emisq- emjsq)**2/s
          pmom2ij = vex1/(1.d0+p2/pmu1)+vex2/(1.d0+p2/pmu2)
          pmom2ji = pmom2ij
          p2i=p2
          p2j=p2

c...distance is measured from the rest-frame of particle j.
          else

          r2i=rs + (rx*p(1,i)+ry*p(2,i)+rz*p(3,i))**2/emisq
          r2j=rs + (rx*p(1,j)+ry*p(2,j)+rz*p(3,j))**2/emjsq
          rhom(i,j)= exp(-r2j*wg) * gfac /fac
          rhom(j,i)= exp(-r2i*wg) * gfac /fac

c...Distance in the two-body cm is used for the argument of MD pot.
          if(mstc(120).eq.1) then
            s=(p(4,i)+p(4,j))**2
     &       -(p(1,i)+p(1,j))**2-(p(2,i)+p(2,j))**2-(p(3,i)+p(3,j))**2
            p2i = ps - ((p(4,i)-p(4,j))**2 - (emisq- emjsq)**2/s)
            p2j=p2i
          else
            pij=p(4,i)*p(4,j)-p(1,i)*p(1,j)-p(2,i)*p(2,j)-p(3,i)*p(3,j)
            p2i=ps - (p(4,i)-p(4,j))**2 + (pij-emisq)**2/emisq
            p2j=ps - (p(4,i)-p(4,j))**2 + (pij-emjsq)**2/emjsq
          endif
          pmom2ij = vex1/(1.d0+p2j/pmu1)+vex2/(1.d0+p2j/pmu2)
          pmom2ji = vex1/(1.d0+p2i/pmu1)+vex2/(1.d0+p2i/pmu2)
          ei = emi
          ej = emj
          gi = p(4,i)/emi
          gj = p(4,j)/emj

          endif

          if(mstc(118).ge.1) then
            rhom(i,j)=rhom(i,j)/(1.0+p2j/parc(105)**2)
            rhom(j,i)=rhom(j,i)/(1.0+p2i/parc(105)**2)
          endif

          rhoi(i) = rhoi(i) + rhom(i,j)*bj*vfac
          rhoi(j) = rhoi(j) + rhom(j,i)*bi*vfac

          if (mstd(101).eq.0) cycle
c         vmom(i) = vmom(i) + pmom2ij*rhom(i,j)
c         vmom(j) = vmom(j) + pmom2ji*rhom(j,i)

          vmom(0,i) = vmom(0,i) + pmom2ij*rhom(i,j)*gj*vfac
          vmom(1,i) = vmom(1,i) + pmom2ij*rhom(i,j)*p(1,j)/ej*vfac
          vmom(2,i) = vmom(2,i) + pmom2ij*rhom(i,j)*p(2,j)/ej*vfac
          vmom(3,i) = vmom(3,i) + pmom2ij*rhom(i,j)*p(3,j)/ej*vfac

          vmom(0,j) = vmom(0,j) + pmom2ji*rhom(j,i)*gi*vfac
          vmom(1,j) = vmom(1,j) + pmom2ji*rhom(j,i)*p(1,i)/ei*vfac
          vmom(2,j) = vmom(2,j) + pmom2ji*rhom(j,i)*p(2,i)/ei*vfac
          vmom(3,j) = vmom(3,j) + pmom2ji*rhom(j,i)*p(3,i)/ei*vfac

 110     continue               ! Loop end of sigma_j(>i)
 100  continue                  ! Loop end of sigma_i
 
c...Compute baryon current J_i^\mu
      do 8 i=1,nv
         if(MF_on(i).eq.0) goto 8

         bi=k(9,i)/3
         pk1(:)=p(:,i)
         if(mstc(112).eq.1) call getkinmom2(i,pk1)
         emi=pk1(4)
         if(mstc(114).eq.3.or.mstc(113).eq.2) then
           emi=sqrt(pk1(4)**2-pk1(1)**2-pk1(2)**2-pk1(3)**2)
         endif
         vfac1=vfacs(i)
         do 9 j=i+1,nv
           if(MF_on(j).eq.0) goto 9
           bj=k(9,j)/3
           pk2(:)=p(:,j)
           if(mstc(112).eq.1) call getkinmom2(j,pk2)
           emj=pk2(4)
           if(mstc(114).eq.3.or.mstc(113).eq.2) then
             emj=sqrt(pk2(4)**2-pk2(1)**2-pk2(2)**2-pk2(3)**2)
           endif
           vfac=vfac1*vfacs(j)
           rhoj(0,i) = rhoj(0,i) + rhom(i,j)*pk2(4)/emj*vfac*bj
           rhoj(1,i) = rhoj(1,i) + rhom(i,j)*pk2(1)/emj*vfac*bj
           rhoj(2,i) = rhoj(2,i) + rhom(i,j)*pk2(2)/emj*vfac*bj
           rhoj(3,i) = rhoj(3,i) + rhom(i,j)*pk2(3)/emj*vfac*bj
           rhoj(0,j) = rhoj(0,j) + rhom(j,i)*pk1(4)/emi*vfac*bi
           rhoj(1,j) = rhoj(1,j) + rhom(j,i)*pk1(1)/emi*vfac*bi
           rhoj(2,j) = rhoj(2,j) + rhom(j,i)*pk1(2)/emi*vfac*bi
           rhoj(3,j) = rhoj(3,j) + rhom(j,i)*pk1(3)/emi*vfac*bi
 9       continue
 8    continue

c...Compute invariant baryon density.
      if(mstc(114).ne.0) then
        do i=1,nv
          if(MF_on(i).eq.0) cycle
          rhoi(i)=sqrt(max(0d0,rhoj(0,i)**2
     &        -rhoj(1,i)**2-rhoj(2,i)**2-rhoj(3,i)**2))
        enddo
      endif

      end
************************************************************************
c....Pre-factor from the baryon current.
      subroutine fengd(i,pkin,Ai)
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
        dv = (t1 + t3f*rhog(i))/rhoi(i)
        Ai(:)=dv*rhoj(:,i)
        return
      endif

      vj=rhoj(0,i)-bi(1)*rhoj(1,i)-bi(2)*rhoj(2,i)-bi(3)*rhoj(3,i)
      vv = t1 + t3*rhog(i)  ! V/rho_B
      dv = (pard(103)-1)*t3*rhoi(i)**(pard(103)-3) ! del(V/rho)/rho
      Ai(:)=vj*dv*rhoj(:,i) + vv*bi(:)

      end

**********************************************************************

      subroutine rqmddev0e(i,j,Ai,Aj,pkin1,pkin2,vfac)

c...Purpose: to compute derivatives of the non-relativistic distance
c...The Skyrme potential is treated as Lorentz vector.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,n
      real*8 pc1(5),pc2(5),bi(0:3),bj(0:3),Ai(0:3),Aj(0:3)
      real*8 vfac,fengi,fengj
      real*8 fsky3,fsky4,fsky
      real*8 fsgam2i,fsgam2j,fsgam3i,fsgam3j
      real*8 p2,fac1,fac2,fengij,fmome,fmomd,facmom
      real*8 pkin1(5),pkin2(5),vf1,vf2,fm1,fm2,fm3,fm4
      real*8 A1(0:3),A2(0:3)

      pc1(:)=p(:,i)
      pc2(:)=p(:,j)

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

c...Vector part.
      fsky3= fengi*rhom(i,j)*vfac*k(9,i)*k(9,j)/9
      fsky4= fengj*rhom(j,i)*vfac*k(9,i)*k(9,j)/9

c...Total.
      fsky= -pard(109)*(fsky3 + fsky4)

      fsgam2i=0
      fsgam2j=0
      fsgam3i=0
      fsgam3j=0
      A1=0.0
      A2=0.0
      if(mstc(111).eq.1) then

c....Derivatives of p^\mu/_i/p^0_i term in the vector part.
        
        call fengd(i,pkin1,A1)
        call fengd(j,pkin2,A2)
        fsgam2i=fsgam2i + mstc(116)*rhom(j,i)*
     &         ( A2(1)*bi(1) + A2(2)*bi(2) + A2(3)*bi(3) )/pc1(4)
        fsgam3i=-rhom(j,i)/pc1(4)

        fsgam2j=fsgam2j + mstc(116)*rhom(i,j)*
     &         ( A1(1)*bj(1) + A1(2)*bj(2) + A1(3)*bj(3) )/pc2(4)
        fsgam3j=-rhom(i,j)/pc2(4)
      endif

      do n=1,3
        forcer(n,i) = forcer(n,i) + bi(n)*fsgam2i + A2(n)*fsgam3i
        forcer(n,j) = forcer(n,j) + bj(n)*fsgam2j + A1(n)*fsgam3j
        force(n,i)  = force(n,i)  - 2*fsky*(r(n,i)-r(n,j)) ! dP_i/dt
        force(n,j)  = force(n,j)  - 2*fsky*(r(n,j)-r(n,i)) ! dP_j/dt
      end do

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
 
c....Momentum dependent potential part start.
      p2 = (p(1,i)-p(1,j))**2 + (p(2,i)-p(2,j))**2 + (p(3,i)-p(3,j))**2
      fac1 = 1.d0 + p2/pmu1
      fac2 = 1.d0 + p2/pmu2

c...p Derivative term.
      fengij=vf1*rhom(i,j) + vf2*rhom(j,i)
      fmome=-fengij*(vex1/pmu1/fac1**2 + vex2/pmu2/fac2**2)

c...r Derivative term.
      facmom=vex1/fac1 + vex2/fac2
      fmomd=-pard(109)*fengij*facmom

      do n=1,3
        force(n,i)  = force(n,i)  - 2*fmomd*(r(n,i)-r(n,j))
        force(n,j)  = force(n,j)  - 2*fmomd*(r(n,j)-r(n,i))
        forcer(n,i) = forcer(n,i) + 2*fmome*(p(n,i)-p(n,j))
        forcer(n,j) = forcer(n,j) + 2*fmome*(p(n,j)-p(n,i))
      enddo

c...Derivative of  p_i/p^0_i part.
      if(mstc(111).eq.1.and.mstc(114).ne.0.and.mstc(115).eq.2) then
        fm1=-rhom(j,i)*facmom/p(4,i)/pkin2(4)
        fm2=-rhom(i,j)*facmom/p(4,j)/pkin1(4)
c       vf2 = (pkin2(1)*bi(1)+pkin2(2)*bi(2)+pkin2(3)*bi(3))/pkin2(4)
c       vf1 = (pkin1(1)*bj(1)+pkin1(2)*bj(2)+pkin1(3)*bj(3))/pkin1(4)
        vf1 = 1.0 - vf1
        vf2 = 1.0 - vf2
        fm3=mstc(116)*facmom/p(4,i)*rhom(j,i)*vf1
        fm4=mstc(116)*facmom/p(4,j)*rhom(i,j)*vf2
        do n=1,3
          forcer(n,i) = forcer(n,i) + fm1*pkin2(n) + bi(n)*fm3
          forcer(n,j) = forcer(n,j) + fm2*pkin1(n) + bj(n)*fm4
        enddo
      endif

      end

**********************************************************************

      subroutine rqmddev1e(i,j,Ai,Aj,emi,emj,pkin1,pkin2,vfac)

c...Purpose: to compute derivatives of the squared four-vector distance
c...q_{Tij}^2 and p_{Tij}^2 in the two body c.m. of particle i and j.
c...The Skyrme potential is treated as Lorentz vector.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,n
      real*8 bi(0:3),bj(0:3),vi(0:3),vj(0:3)
      real*8 dr2ri(3),dr2pi(3),dp2pi(3)
      real*8 dr2rj(3),dr2pj(3),dp2pj(3)
      real*8 pc1(5),pc2(5),pcm(4),bet(3),rk(3),pk(4)
      real*8 s,rbij,pma,fengi,fengj,fengij,fsky
      real*8 fsky3,fsky4,vfac
      real*8 fsgam1,fsgam2i,fsgam2j,fmgam1,facsk
      real*8 fsgam3i,fsgam3j,fmgam2i,fmgam2j,fmgam3i,fmgam3j
      real*8 p2,fac1,fac2,facmom,fmomd,fmome
      real*8 bbi(3),bbj(3),wg,emi,emj,pe1,pe2
      real*8 Ai(0:3),Aj(0:3), A1(0:3),A2(0:3)
      real*8 pkin1(5),pkin2(5),fm1,fm2,fm3,fm4,vf1,vf2,facm1,facm2

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
        vi(n)=p(n,i)/pe1
        vj(n)=p(n,j)/pe2
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
c       fengi=0.0
c       fengj=0.0
c       do n=0,3
c         fengi = fengi + Ai(n)*vi(n)
c         fengj = fengj + Aj(n)*vj(n)
c       end do

        fengi=Ai(0)*vi(0)
        fengj=Aj(0)*vj(0)
        do n=1,3
          fengi = fengi - Ai(n)*vi(n)
          fengj = fengj - Aj(n)*vj(n)
        end do

c     fengi=fengi*k(9,i)/3
c     fengj=fengj*k(9,j)/3

c...Width parameter.
      wg = 1.0/(4.d0*pard(104))

c...Vector part.
      fsky3= fengi*rhom(i,j)*vfac*k(9,i)*k(9,j)/9
      fsky4= fengj*rhom(j,i)*vfac*k(9,i)*k(9,j)/9
      fsky= -wg*(fsky3 + fsky4)

c...Coulomb, Symmetry e, and Yukawa potentials are treated as vector.
      fsky=fsky + rhc(i,j)+rhc(j,i)  ! Add Coulomb potential
      fsky=fsky + rhs(i,j)+rhs(j,i)  ! Add symmetry energy
      fsky=fsky + rhy(i,j)+rhy(j,i)  ! Add Yukawa potential

      fsgam1=0.0  ! p_i(n) + p_j(n) term
      fsgam2i=0.0 ! bi(n) term
      fsgam2j=0.0 ! bi(n) term
      fsgam3i=0.0 ! j_i(n) term
      fsgam3j=0.0 ! j_i(n) term
      A1=0.0
      A2=0.0
c...for now option mstc(114)=2,3 is only for Skyrme + mom. dep. force.

c.....Derivatives of gamma_{ij}
      if(mstc(114).ne.3) then
         facsk = fsky3 + fsky4
         fsgam1= facsk/s
         fsgam2i=mstc(116)*facsk*(1/pcm(4) - pcm(4)/s)
         fsgam2j=fsgam2i
      endif

      if(mstc(111).eq.1) then
        call fengd(i,pkin1,A1)
        call fengd(j,pkin2,A2)
      if(mstc(114).eq.2) then

c....Derivatives of p^\mu/_i/p^0_i term in the vector part.
         fsgam2i=fsgam2i + mstc(116)*rhom(j,i)*
     &         ( A2(1)*bi(1) + A2(2)*bi(2) + A2(3)*bi(3) )/p(4,i)
         fsgam3i=-rhom(j,i)/p(4,i)

         fsgam2j=fsgam2j + mstc(116)*rhom(i,j)*
     &         ( A1(1)*bj(1) + A1(2)*bj(2) + A1(3)*bj(3) )/p(4,j)
         fsgam3j=-rhom(i,j)/p(4,j)

c...Derivative of p^\mu_j/m_j
      else if(mstc(114).eq.3) then

         fsgam2i=fsgam2i + mstc(116)*rhom(j,i)*A2(0)/emi
         fsgam3i=fsgam3i -rhom(j,i)/emi/p(4,i)
         fsgam2j=fsgam2j + mstc(116)*rhom(j,i)*A1(0)/emj
         fsgam3j=fsgam3j -rhom(i,j)/emj/p(4,j)

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
 
c....Momentum dependent potential part start.
      pma = mstc(113) * (emi**2- emj**2)**2 / s
      pk(4)=mstc(113)*(pc1(4)-pc2(4))
      p2 = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + pma
      fac1 = 1.d0 + p2/pmu1
      fac2 = 1.d0 + p2/pmu2

c...Derivative of p_{Tij} term.
      fengij=vf1*rhom(i,j) + vf2*rhom(j,i)
      fmome=-fengij*(vex1/pmu1/fac1**2+vex2/pmu2/fac2**2)

c...Derivative of rho_{ij} term
      facmom=(vex1/fac1+vex2/fac2)
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
      fmgam3i=0.0
      fmgam3j=0.0
c...Derivatives of gamma_{ij} part.
      if(mstc(114).eq.1) then
        fmgam1=fengij*facmom/s
        fmgam2i=mstc(116)*fengij*facmom*(1/pcm(4) - pcm(4)/s)
        fmgam2j=fmgam2i
      do n=1,3
        forcer(n,i) = forcer(n,i) + pcm(n)*fmgam1
     &                            + bi(n)*fmgam2i
 
        forcer(n,j) = forcer(n,j) + pcm(n)*fmgam1
     &                            + bj(n)*fmgam2j
      enddo

      endif

c...Derivative of  p_i/p^0_i part.
      if(mstc(111).eq.1.and.mstc(114).ne.0.and.mstc(115).eq.2) then
        facm1=facmom*rhom(j,i)
        facm2=facmom*rhom(i,j)
        if(mstc(114).eq.3) then
          fm1=-1.0/emi/pkin2(4)*facm1
          fm2=-1.0/emj/pkin1(4)*facm2
          fm3=mstc(116)*facm1/emi
          fm4=mstc(116)*facm2/emj
        else
          fm1=-facm1/p(4,i)/pkin2(4)
          fm2=-facm2/p(4,j)/pkin1(4)
c         vf2 = (pkin2(1)*bi(1)+pkin2(2)*bi(2)+pkin2(3)*bi(3))/pkin2(4)
c         vf1 = (pkin1(1)*bj(1)+pkin1(2)*bj(2)+pkin1(3)*bj(3))/pkin1(4)
          vf1 = 1.0 - vf1
          vf2 = 1.0 - vf2
          fm3=mstc(116)*facm1/p(4,i)*vf1
          fm4=mstc(116)*facm2/p(4,j)*vf2
        endif
        do n=1,3
          forcer(n,i) = forcer(n,i) + fm1*pkin2(n) + bi(n)*fm3
          forcer(n,j) = forcer(n,j) + fm2*pkin1(n) + bj(n)*fm4
        enddo
      endif

      end

**********************************************************************

      subroutine rqmddev2e(i,j,Ai,Aj,emi,emj,pkin1,pkin2,vfac)

c...Purpose: to compute derivatives of the squared four-vector distance
c...q_{Tij}^2 and p_{Tij}^2 in the rest frame of a particle i or j.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,n
      real*8 fengi,fengj,emi,emj
      real*8 dr2ri(3),dr2rj(3),dr2pi(3),dr2pj(3)
      real*8 dpij2pi(3),dpji2pi(3),dpij2pj(3),dpji2pj(3)
      real*8 fskyi,fskyj,fmomdi,fmomdj,fmomei,fmomej,pp2i,pp2j
      real*8 rk(3),pk(4),bi(3),bj(3),vfac
      real*8 rbi,rbj,pij,bb,bb1,bb2,wg,fac1,fac2,bij
      real*8 vi(0:4),vj(0:4),dvi(3),dvj(3)
      real*8 Ai(0:3),Aj(0:3),A1(0:3),A2(0:3)
      real*8 pkin1(5),pkin2(5),vf1,vf2
      real*8 fm1,fm2,fm3,fm4,facp1,facp2
      real*8 pcm(5),bet(3),p1(5),p2(5),bbi(3),bbj(3),s,pma,ppk

      p1(:)=p(:,i)
      p2(:)=p(:,j)
c...Factors for vector potential
      vi(0)=p1(4)/emi
      vj(0)=p2(4)/emj
      do n=1,3
        vi(n)=p1(n)/emi
        vj(n)=p2(n)/emj
        bi(n)=p1(n)/p1(4)
        bj(n)=p2(n)/p2(4)
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
      bij=k(9,i)*k(9,j)/9
      fskyi=-wg*vfac*fengi*rhom(i,j)*bij
      fskyj=-wg*vfac*fengj*rhom(j,i)*bij

c     fskyi=fskyi + rhc(i,j)  ! add Coulomb potential
c     fskyj=fskyj + rhc(j,i)  ! add Coulomb potential
c     fskyi=fskyi + rhs(i,j)  ! add symmetry energy
c     fskyj=fskyj + rhs(j,i)  ! add symmetry energy
c     fskyi=fskyi + rhy(i,j)  ! add Yukawa potential
c     fskyj=fskyj + rhy(j,i)  ! add Yukawa potential

      rk(1)=r(1,i)-r(1,j)
      rk(2)=r(2,i)-r(2,j)
      rk(3)=r(3,i)-r(3,j)
      rbi=rk(1)*vi(1)+rk(2)*vi(2)+rk(3)*vi(3)
      rbj=rk(1)*vj(1)+rk(2)*vj(2)+rk(3)*vj(3)
      do n=1,3
        dr2ri(n) =  2*(rk(n)+rbj*vj(n))   ! dR~^2_ij/dR_i
        dr2rj(n) =  2*(rk(n)+rbi*vi(n))   ! dR~^2_ji/dR_i
        dr2pi(n) =  2*rk(n)*rbi/emi       ! dR^2_ji/dP_i
        dr2pj(n) =  2*rk(n)*rbj/emj       ! dR^2_ij/dP_j
      end do

      do n=1,3
        force(n,i)  = force(n,i)  - fskyi*dr2ri(n) - fskyj*dr2rj(n)
        force(n,j)  = force(n,j)  + fskyi*dr2ri(n) + fskyj*dr2rj(n)
        forcer(n,i) = forcer(n,i) + fskyj*dr2pi(n)
        forcer(n,j) = forcer(n,j) + fskyi*dr2pj(n)
      end do

c...Derivative of p^\mu_j/m_j in the vector potential.
      if(mstc(111).eq.1) then
        call fengd(i,pkin1,A1)
        call fengd(j,pkin2,A2)
        do n=1,3
          dvi(n)=(mstc(116)*A2(0)*vi(n)/p(4,i)-A2(n)/emi)*rhom(j,i)*bij
          dvj(n)=(mstc(116)*A1(0)*vj(n)/p(4,j)-A1(n)/emj)*rhom(i,j)*bij
          forcer(n,i) = forcer(n,i) + dvi(n)
          forcer(n,j) = forcer(n,j) + dvj(n)
        end do
      endif

 
      ! Momentum dependent potential on.
      if(mstd(101).eq.0) return

c     vf1=1.0
c     vf2=1.0
      vf1=vj(0)
      vf2=vi(0)
c...V^mu(p) is fully included.
      if(mstc(115).eq.2) then
c       vf1=vj(0)
c       vf2=vi(0)
        do n=1,3
          vf1 = vf1 - pkin1(n)/pkin1(4)*vj(n)
          vf2 = vf2 - pkin2(n)/pkin2(4)*vi(n)
        end do
      endif

      pk(1)=p(1,i)-p(1,j)
      pk(2)=p(2,i)-p(2,j)
      pk(3)=p(3,i)-p(3,j)
      pk(4)=p(4,i)-p(4,j)
      ppk=pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2
      if(mstc(120).eq.1) then
        pcm(:)=p1(:)+p2(:)
        s=pcm(4)**2-pcm(1)**2-pcm(2)**2-pcm(3)**2
        pma =  (emi**2- emj**2)**2 / s
        pp2i = ppk + pma
        pp2j = pp2i
      else
        pij=p(4,i)*p(4,j)-p(1,i)*p(1,j)-p(2,i)*p(2,j)-p(3,i)*p(3,j)
        pp2j = ppk + (pij-emj**2)**2/emj**2
        pp2i = ppk + (pij-emi**2)**2/emi**2
      endif
c     facm1 = vf1*rhom(i,j)
c     facm2 = vf2*rhom(j,i)

      fac1 = 1.d0 + pp2j/pmu1
      fac2 = 1.d0 + pp2j/pmu2
      fmomei=-(vex1/pmu1/fac1**2+vex2/pmu2/fac2**2)*rhom(i,j)*vf1
      fmomdi=-(vex1/fac1+vex2/fac2)*wg*rhom(i,j)*vf1
      facp2 = (vex1/fac1+vex2/fac2)*rhom(i,j)/emj

      fac1 = 1.d0 + pp2i/pmu1
      fac2 = 1.d0 + pp2i/pmu2
      fmomej=-(vex1/pmu1/fac1**2+vex2/pmu2/fac2**2)*rhom(j,i)*vf2
      fmomdj=-(vex1/fac1+vex2/fac2)*wg*rhom(j,i)*vf2
      facp1 = (vex1/fac1+vex2/fac2)*rhom(j,i)/emi

c...Distance in the two-body cm is used for the argument of MD pot.
      if(mstc(120).eq.1) then

      do n=1,3
        bet(n)=pcm(n)/pcm(4)
        bbi(n)=bet(n)-mstc(116)*bi(n)
        bbj(n)=bet(n)-mstc(116)*bj(n)
      end do
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
        bb=mstc(116)*p(4,j)*bi(n)-p(n,j)
        dpij2pi(n) = 2*pk(n) - 2*mstc(116)*pk(4)*bi(n) + bb*bb2 ! dpij^2/dp_i
        dpji2pi(n) = 2*pk(n) - 2*mstc(116)*pk(4)*bi(n) + bb*bb1 ! dpji^2/pd_i

        bb=mstc(116)*p(4,i)*bj(n)-p(n,i)
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
c       facm1=-fmomdj/wg
c       facm2=-fmomdi/wg
        fm1=-1.0/pkin2(4)*facp1
        fm2=-1.0/pkin1(4)*facp2
        fm3=mstc(116)*facp1
        fm4=mstc(116)*facp2
        do n=1,3
          forcer(n,i) = forcer(n,i) + fm1*pkin2(n) + bi(n)*fm3
          forcer(n,j) = forcer(n,j) + fm2*pkin1(n) + bj(n)*fm4
        enddo
      endif

      end

**********************************************************************
      subroutine jamvpotdot4

c...Compute time-derivatives of the baryon current.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      real* 8 Ai(0:3),Aj(0:3),Bi(0:3),Bj(0:3),vi(0:3),vj(0:3)
      integer i,j,n,nn
      real*8 pcm(4),bet(3),rk(3),pk(5),bbi(3),bbj(3),emi,emj
      real*8 s,rbij,dr2ri(3),dr2rj(3),dr2pi(3),dr2pj(3)
      real*8 drji2ri(3),drji2rj(3),drji2pi(3),drji2pj(3)
      real*8 wg,dvi,dvj,vvi,vvj,vij,vji,doti,dotj,fai,faj,fvi,fvj
      real*8 xdoti,pdoti,xdotj,pdotj
      real*8 vdott(3),fai2,faj2,fai3,faj3,fvi2,fvj2
      real*8 dgami(3),dgamj(3),p1(5),p2(5),rbi,rbj,b1(3),b2(3)
      real*8 fmomdi,fmomdj,dp2ijpi(3),dp2ijpj(3),dp2jipi(3),dp2jipj(3)
      real*8 bb1,bb2,bb,p2i,p2j,psq1,psq2,fac1,fac2,fmomei,fmomej
      real*8 pma,pij

      real*8 vdot
      common/jamvdot/vdot(3,mxv)

c...Width parameter.
      wg = 1.0/(4.d0*pard(104))

      do i=1,nv
        if(MF_on(i).eq.0) cycle

        vvi = t1 + t3*rhog(i)  ! V_i/rho_i
        dvi=0.0
        if(abs(rhoi(i)).gt.1d-7) ! del(V_i/rho_i)/rho_i
     &  dvi = (pard(103)-1)*t3*rhoi(i)**(pard(103)-3)
        Bi(:)=dvi*rhoj(:,i)

c       emi=p(4,i)
c       if(mstc(114).eq.3.or.mstc(113).eq.2) then
c         emi=sqrt(p(4,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2)
c       endif

        p1=p(:,i)
        if(mstc(112).eq.1) call getkinmom2(i,p1)
        p1(5)=sqrt(p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2)
        emi=p1(4)
        if(mstc(114).eq.3.or.mstc(113).eq.2) emi=p1(5)
        
        vi(0)=p1(4)/emi
        vi(1)=p1(1)/emi
        vi(2)=p1(2)/emi
        vi(3)=p1(3)/emi
        b1(1)=p1(1)/p1(4)
        b1(2)=p1(2)/p1(4)
        b1(3)=p1(3)/p1(4)

      do j=i+1,nv
        if(MF_on(j)==0) cycle

        vvj = t1 + t3*rhog(j)
        dvj=0.0
        if(abs(rhoi(j)).gt.1d-7)
     &  dvj = (pard(103)-1)*t3*rhoi(j)**(pard(103)-3)
        Bj(:)=dvj*rhoj(:,j)

c       emj=p(4,j)
c       if(mstc(114).eq.3.or.mstc(113).eq.2) then
c         emj=sqrt(p(4,j)**2-p(1,j)**2-p(2,j)**2-p(3,j)**2)
c       endif

        p2=p(:,j)
        if(mstc(112).eq.1) call getmom5(j,p2)
        p2(5)=sqrt(p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2)
        emj=p2(4)
        if(mstc(114).eq.3.or.mstc(113).eq.2) emj=p2(5)

        vj(0)=p2(4)/emj
        vj(1)=p2(1)/emj
        vj(2)=p2(2)/emj
        vj(3)=p2(3)/emj
        b2(1)=p2(1)/p2(4)
        b2(2)=p2(2)/p2(4)
        b2(3)=p2(3)/p2(4)

        vji=vj(0)*rhoj(0,i)
        vij=vi(0)*rhoj(0,j)
        do n=1,3
          vji=vji-vj(n)*rhoj(n,i)
          vij=vij-vi(n)*rhoj(n,j)
        end do
        Ai(:)=vji*Bi(:)
        Aj(:)=vij*Bj(:)

       pk(:)=p1(:)-p2(:)

c...Compute derivatives of interaction densities.
      if(mstc(113).eq.1) then

      pcm(4)=p1(4)+p2(4)
      do n=1,3
        pcm(n)=p1(n)+p2(n)
        bet(n)=pcm(n)/pcm(4)
        rk(n)=r(n,i)-r(n,j)
        bbi(n)=bet(n)-mstc(116)*p1(n)/p1(4)
        bbj(n)=bet(n)-mstc(116)*p2(n)/p2(4)
      end do
      s=pcm(4)**2-pcm(1)**2-pcm(2)**2-pcm(3)**2
      rbij=(rk(1)*pcm(1)+rk(2)*pcm(2)+rk(3)*pcm(3))/s

      do n=1,3
        dr2ri(n)= 2*(rk(n)+rbij*pcm(n))               ! dR~^2_ij/dR_i
        dr2rj(n)= -dr2ri(n)
        dr2pi(n)=2*(rk(n) + pcm(4)*rbij*bbi(n))*rbij  ! dR~^2_ij/dP_i 
        dr2pj(n)=2*(rk(n) + pcm(4)*rbij*bbj(n))*rbij

        drji2ri(n)= 2*(rk(n)+rbij*pcm(n))
        drji2rj(n)= -drji2ri(n)
        drji2pi(n)=2*(rk(n) + pcm(4)*rbij*bbi(n))*rbij  ! dR~^2_ij/dP_i 
        drji2pj(n)=2*(rk(n) + pcm(4)*rbij*bbj(n))*rbij

      enddo

c...derivatives from gamma_{ij}.
      do n=1,3
        dgami(n)=mstc(116)*(1.0/pcm(4)-pcm(4)/s)*vi(n)+pcm(n)/s
        dgamj(n)=mstc(116)*(1.0/pcm(4)-pcm(4)/s)*vj(n)+pcm(n)/s
      end do

c...Non-relativistic.
      else if(mstc(113).eq.0) then

      dgami(:)=0.0
      dgamj(:)=0.0
      do n=1,3
        dr2ri(n)= 2*(r(n,i)-r(n,j))
        dr2rj(n)= -dr2ri(n)
        dr2pi(n)=0.0
        dr2pj(n)=0.0

        drji2ri(n)= dr2ri(n)
        drji2rj(n)= dr2rj(n)
        drji2pi(n)=0.0
        drji2pj(n)=0.0

      enddo

      else

        dgami(:)=0.0
        dgamj(:)=0.0
        rk(1)=r(1,i)-r(1,j)
        rk(2)=r(2,i)-r(2,j)
        rk(3)=r(3,i)-r(3,j)
        rbi=(rk(1)*vi(1)+rk(2)*vi(2)+rk(3)*vi(3))
        rbj=(rk(1)*vj(1)+rk(2)*vj(2)+rk(3)*vj(3))
        do n=1,3
          dr2ri(n) =  2*(rk(n)+rbj*vj(n))   ! 1/2*dR~^2_ij/dR_i
          dr2rj(n) =  -dr2ri(n)              ! 1/2*dR~^2_ij/dR_j
          dr2pi(n) =  0.0                    ! 1/2*dR~^2_ij/dP_i
          dr2pj(n) =  2*rk(n)*rbj/emj        ! 1/2*dR~^2_ij/dP_j

          drji2rj(n) = 2*(rk(n)+rbi*vi(n))
          drji2rj(n) =  -drji2ri(n)
          drji2pi(n) =  2*rk(n)*rbi/emi
          drji2pj(n) =  0.0
        end do

      endif

      xdoti=0.0
      pdoti=0.0
      xdotj=0.0
      pdotj=0.0
      do n=1,3
        xdoti=xdoti + (b1(n)+forcer(n,i))*dr2ri(n)
     &              + (b2(n)+forcer(n,j))*dr2rj(n)
        xdotj=xdotj + (b2(n)+forcer(n,j))*drji2rj(n)
     &              + (b1(n)+forcer(n,i))*drji2ri(n)
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
      else if(mstc(113).eq.1.or.mstc(120).eq.1) then
        pcm(4)=p1(4)+p2(4)
        do n=1,3
          pcm(n)=p1(n)+p2(n)
          bet(n)=pcm(n)/pcm(4)
          bbi(n)=bet(n)-mstc(116)*p1(n)/p1(4)
          bbj(n)=bet(n)-mstc(116)*p2(n)/p2(4)
        end do
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
      else
        pij=p1(4)*p2(4) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
c       pij=p(4,i)*p(4,j)-p(1,i)*p(1,j)-p(2,i)*p(2,j)-p(3,i)*p(3,j)
        p2i=(pij-emi**2)**2/emi**2
        p2j=(pij-emj**2)**2/emj**2
        psq2 = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + p2j
        psq1 = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + p2i
        bb1=2*(pij-emi**2)/emi**2
        bb2=2*(pij-emj**2)/emj**2
        do n=1,3
          bb=mstc(116)*p2(4)*b1(n)-p2(n)
          dp2ijpi(n) = 2*pk(n) - 2*mstc(116)*pk(4)*b1(n) + bb*bb2 ! dpij^2/dp_i
          dp2jipi(n) = 2*pk(n) - 2*mstc(116)*pk(4)*b1(n) + bb*bb1 ! dpji^2/pd_i

          bb=mstc(116)*p1(4)*b2(n)-p1(n)
          dp2ijpj(n) = -2*pk(n) + 2*mstc(116)*pk(4)*b2(n) + bb*bb2 ! dpij^2/dp_j
          dp2jipj(n) = -2*pk(n) + 2*mstc(116)*pk(4)*b2(n) + bb*bb1 ! dpji^2/pd_j
        end do
      endif
      fac1 = 1.d0 + psq1/pmu2
      fac2=  1.d0 + psq2/pmu2
      fmomdi = vex1/fac2 + vex2/fac2
      fmomdj = vex1/fac1 + vex2/fac1
      do n=1,3
        vdot(n,i) = vdot(n,i) + doti*fmomdi*vj(n)
        vdot(n,j) = vdot(n,j) + dotj*fmomdj*vi(n)
      end do
      pdoti=0.0
      pdotj=0.0
      do n=1,3
        pdoti=pdoti + force(n,i)*dp2ijpi(n) + force(n,j)*dp2ijpj(n)
        pdotj=pdotj + force(n,j)*dp2jipj(n) + force(n,i)*dp2jipi(n)
      end do
      fmomei = -pdoti*(vex1/pmu1/fac2**2 + vex2/pmu2/fac2**2)*rhom(i,j)
      fmomej = -pdotj*(vex1/pmu1/fac1**2 + vex2/pmu2/fac1**2)*rhom(j,i)
      do n=1,3
        vdot(n,i) = vdot(n,i) + fmomei*vj(n)
        vdot(n,j) = vdot(n,j) + fmomej*vi(n)
      end do

      endif  ! End momentum dependent potential part.


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

      end do
      end do

c     return

      vdott=0.0
      nn=0
      do i=1,nv
       if(MF_on(i).eq.0) cycle
       vdott(:) = vdott(:) + vdot(:,i)
       nn=nn+1
      end do
      do i=1,nv
       if(MF_on(i).eq.0) cycle
       vdot(:,i) = vdot(:,i) - vdott(:)/nn
      end do

      end

**********************************************************************
      subroutine kineticmomit4
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

      call jamrqmm4

      diff=0.0
      do i=1,nv
        if(MF_on(i).eq.0) cycle
        vv=k(9,i)/3*(t1+t3*rhog(i))
        pots(i) = 0.0
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

      subroutine jamsvpot4(difm,difv)

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
c     call jamrqmm4
c     call jamepart4(difm,difv)

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

         call jamrqmm4
         call jamepart4(difm,difv)
         call setenergy(iopt)
      end do
      write(6,800)it,difm,difv,difm+difv
3000  continue

800   format('jamsvpot4: not conserved diff',i3,3(e15.8,1x))
      end

**********************************************************************
      subroutine getkinmom2(i,pk)
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      integer i
      real*8 pk(5)

      pk(:)=p(:,i)
      pk(5)=p(5,i)
c     if(mstc(115).ne.2) return

      pk(1)=p(1,i)+potv(1,i)
      pk(2)=p(2,i)+potv(2,i)
      pk(3)=p(3,i)+potv(3,i)
      pk(4)=sqrt(pk(5)**2+pk(1)**2+pk(2)**2+pk(3)**2)

      end

************************************************************************
      subroutine rqmdpotparam4
      implicit none
      real*8 hc,rho0,conv
      real*8 alpha,beta,gam,C1,C2,mu1,mu2
      include 'jam2.inc'

      hc=paru(3)
      rho0=parc(21)
      conv=hc*hc*hc

      mstd(101)=1  ! Momentum dependent potential is used.

      if (mstc(106).eq.1) then ! Hard EoS K=380MeV
        alpha = -0.124838736062805
        beta  =  0.0707373967100532
        gam   =  2.00261915156202
        mstd(101)=0  ! momentum dependent potential is not used.
        C1=0.0
        C2=0.0
        mu1=1.0
        mu2=1.0

      elseif (mstc(106).eq.2) then ! Soft EoS K=210MeV
        alpha = -0.310514700691548
        beta  =  0.256413361338797
        gam   =  1.20292928384297
        C1=0.0
        C2=0.0
        mu1=1.0
        mu2=1.0
        mstd(101)=0  ! momentum dependent potential is not used.

      elseif (mstc(106).eq.3) then ! MH2 Skyrme Hard Nara2019 K=380
c    U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
c    U( 0.65 )= 0.0 Elab= 0.20320287416392357

        alpha = -0.00595607567906286
        beta = 0.08157802422694654
        gam = 1.718563134604995
        C1 = -0.38616561269531274
        mu1 = 2.02*hc
        C2= 0.34267680471490275
        mu2= 1.0*hc

      elseif (mstc(106).eq.4) then ! MS2 Skyrme Hard Nara2019 K=210

c   // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
c   // U( 0.65 )= 0.0 Elab= 0.20320287416392357
        alpha = -0.3085250472810656
        beta = 0.3804296791388347
        gam = 1.1142288606816626
        C1 = -0.38616561269531274
        C2= 0.34267680471490275
        mu1 = 2.02*hc
        mu2= 1.0*hc

      elseif (mstc(106).eq.5) then ! Nara2019 MH3 Hard EoS K=380MeV
        alpha = 0.031029748068250693
        beta = 0.04755488958383087
        gam = 2.1354953844773155
        C1 = -0.20192630016292473
        C2= 0.05694208245922472
        mu1 = 2.8*hc
        mu2= 1.0*hc
      
      ! 2020/1/14
      elseif (mstc(106).eq.6) then ! Nara2019 MS3 Soft EoS K=210MeV
        alpha = -0.8290165881158311
        beta = 0.907601225767913
        gam = 1.0386838048982703
        C1 = -0.20192630016293525
        C2= 0.05694208245927023

       !alpha = -0.6285006244768873
       !beta = 0.7070852697952635
       !gam = 1.0496537865787445
       !C1 = -0.20192630016292473
       !C2= 0.05694208245922472

       mu1 = 2.8*hc
       mu2= 1.0*hc

      elseif (mstc(106).eq.7) then ! MH4 Skyrme Hard Nara2019 K=380
        alpha = 0.03953021450755792
        beta = 0.040627668934799493
        gam = 2.2928923888385695
        C1 = -0.17086291974074935
        C2 = 0.0
        mu1 = 3.1466990715061636*hc
        mu2 = 1.0*hc

      elseif (mstc(106).eq.8) then ! MS4 Skyrme Soft Nara2019 K=210
        alpha = -0.22912691249560302
        beta = 0.30928479593796043
        gam = 1.1087616187247877
        C1 = -0.17086291974074935
        mu1 = 3.1466990715061636*hc
        C2 = 0.0
        mu2 = 1.0*hc

cc optical potential is defined by sqrt{(m_N+S)^2+p^2} - sqrt{m_N^2 + p^2}
      elseif (mstc(106).eq.11) then ! MH2 Skyrme Nara2021 K=380
c   // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
c   // U( 0.65 )= 0.0 Elab= 0.20320287416392357
        alpha = -0.013119515535259911
        beta = 0.08885442751779084
        gam = 1.674140687543709
        C1 = -0.3989613178044121
        C2= 0.36728513480692454
        mu1 = 2.02*hc
        mu2= 1.0*hc

      elseif (mstc(106).eq.12) then ! MS2 Skyrme Nara2021 K=210
c   // U( Pinf= 1.7 )= 0.06  Einf= 1.0036086114353737 Pinf= 1.7
c   // U( 0.65 )= 0.0 Elab= 0.20320287416392357
        alpha = -0.5157475041588349
        beta = 0.5906455475692723
        gam = 1.0708570434690778
        C1 = -0.3989613178044121
        C2= 0.36728513480692454
        mu1 = 2.02*hc
        mu2= 1.0*hc


      else
        write(6,*)'RQMDv mode wrong number mstc(106)=',mstc(106)
        stop
      endif

      if(mstc(109).ne.0) then
        write(6,*)' mstc(109) must be 0 in the RQMDv mode'
        mstc(109)=0
      endif

      pard(101)=alpha
      pard(102)=beta
      pard(103)=gam
      pard(105)= mu1
      pard(106)= mu2
      pard(107)= C1
      pard(108)= C2

      if(mstc(6).ge.2.and.mstc(119).eq.0.and.mstc(115).eq.2) then
        mstd(91)=1
      endif

      end
