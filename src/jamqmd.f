c***********************************************************************
c 2019/9/21                                                            *
c        PART  : Quantum Molecular Dynamics                            *
c                                                                      *
c   List of Subprograms in rough order of relevance with main purpose  *
c      (s = subroutine, f = function, b = block data, e = entry)       *
c  s jamqmd     to compute QMD force                                   *
c  s qmddev     to calculate derivatives                               *
c                                                                      *
************************************************************************

      subroutine jamqmd

c...Purpose: to calculate force in QMD
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,NrqmdCount,icheckMF
      real*8 dt,emfsq

c...Check which particle feels potentials.
      NrqmdCount=icheckMF(0,pard(1))

      call jamqmm
      call jamepart0
c     call setenergy(1)

c...Loop over all particles.
      forcer(:,:) = 0.0d0
      force(:,:) = 0.0d0
      rhog(:)=0.0
      do i=1,nv
        if(MF_on(i)==0) cycle
        if(rhoi(i).lt.1d-8) cycle
        rhog(i)=rhoi(i)**(pard(103)-1.0d0)
      end do
 
      do i=1,nv !...Loop over all particles i.
        if(MF_on(i).eq.0) cycle
        do j=i+1,nv !...Loop over particles j.
          if(MF_on(j)==0) cycle
          call qmddev(i,j)
        end do
      end do

c...Compute the time derivative of the vector potential to compute the kinetic
c...momentum in case vector potential is fully included.
      dt=parc(2)
      do i=1,nv
        if(MF_on(i)==0) cycle
        emfsq=p(4,i)**2 - p(1,i)**2 - p(2,i)**2 - p(3,i)**2
        p(1,i)=p(1,i)+dt*force(1,i)
        p(2,i)=p(2,i)+dt*force(2,i)
        p(3,i)=p(3,i)+dt*force(3,i)
        p(4,i)=sqrt(emfsq+p(1,i)**2+p(2,i)**2+p(3,i)**2)
c displacement only with interaction
        r(1,i)=r(1,i)+dt*forcer(1,i)
        r(2,i)=r(2,i)+dt*forcer(2,i)
        r(3,i)=r(3,i)+dt*forcer(3,i)
      enddo

c     call jamqmm
c     call jamepart0
c     call setenergy(1)

      end

**********************************************************************

      subroutine qmddev(i,j)

c...Purpose: to compute derivatives of the non-relativistic distance for
c...scalar Skyrme.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,n
      real*8 wg
      real*8 fsky1,fsky2,fsky
      real*8 p2,fac1,fac2,fengij,fmome,fmomd,facmom

c...Width parameter.
      wg = 1.0/(4.d0*pard(104))

c...Scalar potential.
      fsky1=(t1+t3f*rhog(i))*rhom(i,j)
      fsky2=(t1+t3f*rhog(j))*rhom(j,i)

c...Total.
      fsky= -wg*(fsky1 + fsky2)

c...Coulomb, Symmetry e, and Yukawa potentials are treated as vector.
c     fsky=fsky + rhc(i,j)+rhc(j,i)  ! Add Coulomb potential
c     fsky=fsky + rhs(i,j)+rhs(j,i)  ! Add symmetry energy
c     fsky=fsky + rhy(i,j)+rhy(j,i)  ! Add Yukawa potential

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
      fengij=rhom(i,j) + rhom(j,i)
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

      subroutine jamqmm

c...Purpose: to prepare matrix in calculating force
c...Option mstc(113)=2 is not implemented for Coulomb, Yukawa Symmetry force 

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j
      real*8 fac,wg,bi,bj,rs,den,ps,pmom2

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
        do 110 j = i+1 , nv      ! Sum for j (sigma_j(>i))
          if(MF_on(j).eq.0) goto 110
          bj=k(9,j)/3
          rs=(r(1,i)-r(1,j))**2+(r(2,i)-r(2,j))**2+(r(3,i)-r(3,j))**2
          den = exp(-rs*wg) /fac
          rhom(i,j) = den
          rhom(j,i) = den
          ps=(p(1,i)-p(1,j))**2+(p(2,i)-p(2,j))**2+(p(3,i)-p(3,j))**2
          pmom2 = vex1/(1.d0+ps/pmu1)+vex2/(1.d0+ps/pmu2)

          rhoi(i) = rhoi(i) + rhom(i,j)
          rhoi(j) = rhoi(j) + rhom(j,i)

          if (mstd(101).eq.0) goto 110
          vmoms(i) = vmoms(i) + pmom2*rhom(i,j)
          vmoms(j) = vmoms(j) + pmom2*rhom(j,i)

 110     continue               ! Loop end of sigma_j(>i)
 100  continue                  ! Loop end of sigma_i
 
      end

c***********************************************************************

      subroutine jamepart0 ! energy

c...Purpose: to calculate single particle potential energy in RQMD/S
c H = sigma_i=1^N  1/(2*E_i) *[E_i^2 -vec{p}_i^2 - m_i^2 - 2m_i*V_i ]
c here V_i = t1 * <rho_i> + t3 *<rho_i>^gamma : Skyrme Potential

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
c...Local Variable
      integer i
      real*8 pp2,vsky

      do 100 i=1,nv

       if(MF_on(i)==0) then
         pots(i)=0.0d0
         potv(:,i)=0.0d0
         goto 100
       endif

       potv(1,i)=0.0
       potv(2,i)=0.0
       potv(3,i)=0.0
       pots(i)=0.0
       vsky=k(9,i)/3*(t1*rhoi(i) + t3*rhoi(i)**pard(103))
       potv(0,i) = vsky + rhoc(i) + rhos(i) + rhoy(i) + vmom(0,i)

100   continue

      mste(44)=1
      end
