c***********************************************************************
c                                                                      *
c        PART  : RQMD/S Evolution                                      *
c                                                                      *
c   List of Subprograms in rough order of relevance with main purpose  *
c      (s = subroutine, f = function, b = block data, e = entry)       *
c  s jamrqmd1   to calculate force in RQMD/S                           *
c  s jamrqmm1   to prepare matrix in calculating force                 *
c                                                                      *
c***********************************************************************

**********************************************************************

      subroutine jamrqmd1(ind)

c...Purpose: to calculate force in RQMD/S
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      real*8 vdot
      common/jamvdot/vdot(3,mxv)

c...mstc(116)=1: Include the contribution of
c...{\partial p_i^0}{\partial p_i} = p_i/p_i^0
c...neglecting the potential term in p_i^0 in RQMD/S.

      integer ind,n0,nn,i,j,NrqmdCount,icheckMF,it
      real*8 diff,fengi,fengj,emisq,emjsq,feng
      real*8 dt,emfsq,diffv,difm,difv

c     t1=pard(101)/2.d0/parc(21)
c     t3f=pard(103)*pard(102)/(pard(103)+1.d0)/(parc(21)**pard(103))
c     wg = 1.0/(4.d0*pard(104))

      if (ind.eq.1) then         ! for target baryon
         n0 = 1
         nn = mstd(5)
      elseif (ind.eq.2) then     ! for projectile baryon
         n0 = mstd(5)+1
         nn = mstd(2)+mstd(5)
      elseif (ind.eq.0) then
         n0 = 1
         nn = nv
      endif

c...Check which particle feels potentials.
      NrqmdCount=icheckMF(ind,pard(1))

c....set p0 to be the free kinetic energy for the calculation of
c....potential.
      if(mstc(112).ge.1.and.mstc(109).ge.2) call setenergy(0)

      call jamrqmm1(ind)
      call jamepart(diff,diffv)
      if(mstc(112).eq.0) call setenergy(1)

c...Compute scalar and vector potential self-consistently.
      if(mstc(109).ge.3) then

        if(mstc(108).ge.1.and.mstc(115).ge.2) then
           call jamsvpot1(difm,difv)
        else

c....Effective mass is determined self-consistently.
        do it=1,30
           if(diff.le.1d-5) goto 3000
           call jamrqmm1(ind)
           call jamepart(diff,diffv)
           call setenergy(1)
        end do
           print *,'effective mass does not converge diff=',diff
3000  continue

      endif

      endif

c...Loop over all particles.
      forcer(:,:) = 0.0d0
      force(:,:) = 0.0d0
      vdot(:,:)=0.0d0
      rhog(:)=0.0
      if(t3.ne.0.0d0) then
      do 100 i=n0,nn
        if(MF_on(i)==0) goto 100
        if(rhoi(i).lt.1d-8) goto 100
        rhog(i)=rhoi(i)**(pard(103)-1.0d0)
 100  continue
      endif
 
c...Loop over all particles i.
      do 10 i=n0,nn  ! Sum for ith particle (sigma_i) 

         if(MF_on(i).eq.0) goto 10
         fengi=feng(i) ! Pre-factor from the scalar potential.
         emisq=p(4,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2

c...Loop over particles j.
         do 11 j=i+1,nn     ! Sum for j (sigma_j(not=i))
           if(MF_on(j)==0) goto 11
           fengj= feng(j)
           emjsq=p(4,j)**2-p(1,j)**2-p(2,j)**2-p(3,j)**2

          if(mstc(113).eq.0) then
            call rqmddev0(i,j,fengi,fengj)
          else if(mstc(113).le.1) then
            if(mstc(108).ge.1) then
              call rqmddev1b(i,j,fengi,fengj,emisq,emjsq)
            else
              call rqmddev1(i,j,fengi,fengj,emisq,emjsq)
            endif
          else
            call rqmddev2(i,j,fengi,fengj,emisq,emjsq)
          endif

 11      continue               ! Loop end of sigma_j [1,nn](j not= i)
 10   continue                  ! Loop end of sigma_i [1,nn]


c...Compute the time derivative of the vector potential to compute the kinetic
c...momentum in case vector potential is fully included.
      if(mstc(115).ge.2) call jamvpotdot4

c....set p0 using effective mass.
      if(mstc(112).ge.1) call setenergy(1)

c...[AO:030126
      if(ind.eq.0.and.mstc(103).eq.0) then
        dt=parc(2)
        do i=1,nv
          if(MF_on(i).eq.1) then
            emfsq=p(4,i)**2 - p(1,i)**2 - p(2,i)**2 - p(3,i)**2
            p(1,i)=p(1,i)+dt*(force(1,i) - vdot(1,i))
            p(2,i)=p(2,i)+dt*(force(2,i) - vdot(2,i))
            p(3,i)=p(3,i)+dt*(force(3,i) - vdot(3,i))
            p(4,i)=sqrt(emfsq+p(1,i)**2+p(2,i)**2+p(3,i)**2)
c displacement only with interaction
            r(1,i)=r(1,i)+dt*forcer(1,i)
            r(2,i)=r(2,i)+dt*forcer(2,i)
            r(3,i)=r(3,i)+dt*forcer(3,i)
          endif
        enddo

c       if(mstc(108).eq.1) call kineticmomit

        call jamrqmm1(ind)
        call jamepart(diff,diffv)

c....set p0 using effective mass.
      if(mstc(109).ge.2) call setenergy(1)

      endif
c...]AO:030126

c...[AO:050728
      if(ind.eq.0.and.NrqmdCount.gt.0.and.mstc(103).gt.0) then
        call jamdtfree
      endif
c...]AO:050728

      end

**********************************************************************

      subroutine jamrqmm1(ind)  ! makemat

c...Purpose: to prepare matrix in calculating force
c...Option mstc(113)=2 is not implemented for Coulomb, Yukawa Symmetry force 

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
c     logical jamrqpb
      integer ind,n0,nn,i,j,kc1,kc2,iz1,iz2,ic1,ic2,jamcomp
      real*8 fac,wg,wc,clw,gamy,yex,emisq,emjsq,emi,emj,ef1,ef2,bi,bj
      real*8 px,py,pz,pe,rx,ry,rz,s,r2,ps,p2,rs,erfr,dev1,dev2
      real*8 rhcpij,rhoyij,qfacr,cls,clsd,sym
      real*8 r2i,r2j,pij,p2i,p2j,gfac,den1,den2,pmom2ij,pmom2ji
      real*8 den,gam,gam2

      if (ind.eq.1) then         ! for target baryon
         n0 = 1
         nn = mstd(5)
      elseif (ind.eq.2) then     ! for projectile baryon
         n0 = mstd(5)+1
         nn = mstd(2)+mstd(5)
      elseif (ind.eq.0) then
         n0 = 1
         nn = nv
      endif

c     do i=n0,nn
c        MF_on(i)=1
c        if(jamrqpb(i)) MF_on(i)=0
c     enddo

      fac  =(4.d0*paru(1)*pard(104))**1.5d0 !    [(4*pi*L)^3/2]
      wg   = 1.0/4.d0/pard(104)
      wc   = sqrt(wg)                   ! for Coulomb potential
      clw  = 2.0/sqrt(4.0*paru(1)*pard(104)) ! for Coulomb potential
      gamy = parc(99)
      yex  = exp(pard(104)/gamy**2)     ! for Yukawa potential
      cls  = parc(100)/(2.0*parc(21))   ! for symmetry potential
      clsd = cls*wg                     ! for symmetry potential

      rho=0.d0   !...  <rho_i> = sigma_j(not= i) rho_ij
      rhoi=0.d0
      vmom=0.0
      rhoj=0.d0
      rhom=0.d0
      rhoc=0.d0  ! Coulomb
      rhc=0.d0   ! Coulomb
      rhos=0.d0  ! Symmetry energy
      rhs=0.d0   ! Symmetry energy
      rhoy=0.d0  ! Yukawa potential
      rhy=0.d0   ! Yukawa potential

      do 100 i=n0,nn              ! sum for i-th particle (sigma_i) 

        if(MF_on(i).eq.0) goto 100
        emisq=p(4,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2
        emi=1.0
        bi=k(9,i)/3
c       if(mstc(106).ge.201) emi=p(5,i)
        if(mstc(114).eq.4) emi=p(5,i)/p(4,i) ! Factor in the scalar density
c       if(mstc(108).eq.1.and.mstc(114).ne.0) emi=p(5,i)/p(4,i)
        if(mstc(108).eq.1.and.mstc(114).ne.0) emi=sqrt(emisq)/p(4,i)

c...For Coulomb potential.
        if(mstc(101).ge.1) then
          kc1=jamcomp(k(2,i))
          iz1=kchg(kc1,1)*isign(1,k(2,i))/3
        endif
        if(mstc(100).ge.1) then  ! Symmetry energy
          ic1=-1
          if(k(2,i).eq.2212) ic1=1
          if(k(2,i).eq.2112) ic1=0
        endif

        do 110 j = i+1 , nn      ! Sum for j (sigma_j(>i))
          if(MF_on(j).eq.0) goto 110
          emj=1.0
          bj=k(9,j)/3
          emjsq = p(4,j)**2 - p(1,j)**2 - p(2,j)**2 - p(3,j)**2
c         if(mstc(106).ge.201) emj=p(5,j)
          if(mstc(114).eq.4) emj=p(5,j)/p(4,j)
c         if(mstc(108).eq.1.and.mstc(114).ne.0) emj=p(5,j)/p(4,j)
          if(mstc(108).eq.1.and.mstc(114).ne.0) emj=sqrt(emjsq)/p(4,j)
          rx = r(1,i) - r(1,j)
          ry = r(2,i) - r(2,j)
          rz = r(3,i) - r(3,j)
          ps=(p(1,i)-p(1,j))**2+(p(2,i)-p(2,j))**2+(p(3,i)-p(3,j))**2
          rs= rx**2 + ry**2 + rz**2
          gfac=qfacr(i)*qfacr(j)
          gam2=1.0

c...Non-relativistic.
          if(mstc(113).eq.0) then

          den1 = exp(-rs*wg) * gfac
          den2 = den1
          rhom(i,j) = den1 /fac
          rhom(j,i) = den2 /fac
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
          r2= rs + mstc(113)*(rx*px + ry*py + rz*pz)**2/s
          den  = exp(-r2*wg) * gfac
          gam=1.0
          if(mstc(114).eq.1.or.mstc(114).eq.2.or.mstc(114).eq.4)
     &          gam = pe/sqrt(s)
          den1=den*gam
          den2=den*gam
          rhom(i,j) = den1 /fac !* emj*bj
          rhom(j,i) = den2 /fac !* emi*bi
          if(mstc(114).eq.2) gam2=gam

          p2 = ps - mstc(113)*((p(4,i)-p(4,j))**2 - (emisq- emjsq)**2/s)
          pmom2ij = vex1/(1.d0+p2/pmu1)+vex2/(1.d0+p2/pmu2)
          pmom2ji = pmom2ij
          p2i=p2
          p2j=p2

c...Two body distance is defined by the rest frame of i-th or j-th particle.
          else

          if(mstc(108).ge.1) then
            emi=1.0
            emj=1.0
          endif
          r2i=rs + (rx*p(1,i)+ry*p(2,i)+rz*p(3,i))**2/emisq
          r2j=rs + (rx*p(1,j)+ry*p(2,j)+rz*p(3,j))**2/emjsq
          den1= exp(-r2j*wg) * gfac
          den2= exp(-r2i*wg) * gfac
          if(mstc(114).eq.1.or.mstc(114).eq.2) then
           den1 = den1 * p(4,j)/sqrt(emjsq)
           den2=  den2 * p(4,i)/sqrt(emisq)
          endif
          rhom(i,j)= den1 /fac ! * emj*bj
          rhom(j,i)= den2 /fac ! * emi*bi
          pij=p(4,i)*p(4,j)-p(1,i)*p(1,j)-p(2,i)*p(2,j)-p(3,i)*p(3,j)
          p2i=ps - (p(4,i)-p(4,j))**2 + (pij-emisq)**2/emisq
          p2j=ps - (p(4,i)-p(4,j))**2 + (pij-emjsq)**2/emjsq
          pmom2ij = vex1/(1.d0+p2j/pmu1)+vex2/(1.d0+p2j/pmu2)
          pmom2ji = vex1/(1.d0+p2i/pmu1)+vex2/(1.d0+p2i/pmu2)

          endif

          if(mstc(118).ge.1) then
            rhom(i,j)=rhom(i,j)/(1.0+p2i/parc(105)**2)
            rhom(j,i)=rhom(j,i)/(1.0+p2j/parc(105)**2)
          endif

          rho(i) = rho(i) + rhom(i,j)*emj
          rho(j) = rho(j) + rhom(j,i)*emi

          rhoi(i) = rhoi(i) + rhom(i,j)*bj
          rhoi(j) = rhoi(j) + rhom(j,i)*bi

          vmom(0,i) = vmom(0,i) + pmom2ij*rhom(i,j)/gam2
          vmom(0,j) = vmom(0,j) + pmom2ji*rhom(j,i)/gam2

c....Coulomb potential
          if(mstc(101).ge.1) then
            kc2=jamcomp(k(2,j))
            iz2=kchg(kc2,1)*isign(1,k(2,j))/3  ! charge

c...Note! current version, r2 id not defined for mstc(113)=2
            rs   = sqrt( r2 )
            erfr  = erf( rs*wc ) / rs
            ! derivative
            dev1 = 0.25*parc(23)*iz1*iz2*( - erfr + clw * den1)/r2
            dev2 = 0.25*parc(23)*iz1*iz2*( - erfr + clw * den2)/r2
            rhc(i,j) = dev1 * emj
            rhc(j,i) = dev2 * emi
            rhcpij = 0.5*parc(23)*iz1*iz2*erfr
            rhoc(i) = rhoc(i) +  rhcpij
            rhoc(j) = rhoc(j) +  rhcpij
          endif

          if(mstc(100).ge.1) then  ! Symmetry energy
              ic2=-1
              if(k(2,j).eq.2212) ic2=1
              if(k(2,j).eq.2112) ic2=0
              if(ic1.ge.0.and.ic2.ge.0) then
                sym=(1.0-2.0*abs(ic1-ic2))*rhom(i,j)
                rhos(i)=rhos(i)+cls*sym
                rhos(j)=rhos(j)+cls*sym
                rhs(i,j)= -clsd*sym * emj    ! Derivative
                rhs(j,i)= -clsd*sym * emi
              endif
          endif

          if(mstc(99).ge.1) then  ! Yukawa potential
            rs   = sqrt( r2 )
            ef1=exp(-rs/gamy)*(1.0-erf(0.5/(wc*gamy)-rs*wc))
            ef2=exp( rs/gamy)*(1.0-erf(0.5/(wc*gamy)+rs*wc))
            rhoyij = 0.5*parc(98)*yex/rs*(ef1-ef2)
            rhoy(i) = rhoy(i) +  rhoyij
            rhoy(j) = rhoy(j) +  rhoyij
            ! Derivative
            dev1 = 0.25*parc(98)/r2*(
     &              yex*( (-1.0/gamy - 1.0/rs )*ef1
     &                   +(-1.0/gamy + 1.0/rs )*ef2 )
     &                + 2.0*clw*den1 )
            dev2 = 0.25*parc(98)/r2*(
     &              yex*( (-1.0/gamy - 1.0/rs )*ef1
     &                   +(-1.0/gamy + 1.0/rs )*ef2 )
     &                + 2.0*clw*den2 )
            rhy(i,j) = dev1 * emj
            rhy(j,i) = dev2 * emi
          endif

 110     continue               ! Loop end of sigma_j(>i)
 100  continue                  ! Loop end of sigma_i
 
c...Compute baryon current J_i^\mu
c     if(mstc(108).ge.1.or.(mstc(114).ge.2.and.mstc(114).le.4)) then
      if(mstc(108).ge.1.or.mstc(109).eq.0) then

      do 8 i=n0,nn
         if(MF_on(i).eq.0) goto 8
         emi=p(4,i)
         if(mstc(114).eq.3) ! emi=p(5,i)
     &   emi=sqrt(p(4,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2)
         do 9 j=i+1,nn
           if(MF_on(j).eq.0) goto 9
           emj=p(4,j)
           if(mstc(114).eq.3) !  emj=p(5,j)
     &     emj=sqrt(p(4,j)**2-p(1,j)**2-p(2,j)**2-p(3,j)**2)
            rhoj(0,i) = rhoj(0,i) + rhom(i,j)*p(4,j)/emj
            rhoj(1,i) = rhoj(1,i) + rhom(i,j)*p(1,j)/emj
            rhoj(2,i) = rhoj(2,i) + rhom(i,j)*p(2,j)/emj
            rhoj(3,i) = rhoj(3,i) + rhom(i,j)*p(3,j)/emj
            rhoj(0,j) = rhoj(0,j) + rhom(j,i)*p(4,i)/emi
            rhoj(1,j) = rhoj(1,j) + rhom(j,i)*p(1,i)/emi
            rhoj(2,j) = rhoj(2,j) + rhom(j,i)*p(2,i)/emi
            rhoj(3,j) = rhoj(3,j) + rhom(j,i)*p(3,i)/emi
 9       continue
 8    continue

      if(mstc(114).ne.0) then

c...Compute invariant baryon density.
        do i=n0,nn
          if(MF_on(i).eq.0) cycle
          rhoi(i)=sqrt(max(0d0,rhoj(0,i)**2
     &        -rhoj(1,i)**2-rhoj(2,i)**2-rhoj(3,i)**2))
c         if(rho(i).le.0d0) MF_on(i)=0
        enddo

      endif

      else

        rhoi=rho

      endif

c     call checkdens
c     call bcurrent

      end

**********************************************************************

      subroutine jamrqmd_simple(ind)
c     subroutine jamrqmd(ind)

c...Purpose: to calculate force in RQMD/S

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      logical jamrqpb
      real*8 bi(3),rk(3)
      real*8 dr2ri(3),dr2pi(3),dp2pi(3)
      real*8 dr2rj(3),dr2pj(3),dp2pj(3)
      real*8 pk(3),bet(3)

c...Loop over all particles.
      force=0.0
      forcer=0.0
      rho=0.0
      vmom=0.0
      rhom=0.0
      do i=1,nv
      MF_on(i)=1
      if(jamrqpb(i)) MF_on(i)=0
      end do

      fac =(4.d0*paru(1)*pard(104))**1.5d0 !    [(4*pi*L)^3/2]
      wg  = 1.0/4.d0/pard(104)
      do 200 i=1,nv              ! sum for i-th particle (sigma_i) 
        if(MF_on(i).eq.0) goto 200
        emisq = p(4,i)**2 - p(1,i)**2 - p(2,i)**2 - p(3,i)**2
        emi=1.0
        do 210 j = i+1 , nv      ! Sum for j (sigma_j(>i))
          if(MF_on(j).eq.0) goto 210
          emj=1.0
          px = p(1,i) + p(1,j)
          py = p(2,i) + p(2,j)
          pz = p(3,i) + p(3,j)
          pe = p(4,i) + p(4,j)
          rx = r(1,i) - r(1,j)
          ry = r(2,i) - r(2,j)
          rz = r(3,i) - r(3,j)
          s = pe**2 - px**2 - py**2 - pz**2
          r2 = rx**2+ry**2+rz**2 + mstc(113)*(rx*px+ry*py+rz*pz)**2/s
          den  = exp(-r2*wg)
          rhom(i,j) = den / fac * qfacr(i) * qfacr(j) 
          rhom(j,i) = rhom(i,j)
          rho(i) = rho(i) + rhom(i,j)
          rho(j) = rho(j) + rhom(i,j)
          emjsq = p(4,j)**2 - p(1,j)**2 - p(2,j)**2 - p(3,j)**2
          ps=(p(1,i)-p(1,j))**2+(p(2,i)-p(2,j))**2+(p(3,i)-p(3,j))**2
          p2 = ps - mstc(113)*((p(4,i)-p(4,j))**2 - (emisq- emjsq)**2/s)
          pmom2ij = vex1/(1.d0+p2/pmu1)+vex2/(1.d0+p2/pmu2)
          pmom2ji = pmom2ij
          vmom(0,i) = vmom(0,i) + pmom2ij*rhom(i,j)*emj
          vmom(0,j) = vmom(0,j) + pmom2ji*rhom(j,i)*emi
 210     continue               ! loopend of sigma_j(>i)
 200  continue                  ! loopend of sigma_i

c...Loop over all particles.
      rhog=0.0
      if(t3.ne.0d0) then
      do 300 i=1,nv
      if(MF_on(i)==1) then
        rhog(i)=max(0d0,rho(i))**(pard(103)-1.0d0)
      endif
 300  continue
      endif
 
c...Loop over all particles i.
      do 10 i=1,nv  ! sum for ith particle (sigma_i) 

         if(MF_on(i).eq.0) goto 10
         emisq=p(4,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2
         fengi= 1.0d0   ! In the case of vector potential.
         bi(1) = p(1,i)/p(4,i) ! P_xi/E_i
         bi(2) = p(2,i)/p(4,i) ! P_yi/E_i
         bi(3) = p(3,i)/p(4,i) ! P_zi/E_i

c...Loop over all particles j.
         do 11 j=i+1,nv     ! sum for j (sigma_j(not=i))
           if(MF_on(j)==0) goto 11
           emjsq=p(4,j)**2-p(1,j)**2-p(2,j)**2-p(3,j)**2
           fengj= 1.0d0

c...Compute derivatives of the squared four-vector distance q_{Tij}^2
c...and p_{Tij}^2.  mstc(113)=0: non-relativistic distance.
           feng = fengi+fengj    ! 2m_i/(2E_i)+2m_j/(2E_j)
           pe=p(4,i)+p(4,j)
           px=p(1,i)+p(1,j)
           py=p(2,i)+p(2,j)
           pz=p(3,i)+p(3,j)
           s=pe**2-px**2-py**2-pz**2
           bet(1)=px/pe
           bet(2)=py/pe
           bet(3)=pz/pe
           rk(1)=r(1,i)-r(1,j)
           rk(2)=r(2,i)-r(2,j)
           rk(3)=r(3,i)-r(3,j)
           rbij=mstc(113)*(rk(1)*bet(1)+rk(2)*bet(2)+rk(3)*bet(3))/s
           pma  = mstc(113) * (emisq- emjsq)**2 / s
           deng=mstc(113)*(p(4,i)-p(4,j)) ! E_i - E_j !m
           do n=1,3
             dr2ri(n)= rk(n)+rbij*bet(n)         ! 1/2*dR~^2_ij/dR_i
             dr2rj(n)=-dr2ri(n)
 
             bj=p(n,j)/p(4,j)                 ! P_j/E_j
             bbi=bet(n)-mstc(116)*bi(n)         ! beta_ij - P_i/E_i
             bbj=bet(n)-mstc(116)*bj            ! beta_ij - P_j/E_j

             ! 1/2*dR~^2_ij/dP_i 
             dr2pi(n)=rbij*(rk(n)+pe*rbij*bbi)
             dr2pj(n)=rbij*(rk(n)+pe*rbij*bbj)

             pk(n) = p(n,i)-p(n,j)
             dp2pi(n) = pk(n) - mstc(116)*deng*bi(n) + pe/s*pma*bbi
             dp2pj(n) =-pk(n) + mstc(116)*deng*bj    + pe/s*pma*bbj
           enddo

c... rho_ij = exp[-R~_ij^2/(4*L)]/[(4*pi*L)^3/2]
       fsky=feng*t1+t3f*(fengi*rhog(i) + fengj*rhog(j))
       fsky=(-0.25d0/pard(104))*rhom(i,j)*fsky

c....Momentum dependent potential.
       fmome=0.0 ! E_ij
       fmomd=0.0 ! D_ij
       if (mstd(101).eq.1) then
         ps = (p(1,i)-p(1,j))**2+(p(2,i)-p(2,j))**2+(p(3,i)-p(3,j))**2
         p2 = ps - deng**2 + pma
         fac1 = 1.d0 + p2/pmu1
         fac2 = 1.d0 + p2/pmu2
         fmome=-(vex1/pmu1/fac1**2+vex2/pmu2/fac2**2)*feng*rhom(i,j)
         fmomd=-(vex1/fac1+vex2/fac2)*feng*rhom(i,j)/(4*pard(104))
       endif

           do n=1,3
             forcer(n,i) = forcer(n,i) + 2*fsky*dr2pi(n) ! dR_i/dt
             force(n,i)  = force(n,i)  - 2*fsky*dr2ri(n) ! dP_i/dt
             forcer(n,j) = forcer(n,j) + 2*fsky*dr2pj(n) ! dR_j/dt
             force(n,j)  = force(n,j)  - 2*fsky*dr2rj(n) ! dP_j/dt
             forcer(n,i) = forcer(n,i) + 2*fmomd*dr2pi(n)
     &                                 + 2*fmome*dp2pi(n)
             force(n,i)  = force(n,i)  - 2*fmomd*dr2ri(n)
             forcer(n,j) = forcer(n,j) + 2*fmomd*dr2pj(n)
     &                                 + 2*fmome*dp2pj(n)
             force(n,j)  = force(n,j)  - 2*fmomd*dr2rj(n)
           enddo

 11      continue               ! loopend of sigma_j [1,nn](j not= i)
 10   continue                  ! loopend of sigma_i [1,nn]

      dt=parc(2)
      do i=1,nv
        if(MF_on(i).eq.1) then
          emfsq=p(4,i)**2 - p(1,i)**2 - p(2,i)**2 - p(3,i)**2
          p(1,i)=p(1,i)+dt*force(1,i)
          p(2,i)=p(2,i)+dt*force(2,i)
          p(3,i)=p(3,i)+dt*force(3,i)
          p(4,i)=sqrt(emfsq+p(1,i)**2+p(2,i)**2+p(3,i)**2)
c displacement only with interaction
          r(1,i)=r(1,i)+dt*forcer(1,i)
          r(2,i)=r(2,i)+dt*forcer(2,i)
          r(3,i)=r(3,i)+dt*forcer(3,i)
        endif
      enddo

      call jamrqmm1(ind)
      call jamepart(difm,difv)

      end

**********************************************************************
      real*8 function feng(i)
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      integer i
      real*8 effmsq,emff,emfsq

c...This is for jamrqmd

      feng= 1.0d0   ! In the case of vector potential.
      if(mstc(109).ge.1.or.mstc(108).ge.1) then
        feng= p(5,i)/p(4,i)  ! 2m_i/(2E_i)
        ! Effective mass is included only in this part:factor in
        ! front of derivative.
        if(mstc(109).ge.2) 
     &    feng=p(5,i)/sqrt(effmsq(i)+p(1,i)**2+p(2,i)**2+p(3,i)**2)
      endif
      if(mstc(106).ge.201) then
           emff=p(5,i)*(1d0+pots(i))
           feng= emff*p(5,i)/p(4,i)
      else if(mstc(106).ge.31) then
         emfsq=effmsq(i)
         feng=sqrt(emfsq)/sqrt(emfsq+p(1,i)**2+p(2,i)**2+p(3,i)**2)
      endif

      end


**********************************************************************
c....Pre-factor from the baryon current.
      real*8 function fengb(i,j)
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      integer i,j

      fengb=1.0
      if(rhoi(i).gt.0d0) then
        fengb=(rhoj(0,i)*p(4,j)
     &         -rhoj(1,i)*p(1,j)-rhoj(2,i)*p(2,j)-rhoj(3,i)*p(3,j))
     &          /(rhoi(i)*p(4,j))
      endif

      end

**********************************************************************
c....Pre-factor from the baryon current.
      real*8 function fengm(i,j,emj)
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      integer i,j
      real*8 emj

      fengm=1.0
      if(rhoi(i).gt.0d0) then
        fengm=(rhoj(0,i)*p(4,j)
     &         -rhoj(1,i)*p(1,j)-rhoj(2,i)*p(2,j)-rhoj(3,i)*p(3,j))
     &          /(rhoi(i)*emj)
      endif

      end

**********************************************************************

      subroutine rqmddev2(i,j,feng1,feng2,emisq,emjsq)

c...Purpose: to compute derivatives of the squared four-vector distance
c...q_{Tij}^2 and p_{Tij}^2 in the rest frame of a particle i or j.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i,j,n
      real*8 fengb,feng1,feng2,fengi,fengj,emisq,emjsq
      real*8 dr2ri(3),dr2rj(3),dp2pi(3),dp2pj(3)
      real*8 fskyi,fskyj,fmomdi,fmomdj,fmomei,fmomej,pp2i,pp2j
      real*8 rk(3),pk(4),bi(3),bj(3),emi,emj,fengm
      real*8 rbi,rbj,pij,p2i,p2j,bbi,wg,fac1,fac2,facmom

      fengi=feng1
      fengj=feng2
      if(mstc(114).eq.2) then
        fengi=feng1*fengb(i,j)
        fengj=feng2*fengb(j,i)
      else if(mstc(114).eq.3) then
        emi=sqrt(emisq)
        emj=sqrt(emjsq)
        fengi=feng1*fengm(i,j,emj)
        fengj=feng2*fengm(j,i,emi)
      endif

c...Skyrme potential.
      wg = 1.0/(4.d0*pard(104))
      fskyi= -wg*fengi*(t1+t3f*rhog(i))*rhom(i,j)
      fskyj= -wg*fengj*(t1+t3f*rhog(j))*rhom(j,i)

      fskyi=fskyi + fengi*rhc(i,j)  ! add Coulomb potential
      fskyj=fskyj + fengj*rhc(j,i)  ! add Coulomb potential
      fskyi=fskyi + fengi*rhs(i,j)  ! add symmetry energy
      fskyj=fskyj + fengj*rhs(j,i)  ! add symmetry energy
      fskyi=fskyi + fengi*rhy(i,j)  ! add Yukawa potential
      fskyj=fskyj + fengj*rhy(j,i)  ! add Yukawa potential

      rk(1)=r(1,i)-r(1,j)
      rk(2)=r(2,i)-r(2,j)
      rk(3)=r(3,i)-r(3,j)
      rbi=(rk(1)*p(1,i)+rk(2)*p(2,i)+rk(3)*p(3,i))/emisq
      rbj=(rk(1)*p(1,j)+rk(2)*p(2,j)+rk(3)*p(3,j))/emjsq
      do n=1,3
        dr2ri(n) =  rk(n)+rbj*p(n,j)    ! 1/2*dR~^2_ij/dR_i
        dr2rj(n) =  rk(n)+rbi*p(n,i)    ! 1/2*dR~^2_ij/dR_i
        force(n,i) = force(n,i) - 2*fskyi*dr2ri(n) - 2*fskyj*dr2rj(n) ! dP_i/dt
        force(n,j) = force(n,j) + 2*fskyi*dr2ri(n) + 2*fskyj*dr2rj(n) ! dP_j/dt
        forcer(n,i) = forcer(n,i) - 2*fskyj*rk(n)*rbi
        forcer(n,j) = forcer(n,j) - 2*fskyi*rk(n)*rbj
      end do

      ! Momentum dependent potential on.
      if(mstd(101).eq.0) return

      pk(1)=p(1,i)-p(1,j)
      pk(2)=p(2,i)-p(2,j)
      pk(3)=p(3,i)-p(3,j)
      pk(4)=p(4,i)-p(4,j)
      pij=p(4,i)*p(4,j)-p(1,i)*p(1,j)-p(2,i)*p(2,j)-p(3,i)*p(3,j)
      p2i=(pij-emisq)**2/emisq
      p2j=(pij-emjsq)**2/emjsq

      pp2j = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + p2j
      fac1 = 1.d0 + pp2j/pmu1
      fac2 = 1.d0 + pp2j/pmu2
      facmom = -fengi*rhom(i,j)
      fmomei=(vex1/pmu1/fac1**2+vex2/pmu2/fac2**2)*facmom
      fmomdi=(vex1/fac1+vex2/fac2)*wg*facmom

      pp2i = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + p2i
      fac1 = 1.d0 + pp2i/pmu1
      fac2 = 1.d0 + pp2i/pmu2
      facmom = -fengj*rhom(j,i)
      fmomej=(vex1/pmu1/fac1**2+vex2/pmu2/fac2**2)*facmom
      fmomdj=(vex1/fac1+vex2/fac2)*wg*facmom

      do n=1,3
        bi(n)=p(n,i)/p(4,i)
        bj(n)=p(n,j)/p(4,j)
        bbi=mstc(116)*bi(n)-bj(n)
c       bbj=mstc(116)*bj(n)-bi(n)
        dp2pi(n) = pk(n) - mstc(116)*pk(4)*bi(n) + p(4,j)*p2j*bbi
        dp2pj(n) =-pk(n) - mstc(116)*pk(4)*bi(n) + p(4,j)*p2i*bbi
        forcer(n,i) = forcer(n,i) + 2*fmomei*dp2pi(n)
        forcer(n,j) = forcer(n,j) + 2*fmomej*dp2pj(n)
        force(n,i) = force(n,i) - 2*fmomdi*dr2ri(n) - 2*fmomdj*dr2rj(n)
        force(n,j) = force(n,j) + 2*fmomdi*dr2ri(n) + 2*fmomdj*dr2rj(n)
      enddo

      end

**********************************************************************

      subroutine rqmddev1(i,j,feng1,feng2,emisq,emjsq)

c...Purpose: to compute derivatives of the squared four-vector distance
c...q_{Tij}^2 and p_{Tij}^2 in the two body c.m. of particle i and j.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      real*8 feng1,feng2,fengi,fengj,emisq,emjsq
      integer i,j,n
      real*8 bi(3),bj(3)
      real*8 dr2ri(3),dr2pi(3),dp2pi(3)
      real*8 dr2rj(3),dr2pj(3),dp2pj(3)
      real*8 pc(4),bet(3),rk(3),pk(4)
      real*8 s,rbij,pma,fengb,fengij,fengip,fengjp,fsky1,fsky2,fsky
      real*8 fsgam1,fsgam2i,fsgam2j,fmgam1,facsk
      real*8 fsgam3i,fsgam3j,fmgam2i,fmgam2j,fmgam3i,fmgam3j
      real*8 p2,fac1,fac2,facmom,fmomd,fmome
      real*8 bbi(3),bbj(3),wg,gam,emi,emj,fengm


      pc(4)=p(4,i)+p(4,j)
      do n=1,3
        pc(n)=p(n,i)+p(n,j)
        bet(n)=pc(n)/pc(4)
        rk(n)=r(n,i)-r(n,j)
        pk(n)=p(n,i)-p(n,j)
        bi(n)=p(n,i)/p(4,i)
        bj(n)=p(n,j)/p(4,j)
        bbi(n)=bet(n)-mstc(116)*bi(n)
        bbj(n)=bet(n)-mstc(116)*bj(n)
      end do

      s=pc(4)**2-pc(1)**2-pc(2)**2-pc(3)**2
      rbij=mstc(113)*(rk(1)*pc(1)+rk(2)*pc(2)+rk(3)*pc(3))/s
      do n=1,3
        dr2ri(n)= rk(n)+rbij*pc(n)         ! 1/2*dR~^2_ij/dR_i
        dr2rj(n)= -dr2ri(n)
        dr2pi(n)=(rk(n) + pc(4)*rbij*bbi(n))*rbij  ! 1/2*dR~^2_ij/dP_i 
        dr2pj(n)=(rk(n) + pc(4)*rbij*bbj(n))*rbij
      enddo

      fengi=feng1
      fengj=feng2
      fengip=feng1
      fengjp=feng2
      gam=1.0
      if(mstc(114).eq.2) then
        fengi=feng1*fengb(i,j)
        fengj=feng2*fengb(j,i)
c...Baryon current is not used in the momentum dependent potential.
        gam=pc(4)/sqrt(s)
      else if(mstc(114).eq.3) then
        emi=sqrt(emisq)
        emj=sqrt(emjsq)
        fengi=feng1*fengm(i,j,emj)
        fengj=feng2*fengm(j,i,emi)
      endif

c...Skyrme potential.
      wg = 1.0/(4.d0*pard(104))
      fsky1= (t1+t3f*rhog(i))*rhom(i,j)
      fsky2= (t1+t3f*rhog(j))*rhom(j,i)
      facsk= fengi*fsky1+fengj*fsky2
      fsky= -wg*facsk

      fsky=fsky + fengi*rhc(i,j)+fengj*rhc(j,i)  ! Add Coulomb potential
      fsky=fsky + fengi*rhs(i,j)+fengj*rhs(j,i)  ! Add symmetry energy
      fsky=fsky + fengi*rhy(i,j)+fengj*rhy(j,i)  ! Add Yukawa potential

      fsgam1=0.0  ! p_i(n) + p_j(n) term
      fsgam2i=0.0 ! bi(n) term
      fsgam2j=0.0 ! bi(n) term
      fsgam3i=0.0 ! j_i(n) term
      fsgam3j=0.0 ! j_i(n) term
c...for now option mstc(114)=1,4 is only for Skyrme + mom. dep. force.

c.....Derivatives of gamma_{ij}
      if(mstc(114).ge.1.and.mstc(114).le.2) then
         fsgam1= facsk/s
         fsgam2i=mstc(116)*facsk*(1/pc(4) - pc(4)/s)
         fsgam2j=fsgam2i
      endif

c....Derivatives of p^\mu/_i/p^0_i term.
      if(mstc(114).eq.2) then

         if(rho(j).gt.0d0) then
         fsgam2i=fsgam2i + mstc(116)*fsky2*
     &         (rhoj(0,j)/rho(j) - fengj)/p(4,i)
         fsgam3i=-fsky2/(rho(j)*p(4,i))
         endif
         if(rho(i).gt.0d0) then
         fsgam2j=fsgam2j + mstc(116)*fsky1*
     &         ( rhoj(0,i)/rho(i) - fengi )/p(4,j)
         fsgam3j=-fsky1/(rho(i)*p(4,j))
         endif

c...Derivative of p^\mu_j/m_j
      else if(mstc(114).eq.3) then

         if(rho(j).gt.0d0) then
         fsgam2i=fsgam2i + mstc(116)*fsky2*rhoj(0,j)/rho(j)/emi
         fsgam3i=-fsky2/(rho(j)*emi)
         endif
         if(rho(i).gt.0d0) then
         fsgam2j=fsgam2j + mstc(116)*fsky1*rhoj(0,i)/rho(i)/emj
         fsgam3j=-fsky1/(rho(i)*emj)
         endif

c...Derivative of m_j/p^0_j
      else if(mstc(114).eq.4) then
         fsgam2i=fsgam2i-mstc(116)*fsky2*fengj/p(4,i)
         fsgam2j=fsgam2j-mstc(116)*fsky1*fengi/p(4,j)
      endif

      do n=1,3
        forcer(n,i) = forcer(n,i) + 2*fsky*dr2pi(n) ! dR_i/dt
     &           + pc(n)*fsgam1 + bi(n)*fsgam2i + rhoj(n,j)*fsgam3i
        forcer(n,j) = forcer(n,j) + 2*fsky*dr2pj(n) ! dR_j/dt
     &           + pc(n)*fsgam1 + bj(n)*fsgam2j + rhoj(n,i)*fsgam3j
        force(n,i) = force(n,i) - 2*fsky*dr2ri(n) ! dP_i/dt
        force(n,j) = force(n,j) - 2*fsky*dr2rj(n) ! dP_j/dt
      end do

c...Done if momentum dependent potential is off.
      if (mstd(101).eq.0) return
 
c....Momentum dependent potential part start.
      pma = mstc(113) * (emisq- emjsq)**2 / s
      pk(4)=mstc(113)*(p(4,i)-p(4,j))
      p2 = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + pma
      fac1 = 1.d0 + p2/pmu1
      fac2 = 1.d0 + p2/pmu2

c...Derivative of p_{Tij} term.
      fengij=(fengip*rhom(i,j)+fengjp*rhom(j,i))/gam
      fmome=-fengij*(vex1/pmu1/fac1**2+vex2/pmu2/fac2**2)

c...Derivative of rho_{ij} term
      facmom=(vex1/fac1+vex2/fac2)
      fmomd=-wg*facmom*fengij

      do n=1,3
        dp2pi(n) = pk(n) - mstc(116)*pk(4)*bi(n) + pc(4)/s*pma*bbi(n)
        dp2pj(n) =-pk(n) + mstc(116)*pk(4)*bj(n) + pc(4)/s*pma*bbj(n)
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
        fmgam2i=mstc(116)*fengij*facmom*(1/pc(4) - pc(4)/s)
        fmgam2j=fmgam2i

c...Invariant baryon density.
c     else if(mstc(114).eq.2) then
c       fmgam1=fengij*facmom/s
c       fmgam2i=mstc(116)*fengij*facmom*(1/pc(4) - pc(4)/s)
c       fmgam2j=fmgam2i

c...Derivative of  m_j/p^0_j part.
      else if(mstc(114).eq.4) then
        fmgam2i=fmgam2i-mstc(116)*facmom*fengjp/p(4,i)*rhom(j,i)
        fmgam2j=fmgam2j-mstc(116)*facmom*fengip/p(4,j)*rhom(i,j)
      endif

      do n=1,3
        forcer(n,i) = forcer(n,i) + pc(n)*fmgam1
     &                            + bi(n)*fmgam2i
c    &                            + rhoj(n,j)*fmgam3i
 
        forcer(n,j) = forcer(n,j) + pc(n)*fmgam1
     &                            + bj(n)*fmgam2j
c    &                            + rhoj(n,i)*fmgam3j
      enddo

      end

**********************************************************************

      subroutine rqmddev1b(i,j,feng1,feng2,emisq,emjsq)

c...Purpose: to compute derivatives of the squared four-vector distance
c...q_{Tij}^2 and p_{Tij}^2 in the two body c.m. of particle i and j.
c...Here attractive part of the Skyrme potential is treated as Lorentz scalar,
c...repulsive part is treated as the Lorentz vector.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      real*8 feng1,feng2,fengi,fengj,emisq,emjsq
      integer i,j,n
      real*8 bi(3),bj(3)
      real*8 dr2ri(3),dr2pi(3),dp2pi(3)
      real*8 dr2rj(3),dr2pj(3),dp2pj(3)
      real*8 pc(4),bet(3),rk(3),pk(4)
      real*8 s,rbij,pma,fengb,fengij,fengip,fengjp,fsky1,fsky2,fsky
      real*8 fsgam1,fsgam2i,fsgam2j,fmgam1,facsk
      real*8 fsgam3i,fsgam3j,fmgam2i,fmgam2j,fmgam3i,fmgam3j
      real*8 p2,fac1,fac2,facmom,fmomd,fmome
      real*8 bbi(3),bbj(3),wg,gam,emi,emj,fengm
      real*8 fsky1b,fsky2b

      pc(4)=p(4,i)+p(4,j)
      do n=1,3
        pc(n)=p(n,i)+p(n,j)
        bet(n)=pc(n)/pc(4)
        rk(n)=r(n,i)-r(n,j)
        pk(n)=p(n,i)-p(n,j)
        bi(n)=p(n,i)/p(4,i)
        bj(n)=p(n,j)/p(4,j)
        bbi(n)=bet(n)-mstc(116)*bi(n)
        bbj(n)=bet(n)-mstc(116)*bj(n)
      end do

      s=pc(4)**2-pc(1)**2-pc(2)**2-pc(3)**2
      rbij=mstc(113)*(rk(1)*pc(1)+rk(2)*pc(2)+rk(3)*pc(3))/s
      do n=1,3
        dr2ri(n)= rk(n)+rbij*pc(n)         ! 1/2*dR~^2_ij/dR_i
        dr2rj(n)= -dr2ri(n)
        dr2pi(n)=(rk(n) + pc(4)*rbij*bbi(n))*rbij  ! 1/2*dR~^2_ij/dP_i 
        dr2pj(n)=(rk(n) + pc(4)*rbij*bbj(n))*rbij
      enddo

c...Factors for vector potential
      fengi=1.0
      fengj=1.0
      if(mstc(114).eq.2) then
        fengi=fengb(i,j)
        fengj=fengb(j,i)
      else if(mstc(114).eq.3) then
        emi=sqrt(emisq)
        emj=sqrt(emjsq)
        fengi=fengm(i,j,emj)
        fengj=fengm(j,i,emi)
      endif

      gam=1.0
      if(mstc(114).eq.1.or.mstc(114).eq.2) then
        gam=pc(4)/sqrt(s)
      endif

c...Skyrme potential.
      wg = 1.0/(4.d0*pard(104))
      fsky1= t3f*rhog(i)*rhom(i,j)
      fsky2= t3f*rhog(j)*rhom(j,i)
      fsky1b= feng1*t1*rhom(i,j)/gam + fengi*fsky1
      fsky2b= feng2*t1*rhom(j,i)/gam + fengj*fsky2
      fsky= -wg*(fsky1b+fsky2b)

c...Coulomb, Symmetry e, and Yukawa potentials are treated as vector.
      if(mstc(109).eq.0) then
        fsky=fsky + rhc(i,j)+rhc(j,i)  ! Add Coulomb potential
        fsky=fsky + rhs(i,j)+rhs(j,i)  ! Add symmetry energy
        fsky=fsky + rhy(i,j)+rhy(j,i)  ! Add Yukawa potential
c...Coulomb, Symmetry e, and Yukawa potentials are treated as scalar.
      else
        fsky=fsky + feng1*rhc(i,j)+feng2*rhc(j,i)  ! Add Coulomb potential
        fsky=fsky + feng1*rhs(i,j)+feng2*rhs(j,i)  ! Add symmetry energy
        fsky=fsky + feng1*rhy(i,j)+feng2*rhy(j,i)  ! Add Yukawa potential
      endif

      fsgam1=0.0  ! p_i(n) + p_j(n) term
      fsgam2i=0.0 ! bi(n) term
      fsgam2j=0.0 ! bi(n) term
      fsgam3i=0.0 ! j_i(n) term
      fsgam3j=0.0 ! j_i(n) term
c...for now option mstc(114)=2,3 is only for Skyrme + mom. dep. force.

c.....Derivatives of gamma_{ij}
      if(mstc(114).ge.1.and.mstc(114).le.2) then
         facsk = fengi*fsky1 + fengj*fsky2
         fsgam1= facsk/s
         fsgam2i=mstc(116)*facsk*(1/pc(4) - pc(4)/s)
         fsgam2j=fsgam2i
      endif

c....Derivatives of p^\mu/_i/p^0_i term.
      if(mstc(114).eq.2) then

         if(rhoi(j).gt.0d0) then
         fsgam2i=fsgam2i + mstc(116)*fsky2*
     &         (rhoj(0,j)/rhoi(j) - fengj)/p(4,i)
         fsgam3i=-fsky2/(rhoi(j)*p(4,i))
         endif
         if(rhoi(i).gt.0d0) then
         fsgam2j=fsgam2j + mstc(116)*fsky1*
     &         ( rhoj(0,i)/rhoi(i) - fengi )/p(4,j)
         fsgam3j=-fsky1/(rhoi(i)*p(4,j))
         endif

c...Derivative of p^\mu_j/m_j
      else if(mstc(114).eq.3) then

         if(rhoi(j).gt.0d0) then
         fsgam2i=fsgam2i + mstc(116)*fsky2*rhoj(0,j)/rhoi(j)/emi
         fsgam3i=-fsky2/(rhoi(j)*emi)
         endif
         if(rhoi(i).gt.0d0) then
         fsgam2j=fsgam2j + mstc(116)*fsky1*rhoj(0,i)/rhoi(i)/emj
         fsgam3j=-fsky1/(rhoi(i)*emj)
         endif

c...Derivative of m_j/p^0_j
      else if(mstc(114).eq.4) then
         fsgam2i=fsgam2i-mstc(116)*fsky2*fengj/p(4,i)
         fsgam2j=fsgam2j-mstc(116)*fsky1*fengi/p(4,j)
      endif

      do n=1,3
        forcer(n,i) = forcer(n,i) + 2*fsky*dr2pi(n) ! dR_i/dt
     &           + pc(n)*fsgam1 + bi(n)*fsgam2i + rhoj(n,j)*fsgam3i
        forcer(n,j) = forcer(n,j) + 2*fsky*dr2pj(n) ! dR_j/dt
     &           + pc(n)*fsgam1 + bj(n)*fsgam2j + rhoj(n,i)*fsgam3j
        force(n,i) = force(n,i) - 2*fsky*dr2ri(n) ! dP_i/dt
        force(n,j) = force(n,j) - 2*fsky*dr2rj(n) ! dP_j/dt
      end do

c...Done if momentum dependent potential is off.
      if (mstd(101).eq.0) return
 
c....Momentum dependent potential part start.
      pma = mstc(113) * (emisq- emjsq)**2 / s
      pk(4)=mstc(113)*(p(4,i)-p(4,j))
      p2 = pk(1)**2 + pk(2)**2 + pk(3)**2 - pk(4)**2 + pma
      fac1 = 1.d0 + p2/pmu1
      fac2 = 1.d0 + p2/pmu2

c...Factors for momentum dependent potential.
      if(mstc(109).eq.0) then
        fengip=1.0
        fengjp=1.0
      else
        fengip=feng1
        fengjp=feng2
      endif

c...Derivative of p_{Tij} term.
c...Baryon current is not used in the momentum dependent potential:
c...cancel the gamma factor in front of Gaussian rhom(i,j).
      fengij=(fengip*rhom(i,j)+fengjp*rhom(j,i))/gam
      fmome=-fengij*(vex1/pmu1/fac1**2+vex2/pmu2/fac2**2)

c...Derivative of rho_{ij} term
      facmom=(vex1/fac1+vex2/fac2)
      fmomd=-wg*facmom*fengij

      do n=1,3
        dp2pi(n) = pk(n) - mstc(116)*pk(4)*bi(n) + pc(4)/s*pma*bbi(n)
        dp2pj(n) =-pk(n) + mstc(116)*pk(4)*bj(n) + pc(4)/s*pma*bbj(n)
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
        fmgam2i=mstc(116)*fengij*facmom*(1/pc(4) - pc(4)/s)
        fmgam2j=fmgam2i

c...Invariant baryon density.
c     else if(mstc(114).eq.2) then
c       fmgam1=fengij*facmom/s
c       fmgam2i=mstc(116)*fengij*facmom*(1/pc(4) - pc(4)/s)
c       fmgam2j=fmgam2i

c...Derivative of  m_j/p^0_j part.
      else if(mstc(114).eq.4) then
        fmgam2i=fmgam2i-mstc(116)*facmom*fengjp/p(4,i)*rhom(j,i)
        fmgam2j=fmgam2j-mstc(116)*facmom*fengip/p(4,j)*rhom(i,j)
      endif

      do n=1,3
        forcer(n,i) = forcer(n,i) + pc(n)*fmgam1
     &                            + bi(n)*fmgam2i
c    &                            + rhoj(n,j)*fmgam3i
 
        forcer(n,j) = forcer(n,j) + pc(n)*fmgam1
     &                            + bj(n)*fmgam2j
c    &                            + rhoj(n,i)*fmgam3j
      enddo

      end


**********************************************************************

      subroutine rqmddev0(i,j,feng1,feng2)

c...Purpose: to compute the force in the non-rel. QMD model.
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      real*8 feng1,feng2,fengi,fengj
      integer i,j,n
      real*8 bi(3),bj(3)
      real*8 fsky1,fsky2,fsky
      real*8 fengb,fengij,fengip,fengjp,fsgam2i,fsgam2j
      real*8 fmgam2i,fmgam2j,fsgam3i,fsgam3j
      real*8 p2,fac1,fac2,facmom,fmomd,fmome
      real*8 wg,emi,emj,feng3,feng4,fengm

      fengi=feng1
      fengj=feng2
      fengip=feng1
      fengjp=feng2
c...Pre-factors from the baryon current.
      if(mstc(114).eq.2) then
        fengi=feng1*fengb(i,j)
        fengj=feng2*fengb(j,i)
      else if(mstc(114).eq.3) then
        fengi=feng1*fengm(i,j,p(5,j))
        fengj=feng2*fengm(j,i,p(5,i))
      endif

c...Skyrme potential.
      wg = 1.0/(4.d0*pard(104))
      fsky1= (t1+t3f*rhog(i))*rhom(i,j)
      fsky2= (t1+t3f*rhog(j))*rhom(j,i)
      fsky= -wg*(fengi*fsky1+fengj*fsky2)
      fsky=fsky + fengi*rhc(i,j)+fengj*rhc(j,i)  ! Add Coulomb potential
      fsky=fsky + fengi*rhs(i,j)+fengj*rhs(j,i)  ! Add symmetry energy
      fsky=fsky + fengi*rhy(i,j)+fengj*rhy(j,i)  ! Add Yukawa potential
      do n=1,3
        force(n,i) = force(n,i) - 2*fsky*(r(n,i)-r(n,j))
        force(n,j) = force(n,j) - 2*fsky*(r(n,j)-r(n,i))
        bi(n)=p(n,i)/p(4,i)
        bj(n)=p(n,j)/p(4,j)
      end do

c....Derivatives of p^\mu/_i/p^0_i term in the baryon current.
      if(mstc(114).eq.2.or.mstc(114).eq.3) then
        if(mstc(114).eq.2) then
          emi=p(4,i)
          emj=p(4,j)
          feng3=fengi
          feng4=fengj
        else
          emi=p(5,i)
          emj=p(5,j)
          feng3=0.0
          feng4=0.0
        endif
        if(rho(j).gt.0d0) then
        fsgam2i = mstc(116)*fsky2*(rhoj(0,j)/rho(j) - feng4)/emi
        fsgam3i = -fsky2/(rho(j)*emi)
        endif
        if(rho(i).gt.0d0) then
        fsgam2j = mstc(116)*fsky1*(rhoj(0,i)/rho(i) - feng3)/emj
        fsgam3j = -fsky1/(rho(i)*emj)
        endif
        do n=1,3
          forcer(n,i) = forcer(n,i) + bi(n)*fsgam2i + rhoj(n,j)*fsgam3i
          forcer(n,j) = forcer(n,j) + bj(n)*fsgam2j + rhoj(n,i)*fsgam3j
        end do

c...Derivative of m_j/p^0_j in the scalar density.
      else if(mstc(114).eq.4) then
        fsgam2i = -mstc(116)*fsky2*fengj/p(4,i)
        fsgam2j = -mstc(116)*fsky1*fengi/p(4,j)
        do n=1,3
          forcer(n,i) = forcer(n,i) + bi(n)*fsgam2i
          forcer(n,j) = forcer(n,j) + bj(n)*fsgam2j
        end do
      endif

c...Done if momentum dependent potential is off.
      if (mstd(101).eq.0) return
 
c....Momentum dependent potential part start.
      p2 = (p(1,i)-p(1,j))**2 + (p(2,i)-p(2,j))**2 + (p(3,i)-p(3,j))**2
      fac1 = 1.d0 + p2/pmu1
      fac2 = 1.d0 + p2/pmu2

c...p Derivative term.
      fengij=fengip*rhom(i,j) + fengjp*rhom(j,i)
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
      if(mstc(114).eq.4) then
        fmgam2i=-mstc(116)*facmom*fengjp/p(4,i)*rhom(j,i)
        fmgam2j=-mstc(116)*facmom*fengip/p(4,j)*rhom(i,j)
        do n=1,3
          forcer(n,i) = forcer(n,i) + p(n,i)/p(4,i)*fmgam2i
          forcer(n,j) = forcer(n,j) + p(n,j)/p(4,j)*fmgam2j
        enddo
      endif

      end

**********************************************************************

      subroutine jamsvpot1(difm,difv)

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

      difm=0.0
      difv=0.0
      do it=1,50
        if(it.ge.2.and.difm+difv.le.1d-5) goto 3000

        if(mstc(112).eq.1) then
          do i=1,nv
            if(MF_on(i)==0) cycle
             p(1,i)=pkin(1,i) + potv(1,i)
             p(2,i)=pkin(2,i) + potv(2,i)
             p(3,i)=pkin(3,i) + potv(3,i)
             p(4,i)=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
          end do
         endif

         call jamrqmm1(0)
         call jamepart(difm,difv)
         call setenergy(iopt)
      end do
      write(6,800)it,difm,difv,difm+difv
3000  continue

800   format('jamsvpot1: not conserved diff',i3,3(e15.8,1x))
      end

c***********************************************************************

      subroutine checkdens

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      integer i,j
      real*8 den1,den2,den3,den4,den5,den6,den7
      real*8 rx,ry,rz,rs,r2,px,py,pz,pe,s,fac,wg

      fac  =(4.d0*paru(1)*pard(104))**1.5d0 !    [(4*pi*L)^3/2]
      wg   = 1.0/4.d0/pard(104)

      do i=1,3
        den1=0.0
        den2=0.0
        den3=0.0
        den4=0.0
        den5=0.0
        den6=0.0
        den7=0.0
        do j=1,nv
          if(i.eq.j) cycle
          rx = r(1,i) - r(1,j)
          ry = r(2,i) - r(2,j)
          rz = r(3,i) - r(3,j)
          rs= rx**2 + ry**2 + rz**2
          den4=den4+exp(-rs*wg) /fac

          r2= rs+(rx*rhoj(1,i)+ry*rhoj(2,i)+rz*rhoj(3,i))**2/rho(i)**2
          den5=den5+exp(-r2*wg) /fac * rhoj(0,i)/rho(i)
          den1=den1+exp(-r2*wg) /fac

          px = p(1,i) + p(1,j)
          py = p(2,i) + p(2,j)
          pz = p(3,i) + p(3,j)
          pe = p(4,i) + p(4,j)
          s = pe**2 - px**2 - py**2 - pz**2
          r2= rs+(rx*px+ry*py+rz*pz)**2/s
          den2=den2+exp(-r2*wg) /fac
          den6=den6+exp(-r2*wg) /fac * pe/sqrt(s)

          px =  p(1,j)
          py =  p(2,j)
          pz =  p(3,j)
          pe =  p(4,j)
          s = pe**2 - px**2 - py**2 - pz**2

          r2= rs+(rx*px+ry*py+rz*pz)**2/s
          den3=den3+exp(-r2*wg) /fac
          den7=den7+exp(-r2*wg) /fac * pe/sqrt(s)
        end do
        print *,i,'rho_inv=',rho(i),'rho_0=',rhoj(0,i)
        print *,'rest frame         =',den1,'w gam=',den5
        print *,'two-body c.m.farame=',den2,'w gam=',den6
        print *,'rest frame of part =',den3,'w gam=',den7
        print *,'nonrel             =',den4
      end do
c       read(5,*)

      end

************************************************************************
      subroutine rqmdpotparam1
      implicit none
      real*8 hc,rho0,conv,ck,be,emnucl,pf,ef
      real*8 alpha,beta,gam,C1,C2,mu1,mu2
      include 'jam2.inc'

      hc=paru(3)
      rho0=parc(21)
      conv=hc*hc*hc

      pard(105)=1.0  ! pmu1 dummy
      pard(106)=1.0  ! pmu2 dummy
      pard(107)=0.0d0   ! vex1
      pard(108)=0.0d0   ! vex2

      mstd(101)=1  ! Momentum dependent potential is used.

      if (mstc(106).eq.-3) then !
c       CK=mstc(106)*0.001       ! Compressibility 
c       BE = -0.0158   ! Binding energy
        CK=parc(104)             ! Compressibility 
        BE = parc(22)  ! Binding energy
        emnucl=0.938   ! Nucleon mass
        pf = (1.5*paru(1)**2*rho0)**(1.0/3.0)*hc
        Ef=3.0/5.0*pf*pf/(2*emnucl)
        pard(103)=(CK/9.0 + 2.0/9.0*Ef)/(Ef/3.0 - BE)
        pard(102)=(Ef/3.0 - BE)*(pard(103)+1.0)/(pard(103)-1.0)
        pard(101)=2*( BE - Ef - pard(102)/(pard(103)+1.0) )
        mstd(101)=0  ! Momentum dependent potential is not used.
        print *,'K alpha beta gamma=',ck,pard(101),pard(102),pard(103)

      elseif (mstc(106).eq.-2) then ! test momentum dependent part.
        pard(101)= 0.0
        pard(102)= 0.0
        pard(103)= 1.0             ! gamma
        pard(105)= 2.02*hc             ! mu1
        pard(106)= 1.0*hc              ! mu2
        pard(107)= -0.38314            ! C1
        pard(108)=  0.33741            ! C2

      else if (mstc(106).eq.-1) then ! Free hadron gas
        pard(101)=0.0
        pard(102)= 0.0
        pard(103)= 1.d0
        pard(105)=1.0  ! pmu1 dummy
        pard(106)=1.0  ! pmu2 dummy
        pard(107)= 0.0d0   ! vex1
        pard(108)= 0.0d0   ! vex2
        mstd(101)=0  ! momentum dependent potential is not used.
      else if (mstc(106).eq.0) then
c       pard(101)=0.45*rho0 ! EoS-Q
        pard(101)=0.54*rho0*2 ! EoS-Q
        pard(102)= 0.0
        pard(103)= 1.d0
        pard(105)=1.0  ! pmu1 dummy
        pard(106)=1.0  ! pmu2 dummy
        pard(107)= 0.0d0   ! vex1
        pard(108)= 0.0d0   ! vex2
        mstd(101)=0  ! momentum dependent potential is not used.

      elseif (mstc(106).eq.1) then
        pard(101)=- 33.d-3 ! New Hard 2005
        pard(102)= 110.d-3
        pard(103)=   1.66667d0
c       pard(104)=   2.05d0
        pard(105)=2.35d0*0.19732705d0  ! pmu1
        pard(106)=0.4d0*0.19732705d0    ! pmu2
        pard(107)= -0.277d0   ! vex1
        pard(108)= 0.663d0    ! vex2
      elseif (mstc(106).eq.2) then  
        pard(101)=-116.d-3 ! New Medium 2005
        pard(102)= 193.d-3
        pard(103)=   1.33333d0
c       pard(104)=   2.1d0
        pard(105)=2.35d0*0.19732705d0  ! pmu1
        pard(106)=0.4d0*0.19732705d0    ! pmu2
        pard(107)= -0.277d0   ! vex1
        pard(108)= 0.663d0    ! vex2
      elseif (mstc(106).eq.3) then
        pard(101)=-268.d-3 ! New Soft  2005
        pard(102)= 345.d-3
        pard(103)=   1.16667d0
c       pard(104)=   2.1d0
        pard(105)=2.35d0*0.19732705d0  ! pmu1
        pard(106)=0.4d0*0.19732705d0    ! pmu2
        pard(107)= -0.277d0   ! vex1
        pard(108)= 0.663d0    ! vex2

      elseif (mstc(106).eq.4) then
        pard(101)=-124.d-3 ! Aich Hard
        pard(102)=  70.5d-3
        pard(103)=   2.d0
c       pard(104)=   1.08d0
        pard(105)=1.0  ! pmu1 dummy
        pard(106)=1.0  ! pmu2 dummy
        pard(107)= 0.0d0   ! vex1
        pard(108)= 0.0d0   ! vex2
        mstd(101)=0  ! momentum dependent potential is not used.

      elseif (mstc(106).eq.5) then
        pard(101)=-356.d-3 ! Aich Soft
        pard(102)= 303.d-3
        pard(103)=   1.16667d0
c       pard(104)=   1.08d0
        pard(105)=1.0  ! pmu1 dummy
        pard(106)=1.0  ! pmu2 dummy
        pard(107)= 0.0d0   ! vex1
        pard(108)= 0.0d0   ! vex2
        mstd(101)=0  ! momentum dependent potential is not used.

      elseif (mstc(106).eq.6) then ! GBD
        pard(101)= -0.1449
        pard(102)=  0.2033
        pard(103)=   7.0d0/6.0d0
c       pard(104)=   2.1d0
        pard(105)= 0.267*1.5*0.19732705d0  ! pmu1
        pard(106)= 1.0
        pard(107)= -0.075
        pard(108)= 0.0d0
      elseif (mstc(106).eq.7) then ! Maruyama soft K=210MeV
        pard(101)= -0.22356
        pard(102)=  0.29878
        pard(103)=   1.16667
c       pard(104)=   2.1d0
        pard(105)=2.35d0*0.19732705d0   ! pmu1
        pard(106)=0.4d0*0.19732705d0    ! pmu2
        pard(107)= -0.25854
        pard(108)= 0.3756
      elseif (mstc(106).eq.8) then ! GiBUU soft K=215MeV
        pard(101)= -0.1086
        pard(102)=  0.1368
        pard(103)=   1.26
c       pard(104)=   2.1d0
        pard(105)=2.13d0*0.19732705d0   ! pmu1
        pard(106)=1.0d0                 ! pmu2
        pard(107)= -0.0636*2
        pard(108)= 0.0

      elseif (mstc(106).eq.9) then
c       A0=-87.67
c       chi = 0.93
        pard(101)= -0.93*2.0*87.67*1d-3 ! urqmd Hard
        pard(102)=  125.98*1d-3
        pard(103)=  1.676
        pard(105)=1.0  ! pmu1 dummy
        pard(106)=1.0  ! pmu2 dummy
        pard(107)= 0.0d0   ! vex1
        pard(108)= 0.0d0   ! vex2
        mstd(101)=0  ! momentum dependent potential is not used.


c.....2015 A.Ohnishi new fit1 mstc(106)=11-14
      elseif (mstc(106).eq.11) then ! MEH1 K=436.5
        pard(101)= 0.00960346793263952 ! alpha
        pard(102)= 0.0655471283733455  ! beta
        pard(103)= 2.0d0               ! gamma
c       pard(104)= 2.05d0              ! The width of Gaussian L (fm^2)
        pard(105)= 2.02*hc             ! mu1
        pard(106)= 1.0*hc              ! mu2
        pard(107)= -0.38314            ! C1
        pard(108)=  0.33741            ! C2
      elseif (mstc(106).eq.12) then ! MH1 K=370.919367532666
        pard(101)= -0.0122455748584756 ! alpha
        pard(102)=  0.0873961711644606 ! beta
        pard(103)= 5.0/3.0             ! gamma
c       pard(104)= 2.05d0              ! The width of Gaussian L (fm^2)
        pard(105)= 2.02*hc             ! mu1
        pard(106)= 1.0*hc              ! mu2
        pard(107)= -0.38314            ! C1
        pard(108)=  0.33741            ! C2
      elseif (mstc(106).eq.13) then ! MM1 K= 305.37223915932
        pard(101)=-0.0777927032318212  ! alpha
        pard(102)= 0.152943299537806   ! beta
        pard(103)= 4.0/3.0             ! gamma
c       pard(104)= 2.05d0              ! The width of Gaussian L (fm^2)
        pard(105)= 2.1*hc              ! mu1
        pard(106)= 1.0*hc              ! mu2
        pard(107)= -0.38314            ! C1
        pard(108)=  0.33741            ! C2
      elseif (mstc(106).eq.14) then ! MS1 K=272.598674972647
        pard(101)= -0.208886959978512  ! alpha
        pard(102)=  0.284037556284497  ! beta
        pard(103)= 7.0/6.0             ! gamma

c       pard(104)= 2.1d0               ! The width of Gaussian L (fm^2)
c       pard(104)= 1.08d0               ! The width of Gaussian L (fm^2)

        pard(105)= 2.02*hc             ! mu1
        pard(106)= 1.0*hc              ! mu2
        pard(107)= -0.38314            ! C1
        pard(108)=  0.33741            ! C2
c.....2015 A.Ohnishi new fit2 mstc(106)=21-24
      elseif (mstc(106).eq.21) then ! MEH2 K=466.821092509211
        pard(101)= -0.00353218032039811 ! alpha
        pard(102)=  0.0693982863782082  ! beta
        pard(103)= 2.0d0               ! gamma
c       pard(104)= 2.05d0              ! The width of Gaussian L (fm^2)
        pard(105)= 1.55*hc             ! mu1
        pard(106)= 1.0*hc              ! mu2
        pard(107)= -0.70388            ! C1
        pard(108)=  0.73311            ! C2
      elseif (mstc(106).eq.22) then ! MH1 K=397.422806131003
        pard(101)= -0.0266649424464675 ! alpha
        pard(102)= 0.0925310485042777  ! beta
        pard(103)= 5.0/3.0             ! gamma
c       pard(104)= 2.05d0              ! The width of Gaussian L (fm^2)
        pard(105)= 1.55*hc             ! mu1
        pard(106)= 1.0*hc              ! mu2
        pard(107)= -0.70388            ! C1
        pard(108)=  0.73311            ! C2
      elseif (mstc(106).eq.23) then ! MM1 K=328.024519752795
        pard(101)= -0.0960632288246758 ! alpha
        pard(102)=  0.161929334882486  ! beta
        pard(103)= 4.0/3.0             ! gamma
c       pard(104)= 2.05d0              ! The width of Gaussian L (fm^2)
        pard(105)= 1.55*hc             ! mu1
        pard(106)= 1.0*hc              ! mu2
        pard(107)= -0.70388            ! C1
        pard(108)=  0.73311            ! C2
      elseif (mstc(106).eq.24) then ! MS1 K=293.325376563691 
        pard(101)= -0.234859801581092  ! alpha
        pard(102)= 0.300725907638902   ! beta
        pard(103)= 7.0/6.0             ! gamma
c       pard(104)= 2.1d0               ! The width of Gaussian L (fm^2)
        pard(105)= 1.55*hc             ! mu1
        pard(106)= 1.0*hc              ! mu2
        pard(107)= -0.70388            ! C1
        pard(108)=  0.73311            ! C2

c...Danielewicz NPA673(2000)375
      elseif (mstc(106).eq.51) then ! Soft K=380 MeV
        pard(101)= -0.12126
        pard(102)=  0.0521
        pard(103)=  2.4624
c       pard(104)=   1.08d0
        pard(105)=1.0  ! pmu1 dummy
        pard(106)=1.0  ! pmu2 dummy
        pard(107)= 0.0d0   ! vex1
        pard(108)= 0.0d0   ! vex2
        mstd(101)=0  ! momentum dependent potential is not used.
      elseif (mstc(106).eq.52) then ! Soft K=210 MeV
        pard(101)= -0.18724
        pard(102)=  0.10262
        pard(103)=  1.6339
c       pard(104)=   1.08d0
        pard(105)=1.0  ! pmu1 dummy
        pard(106)=1.0  ! pmu2 dummy
        pard(107)= 0.0d0   ! vex1
        pard(108)= 0.0d0   ! vex2
        mstd(101)=0  ! momentum dependent potential is not used.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c...Under construction. Mean field with phase transition. 
      elseif (mstc(106).eq.201) then
        pard(101)=parc(101)*2.0d0*rho0
        pard(103)=parc(103)
        pard(102)=parc(102)*(rho0**pard(103))*(pard(103)+1.0d0)
        mstc(109)=2
      elseif (mstc(106).eq.202) then
        pard(101)=0.0d0
c       pard(103)=1.0/3.0
c       pard(102)=-0.3*(rho0**pard(103))*(pard(103)+1.0d0)
        pard(103)=parc(103)
        pard(102)=parc(101)*(rho0**pard(103))*(pard(103)+1.0d0)
        mstd(101)=0  ! momentum dependent potential is not used.
        mstc(109)=2
      else
        write(6,*)'RQMD/S wrong mode number mstc(106)=',mstc(106)
        stop
      endif

      end


