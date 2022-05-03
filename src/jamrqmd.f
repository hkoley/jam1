**********************************************************************
! YN  last modified 2017.01.26 bug fix in fsky
! YN  last modified 2015.08.11
! YN new param mstd(101)=1:mom.dep.pot. 2015.04.17
! YN revised mstc(109) 2015.02.20
! YN mstc(198)->mstc(106) mstc(199)->mstc(107) 2014.11.23 YN
! rel. mom-dep incl 2002.09.12
! AO revised        2002.08.30
! introduced mom-dep2002.08.23
! attempt mom-dep   2002.07.24 
! Otuka-san revised 2002.02.14
! Cooling succeed   2002.01.06
! Compile succeed   2001.10.29
! Included V        2001.10.26 
! Otuka-san revised 2001.10.10
**********************************************************************
c***********************************************************************
c                                                                      *
c        PART  : RQMD/S Evolution                                      *
c                                                                      *
c   List of Subprograms in rough order of relevance with main purpose  *
c      (s = subroutine, f = function, b = block data, e = entry)       *
c  f jamrqpb    to judge potential act on it or not                    *
c  s jamrqmd    to calculate force in RQMD/S                           *
c  s jamepart   to calculate single particle energy
c  s jamrqen    to calculate energy in RQMD/S                          *
c                                                                      *
c  s jamrqch    cool or heat proj./targ. to fit binding energy         *
c  s jamrqpt                                                           *
c  s findbeng                                                          *
c  f beld                                                              *
c  f clust                                                             *
c  s jamcluster                                                        *
c***********************************************************************

      subroutine jamrqmd(ind)

c...Purpose: to calculate force in RQMD
      implicit none
      include 'jam2.inc'
      integer ind


      if(mstc(108).eq.0) then      ! non-relativistic QMD
        call jamqmd
      else if(mstc(108).eq.1) then ! Skyrme force RQMD/S original
        call jamrqmd1(ind)
      else if(mstc(108).eq.2) then ! Skyrme force RQMD/S rewritten
        call jamrqmd2
      else if(mstc(108).eq.3) then !  Scalar Skyrme force
        call jamrqmd3
      else if(mstc(108).eq.4) then !  Vector Skyrme force
        call jamrqmd4
      else if(mstc(108).eq.5) then ! Sigma-omega model new
        call jamrqmd5
      else
        print *,'jamrqmd::invalid value mstc(108)=',mstc(108)
        stop
      endif

      end

**********************************************************************
      integer function icheckMF(iopt,gtime)

      implicit none
      include 'jam1.inc'
      include 'jam3.inc'
      logical jamrqpb
      integer NrqmdCount,i,iopt,MFon
      real*8 gtime
      save NrqmdCount
      data NrqmdCount/0/

      if(gtime.eq.0.0d0) then
        NrqmdCount=0
        do i=1,mxv
          dtfree(i)=0.0d0
        enddo
      else
        NrqmdCount=NrqmdCount+1
      endif

      if(iopt.eq.0.and.NrqmdCount.gt.0) then
        do i=1,nv
           MFon=MF_on(i)
           MF_on(i)=1
           if(jamrqpb(i)) MF_on(i)=0
           if(MFon.eq.1.or.MF_on(i).eq.0) then
             dtfree(i)=0.0d0
           endif
           if(MF_on(i).eq.1.and.i.gt.mxw) then
             print *,'mean-field particle index exceed mxw',i
             print *,'k1=',k(1,i),'kf=',k(2,i),'m=',p(5,i)
             MF_on(i)=0
           endif
        enddo
      endif

      icheckMF=NrqmdCount

      end

c***********************************************************************
      subroutine jamdtfree

      implicit none
      include 'jam1.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i
      real*8 dt,emfsq

        do i=1,nv
          if(MF_on(i).eq.1.and.dtfree(i).gt.0) then
            dt=dtfree(i)
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
        do i=1,mxv
          dtfree(i)=0.0d0
        enddo

      end

c***********************************************************************
c...Note used
      logical function jamrqpb2(kf)   ! judge

c...Purpose: to judge potential should
c              act on i-th particle (false)
c       or not act on i-th particle (true)
c mstc(104) :(D=1) option for RQMD/S transport.
c =0 : Potential effects are counted only for formed nucleons.
c =1 : Potential effects are counted only for formed baryons.
c =10: Potential effects are counted for nucleons (even before formation).
c =11: Potential effects are counted for baryons (even before formation).
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      integer kf

      jamrqpb2=.true.
      if(abs(kf).le.100) return
      if(mstc(104).lt.13) then
          if(mod(iabs(kf),10000)/1000.eq.0) return ! Non-Baryons
      endif
      if(mod(mstc(104),10).eq.0
     &  .and.kf.ne.2212.and.kf.ne.2112) return  ! Non-Nucleons
      jamrqpb2=.false.

      end


c***********************************************************************

      logical function jamrqpb(i)   ! judge

c...Purpose: to judge potential should
c              act on i-th particle (false)
c       or not act on i-th particle (true)
c mstc(104) :(D=1) option for RQMD/S transport.
c =0 : Potential effects are counted only for formed nucleons.
c =1 : Potential effects are counted only for formed baryons.
c =10: Potential effects are counted for nucleons (even before formation).
c =11: Potential effects are counted for baryons (even before formation).
c =12: Potential effects are counted for baryons which has original
c......constituent quarks (even before formation).
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      integer i,k1,kf,ibar,iq
      real*8 tf,t,teps
      parameter(teps=1.0d-8)
 
      jamrqpb=.true.
      k1=k(1,i)
      if(k1.gt.10) return
      kf=k(2,i)
      if(abs(kf).le.100) return  ! Exclude leptons etc.
      ibar=k(9,i)
      tf=r(5,i)
      t=pard(1)+teps
c     t=mstd(23)*parc(2)+teps

c     if(mod(mstc(104),10).eq.0) then
      if(mstc(104).eq.0 .or. mstc(104).eq.10) then

        if(kf.ne.2212.and.kf.ne.2112) return  ! Only nucleons feel pot.
        if(mstc(104).eq.0 .and. tf.ge.t) return

c....Baryons feel potentials.
      else if(mstc(104).eq.1 .or. mstc(104).eq.11) then

        if(ibar.eq.0) return ! Exclude mesons.

c.....only for formed baryons.
        if(mstc(104).eq.1 .and. (tf.ge.t)) return

c...only formed delta and nucleon
      else if(mstc(104).eq.2) then
        if(tf.ge.t) return
        if(kf.ne.2212.and.kf.ne.2112.and.
     &  kf.ne.1114.and.kf.ne.2114.and.kf.ne.2214.and.kf.ne.2224) return

c...Baryons that have const.quarks feel potentials.
      else if(mstc(104).eq.12) then
        if(ibar.eq.0) return           ! Exclude mesons
        if(tf.ge.t) then
          iq=mod(abs(k1)/10,10)
          if(iq.eq.0) return
        endif

c...Formed baryon + hadrons that have const.quarks feel potentials.
      else if(mstc(104).eq.13) then
        if(tf.ge.t) then
          iq=mod(abs(k1)/10,10)
          if(iq.eq.0) return
        else
          if(ibar.eq.0) return  ! Exclude formed meson
        endif

c...Formed hadrons + hadrons that have const.quarks feel potentials.
      else if(mstc(104).eq.14) then
        if(tf.ge.t) then
          iq=mod(abs(k1)/10,10)
          if(iq.eq.0) return
        endif
      endif

      jamrqpb=.false.
      end

**********************************************************************
      real*8 function gfacs(i)
c...option for the scalar coupling
      implicit none
      integer i,kf,kc,str,jamcomp
      include 'jam1.inc'
      include 'jam2.inc'

      gfacs=1.0
      if(mstc(121).eq.0) return

      if(mstc(121).eq.1) then ! all resonance
        if(abs(k(2,i)).ne.2212.and.abs(k(2,i)).ne.2112) then
          gfacs=p(5,i)/parc(25) 
        endif
      else if(mstc(121).eq.2) then ! only delta
          kf=abs(k(2,i))
          if(kf.eq.1114.or.kf.eq.2114.or.kf.eq.2214.or.kf.eq.2224)
     &          gfacs=p(5,i)/parc(25) 

      else if(mstc(121).eq.3) then
        if(abs(k(2,i)).ne.2212.and.abs(k(2,i)).ne.2112) then
          gfacs=parc(106) 
        endif

      else if(mstc(121).eq.4) then
          gfacs=parc(106) 

      endif

c....2020/10/26 strange baryons.
      if(mstc(125).ne.0) then
        kc=jamcomp(k(2,i))
c       str=kchg(kc,7)*isign(1,k(2,i))
        str=kchg(kc,7)
        if(str.ne.0) gfacs = parc(108)*abs(str)*gfacs
      endif

      end

**********************************************************************
      real*8 function vfacs(i)
c...option for the vector coupling
      implicit none
      integer i,kc,jamcomp,str
      include 'jam1.inc'
      include 'jam2.inc'

      vfacs=1.0
      if(mstc(122).eq.1) then
        if(abs(k(2,i)).ne.2212.and.abs(k(2,i)).ne.2112) then
          vfacs=parc(107)
        endif
      else if(mstc(122).eq.2) then
          vfacs=parc(107)
      endif

c....2020/10/26 strange baryons.
      if(mstc(125).ne.0) then
        kc=jamcomp(k(2,i))
        str=kchg(kc,7)
        if(str.ne.0) vfacs = parc(109)*abs(str)*vfacs
      endif

      end

**********************************************************************
      real*8 function qfacr(i)

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      real*8 qnum1
      integer i,iqcnum

      qfacr=1.0d0

      if(mstc(117).ge.1) then
       if(abs(k(2,i)).eq.2212.or.abs(k(2,i)).eq.2112) then
       else
         qfacr=2.0/3.0
       endif
      endif

      if(mstc(104).lt.12) return
      if(mstc(104).eq.15) return
      if(r(5,i).lt.pard(1)+1d-8) return  ! formed hadron

      qnum1=3.0d0
      if(k(9,i).eq.0) qnum1=2.d0
      iqcnum=mod(abs(k(1,i))/10,10)
      if(iqcnum.eq.3) iqcnum=2
      qfacr=qfacr*iqcnum/qnum1

      end


************************************************************************
      real*8 function restmass(meff,pot)
      implicit none
      include 'jam2.inc'
      real*8 meff,pot

c     if(mstc(106).ge.201) then
c       restmass=meff/(1+pot)
c     else if(mstc(106).ge.101) then
      if(mstc(108).ge.3) then
        restmass=meff-pot
      else
        restmass=sqrt(meff**2+pot**2) - pot
      endif

      end

************************************************************************
      real*8 function effmsq(i)

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      integer i

      effmsq=p(5,i)**2
      if(MF_on(i).eq.0) return

      if(mstc(109).le.1) return

c     if(mstc(106).ge.201) then
c       effmsq=p(5,i)*(1.0+pots(i))
c     else if(mstc(106).ge.31) then
      if(mstc(108).ge.3) then
        effmsq=(p(5,i)+pots(i))**2
      else
        effmsq=max(0d0,p(5,i)**2+2*p(5,i)*pots(i))
      endif

      end

************************************************************************
      subroutine setenergy(iopt)

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      integer i,iopt
      real*8 effmsq,emf

      do i=1,nv
        emf=p(5,i)**2
        if(iopt.eq.1) emf=effmsq(i)
        p(4,i)=sqrt(emf+p(1,i)**2+p(2,i)**2+p(3,i)**2)
      end do

      end

c***********************************************************************

      subroutine jamepart(difm,difv)  ! energy

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
      real*8 difm,difv,emf,pp2,em0,effmsq,potv0(0:3),vsky

      difm=0.d0   ! Difference of effective mass
      difv=0.d0   ! Difference of vector potential

c....Sigma-omega model.
      if(mstc(108).eq.2) then
        call jamepart2(difm,difv)
        return
      else if(mstc(108).eq.3) then
        call jamepart3(difm,difv)
        return
      else if(mstc(108).eq.4) then
        call jamepart4(difm,difv)
        return
      else if(mstc(108).eq.5) then
        call jamepart5(difm,difv)
        return
      else if(mstc(108).eq.0) then
        call jamepart0
        return
      endif

      do 100 i=1,nv

       if(MF_on(i)==0) then
         pots(i)=0.0d0
         potv(:,i)=0.0d0
         goto 100
       endif

       emf=p(5,i)
       potv0(:)=potv(:,i)
       potv(1,i)=0.0
       potv(2,i)=0.0
       potv(3,i)=0.0

c...Skyrme potential
         if(mstc(109).ge.1) then
c          potv(0,i)=0.0
           potv(:,i)=0.0
           pots(i) = t1*rho(i) + t3*rho(i)**pard(103)
     $         +   rhoc(i) + rhos(i) + rhoy(i) + vmoms(i)
         else
           pots(i)=0.0
           vsky=k(9,i)/3*(t1*rhoi(i) + t3*rhoi(i)**pard(103))
           potv(0,i) = vsky + rhoc(i) + rhos(i) + rhoy(i) + vmom(0,i)
           if(mstc(115).eq.2.and.rhoi(i).gt.1d-8) then
             potv(:,i)=vsky*rhoj(:,i)/rhoi(i) + vmom(:,i)
           endif
         endif

       if(mstc(112).eq.0) emf=sqrt(effmsq(i))

       pp2 = p(1,i)**2 + p(2,i)**2 + p(3,i)**2
       em0=sqrt(max(0d0,p(4,i)**2 - pp2))
       difm=difm+abs(em0-emf)
       difv=difv+sqrt((potv0(0)-potv(0,i))**2+(potv0(1)-potv(1,i))**2
     &    +(potv0(2)-potv(2,i))**2+(potv0(3)-potv(3,i))**2)
c      p(4,i) = sqrt(emf**2 + pp2)

100   continue

      mste(44)=1
      end

************************************************************************
      subroutine jamrqmde(ekin,epot,etot)
************************************************************************
c...Purpose: to calculate energy in RQMD/S
c H = sigma_i=1^N  1/(2*E_i) *[E_i^2 -vec{p}_i^2 - m_i^2 - 2m_i*V_i ]
c here V_i = t1 * <rho_i> + t3 *<rho_i>^gamma : Skyrme Potential

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer i
      real*8 ekin,epot,etot,em,pp2,p4,ee,emf2,diff,diffv
      logical jamrqpb

      if(mste(44).eq.0) then
        do i=1,nv
           MF_on(i)=1
           if(jamrqpb(i)) MF_on(i)=0
        enddo
        if(mstc(108).eq.2) then
          call jamrqmm2
        else if(mstc(108).eq.3) then
          call jamrqmm3
        else if(mstc(108).eq.4) then
          call jamrqmm4
        else if(mstc(108).eq.5) then
          call jamrqmm5
        else if(mstc(108).eq.0) then
          call jamqmm
        else
          call jamrqmm1(0)
        endif
        call jamepart(diff,diffv)
      endif

c...Local Variable
      ekin=0.0d0
      epot=0.0d0

      do 100 i=1,nv
        if(k(1,i).gt.10) goto 100

        em=p(5,i)
        pp2=p(1,i)**2+p(2,i)**2+p(3,i)**2
        if(MF_on(i)==0) then
         ekin = ekin + sqrt(pp2+em**2)
         goto 100
        endif

        if(mstc(119).eq.0) then
          pp2=(p(1,i)-potv(1,i))**2
     &      + (p(2,i)-potv(2,i))**2 + (p(3,i)-potv(3,i))**2
        endif

        p4=sqrt(pp2+em**2)
        ekin=ekin+p4

        if(mstc(109).eq.2) then
          ee=p(4,i) + potv(0,i)
          if(mstc(119).eq.0) then
            em=p(5,i)+pots(i)
            ee=sqrt(em**2 + pp2) + potv(0,i)
          endif

        else
c         if(mstc(106).ge.201) then
c           emf2 = (em*(1+pots(i)))**2
c         else if(mstc(106).ge.31) then
          if(mstc(108).ge.3) then
            emf2 = (em+pots(i))**2
          else
            emf2 = em**2 + 2*em*pots(i)
          endif
          ee=sqrt(pp2+emf2)+potv(0,i)
        endif
        epot=epot+ee-p4
100   continue
      ekin=ekin/mstd(11)
      epot=epot/mstd(11)
      etot=ekin+epot

      end

************************************************************************

      subroutine jamrqen(ind)  ! energy

c...Purpose: to calculate energy in RQMD/S
c H = sigma_i=1^N  1/(2*E_i) *[E_i^2 -vec{p}_i^2 - m_i^2 - 2m_i*V_i ]
c here V_i = t1 * <rho_i> + t3 *<rho_i>^gamma : Skyrme Potential

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      integer ind,n0,nn,i
      real*8 vfac

      call jamrqmm1(ind)

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

      do 13 i=n0,nn
       pots(i)=0.0
       potv(:,i)=0.0
       if(MF_on(i)==0) goto 13

        vfac=1.0d0
        pots(i) = vfac*(t1*rho(i) + t3*(rho(i)**pard(103)))
     $         + rhoc(i) + rhos(i) + rhoy(i) + vmoms(i)

13    continue

      end
************************************************************************
      subroutine jamrqch(dt2,ic,ind) ! coolheat

c...Purpose: Cool and Heat Projectile(1)/Target(2) Nuclei to fit 
c            the binding energy
c
c  d P_i/dt = - del<H>/del R_i - mu/b *  del <H>/del P_i 
c  d R_i/dt =   del<H>/del P_i - mu   *  del <H>/del R_i non-rel 
c  force(n,i)  = +del <H>/del R_i <= dP_i/dt   <H>is rel. 
c  forcer(n,i)  = -del <H>/del P_i <= dR_i/dt
c  sp =  mu/b > 0
c  sr =  mu   > 0

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
c     dimension force(3,mx),forcer(3,mx)
      data sp,sr/-1.d-1,-1.d-1/
      call jamrqmd(ind)
      call jamrqen(ind)

      if (ind.eq.1) then         ! for target baryon
         n0 = 1
         nn = mstd(5)
      elseif (ind.eq.2) then     ! for projectile baryon
         n0 = mstd(5)+1
         nn = mstd(2)+mstd(5)
      endif

      do i=n0,nn
        do n=1,3
c...[Modified by AO 2002/08/30
            vn=p(n,i)/p(4,i)
c forcer(n,i)+vn=   dH/dp(n,i)
c force(n,i)   = - dH/dr(n,i)
c cooling Eq.: dp/dt=-dH/dp =  fr+vn
c              dr/dt=-dH/dr = -fp
          if (ic.eq.1) then      ! cool
            p(n,i)=p(n,i) + dt2*(force(n,i)     + sp*(forcer(n,i)+vn))
            r(n,i)=r(n,i) + dt2*(forcer(n,i)+vn - sr* force(n,i)    )
          elseif (ic.eq.2) then  ! heat
            p(n,i)=p(n,i) + dt2*(force(n,i)     - sp*(forcer(n,i)+vn))
            r(n,i)=r(n,i) + dt2*(forcer(n,i)+vn + sr* force(n,i)    )
          endif 
c...]Modified by AO 2002/08/30
        enddo
        e2 = p(5,i)**2 + p(1,i)**2 + p(2,i)**2 + p(3,i)**2
        p(4,i) = sqrt( e2 +2.d0*p(5,i)*pots(i)) !m
c       p(4,i) = sqrt( e2 )
      enddo

      return 
      end

************************************************************************
      subroutine jamrqpt(ind) ! trans
************************************************************************
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'
      parameter (de=0.1d-3)
      data dt2/0.05d0/
      data sau/-7.9d-3/         ! binding energy per nucleon for 197Au
      if (ind.eq.1) then         ! for target baryon
         n0 = 1
         nn = mstd(5)
         ia = mstd(5)
         iz = mstd(6)
      elseif (ind.eq.2) then     ! for projectile baryon
         n0 = mstd(5)+1
         nn = mstd(2)+mstd(5)
         ia = mstd(2)
         iz = mstd(3)
      endif
        call findbeng(ia,iz,beng)
        sau = -beng/dble(nn-n0+1)
      do 1000 it=1,10000
        call jamrqen(ind)
        esky = 0.d0
        emom = 0.d0
        eqmas= 0.d0
        ekin = 0.d0
        do i = n0,nn
          esky  = esky  + pots(i)
c         esky  = esky  + vsky(i)
          emom  = emom  + vmom(0,i) !m
          eqmas = eqmas + p(5,i)
          ekin  = ekin  + sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+p(5,i)**2)
        enddo
        ekin= ekin - eqmas
        if(mstd(101).eq.1) then
          eng = esky + emom + ekin !m
        else
          eng = esky + ekin 
        endif
        esky= esky/dble(nn-n0+1)
        emom= emom/dble(nn-n0+1)
        ekin= ekin/dble(nn-n0+1)
        sine= eng /dble(nn-n0+1)  ! Binding energy (AGeV)
        if(mstc(8).ge.2.or.(mstc(8).ge.1.and.mod(it,10).eq.0)) then
          write(*,820) it, sau-de , sine , sau+de , eng, ekin ! check
        endif
820     format (i5,1x,6(f13.7,1x))

        if ((sine.ge.sau-de) .and. (sine.le.sau+de))  then
          return
        else
c          dt2=abs(sau-sine)*10.d0
c         dt2=abs(sau-sine)
          if (sine.gt.sau) call jamrqch(dt2,1,ind) ! cool (underbinding)
          if (sine.lt.sau) call jamrqch(dt2,2,ind) ! heat (overbinding)
        endif
          
1000  continue

      return
      end 

************************************************************************
        subroutine findbeng(ia,iz,beng) ! findbeng
************************************************************************
      implicit double precision(a-h, o-z)
       include 'jam1.inc'
       character symbol*5
       open(410,file='beng.dat',status='old',err=1000)
c     print *, ia,iz
10     read(410,*,end=1000,err=1000) kz,kn,beng,symbol
c     print *, kz,kn,beng,symbol
       if (kn+kz.eq.ia.and.kz.eq.iz) then
        write(*,*)"B.E of",ia,"-",symbol,"-",iz,"=",beng,"MeV"
        beng=beng/1000.d0
        close(410)
        return
       endif
       goto 10

1000   continue
       beng=beld(ia,iz)/1000.d0
c20    continue
c      write(*,*)"(findbeng) error: cannot find B.E."
c30     continue
c       write(*,*)"(findbeng) error: read error occur in subroutine"
      end
************************************************************************
        function beld(ia,iz)
************************************************************************
        implicit real*8(a-h,o-z)
c a1              = -15.4602          #       +/- 0.105768
c a2              = 18.4107           #       +/- 0.256056
        parameter(av=-15.68d0,as=18.56d0,ac=0.717d0,ai=28.1d0)
        beld=0.0d0
        if(ia.le.1) return
        aa=dble(ia)
        aa3=aa**(1.0d0/3.0d0)
        beld=av*ia+as*aa3*aa3+ac*iz**2/aa3+ai*(2*iz-ia)**2/ia
        end

************************************************************************
        integer function zbstable(ia)
************************************************************************
      implicit none
      integer ia
      real*8 ac,ai
      parameter(ac=0.717d0,ai=28.1d0)
      zbstable=0.5d0*ia/(1d0+ac/ai*0.25*dble(ia)**(2.0d0/3.0d0))
      end

************************************************************************
      subroutine qmdpotparam
      implicit none
      real*8 hc,rho0,conv
      integer i
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      include 'jam4.inc'

c...Width parameter.
      pard(104)= parc(15)            ! The width of Gaussian L (fm^2)

c...Option for the initial nucleon position.
      if(mstc(6).ge.100.and.mstc(95).eq.1) then 
        ! WS for L=1fm^2 for Au and Pb nculeus.
        mstc(98)=4
        parc(13)=0.9
        parc(14)=0.9
        parc(20)=0.03
        print *,'Initial condition for QMD mode'
        if((mstd(2).lt.100.and.mstd(2).ge.3)
     &        .or.  (mstd(5).lt.100.and.mstd(5).ge.3)) then
          write(6,*)'For now, this initialization mstc(95)=',mstc(95),
     &          'cannot be used for small collision system'
        endif
      endif

c...Initialize potentials.
      do i=1,mxv
        pots(i)=0.0d0
        potv(:,i)=0.0d0
        MF_on(i)=0
      end do

      if(mstc(108).eq.2) then ! RQMD/S old
        call rqmdpotparam2
      else if(mstc(108).eq.5) then ! sigma-omega interaction
        call rqmdpotparam5
        return
      else if(mstc(108).eq.3) then ! Skyrme scalar
        call rqmdpotparam3
      else if(mstc(108).eq.4) then ! Skyrme vector
        call rqmdpotparam4
      else if(mstc(108).eq.1) then ! RQMD/S original (obsolete)
        call rqmdpotparam1
      else if(mstc(108).eq.0) then ! non-relativistic QMD
        call rqmdpotparam1
      else
        print *,'qmdpotparm wrong mstc(108) = ',mstc(108)
        stop
      endif

      if(mstc(105).eq.0.and.mstc(119).eq.0) then
        print *,'inconsistent mstc(105)=',mstc(105),
     &        'mstc(119)=',mstc(119)
        stop
      endif

      hc=paru(3)
      rho0=parc(21)
      conv=hc*hc*hc

      t1=pard(101)/2.d0/rho0
      t3=pard(102)/(pard(103)+1.d0)/(rho0**pard(103))
      t3f=pard(103)*t3
      t2=0.0
      t4=0.0
      t4f=0.0
      pmu1=pard(105)**2
      pmu2=pard(106)**2
      vex1=pard(107)/(2*rho0)
      vex2=pard(108)/(2*rho0)

      pard(109)= 1.0/(4.d0*pard(104))
      pard(110)=1.0d0/(4.d0*paru(1)*pard(104))**1.5d0 ! 1/[(4*pi*L)^3/2]

      if(mstc(6).ge.100) then
      if(mstc(108).eq.3) then
        print *,'RQMD Skyrme scalar mode'
      else if(mstc(108).eq.4) then
        print *,'RQMD Skyrme vector mode'
      endif
      print *,'t1=',t1,'t3=',t3,'t3f=',t3f
      print *,'pard(109)=',pard(109),pard(110)
      endif

c....Vector potential on.
c     if(mstc(111).ge.1) then
c....This value is for mstc(106)=201 and pard(101)=0,pard(102)=-0.3, gamma=0.3
c       pard(111)= -0.459272   ! alpha
c       pard(112)=  0.398717   ! beta
c       pard(113)=  1.13482    ! gamma

c....K=200 for a=-0.8 g=0.6
c       pard(111)= -0.555522
c       pard(112)=  0.489625
c       pard(113)=  1.11408

c....K=380 for a=-0.8 g=0.6
c       pard(111)= -0.152978
c       pard(112)=  0.0870809
c       pard(113)=  1.87108
c     endif

      end

