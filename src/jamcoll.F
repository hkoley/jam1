c#define _CUDA 1
***********************************************************************
c***********************************************************************
c                                                                      *
c        PART 2: collision part                                        *
c                                                                      *
c   List of Subprograms in rough order of relevance with main purpose  *
c      (s = subroutine, f = function, b = block data, e = entry)       *
c                                                                      *
c  s  jamcoll  to simulate collision and decay                         *
c  s  jamhit1  to determine whether the particle can collide or not    *
c  f  jamhit   to determine whether collision occurs                   *
c  s  jamclist to load initial values into collision predictor arrays  *
c  s  jamcfind to search decay and collision arrays to find next event *
c  s  jamcupda to update collision array                               *
c  s  jamupdcB to update collision array for particles feeling MF      *
c  s  jamtrspt to transport all particles ahead in time                *
c  f  jamcltyp to give collision type                                  *
c  s  jamscatt to perform two-body collisions                          *
c  s  jamchanl to select final scattering channel                      *
c  s  jamabsrb to treat annihilation scattering                        *
c  s  jamabbann to treat annihilation for antibaryon-baryon collision  *
c  s  jamivide to subdivide hadron into quarks                        *
c  s  jamscat2 to do kinematics for two-body final scattering          *
c  s  jamscatd to handle deuteron pion final scattering                *
c  s  jamscatd2to handle double pionic fusion scattering               *
c  s  jamn2dec to perform two-body decay                               *
c  s  jamscatm to handle multi-particle final scattering               *
c  s  jamscat3 to do kinematics for three-body final scattering        *
c  s  jamscatp to scatter the particle according to phase space        *
c  s  jamkupda to update particle status code k() after collision      *
c  s  jamvupda to update vertex of the particles                       *
c  s  jamtupda to update time of particle                              *
c  s  jamsave  to save or reset particle information                   *
c  s  jamangel to calculate elastic angular distribution               *
c  f  jamslope to give slope parameter in two-body collisions          *
c  s  jamangin to calculate inelastic angular distribution             *
c  s  jamangrr to generate Angular distribution of RR                  *
c  s  jamangdn to determine inelastic scattering angle for nn->nd      *
c  s  jamedit  to remove unwanted entries from record                  *
c  s  jamexch  to exchange position of event record vectors            *
c  s  jamzero  to zero the vectors r(),p(),v(),k(),kq(),vq()           *
c  s  jamprcl  to print collision information                          *
c  s  jamcheck to check possible errors after collision or decay       *
c  s  jampauli to account for Pauli-blocking in collision              *
c  f  jamisjet to determine number of jet                              *
c  s  jamcmom  to construct energy-momenta of proj. and targ.          *
c  s  jamdeut  to determine deuteron nuclear cluster by coalescence    *
c  s  jamglaub to do AA collision according to the Glauber theory      *
c  s  jamfpath to hA collision according to the mean free path argument*
c  s  jaminil  to initialize some values for leading particle cascade  *
c  s  jamfgas  to calculate Fermi momentum using Fermi gas model       *
c  s  jamdnsf  to give parameters for nuclear Fermi density            *
c                                                                      *
c***********************************************************************
c                                                                      *
c***********************************************************************
c***********************************************************************

      subroutine jamcoll

c...Purpose: to simulate collision and decay.
c...Modify for box.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      character rtype*3
      dimension indd(100)
      real*8 u(4),pk1(5),pk2(5)

c...Set max. time for cascading.
c     pard(1)=(mstd(23)-1)*parc(2)
c     pare(1)=mstd(23)*parc(2)
      ftime=pare(1)
      gtime=pard(1)
crqmd2014/11/12 
      dtRQMD=parc(16)
      if(dtRQMD.eq.0) dtRQMD=0.5d0

c...Initialize time dependent analysis.
      if(mstc(161).ne.0.and.mstc(3).eq.1.and.mstd(23).eq.0)
     $ call jamanat(1)

crqmd...Evaluate hadron density and force in RQMD/S.
      if(mstc(6).ge.100) then
#ifdef _CUDA
           if(mstc(129).ne.0) then
              NrqmdCount=icheckMF(0,pard(1))
              call jamrqmdcuda(0)
              if(NrqmdCount.gt.0.and.mstc(103).gt.0) call jamdtfree
           else
              call jamrqmd(0)
           endif
#else
        call jamrqmd(0)
#endif
c...Evaluate hadron density and force in BUU.
      else if(mstc(6).ge.2) then
        call jambuud
        call jambuuf
      endif

      nsimul=mstc(5)
      if(mstc(46).eq.2) nsimul=1

      do isimul=1,nsimul

       mstd(8)=isimul
       pare(1)=ftime
       pard(1)=gtime

crqmd
       istepRQMD=0
       ftimeRQMD=gtime+dtRQMD

c....No collision and no mean field.
      if(mstc(6).eq.1) then 
        goto 1999
      else if(mstc(6).ge.2.and.mstc(6).le.10) then
        goto 2000
      endif

      if(mstc(7).eq.0) goto 1999 ! No collision.

c...Load collision arrays.
10000 continue

c...Collision list will not be update for hydro+cascade simulation.
c     if(mstd(23).eq.0.or.mstc(140).eq.0) then
c     if(mstc(140).eq.0) then
c     if(mstc(140).eq.0.or.mstc(6).gt.0) then
      if(mstd(20).eq.0) then

c...2019/4/2
      if(mstc(18).eq.0.and.mstd(30).ne.0) call jamedit

#ifdef _CUDA
      if(mstc(129).ne.0) then
        call jamclistcuda
      else
        call jamclist
      endif
#else
      call jamclist
#endif
      endif

c...Search for next collision event.
      ifnd0=0   ! rqmd
      ifnd=0
 1000 continue
      call jamcfind(ind,rtype)

crqmd...[ AO: 050727
      if(mstc(6).ge.100.and.mstc(103).eq.2) then
        tnext=pard(1)
        if(rtype.eq.'col') then
          tnext=coll(1,ind)
        elseif(rtype.eq.'dec') then
          tnext=v(5,ind)
        endif
        if(tnext.gt.ftimeRQMD) then
          istepRQMD=istepRQMD+1
          if(mstc(8).ge.10) then
            write(mstc(38),'(a,3(1x,i4),1x,a3,3(1x,f7.2))')
     &        'step,ment,ifnd,typ,T(col,cur,stp)=',
     &       istepRQMD,mentry,ifnd-ifnd0,rtype,tnext,pard(1),ftimeRQMD
          endif
c            write(*,'(a,3(1x,i4),1x,a3,3(1x,f7.2))')
c     &	      'step,ment,ifnd,typ,T(col,cur,stp)=',
c     &	      istepRQMD,mentry,ifnd-ifnd0,rtype,tnext,pard(1),ftimeRQMD
          ifnd0=ifnd
          pard(1)=ftimeRQMD
          call jamtrspt(0,pard(1))
#ifdef _CUDA
           if(mstc(129).ne.0) then
              NrqmdCount=icheckMF(0,pard(1))
              call jamrqmdcuda(0)
              if(NrqmdCount.gt.0.and.mstc(103).gt.0) call jamdtfree
           else
              call jamrqmd(0)
           endif
#else
          call jamrqmd(0)
#endif
          ftimeRQMD=ftimeRQMD+dtRQMD
          if(mstc(161).ne.0.and.mstc(3).eq.1) call jamanat(2)
          if(mod(istepRQMD,mstc(105)).eq.0) then
           call jamupdcB
          endif
          goto 1000
        endif
      endif
c...] AO: 050727

c...Finish collisions.
      if(rtype.eq.'end') goto 1999
      ifnd=ifnd+1
      if(ifnd.ge.mxcoll-30) then
        if(ifnd.le.mxcoll) then
          write(mstc(38),*)'ind ctime=',ind,coll(1,ind),' ',rtype
          j1=icoll(1,ind)
          write(mstc(38),*)'i1 i2',icoll(1,ind),icoll(2,ind),
     &        'kf1=',k(2,j1),' m= ',p(5,j1)
        else
         write(check(1),'(''ifnd='',i10)')ifnd
         goto 1999
c        call jamerrm(30,1,'infinit loop ? (jamcoll) ifnd')
        endif
      endif

c...Perform time dependent analysis. 2011/12/9
      if(mstc(161).ne.0.and.mstc(3).eq.1) then
        if (rtype .eq. 'col') then
          pard(1)=coll(1,ind)    ! Collision time.
        else if (rtype .eq. 'dec') then
          pard(1)=v(5,ind)
        else
          print *,'jamcoll unknown rtype ',rtype
          stop
        endif
        call jamanat(2)
      endif


c       tnext=pard(1)
c       if(rtype.eq.'col') then
c         tnext=coll(1,ind)
c         t1=coll(2,ind)
c         t2=coll(3,ind)
c         i1=icoll(1,ind)
c         i2=icoll(2,ind)
c       elseif(rtype.eq.'dec') then
c         tnext=v(5,ind)
c         t1=0.0
c         t2=0.0
c         i1=ind
c         i2=0
c       endif
c     print *,'rypte=',rtype,' t=',pard(1),pare(1),mentry
c     print *,'i1=',i1,' i2=',i2
c     print *,'tcol=',tnext,t1,t2


      mste(1)=0
c...Collision will occur.
c===================================================================
      if (rtype .eq. 'col') then
c===================================================================

c...Set time of collision or decay.
        i1=icoll(1,ind)        ! particle 1 line number 
        i2=icoll(2,ind)        ! particle 2 line number
        mste(6)=icoll(3,ind)   ! box: by maru.
        pard(1)=coll(1,ind)    ! Collision time.
        pare(8)=coll(2,ind)    ! collision time
        pare(9)=coll(3,ind)    ! collision time
        pare(4)=coll(4,ind)    ! sig
        pare(5)=coll(5,ind)    ! sigel
        pare(6)=coll(6,ind)    ! impact parameter sq
        t1=pare(8)
        t2=pare(9)

c       if(pard(1).gt.pare(1)+parc(2)) then
        if(pard(1).gt.pare(10)) then
          write(mstc(38),*)'??? time>timmax',pard(1),pare(10),rtype
          stop
          goto 1999
        end if

c     call checkCDA(i1,i2,coll(2,ind),coll(3,ind))

      if(mstc(8).ge.10) then
      write(mstc(38),*)'bef jamtrspt r4=',i1,r(4,i1),v(4,i1),coll(2,ind)
     &  ,k(1,i1),k(8,i1),mstd(8)
      write(mstc(38),*)'bef jamtrspt r4=',i2,v(4,i2),v(4,i2),coll(3,ind)
     &       ,k(1,i2),k(8,i2),mstd(8)
      endif

c...Save closest points in the collision.
      pk1(:)=p(:,i1)
      pk2(:)=p(:,i2)
      em1sq=pk1(4)**2-pk1(1)**2-pk1(2)**2-pk1(3)**2
      em2sq=pk2(4)**2-pk2(1)**2-pk2(2)**2-pk2(3)**2

c...2020/2/1
      pk1(5)=sqrt(em1sq)
      pk2(5)=sqrt(em2sq)
      if(mstd(91).ge.1) then
        if(MF_on(i1).eq.1) then
          do n=1,3
          pk1(n)=pk1(n)-potv(n,i1)
          end do
          pk1(4)=sqrt(em1sq+pk1(1)**2+pk1(2)**2+pk1(3)**2)
        endif
        if(MF_on(i2).eq.1) then
          do n=1,3
          pk2(n)=pk2(n)-potv(n,i2)
          end do
          pk2(4)=sqrt(em2sq+pk2(1)**2+pk2(2)**2+pk2(3)**2)
        endif

      endif

      dt1=t1-r(4,i1)
      pare(43)=r(1,i1)+dt1*pk1(1)/pk1(4)
      pare(44)=r(2,i1)+dt1*pk1(2)/pk1(4)
      pare(45)=r(3,i1)+dt1*pk1(3)/pk1(4)

      dt2=t2-r(4,i2)
      pare(46)=r(1,i2)+dt2*pk2(1)/pk2(4)
      pare(47)=r(2,i2)+dt2*pk2(2)/pk2(4)
      pare(48)=r(3,i2)+dt2*pk2(3)/pk2(4)

c     if(dt1.lt.0.0 .or.dt2.lt.0.0) then
c       print *,'dt1=',dt1,t1,r(4,i1),r(5,i1)
c       print *,'dt2=',dt2,t2,r(4,i2),r(5,i2)
c       read(5,*)
c     endif

      gam1=pk1(4)/pk1(5)
      gam2=pk2(4)/pk2(5)
      pare(49)=(t1-v(4,i1))/gam1
      pare(50)=(t2-v(4,i2))/gam2
      pare(56)=gam1
      pare(57)=gam2

c...Transport to time of collision.
        if(mstc(52).le.11) then
c         if(mstc(52).eq.6) then
c           t1=pard(1)
c           t2=pard(1)
c         endif
          call jamtrspt(i1,t1)
          call jamtrspt(i2,t2)

        else
          call jamtrspt(0,pard(1))
        endif

      if(mstc(8).ge.10) then
        write(mstc(38),*)'after jamtrspt r4=',i1,r(4,i1),v(4,i1),t1
        write(mstc(38),*)'after jamtrspt r4=',i2,r(4,i2),v(4,i2),t2
        if(t1 .gt.r(4,i1)) then
          print *,'transport? ',i1,t1,r(4,i1)
        endif
      endif

c...Save particle information.
        call jamsave(1,1,i1)
        call jamsave(1,2,i2)

      if(mstc(8).ge.10) then
        write(mstc(38),*)'after jamsave r4=',i1,r(4,i1),v(4,i1)
        write(mstc(38),*)'after jamsave r4=',i2,r(4,i2),v(4,i2)
      endif


c...qmd:Save initial energy of two particles.
cq      if(mstc().ge.1) then
cq        call epotall(epot0,epotpa)
cq        pare(11)=p(4,i1)+p(4,i2)+epot0
cq      endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxc
c       print *,'before scatt',k(2,i1),k(2,i2),t1,t2
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

c...Scatter two-particles.
        call jamscatt

        mste(49)=kcp(1,3)
        mste(50)=kcp(1,4)
c       pare(49)=pcp(5,3)
c       pare(50)=pcp(5,4)

c...Information.
        if(mstc(8).ge.10) then
          write(mstc(38),*)'after jamscatt channel',mste(1)
     $    ,i1,i2,k(2,i1),k(2,i2),p(5,i1),p(5,i2)
          write(mstc(38),*)'r4=',r(4,i1),(v(j,i1),j=1,4)
          write(mstc(38),*)'v4=',r(4,i2),(v(j,i2),j=1,4)
        endif

c...Count collision.
        call jamanacl(2)

c...User defined routine: normally dummy.
        call jamanaus(indd,nadd)


c...Collision forbidden.
        if(mste(1).le.0) then
          call jamsave(2,1,mste(21))
          call jamsave(2,2,mste(23))
c...2017/5/6
c         do j=1,4
c           v(j,mste(21))=r(j,mste(21))
c           v(j,mste(23))=r(j,mste(23))
c         end do
          if(mstc(8).ge.10) then
            write(mstc(38),*)'after jamsave kf1 kf2=',k(2,i1),k(2,i2)
            write(mstc(38),*)'after jamsave r4 v4 i1=',r(4,i1),v(4,i2)
            write(mstc(38),*)'after jamsave r4 v4 i2=',r(4,i1),v(4,i2)
            write(mstc(38),*)'v1=',(v(j,i1),j=1,5)
            write(mstc(38),*)'v2=',(v(j,i2),j=1,5)
            write(mstc(38),*)'i1 i2',i1,i2
            write(mstc(38),*)'i3 i4',mste(21),mste(23)
          endif
          icoll(1,ind)=-1
c         call jamcupda(mste(21),mste(23),mstc(25),mste(27),1)
cqmd      if(mstc().ge.1) call caldis2(mste(21),mste(23))
          goto 1000
        endif

c...Check absorbed particles.
        if(k(1,i1).ge.11.or.k(1,i2).ge.11) then
          mstd(30)=mstd(30)+1
        end if

cTABARA
c       call ttchk(indd,0)
c....Print collision information.
        if(mstc(8).ge.2) call jamprcl(indd,0)

c...Update the collision array.
        call jamcupda(mste(25),mste(27),mste(25),mste(27),1)
        if(mste(29).ge.1) call jamcupda(mste(29),-1,0,0,1)
        if(mste(31).ge.1) call jamcupda(mste(31),-1,0,0,1)

c...String decay.
cc      if(mstc(76).le.1) then
cc        ind=mste(21)
cc        if(k(2,ind).eq.92) then
cc        endif
cc        ind=mste(23)
cc        if(k(2,ind).eq.92) then
cc        endif
cc      endif

c...Decay events
c===================================================================
      else if (rtype .eq. 'dec') then
c===================================================================

c...Decay time.
        pard(1)=v(5,ind)

c...Transport to time of decay.
         if(mstc(52).le.11) then
           call jamtrspt(ind,pard(1))
         else
           call jamtrspt(0,pard(1))
         endif

c...Save particle information.
        call jamsave(1,1,ind)

c...Set collision type.
        mste(2)=-1
        mste(7)=0 ! Resonance decay
c...String fragmentation.
        if(mod(abs(k(1,ind)),10).eq.3.or.k(1,ind).eq.4) mste(7)=1

c...qmd:Save initital energy (resonance decay only)
cq      if(mstc().ge.1.and.k(1,ind).eq.2) then
cq        call epotall(epot,epotpa)
cq        pare(11) = p(4,ind) + epot
cq        pare(11)=p(4,ind)
cq      endif

c...Convert this particle into fluid after formation time.
      if(k(7,ind).eq.-1200.and.k(2,ind).ne.92) then
        call jam2fluid_aformation(ind)
        goto 1000
      endif

c...Compute density before decay for hydro simulation.
      isfluid=0
      if(mstc(140).ge.1.and.mste(41).ge.0.and.mstc(141).eq.3) then
c     if(mstc(142).gt.10.and.mste(7).eq.1) then
         rx=r(1,ind)
         ry=r(2,ind)
         rz=r(3,ind)
        if(k(1,ind).eq.4) then
          rx=0.0d0
          ry=0.0d0
          rz=0.0d0
          mpa=k(11,ind)-k(10,ind)+1
          do i=k(10,ind),k(11,ind)
            rx=rx+r(1,i)/mpa
            ry=ry+r(2,i)/mpa
            rz=rz+r(3,i)/mpa
          end do
        endif
        isfluid=isinside(rx,ry,rz,ix,iy,iz)
        if(isfluid.eq.1) then
          call jamdens(icon,ind,cc,rhop,rhobar,einv,u,1)
          eden=eden_fluid(ix,iy,iz)+einv
          bden=bden_fluid(ix,iy,iz)+rhobar
          if(icon.ne.0.or.eden.lt.edminfluid(bden)) isfluid=0
        endif
      endif
c     endif

        call jamdec(ind,indd,nadd,icon)

        if(icon.ne.0) then
          call jamcupda(ind,-1,ind,0,1)
          goto 1000
        endif


cTABARA
c       call ttchk(indd,nadd)
c...Print information after decay.
        if(mstc(8).ge.2) call jamprcl(indd,nadd)

c...User defined routine: normally dummy.
        call jamanaus(indd,nadd)

c...Update the collision array
c       if(mstc(140).eq.0.or.mstc(141).ne.3) then
c       if(mste(41).ge.0.and.mstc(141).eq.3) then
        if(isfluid.eq.1) then
            call jam2fluid_decay(einv,efluid,indd,nadd,ind)
        else

        kf0=kcp(2,1)
        ist=0    ! rqmd
        ist0=0
        if(kf0.ne.92) then
          call attflv2(kf0,kfl1,kfl2,kfl3)
          if(abs(kfl1).eq.3) ist0=ist0+1
          if(abs(kfl2).eq.3) ist0=ist0+1
          if(abs(kfl3).eq.3) ist0=ist0+1
        endif
        do i=1,nadd
         j=indd(i)
         call jamcupda(indd(i),-1,ind,0,1)
         call attflv2(k(2,j),kfl1,kfl2,kfl3)
         ist=0  ! 2010/8/24
         if(abs(kfl1).eq.3) ist=ist+1
         if(abs(kfl2).eq.3) ist=ist+1
         if(abs(kfl3).eq.3) ist=ist+1
         if(ist.gt.ist0) then
           if(kf0.eq.92) then
             mstd(63)=mstd(63)+1
           else
             mstd(64)=mstd(64)+1
           endif
         endif
        end do

c...hydro: convert produced mesons to fluid elements.
c       else
c         if(mste(41).ge.0.and.mstc(141).eq.3) then
c           call jam2fluid_decay(indd,nadd,ind)
c         endif
        endif


 
      else
c===================================================================
      endif
c===================================================================

c...Perform time dependent analysis.
 2000 continue
c     if(mstc(161).ne.0.and.mstc(3).eq.1) call jamanat(2)

c....Save max. number used in vector.
      mstd(81)=max(mstd(81),nv)
      if(nv.ge.mxv*0.8.and.mstd(30).ne.0) then
          call jamedit
          goto 10000
      endif

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c     idbx=0
c     if(idbx.eq.1) then
c     do i=1,nv
c       if(k(1,i).gt.10) cycle
c       kc1=jamcomp(k(2,i))
c       id1=kchg(kc1,5)
c       pp=p(1,i)**2+p(2,i)**2+p(3,i)**2
c       emf=sqrt(p(4,i)**2-pp)
c       if(abs(k(9,i)).eq.0.and.emf.le.0.1) then
c       print *,'rtype=',rtype,' mste',mste(1),' mass? ',i,k(2,i),p(5,i)
c       print *,'i1 i2',i1,i2,emf
c       emf1=sqrt(pcp(4,1)**2-pcp(1,1)**2-pcp(2,1)**2-pcp(3,1)**2)
c       emf2=sqrt(pcp(4,2)**2-pcp(1,2)**2-pcp(2,2)**2-pcp(3,2)**2)
c       print *,kcp(2,1),pcp(5,1),emf1
c       print *,kcp(2,2),pcp(5,2),emf2
c       endif
c     end do

c     do i=1,mentry
c        l1=icoll(1,i)
c        l2=icoll(2,i)
c        dt1=coll(2,i)-r(4,l1)
c        dt2=coll(3,i)-r(4,l2)
c        if(dt1.lt.0d0) then
c             print *,'after i1 i2',rtype,i1,i2,j1,j2,ifd
c             print *,'mste ',mste(1)
c             print *,'i1 i2 i3 i4',mste(21),mste(23),mste(25),mste(27)
c             print *,'coload',l1,dt1,k(1,l1),k(2,l1)
c             print *,'coll',coll(1,i),r(4,l1),k(4,l1)
c        endif
c        if(dt2.lt.0d0) then
c             print *,'after decay i2',i1,i2,j1,j2,ifd
c             print *,'mste ',mste(1)
c             print *,'i1 i2 i3 i4',mste(21),mste(23),mste(25),mste(27)
c             print *,'coload',l2,dt2,k(1,l2),k(2,l2)
c             print *,'coll',coll(2,i),r(4,l2),k(4,l2)
c        endif
c     end do
c     endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

c...Go back to the next collision.
      goto 1000
c===================================================================
 1999 continue                              ! collision loop closed
c===================================================================

c2019/4/2 move to the upper: before "call jamclist"
c...Delete absorbed particles.
c...2014?/1/27 does not want jamedit in hydro  since after that,
c...collision list has to be updated.
c     if(mstc(140).eq.0) then
c     if(mstc(18).eq.0.and.mstd(30).ne.0) call jamedit
c     endif

c...Final time dependent analysis.
      if(mstc(161).ne.0.and.mstc(3).eq.1) call jamanat(3)

c...Transport to end of time-slice if time step is not one.
      if (mstc(3).gt.1.and.pard(1).lt.pare(10)) then
        call jamtrspt(0,pare(10))
        pard(1)=pare(10)
      endif

c...End loop over parallel ensemble.
      end do

c     do i=1,nv
c      if(k(2,i).eq.92) then
c        print *,'end of collision string?',k(2,i),k(1,i),
c    &   v(5,i),pard(1),pare(1)
c      endif
c     end do

      end

crqmd
c***********************************************************************

      subroutine jamhit1(i1,ICC)

c...Purpose: to determine whether collision occurs by closest distance
c...approach.
c...Modify for box.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      ICC=0

c...Different Simulation
      if(mstc(46).ne.2.and.k(8,i1).ne.mstd(8)) return

c...This particle is already ahead in time
      nop=1
      if(r(4,i1).gt.pare(1)) return

c...This particle was already dead.
      nop=2
      kp1=k(1,i1)
      if(kp1.ge.11) return

c...Skip if switch of jet system having no life time on.
      nop=3
      if(mstc(76).le.1.and.mod(abs(kp1),10).eq.3) return

c....Skip No Flavor/lepton, gamma, etc.
      nop=4
      nop=5
      kf1=k(2,i1)
      kf1a=abs(kf1)
      if(kf1.eq.0) return
      if((kf1a.gt.10.and.kf1a.le.100).and.kf1.ne.21) return

c...Option for only hadron-hadron collisions.
      nop=7
      kfl1c=mod(kf1a/10,10)
      if(abs(kf1).le.10) then
        iq1=1  ! Quark
      else if(kfl1c.eq.0) then
        iq1=2  ! Diquark
      else if(kf1.eq.21) then
        iq1=3  ! Gluon
      else
        iq1=0  ! Hadron
      endif
      if(mstc(51).eq.2.and.iq1.ne.0) return

c...Partons
      nop=8
      if(iq1.ne.0) then
        if(i1.ne.k(10,i1).and.i1.ne.k(11,i1)) return
        if(k(2,k(10,i1)).eq.21.or.k(2,k(11,i1)).eq.21) return
      endif

c...This particle can collide.
      ICC=1
c
      end

c***********************************************************************

      function jamhit(i1,i2,icell,ctime,tcol1,tcol2,bsq,sig,sigel)

c...Purpose: to determine whether collision occurs by closest distance
c...approach.
c...Modify for box.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      logical jamhit
      parameter(mxchan=30)
      dimension sigin(mxchan)
      dimension pv(2,5),kfv(2),k9v(2),kfq(2,2)
      real*8 pc1(5),pc2(5)

      real*8 pcm(5),prl(5),q(5),ppcm,pcmx,qx,bbb

      ctime=0.0d0
      tcol1=0.0d0
      tcol2=0.0d0
      bsq=0.0d0
      sig=0.0d0
      sigel=0.0d0
      nop=1

c....Collision are allowed between same test particles.
c...This is already checked.
c     if(mstc(46).eq.1.and.(k(8,i1).ne.k(8,i2))) goto 100

c...Option for only BB collisions.
c     if(mstc(51).eq.0.and.(k(9,i1).ne.3.or.k(9,i2).ne.3)) goto 100

c...Option for only BB MB collisions, i.e. no mm collisions.
c     if(mstc(51).eq.1.and.(k(9,i1).eq.0.and.k(9,i2).eq.0)) goto 100

      nop=2
      kp1=k(1,i1)
      kt1=k(1,i2)
c...This particle was already dead.
      if(kp1.ge.11) go to 100
      if(kt1.ge.11) go to 100

      nop=3
c...Skip if switch of jet system having no life time on.
      if(mstc(76).le.1.and.
     $   (mod(abs(kp1),10).eq.3.or.mod(abs(kt1),10).eq.3) ) goto 100

      nop=4
      kf1=k(2,i1)
      kf1a=abs(kf1)
      if(kf1.eq.0) goto 100

      nop=5
c....Skip if lepton, gamma, etc.
      if((kf1a.gt.10.and.kf1a.le.100).and.kf1.ne.21) goto 100
      kfl1c=mod(kf1a/10,10)
c...Quark
      if(abs(kf1).le.10) then
        iq1=1
c...Diquark
      else if(kfl1c.eq.0) then
        iq1=2
c...Gluon
      else if(kf1.eq.21) then
        iq1=3
c...Hadron
      else
        iq1=0
      endif

      nop=6
      kf2=k(2,i2)
      if(kf2.eq.0) goto 100
      kf2a=abs(kf2)
      if((kf2a.gt.10.and.kf2a.le.100).and.kf2.ne.21) goto 100

      kfl2c=mod(kf2a/10,10)
      if(kf2.eq.21) then
        iq2=3
      else if(abs(kf2).le.10) then
        iq2=1
      else if(kfl2c.eq.0) then
        iq2=2
      else
        iq2=0
      endif

      nop=7
c...Option for only hadron-hadron collisions.
      if(mstc(51).eq.2.and.(iq1.ne.0.or.iq2.ne.0)) goto 100

c...Option for only h-h and h-p
cc    if(mstc(51).eq.3.and.(iq1.ne.0.and.iq2.ne.0)) goto 100


cc    if(mstc(51).eq.3) then
        nop=8
        if(iq1.ne.0) then
          if(i1.ne.k(10,i1).and.i1.ne.k(11,i1)) goto 100
          if(k(2,k(10,i1)).eq.21.or.k(2,k(11,i1)).eq.21) goto 100
        endif
        nop=9
        if(iq2.ne.0) then
          if(i2.ne.k(10,i2).and.i2.ne.k(11,i2)) goto 100
          if(k(2,k(10,i2)).eq.21.or.k(2,k(11,i2)).eq.21) goto 100
        endif
cc    endif

c...Avoid collision of same string system.
      nop=10
      if(iq1.ne.0.and.iq2.ne.0) then
        if(i1.ge.k(10,i2).and.i1.le.k(11,i2)) goto 100
        if(i2.ge.k(10,i1).and.i2.le.k(11,i1)) goto 100
      endif

      nop=11
c...Avoid first collisions within the same nucleus
      if(k(7,i1)*k(7,i2).eq.mstc(54)) go to 100

      nop=12
c...Avoid second collisions for the same pairs
      if((k(5,i1).eq.k(5,i2)).and.(k(5,i2).ne.-1)) go to 100

      kf1=k(2,i1)
      kf2=k(2,i2)
      ibar1=k(9,i1)
      ibar2=k(9,i2)
c     srt=sqrt(max(0d0,(p(4,i1)+p(4,i2))**2
c    $ -(p(1,i1)+p(1,i2))**2-(p(2,i1)+p(2,i2))**2-(p(3,i1)+p(3,i2))**2))

c...Check if the collision between particles which is supposed to be
c...fluid after formation time can happen or not.
      if(k(7,i1).eq.-1200.or.k(7,i2).eq.-1200) then
        goto 100 ! 2019/8/29
c       if(srt.lt.parc(145)) goto 100 ???
      endif

      call jamcmom(i1,i2,pv,kfv,k9v,kfq,srt)

      nop=13
c...Option for only BB collisions.
      if(mstc(51).eq.0.and.(k9v(1).ne.3.or.k9v(2).ne.3)) goto 100

      nop=14
c...Option for only BB MB collisions, i.e. no mm collisions.
      if(mstc(51).eq.1.and.(k9v(1).eq.0.and.k9v(2).eq.0)) goto 100

c...Box by maru
      dx=r(1,i2)-r(1,i1)+rcell(1,icell)
      dy=r(2,i2)-r(2,i1)+rcell(2,icell)
      dz=r(3,i2)-r(3,i1)+rcell(3,icell)
      rsqare=dx**2+dy**2+dz**2

c...Skip pair if separation of the two partons is too large.
c     r2=rsqare
c    $    +(dx*(px1+py1)+dy*(py1+py2)+dz*(pz1+pz2))**2/(srt*(e1+e2))
c     if (r2.gt.parc(31)**2) go to 100


c...Does this particle have constituent quarks within a formation time?
c     iqc=0
c     if(kp1.le.-11.or.kt1.le.-11) iqc=1

c...Get collision type.
      icltyp=0
      icltyp=jamcltyp(kf1,kf2,ibar1,ibar2)
      sigfac=1.0d0

c...Baryon-baryon(antiB-antiB)
      if(icltyp.eq.1) then  ! b-b/antiB-antiB
c       if(iqc.eq.1.and.srt.lt.parc(61)) goto 100
        rcut=pard(53)
        sigfac=parc(58)

c...Meson-meson
      else if(icltyp.eq.3) then
c       if(iqc.eq.1.and.srt.lt.parc(63)) goto 100
        rcut=pard(55)
        sigfac=parc(60)

c...Meson-baryon
      else if(icltyp.eq.2) then ! m-b

c       kc1=jamcomp(k(2,i1))
c       kc2=jamcomp(k(2,i2))
c       id1=kchg(kc1,5)
c       id2=kchg(kc2,5)
c       if(id1.eq.id_lamb.or.id1.eq.id_lambs
c    $ .or.id1.eq.id_sigm.or.id1.eq.id_sigms) sigfac=parc(59)
c       if(id2.eq.id_lamb.or.id2.eq.id_lambs
c    $ .or.id2.eq.id_sigm.or.id2.eq.id_sigms) sigfac=parc(59)

        sigfac=parc(59)

c       if(iqc.eq.1.and.srt.lt.parc(62)) goto 100
        rcut=pard(54)
      else if(icltyp.eq.5) then  ! hadron-quark
        rcut=2.0d0
      else if(icltyp.eq.6) then ! partion-parton
        rcut=2.0d0
      else if(icltyp.eq.4) then ! antib-b
        rcut=pard(56)
      else
        write(check(1),'(''ibar1 ibar2='',2i7)')ibar1,ibar2
        call jamerrm(30,1,'(jamhit:)fatal error invalid baryon number')
      endif

c....Determine max. cross section and max. impact par.
c....as well as low energy cutoff

c...2020/2/1
      pc1(:)=p(:,i1)
      pc2(:)=p(:,i2)
      em1sq=pc1(4)**2-pc1(1)**2-pc1(2)**2-pc1(3)**2
      em2sq=pc2(4)**2-pc2(1)**2-pc2(2)**2-pc2(3)**2
      if(mstd(91).ge.1) then
        if(MF_on(i1).eq.1) then
        pc1(1)=pc1(1)-potv(1,i1)
        pc1(2)=pc1(2)-potv(2,i1)
        pc1(3)=pc1(3)-potv(3,i1)
        pc1(4)=sqrt(em1sq+pc1(1)**2+pc1(2)**2+pc1(3)**2)
        endif
        if(MF_on(i2).eq.1) then
        pc2(1)=pc2(1)-potv(1,i2)
        pc2(2)=pc2(2)-potv(2,i2)
        pc2(3)=pc2(3)-potv(3,i2)
        pc2(4)=sqrt(em2sq+pc2(1)**2+pc2(2)**2+pc2(3)**2)
        endif
      endif

      em1=sqrt(max(0d0,em1sq))
      em2=sqrt(max(0d0,em2sq))

      cutoff=em1+em2
      if(ibar1*ibar2.eq.9) then
        inucl1=0
        inucl2=0
        if(kf1.eq.2112.or.kf1.eq.2212.or.kf1.eq.3122)inucl1=1
        if(kf1.eq.-2112.or.kf1.eq.-2212.or.kf1.eq.-3122)inucl1=-1
        if(kf2.eq.2112.or.kf2.eq.2212.or.kf2.eq.3122)inucl2=1
        if(kf2.eq.-2112.or.kf2.eq.-2212.or.kf2.eq.-3122)inucl2=-1
        if(inucl1*inucl2.eq.1)  then
          cutoff=em1+em2+parc(38)
          rcut=pard(51)
        endif
      endif

      if( (kf1.eq.111.or.abs(kf1).eq.211)
     $   .and.(kf2.eq.111.or.abs(kf2).eq.211) ) then
         cutoff=cutoff+parc(38)
          rcut=pard(55)
      endif

      if(srt.ge.4.9d0) rcut=1.6d0
      if(srt.ge.50.d0) rcut=parc(31)

c...Low energy cut off  i.e. Pauli block
      if(srt.lt.cutoff) then
        mstd(52)=mstd(52)+1
        go to 100
      endif

       nop=21
       if(srt.le.em1+em2+parc(41)) then
         goto 100
       endif
       s=srt*srt
       prsq=(s-(em1+em2)**2)*(s-(em1-em2)**2)/(4*s)
       nop=22
c...Too low relative momentum.
       if(prsq.lt.0.000001d0) then  ! pr<0.001GeV/c
         mstd(52)=mstd(52)+1
         goto 100
       endif
       pr=sqrt(prsq)

c...Will particles get closest point in this time interval ?
      t01=r(4,i1)
      t02=r(4,i2)
      dt=t02-t01
      dx12=dt**2-rsqare
      dxp1=dt*pc1(4)-dx*pc1(1)-dy*pc1(2)-dz*pc1(3)
      dxp2=dt*pc2(4)-dx*pc2(1)-dy*pc2(2)-dz*pc2(3)
      dp12=pc1(4)*pc2(4)-pc1(1)*pc2(1)-pc1(2)*pc2(2)-pc1(3)*pc2(3)

      dn=dp12*dp12-em1sq*em2sq
      nop=31
      if(dn.lt.1d-5) goto 100
      b12=dxp1**2*em2sq+dxp2**2*em1sq-2.d0*dxp1*dxp2*dp12
      bsq=-dx12-b12/dn
      brel=sqrt(max(0.0d0,bsq))
      if(brel.gt.rcut) goto 100

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c     pcm(:)=pc1(:)+pc2(:)
c     prl(:)=pc1(:)-pc2(:)
c     pcm(5)=pcm(4)**2 - pcm(1)**2-pcm(2)**2-pcm(3)**2
c     prl(5)=prl(4)**2 - prl(1)**2-prl(2)**2-prl(3)**2
c     ppcm = pcm(4)*prl(4) - pcm(1)*prl(1) - pcm(2)*prl(2)
c    &       - pcm(3)*prl(3)
c     q(:) = prl(:) - ppcm/pcm(5)*pcm(:)
c     q(5) = q(4)**2-q(1)**2-q(2)**2-q(3)**2
c     pcmx = (pcm(4)*dt -pcm(1)*dx - pcm(2)*dy - pcm(3)*dz)**2
c     qx   = (q(4)*dt - q(1)*dx - q(2)*dy - q(3)*dz)**2
c     bbb=-dx12  +  pcmx/pcm(5) + qx / q(5)
c     if(abs(bsq-bbb).gt.1e-7) then
c     print *,'b=',bsq,bbb,p(5,i1),p(5,i2),em1,em2
c     print *,'x=',-b12/dn, pcmx/pcm(5) + qx / q(5)
c     endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

c...Closest distance approach at the two-body c.m. frame.
      if(mstc(52).ne.11) then
        dt1=-pc1(4)*(dxp1*em2sq-dxp2*dp12)/dn
        dt2= pc2(4)*(dxp2*em1sq-dxp1*dp12)/dn
        tcol1=t01+dt1
        tcol2=t02+dt2
c....Time of collision is defined by the time of closet approach
c...in the observational (computational) frame.
      else
       vx1=pc1(1)/pc1(4)
       vy1=pc1(2)/pc1(4)
       vz1=pc1(3)/pc1(4)
       vx2=pc2(1)/pc2(4)
       vy2=pc2(2)/pc2(4)
       vz2=pc2(3)/pc2(4)
       vx=vx2 - vx1
       vy=vy2 - vy1
       vz=vz2 - vz1
       vsq=vx*vx + vy*vy + vz*vz;
       if(vsq.lt.1d-8) goto 100
       dvx=vx*dx + vy*dy + vz*dz;
       dtime=-dvx/vsq
       dt1=dtime + dt*(vx*vx2 + vy*vy2 + vz*vz2)/vsq
       dt2=dtime + dt*(vx*vx1 + vy*vy1 + vz*vz1)/vsq
       tcol1=t01 + dt1
       tcol2=t02 + dt2  ! Actually  tcol1=tcol2
      endif

c...v(4,) is the last collision or decay time.
c...Usually, r(4,)=v(4,) except a previous collision was Pauli-blocked where
c...particle travels without collision but vertex is kept at the last collision
c...point.
c...If time step mstc(3) is greater than 1,
c... r(4,i) is the same as the current time step pare(1)-parc(2)
c...or the time of collision, and different from v(4,).
c...In this case, time step interval is pare(1)-parc(2) < t < pare(1).

      tlast1=v(4,i1)
      tlast2=v(4,i2)

      nop=32
c...Avoid backward collision. Modified by K.Murase May 2013
      if(tcol1.le.tlast1) goto 100
      if(tcol2.le.tlast2) goto 100

c...Define collision ordering time.
      if(mstc(52).eq.2) then
         ctime=0.5d0*(tcol1+tcol2)
      else if(mstc(52).eq.3) then
         ctime=min(tcol1,tcol2)
      else if(mstc(52).eq.4) then
         ctime=max(tcol1,tcol2)
      else if(mstc(52).eq.5) then
         ctime=0.5d0*(tcol1+tcol2)
         tcol1=ctime
         tcol2=ctime
c...2021/8/4
         if(tcol1 .le.tlast1 .or.tcol2.le.tlast2) go to 100
      else if(mstc(52).eq.6) then
         ctime=min(tcol1,tcol2)
         tcol1=ctime
         tcol2=ctime
c...2021/8/4
         if(tcol1 .le.tlast1 .or.tcol2.le.tlast2) go to 100
      else if(mstc(52).eq.11) then
         ctime=tcol1
      else
         ctime=0.5d0*(tcol1+tcol2)
         tcol1=ctime
         tcol2=ctime
      endif

c     if(tcol1.le.tlast1 .or. tcol2.le.tlast2) goto 100
c     if(tcol1.le.t01 .or. tcol2.le.t02) goto 100

      nop=33
      if(mstc(35).eq.1) then
        if(ctime.lt.t01.or.ctime.lt.t02) go to 100
      else if(mstc(35).eq.2) then
        if(ctime.lt.tlast1.or.ctime.lt.tlast2) go to 100
      endif

c...Check max. time.
      nop=34
      if(ctime.gt.pare(1)) go to 100
c...Check min. time.
      nop=35
      if(ctime.lt.pard(20)) go to 100
c     if(ctime.lt.pard(1)) go to 100

c...Avoid collision that will happen after decay.
      if(tcol1.gt.v(5,i1)) goto 100
      if(tcol2.gt.v(5,i2)) goto 100

      facq1=1.0d0
      facq2=1.0d0
c...2017/12/23 in the case of BB collision only, baryons have full cross
c...section to generate Glauber type initial multiple collision during
c...the passage of two nuclei.
c     if(mstc(51).ne.0) then

c....Can const. quark interact within a formation time?
      nop=37
c2017/05/14
      if(t01.lt.r(5,i1).and.tcol1.lt.r(5,i1)) then
c....Additive quark cross section within t_form
          if(abs(ibar1).eq.3) qnum1=3.d0
          if(ibar1.eq.0) qnum1=2.d0
          iqcnum1=mod(abs(kp1)/10,10)
c...skip newly produced hadron which does not have original q.
          if(iqcnum1.eq.0) goto 100
          if(iqcnum1.eq.3) iqcnum1=2
          facq1=iqcnum1/qnum1
      endif

c....Same for targ.
      nop=38
c2017/05/14
      if(t02.lt.r(5,i2).and.tcol2.lt.r(5,i2)) then
        if(abs(ibar2).eq.3) qnum2=3.d0
        if(ibar2.eq.0) qnum2=2.d0
        iqcnum2=mod(abs(kt1)/10,10)
c...skip newly produced hadron which does not have original q.
        if(iqcnum2.eq.0) goto 100
        if(iqcnum2.eq.3) iqcnum2=2
        facq2=iqcnum2/qnum2
      endif

c     endif  ! end 2017/12/23

      pare(54)=facq1
      pare(55)=facq2


c...Get total cross section.
      srt = srt - em1 - em2 + p(5,i1) + p(5,i2)
      em1 = p(5,i1)
      em2 = p(5,i2)
      pr=sqrt((srt**2-(em1+em2)**2)*(srt**2-(em1-em2)**2))/(2*srt)

      call jamcross(1,icltyp,srt,pr,kf1,kf2,em1,em2,
     $                 sig,sigel,sigin,mchanel,mabsrb,isoft,icon)

      nop=40
      if(icon.ne.0) then
          if(mstc(8).ge.2) then
           write(check(1),'(''kf1 em1='',i9,g8.3,i4)')kf1,em1,k(9,i1)
           write(check(2),'(''kf2 em2='',i9,g8.3,i4)')kf2,em2,k(9,i2)
           call jamerrm(1,2,'after cross something was wrong')
           goto 100
          endif
      endif
      if(sigel.gt.sig.and.mstc(8).ge.2) then
        write(mstc(38),*)'kf1 kf2 sigel>sig',kf1,kf2,sigel,sig
      endif
      if(sig.le.0.0d0) then
        if(mstc(8).ge.2)
     $  write(mstc(38),*)'WARNING(cross:1):kf1 kf2 sig<0.0d0'
     $          ,kf1,kf2,kp1,kt1,em1,em2,srt,sig,icltyp
        sig=0.0d0
        sigel=0.0d0
      endif
      nop=41
      if(sig.le.0.0d0) goto 100

      reds=1.0d0
c...Full ensemble method.
      if(mstc(46).eq.2) reds=1d0/mstc(5)

c...All hadrons can collide elastically.
c     if(mstc(59).ge.30.and.(facq1.le.0d0.or.facq2.le.0d0)) then
c       sigel=sig
c       goto 150
c     endif

c....Treatment of the collision involving virtual quarks.
      if(mstc(58).ne.1) reds=reds*facq1*facq2

c...Is their impact parameter small enough?
      bcmax=0.1d0*sig*sigfac*reds/paru(1)
      if(mstc(19).eq.0) then
c       bcmax=sqrt(0.1d0*sig*sigfac*reds/paru(1))
        if(brel.gt.sqrt(bcmax)) go to 100
      else if(mstc(19).eq.1) then
        if(rn(0).gt.parc(39)*exp(-parc(39)*bsq/bcmax)) goto 100
      else if(mstc(19).eq.2) then
        if(rn(0).gt.sigprob(bcmax,bsq)) goto 100
      else
        if(rn(0).gt.sigprob2(bcmax*paru(1),bsq)) goto 100
      endif

c...Additive quark cross section.
      if(mstc(58).eq.1) then

      if(facq1.lt.0.9d0) then
        nop=42
        if(rn(0).gt.facq1) goto 100
      endif

      if(facq2.lt.0.9d0) then
        nop=43
        if(rn(0).gt.facq2) goto 100
      endif

      endif

c...Collision will occur.
150   jamhit=.true.
      mstd(123)=0
      return

c....Collision will not occur.
100   continue
      jamhit=.false.
c     print *,'nop=',nop,i1,i2
c     if(nop.ne.11) write(30,*)nop,i1,i2
      return

      end

c***********************************************************************

      real*8 function sigprob(sig,rr)
      implicit none
      include 'jam2.inc'
      real*8 sig,rr,a,b,sigeff,tpp
c...There sig is defined by the cross section divided by pi.

      sigprob = parc(39)*sig*exp(-parc(39)*rr/(4*parc(15)))/(4*parc(15))

      end

c***********************************************************************

      real*8 function sigprob2(sig,rr)
      implicit none
      include 'jam2.inc'
      real*8 sig,rr,a,b,sigeff,tpp

c...This is the fit of sigma_eff for L=1 fm^2.
      a=0.0651736
      b=1.13922
      sigeff=a*sig**b

      tpp=parc(39)*exp(-parc(39)*rr/(4*parc(15)))/(4*paru(1)*parc(15))
      sigprob2=1.0 - exp(-sigeff*tpp)


      end

c***********************************************************************

      subroutine jamclist

c...Purpose: to make collision list into collision predictor arrays.
c...Modify for box.
 
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      logical jamhit

      ctimemax=0d0
c....Hydro
c     if(mstc(140).ge.1) then
c       pare1=pare(1)
c       parc2=parc(2)
c       pare(1)=mstc(3)*parc(2) ! =nstep*dt: final time
c       parc(2)=pare(1)
c     endif

crqmd
      isw_CanColl=1
      if(isw_CanColl.eq.1) then
        do i1=1,nv
           call jamhit1(i1,ICC)
           ICanColl(i1)=ICC
        enddo
      endif

      mentry=0
c.....Box by maru loop over cell.
      do icell=1,mstd(15)
c...Loop over all particles.
        do 100 i1=2,nv
crqmd...[AO:050727
          if(isw_CanColl.eq.1) then
            if(ICanColl(i1).eq.0) goto 100
          else
c...]AO:050727
          if(mstc(46).ne.2.and.k(8,i1).ne.mstd(8)) goto 100
c...[AO:030128
            if(k(1,i1).ge.11) goto 100
            kf1=k(2,i1)
            kf1a=abs(kf1)
            if(kf1a.eq.0) goto 100
            if((kf1a.gt.10.and.kf1a.le.100).and.kf1.ne.21) goto 100
c...]AO:030128
          endif

          do 200 i2=1,i1-1
c...[AO:050727
          if(isw_CanColl.eq.1) then
            if(ICanColl(i2).eq.0) goto 200
          else
c...]AO:050727
          if(mstc(46).ne.2.and.k(8,i2).ne.mstd(8)) goto 200
c...[AO:030128
            if(k(1,i2).ge.11) goto 200
            kf2=k(2,i2)
            kf2a=abs(kf2)
            if(kf2a.eq.0) goto 200
            if((kf2a.gt.10.and.kf2a.le.100).and.kf2.ne.21) goto 200
c...]AO:030128
          endif

          ctime=0.0d0
          tcol1=0.0d0
          tcol2=0.0d0
          bsq=0.0d0
          sig=0.0d0
          sigel=0.0d0
          if(jamhit(i1,i2,icell,ctime,tcol1,tcol2,bsq,sig,sigel)) then
            mentry=mentry+1
            if(mentry.gt.mxcoll)
     $      call jamerrm(30,0,'(jamclist:), collision storage exceeded')
            icoll(1,mentry)=i1
            icoll(2,mentry)=i2
            icoll(3,mentry)=icell  !....Box by maru
            coll(1,mentry)=ctime
            coll(2,mentry)=tcol1
            coll(3,mentry)=tcol2
            coll(4,mentry)=sig
            coll(5,mentry)=sigel
            coll(6,mentry)=bsq
c           ih=mstc(38)
c           write(22,*)mentry,i1,i2,ctime,tcol1,tcol2,sig,sigel
            ctimemax=max(ctimemax,ctime)
          end if
  200     continue
  100   continue
      end do

c...Save max. number of collision entry.
      mstd(82)=max(mstd(82),mentry)
      pard(23)=ctimemax

c....Hydro
c     if(mstc(140).ge.1) then
c       pare(1)=pare1
c       parc(2)=parc2
c     endif

      end

c***********************************************************************

      subroutine jamcfind(ind,rtype)

c..modified by K.Murase May 2013
c...Purpose: to search decay and collision arrays to find next event
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      character*3 rtype

      integer nmins
      dimension imins(mxcoll)

      rtype='end'
      tmin=pare(10)
      imin=0
      if(mentry.le.0) goto 2000

c...Search collision array
      nmins=0
      do 100 i=1,mentry
        if(icoll(1,i).gt.0) then
c         if(coll(1,i).le.tmin.and.coll(1,i).gt.pard(1)+0.0001d0) then
          if(coll(1,i).eq.tmin) then
            if(k(1,icoll(1,i)).gt.10.or.k(1,icoll(2,i)).gt.10) goto 100
            nmins=nmins+1
            imins(nmins)=i
          else if(coll(1,i).lt.tmin) then
            if(k(1,icoll(1,i)).gt.10.or.k(1,icoll(2,i)).gt.10) goto 100
            tmin=coll(1,i)
            nmins=1
            imins(nmins)=i
          end if
        end if
  100 continue
      if(nmins.gt.0) then
        imin=imins(1+int(rn(0)*nmins))
      end if

c...Search decay array for next decay event.

 2000 tdec=pare(10)
      jmin=0
      if(mstc(55).eq.0)then          ! keeps the resonance stable
      do 101 i=1,nv
        if(k(1,i).gt.10) goto 101    ! already absorbed particle
        if(v(5,i).le.tdec) then
          tdec=v(5,i)
          jmin=i
        end if
  101 continue
      end if
 
      if(tmin.le.tdec) then
        ind=imin
        rtype='col'
      else if(jmin.ne.0) then ! 2014/2/8
        ind=jmin
        rtype='dec'
      end if
      if(ind.eq.0) rtype='end'

      end

c***********************************************************************

      subroutine jamcupda(i1,i2,j1,j2,ifd)

c...Purpose: to update collision array.
c...Modify for box.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      logical jamhit, jamrqpb



c...Cell boundary check.
      if(mstc(4).eq.10) then
        drcell=pard(21)
        do jt=1,2
        if(jt.eq.1)ip=i1
        if(jt.eq.2)ip=i2
        if(ip.gt.0)then
          r(1,ip)=mod(mod(r(1,ip),drcell)+drcell*1.5d0,drcell)-drcell/2
          r(2,ip)=mod(mod(r(2,ip),drcell)+drcell*1.5d0,drcell)-drcell/2
          r(3,ip)=mod(mod(r(3,ip),drcell)+drcell*1.5d0,drcell)-drcell/2
        endif
        end do
      endif
c...end box.

c...Remove collisions including i1,i2.
      ie=0
      do 100 i=1,mentry
        if((icoll(1,i).eq.i1).or.(icoll(2,i).eq.i1).or.
     $     (icoll(1,i).eq.i2).or.(icoll(2,i).eq.i2).or.
     $     (icoll(1,i).eq.-1)) goto 100
        ie=ie+1
        icoll(1,ie)=icoll(1,i)
        icoll(2,ie)=icoll(2,i)
        icoll(3,ie)=icoll(3,i)  ! Box by maru
        coll(1,ie)=coll(1,i)
        coll(2,ie)=coll(2,i)
        coll(3,ie)=coll(3,i)
        coll(4,ie)=coll(4,i)
        coll(5,ie)=coll(5,i)
        coll(6,ie)=coll(6,i)
  100 continue
      mentry=ie


c#ifdef _CUDA
c      if(mstc(129).ne.0) then
c        call jamc_update_cuda(i1,i2,j1,j2,ifd)
c        return
c      endif
c#endif
      if(ifd.eq.0) return

c....Hydro
c     if(mstc(140).ge.1) then
c       pare1=pare(1)
c       parc2=parc(2)
c       pare(1)=mstc(3)*parc(2) ! =nstep*dt: final time
c       parc(2)=pare(1)
c     endif

      do 400 jt=1,2
        if(jt.eq.1) then
          ip=i1
          iq=j1
        endif
        if(jt.eq.2) then 
          ip=i2
          iq=j2
        endif
        if(ip.le.0) goto 400
        if(k(1,ip).gt.10) goto 400


crqmd...[ AO:050727 (MF off for collided particles until next time slice)
        MF_on(ip)=0
c...YN 2015/2/12 option of MF on 
        if(mstc(6).ge.100) then
        if(mstc(107).ge.1) then
           MF_on(ip)=1
           if(jamrqpb(ip)) MF_on(ip)=0
           if(MF_on(ip).eq.1.and.iq.ge.1) then
             do m=1,3
               force(m,ip)=force(m,iq)
               forcer(m,ip)=forcer(m,iq)
             end do
           endif
        endif
        endif
crqmd

c...Box by maru
        do 300 icell3=1,mstd(15)
        do 200 i3=1,nv
          if(mstc(46).ne.2.and.k(8,i3).ne.mstd(8)) goto 200
          if(i3.eq.i1) go to 200
          if(i3.eq.i2) go to 200
crqmd
          call jamhit1(i3,ICanColl3)
          if(ICanColl3.eq.0) goto 200
          ctime=0.0d0
          tcol1=0.0d0
          tcol2=0.0d0
          bsq=0.0d0
          sig=0.0d0
          sigel=0.0d0
crqmd
          if(jamhit(ip,i3,icell3,ctime,tcol1,tcol2,bsq,sig,sigel)) then
            mentry=mentry+1
            if(mentry.gt.mxcoll)
     $      call jamerrm(30,0,'(jamcupda:) collision storage exceeded')
            icoll(1,mentry)=ip
            icoll(2,mentry)=i3
            icoll(3,mentry)=icell3
            coll(1,mentry)=ctime
            coll(2,mentry)=tcol1
            coll(3,mentry)=tcol2
            coll(4,mentry)=sig
            coll(5,mentry)=sigel
            coll(6,mentry)=bsq
          end if
  200   continue
  300   continue
  400 continue

c...Save max. collision entry.
      mstd(82)=max(mstd(82),mentry)

c....Hydro
c     if(mstc(140).ge.1) then
c       pare(1)=pare1
c       parc(2)=parc2
c     endif

      end

#ifdef _CUDA
c***********************************************************************
      subroutine jamc_update_cuda(i1,i2,j1,j2,ifd)

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      logical jamrqpb

      if(i1.ge.1) call jampupdate_cuda(i1)
      if(i2.ge.1) call jampupdate_cuda(i2)

      if(ifd.eq.0) return

      do 400 jt=1,2
        if(jt.eq.1) then
          ip=i1
          iq=j1
        endif
        if(jt.eq.2) then 
          ip=i2
          iq=j2
        endif
        if(ip.le.0) goto 400
        if(k(1,ip).gt.10) goto 400

        MF_on(ip)=0
        if(mstc(107).ge.1) then
           MF_on(ip)=1
           if(jamrqpb(ip)) MF_on(ip)=0
           if(MF_on(ip).eq.1.and.iq.ge.1) then
             do m=1,3
               force(m,ip)=force(m,iq)
               forcer(m,ip)=forcer(m,iq)
             end do
           endif
        endif
  400 continue

      call  jamcoll_update_cuda(i1,i2)

      end
#endif

crqmd
c***********************************************************************

      subroutine jamupdcB

c...Purpose: to update collision array for particles with Mean Field.
c...Modify for box.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      logical jamhit,jamrqpb

      do i=1,nv
        MF_on(i)=1
        if(jamrqpb(i)) MF_on(i)=0
        call jamhit1(i,ICanColl1)
        ICanColl(i)=ICanColl1
      enddo
c...Cell boundary check.
      if(mstc(4).eq.10) then
        drcell=pard(21)
        do ip=1,nv
          if(ip.gt.0.and.MF_on(ip).eq.1)then
           r(1,ip)=mod(mod(r(1,ip),drcell)+drcell*1.5d0,drcell)-drcell/2
           r(2,ip)=mod(mod(r(2,ip),drcell)+drcell*1.5d0,drcell)-drcell/2
           r(3,ip)=mod(mod(r(3,ip),drcell)+drcell*1.5d0,drcell)-drcell/2
          endif
        end do
      endif
c...end box.

c...Remove collisions including particles feeling MF
      mentry0=mentry
      ie=0
      do 100 i=1,mentry
        i1=icoll(1,i)
        i2=icoll(2,i)
        if(i1.lt.0.or.i2.lt.0) goto 100
        if(MF_on(i1).eq.1.or.MF_on(i2).eq.1) goto 100
        if(ICanColl(i1).eq.0.or.ICanColl(i2).eq.0) goto 100
        ie=ie+1
        icoll(1,ie)=icoll(1,i)
        icoll(2,ie)=icoll(2,i)
        icoll(3,ie)=icoll(3,i)  ! Box by maru
        coll(1,ie)=coll(1,i)
        coll(2,ie)=coll(2,i)
        coll(3,ie)=coll(3,i)
        coll(4,ie)=coll(4,i)
        coll(5,ie)=coll(5,i)
        coll(6,ie)=coll(6,i)
  100 continue
      mentry=ie

      do 400 i1=2,nv
        if(ICanColl(i1).eq.0) goto 400
c...Box by maru
        do 300 icell3=1,mstd(15)
        do 200 i2=1,i1-1
          if(ICanColl(i2).eq.0) goto 200
          if(MF_on(i1).eq.0.and.MF_on(i2).eq.0) goto 200
          ctime=0.0d0
          tcol1=0.0d0
          tcol2=0.0d0
          bsq=0.0d0
          sig=0.0d0
          sigel=0.0d0
          if(jamhit(i1,i2,icell3,ctime,tcol1,tcol2,bsq,sig,sigel)) then
            mentry=mentry+1
            if(mentry.gt.mxcoll)
     $      call jamerrm(30,0,'(jamupdc:) collision storage exceeded')
            icoll(1,mentry)=i1
            icoll(2,mentry)=i2
            icoll(3,mentry)=icell3
            coll(1,mentry)=ctime
            coll(2,mentry)=tcol1
            coll(3,mentry)=tcol2
            coll(4,mentry)=sig
            coll(5,mentry)=sigel
            coll(6,mentry)=bsq
          end if
  200   continue
  300   continue
  400 continue
c      write(*,*) mentry0,mentry

c...Save max. collision entry.
      mstd(82)=max(mstd(82),mentry)

      end

c***********************************************************************

      subroutine jamtrspt(ip,tv)

c...Purpose: to transport particles ahead in time.
c...modified for rqmd 2014/11/12
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      integer ip,nn0,nn1,i
      real*8 tv,ee,dt,rv1,rv2,rv3,pvx,pvy,pvz,emfsq,pp

c...transport all particles
      if(ip.eq.0) then
        nn0=1
        nn1=nv
c...transport one particle
      elseif(ip.ge.1.and.ip.le.nv) then
        nn0=ip
        nn1=ip
      else
        print *,'jamtrsprt',ip,nv
        stop
        return
      endif
 
      do 300 i=nn0,nn1
        dtfree(i)=0.0d0
c.....Different Simulation Event when Test Particle Number > 1
        if(mstc(46).ne.2.and.k(8,i).ne.mstd(8)) goto 300
c.....Already dead particle.
        if(k(1,i).ge.11) goto 300
c.....Time to be evolved
        dt=tv-r(4,i)
        if(dt.le.0.0d0) goto 300


        ee=p(4,i)
        rv1=p(1,i)/ee
        rv2=p(2,i)/ee
        rv3=p(3,i)/ee
c...2020/2/1
        if(mstd(91).ge.1.and.MF_on(i).eq.1) then
          pvx=p(1,i)-potv(1,i)
          pvy=p(2,i)-potv(2,i)
          pvz=p(3,i)-potv(3,i)
          ee=sqrt((p(5,i)+pots(i))**2 + pvx**2 + pvy**2 + pvz**2)
          rv1=pvx/ee
          rv2=pvy/ee
          rv3=pvz/ee
        endif

c...[ transport with force
        if(mstc(6).ge.2.and.mstc(103).ge.1) then
cc...[AO:030108
cc force (ix,i)= - dH/dr(ix,i)
cc forcer(ix,i)=   dH/dp(ix,i) - p(ix,i)/p(4,i)
cc We omit force for particles produced or already collided
cc after the time of force evaluation.
c           if(tforce.gt.v(4,i)) then
cc...]AO:030108
c...[AO:050727
          if(MF_on(i).eq.1) then
c...]AO:050727
            emfsq=p(4,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2
            p(1,i)=p(1,i)+dt*force(1,i)
            p(2,i)=p(2,i)+dt*force(2,i)
            p(3,i)=p(3,i)+dt*force(3,i)
            pp=p(1,i)**2+p(2,i)**2+p(3,i)**2
            ee=sqrt(emfsq+pp)
            p(4,i)=ee
c           if(mstc(109).eq.1) then
c             ee=sqrt(max(0d0,ee**2+2*p(5,i)*vpot(i)))
c           endif
            rv1=(rv1+p(1,i)/ee)/2+forcer(1,i)
            rv2=(rv2+p(2,i)/ee)/2+forcer(2,i)
            rv3=(rv3+p(3,i)/ee)/2+forcer(3,i)
          else
            force(1,i)=0.0d0
            force(2,i)=0.0d0
            force(3,i)=0.0d0
            forcer(1,i)=0.0d0
            forcer(2,i)=0.0d0
            forcer(3,i)=0.0d0
            if(ip.eq.0.and.r(5,i).lt.tv) then
              dtfree(i)=tv-r(5,i)
            endif
          endif
        endif
c...] transport with force

            r(1,i)=r(1,i)+dt*rv1
            r(2,i)=r(2,i)+dt*rv2
            r(3,i)=r(3,i)+dt*rv3
            r(4,i)=tv
            call jamtupda(i)

c         if(r(4,i).ge.r(5,i)) then
c           k(1,i)=mod(abs(k(1,i)),10)
c           r(5,i)=r(4,i)
c          endif
c         if(k(1,i).ge.0) r(5,i)=r(4,i)

c           if(ip.ne.0) v(4,i)=tv
 300  continue

c     if(ip.ne.0) then
c      if(r(4,ip).lt.tv) then
c       write(mstc(38),*) 'After transport',ip,r(4,ip),tv
c      endif
c     endif

      end

c***********************************************************************

      subroutine jamtupda(i)

c...Update time of particle.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'

c......has const. quark.
      if(k(1,i).lt.0) then
         if(r(4,i).ge.r(5,i)) then
           k(1,i)=mod(abs(k(1,i)),10)
           v(4,i)=r(5,i)
           r(5,i)=r(4,i)
         endif
      else
        r(5,i)=r(4,i)
      endif

      end

c***********************************************************************

      function jamcltyp(kf1,kf2,ibar1,ibar2)

c...Purpose: to determine type of collision.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'


      kfl1=mod(abs(kf1)/10,10)
      kfl2=mod(abs(kf2)/10,10)
      iq1=0
      iq2=0
      if(kfl1.eq.0.or.kf1.eq.21.or.abs(kf1).le.10) iq1=1
      if(kfl2.eq.0.or.kf2.eq.21.or.abs(kf2).le.10) iq2=1
c     if(iq2.eq.0.and.ibar2.le.2) iq2=1
c     if(iq1.eq.0.and.ibar1.le.2) iq1=1

c...Hadron-hadron collisions.
      if(iq1.eq.0.and.iq2.eq.0) then
         goto 10
c...Parton-parton collisions
      else if(iq1.ne.0.and.iq2.ne.0) then
         jamcltyp=6
         return
c...Parton-hadron collisions.
      else
         jamcltyp=5
         return
      endif

10    continue
      if(ibar1*ibar2.eq.9) then               ! B-B
        jamcltyp=1
      else if(ibar1.eq.0.and.ibar2.eq.0) then ! M-M
          jamcltyp=3
      else if(ibar1*ibar2.eq.0) then          ! M-B
          jamcltyp=2
      else if(ibar1*ibar2.eq.-9) then         ! AntiB-B
        jamcltyp=4
      else
        write(check(1),'(i9,1x,i9,4(i4,1x))')
     $                         kf1,kf2,ibar1,ibar2,iq1,iq2
        call jamerrm(30,1,
     $'(jamcltyp:)fatal error invalid kf1 kf2 ibar1 ibar2 iq1 iq2')
      endif

      end

c***********************************************************************

      subroutine jamscatt

c...Purpose: to perform two-body collision.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'

      dimension pcm(5),pv(2,5),kfv(2),k9v(2),kfq(2,2)
      dimension pd(10,5),kd(10,2)
      logical jamrqpb2

      i1=mste(21)
      i2=mste(23)
      if(i1.le.0.or.i2.le.0) call jamerrm(30,0,'(jamscatt:)i1,i2=0')
      ks1=k(1,i1)
      ks2=k(1,i2)
      kf1=k(2,i1)
      kf2=k(2,i2)
      ibar1=k(9,i1)
      ibar2=k(9,i2)
      mste(32)=ibar1
      mste(33)=ibar2
      kc1=jamcomp(kf1)
      kc2=jamcomp(kf2)
      mste(22)=kc1
      mste(24)=kc2
      nv0=nv
      nmeson0=nmeson
      mste(25)=i1
      mste(27)=i2
      k6a=k(6,i1)
      k6b=k(6,i2)

      facq1=1.0
      facq2=1.0
      if(pare(8).lt.r(5,i1)) then
        if(abs(ibar1).eq.3) qnum1=3.d0
        if(ibar1.eq.0) qnum1=2.d0
        iqcnum1=mod(abs(ks1)/10,10)
        if(iqcnum1.eq.3) iqcnum1=2
        facq1=iqcnum1/qnum1
      endif
      if(pare(9).lt.r(5,i2)) then
        if(abs(ibar2).eq.3) qnum1=3.d0
        if(ibar2.eq.0) qnum1=2.d0
        iqcnum1=mod(abs(ks2)/10,10)
        if(iqcnum1.eq.3) iqcnum1=2
        facq2=iqcnum1/qnum1
      endif
      pare(54)=facq1
      pare(55)=facq2
      
c...Total momentum and energy.
      call jamcmom(i1,i2,pv,kfv,k9v,kfq,srt)

      icltyp=jamcltyp(kfv(1),kfv(2),k9v(1),k9v(2))
      mste(2)=icltyp
      do j=1,4
      pcm(j)=pv(1,j)+pv(2,j)
      end do
      s=pcm(4)**2-(pcm(1)**2+pcm(2)**2+pcm(3)**2)
      if(s<0) then
        print *,'jamscatt s<0? ',s,pcm(1),pcm(2),pcm(3),pcm(4)
        stop
      endif


c....Save initial momenta of colliding particles.
c     pare(31)=pv(1,1)
c     pare(32)=pv(1,2)
c     pare(33)=pv(1,3)
c     pare(34)=pv(2,1)
c     pare(35)=pv(2,2)
c     pare(36)=pv(2,3)

      srt=sqrt(s)
      pcm(5)=srt
      pare(2)=srt
c     pv15=pv(1,4)**2-pv(1,1)**2-pv(1,2)**2-pv(1,3)**2
c     pv25=pv(2,4)**2-pv(2,1)**2-pv(2,2)**2-pv(2,3)**2
c     if(pv15.le.0d0) then
c       write(mstc(38),*)'hadron mass? kf1=',kf1,pv15,p(5,i1)
c       pv15=0d0
c     else
c       pv15=sqrt(pv15)
c     endif
c     if(pv15.le.0d0) then
c       write(mstc(38),*)'hadron mass? kf2=',kf2,pv25,p(5,i2)
c       pv25=0d0
c     else
c       pv25=sqrt(pv25)
c     endif

      pv15=pv(1,5)
      pv25=pv(2,5)

      pr=sqrt((s-(pv15+pv25)**2)*(s-(pv15-pv25)**2))/(2*srt)
      em1=p(5,i1)
      em2=p(5,i2)
      pr0=pr
      pare(7)=pr
      pxr0=pv(1,1)
      pyr0=pv(1,2)
      pzr0=pv(1,3)
      per0=pv(1,4)
      pare(12)=pcm(1)/pcm(4)
      pare(13)=pcm(2)/pcm(4)
      pare(14)=pcm(3)/pcm(4)
      pare(15)=pcm(4)/pcm(5)
      call jamrobo(0.0d0,0.0d0,-pare(12),-pare(13),-pare(14),pare(15)
     $ ,pxr0,pyr0,pzr0,per0)
      pare(37)=pxr0
      pare(38)=pyr0
      pare(39)=pzr0
      pxr1=-pxr0
      pyr1=-pyr0
      pzr1=-pzr0
      pare(16)=pjangl(pxr0,pyr0)
      call jamrobo(0d0,-pare(16),0d0,0d0,0d0,1.0d0,pxr0,pyr0,pzr0,per0)
      pare(17)=pjangl(pzr0,pxr0)
      call jamrobo(-pare(17),0d0,0d0,0d0,0d0,1.0d0,pxr0,pyr0,pzr0,per0)

      tc1=pare(8)
      tc2=pare(9)
      call jamrobo(0.0d0,0.0d0,-pare(12),-pare(13),-pare(14),pare(15),
     &    pare(43),pare(44),pare(45),tc1)
      call jamrobo(0d0,-pare(16),0d0,0d0,0d0,1.0d0,
     &    pare(43),pare(44),pare(45),tc1)
      call jamrobo(-pare(17),0d0,0d0,0d0,0d0,1.0d0,
     &    pare(43),pare(44),pare(45),tc1)
      call jamrobo(0.0d0,0.0d0,-pare(12),-pare(13),-pare(14),pare(15),
     &       pare(46),pare(47),pare(48),tc2)
      call jamrobo(0d0,-pare(16),0d0,0d0,0d0,1.0d0,
     &       pare(46),pare(47),pare(48),tc2)
      call jamrobo(-pare(17),0d0,0d0,0d0,0d0,1.0d0,
     &       pare(46),pare(47),pare(48),tc2)

      pare(25)=pare(43)-pare(46)
      pare(26)=pare(44)-pare(47)
      pare(27)=pare(45)-pare(48)


c     print *,'x1=',pare(43),pare(44),pare(45),tc1
c     print *,'x2=',pare(46),pare(47),pare(48),tc2
c     print *,'dx dy dz=',pare(25),pare(26),pare(27)

c       if(abs(pare(27)).ge.1d-3) then
c         print *,'pare(27)=',pare(27)
c         print *,kf1,kf2
c         print *,'em1=',em1,'em2=',em2
c         print *,'em1*=',pv15,'em2*=',pv25
c         print *,'v1',v(4,i1),'v2=',v(4,i2)
c         print *,'t1',r(4,i1),'t2=',r(4,i2)
c         print *,'p=',pxr0,pyr0,pzr0
c         print *,'r1=',xr0,yr0,zr0,tr0
c         print *,'r2=',xr1,yr1,zr1,tr1
c         print *,'gam=',pare(15)
c         print *,'the phi=',pare(16),pare(17)
c       endif

c     pare(43)=xr0
c     pare(44)=yr0
c     pare(45)=zr0

c     pare(46)=xr1
c     pare(47)=yr1
c     pare(48)=zr1

c     gam1=p(4,i1)/sqrt(p(4,i1)**2-(p(1,i1)**2+p(2,i1)**2+p(3,i1)**2))
c     gam2=p(4,i2)/sqrt(p(4,i2)**2-(p(1,i2)**2+p(2,i2)**2+p(3,i2)**2))
c     pare(49)=(r(4,i1)-v(4,i1))/gam1
c     pare(50)=(r(4,i2)-v(4,i2))/gam2


c...Information.
      if(mstc(8).ge.5) then
        ih=mstc(38)
        write(ih,*)'Collision will occure srt',
     $ srt,k(1,i1),k(2,i2),i1,i2,kf1,kf2,em1,em2
      endif

      isetm=0
      ichanel=1
      mste(1)=ichanel
      pare(3)=0.0d0  ! xsig
      taucol=0.0d0
      mste(29)=0
      mste(31)=0

c2016/6
      pare(21)=pv15/em1
      pare(22)=pv25/em2
c     if(mstc(109).ge.1.and.mstc(6).ge.101) then
c       pare(21)=0d0
c       pare(22)=0d0
c       if(.not.jamrqpb2(kf1)) pare(21)=vpot(i1)
c       if(.not.jamrqpb2(kf2)) pare(22)=vpot(i2)
c       if(mstc(106).ge.201) then
c         pare(21)=1d0+pare(21)
c         pare(22)=1d0+pare(22)
c       else
c         pare(21)=sqrt(max(0.14d0,1d0 + 2*pare(21)/em1))
c         pare(22)=sqrt(max(0.14d0,1d0 + 2*pare(22)/em2))
c       endif
c     endif
      pare(23)=pare(21)
      pare(24)=pare(22)



c     if(abs(pare(21)-1d0).gt.1d-4.or.abs(pare(22)-1d0).gt.1d-4) then
c            print *,'before coll.',mstc(109),mstc(6)
c            print *,'em1ef em1',kf1,em1,pv15,pare(21)
c            print *,'em2ef em2',kf2,em2,pv25,pare(22)
c     endif

c---------------------------------------------------------------------
c...Determine hard/soft/low energy inelastic channels

      call jamchanl(i1,i2,pcm,srt,pr,kf1,kf2
     $ ,kc1,kc2,em1,em2,pv15,pv25,ibar1,ibar2,ichanel,isoft,n_jet)
      mste(1)=ichanel

c...Information.
      if(mstc(8).ge.5) then
        ih=mstc(38)
c        if(mstc(106).ge.201) then
c          emf1=em1*(1+pare(23))
c          emf2=em2*(1+pare(24))
         if(mstc(108).ge.3) then
           emf1=em1+pare(23)
           emf2=em2+pare(24)
         else
           emf1=sqrt(max(0d0,em1**2+2*em1*pare(23)))
           emf2=sqrt(max(0d0,em2**2+2*em2*pare(24)))
         endif
        write(ih,*)'after jamchanl',ichanel,srt,i1,i2,kf1,kf2,em1,em2,
     &    emf1,emf2,'pare23=',pare(23),pare(24)
      endif

c....Compute energy-momentum tensor at the collision point.
      mste(42)=100
      itens=0
      itens2=0
      if(mstc(170).ge.1.or.mstc(59).ge.30) itens=1
      if((ichanel.ge.1.and.ichanel.le.5).and.ichanel.ne.3) itens2=1

      if(itens.ge.1.and.itens2.ge.1) then
        call jamtensor(icon,i1,i2,cc,rhop,rhobar,einv,pfree,pxx,pyy,pzz)
        mste(42)=icon
        pare(31)=einv
        pare(32)=rhop
        pare(33)=rhobar
        pare(34)=pfree
        pare(35)=(pxx+pyy)/2
        pare(36)=pzz
      endif

 1200 continue
      if(ichanel.eq.4.or.ichanel.eq.5) then

        call jamsoft(icltyp,i1,i2,icon)

        if(mstc(8).ge.5) then
          write(ih,*)'after jamsoft',icon,srt,i1,i2,kf1,kf2,em1,em2
        endif

        if(icon.eq.5000) then ! New diffractive scattering.
          if(mstc(170).ge.1) call jameosout(i1,i2,1)
          goto 6000
        endif

        if(icon.ne.0) then
          mste(1)=-77
          return
        endif

        kf1=k(2,i1)
        kf2=k(2,i2)
        kc1=jamcomp(kf1)
        kc2=jamcomp(kf2)
        if(mstc(170).ge.1) call jameosout(i1,i2,1)
        goto 5500

      else if(ichanel.eq.1.or.ichanel.eq.2.or.ichanel.ge.11) then

        if(ichanel.eq.1) then
          kf1=kcp(2,1)
          kf2=kcp(2,2)
          em1=pcp(5,1)
          em2=pcp(5,2)
        endif
   
        if(ichanel.le.2) then

c         if(mstc(109).le.1) then
c         call jamscat2(i1,i2,pcm,kf1,kf2,kc1,kc2,em1,em2,ibar1,ibar2,
c    $      pr0,pr,taucol)
c         else

          call jamscat2m(i1,i2,pcm,kf1,kf2,kc1,kc2,em1,em2,ibar1,ibar2,
     $      pr0,pr,taucol)

c         endif
          if(mste(1).le.0) return

          if(mstc(170).ge.1) call jameosout(i1,i2,0)


          isetm=2
c....NN -> d pi
        else if(ichanel.ge.12.and.ichanel.le.14) then
          call jamscatd(i1,i2,pcm,kf1,kf2,kc1,kc2,em1,em2,pr0,kd,pd)
c....NN -> d 2pi
        else if(ichanel.ge.15.and.ichanel.le.18) then
          call jamscatd2(i1,i2,pcm,kf1,kf2,kc1,kc2,em1,em2,kd,pd)
        else
          call jamscatm(i1,i2,pcm,kf1,kf2,kc1,kc2,em1,em2,kd,pd)
          if(mste(1).le.0) return
        endif

      else if(ichanel.eq.3) then

        call jamabsrb(pcm,kf1,kc1,ks1,ks2,ibar1,isoft,icltyp)

c....test xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c       if(icltyp.ne.4.or.isoft.ne.2) then
c          if(abs(k(7,i1)).ne.1.and.abs(k(7,i2)).ne.1) then
c          print *,'kf1=',kf1,kf2,kcp(2,1),kcp(2,2),pcm(5)
c          print *,'i1=',i1,i2,mste(21),mste(23)
c          print *,'k7=',k(7,i1),k(7,i2)
c          mste(1)=0
c          endif
c       endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


        if(mste(1).eq.0) return
        isetm=1
        if(icltyp.eq.4) goto 6000  ! antiBB annihilation (one string)

      else if(ichanel.eq.6) then

          if(k(9,i1).eq.0.or.k(9,i2).eq.0) then
             print *,'before jamhard meson?',k(9,i1),k(9,i2)
             print *,'kf1=',k(2,i1),k(2,i2)
             stop
          endif

        call jamhard(n_jet,i1,i2,jflg)
        if(jflg.le.0) then
          ichanel=5
          mste(1)=ichanel
          goto 1200 
        endif
        goto 6000

c...antibaryon-baryon annihilation
      else if(ichanel.eq.7) then

c       call jamabbann(srt,pr,kf1,kf2,em1,em2,icon)
        call jamabbann(srt,pr,kf1,kf2,pv15,pv25,icon)
        if(icon.ne.0) then
          mste(1)=-66
          return
        endif
        goto 6000

      else if(ichanel.le.0) then
        mste(1)=ichanel
        return
      else
        write(check(1),'(''ichanel='',i6)')ichanel
        call jamerrm(30,1,'(jamscatt:) invalid ichanel')
      endif

c----------------------------------------------------------------------*

c...Set the new particle information after collision
      mstd(29)=mstd(29)+1
      call jamkupda(1,mste(25),kf1,kc1,ks1,k6b,mstd(29),mste(1))
      call jamkupda(1,mste(27),kf2,kc2,ks2,k6a,mstd(29),mste(1))
c     if(mste(1).ge.11) then
      if(mste(29).ge.1) then
        call jamkupda(1,mste(29),kd(3,2),kd(3,1),1
     $ ,k(6,mste(25)),mstd(29),mste(1))
      endif
      if(mste(31).ge.1) then
        call jamkupda(1,mste(31),kd(4,2),kd(4,1),1
     $ ,k(6,mste(25)),mstd(29),mste(1))
      endif

c...qmd:Recalculate momenta of the colliding particles
c...in order to recover total energy conservation in case of
c...potential forces act.
c     if(mstc(110).ge.1.and.isetm.ge.1) then
c       call jamecon(isetm)
c     end if

      i3=mstc(25)
      i4=mstc(27)
      em3f=sqrt(p(4,i3)**2-p(1,i3)**2-p(2,i3)**2-p(3,i3)**2)
      if(i4.ge.1) then
        em4f=sqrt(p(4,i4)**2-p(1,i4)**2-p(2,i4)**2-p(3,i4)**2)
      else
        em4f=0.0d0
        p(5,i4)=0.0d0
      endif

c     if(abs(em3f-p(5,i3)).gt.1d-4.or.abs(em4f-p(5,i4)).gt.1d-4) then
c      print *,'after scatter ',ichanel,srt
c      print *,'i3',k(2,i3),em3f,p(5,i3)
c      print *,'i4',k(2,i4),em4f,p(5,i4)
c     endif

c...Information.
      if(mstc(8).ge.10) then
        ih=mstc(38)
        j1=mstc(25)
        j2=mstc(27)
        write(ih,*)'before pauli',ichanel,srt,i1,i2,kf1,kf2,em1,em2
        if(j1.ne.0)write(ih,*)j1,k(1,j1),k(2,j1),kq(1,j1)
        if(j2.ne.0)write(ih,*)j2,k(1,j2),k(2,j2),kq(1,j2)
      endif

c----------------------------------------------------------------------*
c
c   Check on Pauli-blocking
c
c----------------------------------------------------------------------*
5500  continue
      ntag=0
c...Check mass of string system. convert into hadron
      if(mste(25).ge.1) then
      if(k(10,mste(25)).ne.0) then
        call str2had(0,mste(25),kf1,kc1,em1,ic1)
        if(ic1.eq.2.or.ic1.eq.1) then
         ntag=1
         mste(1)=-1
         if(mstc(8).ge.5) then
         ih=mstc(38)
         write(ih,*)'parton system mass becomes small after scatt.1',i1
         write(ih,*)'ichanel icltyp',ichanel,icltyp
         write(ih,*)'i1 k1 k2',i1,k(1,i1),k(2,i1),p(5,i1),em1
         write(ih,*)'i2 k1 k2',i2,k(1,i2),k(2,i2),p(5,i2)
         endif
         return
        endif
      endif
      endif
      if(mste(27).ge.1) then
      if(k(10,mste(27)).ne.0) then
        call str2had(0,mste(27),kf2,kc2,em2,ic1)
        if(ic1.eq.2.or.ic1.eq.1) then
         ntag=1
         mste(1)=-1
         if(mstc(8).ge.5) then
         ih=mstc(38)
         write(ih,*)'parton system mass becomes small after scatt.2',i2
         write(ih,*)'ichanel icltyp',ichanel,icltyp
         write(ih,*)'i1 k1 k2',i1,k(1,i1),k(2,i1),p(5,i1)
         write(ih,*)'i2 k1 k2',i2,k(1,i2),k(2,i2),p(5,i2),em2
         endif
         return
        endif
      endif
      endif

c....Check Pauli-blocking for nucleons.
      if(mstc(56).ne.0) then
        if(kf1.eq.2112.or.kf1.eq.2212) then
          call jampauli(mste(25),ntag,phase)
        end if
        if((kf2.eq.2112.or.kf2.eq.2212).and.(ntag.eq.0)) then
          call jampauli(mste(27),ntag,phase)
        end if

c...This collision was Pauli-blocked.
        if(ntag.eq.1) then
          mste(1)=-1
          nv=nv0
          nmeson=nmeson0
          return
        end if
      end if

c...Information.
      if(mstc(8).ge.5) then
        ih=mstc(38)
        i3=mste(25)
        i4=mste(27)
        write(ih,*)'after pauli ichanel srt=',mste(1),srt
        if(i3.ge.1)then
         write(ih,*)'i3 k1 k2 em',i3,k(1,i3),k(2,i3),p(5,i3)
         Write(3,*)'i1',kq(1,i3),kq(2,i3),(vq(j,i3),j=1,10)
        endif
        if(i4.ge.1) then
          write(ih,*)'i4 k1 k2 em',i4,k(1,i4),kf2,em2
          Write(3,*)'i2',kq(1,i4),kq(2,i4),(vq(j,i4),j=1,10)
        endif
      endif

c----------------------------------------------------------------------*
c
c           Collision Was Successful !
c
c----------------------------------------------------------------------*

c...Store momentum and flavor of the jet system.
      if(mste(25).ge.1) then
        if(mod(abs(k(1,mste(25))),10).eq.3) then
c.....Recalculate relative momentum, because of setmom.
          if(mste(1).ge.3.and.mste(1).le.6) then
          else
            j1=mste(25)
            call jamjetm2(mste(25))
c           print *,'1',mste(1),k(2,j1),p(5,j1)
          endif
        else
          j1=mste(25)
          kq(1,j1)=999999
          kq(2,j1)=0
          vq(1,j1)=pare(12)   ! beta_x
          vq(2,j1)=pare(13)   ! beta_y
          vq(3,j1)=pare(14)   ! beta_z
          vq(4,j1)=pare(15)   ! gamma
          vq(5,j1)=pare(16)   ! phi
          vq(6,j1)=pare(17)   ! theta
          vq(7,j1)=pare(25)   ! xr0-xr1
          vq(8,j1)=pare(26)   ! yr0-yr1
          vq(9,j1)=pare(27)   ! zr0-zr1
          vq(10,j1)=pare(49)  ! tau_1

c....Take care of s-channel scattering.
          if(mste(27).eq.0) then
            vq(10,j1)=pare(49)+pare(50)
            if(ibar1.eq.0) then
              vq(7,j1)=-pare(25)
              vq(8,j1)=-pare(26)
              vq(9,j1)=-pare(27)
              pare(16)=pjangl(pxr1,pyr1)
              pare(17)=pjangl(pzr1,sqrt(pxr1**2+pyr1**2))
              vq(5,j1)=pare(16)   ! phi
              vq(6,j1)=pare(17)   ! theta
            endif
          endif


        endif
crqmd[ AO:050727 (MF off for collided particles until next time slice)
c         MF_on(mste(25))=0   ! move to jamcupda
c...] AO:050727
      endif

      if(mste(27).ge.1) then
        if(mod(abs(k(1,mste(27))),10).eq.3) then
          if(mste(1).ge.3.and.mste(1).le.6) then
          else
            call jamjetm2(mste(27))
c           print *,'2',mste(1),k(2,j2),p(5,j2)
          endif
        else
          j2=mste(27)
          kq(1,j2)=999999
          kq(2,j2)=0
          vq(1,j2)=pare(12)
          vq(2,j2)=pare(13)
          vq(3,j2)=pare(14)
          vq(4,j2)=pare(15)
          vq(5,j2)=pare(16)
          vq(6,j2)=pare(17)
          vq(7,j2)=-pare(25)
          vq(8,j2)=-pare(26)
          vq(9,j2)=-pare(27)

          if(abs(ibar1).eq.3) then
            vq(7,j2)=pare(25)
            vq(8,j2)=pare(26)
            vq(9,j2)=pare(27)
          endif

          vq(10,j2)=pare(50)
        endif
c...[ AO:050727 (MF off for collided particles until next time slice)
c         MF_on(mste(27))=0   ! move to jamcupda
c...] AO:050727
      endif

c...Update life time and vertex.
      call jamvupda(mste(25),mste(25),kf1,kc1,taucol,mste(1))
      call jamvupda(mste(27),mste(27),kf2,kc2,taucol,mste(1))
c     if(mste(1).ge.11)
      if(mste(29).ge.1)
     $  call jamvupda(mste(29),0,kd(3,2),kd(3,1),taucol,mste(1))
      if(mste(31).ge.1)
     $  call jamvupda(mste(31),0,kd(4,2),kd(4,1),taucol,mste(1))

c...Information.
      if(mstc(8).ge.5) then
        ih=mstc(38)
        j1=mste(25)
        j2=mste(27)
        write(ih,*)'after jamvupda',mste(1),srt
        if(j1.ge.1) then
        write(ih,*)'i3',j1,k(1,j1),k(2,j1),p(5,j1),(v(j,j1),j=1,5)
        write(ih,*)kq(1,j1),(vq(j,j1),j=1,5)
        write(ih,*)kq(2,j1),(vq(j,j1),j=6,10)
        endif
        if(j2.ge.1) then
        write(ih,*)'i4',j2,k(1,j2),k(2,j2),p(5,j2),(v(j,j2),j=1,5)
        write(ih,*)kq(1,j2),(vq(j,j2),j=1,5)
        write(ih,*)kq(2,j2),(vq(j,j2),j=6,10)
        endif
      endif

 6000 continue

c...Store particle information.
      if(mste(1).eq.3) then
        call jamsave(1,3,mste(25))
      else
        call jamsave(1,3,mste(25))
        if(mste(27).ge.1) call jamsave(1,4,mste(27))
        if(mste(29).ge.1) call jamsave(1,5,mste(29))
        if(mste(31).ge.1) call jamsave(1,6,mste(31))
      endif
       
c...Count collision number.
      if(mstc(162).eq.1) call jamclhst(1,mste(1))

      end

c***********************************************************************

      subroutine jamchanl(i1,i2,pcm,srt,pr,kf1,kf2
     $ ,kc1,kc2,em1,em2,emf1,emf2,ibar1,ibar2,ichanel,isoft,n_jet)

c...Determine hard/soft/resonance.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      parameter(mxchan=30)
      dimension pcm(5),sigin(mxchan)

      n_jet=0 
      sig=pare(4)
      sigel=pare(5)
      bsq=pare(6)
      icltyp=mste(2)

c....In case of AA simulation,
      if(mstd(11).ne.2) then

c....Sorry not implemented m-m/m-B hard scattering.
      if(icltyp.eq.2.or.icltyp.eq.3) goto 2200
      if(icltyp.eq.5.or.icltyp.eq.6) goto 2200
      if(k(2,i1).eq.21) goto 2200
      if(k(2,i2).eq.21) goto 2200

      endif

c....Hard or soft
c....sharp transition
      if(srt.ge.parc(71).and.icltyp.le.4) then
c....smooth transition 2014/2/12
        y=1
        if(parc(72).gt.parc(71)) 
     $        y=log(srt/parc(71))/log(parc(72)/parc(71))

        if(rn(0).le.y) then
          n_jet=jamisjet(pcm,bsq,i1,i2,sig,sigel)
          if(n_jet.eq.0) then
            ichanel=4
          else if(n_jet.ge.1) then
            ichanel=6
          else if(n_jet.eq.-1) then
           ichanel=1
          else
           ichanel=-77
          endif
          return
        endif

      endif
c---------------------------------------------------------------------

2200  continue
c....In case of 1+1 simulation, no elastic if desired.
c....(NOTE: only works for BB collisions).
      if(mstc(17).eq.1) then
        sig=sig-sigel
        pare(4)=sig
        pare(5)=0.0d0
        pare(3)=sig*rn(0)
        goto 2000
      endif

c...String interaction forbid elastic.
      if((kf1.eq.92.or.kf2.eq.92)
     $  .or.(k(1,i1).eq.4.or.k(1,i2).eq.4)) then
        ichanel=4
        return
      endif

      pare(3)=sig*rn(0)-sigel
c     call gover(b2,sig,sigel,gsin,gstot)

c...An elastic collision will happen.
      if(pare(3).le.0.0d0) then
         ichanel=1
         return
      endif

c======================================================================*
c
c     Start Inelastic Collisions
c
c======================================================================*
2000  continue
      ichanel=2

c...String interaction.
      if(kf1.eq.92.or.kf2.eq.92) then
        ichanel=4
        return
      endif

c...Hadron-quark
      if(icltyp.eq.5) then
        ichanel=4
c....g-hadron scatt.
        if(kf1.eq.21.or.kf2.eq.21) then
          ichanel=1
        endif
        return
c...Parton-parton
      else if(icltyp.eq.6) then
        ichanel=1
        return
      endif

c...Save total charge in order to check charge conservation.
      iz10=kchg(kc1,1)*isign(1,kf1)
      iz20=kchg(kc2,1)*isign(1,kf2)
      izt0=iz10+iz20
      kf01=kf1
      kf02=kf2
      em01=em1
      em02=em2

c...Get inel. cross sections at low energy.
c     srtf=sqrt(em1**2+pr**2)+sqrt(em2**2+pr**2)
c     pe1=sqrt(em1**2+p(1,i1)**2+p(2,i1)**2+p(3,i1)**2)
c     pe2=sqrt(em2**2+p(1,i2)**2+p(2,i2)**2+p(3,i2)**2)
c     s=(pe1+pe2)**2
c    $ -(p(1,i1)+p(1,i2))**2-(p(2,i1)+p(2,i2))**2-(p(3,i1)+p(3,i2))**2
c     prf=sqrt((s-(em1+em2)**2)*(s-(em1-em2)**2)/(4*s))

      s0=srt - emf1 - emf2 + em1 + em2 
      pr0=sqrt((s0**2-(em1+em2)**2)*(s0**2-(em1-em2)**2))/(2*s0)

c     print *,'m=',em1,em2
c     print *,'m=',emf1,emf2
c     print *,'s=',srt,s0
c     print *,'p=',pr,pr0
c     read(5,*)

      call jamcross(3,icltyp,srt,pr0,kf1,kf2,em1,em2,
     $           sig,sigel,sigin,mchanel,mabsrb,isoft,icon)

      if(mstc(8).ge.5) then
        ih=mstc(38)
        write(ih,*)'jamchanl after jamcross3',
     $        icon,srt,i1,i2,kf1,kf2,em1,em2
       endif

c...Energetically forbidden.
      if(icon.eq.1) then
        ichanel=-88
        return
      else if(icon.ne.0) then
         if (.true.) then ! ONO
           ichanel=1 ! elastic as fallback
           return
         endif

        write(check(1),8000)icon,srt
        write(check(2),8100)chaf(kc1,(3-isign(1,kf1))/2)
     $                     ,chaf(kc2,(3-isign(1,kf2))/2)
        write(check(3),8200)kf1,em1,k(9,i1)
        write(check(4),8300)kf2,em2,k(9,i2)
        call jamerrm(1,4,'(after jamcross3:)inel. x-section wrong')
        ichanel=-99
        return
      end if

c         if((kf1.eq.2212.or.kf1.eq.2112).and.em1.le.0.9) then
c           print *,'after cross3 kf01=',kf01,em01
c           print *,'kf1=',kf1,em1
c         endif
c         if((kf2.eq.2212.or.kf2.eq.2112).and.em2.le.0.9) then
c           print *,'after cross3 kf02=',kf02,em02
c           print *,'kf2=',kf2,em2
c         endif

      if(pare(3).le.0.0d0) then

c....Soft interaction.
        if(isoft.eq.1) then
          ichanel=4
          return
c...antibaryon-baryon annihilation
        else if(isoft.eq.7) then
          ichanel=7
          return
c....Elastic.
        else if(isoft.eq.-1) then
          ichanel=1
          return
        endif

        if(kf2.eq.0) then
          ichanel=3
        else if(kf1.ne.0.and.kf2.ne.0) then
          ichanel=2
          if(mste(3).ge.1.and.mste(3).le.3) ichanel=11
          if(mste(3).eq.11) ichanel=12
          if(mste(3).eq.12) ichanel=13
          if(mste(3).eq.13) ichanel=14
          if(mste(3).eq.14) ichanel=15
          if(mste(3).eq.15) ichanel=16
          if(mste(3).eq.16) ichanel=17
          if(mste(3).eq.17) ichanel=18
        else
           write(check(1),'(''kf1 kf2='',i9,1x,i9)')kf1,kf2
           call jamerrm(30,1,'(jamscatt:) ???kf1 kf2=')
        endif

c...Nothing happens.
      else
        if(isoft.eq.1) then
          ichanel=4
        else if(isoft.eq.2) then
          ichanel=5
        else
          ichanel=0
        endif
        return
      endif

      kc1=jamcomp(kf1)
      kc2=jamcomp(kf2)
      mste(26)=kc1
      mste(28)=kc2
      pare(41)=em1
      pare(42)=em2

c...Check charge conservation after coll.
      iz1=kchg(kc1,1)*isign(1,kf1)
      if(kc2.gt.0) then
        iz2=kchg(kc2,1)*isign(1,kf2)
      else
        iz2=0
      endif
      if(izt0.ne.iz1+iz2) then
        write(check(1),8500)icltyp,ichanel
        write(check(2),8600)iz10,iz20,iz1,iz2
        write(check(3),8700)kcp(2,1),kcp(2,2),kf1,kf2,em1,em2
        call jamerrm(30,3,
     $   '(jamchanl:) Charge not conserved after jamcross(3)')
      endif

 8000 format('after jamcross(3) icon=',i8,' srt=',g10.3)
 8100 format('collision of ',a16,'+ ',a16)
 8200 format('kf1 em1 k9=',i9,1x,g15.3,1x,i5)
 8300 format('kf2 em2 k9=',i9,1x,g15.3,1x,i5)
 8500 format('icltyp ichanel',i8,1x,i8)
 8600 format('iz01 iz02 =>iz1 iz2',4(i5,1x))
 8700 format(i9,1x,i9,' ==> ',i9,1x,i9,' em1 em2',g10.3,1x,g10.3)

      end

c***********************************************************************

      subroutine jamabsrb(pcm,kf1,kc1,ks1,ks2,ibar1,isoft,icltyp)

c...Treat annihilation scattering.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      dimension pcm(5),pcmf(5)
      real*8 jamemjet

      i1=mste(21)
      i2=mste(23)
      taucol=0d0
      pcmf(5)=pcm(5)

c2015/4/18
c       if(ks1.le.0.or.ks2.le.0) then
c          mste(1)=0
c          return
cc         k(1,i3)=-k(1,i3)
cc         ks0b=mod(abs(ksb)/10,10)
cc         ks0m=mod(abs(ksm)/10,10)
cc         if(ks0b.eq.2.and.ks0m.eq.1) then
cc         else if(ks0b.eq.1) then
cc            if(ks0m.eq.1) k(1,i3)=-20+k(1,i3)
cc            if(ks0m.eq.0) k(1,i3)=-10+k(1,i3)
cc         else if(ks0b.eq.2) then
cc            if(ks0m.eq.0) k(1,i3)=-20+k(1,i3)
cc         endif
c       endif


c.....AntiB B annihilation.
      if(icltyp.eq.4) then
        nmeson=nmeson+1
        nv=nv+1
        if(nv.gt.mxv) then
          call jamerrm(30,0,'(jamabsrb:)Particle too large [mxv]')
        endif
        i3=nv
        i4=0

c...MB/MM absorption.
      else

c       if(mstc(109).ge.1) then

          emf1=sqrt(p(4,i1)**2-p(1,i1)**2-p(2,i1)**2-p(3,i1)**2)
          emf2=sqrt(p(4,i2)**2-p(1,i2)**2-p(2,i2)**2-p(3,i2)**2)
          sm1=emf1/p(5,i1)
          sm2=emf2/p(5,i2)

          do j=1,4
            pcmf(j)=p(j,i1)/sm1+ p(j,i2)/sm2
          end do
          pcmf(5)=sqrt(pcmf(4)**2-pcmf(1)**2-pcmf(2)**2-pcmf(3)**2)
          ss=pcm(5)/pcmf(5)

c         pkx1=p(1,i1)-potv(1,i1)       
c         pky1=p(2,i1)-potv(2,i1)       
c         pkz1=p(3,i1)-potv(3,i1)       
c         pkx2=p(1,i2)-potv(1,i2)       
c         pky2=p(2,i2)-potv(2,i2)       
c         pkz2=p(3,i2)-potv(3,i2)       
c         ek1=sqrt(emf1**2+pkx1**2+pky1**2+pkz1**2)
c         ek2=sqrt(emf2**2+pkx2**2+pky2**2+pkz2**2)
c         ee=ek1+ek2
c         sf=ee**2 - pcm(1)**2 - pcm(2)**2 -pcm(3)**2
c         print *,'sf=',sqrt(sf),pcm(5)

c         eef=sqrt(p(5,i1)**2+p(1,i1)**2+p(2,i1)**2+p(3,i1)**2)
c    &      +sqrt(p(5,i2)**2+p(1,i2)**2+p(2,i2)**2+p(3,i2)**2)
c         do j=1,3
c           ee = ee - (p(j,i1)+p(j,i2))**2
c         end do
c         print *,k(2,i1),k(2,i2),kf1
c         print *,'fmas=',sqrt(max(0.0,ee)),pcmf(5),p(5,i1)+p(5,i2)


          if(abs(pcmf(5)-pcm(5)).gt.1d-5) then
            if(kchg(kc1,5).eq.id_delt) then
              emin=1.08090d0+parc(41)
            else
              call jamdmass(kf1,kfm1,kfd1,emin1,emdn1,ist)
c             emin=pmas(kc1,1)-pmas(kc1,3)
              emin=emdn1
            endif
            pcmf(5)=max(emin,pcmf(5))
          endif

c         if(ss.le.0.0.or.ss.gt.1.0) then
c         if(abs(ss-1d0).gt.1.0d-4) then
c             print *,'absorb S=',ss,pcm(5),pcmf(5)
c             pcmf(5)=pcm(5)
c         endif

c     if(emf1.le.0d0.or.emf2.le.0d0) then
c     print *,'kf1=',k(2,i1),'kf2=',k(2,i2)
c     print *,'sm1=',sm1,'sm2=',sm2
c     print *,'pot1=',vpot(i1),'pot2=',vpot(i2)
c     print *,'m1=',p(5,i1),'m2=',p(5,i2)
c     print *,'m*1=',emf1,'m*2=',emf2
c     print *,'m*1=',p(5,i1)*(1+vpot(i1)),'m*2=',p(5,i2)*(1+vpot(i2))
c     print *,'kf=',kf1,'m=',pcm(5)
c     print *,'absorb kf1',kf1,pcmf(5),pcm(5),ss
c     print *,'M=',pcm(5)/sm1,pcm(5)/sm2,(pcm(5)/sm1+pcm(5)/sm2)/2
c     print *,'pcm2(5)=',pcmf(5)
c     endif


        if(abs(ibar1).eq.3) then
          i3=i1
          i4=0
          call jamzero(i2)
          if(mstc(6).ge.0) call jamcupda(i2,-1,0,0,0)
c         if(mstc(6).ge.0) call jamcupda(i1,i2,0,0,0)
c....2017/3/22 set time of resonance.
c         rcp(4,1)=min(rcp(4,1),rcp(4,2))
c         rcp(5,1)=min(rcp(5,1),rcp(5,2))
        else
          i3=i2
          i4=0
          ks1=k(1,i2)
          call jamzero(i1)
          if(mstc(6).ge.0) call jamcupda(i1,-1,0,0,0)
c         if(mstc(6).ge.0) call jamcupda(i1,i2,0,0,0)
c....2017/3/22 set time of resonance.
c         rcp(4,2)=min(rcp(4,1),rcp(4,2))
c         rcp(5,2)=min(rcp(5,1),rcp(5,2))
        endif

      endif

      kc1=jamcomp(kf1)
      k(2,i3)=kf1
      p(1,i3)=pcm(1)
      p(2,i3)=pcm(2)
      p(3,i3)=pcm(3)
      p(5,i3)=pcmf(5)

c...Canonical momentum.
      if(mstd(91).eq.1) then
        if(MF_on(i1).eq.1) then
        p(1,i3)=p(1,i3)+potv(1,i1)
        p(2,i3)=p(2,i3)+potv(2,i1)
        p(3,i3)=p(3,i3)+potv(3,i1)
        endif
        if(MF_on(i2).eq.1) then
        p(1,i3)=p(1,i3)+potv(1,i2)
        p(2,i3)=p(2,i3)+potv(2,i2)
        p(3,i3)=p(3,i3)+potv(3,i2)
        endif
      endif
      p(4,i3)=sqrt(pcm(5)**2+p(1,i3)**2+p(2,i3)**2+p(3,i3)**2)


      mste(25)=i3
      mste(27)=i4

c===========================================================================
      istrcheck=0
      if(istrcheck.eq.1) then
c...Check string mass. Currently we ignore the mass shift of string.
      call attflv(kf1,jq1,jqq1)
      emjet1=jamemjet(jq1,jqq1)
      if(pcm(5).le.emjet1) then
        isoft=0
        kf01=kf1
        call jamidres(2,kf1,pcm(5),icon)
        kc1=jamcomp(kf1)
        if(pmas(kc1,2).le.1d-5) then
           print *,'absorb stable paticle?',kf01,kf1,pmas(kc1,2)
           print *,'icon=',icon,'pcmf5=',pcmf(5),pcm(5)
           stop
        endif
c       if(icon.ne.0..and.icon.ne.10) then
c         print *,'absorb kf1=',kf1,kf01
c         print *,'icon=',icon,'pcmf5=',pcmf(5),pcm(5)
c       endif
        k(2,i3)=kf1
      endif
      endif
c===========================================================================

c....String formation
      if(isoft.eq.2) then
        call jamjetm2(i3)
        kf1=92
        kc1=92
        p(1,i3)=pare(12)  ! bex
        p(2,i3)=pare(13)  ! bey
        p(3,i3)=pare(14)  ! bez
        p(4,i3)=pare(16)  ! phi
        p(5,i3)=pare(17)  ! the
        v(1,i3)=pare(15)  ! gam

        v(2,i3)=pare(49)+pare(50) ! taucl
        v(3,i3)=pare(25)  ! dxcm
        v(4,i3)=pare(26)  ! dycm
        v(5,i3)=pare(27)  ! dzcm

c       v(1,i3)=pare(16)  ! phi
c       v(2,i3)=pare(17)  ! the
c       v(3,i3)=pare(15)  ! gam
c       v(4,i3)=pare(25)  ! dxcm
c       v(5,i3)=pare(26)  ! dycm

        if(ibar1.eq.0) then
              pare(16)=pjangl(-pare(37),-pare(38))
              pare(17)=pjangl(-pare(39),sqrt(pare(37)**2+pare(38)**2))
              p(4,i3)=pare(16)   ! phi
              p(5,i3)=pare(17)   ! theta

              v(3,i3)=-pare(25)  ! dxcm
              v(4,i3)=-pare(26)  ! dycm
              v(5,i3)=-pare(27)  ! dzcm
c     pare(37)=pxr0
c     pare(38)=pyr0
c     pare(39)=pzr0
        endif

      else
c       if(mstc(109).ge.2) then
c         em0 = pcm(5)/pare(21)
c         p(5,i3)=max(eminres(kc1),em0)
c         p(5,i3)=pcm(5)/pare(21)
c       endif
      endif

      if(icltyp.eq.4) then
        r(1,i3)=(r(1,i1)+r(1,i2))/2.d0
        r(2,i3)=(r(2,i1)+r(2,i2))/2.d0
        r(3,i3)=(r(3,i1)+r(3,i2))/2.d0
        call jamzero(i1)
        call jamzero(i2)
        if(mstc(6).ge.0) then
          call jamcupda(i1,-1,0,0,0)
          call jamcupda(i2,-1,0,0,0)
        endif
        mstd(29)=mstd(29)+1
        call jamkupda(1,i3,kf1,kc1,ks1,k(6,i3),mstd(29),3)
        call jamvupda(i3,0,kf1,kc1,taucol,3)
      endif

      end

c***********************************************************************
      subroutine jamabbann(srt,pr,kf1,kf2,em01,em02,icon)
c...Purpose: to treat annihilation of antibaryon-baryon collision
c     implicit double precision(a-h, o-z)
c     implicit integer(i-n)
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      dimension pp(3,5),pt(3,5),kfl1(3),kfl2(3)
      dimension pstr(3,5),em(3),iann(3,3),kfm(3)
c     parameter(ems=0.5)
      real*8 vpot(0:3)

c     ems=0.9d0*parj(32)
      ems=1.0d0*parj(32)
      ipion=1
      icon=0

c....bbar+b -> 3 string 
c....bbar+b -> 2 string + meson
c....bbar+b -> 2 string 

      ianti1=1
      ianti2=1
      if(kf1.lt.0) ianti1=-1
      if(kf2.lt.0) ianti2=-1

c...Quark contents.
      ntry=1
 1000 continue
      ntry=ntry+1
      ns=3
      de=0.0d0
      do i=1,3
        kfm(i)=0
        em(i)=0.0d0
      end do
      if(ntry.ge.100) then
        write(mstc(38),*) 'jamabbann infinit loop?'
        icon=1
        return
      endif
      xr=rn(0)
      if(xr.le.0.33333) then
        l1=1000
        l2=100
        l3=10
      else if(xr.le.0.66666) then
        l1=100
        l2=10
        l3=1000
      else
        l1=10
        l2=1000
        l3=100
      endif
      kf1a=abs(kf1)
      kfl1(3)=mod(kf1a/l1,10) 
      kfl1(2)=mod(kf1a/l2,10) 
      kfl1(1)=mod(kf1a/l3,10) 
      xr=rn(0)
      if(xr.le.0.33333) then
        l1=1000
        l2=100
        l3=10
      else if(xr.le.0.66666) then
        l1=100
        l2=10
        l3=1000
      else
        l1=10
        l2=1000
        l3=100
      endif
      kf2a=abs(kf2)
      kfl2(3)=mod(kf2a/l1,10) 
      kfl2(2)=mod(kf2a/l2,10) 
      kfl2(1)=mod(kf2a/l3,10) 


      nann=0
      do i=1,3
      do j=1,3
         iann(i,j)=0
         if(kfl1(i).eq.kfl2(j)) then
           nann=nann+1
           iann(i,j)=1
         endif
      end do
      end do

      do i=1,3
        kfl1(i)=kfl1(i)*ianti1
        kfl2(i)=kfl2(i)*ianti2
      end do
      
c     write(mstc(38),*)'kf1=',(kfl1(i),i=1,3)
c     write(mstc(38),*)'kf2=',(kfl2(i),i=1,3)

c....annihilation into two meson strings.
      if(nann.ge.1) then
         if(srt.le.2.3) ns=2
         if(rn(0).le.parc(18)) ns=2
      endif

      if(ns.eq.2) then
        kfla=0
        kflb=0
        ns=2
        ir=1+rn(0)*(nann-1)
        l=0
        do i=1,3
        do j=1,3
         if(iann(i,j).eq.1)  l=l+1
         if(ir.eq.l) then
           if(i.eq.1) then
             kfla=kfl1(1)
             kfl1(1)=kfl1(2)
             kfl1(2)=kfl1(3)
           else if(i.eq.2) then
             kfla=kfl1(2)
             kfl1(1)=kfl1(1)
             kfl1(2)=kfl1(3)
           else
             kfla=kfl1(3)
           endif
           if(j.eq.1) then
             kflb=kfl2(1)
             kfl2(1)=kfl2(2)
             kfl2(2)=kfl2(3)
           else if(j.eq.2) then
             kflb=kfl2(2)
             kfl2(1)=kfl2(1)
             kfl2(2)=kfl2(3)
           else
             kflb=kfl2(3)
           endif
           kfl1(3)=0
           kfl2(3)=0
           goto 100
         endif
      end do
      end do
100   continue
       if(kfla*kflb.eq.0.or.abs(kfla).ne.abs(kflb)) then
         write(mstc(38),*)'error kfla kflb',kfla,kflb,ir,l,nann
         stop
       endif

      else if(ipion.eq.1.and.(srt.le.2.15.or.rn(0).le.parc(19))) then
        ns=3
c....find smallest mass of the hadron
       emm=10.0d0
       ilf=0
       do i=1,3
        kfla=kfl1(i)
        kflb=kfl2(i)
        call kfcnst(kfla,kflb,kfh,0.0d0)
        call jamdmass(kfh,kfm1,kfd,emin,emdn,ist)
        if(emin.lt.emm) then
          kfm(3)=kfm1
          emm=emin
          ilf=i
        endif
        end do
        if(ilf.eq.0) then
          print *,'funny ilf=0'
          stop
        endif
        if(ilf.ne.3) then
          kfl1(ilf)=kfl1(3)
          kfl2(ilf)=kfl2(3)
        endif
        kfl1(3)=0
        kfl2(3)=0

      endif

      i1=mste(21)
      i2=mste(23)

      pxcm=p(1,i1)+p(1,i2)
      pycm=p(2,i1)+p(2,i2)
      pzcm=p(3,i1)+p(3,i2)
      pecm=p(4,i1)+p(4,i2)
      s=pecm**2-pxcm**2-pycm**2-pzcm**2
      srt0=sqrt(s)
      pr0=sqrt((s-(em01+em02)**2)*(s-(em01-em02)**2)/(4*s))
      e1=sqrt(em01**2+pr0**2)
      e2=sqrt(em02**2+pr0**2)

c     e1=sqrt(em01**2+pr**2)
c     e2=sqrt(em02**2+pr**2)
      call jamdivide(ns,e1,pr0,1,kfl1,pp,icon)
      if(icon.ne.0) goto 1000
      call jamdivide(ns,e2,-pr0,2,kfl2,pt,icon)
      if(icon.ne.0) goto 1000

      do j=1,4
        pstr(1,j)=pp(1,j)+pt(1,j)
        pstr(2,j)=pp(2,j)+pt(2,j)
        pstr(3,j)=0.0
        if(ns.eq.3) pstr(3,j)=pp(3,j)+pt(3,j)
      end do


c....convert meson string into hadron
      if(kfm(3).ne.0) then
        em(3)=pstr(3,4)**2-pstr(3,1)**2-pstr(3,2)**2-pstr(3,3)**2
        if(em(3).le.0.0) goto 1000
        em(3)=sqrt(em(3))
        if(em(3).lt.emm) then
          write(mstc(38),*)'em3<emm',kfm,em(3),emm
          goto 1000
        endif
        e=sqrt(emm**2+pstr(3,1)**2+pstr(3,2)**2+pstr(3,3)**2)
        de=pstr(3,4)-e
        if(de.lt.0.0d0) goto 1000
        pstr(3,4)=e
        em(3)=emm
        pstr(1,4)=pstr(1,4)+de
        pp(1,4)=pp(1,4)+0.5d0*de
        pt(1,4)=pt(1,4)+0.5d0*de
      endif

      do i=1,ns
        em(i)=pstr(i,4)**2-pstr(i,1)**2-pstr(i,2)**2-pstr(i,3)**2
        if(em(i).gt.0.0) then
          em(i)=sqrt(em(i))
        else
c         write(mstc(38),*)'em<0',i,em(i)
          goto 1000
        endif
c       qm=pjmass(kfl1(i))+pjmass(kfl2(i))
        qm=pp(i,5)+pt(i,5)
        if(kfm(i).eq.0.and.em(i).lt.qm+ems) then ! string mass?
c         if(em(i).lt.qm+0.5) goto 100 ! resonance mass?
c            write(mstc(38),*)'em qm',em(i),qm+0.01
          if(em(i).lt.qm+0.01) goto 1000 ! resonance mass?
          kfla=kfl1(i)
          kflb=kfl2(i)
          emr=em(i)
          call kfcnst(kfla,kflb,kfm(i),0.0d0)
          call jamidres(0,kfm(i),emr,icon)

        endif
      end do

      e= pstr(1,4)+pstr(2,4)+pstr(3,4)
      px=pstr(1,1)+pstr(2,1)+pstr(3,1)
      py=pstr(1,2)+pstr(2,2)+pstr(3,2)
      pz=pstr(1,3)+pstr(2,3)+pstr(3,3)
      s=e**2 - px**2 - py**2 - pz**2

      iconv=0
      if(abs(sqrt(s)-srt0).gt.1e-5
     & .or.abs(px).gt.0.01.or.abs(py).gt.0.01
     & .or.abs(pz).gt.0.01) iconv=1

       if(iconv.eq.1.or.mstc(8).ge.3) then
       write(mstc(38),*)'jamabban= ',mste(21),mste(23)
       write(mstc(38),*)'ntry kf1 kf2= ',ntry,kf1,kf2
       write(mstc(38),*)'em1 em2 em3',em(1),em(2),em(3)
       write(mstc(38),*)'kf= ',(kfm(j),j=1,3)
       write(mstc(38),*)'ns srt= ',ns,srt,sqrt(s)
       write(mstc(38),*)px,py,pz,e
       if(iconv.eq.1) then
c        stop
         icon=1
         return
       endif
      endif

c     betx=pare(12)
c     bety=pare(13)
c     betz=pare(14)
c     gam=pare(15)
c     phi=pare(16)
c     the=pare(17)

      betx=pxcm/pecm
      bety=pycm/pecm
      betz=pzcm/pecm
      gam=pecm/srt0
      pxr0=p(1,i1)
      pyr0=p(2,i1)
      pzr0=p(3,i1)
      per0=p(4,i1)
      call jamrobo(0d0,0d0,-betx,-bety,-betz,gam,pxr0,pyr0,pzr0,per0)
      phi=pjangl(pxr0,pyr0)
      the=pjangl(pzr0,sqrt(pxr0**2+pyr0**2))
c     call jamrobo(0d0,-phi,0d0,0d0,0d0,1.0d0,pxr0,pyr0,pzr0,per0)
c     the1=pjangl(pzr0,pxr0)

      rx=(r(1,i1)+r(1,i2))/2.d0
      ry=(r(2,i1)+r(2,i2))/2.d0
      rz=(r(3,i1)+r(3,i2))/2.d0

c...Remove collision of baryon and antibaryon
      call jamcupda(i1,-1,0,0,0)
      call jamcupda(i2,-1,0,0,0)
      mstd(29)=mstd(29)+1
c...bug fix 2017/5/7
      k(5,i1)=mstd(29)
      k(5,i2)=mstd(29)
      taucol=0.0d0

      vpot(:)=0.0
c     if(mstd(91).eq.1) then
c       if(MF_on(i1).eq.1) vpot(:)=potv(:,i1)/ns
c       if(MF_on(i2).eq.1) vpot(:)=vpot(:)+potv(:,i2)/ns
c     endif

c...Loop over meson strings
      do i=1,ns
      nv=nv+1
      nmeson=nmeson+1
      ks=1
      k(6,nv)=0
      if(kfm(i).ne.0) then
        call jamrobo(the,phi,betx,bety,betz,gam,
     &             pstr(i,1),pstr(i,2),pstr(i,3),pstr(i,4))
        call jamzero(nv)
        k(1,nv)=1
        kf=kfm(i)
        kc=jamcomp(kf)
        ks=1
        p(1,nv)=pstr(i,1) + vpot(1)
        p(2,nv)=pstr(i,2) + vpot(2)
        p(3,nv)=pstr(i,3) + vpot(3)
c       p(4,nv)=pstr(i,4)
        p(4,nv)=sqrt(em(i)**2+p(1,nv)**2+p(2,nv)**2+p(3,nv)**2)
        p(5,nv)=em(i)
        r(1,nv)=rx
        r(2,nv)=ry
        r(3,nv)=rz
        do j=1,3
         v(j,nv)=r(j,nv)
        end do
        kq(1,nv)=0
        kq(2,nv)=0
        do kk=1,10
          vq(kk,nv)=0.0d0
        end do
        kq(1,nv)=999999
        kq(2,nv)=0
        vq(1,nv)=pare(12)
        vq(2,nv)=pare(13)
        vq(3,nv)=pare(14)
        vq(4,nv)=pare(15)
        vq(5,nv)=pare(16)
        vq(6,nv)=pare(17)
        call jamkupda(1,nv,kf,kc,ks,0,mstd(29),3)
        call jamvupda(nv,0,kf,kc,taucol,3)
c       write(mstc(38),*)'jamabban k1 k2=',i,k(1,nv),k(2,nv),p(5,nv)
c       write(mstc(38),*)'r=',(r(j,nv),j=1,5)
c       write(mstc(38),*)'v=',(v(j,nv),j=1,5)
c...bug fix 2017/5/7
        do j=1,3
          v(j,nv)=r(j,nv)
        end do
        call jamcupda(nv,-1,0, 0, 1)

      else
        call jamrobo(the,phi,betx,bety,betz,gam,
     &             pp(i,1),pp(i,2),pp(i,3),pp(i,4))
        call jamrobo(the,phi,betx,bety,betz,gam,
     &             pt(i,1),pt(i,2),pt(i,3),pt(i,4))

        do j=1,3
          pp(i,j)=pp(i,j)+vpot(j)/2
          pt(i,j)=pt(i,j)+vpot(j)/2
        end do
        kq(1,nv)=kfl1(i)
        kq(2,nv)=kfl2(i)
        do j=1,5
          vq(j,nv)=pp(i,j)
          vq(j+5,nv)=pt(i,j)
        end do

c     print *,'mq=',vq(5,nv),vq(10,nv)
c     print *,'px=',vq(1,nv)+vq(6,nv),pp(i,1)+pt(i,1),pstr(i,1)
c     print *,'py=',vq(2,nv)+vq(7,nv),pp(i,2)+pt(i,2),pstr(i,2)
c     print *,'pz=',vq(3,nv)+vq(8,nv),pp(i,3)+pt(i,3),pstr(i,3)
c     print *,'pe=',vq(4,nv)+vq(9,nv),pp(i,4)+pt(i,4),pstr(i,4)


        r(1,nv)=rx
        r(2,nv)=ry
        r(3,nv)=rz
        p(1,nv)=pare(12)
        p(2,nv)=pare(13)
        p(3,nv)=pare(14)
        v(3,nv)=pare(15)
        v(1,nv)=pare(16)
        v(2,nv)=pare(17)
        kf=92
        kc=92
        ks=1
        call jamkupda(1,nv,kf,kc,ks,0,mstd(29),3)
        call jamvupda(nv,0,kf,kc,taucol,3)
c...baryon and/or anti-baryon has const.quark
        v(5,nv)=r(5,nv)
c       if(r(5,nv).gt.v(5,nv)) then
c       print *,'r4=',r(4,nv),v(5,nv),k(1,i1),k(1,i2)
c       endif

c       emj=(vq(4,nv)+vq(9,nv))**2-(vq(1,nv)+vq(6,nv))**2
c    &   -(vq(2,nv)+vq(7,nv))**2-(vq(3,nv)+vq(8,nv))**2
c       write(mstc(38),*)'string',i,kq(1,nv),kq(2,nv),em(i),sqrt(emj)
c       write(mstc(38),*)'kfl1 kfl2',kfl1(i),kfl2(i),ianti1,ianti2

      endif


      end do

      call jamzero(i1)
      call jamzero(i2)

      end

c***********************************************************************
      subroutine jamdivide(ns,pe,pz,jt,kfl,p,icon)
      implicit double precision(a-h, o-z)
      implicit integer(i-n)
      include 'jam2.inc'
      dimension p(3,5),pmt(3),pmt2(3),kfl(3)

      icon=0
      p(3,5)=0.0d0
      mstj(93)=1
      p(1,5)=pjmass(kfl(1))
      mstj(93)=1
      p(2,5)=pjmass(kfl(2))
      if(kfl(3).ne.0) then
        mstj(93)=1
        p(3,5)=pjmass(kfl(3))
      endif
      em1=parf(100+abs(kfl(1)))
      em2=parf(100+abs(kfl(2)))
      em3=0.0
      if(kfl(3).ne.0) em3=parf(100+abs(kfl(3)))

          loop=0
  370     loop=loop+1
          call pjptdi(1,p(1,1),p(1,2))
          if(ns.eq.2) then
            p(2,1)=-p(1,1)
            p(2,2)=-p(1,2)
            pmt(1)=p(1,5)**2+p(1,1)**2+p(1,2)**2
            pmt(2)=p(2,5)**2+p(2,1)**2+p(2,2)**2
            pmt2(1)=em1**2+p(1,1)**2+p(1,2)**2
            pmt2(2)=em2**2+p(2,1)**2+p(2,2)**2
          else
            call pjptdi(1,p(2,1),p(2,2))
            p(3,1)=-p(1,1)-p(2,1)
            p(3,2)=-p(1,2)-p(2,2)
            pmt(1)=p(1,5)**2+p(1,1)**2+p(1,2)**2
            pmt(2)=p(2,5)**2+p(2,1)**2+p(2,2)**2
            pmt(3)=p(3,5)**2+p(3,1)**2+p(3,2)**2
            pmt2(1)=em1**2+p(1,1)**2+p(1,2)**2
            pmt2(2)=em2**2+p(2,1)**2+p(2,2)**2
            pmt2(3)=em3**2+p(3,1)**2+p(3,2)**2
          endif

         if(ns.eq.2) then
          chi1=1d0-sqrt(rn(0))
          emt=pmt2(1)/chi1+pmt2(2)/(1d0-chi1)
         else
 380      chi1=1d0-sqrt(rn(0))
          chi2=1d0-sqrt(rn(0))
          if(1d0-chi1-chi2.lt.0d0) goto 380
          emt=pmt2(1)/chi1+pmt2(2)/chi2+pmt2(3)/(1d0-chi1-chi2)
        endif

          if(emt.gt.pe**2) then
            if(loop.lt.100) then
              goto 370
            else
              icon=-1
              return
            endif
          endif
 
C...Subdivide longitudinal momentum according to value selected above.
          pw1=chi1*(pe+abs(pz))
          p(1,4)=0.5d0*(pw1+pmt(1)/pw1)
          p(1,3)=0.5d0*(pw1-pmt(1)/pw1)*(-1)**(jt-1)
          if(ns.eq.2) then
            p(2,4)=pe-p(1,4)
            p(2,3)=pz-p(1,3)
          else
            pw2=chi2*(pe+abs(pz))
            p(2,4)=0.5d0*(pw2+pmt(2)/pw2)
            p(2,3)=0.5d0*(pw2-pmt(2)/pw2)*(-1)**(jt-1)
            p(3,4)=pe-p(1,4)-p(2,4)
            p(3,3)=pz-p(1,3)-p(2,3)
          endif

c         print *,'chi=',chi1
c         print *,'e=',p(1,4),p(2,4)
c         print *,'pz=',p(1,3),p(2,3)

      end

c***********************************************************************
      subroutine jamdivide2(srt,pe,pz,jt,kfl,p,icon)
      implicit double precision(a-h, o-z)
      implicit integer(i-n)
      dimension p(2,5),pmt(2),kfl(3)

      icon=0
      p(1,5)=pjmass(kfl(1))
      p(2,5)=pjmass(kfl(2))

c     chik=parp(92+2*imb)
c     cut=2d0*0.3d0/vint(1)
      chik=3.0
      cut=2d0*0.3d0/srt
       
          loop=0
  370     loop=loop+1
          call pjptdi(1,p(1,1),p(1,2))
          p(2,1)=-p(1,1)
          p(2,2)=-p(1,2)
          pmt(1)=p(1,5)**2+p(1,1)**2+p(1,2)**2
          pmt(2)=p(2,5)**2+p(2,1)**2+p(2,2)**2

  380       chi=pjr(0)**2
            if((chi**2/(chi**2+cut**2))**0.25d0*
     &        (1d0-chi)**chik.lt.pjr(0)) goto 380
C...Construct total transverse mass; reject if too large.
          chi=max(1d-8,min(1d0-1d-8,chi))
          emt=pmt(1)/chi+pmt(2)/(1d0-chi)
          if(emt.gt.pe**2) then
            if(loop.lt.100) then
              goto 370
            else
              icon=-1
              return
            endif
          endif
 
C...Subdivide longitudinal momentum according to value selected above.
          pw1=chi*(pe+abs(pz))
          p(1,4)=0.5d0*(pw1+pmt(1)/pw1)
          p(1,3)=0.5d0*(pw1-pmt(2)/pw1)*(-1)**(jt-1)
          p(2,4)=pe-p(1,4)
          p(2,3)=pz-p(1,3)
          print *,'chi=',chi
          print *,'e=',p(1,4),p(2,4)
          print *,'pz=',p(1,3),p(2,3)

      end

c***********************************************************************

      subroutine jamscat2(i1,i2,pcm,kf1,kf2,kc1,kc2,em1,em2,ibar1,ibar2,
     $ pr0,pr,taucol)

c...Calculate kinematics for the two-body final case.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      dimension pcm(5)
 
      sig=pare(4)
      srt=pcm(5)

c...Scattering and azimuthal angles in two-body cm-fame.
      if(mste(1).eq.1) then

c.......Get elastic angular distribution.
          call jamangel(pr,sig,kf1,kf2,ibar1,ibar2,t1,c1) 

      else

c.....Compute magnitude of rel. momentum using the new masses.
        pr2=(srt**2-(em1+em2)**2)*(srt**2-(em1-em2)**2)
        if(pr2.gt.0.00001d0) then
          pr=sqrt(pr2)/(2.d0*srt)
        else
          write(check(1),'(g15.3,1x,i9,1x,i9,g10.3,1x,g10.3)')
     $       srt,kf1,kf2,em1,em2 
          call jamerrm(3,1,'(jamscat2:) rel.mom.pr<0 after scatter')
          mste(1)=-99
          return
        endif

c.....Get inelastic angular distribution.
        call jamangin(srt,pr,pr0,em1,em2,ibar1,ibar2,
     $                kcp(2,1),kcp(2,2),mste(22),mste(24)
     $                  ,kc1,kc2,t1,c1)

      end if

      betx=pare(12)
      bety=pare(13)
      betz=pare(14)
      gam=pare(15)
      phi=pare(16)
      the=pare(17)
      pxr=0.0d0
      pyr=0.0d0
      pzr=pr
      pare(18)=t1
      pare(19)=c1

c.....Rotate momentum vector to new position.
      sin1=sqrt(1.d0-c1**2)
      pxr=pr*sin1*cos(t1)
      pyr=pr*sin1*sin(t1)
      pzr=pr*c1

      pare(28)=pxr
      pare(29)=pyr
      pare(30)=pzr-pr0
      pare(20)=pxr*pare(25)+pyr*pare(26)+(pzr-pr0)*pare(27)

c....EoS modified collision.
      if(mstc(59).ge.101) then
         call jameoscol(i1,i2,pr0,pxr,pyr,pzr,em1,em2)

c...Select repulsive or attractive orbit.
      else if(mstc(59).ge.1) then
        dot=pare(20)
        irev=0
        if(mstc(59).eq.30) then
          pare(18)=paru(1)+atan2(pare(26),pare(25))
          pxr=pr*sin1*cos(pare(18))
          pyr=pr*sin1*sin(pare(18))
        else if(mod(mstc(59),10).eq.1) then ! Repulsive
         if(dot.lt.0.0d0) irev=1
        else if(mod(mstc(59),10).eq.2) then ! Attractive 
         if(dot.gt.0.0d0) irev=1
        endif
        if(irev.eq.1) then
          pxr=-pxr
          pyr=-pyr
        endif

        pare(20)=pxr*pare(25)+pyr*pare(26)+(pzr-pr0)*pare(27)
        pare(28)=pxr
        pare(29)=pyr

      endif

c...Estimate collision time scale from uncertainty principle.
c     pt=sqrt(pxr**2+pyr**2)
c     if(pt.gt.1d-7) then
c       taucol=gam*paru(3)/pt
c       taucol=-tau*log(max(rn(0),1.d-35))
c       if(taucol.gt.1.0d0) then
c         write(33,*)mste(1),srt,taucol
c         mste(1)=-9
c         return
c       endif
c     endif
 
      per=sqrt(em1**2+pr**2)
      pxr2=-pxr
      pyr2=-pyr
      pzr2=-pzr
      per2=sqrt(em2**2+pr**2)

c...Back to the original frame.
      call jamrobo(the,phi,betx,bety,betz,gam,pxr,pyr,pzr,per)
      p(1,i1)=pxr
      p(2,i1)=pyr
      p(3,i1)=pzr
      p(4,i1)=sqrt(em1**2+p(1,i1)**2+p(2,i1)**2+p(3,i1)**2)
      p(5,i1)=em1
 
      call jamrobo(the,phi,betx,bety,betz,gam,pxr2,pyr2,pzr2,per2)
      p(1,i2)=pxr2
      p(2,i2)=pyr2
      p(3,i2)=pzr2
      p(4,i2)=sqrt(em2**2+p(1,i2)**2+p(2,i2)**2+p(3,i2)**2)
      p(5,i2)=em2

c     pare(37)=pxr
c     pare(38)=pyr
c     pare(39)=pzr

c     pare(40)=pxr2
c     pare(41)=pyr2
c     pare(42)=pzr2

 
c....Check energy momentum conservation.
      srt1=sqrt((p(4,i1)+p(4,i2))**2-(p(1,i1)+p(1,i2))**2
     $               -(p(2,i1)+p(2,i2))**2-(p(3,i1)+p(3,i2))**2)

      if( (abs(pcm(1)-(p(1,i1)+p(1,i2))).ge.1)
     a    .or. (abs(pcm(2)-(p(2,i1)+p(2,i2))).ge.1)
     a      .or. (abs(pcm(3)-(p(3,i1)+p(3,i2))).ge.1)
     a      .or. (abs(pcm(4)-(p(4,i1)+p(4,i2))).ge.1)
     $        .or. (abs((srt1-srt)/srt).gt.0.1d0) ) then

        ih=mstc(38)
        write(ih,'(/,''<<jamscat2>>'')')
        write(ih,*)'beta gam',betx,bety,betz,gam
        write(ih,*)'ichanel srt sig=',
     $    mste(1),srt,sig,pcm(5)
        write(ih,*)'c1 t1',c1,t1
        write(ih,*)'pr0',pr0
        write(ih,*)'pr pxr pyr pzr',pr,pxr,pyr,pzr
        write(ih,*)'kf1 em1',kcp(2,1),pcp(5,1),(pcp(i,1),i=1,4)
        write(ih,*)'kf2 em2',kcp(2,2),pcp(5,1),(pcp(i,2),i=1,4)
        write(ih,*)'kf1 em1',kf1,em1,(p(i,i1),i=1,3)
        write(ih,*)'kf2 em2',kf2,em2,(p(i,i2),i=1,3)
        write(ih,*)'pcm1',pcm(1),p(1,i1)+p(1,i2)
        write(ih,*)'pcm2',pcm(2),p(2,i1)+p(2,i2)
        write(ih,*)'pcm3',pcm(3),p(3,i1)+p(3,i2)
        write(ih,*)'e2  ',pcm(4),p(4,i1)+p(4,i2)
        write(ih,*)'srt ',srt,srt1,abs(srt1-srt)/srt
        write(ih,*) pcp(4,1),
     &   sqrt(pcp(5,1)**2+pcp(1,1)**2+pcp(2,1)**2+pcp(3,1)**2)
        write(ih,*) pcp(4,2),
     &   sqrt(pcp(5,2)**2+pcp(1,2)**2+pcp(2,2)**2+pcp(3,2)**2)
        call jamerrm(3,0,'(jamscat2:)energy or momentum not conserved')
        mste(1)=-99
        return
      endif

      end

c***********************************************************************
      function eminres(kc1)
      implicit none
      include 'jam2.inc'
      integer id,kc1
      real*8 eminres

      id=kchg(kc1,5)
      if(id.eq.id_delt.or.id.eq.id_nucls.or.id.eq.id_delts) then
        eminres = 1.08d0
      else
        eminres=pmas(kc1,1)-pmas(kc1,3)
      endif

      end

c***********************************************************************

      subroutine jamscat2m(i1,i2,pcm,kf1,kf2,kc1,kc2,em1,em2,
     $ ibar1,ibar2,pr0,pr,taucol)

c...Calculate kinematics for the two-body final case.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      dimension pcm(5)
      logical jamrqpb2

      sig=pare(4)
      srt=pcm(5)
c...em3, em4 are the out going mass.
      em1ef=em1
      em2ef=em2

c...p^- and p^+ before collision.
      em1sq=pcp(4,1)**2-pcp(1,1)**2-pcp(2,1)**2-pcp(3,1)**2
      em2sq=pcp(4,2)**2-pcp(1,2)**2-pcp(2,2)**2-pcp(3,2)**2 
      e01=sqrt(em1sq+pr0*pr0)
      e02=sqrt(em2sq+pr0*pr0)
      pm01=sqrt(em1sq+pr0*pr0)-pr0
      pp02=sqrt(em2sq+pr0*pr0)-pr0

c...Scattering and azimuthal angles in two-body cm-fame.
      if(mste(1).eq.1) then

          if(em1sq.ge.0d0) then
            em1ef=sqrt(em1sq)
          else
            print *,'jamscat2m em3<0? ',em3sq,kf1,p(5,i1),pare(21)
            em1ef=0.001d0
          endif
          if(em2sq.ge.0d0) then
            em2ef=sqrt(em2sq)
          else
            print *,'jamscat2m em4<0? ',em4sq,kf2,p(5,i2),pare(22)
            em2ef=0.001d0
          endif

c.......Get elastic angular distribution.
          call jamangel(pr,sig,kf1,kf2,ibar1,ibar2,t1,c1) 

      else

        em1ef = em1*pare(21)
        em2ef = em2*pare(22)

c.....Compute magnitude of rel. momentum using the new masses.
        pr2=(srt**2-(em1ef+em2ef)**2)*(srt**2-(em1ef-em2ef)**2)
        if(pr2.gt.0.00001d0) then
          pr=sqrt(pr2)/(2.d0*srt)
        else
          write(check(1),'(g15.3,1x,i9,1x,i9,g10.3,1x,g10.3)')
     $       srt,kf1,kf2,em1ef,em2ef
          call jamerrm(3,1,'(jamscat2m:) rel.mom.pr<0 after scatter')
c         print *,'mste(1)=',mste(1)
c         print *,'kf1 kf2',kf1,kf2
c         print *,'em1ef em2ef',em1ef,em2ef
c         print *,'em1 em2',em1,em2,pare(21),pare(22)

          mste(1)=-99
          return
        endif


c.....Get inelastic angular distribution.
        snew=sqrt(em1**2+pr2)+sqrt(em2**2+pr2)
        call jamangin(snew,pr,pr0,em1,em2,ibar1,ibar2,
     $                kcp(2,1),kcp(2,2),mste(22),mste(24)
     $                  ,kc1,kc2,t1,c1)

      end if


      betx=pare(12)
      bety=pare(13)
      betz=pare(14)
      gam=pare(15)
      phi=pare(16)
      the=pare(17)
      pare(18)=t1
      pare(19)=c1

c.....Rotate momentum vector to new position.
      sin1=sqrt(1.d0-c1**2)
      pxr=pr*sin1*cos(t1)
      pyr=pr*sin1*sin(t1)
      pzr=pr*c1
      pare(28)=pxr
      pare(29)=pyr
      pare(30)=pzr-pr0
      pare(20)=pxr*pare(25)+pyr*pare(26)+(pzr-pr0)*pare(27)


c...Estimate collision time scale from uncertainty principle.
c     pt=sqrt(pxr**2+pyr**2)
c     if(pt.gt.1d-7) then
c       taucol=gam*paru(3)/pt
c       taucol=-tau*log(max(rn(0),1.d-35))
c       if(taucol.gt.1.0d0) then
c         write(33,*)mste(1),srt,taucol
c         mste(1)=-9
c         return
c       endif
c     endif

c....EoS modified collision.
      if(mstc(59).ge.101) then
         call jameoscol(i1,i2,pr0,pxr,pyr,pzr,em1ef,em2ef)
      else if(mstc(59).ge.1) then
        pt=pr*sin1
        call jameosscat(i1,i2,pxr,pyr,pyz,pt,pr)
      endif
 
      pr2=pxr**2+pyr**2+pzr**2
      per=sqrt(em1ef**2+pr2)
      pxr2=-pxr
      pyr2=-pyr
      pzr2=-pzr
      per2=sqrt(em2ef**2+pr2)

c...p^- after the collision.
      pm1=per-pzr
      pp2=per2+pzr2

c...collision time in the rest frame of particle.
      if(mstc(127).eq.1) then
        dtau1=paru(3)/max(1e-8,abs(pm01-pm1))/pare(56)
        dtau2=paru(3)/max(1e-8,abs(pp02-pp2))/pare(57)

      else if(mstc(127).eq.2) then
        pt=sqrt(pxr**2+pyr**2)
        dtau1=0.5*paru(3)/pt/pare(56)
        dtau2=0.5*paru(3)/pt/pare(57)
      else if(mstc(127).eq.3) then
        dtau1=0.5*paru(3)/max(1e-8,abs(e01-per))/pare(56)
        dtau2=0.5*paru(3)/max(1e-8,abs(e02-per2))/pare(57)
      endif

c...Back to the original frame.
      call jamrobo(the,phi,betx,bety,betz,gam,pxr,pyr,pzr,per)
      p(1,i1)=pxr
      p(2,i1)=pyr
      p(3,i1)=pzr
      p(4,i1)=sqrt(em1ef**2+p(1,i1)**2+p(2,i1)**2+p(3,i1)**2)
      p(5,i1)=em1
      per=p(4,i1)
 
      call jamrobo(the,phi,betx,bety,betz,gam,pxr2,pyr2,pzr2,per2)
      p(1,i2)=pxr2
      p(2,i2)=pyr2
      p(3,i2)=pzr2
      p(4,i2)=sqrt(em2ef**2+p(1,i2)**2+p(2,i2)**2+p(3,i2)**2)
      p(5,i2)=em2
      per2=p(4,i2)

c...Quantum mechanical suppression of very soft multiple collision.
      if(mstc(127).ge.1) then
        if(dtau1.gt.100.0.or.dtau2.gt.100.0) then
          mste(1)=-9
          return
        endif
        if(jamcoltim(i1,i2,dtau1,dtau2).eq.0) then
          mste(1)=-9
          return
        endif
      endif

c...2020/2/4
      if(mstd(91).eq.1) then
        if(MF_on(i1).eq.1) then
        p(1,i1)=p(1,i1)+potv(1,i1)
        p(2,i1)=p(2,i1)+potv(2,i1)
        p(3,i1)=p(3,i1)+potv(3,i1)
        p(4,i1)=sqrt(em1ef**2+p(1,i1)**2+p(2,i1)**2+p(3,i1)**2)
        endif
        if(MF_on(i2).eq.1) then
        p(1,i2)=p(1,i2)+potv(1,i2)
        p(2,i2)=p(2,i2)+potv(2,i2)
        p(3,i2)=p(3,i2)+potv(3,i2)
        p(4,i2)=sqrt(em2ef**2+p(1,i2)**2+p(2,i2)**2+p(3,i2)**2)
        endif
      endif

c     emf1=p(4,i1)**2-p(1,i1)**2-p(2,i1)**2-p(3,i1)**2
c     emf2=p(4,i2)**2-p(1,i2)**2-p(2,i2)**2-p(3,i2)**2
c     if(emf1.le.0.0.or.emf2.le.0.0) then
c     print *,'jamscat1m mste(1)=',mste(1)
c     print *,'jamscat1m emf1=',emf1,em1ef,kf1
c     print *,'jamscat1m emf2=',emf2,em2ef,kf2
c     endif

c     pare(37)=pxr
c     pare(38)=pyr
c     pare(39)=pzr

c     pare(40)=pxr2
c     pare(41)=pyr2
c     pare(42)=pzr2

c     dpx=dpx+p(1,i1)
c     dpy=dpy+p(2,i1)
c     dpz=dpz+p(3,i1)
c     dpe=dpe+p(4,i1)
c     dx=r(1,i1)-r(1,i2)
c     dy=r(2,i1)-r(2,i2)
c     dz=r(3,i1)-r(3,i2)
c     dt=r(4,i1)-r(4,i2)
c     dxp=dpx*dx + dpy*dy + dpz*dz  - dpe*dt
c     tcm1=(pcm(4)*r(4,i1)
c    &   - pcm(1)*r(1,i1) - pcm(2)*r(2,i1) - pcm(3)*r(3,i1))/pcm(5)
c     tcm2=(pcm(4)*r(4,i2)
c    &   - pcm(1)*r(1,i2) - pcm(2)*r(2,i2) - pcm(3)*r(3,i2))/pcm(5)

c     tcm=(pcm(4)*dt - pcm(1)*dx - pcm(2)*dy - pcm(3)*dz)/pcm(5)
c     ecm=(pcm(4)*dpe - pcm(1)*dpx - pcm(2)*dpy - pcm(3)*dpz)/pcm(5)
c     dxp=tcm*ecm+dxp

c     print *,'kf1 kf2=',k(2,i1),k(2,i2)
c     print *,'t1 t2=',r(4,i1),r(4,i2)
c     print *,'tcm1 tcm2=',tcm1,tcm2
c     print *,'dxp=',dxp,pare(20)
c     print *,'r10=',pare(43),pare(44),pare(45)
c     print *,'r11=',(r(j,i1),j=1,3)
c     print *,'r20=',pare(46),pare(47),pare(48)
c     print *,'r21=',(r(j,i2),j=1,3)
c     print *,'dxcm=',pare(25),pare(26),pare(27)


c         if(abs(em1ef-em1).gt.1d-4.or.abs(em2ef-em2).gt.1d-4) then
c            print *,'mste(1)',mste(1)
c            print *,'em1ef em1',kf1,em1ef,em1,pare(21)
c            print *,'em2sq em2',kf2,em2ef,em2,pare(22)
c         endif

c....Check energy momentum conservation.
c     srt1=sqrt((p(4,i1)+p(4,i2))**2-(p(1,i1)+p(1,i2))**2
c    $               -(p(2,i1)+p(2,i2))**2-(p(3,i1)+p(3,i2))**2)
      srt1=sqrt((per+per2)**2-(pxr+pxr2)**2
     $               -(pyr+pyr2)**2-(pzr+pzr2)**2)

      if( (abs(pcm(1)-(pxr+pxr2)).ge.1d-5)
     a    .or. (abs(pcm(2)-(pyr+pyr2)).ge.1d-5)
     a      .or. (abs(pcm(3)-(pzr+pzr2)).ge.1d-5)
     a      .or. (abs(pcm(4)-(per+per2)).ge.1d-5)
     $        .or. (abs((srt1-srt)/srt).gt.0.1d0) ) then

        ih=mstc(38)
        write(ih,'(/,''<<jamscat2m>>'')')
        write(ih,*)'beta gam',betx,bety,betz,gam
        write(ih,*)'ichanel srt sig=',
     $    mste(1),srt,sig,pcm(5)
        write(ih,*)'pot1 pot2',pare(21),pare(22)
        write(ih,*)'c1 t1',c1,t1
        write(ih,*)'pr0',pr0
        write(ih,*)'pr pxr pyr pzr',pr,pxr,pyr,pzr
        write(ih,*)'kf1 em1',kcp(2,1),(pcp(i,1),i=5,5)
        write(ih,*)'kf2 em2',kcp(2,2),(pcp(i,2),i=5,5)
        write(ih,*)'kf3 em3',kf1,em1ef,(p(i,i1),i=5,5)
        write(ih,*)'kf4 em4',kf2,em2ef,(p(i,i2),i=5,5)
        write(ih,*)'pcm1',pcm(1),p(1,i1)+p(1,i2)
        write(ih,*)'pcm2',pcm(2),p(2,i1)+p(2,i2)
        write(ih,*)'pcm3',pcm(3),p(3,i1)+p(3,i2)
        write(ih,*)'e2  ',pcm(4),p(4,i1)+p(4,i2)
        write(ih,*)'srt ',srt,srt1,abs(srt1-srt)/srt

        call jamerrm(3,0,'(jamscat2m:)energy or mom. not conserved')
        mste(1)=-99
        return
      endif

      end

c***********************************************************************

      integer function jamcoltim(i1,i2,dtau1,dtau2)

c...Purpose: to evaluate suppression of multiple collisions.

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      integer i1,i2,i,iq1,iq2
      real*8 dtau1,dtau2,gam1,gam2,tc1,tc2

      jamcoltim=1
      iq1=0
      iq2=0
      if(r(4,i1).lt.r(5,i1)) iq1=1
      if(r(4,i2).lt.r(5,i2)) iq2=1
c     if(iq1.eq.1.or.iq2.eq.1) return

      gam1=p(4,i1)/sqrt(p(4,i1)**2+p(1,i1)**2+p(2,i1)**2+p(3,i1)**2)
      gam2=p(4,i2)/sqrt(p(4,i2)**2+p(1,i2)**2+p(2,i2)**2+p(3,i2)**2)
      tc1 = pare(8) + parc(127)*dtau1*gam1
      tc2 = pare(9) + parc(127)*dtau2*gam2

      do i=1,mentry
        if(icoll(1,i).eq.-1) cycle

        if(iq1.eq.0) then
        if( (icoll(1,i).eq.i1.and.icoll(2,i).ne.i2)
     &  .or.(icoll(2,i).eq.i1.and.icoll(1,i).ne.i2)) then

          if(coll(2,i) .lt.tc1) then
            jamcoltim=0
            return
          endif

        endif
        endif

        if(iq2.eq.0) then
        if((icoll(1,i).eq.i2.and.icoll(2,i).ne.i1)
     &  .or.    (icoll(2,i).eq.i2.and.icoll(1,i).ne.i1)) then

          if(coll(3,i) .lt.tc2) then
            jamcoltim=0
            return
          endif

        endif
        endif

      end do

      end

c***********************************************************************

      subroutine jamscatd2(i1,i2,pcm,kf1,kf2,kc1,kc2,em1,em2,kd,pd)

c...Handle NN -> deuteron + pion + pion scattering.
c...Ref. arXiv:1212.2881[nucl-ex]

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      dimension pcm(5),pms(5),p1(5),p2(5)
      dimension pd(10,5),kd(10,2),pde(5)
      logical anti
     
      if(kf1.le.0.and.kf2.le.0) then
        anti=.true.
      else
        anti=.false.
      endif

      do i=1,10
       do j=1,5
       pd(i,j)=0.0d0
       end do
       kd(i,1)=0
       kd(i,2)=0
      end do

c...First create deuteron* and pion.
      emp=0.93827
      emn=0.93957
      empip=0.13957   ! pi+/-
      empi0=0.13498   ! pi0
      emd=emp+emn+0.0001
      if(mstc(78).eq.4) then
        emd=emp+emn-0.0022   ! Form real deuteron
        kd(1,2)=1000000*1+1000*1+1000000000
        k(2,i1)=kd(1,2) 
        k(1,i1)=5
        k(3,i1)=0
        k(4,i1)=0
        if(anti) k(2,i1)=-k(2,i1)
        k(1,i2)=21
        mste(27)=0  ! i2
        kf1=k(2,i1)
        kc1=0
        kf2=0
        kc2=0
      endif

c....pp -> d pi+ pi0  or nn -> d pi- pi0
      if(mste(3).eq.14.or.mste(3).eq.15) then
        nd=3
        pd(1,5)=emd
        pd(2,5)=empip
        pd(3,5)=empi0
        kd(2,2)=211
        if(mste(3).eq.15) kd(2,2)=-211
        kd(3,2)=111
        call jamscatp(nd,pcm,pd)
        pde(1)=pd(1,1)
        pde(2)=pd(1,2)
        pde(3)=pd(1,3)
        pde(4)=sqrt(emd**2+pde(1)**2+pde(2)**2+pde(3)**2)
        pde(5)=emd

        nv=nv+2
        nmeson=nmeson+2
        if(nv.gt.mxv) call jamerrm(30,0,'(jamscatm:)particle too large')
        i3=nv-1
        i4=nv
        mste(29)=i3
        mste(31)=i4
        do j=1,5
         p(j,i3)=pd(2,j)
         p(j,i4)=pd(3,j)
         r(j,i3)=0.5d0*(p(j,i1)+p(j,i2))
         r(j,i4)=0.5d0*(p(j,i1)+p(j,i2))
        end do

        if(anti)then
          kd(1,2)=-kd(1,2)
          if(kchg(kd(2,1),3).ne.0)kd(2,2)=-kd(2,2)
          if(kchg(kd(3,1),3).ne.0)kd(3,2)=-kd(3,2)
        endif
        k(2,i3)=kd(2,2)
        k(2,i4)=kd(3,2)
        kd(4,2)=kd(3,2)
        kd(3,2)=kd(2,2)
        kd(3,1)=jamcomp(kd(3,2))
        kd(4,1)=jamcomp(kd(4,2))

c... pn->d pi0 pi0   or  pn->d pi+ pi-
      else if(mste(3).eq.16.or.mste(3).eq.17) then

c...first d* -> d sigma' decay
      emsigm=0.3d0
      emr=0.3d0     ! mass of sigma-meson like particle
      wid=0.1d0     ! width of sigma-meson like particle
      emmin=empip+empip+parc(41)
      emmax=pcm(5)-emd-0.001d0
      call jambwmas(emmin,emmax,emr,wid,emsigm,icon)
      pde(5)=emd
      pms(5)=emsigm
      call jamn2dec(pcm,pde,pms)
      nv=nv+1
      nmeson=nmeson+1
      if(nv.gt.mxv) call jamerrm(30,0,'(jamscatd2:)particle too large')
      i3=nv
      mste(29)=i3
      do j=1,5
        p(j,i3)=pms(j)
      end do
      k(2,i3)=10220
      r(1,i3)=(r(1,i1)+r(1,i2))/2.d0
      r(2,i3)=(r(2,i1)+r(2,i2))/2.d0
      r(3,i3)=(r(3,i1)+r(3,i2))/2.d0
      kd(3,2)=10220
      kd(3,1)=jamcomp(kd(3,2))
      pd(3,5)=p2(5)

c....now sigma' -> pi pi decay
      if(mste(3).eq.16) then
        kd(3,2)=111
        kd(4,2)=111
        pd(3,5)=empi0
        pd(4,5)=empi0
      else if(mste(3).eq.17) then
        kd(3,2)=211
        kd(4,2)=-211
        pd(3,5)=empip
        pd(4,5)=empip
      endif
      kd(3,1)=jamcomp(kd(3,2))
      kd(4,1)=jamcomp(kd(4,2))
      nv=nv+1
      nmeson=nmeson+1
      if(nv.gt.mxv) call jamerrm(30,0,'(jamscatd2:) particle too large')
      i4=nv
      mste(31)=i4
      p1(5)=pd(3,5)
      p2(5)=pd(4,5)
      call jamn2dec(pms,p1,p2)
      do j=1,5
        p(j,i3)=p1(j)
        p(j,i4)=p2(j)
      end do
      k(2,i3)=kd(3,2)
      k(2,i4)=kd(4,2)
      do j=1,3
        r(j,i3)=(r(j,i1)+r(j,i2))/2.d0
        r(j,i4)=(r(j,i1)+r(j,i2))/2.d0
      end do

      else
        print *,'jamscatd2 wrong mste(3)',mste(3)
        stop
      endif

c....Propagate deuteron.
      if(mstc(78).eq.4) then
        do j=1,5
          p(j,i1)=pde(j)
        end do
        return
      endif

c...finally deuteron decay into p and n.
      p1(5)=emp
      p2(5)=emn
      call jamn2dec(pde,p1,p2)
      do j=1,5
        p(j,i1)=p1(j)
        p(j,i2)=p2(j)
      end do
      k(2,i1)=2212
      k(2,i2)=2112
      kd(1,2)=2212
      kd(2,2)=2112

      if(anti)then
        kd(1,2)=-kd(1,2)
        kd(2,2)=-kd(2,2)
      endif

      kd(1,1)=jamcomp(kd(1,2))
      kd(2,1)=jamcomp(kd(2,2))
      pd(1,5)=emp
      pd(2,5)=emn
      em1=emp
      em2=emn
      kf1=k(2,i1)
      kf2=k(2,i2)
      kc1=kd(1,1)
      kc2=kd(2,1)

      end

c***********************************************************************

      subroutine jamn2dec(pcm,p1,p2)

c...Purpose: to compute 2-body decay.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension pcm(5),p1(5),p2(5)

      s=pcm(5)*pcm(5)
      pr=sqrt((s-(p1(5)+p2(5))**2)*(s-(p1(5)-p2(5))**2))/(2*pcm(5))
      cos1=1.0d0-2*rn(0)
      sin1=sqrt(1.0d0-cos1**2)
      phi1=paru(2)*rn(0)
      prx=pr*sin1*cos(phi1)
      pry=pr*sin1*sin(phi1)
      prz=pr*cos1
      pre=sqrt(p1(5)**2+pr**2)
      px2=-prx
      py2=-pry
      pz2=-prz
      pe2=sqrt(p2(5)**2+pr**2)
      betx=pcm(1)/pcm(4)
      bety=pcm(2)/pcm(4)
      betz=pcm(3)/pcm(4)
      gamm=pcm(4)/pcm(5)
      call jamrobo(0d0,0d0,betx,bety,betz,gamm,prx,pry,prz,pre)
      call jamrobo(0d0,0d0,betx,bety,betz,gamm,px2,py2,pz2,pe2)
      p1(1)=prx
      p1(2)=pry
      p1(3)=prz
      p1(4)=sqrt(p1(5)**2+p1(1)**2+p1(2)**2+p1(3)**2)
      p2(1)=px2
      p2(2)=py2
      p2(3)=pz2
      p2(4)=sqrt(p2(5)**2+p2(1)**2+p2(2)**2+p2(3)**2)

      end

c***********************************************************************

      subroutine jamscatd(i1,i2,pcm,kf1,kf2,kc1,kc2,em1,em2,pr0,kd,pd)

c...Handle NN -> deuteron pion scattering.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      dimension pcm(5)
      dimension pd(10,5),kd(10,2),pde(5)
      logical anti
     
      if(kf1.le.0.and.kf2.le.0) then
        anti=.true.
      else
        anti=.false.
      endif

      do i=1,10
       do j=1,5
       pd(i,j)=0.0d0
       end do
       kd(i,1)=0
       kd(i,2)=0
      end do

c...First create deuteron(*) and pion.
      emp=0.93827
      emn=0.93957
      emd=emp+emn+0.0001
      if(mstc(77).eq.4) then
        emd=emp+emn-0.0022   ! Form real deuteron
        kd(1,2)=1000000*1+1000*1+1000000000
        k(2,i1)=kd(1,2) 
        k(1,i1)=5
        k(3,i1)=0
        k(4,i1)=0
        if(anti) k(2,i1)=-k(2,i1)
        k(1,i2)=21
        mste(27)=0  ! i2
        kf1=k(2,i1)
        kc1=0
        kf2=0
        kc2=0
      endif
      if(mste(3).eq.11.or.mste(3).eq.12) then
        empi=0.13957   ! pi+/-
        kd(3,2)=211
        if(mste(3).eq.12) kd(3,2)=-211
      else if(mste(3).eq.13) then
        empi=0.13498   ! pi0
        kd(3,2)=111
      else
        print *,'wrong mste(3)',mste(3)
        stop
      endif

      if(anti)then
        if(kchg(kd(3,1),3).ne.0) kd(3,2)=-kd(3,2)
      endif
      kd(3,1)=jamcomp(kd(3,2))
      pd(3,5)=empi

      betx=pare(12)
      bety=pare(13)
      betz=pare(14)
      gam=pare(15)
      phi=pare(16)
      the=pare(17)
      srt=pcm(5)
      s=srt*srt
      pr=sqrt((s-(emd+empi)**2)*(s-(emd-empi)**2))/(2*srt)

c.....Rotate momentum vector to new position.
      t1=2.0d0*paru(1)*rn(0)
      call jamangrr(pr0,pr,srt,c1,emd,empi,mstc(67)) 
      sin1=sqrt(1.d0-c1**2)
      pxr=pr*sin1*cos(t1)
      pyr=pr*sin1*sin(t1)
      pzr=pr*c1
      per=sqrt(emd**2+pr**2)
      pxr2=-pxr
      pyr2=-pyr
      pzr2=-pzr
      per2=sqrt(empi**2+pr**2)

c...Back to the original frame.
      call jamrobo(the,phi,betx,bety,betz,gam,pxr,pyr,pzr,per)
      p(1,i1)=pxr
      p(2,i1)=pyr
      p(3,i1)=pzr
      p(4,i1)=sqrt(emd**2+p(1,i1)**2+p(2,i1)**2+p(3,i1)**2)
      p(5,i1)=emd
 
      call jamrobo(the,phi,betx,bety,betz,gam,pxr2,pyr2,pzr2,per2)
      nv=nv+1
      nmeson=nmeson+1
      if(nv.gt.mxv) call jamerrm(30,0,'(jamscatd:)particle too large')
      i3=nv
      p(1,i3)=pxr2
      p(2,i3)=pyr2
      p(3,i3)=pzr2
      p(4,i3)=sqrt(empi**2+p(1,i3)**2+p(2,i3)**2+p(3,i3)**2)
      p(5,i3)=empi
      mste(29)=i3
      r(1,i3)=(r(1,i1)+r(1,i2))/2.d0
      r(2,i3)=(r(2,i1)+r(2,i2))/2.d0
      r(3,i3)=(r(3,i1)+r(3,i2))/2.d0

      if(mstc(77).eq.4) return
   
c...deuteron decay into p and n, (p,n) pair will be deuteron-like
      do j=1,5
        pde(j)=p(j,i1)
      end do
      prsq=(emd*emd-(emp+emn)**2)*(emd*emd-(emp-emn)**2)
      pr=sqrt(prsq)/(2*emd)
      cos1=1.0d0-2*rn(0)
      sin1=sqrt(1.0d0-cos1**2)
      phi1=paru(2)*rn(0)
      prx=pr*sin1*cos(phi1)
      pry=pr*sin1*sin(phi1)
      prz=pr*cos1

      pcs=pde(1)*prx+pde(2)*pry+pde(3)*prz
      ecm1=sqrt(emp**2+prx**2+pry**2+prz**2)
      transf=(pcs/(pde(4)+pde(5))+ecm1)/pde(5)
      p(1,i1)=prx+pde(1)*transf
      p(2,i1)=pry+pde(2)*transf
      p(3,i1)=prz+pde(3)*transf
      p(4,i1)=sqrt(emp**2+p(1,i1)**2+p(2,i1)**2+p(3,i1)**2)
      p(5,i1)=emp
      k(2,i1)=2212
      kf1=2212
      em1=emp
      kd(1,2)=2212

      ecm2=sqrt(emn**2+prx**2+pry**2+prz**2)
      transf=(-pcs/(pde(4)+pde(5))+ecm2)/pde(5)
      p(1,i2)=-prx+pde(1)*transf
      p(2,i2)=-pry+pde(2)*transf
      p(3,i2)=-prz+pde(3)*transf
      p(4,i2)=sqrt(emn**2+p(1,i2)**2+p(2,i2)**2+p(3,i2)**2)
      p(5,i2)=emn
      k(2,i2)=2112
      kf2=2112
      em2=emn
      kd(2,2)=2112

      kd(1,1)=jamcomp(kd(1,2))
      kd(2,1)=jamcomp(kd(2,2))
      pd(1,5)=emp
      pd(2,5)=emn

      if(anti)then
        kd(1,2)=-kd(1,2)
        kd(2,2)=-kd(2,2)
      endif
      kf1=kd(1,2)
      kf2=kd(2,2)
      kc1=kd(1,1)
      kc2=kd(2,1)

      end

c***********************************************************************

      subroutine jamscatm(i1,i2,pcm,kf1,kf2,kc1,kc2,em1,em2,kd,pd)

c...Handle multi-particle final scattering.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      dimension pcm(5)
      dimension pd(10,5),kd(10,2),pr(3)
      logical anti
     
      if(kf1.le.0.and.kf2.le.0) then
        anti=.true.
      else
        anti=.false.
      endif

      do i=1,10
       do j=1,5
       pd(i,j)=0.0d0
       end do
       kd(i,1)=0
       kd(i,2)=0
      end do

c2014/11/6 bug fix  pointed out by A.Ono
      pr(1)=p(1,i1)
      pr(2)=p(2,i1)
      pr(3)=p(3,i1)

c....s-wave pion production.
      if(mste(3).le.3) then
        nd=3
        pd(1,5)=p(5,i1)
        pd(2,5)=p(5,i2)
        if(mste(3).eq.1) then ! pp
          if(rn(0).le.0.3333333d0) then
            kd(1,2)=2212
            kd(2,2)=2212
            kd(3,2)=111
          else
            kd(1,2)=2212
            kd(2,2)=2112
            kd(3,2)=211
          endif
        else if(mste(3).eq.2) then ! nn
          if(rn(0).le.0.3333333d0) then
            kd(1,2)=2112
            kd(2,2)=2112
            kd(3,2)=111
          else
            kd(1,2)=2212
            kd(2,2)=2112
            kd(3,2)=-211
          endif
        else ! pn
         if(rn(0).le.0.5d0) then  ! np->np*
c          if(rn(0).le.0.666666d0) then ! np->n(p pi0)
c   bug fix pointed by A.Ono 2015/2/5
           if(rn(0).le.0.333333d0) then ! np->n(p pi0) 
             kd(1,2)=2112
             kd(2,2)=2212
             kd(3,2)=111
           else              ! np->n(n pi+)
             kd(1,2)=2112
             kd(2,2)=2112
             kd(3,2)=211
           endif
         else  ! np->pn*
           if(rn(0).le.0.333333d0) then ! np->p(n pi0)
             kd(1,2)=2212
             kd(2,2)=2112
             kd(3,2)=111
           else               ! np->p(p pi-)
             kd(1,2)=2212
             kd(2,2)=2212
             kd(3,2)=-211
           endif
         endif
        endif
        kd(1,1)=jamcomp(kd(1,2))
        kd(2,1)=jamcomp(kd(2,2))
        kd(3,1)=jamcomp(kd(3,2))
        pd(1,5)=pjmass(kd(1,2))
        pd(2,5)=pjmass(kd(2,2))
        pd(3,5)=pjmass(kd(3,2))
      endif

c...Kinematics.
      if(mste(3).le.3) then
        if(pcm(5).le.2.9d0) then
          call jamscatp(nd,pcm,pd)
        else
          call jamscat3(nd,pcm,pd,pr)
        endif
      else
        call jamscatp(nd,pcm,pd)
      endif

      if(mste(3).le.3) then
        nv=nv+1
        nmeson=nmeson+1
        if(nv.gt.mxv) call jamerrm(30,0,'(jamscatm:)particle too large')
        i3=nv
        do j=1,5
         p(j,i1)=pd(1,j)
         p(j,i2)=pd(2,j)
         p(j,i3)=pd(3,j)
         r(j,i3)=0.5d0*(p(j,i1)+p(j,i2))
        end do

        if(anti)then
          kd(1,2)=-kd(1,2)
          kd(2,2)=-kd(2,2)
          if(kchg(kd(3,1),3).ne.0)kd(3,2)=-kd(3,2)
        endif

      endif

      kf1=kd(1,2)
      kf2=kd(2,2)
      kf3=kd(3,2)

      kf4=kd(4,2)
      kc1=kd(1,1)

      kc2=kd(2,1)
      em1=pd(1,5)
      em2=pd(2,5)

      mste(29)=i3
      k(2,i3)=kd(3,2)
      p(5,i3)=pd(3,5)
 
      end

c***********************************************************************

      subroutine jamscat3(nd,pcm,pd,pr)

c...Kinematics for three-body final scattering.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension pcm(5),pd(10,5),pr(3)

      if(nd.ne.3) call jamerrm(30,0,'(jamscat3:) invalid nd')

      shp=pcm(5)**2
      shpr=sqrt(shp)
      betx=pcm(1)/pcm(4)
      bety=pcm(2)/pcm(4)
      betz=pcm(3)/pcm(4)
      gam=pcm(4)/pcm(5)

      prr=sqrt(pr(1)**2+pr(2)**2+pr(3)**2)
      pre=sqrt(pd(1,5)**2+prr**2)
      call jamrobo(0.0d0,0.0d0,-betx,-bety,-betz,gam,
     $ pr(1),pr(2),pr(3),pre)
      phi=pjangl(pr(1),pr(2))
      call jamrobo(0.0d0,-phi,0.0d0,0.0d0,0.0d0,1.0d0,pr(1),pr(2),pr(3),
     & pre)
      the=pjangl(pr(3),pr(1))
      call jamrobo(-the,0.d0,0.d0,0.d0,0.d0,1.0d0,pr(1),pr(2),pr(3),pre)

 100  do i=1,2
        pt=min(prr-0.001d0,0.4d0*sqrt(-log(max(1.d0-10,rn(0)))))
        phi1=2*paru(1)*rn(0)
        pd(i,1)=pt*cos(phi1)
        pd(i,2)=pt*sin(phi1)
      end do

      pd(3,1)=-pd(1,1)-pd(2,1)
      pd(3,2)=-pd(1,2)-pd(2,2)
      pms1=pd(1,5)**2+pd(1,1)**2+pd(1,2)**2
      pms2=pd(2,5)**2+pd(2,1)**2+pd(2,2)**2
      pms3=pd(3,5)**2+pd(3,1)**2+pd(3,2)**2
      pmt3=sqrt(pms3)
      pm12=(sqrt(pms1)+sqrt(pms2))**2


C...Select rapidity for particle 3 and check phase space not closed.
      y3max=log((shp+pms3-pm12+sqrt(max(0.d0,(shp-pms3-pm12)**2-
     &  4.d0*pms3*pm12)))/(2.d0*shpr*pmt3))
      if(y3max.lt.1d-6) then
        goto 100
      endif
      y3=(2.d0*rn(0)-1.d0)*0.999999d0*y3max

      eps=1.d0
      if(rn(0).gt.0.5d0) eps=-1.d0

      pd(3,3)=pmt3*sinh(y3)
      pd(3,4)=pmt3*cosh(y3)
      pms12=(pcm(5)-pd(3,4))**2-pd(3,3)**2
      sql12=(pms12-pms1-pms2)**2-4.d0*pms1*pms2
      if(sql12.le.0.d0) then
         goto 100
      endif
      pd(1,3)=(-pd(3,3)*(pms12+pms1-pms2)+
     &  eps*(pcm(5)-pd(3,4))*sqrt(sql12))/(2.d0*pms12)
      pd(2,3)=-pd(1,3)-pd(3,3)
      pd(1,4)=sqrt(pms1+pd(1,3)**2)
      pd(2,4)=sqrt(pms2+pd(2,3)**2)


      end

c***********************************************************************

      subroutine jamscatp(nd,pcm,pd)

c...Scatter the particle according to phase space.

      implicit double precision(a-h, o-z)
      dimension pcm(5)
      dimension pv(10,5),pd(10,5),rord(10),ue(3),be(3),wtcor(10)
      data wtcor/2.d0,5.d0,15.d0,60.d0,250.d0,1500.d0,1.2d4,1.2d5, 
     & 150.d0,16.d0/ 

C...Functions: momentum in two-particle decays.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a)

c...Calculate maximum weight nd-particle decay. 
      do j=1,5
       pv(1,j)=pcm(j)
      end do
      ps=0.0d0
      do i=1,nd
        ps=ps+pd(i,5)
      end do

      pv(nd,5)=pd(nd,5) 
      wtmax=1.d0/wtcor(nd-2) 
      pmax=pv(1,5)-ps+pd(nd,5) 
      pmin=0.d0 
      do il=nd-1,1,-1 
        pmax=pmax+pd(il,5) 
        pmin=pmin+pd(il+1,5) 
        wtmax=wtmax*pawt(pmax,pmin,pd(il,5)) 
      end do

c...M-generator gives weight. If rejected, try again. 
  440 rord(1)=1.d0 
      do 470 il1=2,nd-1 
        rsav=rn(0) 
        do 450 il2=il1-1,1,-1 
          if(rsav.le.rord(il2)) goto 460 
          rord(il2+1)=rord(il2) 
  450   continue 
  460   rord(il2+1)=rsav 
  470 continue 
        rord(nd)=0.d0 
        wt=1.d0 
        do 480 il=nd-1,1,-1 
        pv(il,5)=pv(il+1,5)+pd(il,5)+(rord(il)-rord(il+1))*(pv(1,5)-ps) 
        wt=wt*pawt(pv(il,5),pv(il+1,5),pd(il,5)) 
  480   continue 
        if(wt.lt.rn(0)*wtmax) goto 440 

c...Perform two-particle decays in respective CM frame. 
      do 510 il=1,nd-1 

      pa=pawt(pv(il,5),pv(il+1,5),pd(il,5)) 
      ue(3)=2.d0*rn(0)-1.d0 
      phi=2*3.14159d0*rn(0) 
      ue(1)=sqrt(1.d0-ue(3)**2)*cos(phi) 
      ue(2)=sqrt(1.d0-ue(3)**2)*sin(phi) 

      do j=1,3 
      pd(il,j)=pa*ue(j) 
      pv(il+1,j)=-pa*ue(j) 
      end do
      pd(il,4)=sqrt(pa**2+pd(il,5)**2) 
      pv(il+1,4)=sqrt(pa**2+pv(il+1,5)**2) 

  510 continue 
 
c...Lorentz transform decay products to lab frame. 
      do j=1,4 
      pd(nd,j)=pv(nd,j) 
      end do

      do 560 il=nd-1,1,-1 

        do j=1,3 
        be(j)=pv(il,j)/pv(il,4) 
        end do 
        ga=pv(il,4)/pv(il,5) 

        do i=il,nd 
          bep=be(1)*pd(i,1)+be(2)*pd(i,2)+be(3)*pd(i,3) 
          do j=1,3 
          pd(i,j)=pd(i,j)+ga*(ga*bep/(1.d0+ga)+pd(i,4))*be(j) 
          end do
          pd(i,4)=sqrt(pd(i,5)**2+pd(i,1)**2+pd(i,2)**2+pd(i,3)**2)
        end do

  560 continue 

      end

c***********************************************************************

      subroutine jamkupda(msel,i1,kf,kc,ks,k6,icltag,ichanel)

c...Purpose: to update particle status code k() after collision.
c...msel=1: collision
c...    =2: resonance decay
c...    =3: string decay
c...    =7: partons
c....Particle arrays updated here are:
c        k(1,) : status code
c        k(2,) : flavor code
c        k(3,) : origin of production.
c        k(4,) : where this particle comes from  (old:spin)
c        k(5,) : collision counter
c        k(6,) : multi-step history
c        k(7,) : number of collisions sufferd so far
c        k(8,) : test particle identifer
c        k(9,) : baryon number
c        k(10,): color flow connection
c        k(11,): color flow connection

      implicit double precision(a-h, o-z)
      include 'jam1.inc' 
      include 'jam2.inc'

c....Dead particle.
      if(i1.eq.0) return
      if(ks.eq.0) then
        k7=k(7,i1)
        call jamzero(i1)
        k(7,i1)=k7
        k(1,i1)=12
        return
      endif

c...Set new ID-numbers.
      k(2,i1)=kf
      kfa=abs(kf)

c...Set origin.
      k(3,i1)=1000*mste(1)+mste(2)
      if(kfa.gt.100000000) k(3,i1)=0
c..if you want the index of the parent for resonance decay msel=2
c     k(3,i1)=mste(21)

c...Decay.
      if(msel.ge.2) then
c....spin
c       k(4,i1)=int(max(1,mod(kf,10))*rn(0))
        k(4,i1)=kcp(2,1)
c...Collisions.
      else if(kfa.gt.100000000) then
        k(4,i1)=0
      else
        k(4,i1)=1000*abs(mste(22))+abs(mste(24))
      endif

      k(8,i1)=kcp(8,1)

c...Baryon number.
      iq=0
      if(kf.eq.92) then
        if(abs(kq(1,i1)).lt.10) then
          ibq1=1*isign(1,kq(1,i1))
        else
          ibq1=2*isign(1,kq(1,i1))
        endif
        if(abs(kq(2,i1)).lt.10) then
          ibq2=1*isign(1,kq(2,i1))
        else
          ibq2=2*isign(1,kq(2,i1))
        endif
        k(9,i1)=ibq1+ibq2
      else
          kflc=mod(kfa/10,10)
          if(kfa.lt.10) then
            k(9,i1)=1*isign(1,kf)
            iq=1
          else if(kfa.ge.1000000000) then ! Nucleus
            k(9,i1)=3*(mod(kfa/1000,1000)
     &               +mod(kfa/1000000,1000))*isign(1,kf)
            iq=-5
          else if(kfa.gt.1000.and.kflc.eq.0) then
            k(9,i1)=2*isign(1,kf)
            iq=2
          else if(kfa.eq.21) then
            k(9,i1)=0
            iq=3
          else
            k(9,i1)=kchg(kc,6)*isign(1,kf)
          endif
      endif


c...Color flow.
      if(iq.le.0) then
        k(10,i1)=0
        k(11,i1)=0
      endif

c...Check whether this particle is stable,resonance or jet system.
      if(kf.eq.92) then
        k(1,i1)=3
      else if(iq.eq.0) then
        if(pmas(kc,2).le.1d-7.or.mdcy(kc,1).eq.0
     $              .or.mdcy(kc,2).eq.0.or.mdcy(kc,3).eq.0)then
          k(1,i1)=1

c.....2014/2/15
          if(kf.eq.223) k(1,i1)=2 ! omega
          if(kf.eq.113) k(1,i1)=2 ! rho0
          if(kf.eq.221) k(1,i1)=2 ! eta

        else
          k(1,i1)=2
        endif
      else if(iq.eq.-5) then
        k(1,i1)=5
      else
        k(1,i1)=4
      endif

      if(ks.le.0) then
         k(1,i1)=-k(1,i1)
         ks01=mod(abs(ks)/10,10)
         if(ks01.eq.1) k(1,i1)=-10+k(1,i1)
         if(ks01.eq.2) k(1,i1)=-20+k(1,i1)
         if(ks01.eq.3) k(1,i1)=-30+k(1,i1)
      endif

c...Update collision counter.
      k(5,i1)=icltag

c...Reset multi-step history.
      if(msel.ge.2) then
        k(6,i1)=k6
      else
        ihis1=1
        if(k(6,i1).lt.0) ihis1=-1
        if(ichanel.ne.1) ihis1=-1
        ihis0=abs(k(6,i1))+abs(k6)+1
        k66=ihis1*ihis0
        k(6,i1)=min(1000,abs(k66))*isign(1,k66)
      endif

c...Update number of collision.
      if(k(7,i1).eq.-1200) then
      else if(k(7,i1).eq.0.or.msel.eq.2.or.msel.eq.6) then
        k(7,i1)=kcp(7,1)
      else
        k(7,i1)=k(7,i1)+k(7,i1)/abs(k(7,i1))
      endif

      end

c***********************************************************************

      subroutine jamvupda(i1,i2,kf,kc,taucol,ichanel)

c...Purpose: to update particle status code v() after collision.

c....Particle arrays updated here are:
c        i1    : line number of the particle.
c        i2    : 
c        v(1,)-v(4): production vertex.
c        v(5,) : life time
c        r(5,) : formation time

      implicit double precision(a-h, o-z)
      include 'jam1.inc' 
      include 'jam2.inc'
      real*8 jamdtim

      if(i1.eq.0) return
      if(k(1,i1).gt.10) return

      if(i2.ge.1) then
        if(i2.eq.mste(21)) then
          ip=1 
        else if(i2.eq.mste(23)) then
          ip=2
        else
          write(check(1),'(5(i9,1x))')i1,i2,kf,kc,ichanel
          call jamerrm(30,1,'(jamvupda:)invalid i2')
        endif
      else if(i2.eq.0) then
        ip=1
        if(rn(0).ge.0.5d0)ip=2
      endif

      if(ip.eq.1) then
        r4=rcp(4,1)
        r5=rcp(5,1)
      else
        r4=rcp(4,2)
        r5=rcp(5,2)
      endif

c...Production time
      if(mstc(52).ge.2) then
        r(4,i1)=r4+taucol
        v(4,i1)=r4
c       r(5,i1)=r5
        if(k(1,i1).gt.0) r(5,i1)=r4
      else
        v(4,i1)=pard(1)+taucol
        r(4,i1)=pard(1)+taucol
        r(5,i1)=pard(1)+taucol
      endif

c...Set Life time.
      if(k(10,i1).eq.0) then
          if(k(1,i1).gt.0) then
            tm=r5
          else
            ks0=mod(abs(k(1,i1)),10)
            if(ks0.eq.3.and.mstc(76).le.1) then
              tm=r4
            else
              tm=r5
            endif
          endif
          if(kf.ge.1000000000) then ! Nucleus
            v(5,i1)=1d+34
          else
            v(5,i1)=tm+jamdtim(1,kf,kc,k(1,i1),p(5,i1),p(4,i1))
          endif
      endif

c...Save vertex point.
      if(kf.ne.92) then
        if(ip.eq.1) then
          v(1,i1)=rcp(1,1)
          v(2,i1)=rcp(2,1)
          v(3,i1)=rcp(3,1)
        else
          v(1,i1)=rcp(1,2)
          v(2,i1)=rcp(2,2)
          v(3,i1)=rcp(3,2)
        endif
        vv=v(1,i1)**2+v(2,i1)**2+v(3,i1)**2
      if(vv.le.0.0d0) then
        write(mstc(38),*)i1,i2,kf,kc,ichanel
        write(mstc(38),*)'r1',(rcp(j,1),j=1,3),(r(j,i1),j=1,3)
        write(mstc(38),*)'r2',(rcp(j,2),j=1,3)
        write(mstc(38),*)'(jamvupda:)vv=0 ip',ip,i1,k(1,i1),k(2,i1)
     $ ,p(5,i1),(v(j,i1),j=1,5)
      endif
      endif

      end

c***********************************************************************

      subroutine jamsave(isave,isv,ip)

c...Purpose: to save and reset original informations of particles.
c...isave: =1: save, =2:restore.
c...isv : should be 1,2,3,4,5 or 6.
c...ip  : line number of particle to be saved.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'

c...Save particle types and masses.
      if(isave.eq.1) then

c...mste(21),mste(22) or mste(23),mste(24)
      mste(isv*2+19)=ip
      mste(isv*2+20)=jamcomp(k(2,ip))

      do i=1,11
      kcp(i,isv)=k(i,ip)
      end do
      kqcp(1,isv)=kq(1,ip)
      kqcp(2,isv)=kq(2,ip)
 
      do i=1,5
      pcp(i,isv)=p(i,ip)
      rcp(i,isv)=r(i,ip)
      vcp(i,isv)=v(i,ip)
      vqcp(i,isv)=vq(i,ip)
      vqcp(i+5,isv)=vq(i+5,ip)
      end do

      vpcp(:,isv)=0.0d0
      MFcp(isv)=MF_on(ip)
      dtfcp(isv)=dtfree(ip)
      if(MF_on(ip).eq.1) then
        do j=1,3
          fcp(j,isv)=force(j,ip)
          frcp(j,isv)=force(j,ip)
        end do
c       if(mstc(109).ge.1) vpcp(isv)=vpot(ip)
        spcp(isv)=pots(ip)
        vpcp(:,isv)=potv(:,ip)
      endif

c...Reset outgoing types and masses.
      if(isv.eq.1) then

        do i=25,32
         mste(i)=0
        end do
        do i=1,11
          kcp(i,3)=0
          kcp(i,4)=0
          kcp(i,5)=0
          kcp(i,6)=0
        end do
        do j=3,6
          kqcp(1,j)=0
          kqcp(2,j)=0
          MFcp(j)=0
          vpcp(:,j)=0.0d0
          dtfcp(j)=0.0d0
          do i=1,3
            fcp(i,j)=0d0
            frcp(i,j)=0d0
          end do
        end do
        do i=1,5
          do j=3,6
          pcp(i,j)=0.0d0
          vcp(i,j)=0.0d0
          rcp(i,j)=0.0d0
          end do
        end do


      endif

c...Reset the previous particle information.
      else

      do i=1,11
        k(i,ip)=kcp(i,isv)
      end do
      kq(1,ip)=kqcp(1,isv)
      kq(2,ip)=kqcp(2,isv)
      do i=1,5
        p(i,ip)=pcp(i,isv)
        r(i,ip)=rcp(i,isv)
        v(i,ip)=vcp(i,isv)
        vq(i,ip)=vqcp(i,isv)
        vq(i+5,ip)=vqcp(i+5,isv)
      end do

      pots(ip)=0.0d0
      potv(:,ip)=0.0d0
      MF_on(ip)=MFcp(isv)
      dtfree(ip)=dtfcp(isv)
      if(MF_on(ip).eq.1) then
        do j=1,3
          force(j,ip)=fcp(j,isv)
          forcer(j,ip)=frcp(j,isv)
        end do
c       if(mstc(109).ge.1) vpot(ip)=vpcp(isv)
        pots(ip)=spcp(isv)
        potv(:,ip)=vpcp(:,isv)
      endif

      endif

      end

c******************************************************************

      subroutine jamangel(pr,sig,kf1,kf2,ibar1,ibar2,t1,c1) 

c...Purpose: to determine elastic scattering angule.
c=================================================================*
c...determine the new momentum direction direction
c...(relative to the original direction).
c...calculate the azimuthal angle, t1, and the cosine of the polar
c...angle, c1.
c
c        pr     : magnitude of rel. momentum.      (input)
c        sig    : total cross section.             (input)
c        t1     : azimuthal angle                  (output)
c        c1     : cosine(polar angle)              (output)
c
c=================================================================*

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      double precision jamslope
      data  emnuc/0.938d0/

c     t1=2.0d0*paru(1)*rn(0)
      if(mstc(57).eq.0) then
        t1=2.0d0*paru(1)*rn(0)
      else
        t1=0.0d0
      endif

c     c1=1.0-2.0*rn(0)
c...Elastic scattering of identical particles
      if(kf1.eq.kf2) then
        factor=1
      else
c...Elastic scattering of non-identical particles
        factor=2
      endif

      sh=4.d0*(emnuc**2+pr**2)
      snew=sqrt(sh)
      plab=sqrt(sh*(sh-4.d0*emnuc**2))/(2.d0*emnuc)

c....Get slope parameter
      a=jamslope(kf1,kf2,ibar1,ibar2,plab,sh,snew,sig)
      ta=2.d0*pr*pr
      ata=ta*a
      if(ata.lt.1.d-7) then
        c1=1.d0-2.0*rn(0)
        return
      else
        y3=rn(0)
        tt1=log((1.d0-y3)*exp(-min(50.0d0,factor*ata))+y3)/ata
        c1=1.d0+tt1
      end if

c...Backward scattering prob.
      iback=1
      if(iback.eq.1) then
      if(mste(2).eq.1.and.kf1.ne.kf2.and.snew.le.3.0d0) then
        bprob=1.0d0
c ----  Cugnon
        if(plab.gt.0.8d0) bprob=0.8d0/plab
c ---  Niita
c       if(plab.gt.0.8) bprob=(0.8/plab)**2.0d0
 
        bprob=bprob/(1.d0+bprob)
        if(rn(0).le.bprob) c1=-c1
      end if
      end if
 
      if(abs(c1).gt.1.0d0) c1=sign(1.0d0,c1)

      end

c***********************************************************************

      function jamslope(kf1,kf2,ibar1,ibar2,plab,s,snew,sig)

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      double precision jamslope

      if(sig.le.0.0d0) then
       write(check(1),'(g10.3,1x,i9,1x,i9)')s,kf1,kf2
       call jamerrm(30,1,'(jamslope:)sig=0 s kf1 kf2=')
      endif

      if(snew.ge.10.d0) then
        b1=2.3d0
        b2=2.3d0
        if(ibar1.eq.0) then
          b1=1.4d0
          kfl1=mod(abs(kf1)/100,10)
          kfl2=mod(abs(kf1)/10,10)
          if(10*kfl1+kfl2.eq.44) b1=0.23d0  ! J/psi
        endif
        if(ibar2.eq.0) then
          b2=1.4d0
          kfl1=mod(abs(kf2)/100,10)
          kfl2=mod(abs(kf2)/10,10)
          if(10*kfl1+kfl2.eq.44) b2=0.23d0  ! J/psi
        endif
        a=2.0d0*(b1+b2)+4*s**0.0808d0-4.2d0
        goto 100
      endif


c...Anti-p p
      if(ibar1*ibar2.eq.-9) then

c...M.R.Clover,et al., Phys. Rev. C26 (1982) 2138.
c...This fit is valid from 100MeV to 2GeV lab. eng.
c       a=12.94+39.03*exp(-2.075*plab)
c...J. Cugnon, et al., Phys. Rev. C41 (1990) 1701.
        a=3.34d0

      else
c...p+n
        if(kf1.ne.kf2) then

          if( plab .le. 0.6d0 ) then
            a = 6.2d0 * ( plab - 0.225d0 )/ 0.375d0
          else if( plab .le. 1.6d0 ) then
            a =  - 1.63d0 * plab + 7.16d0
          else if( plab .le. 2.0d0 ) then
            a = 5.5d0 * plab**8 / ( 7.7d0 + plab**8 )
          else if(plab.le.4.2d0) then
            a = 5.34d0 + 0.67d0 * ( plab - 2.0d0 )
          else
c           a = 5.656 + 0.863*log(plab)
            a = 5.656d0 + 1.300d0*log(plab)
          end if

        else

          if( plab .le. 2.0d0 ) then
            a = 5.5d0 * plab**8 / ( 7.7d0 + plab**8 )
          else if(plab.le.4.2d0) then
            a = 5.34d0 + 0.67d0 * ( plab - 2.0d0 )
          else
c           a = 5.656 + 0.863*log(plab)
c......srt=19GeV
c           a = 5.656 + 1.100*log(plab)
c......srt=10GeV
            a = 5.656d0 + 1.300d0*log(plab)
          end if

        endif
      endif

 100  fac=sig/40.d0
      if((kf1.eq.2112.or.kf1.eq.2212).and.
     $         (kf2.eq.2112.or.kf2.eq.2212)) fac=1.d0
      jamslope=a*fac

      end

c******************************************************************

      subroutine jamangin(srt,pr,pr0,em1,em2,ibar1,ibar2,
     $                  kf01,kf02,kc01,kc02,kc1,kc2,t1,c1)

c...Purpose: to determine inelastic scattering angle.

c        srt    : c.m. energy                              (input)*
c        pr0    : magnitude of initial rel. momentum       (input)*
c        pr     : magnitude of final rel. momentum         (input)*
c        em1    : mass of the particle 1                   (input)*
c        em2    : mass of the particle 2                   (input)*
c        ibar1  : baryon number of the particle 1          (input)*
c        ibar2  : baryon number of the particle 2          (input)*
c        kf01   : KF code  of the initial  particle 1      (input)*
c        kf02   : KF code  of the initial  particle 2      (input)*
c        kc01   : KC code  of the initial  particle 1      (input)*
c        kc02   : KC code  of the initial  particle 2      (input)*
c        t1     : azimuthal angle                         (output)*
c        c1     : cosine(polar angle)                     (output)*

      implicit double precision(a-h, o-z)
      include 'jam2.inc'

c...Azimuthal angle.
      if(mstc(57).eq.0) then
        t1=2.0d0*paru(1)*rn(0)
      else
        t1=0.0d0
      endif

c...BB collisions.
      if(ibar1*ibar2.eq.9) then
 
c.....No strange BB collisions.
        if(abs(kchg(kc01,7))+abs(kchg(kc02,7)).eq.0) then
           if(srt.le.2.17d0) then
             call jamangdn(srt,pr,em1,em2,kc1,kc2,c1) 
           else
             call jamangrr(pr0,pr,srt,c1,em1,em2,mstc(67)) 
           endif

        else
          call jamangrr(pr0,pr,srt,c1,em1,em2,mstc(67)) 
        endif

c...MM
      else if(ibar1.eq.0.and.ibar2.eq.0) then
           call jamangrr(pr0,pr,srt,c1,em1,em2,mstc(69)) 

c...MB
      else if(ibar1*ibar2.eq.0) then
           call jamangrr(pr0,pr,srt,c1,em1,em2,mstc(68))
c...Ani-B B
      else if(ibar1*ibar2.eq.-9) then
           call jamangrr(pr0,pr,srt,c1,em1,em2,mstc(70))
      else
        write(check(1),'(i9,1x,i9,4(i4,1x)))')
     $                       kf01,kf02,ibar1,ibar2,kc01,kc02
        call jamerrm(30,1,'(jamangin:) unrecognize baryon number')
      endif

      end

c***********************************************************************

      subroutine jamangrr(pr0,pr,srt,c1,em1,em2,iang)

c...Purpse: to generate angular distribution for hh interactions.
c     input: pr0 : c.m. mom. before collision
c            pr  : c.m. mom. after  collision
c            srt : inv. mass of collision
c     output: c1 : scattering angle in c.m.
c======================================================================c
c     Method:
c (1) Follwoing HIJING, 
c         prob(pt^2)=1.0/(x**2+hip0**2)/(x**2+hic**2)
c    &                  /(1+exp((x-hip0)/0.4))
c     with hic=0.1, hip0=1.0
c     Since it is difficult to generate this pt distribution, 
c     by using Monte-Carlo Method (without Table), 
c     we fit the above function as follows,
c       prob(pt^2)=a*(
c    &             exp(-x**2/b**2)/b/b		! Gauss
c    &           c*exp(-x**2/d**2)/d/d		! Gauss
c    &          +e*(exp(-x/f)-exp(-2*x/f))/x/f)	! 2-Exp.
c     Then, each part can be generated by 
c       pt=b*sqrt(-log(1.e0-rnp*(1.0-exp(-x0**2/b**2)))) ! Gauss
c     or
c       pt=-f*log(1.0-sqrt(rn(0))*(1.0-exp(-x0/f))) ! 2-Exp.
c     Selection of each part is done by Monte-Carlo, therefore,
c     This procedure uses two random numbers
c
c (2) Consideration in small pr
c     The above pt^2 distribution of HIJING can be used around 24 GeV/c.
c     At lower energies, we considered as follows.
c         d(sigma)/d(p_t) = d(sigma)/d(theta)/pr
c     where theta is the CM scattering angle (< pi/2).
c     This is true for large pr.
c     Normalization change due to the finite range of pt is ignored.
c
c     Fitting part is done by A.Ohnishi
c======================================================================c

c     implicit double precision(a-h, o-z)
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
c     real*8 jamrnd2
c...Original Parameter Set
      parameter(ag1=1.0d0,ag2=1.60207d0,ag3=5.63708d0)
      parameter(bg1=0.0944812d0,bg2=0.196705d0,bg3=0.288274d0)
c...Parameter set II:
c     parameter(ag1=1.0d0,ag2=0.78d0,ag3=2.78d0)
c     parameter(bg1=0.11d0,bg2=0.28d0,bg3=0.29d0)

c...Ohnishi
c     parameter(ag1=1.0d0,ag2=1.0d0,ag3=2.5d0)
c     parameter(bg1=0.15d0,bg2=0.32d0,bg3=0.35d0)
c...Artificial Parameter Set: 
c     parameter(ag1=0.0d0,ag2=0.0d0,ag3=1.00d0)
c     parameter(bg1=0.11d0,bg2=0.30d0,bg3=0.32d0)

      if(iang.eq.1) then

c...Select the Functional Form of Transverse Momentum Distribution
c...pt=pr*theta, and theta is limited to be < pi/2
        ptx=pr*paru(1)/2
        rnx=rn(0)*(ag1+ag2+ag3)
        rnp=rn(0)
                if(rnx.le.ag3) then
        expf3=1.0d0-exp(-ptx/bg3)
        pt=-bg3*log(1.0d0-sqrt(rnp)*expf3)
                elseif(rnx.le.ag3+ag2) then
        expf2=1.0d0-exp(-ptx**2/bg2**2)
        pt=bg2*sqrt(-log(1.d0-rnp*expf2))
                else
        expf1=1.0d0-exp(-ptx**2/bg1**2)
        pt=bg1*sqrt(-log(1.d0-rnp*expf1))
                endif

c...In the original (HIJING) treatment, pt is real transverse momentum.
c...However, at lower energies, it easily becomes larger than pr,
c...and the measure d(pt^2) does not coincide with the normal measure
c...pr^2 d(cos(theta)), we re-interpret d(sigma)/d(pt)
c...as d(sigma)/d(theta)/pr at high energies and forward angles,
c...these become the same.
        c1=cos(pt/pr)

      else if(iang.eq.2) then

c...See for example, P.T.P.suppl. 41 and 42(1967) p291
c       a=parc(44)
c...Diffraction dissosiation
c       a=7.72+max(0.0d0,0.6d0*(log(srt*srt)-log(2.d0)))
c...binary diagram
        a=parc(44)+max(0.0d0,parc(45)*(log(srt*srt)-log(parc(46))))

        ata=4.d0*pr0*pr*a
        if(ata.lt.1.d-7) then
          c1=1.d0
          return
        else
          xran=rn(0)
          c1=1.d0+log((1.d0-xran)*exp(-min(50.0d0,ata))+xran)/ata
        end if
        if(abs(c1) .gt. 1.0d0) c1=sign(1.0d0,c1)

      else if(iang.eq.3) then

        b=parc(44)
        c=parc(45)
        sqm1=pcp(5,1)**2
        sqm2=pcp(5,2)**2
        sqm3=em1**2
        sqm4=em2**2
        sh=srt*srt

c...Determine maximum possible t range and coefficients of generation.
        sqla12=(sh-sqm1-sqm2)**2-4.d0*sqm1*sqm2
        sqla34=(sh-sqm3-sqm4)**2-4.d0*sqm3*sqm4
        tha=sh-(sqm1+sqm2+sqm3+sqm4)+(sqm1-sqm2)*(sqm3-sqm4)/sh
        thb=sqrt(max(0d0,sqla12))*sqrt(max(0d0,sqla34))/sh
        thc=(sqm3-sqm1)*(sqm4-sqm2)+(sqm1+sqm4-sqm2-sqm3)*
     &  (sqm1*sqm4-sqm2*sqm3)/sh
        thl=-0.5d0*(tha+thb)
        thu=thc/thl
        if(c.gt.0.0d0) then
          icpos=1
          thl=max(thl,-0.5d0*b/c)
          thrnd=exp(max(-50.d0,0.5d0*b*(thl-thu)))-1.d0
          thwmx=exp(max(-50.d0,0.5d0*b*thu+c*thu**2))
        else
          icpos=0
          thrnd=exp(max(-50.d0,b*(thl-thu)))-1.d0
          thwmx=exp(max(-50.d0,c*thu**2))
        endif

C...Select t according to exp(B*t + C*t^2).
        itry=0
  140   if(icpos.eq.1) then
          itry1=0
  144     th=thu+(2.d0/b)*log(1.+thrnd*rn(0))
          itry1=itry1+1
          if(exp(max(-50.d0,0.5d0*b*th+c*th**2)).lt.rn(0)*thwmx)
     $          goto 144
          if(itry1.gt.100) goto 147
        else
          itry2=0
  146     th=thu+(1.d0/b)*log(1.d0+thrnd*rn(0))
          itry2=itry2+1
          if(exp(max(-50.d0,c*th**2)).lt.rn(0)*thwmx) goto 146
          if(itry2.gt.100) goto 147
        endif

  147    itry=itry+1
         if(itry.ge.200) then
           ata=4.d0*pr0*pr*a
           if(ata.lt.1.d-7) then
             c1=1.d0
             return
           else
             xran=rn(0)
             c1=1.d0+log((1.d0-xran)*exp(-min(50.0d0,ata))+xran)/ata
           end if
           if(abs(c1) .gt. 1.0d0) c1=sign(1.0d0,c1)
           goto 148
         endif

        if(th.lt.thl.or.th.gt.thu) goto 140
  148   c1=(th-sqm1-sqm3+2.d0*sqrt(sqm1+pr0**2)*sqrt(sqm3+pr**2))/
     $        (2.d0*pr0*pr)

c...Gaussian dist.
      else if(iang.eq.4) then
        ptsq=-parc(47)**2*log(1.d0-rn(0)*(1.0d0-exp(-(pr/parc(47))**2)))
        c1=cos(sqrt(ptsq)/pr)

c...Gaussian + exponential dist.
      else if(iang.eq.5) then
        ptx=pr*paru(1)/2
        wg1=parc(47)
        wg2=parc(48)
        rnx=rn(0)*(1.0d0+parc(49))
        if(rnx.le.parc(49)) then
          pt=-wg2*log(1.0d0-rn(0)*(1.0d0-exp(-ptx/wg2)))
        else
          pt=wg1*sqrt(-log(1.d0-rn(0)*(1.0d0-exp(-(ptx/wg1)**2))))
        endif
        c1=cos(pt/pr)

      else if(iang.eq.6) then
        ptmx=pr*paru(1)/2
        ptmx2=ptmx*ptmx
        pt=sqrt(hirnd2(8,0.0d0,ptmx2))
c       pt=jamrnd2(8,0.0d0,min(6.0d0,ptmx))
        if(pt.gt.parc(70)) then
          expf=exp(-(parc(70)/parc(66))**2)
          pt=parc(66)*sqrt(-log(expf-rn(0)
     $             *(expf-exp(-ptmx2/parc(66)**2))))
        endif
        c1=cos(pt/pr)

      endif

      end

c******************************************************************

      subroutine jamangdn(srt,pr,em1,em2,kc1,kc2,c1) 

c...Purpose: to determine inelastic scattering angle for nn->nd.
c-----------------------------------------------------------------*
c   determine the new momentum direction                          *
c   for inelastic scattering (relative to the original direction).*
c   return the azimuthal angle, t1,                               *
c   and the cosine of the polar angle, c1.                        *
c                                                                 *
c        srt    : c.m. energy                              (input)*
c        pr     : magnitude of rel. momentum               (input)*
c        em1    : mass of the particle 1                   (input)*
c        em2    : mass of the particle 2                   (input)*
c        kc1    : KC code of the outgoing particle 1       (input)*
c        kc2    : KC code of the outgoing particle 2       (input)*
c        c1     : cosine(polar angle)                     (output)*
c-----------------------------------------------------------------*

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      data em_nuc,em_del/0.938d0,1.232d0/

      srte=srt-em1-em2
      srtd=srt
      srtde=srte
      srtne=srte
      prad=pr
      id1=kchg(kc1,5) 
      id2=kchg(kc2,5) 

c...Rescale expept for the delta-n final state.
      if(jamcpair(id1,id2).ne.jamcpair(id_delt,id_nucl)) then
        dmasp=pmas(kc1,1)+pmas(kc2,1)-em_nuc-em_del 
        srtd=srt-dmasp
        srtne=srte-dmasp
        if(srtne.gt.0.0d0) then
          srtde=srtne
          prad=sqrt((srtd**2-em1**2-em2**2)**2
     $              -4.0d0*(em1*em2)**2)/(2.0d0*srtd)
        end if
      end if

      if(rn(0).lt.0.5d0) then
        if(srtne.gt.0.0d0) then
          as=(3.65d0*srtde)**6
          a=as/(1.0d0+as)*srtd**4*0.14d0
          ta=-2.0d0*prad**2
          x=rn(0)
          t1=log((1-x)*exp(max(-50.0d0,2.d0*a*ta))+x)/a
          c1=1.0d0-t1/ta
          if(abs(c1).gt.1.0d0) c1=2.0d0*x-1.0d0
        else
          c1 =2.0d0*rn(0)-1.0d0
        end if
      else
        if( srtd .lt. 2.14d0 ) then
          c1=2.0d0*rn(0)-1.0d0
        else
          if(srtd.gt.2.4d0) then
            b1=0.06d0
            b3=0.4d0
          else
            b1= 29.0286d0 - 23.749d0  * srtd + 4.86549d0 * srtd**2
            b3=-30.3283d0 + 25.5257d0 * srtd - 5.30129d0 * srtd**2
          end if
          pp3=b1/(3.d0*b3)
          qq3=0.5d0 * (0.5d0 - rn(0)) / b3
          pq3=sqrt(qq3**2 + pp3**3)
          uu =(-qq3 + pq3 )**(1.d0/3.d0)
          vv =( qq3 + pq3 )**(1.d0/3.d0)
          c1 =uu-vv
          if(abs(c1).gt.1.d0) c1=c1/abs(c1)
        end if
      end if

      end

c********************************************************************* 
 
      subroutine jamedit
 
c...Purpose: to exclude dead partons/particles. 
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
 
      if(mstd(30).eq.0) return

      nbary0=nbary
c...Remove unwanted partons/particles. 
      i1=0
      do 110 i=1,nv

        if(k(1,i).ge.11) then
c         if(abs(k(9,i)).ge.3) then
          if(i.le.nbary0) then ! 2018/1/8
            nbary=nbary-1
          else
            nmeson=nmeson-1
          endif
          goto 110
        endif
 
c...Pack remaining partons/particles. origin no longer known. 
        i1=i1+1 
        idel=i-i1
        do m=1,5
          r(m,i1)=r(m,i) 
          p(m,i1)=p(m,i) 
          v(m,i1)=v(m,i) 
        end do
crqmd[AO:030108
c       if(mstc(6).ge.2.and.mstc(103).ge.1) then
        if(mstc(6).ge.2) then ! 2019/2/14
c...[AO:050727
            MF_on(i1)=MF_on(i)
            dtfree(i1)=dtfree(i)
            pots(i1)=pots(i) ! 2019/2/14
            potv(:,i1)=potv(:,i) ! 2019/2/14
c           if(mstc(109).ge.1) vpot(i1)=vpot(i)
c...]AO:050727
          do m=1,3
            force(m,i1)=force(m,i)
            forcer(m,i1)=forcer(m,i)
          end do
        endif
c...]AO:030108

         do m=1,11
           k(m,i1)=k(m,i)
         end do
         kq(1,i1)=kq(1,i)
         kq(2,i1)=kq(2,i)
         do l=1,10
           vq(l,i1)=vq(l,i)
         end do

c....Color flow information of partons.
         if(k(10,i1).ne.0) then
           k(10,i1)=k(10,i1)-idel
           k(11,i1)=k(11,i1)-idel
         endif

c...qmd:Mean field vectors.
cq       if(mstc(6).ge.1) then
cq         call medit(i1,i)
cq       endif
 
  110 continue 

      if(i1.ne.nbary+nmeson) then
        write(check(1),'(''i1 nv'',i10,1x,i10)')i1,nv
        write(check(2),'(''nbary+nmeson='',i10)')nbary+nmeson
        call jamerrm(1,2,'(jamedit:)??? i1 n=')
      endif
      nv=nbary+nmeson

c...2018/1/8 Recompute nbary and nmeson, since some dead particle
c...does not have information about k(9,i) which is zero.
      nv=i1
      nbary=0
      nmeson=0
      do i=1,nv
        if(k(1,i).ge.11) cycle
        if(abs(k(9,i)).ge.3) then
          nbary=nbary+1
        else
          nmeson=nmeson+1
        endif
      end do

      if(nv.ne.nbary+nmeson) then
        print *,'jamedit nv=',nv,' nbary=',nbary,' nmeson=',nmeson
        stop
      endif

      mstd(30)=0

      end

c***********************************************************************

      subroutine jamexch(i2,i1)
 
c...Purpose: to exchange position of event record vectors.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'

      do m=1,5
        r(m,i2)=r(m,i1) 
        p(m,i2)=p(m,i1) 
        v(m,i2)=v(m,i1) 
      end do

crqmd[AO:030108
c     if(mstc(6).ge.2.and.mstc(103).ge.1) then
      if(mstc(6).ge.2) then
c...[AO:050727
            MF_on(i2)=MF_on(i1)
            dtfree(i2)=dtfree(i1)
            pots(i2)=pots(i1) ! 2019/2/14
            potv(:,i2)=potv(:,i1) ! 2019/2/14
c           if(mstc(109).ge.1) vpot(i2)=vpot(i1)
c...]AO:050727
        do m=1,3
            force(m,i2)=force(m,i1) 
            forcer(m,i2)=forcer(m,i1) 
        end do
      endif
c...]AO:030108

      do m=1,11
        k(m,i2)=k(m,i1)
      end do
      kq(1,i2)=kq(1,i1)
      kq(2,i2)=kq(2,i1)
      do l=1,10
       vq(l,i2)=vq(l,i1)
      end do

c...qmd:Mean field vectors.
cq    if(mstc(6).ge.1) then
cq      call medit(i2,i1)
cq    endif

      end

c***********************************************************************

      subroutine jamzero(i)
 
c...Purpose: to zero the vectors.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam3.inc'

      r(:,i)=0.d0
      p(:,i)=0.d0
      v(:,i)=0.d0
      k(:,i)=0
c     k(1,i)=16
      k(7,i)=2
      vq(:,i)=0.0d0
      kq(1,i)=0
      kq(2,i)=0
      k(1,i)=13

      force(:,i)=0.0d0
      forcer(:,i)=0.0d0
      dtfree(i)=0.0d0
      pots(i)=0.0d0
      potv(:,i)=0.0d0
      MF_on(i)=0
      IcanColl(i)=0

      end

c***********************************************************************

      subroutine jamprcl(indd,nadd)

c...Purpose: to print collision/decay information.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      dimension indd(100)
      character*16 chaf1,chaf2,chaf3,chaf4
      character check1*120
      dimension ptot(4)

      ih=mstc(38)
      icltyp=mste(2)
      if(icltyp.lt.0) goto 1000  ! decay
      srt=pare(2)
      ichanel=mste(1)
      sig=pare(4)
      if(mstc(8).lt.3) goto 100

      if(ichanel.eq.1) then
        write(ih,705)pard(1),ichanel,srt,sig
      else if(ichanel.ge.2) then
        if(icltyp.eq.1) then
          write(ih,700)pard(1),ichanel,mste(4),srt,sig
        else if(icltyp.eq.2) then
          write(ih,701)pard(1),ichanel,srt,sig
        else if(icltyp.eq.3) then
          write(ih,702)pard(1),ichanel,srt,sig
        else if(icltyp.eq.4) then
          write(ih,703)pard(1),ichanel,srt,sig
        else if(icltyp.eq.5) then
          write(ih,704)pard(1),ichanel,srt,sig
        endif
      endif

      call pjname(kcp(2,1),chaf1)
      call pjname(kcp(2,2),chaf2)
      i1=mste(21)
      i2=mste(23)
      i3=mste(25)
      i4=mste(27)

c....Absorption.
      if(ichanel.eq.3) then

        ia=mste(25)
c         write(ih,*)i1,(k(j,i1),j=1,11),p(5,i1)
c         write(ih,*)i2,(k(j,i2),j=1,11),p(5,i2)
          em1sq=pcp(4,1)**2-pcp(1,1)**2-pcp(2,1)**2-pcp(3,1)**2
          em2sq=pcp(4,2)**2-pcp(1,2)**2-pcp(2,2)**2-pcp(3,2)**2 
          em3sq=p(4,ia)**2-p(1,ia)**2-p(2,ia)**2-p(3,ia)**2 
          write(ih,*)i1,kcp(2,1),pcp(5,1),max(0d0,em1sq)
          write(ih,*)i2,kcp(2,2),pcp(5,2),max(0d0,em2sq)
          write(ih,*)ia,(k(j,ia),j=1,2),p(5,ia),max(0d0,em3sq)


        call pjname(k(2,ia),chaf3)
        write(ih,800)i1,chaf1,pcp(5,1),i2,chaf2,pcp(5,2)
     $               ,ia,chaf3,p(5,ia)

        write(ih,812)i1,kcp(1,1),rcp(4,1),rcp(5,1),vcp(5,1),
     $  chaf1,pcp(5,1),kcp(2,1),
     $  i2,kcp(1,2),rcp(4,2),rcp(5,2),vcp(5,2),
     $  chaf2,pcp(5,2),kcp(2,2)

        write(ih,803)ia,k(1,ia),r(4,ia),r(5,ia),v(5,ia),chaf3,p(5,ia)
     $  ,k(2,ia)
        write(ih,*)'r4 r5=',r(4,ia),v(4,ia)

c...Two-body collisions.
      else if(ichanel.ge.1.and.ichanel.le.6) then

        call pjname(kcp(2,3),chaf3)
        call pjname(kcp(2,4),chaf4)
        write(ih,811)chaf1,pcp(5,1)
     $              ,chaf2,pcp(5,2)
     $              ,chaf3,p(5,i3)
     $              ,chaf4,p(5,i4)
         write(ih,812)
     $ i1,kcp(1,1),rcp(4,1),rcp(5,1),vcp(5,1),chaf1,pcp(5,1),kcp(2,1),
     $ i2,kcp(1,2),rcp(4,2),rcp(5,2),vcp(5,2),chaf2,pcp(5,2),kcp(2,2)
         write(ih,812)
     $  i1,k(1,i3),r(4,i3),r(5,i3),v(5,i3),chaf3,p(5,i3),k(2,i3)
     $ ,i2,k(1,i4),r(4,i4),r(5,i4),v(5,i4),chaf4,p(5,i4),k(2,i4)
       write(ih,*)'r4 v4',r(4,i3),v(4,i3),r(4,i4),v(4,i4)
       write(ih,*)'r3',(r(j,i3),j=1,4)
       write(ih,*)'v3',(v(j,i3),j=1,4)
       write(ih,*)'r4',(r(j,i4),j=1,4)
       write(ih,*)'v4',(v(j,i4),j=1,4)

      else if(ichanel.ge.11) then
        write(ih,*)'three body collision'
      else if(ichanel.eq.-1) then

        write(ih,822)chaf1,pcp(5,1),chaf2,pcp(5,2)

      endif

 100  if(mstc(8).eq.3.or.mstc(8).ge.5) then
         write(check1,850)ichanel,i1,kcp(2,1),pcp(5,1)
     $ ,i2,kcp(2,2),pcp(5,2)
     $ ,kcp(2,3),pcp(5,3),kcp(2,4),pcp(5,4)
        call jamcheck('<<after coll>> '//check1)

      write(ih,*)'after coll'
      ptot(1)=0.0d0
      ptot(2)=0.0d0
      ptot(3)=0.0d0
      ptot(4)=0.0d0
      do i=1,nv
        if((k(1,i).ge.1.and.k(1,i).le.10).or.k(1,i).lt.0) then
c         if(k(2,i).ne.92.and.r(1,i).ne.v(1,i)) then
c           write(ih,*)'v4? i=',i,k(1,i),k(2,i)
c           write(ih,*)'r=',(r(j,i),j=1,5)
c           write(ih,*)'v=',(v(j,i),j=1,5)
c         endif
          ptot(1)=ptot(1)+p(1,i)
          ptot(2)=ptot(2)+p(2,i)
          ptot(3)=ptot(3)+p(3,i)
          ptot(4)=ptot(4)+p(4,i)
          kf=k(2,i)
          if(kf.eq.111.or.abs(kf).eq.211) then
             if(p(5,i).ge.0.16) then
               emff=sqrt(p(4,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2)
               print *,'after the collision pion mass',p(5,i),emff,kf,i
               stop
             endif
          endif
        endif
      end do

      write(ih,*)'coll ptot=',(ptot(j),j=1,4)
c     call jamlist(1)

          if(pcp(5,1).le.0.28.and.
     &     (kcp(2,1).eq.113.or.abs(kcp(2,1)).eq.213)) then
           write(6,*)'coll.1 rho mass ? ',pcp(5,1),kcp(2,1)
           stop
          endif

          if(pcp(5,2).le.0.28.and.
     &     (kcp(2,2).eq.113.or.abs(kcp(2,2)).eq.213)) then
           write(6,*)'coll2. rho mass ? ',pcp(5,2),kcp(2,2)
           stop
          endif

      endif

705   format(/'Elastic collision! time=',f8.3,'fm/c ichanel=',i3,
     $ ' srt=',f11.3,'GeV sig=',f9.3,'mb')
700   format(/'B-B collision! time=',f8.3,'fm/c ichanel=',i3,
     $ 'diffra=',i3,
     $ ' srt=',f11.3,'GeV sig=',f9.3,'mb')
701         format(/'m-B collision! time=',f8.3,'fm/c ichanel',i3,
     $ ' srt=',f11.3,'GeV sig=',f9.3,'mb')
702         format(/'m-m collision! time=',f8.3,'fm/c ichanel',i3,
     $ ' srt=',f11.3,'GeV sig=',f9.3,'mb')
703   format(/'B-antiB collision! time=',f8.3,'fm/c ichanel',i3,
     $  ' srt=',f11.3,'GeV sig=',f9.3,'mb')
704   format(/'Hadron-Parton collision! time=',f8.3,'fm/c ichanel',i3,
     $  ' srt=',f11.3,'GeV sig=',f9.3,'mb')
 800  format('absorb: ',i4,1x,a9,1x,f6.3,' + ',i4,1x,a9,f6.3,
     $      '=>',i5,1x,a9,1x,f6.3)
 803  format('(',i5,') k1=',i4,' r4=',f8.4,' r5=',f8.4,' v5=',d9.4,
     $   1x,a9,1x,f9.4,1x,i7)
 811  format('coll.: ',a8,1x,f6.3,' + ',a8,1x,f6.3,
     $          ' => ',a8,1x,f6.3,' + ',a8,1x,f6.3)
 812  format('(',i5,') k1=',i4,' r4=',f8.4,' r5=',f8.4,' v5=',d9.4,
     $    1x,a9,1x,f9.4,1x,i7/,
     $    '(',i5,') k1=',i4,' r4=',f8.4,' r5=',f8.4,' v5=',
     $    d9.4,1x,a9,1x,f9.4,1x,i7)
 822  format('Pauli-blocked: ',a8,1x,f6.3,' + ',a8,1x,f6.3)
 850  format(i2,i4,i9,g10.3,'+',i4,i9,g10.3,'->',i9,g10.3,'+',i9,g10.3)


      return

c....Information of decay
1000  continue

      if(mstc(8).eq.2) goto 200

c...Print mother.
        ind1=indd(1)
        ind2=indd(2)
        if(ind1.eq.0.or.ind2.eq.0) then
          print *,'after decay information ind1=',ind1,ind2,nadd
        endif
        kc1=jamcomp(k(2,ind1))
        kc2=jamcomp(k(2,ind2))
        if(kc1.le.0.or.kc2.le.0) then
          print *,'after decay information ind1=',ind1,ind2,nadd
          print *,'kc1 kc2=',kc1,kc2
          print *,'mste(1)=',mste(1),'mste(2)=',mste(2),pare(2)
          print *,'kf1=',kchg(kc1,4)
          stop
        endif
      write(ih,901)mste(21),pard(1),kcp(1,1),kcp(2,1),pcp(5,1),nadd
     $  ,chaf(mste(22),(3-isign(1,kcp(2,1)))/2)
     $  ,chaf(kc1,(3-isign(1,k(2,ind1)))/2)
     $  ,chaf(kc2,(3-isign(1,k(2,ind2)))/2)

      write(ih,'(''r='',5(g10.3,1x),''v5='',g15.5)')
     $ (rcp(j,1),j=1,5),vcp(5,1)

c...Print daughters.
      do i=1,nadd
        ind1=indd(i)
        kc1=jamcomp(k(2,ind1))
        write(ih,902)ind1,(k(l,ind1),l=1,2),p(5,ind1)
     $          ,r(4,ind1),v(4,ind1),r(5,ind1),v(5,ind1)
     $        ,chaf(kc1,(3-isign(1,k(2,ind1)))/2)

           write(ih,*)'p=',(p(j,ind1),j=1,5)
           write(ih,*)'r=',(r(j,ind1),j=1,5)
           write(ih,*)'v=',(v(j,ind1),j=1,5)

          if(p(5,ind1).le.0.28.and.
     &          (k(2,ind1).eq.113.or.abs(k(2,ind1)).eq.213)) then
           write(6,*)'rho mass ? ',p(5,ind1),k(2,ind1)
           stop
          endif

      end do 
c     if(mstc(8) .ge.3) call jamlist(1)

 200  if(mstc(8).eq.2.or.mstc(8).ge.5) then
        if(icltyp.eq.-1) then
          call jamcheck('<<after decay>>')
        else
          call jamcheck('<<after Final decay>>')
        endif

      ptot(1)=0.0d0
      ptot(2)=0.0d0
      ptot(3)=0.0d0
      ptot(4)=0.0d0
      do i=1,nv
        if((k(1,i).ge.1.and.k(1,i).le.10).or.k(1,i).lt.0) then
c         if(k(2,i).ne.92.and.r(1,i).ne.v(1,i)) then
c           write(ih,*)'v1? i=',i,k(1,i),k(2,i)
c           write(ih,*)'r=',(r(j,i),j=1,5)
c           write(ih,*)'v=',(v(j,i),j=1,5)
c         endif
          ptot(1)=ptot(1)+p(1,i)
          ptot(2)=ptot(2)+p(2,i)
          ptot(3)=ptot(3)+p(3,i)
          ptot(4)=ptot(4)+p(4,i)
          kf=k(2,i)
          if(kf.eq.111.or.abs(kf).eq.211) then
             if(p(5,i).ge.0.16) then
               emff=sqrt(p(4,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2)
               print *,'after the collision pion mass',p(5,i),emff,kf,i
               stop
             endif
          endif
        endif
      end do

      write(ih,*)'decay ptot=',(ptot(j),j=1,4)

      endif


 901  format('<<Decay i=',i6,1x,f8.3,'fm/c ks=',i4,' kf=',i7,
     $  ' em=',g10.3,' nadd=',i4,1x,a8,1x,a8,1x,a8)
 902  format(i5,' k1=',i4,' k2=',i7,
     $ ' m=',f7.3,' r4=',g9.3,' v4',g9.3,' r5=',g9.3,' v5=',g9.3,1x,a8)


      end

c***********************************************************************

      subroutine jamcheck(echar)

c...Purpose: to check possible errors.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      dimension ptotc(3),dps(5)
      character echar*(*)

      ierror=0
      ih=mstc(38)
c     if(mstc(8).ne.2) then
c       write(ih,'(/,''  (jamcheck) nv nbary nmeson'',3(i5,1x))')
c    $  nv,nbary,nmeson
c       write(ih,*)echar
c     endif

c...Check baryon number and charge conservations.
      ntchg=0
      numbary=0 
      nstr=0
      ptotc(1)=0.0d0
      ptotc(2)=0.0d0
      ptotc(3)=0.0d0
      ekin=0.0d0
      nq=0

c...Loop over all particles.
      do 1000 i=1,nv

        k1=k(1,i)
        if(k1.gt.10) goto 1000  ! Skip dead particles

c.....Check k1.
        if(k1.gt.5.or.k1.lt.-33) then
           write(ih,*)'Funny K1 ',k1,k(2,i),p(5,i)
           ierror=1
        endif

c...Check particle code.
        kf=k(2,i)
        kc=jamcomp(kf)
        if(abs(kf).lt.1000000000) then ! Nucleus
        if(kf.eq.0.or.(kc.le.0.or.kc.gt.mstu(6))) then
          write(ih,*)'Invalid KF code',i,kf,kc
          ierror=2
          goto 100
        endif
        endif
        if(kf.eq.11103.or.kf.eq.22101.or.kf.eq.22103) then
          write(ih,*)'Invalid KF code',i,kf,kc
          ierror=2
          goto 100
        endif

c       em0=pmas(kc,1)
c       if(abs(em0-p(5,i)).gt.1d-2) then
c         print *,'p(5,i) is different from the mass0',k(2,i),p(5,i)
c         print *,'pmas',em0
c         print *,'diff',em0-p(5,i)
c       endif

c...Check string.
        if(kf.eq.92) then
          ibb=kfprop(kq(1,i),2)+kfprop(kq(2,i),2)
          if(k(9,i).ne.ibb) then
            write(ih,*)'bad k9',i,k(2,i),p(5,i),k(9,i),kq(1,i),kq(2,i)
            call jamerrm(30,0,
     $             '(jamcheck:) invalid baryon number of string')
          endif

c...Check invariant masses in parton systems.
          do j=1,5
          dps(j)=vq(j,i)+vq(j+5,i)
          end do
          emsq=dps(4)**2-dps(1)**2-dps(2)**2-dps(3)**2

          if(emsq.lt.(0.9d0*parj(32)+dps(5))**2) then
            if(emsq.gt.0.0d0) then
            write(check(1),'(i9,1x,i9,1x,i9)')kf,kq(1,i),kq(2,i)
            write(check(2),'(''mass='',g13.4)')sqrt(emsq)
            call jamerrm(3,2,'(jamcheck:) too small mass in jet system')
            write(ih,*)'too small mass in jet system'
            write(ih,*)check(1)
            write(ih,*)check(2)
c           ierror=101
            else
              write(ih,*)'string mass not defiend i kf emsq',i,kf,emsq
              write(ih,*)'p',k(1,i),(p(j,i),j=1,5)
              write(ih,*)'vq1',kq(1,i),(vq(j,i),j=1,5)
              write(ih,*)'vq2',kq(2,i),(vq(j,i),j=6,10)
              call jamerrm(30,0,'(jamcheck:) funny jet system')
            endif 
          endif

        else if(abs(kf).gt.1000000000) then ! Nucleus
        else if(k(9,i).ne.kchg(kc,6)*isign(1,k(2,i))) then
           write(ih,*)'invalid baryon num i k em',i,k(2,i),p(5,i),k(9,i)
           call jamerrm(30,0,'(jamcheck:) invalid baryon number')
        endif

        if(abs(kf).lt.1000000000) then ! Nucleus
        id=kchg(kc,5)

c....Check delta mass
        emeff=p(4,i)**2 -p(1,i)**2-p(2,i)**2-p(3,i)**2
        if(kf.gt.100.and.emeff.le.0d0) then
           write(ih,*)'jamcheck:emeff<0?',i,emeff,kf,
     &           (p(j,i),j=1,5),pots(i),potv(0,i)
           ierror=2
        endif

        if(id.eq.id_delt.or.id.eq.id_nucls.or.id.eq.id_delts) then
          if(p(5,i).le.1.076d0) then
          write(ih,*)'Invalid Delta mass em=',i,kf,emeff,p(5,i),pots(i)
          ierror=3
          endif
        endif

        if(abs(kf).eq.213.and.p(5,i).lt.0.28) then
          write(ih,*)'Invalid rho mass',p(5,i),emeff
          ierror=4
        endif

c       if(mstc(109).le.1.and.(kf.eq.3122.and.p(5,i).le.1.00d0)) then
        if(kf.eq.3122.and.p(5,i).lt.1.11568d0) then
          write(ih,*)'Invalid lambda mass em=',p(5,i),emeff
          ierror=4
        endif
c       endif

        endif

        if(kf.gt.100.and.p(5,i).le.0d0) then
          write(ih,*)'Invalid hadron mass em=',kf,p(5,i),pots(i)
          ierror=21
        endif
  

c...Check partons.
        if(k1.eq.4) then
          do ii=k(10,i),k(11,i)
            kfl=mod(abs(k(2,ii))/10,10)
            if(abs(k(2,ii)).gt.1000.and.kfl.ne.0) then
              write(ih,*)'(jamcheck:) invalid color',k(1,ii),k(2,ii)
              ierror=1000
            endif
          end do
        endif

c...Check r vector and vertex.
c       if(abs(k(7,i)).ne.1) then
c         rr=r(1,i)**2+r(2,i)**2+r(3,i)**2
c         if(rr.le.0.0d0.or.r(4,i).gt.1d+6) then
c           write(6,*)'r?? r',i,(k(j,i),j=1,4),(r(j,i),j=1,3),p(5,i)
c           ierror=4
c         endif
c         vv=v(1,i)**2+v(2,i)**2+v(3,i)**2
c         if(vv.le.0.0d0.or.v(4,i).gt.1d+15) then
c           write(6,*)'vertex?? v',i,(k(j,i),j=1,4),(v(j,i),j=1,3)
c           ierror=5
c         endif
c       endif

c...Baryon number,charge,strangeness.
        numbary = numbary + k(9,i)
        ntchg   = ntchg   + jamk(1,i)
        nstr    = nstr    + jamk(2,i)

c...Total momentum and energy.
        if(kf.eq.92) then
        ptotc(1)=ptotc(1)+vq(1,i)+vq(6,i)
        ptotc(2)=ptotc(2)+vq(2,i)+vq(7,i)
        ptotc(3)=ptotc(3)+vq(3,i)+vq(8,i)
        ekin=ekin+vq(4,i)+vq(9,i)
        else
        ptotc(1)=ptotc(1)+p(1,i)
        ptotc(2)=ptotc(2)+p(2,i)
        ptotc(3)=ptotc(3)+p(3,i)
        ekin=ekin+p(4,i)
        endif

 1000 continue

        write(ih,*)'current time=',mstd(23),mstc(3),pard(1),pare(10),
     &        parc(2)*mstd(23)

c...Check total baryon number.
      if( numbary/mstc(5) .ne. 3*mstd(12)) then
        write(ih,*)' Baryon number not conserved='
     $                     ,mstd(12)*3,numbary
        ierror=123
c       call jamlist(1)
      end if

c...Check total charge.
      if( ntchg/mstc(5) .ne. mstd(13)*3 ) then
        write(ih,*)'Charge error ini. and tot. charge=',
     $                                        mstd(13),ntchg/3
        ierror=66
      end if

c...Check total strangeness.
      if( nstr/mstc(5) .ne. mstd(14) ) then
        write(ih,*)'Strangeness not conserved=',mstd(14),nstr
        if(mste(2).ge.1.and.mste(1).ne.6) then
          j1=mste(21)
          j2=mste(23)
          j3=mste(25)
          j4=mste(27)
          write(ih,*)'i1',j1,kcp(2,1)
     $                ,chaf(mste(22),(3-isign(1,kcp(2,1)))/2),pcp(5,1)
          write(ih,*)'i2',j2,kcp(2,2)
     $                ,chaf(mste(24),(3-isign(1,kcp(2,2)))/2),pcp(5,2)
          write(ih,*)'i3',j3,kcp(2,3)
     $                ,chaf(mste(26),(3-isign(1,kcp(2,3)))/2),pcp(5,3)
          if(kcp(2,3).eq.92) then
            write(ih,*)' q1 q2',kq(1,j3),kq(2,j3)
          endif
          if(mste(1).ne.3) then
            write(ih,*)'i4',j4,kcp(2,4)
     $                ,chaf(mste(28),(3-isign(1,kcp(2,4)))/2),pcp(5,4)
            if(kcp(2,4).eq.92) then
              write(ih,*)' q1 q2',kq(1,j4),kq(2,j4)
            endif
          endif
        endif
        ierror=77
      end if

100   continue

c...Check total momentum conservation.
      delp1=abs(pard(9)-ptotc(1)/mstc(5))
      delp2=abs(pard(10)-ptotc(2)/mstc(5))
      delp3=abs(pard(11)-ptotc(3)/mstc(5))
      if(   ( delp1 .ge. 1.0d-5 )
     $  .or.( delp2 .ge. 1.0d-5 )
     $  .or.( delp3 .ge. 1.0d-5 ) ) then
          write(ih,*)' momentum is not conserved'
c         write(ih,*)'ptot_x ',pard(9),ptotc(1)/mstc(5)
c         write(ih,*)'ptot_y ',pard(10),ptotc(2)/mstc(5)
c         write(ih,*)'ptot_z ',pard(11),ptotc(3)/mstc(5)
          ierror=20
      end if
          write(ih,*)'ptot_x ',pard(9),ptotc(1)/mstc(5)
          write(ih,*)'ptot_y ',pard(10),ptotc(2)/mstc(5)
          write(ih,*)'ptot_z ',pard(11),ptotc(3)/mstc(5)

c...Check total energy conservation.
      etot=ekin/mstd(11)/mstc(5)

c...Mean field potential.
      if(mstc(6).ge.2.and.mstc(6).lt.100) then ! rqmd
        call jambuue(ekin,epot,etot)
      endif

c     if(abs(etot-pard(13)).gt.0.1d0*pard(13)) then
c       write(ih,*)'Total energy not conserved',pard(13),etot
c       ierror=20
c     endif

      if(mstc(8).ne.2)
     $ write(ih,*)'***(jamcheck:) e.con(MeV/A)=',
     $  (pard(13)-etot)*1000,delp1*1000,delp2*1000,delp3*1000
        
      if(mstc(8).ge.5.and.ierror.eq.20) then
        write(ih,*)echar
c       if(ierror.eq.20) return
        call jamerrm(30,0,'(jamcheck:)Total mom.or energy funny')
      endif

      if(mstc(8).ge.5.and.ierror.ne.0) then
        write(ih,*)echar
        write(check(1),'(''ierror='',i4)')ierror
        call jamerrm(30,1,'(jamcheck:) something was wrong')
      endif

      end

c***********************************************************************

      subroutine jampauli(ip,ntag,phase)

c...Purpose: to calculate the Pauli blocking factor
c----------------------------------------------------------------------*
c...ntag: Flag which tells if phase-space is  Pauli-blocked.
c         ntag =  0 => phase space open
c         ntag =  1 => phase space blocked
c...phase  -  phase space factor
c----------------------------------------------------------------------*

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
c...BUU:Gauss smearing.
      common/pq/cm(130,130),dm(130,350),minner,mouter,moutep
c...Commonblock for LPC.
      common /jamlpc1/wa0,wr0,ww,rho0,ext,idens,icscde

      parameter (epsx=-20.0d0)
      logical first
      save prad2,pcount,numpau,first
      data first/.true./

      ntag=0
      phase=0.0d0

      if(mstc(6).ge.2.and.mstc(6).lt.100) goto 1000   ! BUU/rqmd
      if(mstc(6).lt.-100) goto 2000     ! LPC

c.....One-body phase space function.
      if(mstc(56).eq.1) then
          cpw=0.5d0/parc(15)
          cph=2.0d0*parc(15)/paru(3)**2
          cpc=4.0d0
c....Fusimi function.
      else if(mstc(56).eq.2) then
          cpw=0.25d0/parc(15)
          cph=parc(15)/paru(3)**2
          cpc=0.5d0
      else
         call jamerrm(30,0,'(jampauli:)invalid mstc(56)')
      endif

      kf1=k(2,ip)
      rx=r(1,ip)
      ry=r(2,ip)
      rz=r(3,ip)
      px=p(1,ip)
      py=p(2,ip)
      pz=p(3,ip)
c     em1=p(5,ip)
      e1=p(4,ip)
      em1sq=e1*e1 - px*px - py*py - pz*pz

      do 100 i=1,nbary

        if(i.eq.ip) goto 100
        if(k(2,i).ne.kf1) goto 100
        if(k(1,i).gt.10.or.k(1,i).le.0) goto 100 ! dead particle


        if(mstc(52).ge.2.and.mstc(52).le.11) then
          dt=pard(1)-r(4,i)
          if(dt.lt.0.0d0) goto 100  ! not formed.
          vx=p(1,i)/p(4,i)
          vy=p(2,i)/p(4,i)
          vz=p(3,i)/p(4,i)
          if(mstd(91).ge.1.and.MF_on(i).eq.1) then
            pvx=p(1,i)-potv(1,i)
            pvy=p(2,i)-potv(2,i)
            pvz=p(3,i)-potv(3,i)
            ee=sqrt((p(5,i)+pots(i))**2+pvx**2+pvy**2+pvz**2)
            vx=pvx/ee
            vy=pvy/ee
            vz=pvz/ee
          endif

          drx=rx-r(1,i)-dt*vx
          dry=ry-r(2,i)-dt*vy
          drz=rz-r(3,i)-dt*vz

        else
          drx=rx-r(1,i)
          dry=ry-r(2,i)
          drz=rz-r(3,i)
        endif

        dpx=px-p(1,i)
        dpy=py-p(2,i)
        dpz=pz-p(3,i)
        eij=e1+p(4,i)
        dbx=(px+p(1,i))/eij
        dby=(py+p(2,i))/eij
        dbz=(pz+p(3,i))/eij
        bij2=dbx*dbx+dby*dby+dbz*dbz
        em2sq=p(4,i)**2 - p(1,i)**2 - p(2,i)**2 - p(3,i)**2
        gam2=eij**2
     $         /(eij**2-((px+p(1,i))**2+(py+p(2,i))**2+(pz+p(3,i))**2))
        r2=drx*drx+dry*dry+drz*drz+gam2*(drx*dbx+dry*dby+drz*dbz)**2
        p2=dpx*dpx+dpy*dpy+dpz*dpz+(-(e1-p(4,i))**2
     $      +gam2*((em1sq-em2sq)/eij)**2)
        expap=-cpw*r2-cph*p2
        if(expap.gt.epsx) phase=phase+exp(expap)

 100  continue
      phase=cpc*phase/mstc(5)
      goto 9000


c...BUU:
 1000 continue
      if(first) then
        numpau=20
        rad=3.0d0
        zahl=1.d0/6.d0
        prad=(3d0*zahl/16.0d0/paru(1))**(1d0/3d0)*(paru(2)*paru(3))/rad
        prad2=prad**2
        pcount=zahl*dble(mstc(5))
        first=.false.
      endif

      sek=0.0d0
      in=minner
      iq=mouter
      ir=moutep
c     i3=mste(25)
c     i4=mste(27)
c     ix0 = nint((r(1,i3)+r(1,i4))/2d0)
c     iy0 = nint((r(2,i3)+r(2,i4))/2d0)
c     iz0 = nint((r(3,i3)+r(3,i4))/2d0)
      ix0 = nint(r(1,ip))
      iy0 = nint(r(2,ip))
      iz0 = nint(r(3,ip))
 
      px=p(1,ip)
      py=p(2,ip)
      pz=p(3,ip)
      do 200 j=1,nbary

c....dead particles.
        if(k(1,j).gt.10.or.k(1,j).le.0) goto 200
        if(mod(k(8,j)-1,numpau).ne.0) goto 200
c....not nucleons.
        if(k(2,j).ne.2112.and.k(2,j).ne.2212) goto 200

        if((px-p(1,j))**2+(py-p(2,j))**2+(pz-p(3,j))**2.gt.prad2)
     $ goto 200

        jx=nint(r(1,j))
        kx=ix0-jx
        if(abs(kx).gt.ir) goto 200
        jy=nint(r(2,j))
        ky=iy0-jy
        if(abs(ky).gt.ir) goto 200
        jz=nint(r(3,j))
        kz=iz0-jz
        if(abs(kz).gt.ir) goto 200

        lx=nint(dble(2*in+1)*(r(1,j)-dble(jx)))
        if(abs(lx) .eq. in+1) lx = lx/abs(lx) * in
        ly=nint(dble(2*in+1)*(r(2,j)-dble(jy)))
        if(abs(ly) .eq. in+1) ly = ly/abs(ly) * in
        lz=nint(dble(2*in+1)*(r(3,j)-dble(jz)))
        if(abs(lz) .eq. in+1) lz = lz/abs(lz) * in

        ic=1+(lz+in)+(ly+in)*(2*in+1)+(lx+in)*(2*in+1)**2
        ie=1+(kz+ir)+(ky+ir)*(2*ir+1)+(kx+ir)*(2*ir+1)**2

        if(dm(ic,ie).le.1.0d-12) dm(ic,ie) = 0.0d0
        sek=sek+dm(ic,ie)*dble(mstc(5))
  200 continue
      phase=min(1.0d0,sek/pcount*dble(numpau))

 9000 continue
      if(phase.gt.rn(0)) ntag=1
      return

c...LPC:
 2000 continue
      ppn=p(1,ip)**2+p(2,ip)**2+p(3,ip)**2
      rr=sqrt(r(1,ip)**2+r(2,ip)**2+r(3,ip)**2)

      if(idens.eq.1) then
        rhop=0.0d0
c...To help Just On the Surface !
        if(rr.lt.wr0+1.0d-3) rhop=rho0
      elseif(rr.lt.wr0) then
        rhop=max(0d0,rho0*(1+ww*(rr/wr0)**2)/(1+exp((rr-wr0)/wa0)))
      else
        expx=exp(-(rr-wr0)/wa0)
        rhop=max(0d0,rho0*expx*(1+ww*(rr/wr0)**2)/(1+expx))
      endif
      pff=0.6104643d0   !...pff=hc*(3*pi*pi)**(1.0/3.0)
      pf=pff*rhop**(1.0d0/3.0d0)
      if(ppn.lt.pf*pf) ntag=1
  

      end

c***********************************************************************

      function jamisjet(pcm,bsq,jp,jt,sig,sigel)

c...Purpose: to determine number of jet.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common/hijdat/hidat0(10,10),hidat(10)
      save  /hiparnt/,/hijdat/
      dimension pcm(5)

c....Fit parameters for jet cross section:
c....Duke-Owens set 1 structure functions with pt cut pt0=2.0GeV
      data a1,b1,c1,d1,e1/0.0704853d0,7.79423d0,2.56045d0,1.51995d0,
     $     137.426d0/
      data a3,b3,c3/ 0.556925d0,0.315648d0,49251.8d0/
      data a2,b2,c2/0.0912787d0,123.004d0,2.75874d0/
      data ifit/1/

c...jamisjet=-2: no interaction
c...        =-1: elastic
c...        = 0: soft
c...        = n: number of jet
c...mb=0.1*fm, YP is in fm,HIPR1(31) is in mb
c...hipr1(30): cross section 57.0mb sigma_{soft}=2*sigma_0
c...hipr1(31): cross section 28.5mb sigma_0
c...hipr1(40): value of pi=3.14159
c...hint1(14): the jet production cross section without
c...           nuclear shadowing effect.
c...hint1(18): the effective cross section for jet production
c...bsq is a impact parameter squard

      hint1(1)=pcm(5)
      ihnt2(5)=k(2,jp)
      ihnt2(6)=k(2,jt)

c...Soft cross section
        i=0
 20     i=i+1
        if(i.eq.10) go to 30
        if(hidat0(10,i).le.hint1(1)) go to 20
 30     if(i.eq.1) i=2
         hidat(5)=hidat0(5,i-1)+(hidat0(5,i)-hidat0(5,i-1))
     &     *(hint1(1)-hidat0(10,i-1))/(hidat0(10,i)-hidat0(10,i-1))
        hipr1(31)=hidat(5)
        hipr1(30)=2.0d0*hidat(5)

      if(ifit.eq.1) then
        srt=hint1(1)
        if(srt.le.25.0d0) then
          hint1(14)=a1*(srt-b1)**c1/(srt**d1+e1)
        else if(srt.le.300.0d0) then
          hint1(14)=a2*log(srt*srt/b2)**c2
        else
          hint1(14)=a3*(srt*srt-c3)**b3
        endif
        goto 2000
      endif

      if(mstd(11).ne.2) then
 
c....Calculate jet cross section

        if(mod(abs(k(1,jp)),10).eq.4) then
          jq1=k(2,k(10,jp))
          jqq1=k(2,k(11,jp))
          call kfcnst(jq1,jqq1,kf01,0.0d0)
          ihnt2(5)=kf01
        endif
        if(mod(abs(k(1,jt)),10).eq.4) then
          jq1=k(2,k(10,jt))
          jqq1=k(2,k(11,jt))
          call kfcnst(jq1,jqq1,kf02,0.0d0)
          ihnt2(6)=kf02
        endif

        call crsjet
      endif

c...Account for nuclear shadowing effect.
 2000 continue
      if(mstc(84).eq.1) then
           rrb1=min((r(1,jp)**2+r(2,jp)**2)
     &         /1.2d0**2/dble(ihnt2(1))**0.6666667d0,1.0d0)
           rrb2=min((r(1,jt)**2+r(2,jt)**2)
     &          /1.2d0**2/dble(ihnt2(3))**0.6666667d0,1.0d0)
           aphx1=hipr1(6)*4.0d0/3.0d0*(ihnt2(1)**0.3333333d0-1.0d0)
     &           *sqrt(1.0d0-rrb1)
           aphx2=hipr1(6)*4.0d0/3.0d0*(ihnt2(3)**0.3333333d0-1.0d0)
     &           *sqrt(1.0d0-rrb2)

        hint1(18)=hint1(14)-aphx1*hint1(15)
     &                  -aphx2*hint1(16)+aphx1*aphx2*hint1(17)
      else
        hint1(18)=hint1(14)
      endif

      sjet=hint1(18)
      s0=hipr1(31)
      ssoft=2*s0
      r2=bsq*hipr1(40)/hipr1(31)/0.1d0

c...1+1 simulation
      is11=0
      if(mstc(17).eq.1) is11=1
      if(k(2,jp).eq.92.or.k(2,jt).eq.92)is11=2
      if(k(1,jp).eq.4.or.k(1,jt).eq.4)is11=3

      if(is11.ne.0) then
        gs=1.0d0-exp(-(ssoft+sjet)/s0*romg(r2))
        rantot=rn(0)*gs

      else

c...Total at b=0
c       gstot_0=2.0d0*(1.0d0-exp(-(ssoft+sjet)/s0/2.0d0*romg(0.0d0)))
c       r2=r2/gstot_0
c...Inel.
        gs=1.0d0-exp(-(ssoft+sjet)/s0*romg(r2))
c...Total
        gstot=2.0d0*(1.0d0-sqrt(1.0d0-gs))
        rantot=rn(0)*gstot
c...Elastic
        if(rantot.gt.gs) then
          jamisjet=-1
          return
        endif

      endif

      jamisjet=0

      nop=1
c...When IHPR2(8)=0 no jets are produced
      if(ihpr2(8).eq.0 .and. ihpr2(3).eq.0) return

c...Energy too small.
      nop=3
      if(pcm(5).lt.parc(71)) return

c...Hard eikonal
      tt=sjet*romg(r2)/s0
c...Soft eikonal
      tts=ssoft*romg(r2)/s0

c...This is the probability for no jet production
      nop=5
      if(ihpr2(8).gt.0) then
        probj=exp(-tt)*(1.0d0-exp(-tts))
        if(rantot.lt.probj) return
      endif

c...Determine number of mini jet production
      n_jet=0
      xr=-dlog(exp(-tt)+rn(0)*(1.0d0-exp(-tt)))
111   n_jet=n_jet+1
      xr=xr-dlog(rn(0))
      if(xr.lt.tt) go to 111

      n_jet=min(n_jet,ihpr2(8))
      if(ihpr2(8).lt.0)  n_jet=abs(ihpr2(8))
      jamisjet=n_jet

      end



c***********************************************************************

      subroutine jamcmom(i1,i2,pv,kfv,k9v,kfq,srt)

c...Construct energy-momentum for collisions.
c...pv:momenta
c...kfv:kf flavor codes
c...k9v: 3times baryon numbers
c..srt: invariant mass of the system

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      dimension pv(2,5),kfv(2),k9v(2),kfq(2,2)

c...Loop over proj. and targ.
      do jt=1,2
        if(jt.eq.1) ii=i1
        if(jt.eq.2) ii=i2
        k1=k(1,ii)
        kf=k(2,ii)
        kfv(jt)=kf
        k9v(jt)=k(9,ii)
        if(mod(abs(k1),10).ne.4) then
c......String from soft scattering.
          if(kf.eq.92) then
            do j=1,4
            pv(jt,j)=vq(j,ii)+vq(j+5,ii)
            end do
            kfq(jt,1)=kq(1,ii)
            kfq(jt,2)=kq(2,ii)
            call kfcnst(kq(1,ii),kq(2,ii),kft,0d0)
            if(kft.eq.0)then
              call jamerrm(30,0,'(jamcmom:)kf=0 at kf=92')
            endif
            kfv(jt)=kft
          else
            do j=1,4
            pv(jt,j)=p(j,ii)
            end do
            call attflv(kf,kfq(jt,1),kfq(jt,2))
          endif

c...Remnant from hard scattering.
        else
          jp1=k(10,ii)
          jp2=k(11,ii)
          kfq(jt,1)=k(2,jp1)
          kfq(jt,2)=k(2,jp2)

c....In the case of gluon-gluon system.
          if(kfq(jt,1).eq.21.or.kfq(jt,2).eq.21) then
            kf01=21
          else
            call kfcnst(kfq(jt,1),kfq(jt,2),kf01,0.0d0)
            if(kf01.eq.0)then
              call jamerrm(30,0,'(jamcmom:)kf01=0 at k1=4')
            endif
          endif
          do j=1,4
          pv(jt,j)=p(j,jp1)+p(j,jp2)
          end do
          kqc1=jamcomp(kfq(jt,1))
          kqc2=jamcomp(kfq(jt,2))
          kfv(jt)=kf01
          k9v(jt)=kchg(kqc1,6)*isign(1,kfq(jt,1))
     $           +kchg(kqc2,6)*isign(1,kfq(jt,2))
        endif

        pv(jt,5)=sqrt(max(0d0,pv(jt,4)**2-pv(jt,1)**2-pv(jt,2)**2-
     $  pv(jt,3)**2))

c...2020/2/4
        if(mstd(91).eq.1.and.MF_on(ii).eq.1) then
          pv(jt,1)=pv(jt,1)-potv(1,ii)
          pv(jt,2)=pv(jt,2)-potv(2,ii)
          pv(jt,3)=pv(jt,3)-potv(3,ii)
          pv(jt,4)=sqrt(pv(jt,5)**2+pv(jt,1)**2+pv(jt,2)**2+pv(jt,3)**2)
        endif

      end do

      srt=sqrt(max(0d0,(pv(1,4)+pv(2,4))**2-(pv(1,1)+pv(2,1))**2
     $-(pv(1,2)+pv(2,2))**2-(pv(1,3)+pv(2,3))**2))


      end

c***********************************************************************

      subroutine jamdeut

c...Purpose: to determine deuteron nuclear cluster by coalescence.
c...Last modified by Yuichi Hirata
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      dimension p1(5),p2(5),r1(5),r2(5)

c..8 choices are available
c..See Readme in Data-From-Phase/Readme
c..mcoal=5 and 8 seems to be good( give almost same results)
c..****************
c..need to multiply weid(i) on multiplicity of deuteron
c..****************
c     data mcoal/5/
      mcoal=mstc(158)
      mstd(65)=mstd(65)+1
      nclust=0

c...Define cluster distance.
      if((mcoal.eq.1).or.(mcoal.eq.2)) then
        rcc2=parc(151)**2
        pcc2=parc(152)**2
      else
        wdeut=parc(153)
c       wdeut=1.71d0
c       wdeut=1.76d0
c       wdeut=1.56d0
      endif
 
      mentry=0
      do 100 i1=2,nv
        kf1=k(2,i1)
        if(kf1.eq.1001001000) k(1,i1)=11
        if(kf1.eq.2112.or.kf1.eq.2212) k(1,i1)=1
        if(kf1.ne.2212.and.kf1.ne.2112) goto 100
      do 110 i2=1,i1-1
        kf2=k(2,i2)
        if((kf1.eq.2212.and.kf2.eq.2112)
     $    .or.(kf1.eq.2112.and.kf2.eq.2212))then
        else
          goto 110
        endif
        pcm1=p(1,i1)+p(1,i2)
        pcm2=p(2,i1)+p(2,i2)
        pcm3=p(3,i1)+p(3,i2)
        pcm4=p(4,i1)+p(4,i2)
        srt=sqrt(pcm4**2-pcm1**2-pcm2**2-pcm3**2)
        bex=pcm1/pcm4
        bey=pcm2/pcm4
        bez=pcm3/pcm4
        gam=pcm4/srt
        do j=1,5
         p1(j)=p(j,i1)
         p2(j)=p(j,i2)
         r1(j)=r(j,i1)
         r2(j)=r(j,i2)
        end do
        call jamrobo(0d0,0d0,-bex,-bey,-bez,gam,p1(1),p1(2),p1(3),p1(4))
        call jamrobo(0d0,0d0,-bex,-bey,-bez,gam,p2(1),p2(2),p2(3),p2(4))
        call jamrobo(0d0,0d0,-bex,-bey,-bez,gam,r1(1),r1(2),r1(3),r1(4))
        call jamrobo(0d0,0d0,-bex,-bey,-bez,gam,r2(1),r2(2),r2(3),r2(4))
        t=max(r1(4),r2(4))
        t1=t-r1(4)
        t2=t-r2(4)
        r1(4)=t
        r2(4)=t
        do j=1,3
         r1(j)=r1(j)+t1*p1(j)/p1(4)
         r2(j)=r2(j)+t2*p2(j)/p2(4)
        end do
        rr2=(r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2
        pp2=(p2(1)-p1(1))**2+(p2(2)-p1(2))**2+(p2(3)-p1(3))**2


        ideut=0
c..spin factor
        fac1=3.0d0/4.0d0
c..isospin factor
        fac2=1.0d0/2.0d0
        if(mcoal.eq.1)    then
          if((rr2.le.rcc2.and.pp2.le.pcc2).and.rn(0).le.fac1)then
            ideut=1
            weid(i1)=1.0d0
          endif
        elseif(mcoal.eq.2)then
          if(rr2.le.rcc2.and.pp2.le.pcc2)then
            ideut=1
            weid(i1)=fac1
          endif
        elseif(mcoal.eq.3)then
          prob1=8d0*exp(-rr2/wdeut**2-pp2*(wdeut/paru(3))**2)
          if(rn(0).le.fac1*prob1)then
            ideut=1
            weid(i1)=1.0d0
          endif
        elseif(mcoal.eq.4)then
          prob1=8d0*exp(-rr2/wdeut**2-pp2*(wdeut/paru(3))**2)
          if(rn(0).le.prob1)then
            ideut=1
            weid(i1)=fac1
          endif
        elseif(mcoal.eq.5)then
          prob2=8d0*exp(-rr2/wdeut**2-pp2*(wdeut/paru(3)/2d0)**2)
          if(rn(0).le.prob2)then
            ideut=1
            weid(i1)=fac1
          endif
        elseif(mcoal.eq.6)then
          prob2=8d0*exp(-rr2/wdeut**2-pp2*(wdeut/paru(3)/2d0)**2)
          if(rn(0).le.prob2)then
            ideut=1
            weid(i1)=fac1*fac2
          endif
        elseif(mcoal.eq.7)then
          prob2=8d0*exp(-rr2/wdeut**2-pp2*(wdeut/paru(3)/2d0)**2)
          if(rn(0).le.fac1*prob2)then
            ideut=1
            weid(i1)=1.0d0
          endif
        elseif(mcoal.eq.8)then
          prob2=8d0*exp(-rr2/wdeut**2-pp2*(wdeut/paru(3)/2d0)**2)
          rkx=(p2(1)-p1(1))/2.0d0/paru(3)
          rky=(p2(2)-p1(2))/2.0d0/paru(3)
          rkz=(p2(3)-p1(3))/2.0d0/paru(3)
          rrx=r2(1)-r1(1)
          rry=r2(2)-r1(2)
          rrz=r2(3)-r1(3)
          cs=cos(rkx*rrx+rky*rry+rkz*rrz)
          if(rn(0).le.prob2)then
            ideut=1
            weid(i1)=fac1*fac2*2.0d0*cs*cs
          endif
        endif

        if(ideut.eq.1) then
          call jamrobo(0d0,0d0,bex,bey,bez,gam,r1(1),r1(2),r1(3),r1(4))
          call jamrobo(0d0,0d0,bex,bey,bez,gam,r2(1),r2(2),r2(3),r2(4))
          mentry=mentry+1 
          icoll(1,mentry)=i1
          icoll(2,mentry)=i2
          coll(1,mentry)=r1(4)
          coll(2,mentry)=r2(4)
          coll(3,mentry)=max(r1(4),r2(4))
        endif
 110  continue
 100  continue

 1000 continue
      if(mentry.le.0) goto 2000

c...Find deuteron.
      tcol=100000.0d0 
      jent=0
      do ient=1,mentry
        t=coll(3,ient)
        if(t.le.tcol) then
           tcol=t
           jent=ient
        endif
      end do
      if(jent.eq.0)then
        goto 2000
      endif 

c...Form deuteron.
      nclust=nclust+1
      mstd(30)=mstd(30)+1
      i1=icoll(1,jent)
      i2=icoll(2,jent)
      t=coll(3,jent)
      dt1=t-coll(1,jent)
      dt2=t-coll(2,jent)
      do j=1,3
       r(j,i1)=r(j,i1)+dt1*p(j,i1)/p(4,i1)
       r(j,i2)=r(j,i2)+dt2*p(j,i2)/p(4,i2)
      end do

      k(1,i2)=21
      k(1,i1)=5
      k(2,i1)=1000000*1+1000*1+1000000000
      k(3,i1)=0
      k(4,i1)=0
      k(5,i1)=0
      k(6,i1)=0
      k(7,i1)=abs(k(7,i1))+abs(k(7,i2))
      k(9,i1)=6
      p(1,i1)=p(1,i1)+p(1,i2)
      p(2,i1)=p(2,i1)+p(2,i2)
      p(3,i1)=p(3,i1)+p(3,i2)
      p(4,i1)=p(4,i1)+p(4,i2)
      p(5,i1)=sqrt(p(4,i1)**2-p(1,i1)**2-p(2,i1)**2-p(3,i1)**2)
      r(1,i1)=(r(1,i1)+r(1,i2))/2.0d0
      r(2,i1)=(r(2,i1)+r(2,i2))/2.0d0
      r(3,i1)=(r(3,i1)+r(3,i2))/2.0d0
      r(4,i1)=t

      
      mstd(30)=mstd(30)+1

c...Remove collisions including i1,i2.
      idelte=0
      do 300 jjj=1,mentry
        if((icoll(1,jjj).eq.i1).or.(icoll(2,jjj).eq.i1).or.
     $     (icoll(1,jjj).eq.i2).or.(icoll(2,jjj).eq.i2)) then
           idelte = idelte + 1
           go to 300
        end if
        jj=jjj-idelte
        if(idelte .gt. 0) then
          icoll(1,jj)=icoll(1,jjj)
          icoll(2,jj)=icoll(2,jjj)
          coll(1,jj)=coll(1,jjj)
          coll(2,jj)=coll(2,jjj)
          coll(3,jj)=coll(3,jjj)
        end if
  300 continue
      mentry=mentry-idelte
      goto 1000
 2000 continue


c...Save cluster number.
      mstd(66)=nclust

c...Remove unwanted entries.
      call jamedit
 
      end

c***********************************************************************

      subroutine jamglaub(icon,msel)

c....Perform AA collision according to the Glauber theory.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      parameter(mxchan=30)
      dimension pv(2,5),kfv(2),k9v(2),kfq(2,2),sigin(mxchan),indd(100)
      common/jamhrdev1/lead(2)
      save /jamhrdev1/

      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 

c...Get collision number.
      icon=0
      mentry=0
      do 100 jt=1,mstd(5)
      do 110 jp=1+mstd(5),mstd(2)+mstd(5)
        kf1=k(2,jp)
        kf2=k(2,jt)
        ibar1=k(9,jp)
        ibar2=k(9,jt)
        em1=p(5,jp)
        em2=p(5,jt)
c       b2=(r(1,jp)-r(1,jt))**2+(r(2,jp)-r(2,jt))**2
        t01=r(4,jp)
        t02=r(4,jt)
        em1sq=p(4,jp)**2-p(1,jp)**2-p(2,jp)**2-p(3,jp)**2
        em2sq=p(4,jt)**2-p(1,jt)**2-p(2,jt)**2-p(3,jt)**2
        dt=t02-t01
        dx=r(1,jt)-r(1,jp)
        dy=r(2,jt)-r(2,jp)
        dz=r(3,jt)-r(3,jp)
        dx12=dt**2-dx**2-dy**2-dz**2
        dxp1=dt*p(4,jp)-dx*p(1,jp)-dy*p(2,jp)-dz*p(3,jp)
        dxp2=dt*p(4,jt)-dx*p(1,jt)-dy*p(2,jt)-dz*p(3,jt)
        dp12=p(4,jp)*p(4,jt)
     $   -p(1,jp)*p(1,jt)-p(2,jp)*p(2,jt)-p(3,jp)*p(3,jt)
        dn=dp12*dp12-em1sq*em2sq
        if(dn.lt.1d-5) goto 100
        tcol1=t01-p(4,jp)*(dxp1*em2sq-dxp2*dp12)/dn
        tcol2=t02+p(4,jt)*(dxp2*em1sq-dxp1*dp12)/dn
        b2=-dx12-(dxp1**2*em2sq+dxp2**2*em1sq-2.d0*dxp1*dxp2*dp12)/dn
        icltyp=jamcltyp(kf1,kf2,ibar1,ibar2)
        srt=sqrt((p(4,jp)+p(4,jt))**2-(p(1,jp)+p(1,jt))**2
     $    -(p(2,jp)+p(2,jt))**2-(p(3,jp)+p(3,jt))**2)
        if(srt.le.em1+em2+0.0001) goto 110
        pr=pawt(srt,em1,em2)
        call jamcross(1,icltyp,srt,pr,kf1,kf2,em1,em2,
     $                 sig,sigel,sigin,mchanel,mabsrb,isoft,icon)
        if(mstd(11).eq.2) then
          bcmax2=0.1d0*(sig-sigel)/paru(1)
          if(b2.gt.bcmax2) goto 110
        else
          bcmax2=0.1d0*sig/paru(1)
          if(b2.gt.bcmax2) goto 110
        endif
        mentry=mentry+1
        icoll(1,mentry)=jp
        icoll(2,mentry)=jt
        coll(1,mentry)=(tcol1+tcol2)/2d0
        coll(2,mentry)=tcol1
        coll(3,mentry)=tcol2
        coll(4,mentry)=sig
        coll(5,mentry)=sigel
        coll(6,mentry)=b2
        kq(2,jp)=1
        kq(2,jt)=1
  110 continue
  100 continue

      if(mstd(11).eq.2.and.mentry.eq.0) then
        write(6,*)'(jamglaub:) mentry=0:no collision? b=',pard(2)
        write(6,*)'parc(4) bcmax2',parc(4),sqrt(bcmax2),sqrt(b2)
        icon=1
        return
      endif

c...Save predicted initial collision number and participants.
        mstd(79)=mentry
        mstd(80)=0
        do i=1,nv
          if(kq(2,i).eq.1) then
            mstd(80)=mstd(80)+1
            kq(2,i)=0
          endif
        end do

      if(msel.eq.0) return

c...Do multiple collision.
      do 120 ient=1,mentry

        i1=icoll(1,ient)
        i2=icoll(2,ient)
        if(i1.le.0.or.i2.le.0) goto 120
        if(k(1,i1).eq.0.or.k(1,i1).gt.10) goto 120
        if(k(1,i2).eq.0.or.k(1,i2).gt.10) goto 120

        if(mstc(6).eq.-2) then
          if(mod(abs(k(1,i1))/10,10).ne.2) goto 120
          if(mod(abs(k(1,i2))/10,10).ne.2) goto 120
        else if(mstc(6).eq.-3) then
          if(k(1,i1).lt.0) then
            qcnum=mod(abs(k(1,i1))/10,10)
            ibar1=k(9,i1)
            if(abs(ibar1).eq.3) qnum1=3.d0
            if(ibar1.eq.0) qnum1=2.d0
            if(qcnum.eq.3) qcnum=2
            facq1=qcnum/qnum1
            if(rn(0).gt.facq1) goto 120
          endif
          if(k(1,i2).lt.0) then
            qcnum=mod(abs(k(1,i2))/10,10)
            ibar2=k(9,i2)
            if(abs(ibar2).eq.3) qnum2=3.d0
            if(ibar2.eq.0) qnum2=2.d0
            if(qcnum.eq.3) qcnum=2
            facq2=qcnum/qnum2
            if(rn(0).gt.facq2) goto 120
          endif
        endif

        kf1=k(2,i1)
        kf2=k(2,i2)
        em1=p(5,i1)
        em2=p(5,i2)

        call jamcmom(i1,i2,pv,kfv,k9v,kfq,srt)


        if(srt.lt.pv(1,5)+pv(2,5)+0.00001d0) goto 120

        inucl1=0
        inucl2=0
        cutoff=pv(1,5)+pv(2,5)

        kfa1=abs(kfv(1))
        if(kfa1.eq.2112 .or.kfa1.eq.2212 .or.kfa1.eq.3122)inucl1=1
        if(inucl1.eq.1.and.kfv(1).lt.0) inucl1=-1
        kfa2=abs(kfv(2))
        if(kfa2.eq.2112 .or.kfa2.eq.2212.or.kfa2.eq.3122)inucl2=1
        if(inucl2.eq.1.and.kfv(2).lt.0) inucl2=-1
        if(inucl1*inucl2.eq.1) then
          cutoff=cutoff+parc(38)
        endif

        if(srt.le.cutoff) goto 120

        pare(4)=coll(4,ient)
        pare(5)=coll(5,ient)
        pare(6)=coll(6,ient)    ! bsq

        dt1=coll(2,ient)-r(4,i1)
        dt2=coll(3,ient)-r(4,i2)
        do j=1,3
        if(abs(k(7,i1)).ne.1) r(j,i1)=r(j,i1)+dt1*p(j,i1)/p(4,i1)
        if(abs(k(7,i2)).ne.1) r(j,i2)=r(j,i2)+dt2*p(j,i2)/p(4,i2)
        end do

        pard(1)=coll(1,ient)    ! Collision time.
        r(4,i1)=coll(2,ient)
        r(4,i2)=coll(3,ient)
        r(5,i1)=r(4,i1)
        r(5,i2)=r(4,i2)

c...Save particle information.
        call jamsave(1,1,i1)
        call jamsave(1,2,i2)

c...Scatter two-particles.
        call jamscatt


c...Count collision.
        call jamanacl(2)

c...Collision forbidden.
        if(mste(1).le.0) then
          call jamsave(2,1,mste(21))
          call jamsave(2,2,mste(23))

c...Hard collision. i1 and i2 particles are moved to the different
c...slot lead(jt), so we need to change the collision array.
        else if(mste(1).eq.6) then
          do jent=ient+1,mentry
            do jt=1,2
            do jj=1,2
            if(jj.eq.1) ii=i1
            if(jj.eq.2) ii=i2
            if(icoll(jt,jent).eq.ii) then
              icoll(jt,jent)=lead(jt)
            endif
            end do
            end do
          end do
        endif


cTABARA
c       call ttchk(indd,0)
c....Print collision information.
        if(mstc(8).ge.2) call jamprcl(indd,0)
        if(mste(1).le.0) goto 120

c...Fragment strings.
        if(mod(abs(mstc(6)),10).eq.2) then
          if(mod(abs(k(1,i1)),10).eq.3) then
            call jamsave(1,1,i1)
            call jamdec(i1,indd,nadd,icon)
          endif
          if(mod(abs(k(1,i2)),10).eq.3) then
            call jamsave(1,1,i2)
            call jamdec(i2,indd,nadd,icon)
          endif
        endif


 120  continue

c...Fragment strings and resonances.
      if(mod(abs(mstc(6))/10,10).eq.1) call jamfdec

      end

c***********************************************************************
c...LPC
      subroutine jamfpath

c....Perform hA collision according to the mean free path argument.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'

c...Commonblock for LPC.
      common /jamlpc1/wa0,wr0,ww,rho0,ext,idens,icscde

      parameter(mxchan=30)
      real*8 jamdtim
      dimension sigin(mxchan),indd(100)
      dimension r1(5),p1(5),p2(2,5),sigtp(2),sigelp(2)
      logical docoll
      data pnorm,dt,bnde/20d0,0.2d0,0.0d0/
c...Low velocity cut below which leading particle quit scattering.
      data vlcut/0.00001d0/
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a) 


c...Initial values: global time,collision counter,etc.
      pard(1)=0.0d0
      ncolp=0
      ncoln=0
      k(5,1)=0
      r(1,1)=pard(2)
      r(2,1)=0d0
      r(3,1)=-sqrt(abs(pard(4)**2-pard(2)**2))
      r(4,1)=0d0
      r(5,1)=0d0
      do j=1,4
        v(j,1)=r(j,1)
      end do
      v(5,1)=1d+35
      ip=1
      ncoll=0

c....Leading cascade starts.
4000  continue
      ncoll=ncoll+1
      if(ncoll.ge.1000) then
        call jamerrm(30,0,'(jamfpath:) infinit collision???')
      endif

      kf1=k(2,ip)
      kf1a=abs(kf1)
c...Exclude photons, leptons...
      if((kf1a.gt.10.and.kf1a.le.100).and.kf1.ne.21) goto 6000
      kc1=jamcomp(kf1)
      ibar1=isign(kchg(jamcomp(kf1),6),kf1)
      do j=1,5
        p1(j)=p(j,ip)
        r1(j)=r(j,ip)
      end do
      pp1=sqrt(p1(1)**2+p1(2)**2+p1(3)**2)
      v1=pp1/p1(4)
c...Too low velocity: quit cascading.
      if(v1.lt.vlcut) goto 6000
      rr1=sqrt(r1(1)**2+r1(2)**2+r1(3)**2)
      docoll=.true.

c...In case of |r1| > rad. Determine cascade is continue?
      rr1=sqrt(r1(1)**2+r1(2)**2+r1(3)**2)
      if(rr1.gt.pard(4)*ext) goto 6000

c...Calculate proton and neutron density.
      if(idens.eq.1) then
        vol=4*paru(1)*pard(4)**3/3
        rho=dble(mstd(5))/vol
        if(rr1.gt.pard(4)) then
           rho=0.0d0
        endif
      else
        rho=rho0*(1+ww*(rr1/wr0)**2)/(1+exp((rr1-wr0)/wa0))
      endif
      rhop=rho*mstd(6)/mstd(5)
      rhon=rho-rhop

c...Loop over neuton and proton, picking up nucleons from target nucleus.
      do jt=1,2

      if(jt.eq.1) then
        kf2=2112
        p2(1,5)=parc(24)
      else
        kf2=2212
        p2(2,5)=parc(25)
      endif
      ibar2=3
      icltyp=jamcltyp(kf1,kf2,ibar1,ibar2)

c....Generate Fermi momentum.
      call jamfgas(p2(jt,1),p2(jt,2),p2(jt,3),rhon,pf)

      p2(jt,4)=sqrt(p2(jt,5)**2+p2(jt,1)**2+p2(jt,2)**2+p2(jt,3)**2)
      ss=(p1(4)+p2(jt,4))**2-
     a       (p1(1)+p2(jt,1))**2-(p1(2)+p2(jt,2))**2-(p1(3)+p2(jt,3))**2
      srt=sqrt(ss)-bnde
      pr=pawt(srt,p1(5),p2(jt,5))

c....Determine max. cross section and max. impact par.
c....as well as low energy cutoff
      cutoff=p1(5)+p2(jt,5)
      if(ibar1*ibar2.eq.9) then
        inucl1=0
        inucl2=0
        if(kf1.eq.2112.or.kf1.eq.2212.or.kf1.eq.3122)inucl1=1
        if(kf1.eq.-2112.or.kf1.eq.-2212.or.kf1.eq.-3122)inucl1=-1
        if(kf2.eq.2112.or.kf2.eq.2212.or.kf2.eq.3122)inucl2=1
        if(kf2.eq.-2112.or.kf2.eq.-2212.or.kf2.eq.-3122)inucl2=-1
        if(inucl1*inucl2.eq.1) cutoff=p1(5)+p2(jt,5)+parc(38)
      endif

      if( (kf1.eq.111.or.abs(kf1).eq.211)
     $   .and.(kf2.eq.111.or.abs(kf2).eq.211) ) cutoff=cutoff+parc(38)

c...Low energy cutt off  i.e. Pauli block.
      if(srt.lt.cutoff.or.pr.lt.0.001d0) then
        mstd(52)=mstd(52)+1
        go to 6000
      endif

c...Get total cross section between leading particle + neutron(proton).
      pare(3)=0.0d0
      call jamcross(1,icltyp,srt,pr,kf1,kf2,p1(5),p2(jt,5),
     $          sigtp(jt),sigelp(jt),sigin,mchanel,mabsrb,isoft,icon)

      end do

c...Get decay width of leading particle.
      wid=jamdtim(0,kf1,kc1,1,p1(5),p1(4))


c...According to mean free path, estimate collision prob..
500   continue

c...Random number for collision.
      ran1=rn(0)
      prob1=rhon*sigtp(1)/10
      prob2=rhop*sigtp(2)/10
      gam1=p1(4)/p1(5)
      probd=wid/(v1*paru(3)*gam1)

c...Reduce cross section for hadron that has virtual quarks.
      if(k(1,ip).lt.0) then
        qcnum=mod(abs(k(1,ip))/10,10)
        if(abs(ibar1).eq.3) qnum1=3.d0
        if(ibar1.eq.0) qnum1=2.d0
        if(qcnum.eq.3) qcnum=2
        facq1=qcnum/qnum1
c       if(rn(0).gt.facq1) goto 120
        prob1=prob1*facq1
        prob2=prob2*facq1
      endif

c========= cascade procedure 1 =======================================
c....Uniform matter case.
      if(icscde.eq.1) then
        probt=prob1+prob2+probd

c...Normal end.
        if(probt.lt.1.0d-10) goto 6000

        patht=1.0d0/probt  !Mean free path
        path=-patht*log(max(1.0d-15,1.0d0-rn(0)))
        r1(1)=r1(1)+path*p1(1)/pp1
        r1(2)=r1(2)+path*p1(2)/pp1
        r1(3)=r1(3)+path*p1(3)/pp1
        r1(4)=r1(4)+path*p1(4)/pp1
        r1(5)=r1(4)

c========= cascade procedure 2 =======================================
      else if(icscde.eq.2) then

       probt=prob1+prob2+probd
       if(probt.lt.1.0d-10) goto 6000

       patht=1.0d0/probt
       dz=patht/pnorm
       prob1=prob1*dz
       prob2=prob2*dz
       probd=probd*dz
       probt=probt*dz

       if(ran1.gt.probt) then
	 docoll=.false.
       else
         dz=ran1*probt
       endif
       r1(1)=r1(1)+dz*p1(1)/pp1
       r1(2)=r1(2)+dz*p1(2)/pp1
       r1(3)=r1(3)+dz*p1(3)/pp1
       r1(4)=r1(4)+dz*p1(4)/pp1
       r1(5)=r1(4)
       probt=1

c== cascade procedure 3 by Ohnishi ============================
      else

        prob1=prob1*v1
        prob2=prob2*v1
        probd=probd*v1
        probt=prob1+prob2+probd
        if(probt*dt/0.5d0.ge.1.0d0) then
          write(6,*)'time step is large probt*dt=',probt*dt
          write(6,*)'sig1 sig2=',sigtp(1),sigtp(2)
          write(6,*)'kf1',kf1
        endif

c...No Reaction in this step.
        if(exp(-probt*dt).gt.rn(0)) then
          r1(1)=r1(1)+dt*p1(1)/p1(4)
          r1(2)=r1(2)+dt*p1(2)/p1(4)
          r1(3)=r1(3)+dt*p1(3)/p1(4)
          r1(4)=r1(4)+dt
          r1(5)=r1(4)
          docoll=.false.
        endif


c...Reaction occurs in this step. When will it occur ?
c... --> rn(0)=  \int_0^t          dt \exp(-\Gamma t)
c...             / \int_0^{\Delta t} dt \exp(-\Gamma t)

        dt1=-log(1.0d0-rn(0)*(1.0d0-exp(-probt*dt)))/probt
        r1(1)=r1(1)+dt1*p1(1)/p1(4)
        r1(2)=r1(2)+dt1*p1(2)/p1(4)
        r1(3)=r1(3)+dt1*p1(3)/p1(4)
        r1(4)=r1(4)+dt
        r1(5)=r1(4)

      endif
c==============================================================

c...Update position of leading particle.
      do j=1,5
        r(j,ip)=r1(j)
      end do

c...No collision in this step.
      if(.not.docoll) goto 4000

c...In case of |r1| > rad. Determine cascade is continue?
      rr1=sqrt(r1(1)**2+r1(2)**2+r1(3)**2)
      if(rr1.gt.pard(4)*ext) goto 6000

c...Now we determine whether collision partner is proton or neutron?
c...if particle is resonance, decay will be occurred. 

      ireac=0
c...Collision with neutron
      if(prob1.gt.probt*ran1) then
        ireac=1
        kf2=2112
c       rho1=rho*(mstd(5)-mstd(6))/mstd(5)

c...Collision with proton
      elseif(prob1+prob2.gt.probt*ran1) then
        ireac=2
        kf2=2212
c       rho1=rho*mstd(6)/mstd(5)

c...decay event
      elseif(prob1+prob2+probd.gt.probt*ran1) then
        ireac=3
      else
        write(6,*) '* Mysterious Branch'
        write(6,*) 'itypes=',kf1,kf2
        write(6,*) 'prb1,prb2,prb1+prob2,prbt='
        write(6,*) prob1,prob2,prob1+prob2,probt
        write(6,*)'sigtot1 sigtot2= ', sigtp(1),sigtp(2)
        goto 6000
      endif


      if(ireac.eq.1.or.ireac.eq.2) then

        nv=nv+1
        nbary=nbary+1
        if(nv.gt.mxv) call jamerrm(30,0,'(jamfpath:)no more memory mxv')
        i2=nv
        call jamzero(i2)
        k(1,i2)=1
        k(2,i2)=kf2
        k(9,i2)=3
        do j=1,5
          p(j,i2)=p2(ireac,j)
          r(j,i2)=r1(j)
          v(j,i2)=r1(j)
        end do
        v(5,i2)=1d+35

c...Save particle information.
        call jamsave(1,1,ip)
        call jamsave(1,2,i2)

c...Scatter two-particles.
        pare(4)=sigtp(ireac)
        pare(5)=sigelp(ireac)
        pare(6)=0.0d0
        call jamscatt

c...Count collision.
        call jamanacl(2)

c...Collision forbidden.
        if(mste(1).le.0) then
          call jamsave(2,1,mste(21))
          call jamsave(2,2,mste(23))
        endif

cTABARA
c       call ttchk(indd,0)
c....Print collision information.
        if(mstc(8).ge.2) call jamprcl(indd,0)


c...In the case of absorption.
        if(mste(25).eq.0.and.mste(27).eq.0) goto 6100
        ip=mste(25)
        if(mste(25).eq.0.and.mste(27).ge.1) ip=mste(27)

        if(mste(25).ge.1) k(5,mste(25))=0
        if(mste(27).ge.1) k(5,mste(27))=0
        if(ireac.eq.1) ncoln=ncoln+1
        if(ireac.eq.2) ncolp=ncolp+1

      else

c...Set collision type.
        mste(2)=-1
        call jamsave(1,1,ip)
        call jamdec(ip,indd,nadd,icon)
c...Kinetically can not decay, so should be scattered.
        if(icon.eq.1) then
          goto 500
        endif

cTABARA
c       call ttchk(indd,nadd)
c...Print information after decay.
        if(mstc(8).ge.2) call jamprcl(indd,nadd)

        do i=1,nadd
          k(5,indd(i))=0
        end do
      endif

      pard(1)=r1(4)
c...Go to next collision of ip.
      goto 4000

c...Search for next particle which will collied.
 6000 continue
      k(5,ip)=1
 6100 continue
      if(ncolp.ge.mstd(6)) return
      if(ncoln.ge.mstd(5)-mstd(6)) return
      do i=1,nv
        if(k(1,i).le.10.and.k(5,i).eq.0) then
          ip=i
          goto 4000
        endif
      end do

      end

c***********************************************************************
c...LPC

      subroutine jaminil

c...Initialize some values for leading particle cascade.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      parameter(mxchan=30)
      dimension sigin(mxchan)
      parameter (mxa=5)
c...Commonblock for LPC.
      common /jamlpc1/wa0,wr0,ww,rho0,ext,idens,icscde

c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a)

c...  icscde        icscde=1: mean free path prescription
c...=mod(mstc(6)/10,10)    2: mean free path  1 dz
c...                       3: mean free path  2 dt
c...  idens         idens= 1: nuclear matter density
c...=mod(mstc(6),10)       2: Fermi-type density for nucleus

      icscde=mod(abs(mstc(6))/10,10)
      idens=mod(abs(mstc(6)),10)
      rho0=0.165d0

c...proj. + proton cross section.
      srt=pard(16)
      kf1=mstd(1)
      kf2=2212
      em1=pard(34)
      em2=parc(25)
      pr=pawt(srt,em1,em2)
      pare(3)=0.0d0
      ibar1=isign(kchg(jamcomp(kf1),6),kf1)
      ibar2=kchg(jamcomp(kf2),6) 
      icltyp=jamcltyp(kf1,kf2,ibar1,ibar2)
      call jamcross(1,icltyp,srt,pr,kf1,kf2,em1,em2,
     $              sigtot1,sigel,sigin,mchanel,mabsrb,ijet,icon)
c...proj. + neutron cross section.
      srt=pard(16)
      kf1=mstd(1)
      kf2=2112
      em1=pard(34)
      em2=parc(24)
      pr=pawt(srt,em1,em2)
      pare(3)=0.0d0
      ibar1=isign(kchg(jamcomp(kf1),6),kf1)
      ibar2=kchg(jamcomp(kf2),6) 
      icltyp=jamcltyp(kf1,kf2,ibar1,ibar2)
      call jamcross(1,icltyp,srt,pr,kf1,kf2,em1,em2,
     $              sigtot2,sigel,sigin,mchanel,mabsrb,ijet,icon)

c...Average Elementary Cross Section.
      sigav=(sigtot2*(mstd(5)-mstd(6))+sigtot1*mstd(6))/mstd(5)

c...Get min. bias impact paramter.
c...rho=rho0*(1+ww*(r/wr0)**2)/(1+exp((r-wr0)/wa0))
      if(idens.eq.1) then   ! uniform matter density
	wr0=0.0d0
	wa0=0.0d0
	ww=0.0d0
	ext=1.0d0
        if(mstd(5).le.16) then
          rad=1.2d0*dble(mstd(5))**(1.0d0/3.0d0)
        else
c         rad=1.12d0*dble(mstd(5))**(1.0d0/3.0d0)
          rad=1.2d0*dble(mstd(5))**(1.0d0/3.0d0)
        endif
      else
	ext=1.5d0
        call jamdnsf(mstd(5),mstd(6),wa0,wr0,ww,idens)
	rrmax=0.001d0
        rad=wa0*log(rho0/rrmax-1)+wr0
	if(mstd(5).le.16) rad=3.5d0
        if(icscde.ge.4)
     $   rad=1.19d0*mstd(5)**0.3333d0-1.61d0*mstd(5)**(-0.3333d0)
     $  +sqrt(sigav/paru(1)/10)
      endif

      if(icscde.eq.1) then
	drad=0.0d0
      else if(icscde.le.3) then
        if(idens.eq.1) then
	   drad=0.0d0
	else
	   rrmax=1.d-4*0.165d0
	   rho0=0.165d0
           rad2=wa0*log(rho0/rrmax-1)+wr0
	   if(mstd(5).eq.12) rad2=5.5d0
	   drad=rad2-rad
	endif
      else
        drad=sqrt(sigav/paru(1)/10)+0.5d0
      endif

      pard(4)=rad+drad
c     write(6,*)'density radii drad=',rad,drad,pard(4)

c...Set leading particle.
      kc1=jamcomp(mstd(1))
      nv=1
      ibar=kchg(kc1,6)
      if(ibar.eq.0) then
        nmeson=mstc(5)
        nbary=0
      else
        nbary=mstc(5)
        nmeson=0
      endif

      call jamzero(1)
      k(1,1)=1
      k(2,1)=mstd(1)
      k(7,1)=1
      k(9,1)=ibar

      p(1,1)=pard(31)
      p(2,1)=pard(32)
      p(3,1)=pard(33)
      p(4,1)=sqrt(pard(34)**2+p(1,1)**2+p(2,1)**2+p(3,1)**2)
      p(5,1)=pard(34)


      end

c**********************************************************************
c...LPC

      subroutine jamfgas(px,py,pz,rhop,pf)

c...Calculate Fermi momentum using Fermi gas model.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

c.....pff=hc*(3*pi*pi)**(1.0/3.0)
      data pff/0.6104643d0/

      ifrm=mod(abs(mstc(6)),10)

      if(rhop.le.0.0) then
       pf=0.0
      else
       pf=pff*rhop**(1.0/3.0)
      endif
      pfa=pf/sqrt(5.0)

c...Local Fermi Gas
      if(ifrm.eq.1) then
c       pf=hc*(1.5*paru(1)*paru(1)*rhop)**(1.0/3.0)
        pn=pf*(rn(0))**(1.0/3.0)
        cos1=1.0-2*rn(0)
        sin1=sqrt(1.0-cos1**2)
        phi=paru(2)*rn(0)
        px=pn*sin1*cos(phi)
        py=pn*sin1*sin(phi)
        pz=pn*cos1

c...Fermi Gas
      elseif(ifrm.eq.2) then
        pf=0.27d0
        pn=pf*(rn(0))**(1.0/3.0)
        cos1=1.0-2*rn(0)
        sin1=sqrt(1.0-cos1**2)
        phi=paru(2)*rn(0)
        px=pn*sin1*cos(phi)
        py=pn*sin1*sin(phi)
        pz=pn*cos1

c...Gaussian Distribution
c     else
c       px=pfa*rann(0)
c       py=pfa*rann(0)
c       pz=pfa*rann(0)
      endif

      end

c***********************************************************************
c...LPC

      subroutine jamdnsf(n1,iz1,a,r0,w,idens)

      implicit double precision(a-h, o-z)
c...ref. Numerical Data and Functional Reaction ....
c....C.W.DeJager, et al., Nucl.Data Tables 14,479(1974)

      if(idens.eq.0) then
         w=0.0d0
         a=0.53d0
         r0=1.19d0*n1**(0.333333333d0)-1.61d0*n1**(-0.33333333d0)
      else
         if(n1.eq.12.and.iz1.eq.6) then
            w=-0.149d0
            a=0.5224d0
            r0=2.355d0
         else if(n1.eq.27.and.iz1.eq.13) then
            w=0.0d0
            a=0.52d0
            r0=3.07d0
         else if(n1.eq.63.and.iz1.eq.29) then
            w=0.0d0
            a=0.55d0
            r0=4.20d0
         else if(n1.eq.107.and.iz1.eq.47) then
            w=0.0d0
            a=0.52d0
            r0=5.12d0
         else if(n1.eq.208.and.iz1.eq.82) then
            w=0.32d0
            a=0.54d0
            r0=6.40d0
         else
            w=0.0d0
            a=0.52d0
            r0=1.19d0*n1**(0.3333d0)-1.61d0*n1**(-0.3333d0)
         endif
      endif

      end

c***********************************************************************
      subroutine jamecon(isetm)

c...isetm=1: absorption
c...isetm=2: 2-body scattering

cq        eini=pare(11)
cq        call setmom(isetm,eini,pcm,gam,em1,em2
cq   $                    ,pr,pxr,pyr,pzr,i1,i2,icon )
c.......Failed to recover energy conservation.
cq        if( icon.ne.0 .and. mstc(8).ge.3) then
cq            write(mstc(38),*)'ichanel',ichanel
cq            write(mstc(38),*)chaf(mste(22),(3-isign(1,kcp(2,1)))/2)
cq   $       ,' ',pcp(5,1),' ',chaf(mste(24),(3-isign(1,kcp(2,2)))/2)
cq   $       ,' ',pcp(5,2)
cq   $       ,' ',chaf(kc1,(3-isign(1,kf1))/2),' ',p(5,i1)
cq   $       ,' ',chaf(kc2,(3-isign(1,kf2))/2),' ',p(5,i2)
cq          mste(1)=-isetm-4
cq          return
cq        endif

      end


cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      subroutine checkCDA(i1,i2,tcol1,tcol2)
      implicit double precision(a-h, o-z)
      include 'jam1.inc'

      if(k(1,i1).ge.11.or.k(1,i2).ge.11) then
        print *,'checkCDA dead pa?',k(1,i1),k(1,i2)
        stop
      endif

      pcmx=p(1,i1)+p(1,i2)
      pcmy=p(2,i1)+p(2,i2)
      pcmz=p(3,i1)+p(3,i2)
      pcme=p(4,i1)+p(4,i2)
      bex=pcmx/pcme
      bey=pcmy/pcme
      bez=pcmz/pcme
      srt=sqrt(pcme**2 - pcmx**2 - pcmy**2 - pcmz**2)
      ga=pcme/srt
c     print *,'beta=',bex,bey,bez,ga
      pxr0=p(1,i1)
      pyr0=p(2,i1)
      pzr0=p(3,i1)
      per0=p(4,i1)
      call jamrobo(0.0d0,0.0d0,-bex,-bey,-bez,ga,pxr0,pyr0,pzr0,per0)
      phi=pjangl(pxr0,pyr0)
      call jamrobo(0d0,-phi,0d0,0d0,0d0,1.0d0,pxr0,pyr0,pzr0,per0)
      the=pjangl(pzr0,pxr0)
      call jamrobo(-the,0d0,0d0,0d0,0d0,1.0d0,pxr0,pyr0,pzr0,per0)


      dt1=tcol1-r(4,i1)
      if(dt1.le.0d0) then
         print *,'dt1=0? ',dt1
      endif
      xr0=r(1,i1)+ dt1*p(1,i1)/p(4,i1)
      yr0=r(2,i1)+ dt1*p(2,i1)/p(4,i1)
      zr0=r(3,i1)+ dt1*p(3,i1)/p(4,i1)
      tr0=tcol1

      call jamrobo(0.0d0,0.0d0,-bex,-bey,-bez,ga,xr0,yr0,zr0,tr0)
      call jamrobo(0d0,-phi,0d0,0d0,0d0,1.0d0,xr0,yr0,zr0,tr0)
      call jamrobo(-the,0d0,0d0,0d0,0d0,1.0d0,xr0,yr0,zr0,tr0)

      dt2=tcol2-r(4,i2)
      if(dt2.le.0d0) then
         print *,'dt2=0? ',dt2
      endif
      xr1=r(1,i2)+ dt2*p(1,i2)/p(4,i2)
      yr1=r(2,i2)+ dt2*p(2,i2)/p(4,i2)
      zr1=r(3,i2)+ dt2*p(3,i2)/p(4,i2)
      tr1=tcol2

      call jamrobo(0.0d0,0.0d0,-bex,-bey,-bez,ga,xr1,yr1,zr1,tr1)
      call jamrobo(0d0,-phi,0d0,0d0,0d0,1.0d0,xr1,yr1,zr1,tr1)
      call jamrobo(-the,0d0,0d0,0d0,0d0,1.0d0,xr1,yr1,zr1,tr1)

      if(abs(tr1-tr0).ge.1d-3.or.abs(zr1-zr0).ge.1d-3) then
        print *,'jamhit t1 != t2?',tr0,tr0
        print *,'zr1=',zr0,zr1
        print *,'p1=',pxr0,pyr0,pzr0
        print *,'b=',sqrt((xr0-xr1)**2+(yr0-yr1)**2+(zr1-zr0)**2)
      endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      end

