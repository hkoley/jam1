c...A main program for checking p+p reactions: energy dependence of
c...exclusive cross sections.

      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15
      character plabc*11,chap*16
      common/pppp/plab
c...Function:lab. momentum.
      plabsr(a,b,c)=sqrt((a**2-b**2-c**2)**2/(4.d0*c**2)-b**2)

c     call jamnewparam

c....K0 decay allowed
c       mdcy(jamcomp(311),1)=1   ! k0
c       mdcy(jamcomp(-311),1)=1 ! ak0
c       mstc(42)=1
c       parc(52)=2.1 ! mass cut off for S=0 baryon resonance

c     parc(64)=1.09d0 ! min. string mass
c     parc(71)=4.0d0
c     parc(72)=10.0d0
c     mstc(66)=1

      mstc(80)=1
      mstc(77)=1

      bmin=0.0D0          ! minimum impact parameter
      bmax=-1.0D0         ! maximum impact parameter
      mstc(1) =254777     ! random seed.
      mstc(1) =0          ! random seed.
c     mstc(1) =254371     ! random seed.

      dt=10.0D0           ! collision time(fm/c)
      mevent=30000        ! total simulation event
      nstep=1
      frame='nn'          ! comp. frame
      proj='p'            ! projectile
      targ='p'            ! target

c...for analysis of resoance production.
      mstc(8)=0        ! job mode.
      mstc(55)= 1      ! 1:frozen resonance.
      mstc(76)= 1      ! string has a lifetime.
      mstc(41)=0       ! forbid unstable particle decay.
      mstc(13)=2       ! all warning are printed.
      mstc(17)=0       ! only inelastic collisions or not.
      parc(7)=100.0d0  ! output interval.
c     mstc(71)=1

      em1=0.938D0
      em2=0.938D0

c     pmin=2.3D0
c     pmax=4.0D0
c     npmax=10

      pmin=0.90D0
      pmax=14D0
      pmax=2000D0
c     pmax=15D0
      npmax=50


c     pmin=20.0
c     pmax=21.0
c     npmax=2

      delp=(pmax-pmin)/npmax
c     delp=1.d0/npmax*log10(pmax/pmin)

c...Initialize anal for overall run.
      call anal0

c...Loop over incident momentum.
      do ip=1,npmax
        plab=pmin*(pmax/pmin)**((ip-1)/dble(npmax-1))
c       plab=pmin+(ip-1)*delp
  
c       srt=pmin*(pmax/pmin)**((ip-1)/DBLE(npmax-1))

c       srt=pmin+(ip-1)*delp
c       srt=pmin*10**(ip*delp)

c       write(66,*)'srt=',srt

c       plab=plabsr(srt,em1,em2)


c....Initialize JAM.
        write(plabc,'(f11.3)')plab
        cwin=plabc//'gevc'  ! incident energy
        call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)
        nevent=mstc(2)

c...Initialize analysis.
      call anal1

c...Simulation start.
      do iev=1,nevent

c...Simulate one event.
        call jamevt(iev)
        if(mod(iev,10000).eq.0) write(6,*)'event=',iev

c...Data analysis.
        call anal2

      end do

c...Final output.
c     call jamfin

c...Print analysis results.
      call anal3

      end do  ! end loop over momentum

      call anal4

      end


c***********************************************************************

      subroutine anal0

c...Initialize analysis for overall run.

      include 'jam1.inc'
      include 'jam2.inc'
      parameter(ncdet=16,npdet=14) ! stable particles
      parameter(npdet3=16)
      parameter(ncdet2=3)          ! unstable particles
      parameter(mstu6=500)         ! mstu(6)
      character simfile(ncdet)*15,simfile2(ncdet2)*15,
     &  simfile3(npdet3)*15
      character cfile(ncdet)*15,cfile2(ncdet2)*15,cfile3(npdet3)*15
      dimension ncount(ncdet),ncount2(ncdet2),ncount3(npdet3)
      dimension ntyp(npdet), mchan(npdet,ncdet), ktyp(npdet3)
      character chaf1*16,chaf2*16
      character chap*16,ckf*3

c...For resonance production.
      parameter(nchnl=14)
      character simfile4(nchnl)*15,cfile4(nchnl)*15
      dimension numf4(15),numf41(15),ncount4(nchnl),ncount7(nchnl)
      dimension ncount5(mstu6)

      parameter(nres=21)
      dimension ncount6(nres),ktypr(2,nres)
      character simfile6(nres)*15,cfile6(nres)*15
      dimension kfsave(100),kfsave4(100)

      common/pppp/plab

      save ncount,ncount2,ncount3,ncount4,ncount5,ncount6,ncount7
      save sigel,siginel

      data simfile6/
     &     'pp2pd+.dat',    'pp2nd++.dat',
     &     'pp2d++d0.dat',   'pp2d+d+.dat',
     &     'pp2pp1440.dat',  'pp2pp1525.dat',
     &    'pp2pp1688.dat',  'pp2pp1700.dat','pp2pp1720.dat',
     &    'pp2nd1920.dat',  'pp2n1525d++.dat',
     &    'pp2n1688d++.dat',
     &    'pp2pd+1600.dat','pp2pd+1620.dat','pp2pd+1700.dat',
     &    'pp2pd+1900.dat','pp2pd+1905.dat','pp2pd+1910.dat',
     &    'pp2pd+1920.dat','pp2pd+1930.dat','pp2pd+1950.dat'/
      data cfile6/
     &             'pp->pd+',    'pp->nd++',
     &             'pp->d++d0',   'pp->d+d+',
     &             'pp->pp1440',  'pp->pp1525',
     &             'pp->pp1688',  'pp->pp1700','pp-pp1720',
     &             'pp->nd1920++', 'pp->n1525 d++',
     &             'pp->n1688 d++',
     &             'pp->pd+1600','pp->pd+1620','pp->pd+1700',
     &             'pp->pd+1900','pp->pd+1905','pp->pd+1910',
     &             'pp->pd+1920','pp->pd+1930','pp->pd+1950'/

      data simfile/
     &             'pi1a.dat',    'pi1b.dat'
     &            ,'pi2a.dat',   'pi2b.dat', 'pi2c.dat'
     &            ,'pi3a.dat',  'pi3b.dat'
     &            ,'pi4.dat', 'pi5.dat'
     &            ,'lambda1.dat', 'sigma1.dat'
     &            ,'sigma2.dat',  'lambda2.dat'
     &            ,'lambda3.dat',  'eta.dat'
     $            ,'lambda4.dat'/

      data cfile/
     &             'p n pi+',    'p p pion0'
     &            ,'p n pi^+ pi^0', 'p p pi+ pi-','p p pi0 pi0'
     &            ,'pn 2pi+ pi-',  'pp pi+ pi0 pi-' 
     &            ,'pp 2pi+ 2pi-', 'pn 3pi+ 2pi-' 
     &            ,'LpK+',   'S+nK+' 
     &            ,'S0pK+',  'Ln+K+' 
     &            ,'Lp0K+',  'pp eta' 
     $            ,'Lambda p pi+ K0'/

      data simfile2/'pp-omega.dat','pp-rho0.dat','pp-eta.dat'/
      data cfile2/'pp-omega','pp-rho0','pp-eta'/

      data simfile3/'pp-pX.dat','pp-nX.dat','pp-LambdaX.dat',
     $              'pp-S-X.dat','pp-S0X.dat','pp-S+X.dat',
     $              'pp-pi-X.dat','pp-pi0X.dat','pp-pi+X.dat',
     $              'pp-K-X.dat','pp-Kb0X.dat','pp-K0X.dat',
     $              'pp-K+X.dat','pp-etaX.dat',
     $              'pp-omega+X.dat','pp-rho0X.dat'/

      data cfile3/'p + X','pp-n X','Lambda + X',
     $              'S^- +  X','S^0 + X.dat','S^+ + X',
     $              'pi^- X','pi^0 + X','pi^+ + X',
     $              'K^- X','Kb0 + X','K0 + X',
     $              'K^+ X','eta  + X',
     $              'omega + X','rho0  + X'/

      data simfile4/'pp-ND.dat','pp-NNs.dat','pp-DD.dat','pp-NDs.dat',
     $              'pp-NsD.dat','pp-DDs.dat','pp-NsNs.dat',
     $              'pp-NsDs.dat','pp-DsDs.dat',
     $              'pp-NR.dat','pp-RR.dat',
     $              'pp-NStr.dat','pp-RStr.dat','pp-StrStr.dat'/
      data cfile4/'N+D(1232)','NN*','DD','ND*',
     $            'N*D','DD*','N*N*','N*D*','D*D*',
     $            'NR','RR',
     $            'N+String','R+String','String+String'/

      data numf4/0,1,3,2,5,7,4,6,8,9,0,0,0,0,0/
      data numf41/0,10,11,10,2*11,10,3*11,12,3*13,14/
c1 n 2 d 3 n* 4 d* 5 str

c2  1 2  nd    1  10
c3  2 2  dd    3  11
c4  3 1  nn*   2  10
c5  3 2  dn*   5  11
c6  3 3  n*n*  7  11
c7  1 4  nd*   4  10
c8  2 4  dd*   6  11
c9  4 3  n*d*  8  11
c10 4 4  d*d*  9  11
c11 1 5  ns       12
c12 2 5           13
c13 3 5           13
c14 4 5           13
c15 5 5  ss       14

c........... n  p  L  S- S0 S+ p- p0 p+ K- Kb K0 K+ eta
      data mchan/
     &       1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   ! np pi+
     &       0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,   ! pp pi0
     &       1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,   ! np pi0 pi+
     &       0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,   ! pp pi- pi+
     &       0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0,   ! pp pi0 pi0
     &       1, 1, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0,   ! np pi- 2pi+
     &       0, 2, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0,   ! pp pi- pi0 pi+
     &       0, 2, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0,   ! pp 2pi- 2pi+
     &       1, 1, 0, 0, 0, 0, 2, 0, 3, 0, 0, 0, 0, 0,   ! np 2pi- 3pi+
     &       0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,   ! pL  K+
     &       1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,   ! nS+ K+
     &       0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0,   ! pS0 K+
     &       1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0,   ! nL  pi+ K+
     &       0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,   ! pL pi0 K+
     &       0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,   ! pp eta
     &       0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 /  ! Lp pi+ K0
      data ktyp
     &      /2112,2212,3122,3112,3212,3222,
     $              -211,111,211,-321,-311,311,321,221,223,113/

      data ktypr/
     &    2212,2214,  ! p d+
     &    2112,2224,  ! n d++
     &    2224,2114,  ! d++ d0
     &    2214,2214,  ! d+d+
     &    2212,12212, ! pp1440
     &    2212,2124,  ! pp1525
     &    2212,12216, ! pp1688
     &    2212,22124, ! pp(1700)
     &    2212,32124, ! pp(1720)
     &    2112,22224, ! nD(1920)++
     &    1214,2224,  ! n(1525)D++
     &    12116,2224, ! n(1688)D++
c
     &    2212,32214, ! pD(1600)+
     &    2212, 2122, ! pD(1620)+
     &    2212,12214, ! pD(1700)+
     &    2212,12122, ! pD(1900)+
     &    2212, 2126, ! pD(1905)+
     &    2212,22122, ! pD(1910)+
     &    2212,22214, ! pD(1920)+
     &    2212,12126, ! pD(1930)+
     &    2212, 2218/ ! pD(1950)+

      ifile=10
      do i=1,ncdet
        ifile=ifile+1
        open(ifile,file=simfile(i))
        write(ifile,'(a)')'# '//cfile(i)
        write(ifile,'(''# plab(gev/c) srt(GeV)  sigma(mb)  mult.'')')
      enddo

      do i=1,ncdet2
        ifile=ifile+1
        open(ifile,file=simfile2(i))
        write(ifile,'(a)')'# '//cfile2(i)
        write(ifile,'(''# plab(gev/c) srt(GeV)  sigma(mb)  mult.'')')
      enddo

      do i=1,npdet3
        ifile=ifile+1
        open(ifile,file=simfile3(i))
        write(ifile,'(a)')'# '//cfile3(i)
        write(ifile,'(''# plab(gev/c) srt(GeV)  sigma(mb)  mult.'')')
      enddo

      do i=1,nchnl
        ifile=ifile+1
        open(ifile,file=simfile4(i))
        write(ifile,'(a)')'# '//cfile4(i)
        write(ifile,'(''# plab(gev/c) srt(GeV)  sigma(mb)  mult.'')')
      enddo

      do i=1,nchnl
        ifile=ifile+1
        open(ifile,file='ch3_'//simfile4(i))
        write(ifile,'(a)')'# '//cfile4(i)
        write(ifile,'(''# plab(gev/c) srt(GeV)  sigma(mb)  mult.'')')
      enddo

      do i=1,nres
        ifile=ifile+1
        open(ifile,file=simfile6(i))
        write(ifile,'(a)')'# '//cfile6(i)
        write(ifile,'(''# plab(gev/c) srt(GeV)  sigma(mb)  mult.'')')
      enddo


c      ifile=100
c      do kc=1,mstu6
c        kf=kchg(kc,4)
c        if(kf.eq.0) goto 300
c        call pjname(kf,chap)
c        if(chap.eq.'') goto 300
c        ifile=ifile+1
c        leng=index(chap,' ')-1
c      print *,ifile,chap,ckf
cc       write(ckf,*) ifile
cc       leng=index(ckf,' ')-1
cc       open(ifile,file=ckf//'.dat')
cc       open(ifile,file=chap(1:leng)//'.dat')
c        write(ifile,*)'# ', chap
c        write(ifile,'(''# plab(gev/c) srt(GeV)  sigma(mb)  mult.'')')
c300   end do

      open(1,file='sig.dat',status='unknown')

      return

c***********************************************************************

      entry anal1

      do i=1,ncdet
        ncount(i)=0
      enddo
      do i=1,ncdet2
        ncount2(i)=0
      enddo
      do i=1,npdet3
        ncount3(i)=0
      enddo
      do i=1,nchnl
        ncount4(i)=0
        ncount7(i)=0
      enddo

      do i=1,nres
        ncount6(i)=0
      enddo

      do i=1,mstu6
        ncount5(i)=0
      enddo

      sigel=0d0
      siginel=0d0

      return

c***********************************************************************

      entry anal2

c...Count data.
c----------------------------------------------------------------
c...Count resonance production.
c     if(nv.ne.2) print *,'nv=',nv,(k(2,i),i=1,nv)
      em1=0.0d0
      em2=0.0d0
      kf1=0
      kf2=0
      kc1=0
      kc2=0

      if(nv.eq.2) then

        kf1=k(2,1)
        kf2=k(2,2)
        k1=k(1,1)
        k2=k(1,2)
        em1=p(5,1)
        em2=p(5,2)
        kc1=jamcomp(kf1)
        kc2=jamcomp(kf2)
c       if(kc1.ne.0) ncount5(kc1) = ncount5(kc1)+1
c       if(kc2.ne.0) ncount5(kc2) = ncount5(kc2)+1
        call pjname(kf1,chaf1)
        call pjname(kf2,chaf2)


c.....pp ->  resonance production
        do ip=1,nres
          if (ktypr(1,ip).eq.kf1.and.ktypr(2,ip).eq.kf2)
     &            ncount6(ip)=ncount6(ip)+1
          if (ktypr(1,ip).eq.kf2.and.ktypr(2,ip).eq.kf1)
     &            ncount6(ip)=ncount6(ip)+1
        enddo

        if(k1.eq.1) then
          id1=1
        else if(k1.eq.2) then
          id=kchg(kc1,5)
          if(id.eq.id_delt) id1=2
          if(id.eq.id_nucls) id1=3
          if(id.eq.id_delts) id1=4
        else
          id1=5
        endif

        if(k2.eq.1) then
          id2=1
        else if(k2.eq.2) then
          id=kchg(kc2,5)
          if(id.eq.id_delt)  id2=2
          if(id.eq.id_nucls) id2=3
          if(id.eq.id_delts) id2=4
        else
          id2=5
        endif

c     if((id1.eq.3.and.id2.eq.1).or.(id1.eq.1.and.id2.eq.3)) then
c        print *,'NN* collision ',kf1,kf2,chaf1,chaf2
c        pause 
c     endif

        idmin=min(id1,id2)
        idmax=max(id1,id2)
        icpair=(idmax*(idmax-1))/2+idmin
        i1=numf4(icpair)
        i2=numf41(icpair)
        if(i1.ne.0) ncount4(i1)=ncount4(i1)+1
        if(i2.ne.0) ncount4(i2)=ncount4(i2)+1
c       write(66,*)kf1,kf2,' ',chaf1,' ',chaf2

      endif
c----------------------------------------------------------------

      mstc(41)=1       ! unstable particle decay
      mdcy(jamcomp(113),1)=0  ! rho0
      mdcy(jamcomp(223),1)=0  ! omega
c     mdcy(jamcomp(221),1)=0  ! eta
      mdcy(jamcomp(2224),1)=0  ! Dleta++
      mdcy(jamcomp(2214),1)=0  ! Dleta+
      mdcy(jamcomp(2114),1)=0  ! Dleta0
      mdcy(jamcomp(1114),1)=0  ! Dleta-

c....Decay string
      call jamfdec

c...save information
      mv=nv
      do i=1,nv
        kfsave(i) = k(2,i)
        kfsave4(i) = k(4,i)
      end do

c...pp omega (rho0) final
      if(nv.eq.3) then
         if(k(2,1).eq.2212.and.k(2,2).eq.2212.and.k(2,3).eq.223) ! omega
     $      ncount2(1)=ncount2(1)+1
         if(k(2,1).eq.2212.and.k(2,2).eq.2212.and.k(2,3).eq.113) then ! rho
            ncount2(2)=ncount2(2)+1
            write(170,630)pard(16),kf1,kf2,em1,em2,chaf1,chaf2
c           io=mstc(38)
c           mstc(38)=170
c           call jamlist(1)
c           mstc(38)=io
           if(i1.ne.0) ncount7(i1)=ncount7(i1)+1
           if(i2.ne.0) ncount7(i2)=ncount7(i2)+1
         endif
         if(k(2,1).eq.2212.and.k(2,2).eq.2212.and.k(2,3).eq.221) ! eta
     $      ncount2(3)=ncount2(3)+1
      endif

c...pp omega + X  or rho0 + X
      do i=1,nv
        if(k(2,i).eq.223) ncount3(15)=ncount3(15)+1
        if(k(2,i).eq.113) ncount3(16)=ncount3(16)+1
      enddo

      mdcy(jamcomp(113),1)=1  ! rho0
      mdcy(jamcomp(223),1)=1  ! omega
c     mdcy(jamcomp(221),1)=1  ! eta
      mdcy(jamcomp(2224),1)=1  ! Dleta++
      mdcy(jamcomp(2214),1)=1  ! Dleta+
      mdcy(jamcomp(2114),1)=1  ! Dleta0
      mdcy(jamcomp(1114),1)=1  ! Dleta-

      mstc(41)=1       ! unstable particle decay
      call jamfdec
      mstc(41)=0       ! forbid unstable particle decay

      if(mstc(17).eq.1.and.nv.eq.2) then
        print *,'elastic??? nv=',nv
        do i=1,nv
          write(6,*)k(1,i),k(2,i),p(5,i)
        end do
        stop
      endif

      do i=1,npdet
        ntyp(i)=0
      enddo

c...Loop over all particles.
      notdef=0
      do i=1,nv
        kf=k(2,i)
c                     if(kf.eq.221) then
c                        print *,'eta?'
c                     endif
c            if(kf.eq.abs(311)) then
c              print *,'K0 ?',kf
c              stop
c            endif
                      if(kf.eq.223) then
                         print *,'omega?'
                      endif
                      if(kf.eq.113) then
                         print *,'rho0?'
                      endif
        itag=0
        do ipdet=1,npdet
          if (kf.eq.ktyp(ipdet)) then
            itag=1
            ntyp(ipdet)=ntyp(ipdet)+1
            ncount3(ipdet)=ncount3(ipdet)+1
          endif
        enddo
        if(itag.eq.0) notdef=1
3000  end do

      ichanel=-99
c...Loop over reaction channel.
      do icdet=1,ncdet
c....Loop over particle involed this reaction channel.
        do ipdet=1,npdet
          if (ntyp(ipdet).ne.mchan(ipdet,icdet)) goto 10
        enddo
        ichanel=icdet
        goto 4000
10      continue
      enddo

4000  continue
c     if(notdef.eq.0.and.(ichanel.ge.1.and.ichanel.le.ncdet))
      if(ichanel.ge.1.and.ichanel.le.ncdet)
     $   ncount(ichanel)=ncount(ichanel)+1

c...Count elastic collision.
      if(nv.eq.2.and.k(2,1).eq.2212.and.k(2,2).eq.2212) then
        sigel=sigel+1d0
      else
        siginel=siginel+1d0
      endif

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c     if(ichanel.eq.2.or.ichanel.eq.3)then
c      if(nv.ne.3) then
c      write(90,*)' '
c      write(90,*)ichanel,mste(1)
c    $   ,kf1,kf2,em1,em2,chaf(kc1,1),' ',chaf(kc2,1)
c      do i=1,nv
c        kc=jamcomp(k(2,i))
c        write(90,*)k(1,i),k(2,i),p(5,i),' ',chaf(kc,1)
c      end do
c     endif
c     endif
      if(ichanel.eq.6) then
        write(166,630)pard(16),kf1,kf2,em1,em2,chaf1,chaf2
        io=mstc(38)
        mstc(38)=166
        call jamlist(1)
        mstc(38)=io
      endif

c....pp->pn pi+ pi0
      if(ichanel.eq.3) then
        write(167,630)pard(16),kf1,kf2,em1,em2,chaf1,chaf2
        do i=1,mv
         call pjname(kfsave(i),chap)
         write(167,*) kfsave4(i),kfsave(i),chap
        end do
        io=mstc(38)
        mstc(38)=167
        call jamlist(1)
        mstc(38)=io
c       if(i1.ne.0) ncount7(i1)=ncount7(i1)+1
c       if(i2.ne.0) ncount7(i2)=ncount7(i2)+1
      endif

c....pp->pp pi+ pi-
      if(ichanel.eq.4) then
        write(168,630)pard(16),kf1,kf2,em1,em2,chaf1,chaf2
        do i=1,mv
         call pjname(kfsave(i),chap)
         write(168,*)kfsave4(i), kfsave(i),chap
        end do
        io=mstc(38)
        mstc(38)=168
        call jamlist(1)
        mstc(38)=io
      endif

      if(ichanel.eq.8) then
        write(169,630)pard(16),kf1,kf2,em1,em2,chaf1,chaf2
        io=mstc(38)
        mstc(38)=169
        call jamlist(1)
        mstc(38)=io
      endif
630   format(f8.3,1x,2(i6,1x),2(f12.7,1x),2a16)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      return

c***********************************************************************

      entry anal3

c...Output results.

c...Event weight
      wei=1.D0/dble(mstc(2))
      fac=parc(4)**2*paru(1)*10*wei
      srt=pard(16)
      tlab=pard(14)

c...Inclusive data
      ifn=10
      do i=1,ncdet
        ifn=ifn+1
        e=0.0
        if(ncount(i).ge.1) e=1.0/sqrt(dble(ncount(i)))
        write(ifn,820)plab,srt,fac*ncount(i),
     &         fac*ncount(i)*e,ncount(i)*wei,ncount(i)
      enddo

c...Final unstable particle
      do i=1,ncdet2
        ifn=ifn+1
        e=0.0
        if(ncount2(i).ge.1) e=1.0/sqrt(dble(ncount2(i)))
        write(ifn,820)plab,srt,fac*ncount2(i),
     &         fac*ncount2(i)*e,ncount2(i)*wei,ncount2(i)
      enddo

c...Exclusive data
      do i=1,npdet3
        ifn=ifn+1
        e=0.0
        if(ncount3(i).ge.1) e=1.0/sqrt(dble(ncount3(i)))
        write(ifn,820)plab,srt,fac*ncount3(i),
     &         fac*ncount3(i)*e,ncount3(i)*wei,ncount3(i)
      enddo

c...Resonance productions.
      do i=1,nchnl
        ifn=ifn+1
        write(ifn,800)plab,srt,fac*ncount4(i),ncount4(i)*wei
      enddo
      do i=1,nchnl
        ifn=ifn+1
        write(ifn,800)plab,srt,fac*ncount7(i),ncount7(i)*wei
      enddo

c...Specific resonance productions.
      do i=1,nres
        ifn=ifn+1
        e=0.0
        if(ncount6(i).ge.1) e=1.0/sqrt(dble(ncount6(i)))
        write(ifn,830)plab,srt,fac*ncount6(i),
     &         fac*ncount6(i)*e,ncount6(i)*wei,ncount6(i)
      enddo

c      ifile=100
c      do kc=1,mstu6
c        if(kf.eq.0) goto 310
c        ifile=ifile+1
c        kf=kchg(kc,4)
c        call pjname(kf,chap)
c        write(ifn,800)plab,srt,fac*ncount5(i),ncount5(i)*wei
c310   end do

      write(1,810)plab,srt,pare(4),pare(5),sigel*fac
     $ ,pare(4)-pare(5),siginel*fac


830   format(2(f8.3,1x),2(f12.7,1x),f9.4,i8)
820   format(2(f8.3,1x),2(f12.7,1x),f9.4,i5)
800   format(2(f8.3,1x),f12.7,1x,f9.4)
810   format(2(f8.3,1x),10(f10.7,1x))

      return

c***********************************************************************

      entry anal4

      ifn=10
      do i=1,ncdet+ncdet2+npdet3+nchnl+nchnl+nres
        ifn=ifn+1
        close(ifn)
      enddo

      end

c***********************************************************************

      subroutine jamnewparam

c...Initialize parameters for JAM

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension rdec(4),r1520(5),r1720(8)

c...transition range from JAM to HIJING in srt
c     parc(71)=8.0d0
c     parc(72)=30.0d0

      parc(71)=10.0d0
      parc(72)=30.0d0

      parc(64)=1.5d0 ! min. string mass

      call nstar
      call deltastar


c...allow weak decay after simulation
c     mstc(42)=0
c     call jamdecsw(1)
c     mdcy(jamcomp(111),1)=0  ! no decay pi0
c     mdcy(jamcomp(3122),1)=1

c     kc=jamcomp(221)
c     print *,pmas(kc,1),pmas(kc,2),pmas(kc,3)
c     pmas(kc,2)=0.3   ! decay width
c     pmas(kc,3)=0.3   ! decay width

      n1440=0
      n1520=0
      n1720=0

      if(n1440.eq.1) then
c....original N1440 branch
      rdec(1)=0.05  ! n+sigma
      rdec(2)=0.65  ! n+pi
      rdec(3)=0.25  ! d+pi
      rdec(4)=0.05  ! n+rho

      rdec(1)=0.10  ! n+sigma
      rdec(2)=0.55  ! n+pi
      rdec(3)=0.20  ! d+pi
      rdec(4)=0.15  ! n+rho

c     rdec(1)=0.01  ! n+sigma
c     rdec(2)=0.01  ! n+pi
c     rdec(3)=0.97  ! d+pi
c     rdec(4)=0.01  ! n+rho

c...N(1440)0  12112
      kc=jamcomp(12112)
      pmas(kc,1)=1.44  ! mass
      pmas(kc,2)=0.3   ! decay width
      pmas(kc,3)=0.3   ! decay width
      idc=mdcy(kc,2)
      brat(idc)=rdec(1)          ! n0  sigma
      brat(idc+1)=rdec(2)*1./3.  ! n0  pi0
      brat(idc+2)=rdec(2)*2./3.  ! p+  pi-
      brat(idc+3)=rdec(3)*1./2.  ! Delta-   pi+
      brat(idc+4)=rdec(3)*1./3.  ! Delta0   pi0
      brat(idc+5)=rdec(3)*1./6.  ! Delta+   pi-
      brat(idc+6)=rdec(4)*1./3.  ! n0   rho0
      brat(idc+7)=rdec(4)*2./3.  ! p    rho-


c...N(1440)+  12212
      kc=jamcomp(12212)
      pmas(kc,1)=1.44  ! mass
      pmas(kc,2)=0.3   ! decay width
      pmas(kc,3)=0.3   ! decay width
      idc=mdcy(kc,2)
      brat(idc)=rdec(1)          ! p0  sigma
      brat(idc+1)=rdec(2)*1./3.  ! p+  pi0
      brat(idc+2)=rdec(2)*2./3.  ! n0  pi+
      brat(idc+3)=rdec(3)*1./2.  ! Delta++   pi-
      brat(idc+4)=rdec(3)*1./3.  ! Delta+   pi0
      brat(idc+5)=rdec(3)*1./6.  ! Delta0   pi+
      brat(idc+6)=rdec(4)*1./3.  ! p+   rho0
      brat(idc+7)=rdec(4)*2./3.  ! n0   rho+

      endif

      if(n1520.eq.1) then
c...N(1520)0  1214
c     r1520(1)=0.045  ! n+sigma
c     r1520(2)=0.55   ! n+pi
c     r1520(3)=0.225  ! d+pi
c     r1520(4)=0.175  ! n+rho
c     r1520(5)=0.001  ! n+eta

      r1520(1)=0.020  ! n+sigma
      r1520(2)=0.55   ! n+pi
      r1520(3)=0.179  ! d+pi
      r1520(4)=0.250  ! n+rho
      r1520(5)=0.001  ! n+eta

c     r1520(1)=0.001  ! n+sigma
c     r1520(2)=0.001  ! n+pi
c     r1520(3)=0.001  ! d+pi
c     r1520(4)=0.996  ! n+rho
c     r1520(5)=0.001  ! n+eta

      kc=jamcomp(1214)
c     pmas(kc,1)=1.52  ! mass
c     pmas(kc,2)=0.115 ! decay width
c     pmas(kc,3)=0.115 ! decay width
      idc=mdcy(kc,2)
      brat(idc)=r1520(1)          ! n0  sigma
      brat(idc+1)=r1520(2)*1./3.  ! n0  pi0
      brat(idc+2)=r1520(2)*2./3.  ! p+  pi-
      brat(idc+3)=r1520(3)*1./2.  ! Delta-   pi+
      brat(idc+4)=r1520(3)*1./3.  ! Delta0   pi0
      brat(idc+5)=r1520(3)*1./6.  ! Delta+   pi-
      brat(idc+6)=r1520(4)*1./3.  ! n0   rho0
      brat(idc+7)=r1520(4)*2./3.  ! p    rho-
      brat(idc+8)=r1520(5)        ! n    eta



c...N(1520)0  2124
      kc=jamcomp(2124)
c     pmas(kc,1)=1.52  ! mass
c     pmas(kc,2)=0.115 ! decay width
c     pmas(kc,3)=0.115 ! decay width
      idc=mdcy(kc,2)
      brat(idc)=r1520(1)          ! p+  sigma
      brat(idc+1)=r1520(2)*1./3.  ! p+  pi0
      brat(idc+2)=r1520(2)*2./3.  ! n0  pi+
      brat(idc+3)=r1520(3)*1./2.  ! Delta++   pi-
      brat(idc+4)=r1520(3)*1./3.  ! Delta+   pi0
      brat(idc+5)=r1520(3)*1./6.  ! Delta0   pi+
      brat(idc+6)=r1520(4)*1./3.  ! p+   rho0
      brat(idc+7)=r1520(4)*2./3.  ! n0   rho+
      brat(idc+8)=r1520(5)        ! p+   eta

      endif


      if(n1720.eq.1) then
c...N(1720)0  31214
c     r1720(1)=0.1    ! n+sigma
c     r1720(2)=0.15   ! n+pi
c     r1720(3)=0.100  ! d+pi
c     r1720(4)=0.200  ! n+rho
c     r1720(5)=0.04   ! n+eta
c     r1720(6)=0.05   ! L+K
c     r1720(7)=0.06   ! S+K
c     r1720(8)=0.3    ! n+omega

      r1720(1)=0.00   ! n+sigma
      r1720(2)=0.10   ! n+pi
      r1720(3)=0.03   ! d+pi
      r1720(4)=0.85   ! n+rho
      r1720(5)=0.01   ! n+eta
      r1720(6)=0.01   ! L+K
      r1720(7)=0.00   ! S+K
      r1720(8)=0.00   ! n+omega

      kc=jamcomp(31214)
c     pmas(kc,1)=1.52  ! mass
c     pmas(kc,2)=0.115 ! decay width
c     pmas(kc,3)=0.115 ! decay width
      idc=mdcy(kc,2)
      brat(idc)=r1720(1)          ! n0  sigma
      mdme(idc,1)=0
      brat(idc+1)=r1720(2)*1./3.  ! n0  pi0
      brat(idc+2)=r1720(2)*2./3.  ! p+  pi-
      brat(idc+3)=r1720(3)*1./2.  ! Delta-   pi+
      brat(idc+4)=r1720(3)*1./3.  ! Delta0   pi0
      brat(idc+5)=r1720(3)*1./6.  ! Delta+   pi-
      brat(idc+6)=r1720(4)*1./3.  ! n0   rho0
      brat(idc+7)=r1720(4)*2./3.  ! p    rho-
      brat(idc+8)=r1720(5)        ! n    eta
      brat(idc+9)=r1720(6)        ! L   K
      brat(idc+10)=r1720(7)*1./3.  ! S0   K0
      brat(idc+11)=r1720(7)*2./3.  ! S- K+
      mdme(idc+11,1)=0
      brat(idc+12)=r1720(7)       ! p omega
      mdme(idc+12,1)=0


c...N(1720)+  32124
      kc=jamcomp(32124)
c     pmas(kc,1)=1.52  ! mass
c     pmas(kc,2)=0.115 ! decay width
c     pmas(kc,3)=0.115 ! decay width
      idc=mdcy(kc,2)
      brat(idc)=r1720(1)          ! p+  sigma
      mdme(idc,1)=0
      brat(idc+1)=r1720(2)*1./3.  ! p+  pi0
      brat(idc+2)=r1720(2)*2./3.  ! n0  pi+
      brat(idc+3)=r1720(3)*1./2.  ! Delta++   pi-
      brat(idc+4)=r1720(3)*1./3.  ! Delta+   pi0
      brat(idc+5)=r1720(3)*1./6.  ! Delta0   pi+
      brat(idc+6)=r1720(4)*1./3.  ! p+   rho0
      brat(idc+7)=r1720(4)*2./3.  ! n0   rho+
      brat(idc+8)=r1720(5)        ! p+   eta
      brat(idc+9)=r1720(6)        ! L   K
      brat(idc+10)=r1720(7)*1./3.  ! S0   K+
      brat(idc+11)=r1720(7)*2./3.  ! S+ K0
      mdme(idc+11,1)=0
      brat(idc+12)=r1720(7)       ! p omega
      mdme(idc+12,1)=0

      endif

      end


c***********************************************************************

      subroutine nstar

c...Initialize parameters for JAM

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
c...Nucleon*
      parameter(mxbar1=11,maxb1=9)
      common/resn1/emdl1(mxbar1),gamdl1(mxbar1),lres1(maxb1,mxbar1)
      common/resn2/brnch1(maxb1,mxbar1),emmes1(maxb1),embar1(maxb1)
      common/resn2o/brnch1o(maxb1,mxbar1)
      common/resn3/charn(mxbar1),ipdg1(2,mxbar1)
      character charn*8

c...Mass.
      data emdl1/
     & 1.44,1.52,1.535,1.65,1.675,1.68,1.70,1.71,1.72,1.99,2.19/
c...Width.
c     data gamdl1/
c    & 0.35,0.12,0.150,0.15,0.150,0.13,0.10,0.10,0.15,0.35,0.45/
      data gamdl1/
     & 0.30,0.115,0.15,0.15,0.15,0.13,0.15,0.10,0.25,0.35,0.5/

      data charn/
     $ 'N(1440)', 'N(1520)', 'N(1535)', 'N(1650)', 'N(1675)','N(1680)',
     $ 'N(1700)', 'N(1710)', 'N(1720)', 'N(1990)','N(2190)'/
      data (ipdg1(1,i),i=1,mxbar1)/
     $ 12112, 1214, 22112, 32112, 2116, 12116, 21214, 42112,
     $ 31214, 11218, 1218/
      data (ipdg1(2,i),i=1,mxbar1)/
     $  12212, 2124, 22212, 32212, 2216, 12216, 22124, 42212,
     $  32124, 12128, 2128/

c... Branching ratios for N* decay (JAM 1.1x)
c...    Ns   Npi     Dpi   Nrho   Neta   LK     SK   Nomg  Neta'
      data brnch1o /
     1 0.050 ,0.65 ,0.250 ,0.050 ,0.000 ,0.000 ,0.00 ,0.00 ,0.00 !N(1440)
     2,0.049 ,0.55 ,0.225 ,0.175 ,0.001 ,0.000 ,0.00 ,0.00 ,0.00 !N(1520)
     3,0.050 ,0.45 ,0.050 ,0.050 ,0.400 ,0.000 ,0.00 ,0.00 ,0.00 !N(1535)
     4,0.050 ,0.70 ,0.050 ,0.130 ,0.010 ,0.060 ,0.00 ,0.00 ,0.00 !N(1650)
     5,0.009 ,0.34 ,0.550 ,0.090 ,0.010 ,0.001 ,0.00 ,0.00 ,0.00 !N(1675)
     6,0.150 ,0.65 ,0.100 ,0.100 ,0.000 ,0.000 ,0.00 ,0.00 ,0.00 !N(1680)
     7,0.399 ,0.10 ,0.399 ,0.100 ,0.000 ,0.002 ,0.00 ,0.00 ,0.00 !N(1700)
     8,0.015 ,0.15 ,0.150 ,0.125 ,0.300 ,0.200 ,0.06 ,0.00 ,0.00 !N(1710)
     9,0.100 ,0.15 ,0.100 ,0.200 ,0.040 ,0.050 ,0.06 ,0.30 ,0.00 !N(1720)
     $,0.030 ,0.15 ,0.100 ,0.200 ,0.040 ,0.050 ,0.10 ,0.30 ,0.03 !N(1990)
     1,0.157 ,0.15 ,0.110 ,0.300 ,0.030 ,0.003 ,0.02 ,0.20 ,0.02 !N(2190)
     &                                                          /
c... Branching ratios for N* decay
c...    Ns   Npi     Dpi   Nrho   Neta   LK     SK   Nomg  Neta'
      data brnch1 /
     1 0.100 ,0.55 ,0.300 ,0.050 ,0.000 ,0.000 ,0.00 ,0.00 ,0.00 !N(1440)
     2,0.001 ,0.55 ,0.249 ,0.200 ,0.000 ,0.000 ,0.00 ,0.00 ,0.00 !N(1520)
     3,0.090 ,0.45 ,0.040 ,0.000 ,0.420 ,0.000 ,0.00 ,0.00 ,0.00 !N(1535)
     4,0.050 ,0.78 ,0.050 ,0.050 ,0.070 ,0.000 ,0.00 ,0.00 ,0.00 !N(1650)
     5,0.070 ,0.40 ,0.525 ,0.004 ,0.000 ,0.001 ,0.00 ,0.00 ,0.00 !N(1675)
     6,0.110 ,0.69 ,0.150 ,0.050 ,0.000 ,0.000 ,0.00 ,0.00 ,0.00 !N(1680)
     7,0.008 ,0.12 ,0.800 ,0.070 ,0.000 ,0.002 ,0.00 ,0.00 ,0.00 !N(1700)
c    8,0.150 ,0.10 ,0.150 ,0.050 ,0.200 ,0.200 ,0.02 ,0.13 ,0.00 !N(1710)
     8,0.150 ,0.10 ,0.150 ,0.050 ,0.330 ,0.200 ,0.02 ,0.00 ,0.00 !N(1710)
     9,0.000 ,0.11 ,0.020 ,0.700 ,0.040 ,0.070 ,0.06 ,0.00 ,0.00 !N(1720)
     $,0.030 ,0.15 ,0.300 ,0.000 ,0.040 ,0.050 ,0.10 ,0.30 ,0.03 !N(1990)
     1,0.157 ,0.15 ,0.110 ,0.300 ,0.030 ,0.003 ,0.02 ,0.20 ,0.02 !N(2190)
     &                                                          /

c....Loop over all N*
      do ip=1,mxbar1-1
        kf1=ipdg1(1,ip)    ! N*0
        kf2=ipdg1(2,ip)    ! N*+
        kc1=jamcomp(kf1)
        kc2=jamcomp(kf2)

        pmas(kc1,1)=emdl1(ip)  ! mass
        pmas(kc1,2)=gamdl1(ip) ! decay width
        pmas(kc1,3)=pmas(kc1,2)

        pmas(kc2,1)=emdl1(ip)  ! mass
        pmas(kc2,2)=gamdl1(ip) ! decay width
        pmas(kc2,3)=pmas(kc2,2)
        idc1=mdcy(kc1,2)-1
        idc2=mdcy(kc2,2)-1

        brsum=0.0
        do idec=1,maxb1  ! loop over all decay branch
          wd=brnch1(idec,ip)
          wdo=brnch1o(idec,ip)
          brsum=brsum+wd
          if(wd.gt.0.0d0 .or. wdo.gt.0.0d0) then
            idc1=idc1+1
            idc2=idc2+1
            if(idec.eq.1) then ! N sigma-meson
              brat(idc1)=wd     ! n0 + sigma
              brat(idc2)=wd     ! p+ + sigma
              if(wd.eq.0.0d0) then
                 mdme(idc1,1)=0
                 mdme(idc2,1)=0
              endif
            else if(idec.eq.2) then ! N pi
              brat(idc1)=wd*1./3.  ! n0  pi0
              brat(idc2)=wd*1./3.  ! p+  pi0
              idc1=idc1+1
              idc2=idc2+1
              brat(idc1)=wd*2./3.  ! p+  pi-
              brat(idc2)=wd*2./3.  ! n0  pi+
            else if(idec.eq.3) then

              brat(idc1)=wd*1./2.  ! Delta-   pi+
              brat(idc2)=wd*1./2.  ! Delta++   pi-

              idc1=idc1+1
              idc2=idc2+1
              brat(idc1)=wd*1./3.  ! Delta0   pi0
              brat(idc2)=wd*1./3.  ! Delta+   pi0

              idc1=idc1+1
              idc2=idc2+1
              brat(idc1)=wd*1./6.  ! Delta+   pi-
              brat(idc2)=wd*1./6.  ! Delta0   pi+

            else if(idec.eq.4) then

              brat(idc1)=wd*1./3.  ! n0   rho0
              brat(idc2)=wd*1./3.  ! p+   rho0
c                mdme(idc1,1)=0
c                mdme(idc2,1)=0
              if(wd.eq.0.0d0) then
                 mdme(idc1,1)=0
                 mdme(idc2,1)=0
              endif

              idc1=idc1+1
              idc2=idc2+1
              brat(idc1)=wd*2./3.  ! p    rho-
              brat(idc2)=wd*2./3.  ! n0   rho0
c                mdme(idc1,1)=0
c                mdme(idc2,1)=0
              if(wd.eq.0.0d0) then
                 mdme(idc1,1)=0
                 mdme(idc2,1)=0
              endif

            else if(idec.eq.5) then
              brat(idc1)=wd        ! n    eta
              brat(idc2)=wd        ! p    eta
              if(wd.eq.0.0d0) then
                 mdme(idc1,1)=0
                 mdme(idc2,1)=0
              endif
            else if(idec.eq.6) then
              brat(idc1)=wd        ! L   K
              brat(idc2)=wd        ! L   K
              if(wd.eq.0.0d0) then
                 mdme(idc1,1)=0
                 mdme(idc2,1)=0
              endif
            else if(idec.eq.7) then

              brat(idc1)=wd*1./3.  ! S0   K0
              brat(idc2)=wd*1./3.  ! S0   K+

              idc1=idc1+1
              idc2=idc2+1
              brat(idc1)=wd*2./3.  ! S- K+
              brat(idc2)=wd*2./3.  ! S+ K0

            else if(idec.eq.8) then
              brat(idc1)=wd       ! n omega
              brat(idc2)=wd       ! p omega
              if(wd.eq.0.0d0) then
                 mdme(idc1,1)=0
                 mdme(idc2,1)=0
              endif
            else if(idec.eq.9) then
              brat(idc1)=wd        ! n    eta'
              brat(idc2)=wd        ! p    eta'
            endif
          endif

       enddo

        if(abs(brsum).gt.0.0005.and.abs(brsum-1.).gt.0.0005) then
           print *,'Sum of branching ratio not 1.00 N*',ip,brsum
           stop
        endif

      enddo

      end

c***********************************************************************

      subroutine deltastar

c...Set new parameters for Delta*

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
c...Delta*
      parameter(mxbar2=9,maxb2=5)
      common/resd1/emdl2(mxbar2),gamdl2(mxbar2),lres2(maxb2,mxbar2)
      common/resd2/brnch2(maxb2,mxbar2),emmes2(maxb2),embar2(maxb2)
      common/resd2o/brnch2o(maxb2,mxbar2)
      common/resd3/chard(mxbar2),ipdg2(4,mxbar2)
      character chard*8

c... Angular Momentum J*2+1
c     data jspd   /4 ,2 ,4 ,2 ,6 ,2 ,4 ,6 ,8/
      data emdl2  /1.60,1.62,1.70,1.90,1.905,1.91,1.92,1.93, 1.95/
c     data gamdl2 /0.35,0.15,0.30,0.20,0.350,0.25,0.20,0.35, 0.30/
      data gamdl2 /0.32,0.14,0.30,0.20,0.330,0.28,0.26,0.36, 0.285/
      data chard/'D(1600)','D(1620)','D(1700)','D(1900)','D(1905)',
     $           'D(1910)','D(1920)','D(1930)','D(1950)'/

      data (ipdg2(1,i),i=1,mxbar2)/
     $   31114,1112,11114,11112,1116,21112,21114,11116,1118/
      data (ipdg2(2,i),i=1,mxbar2)/
     $   32114,1212,12114,11212,1216,21212,22114,11216,2118/
      data (ipdg2(3,i),i=1,mxbar2)/
     $   32214,2122,12214,12122,2126,22122,22214,12126,2218/
      data (ipdg2(4,i),i=1,mxbar2)/
     $   32224,2222,12224,12222,2226,22222,22224,12226,2228/

c... Branching ratios for delta* decay for JAM1.1xx
c...               N*pi  Npi    Dpi       Nrho    SK
      data brnch2o/
     1             0.15  ,0.175 ,0.550   ,0.125  ,0.00
     2            ,0.05  ,0.250 ,0.500   ,0.200  ,0.00
     3            ,0.00  ,0.200 ,0.400   ,0.3984 ,0.0016
     4            ,0.40  ,0.100 ,0.050   ,0.450  ,0.00
     5            ,0.00  ,0.100 ,0.14998 ,0.750  ,0.00002
     6            ,0.60  ,0.225 ,0.025   ,0.150  ,0.00
     7            ,0.37  ,0.200 ,0.400   ,0.000  ,0.03
     8            ,0.00  ,0.150 ,0.000   ,0.850  ,0.00
     9            ,0.19  ,0.400 ,0.300   ,0.100  ,0.01
     &                                                /

c... Branching ratios for delta* decay
c...               N*pi  Npi    Dpi       Nrho    SK
      data brnch2/
     1             0.15  ,0.175 ,0.550   ,0.125  ,0.00    !D(1600)
     2            ,0.05  ,0.300 ,0.550   ,0.100  ,0.00    !D(1620)
     3            ,0.00  ,0.200 ,0.400   ,0.3984 ,0.0016  !D(1700)
     4            ,0.40  ,0.100 ,0.050   ,0.450  ,0.00    !D(1900)
     5            ,0.00  ,0.100 ,0.14998 ,0.750  ,0.00002 !D(1905)
     6            ,0.174 ,0.225 ,0.600   ,0.001  ,0.00    !D(1910)
     7            ,0.37  ,0.200 ,0.400   ,0.000  ,0.03    !D(1920)
     8            ,0.00  ,0.999 ,0.000   ,0.001  ,0.00    !D(1930)
     9            ,0.19  ,0.450 ,0.300   ,0.050  ,0.01    !D(1950)
     &                                                /

c....Loop over all D*
      do ip=1,mxbar2
        kf1=ipdg2(1,ip)    ! D*-
        kf2=ipdg2(2,ip)    ! D*0
        kf3=ipdg2(3,ip)    ! D*+
        kf4=ipdg2(4,ip)    ! D*++
        kc1=jamcomp(kf1)
        kc2=jamcomp(kf2)
        kc3=jamcomp(kf3)
        kc4=jamcomp(kf4)

c       pmas(kc1,1)=emdl1(ip)  ! mass
c       pmas(kc1,2)=gamdl1(ip) ! decay width
c       pmas(kc1,3)=pmas(kc1,2)


        idc1=mdcy(kc1,2)-1
        idc2=mdcy(kc2,2)-1
        idc3=mdcy(kc3,2)-1
        idc4=mdcy(kc4,2)-1

        brsum=0.0
        do idec=1,maxb2  ! loop over all decay branch
          wd=brnch2(idec,ip)
          wdo=brnch2o(idec,ip)
          brsum=brsum+wd
          if(wd.gt.0.0d0 .or. wdo.gt.0.0d0) then
            idc1=idc1+1
            idc2=idc2+1
            if(idec.eq.1) then ! N* pi

              brat(idc1)=wd        ! D- -> N* pi-

              brat(idc2)=wd*2./3.  ! D0 -> n* pi0
              idc2=idc2+1
              brat(idc2)=wd*1./3.  ! D0 -> p* pi-

              brat(idc3)=wd*2./3.  ! D+ -> p* pi0
              idc3=idc3+1
              brat(idc3)=wd*1./3.  ! D+ -> n* pi+

              brat(idc4)=wd        ! D++ -> p* pi+

            else if(idec.eq.2) then ! N pi

              brat(idc1)=wd  ! D-  -> n0  pi-
              brat(idc4)=wd  ! D++ -> p+  pi+

              brat(idc2)=wd*2./3.  ! D0 -> n0  pi0
              idc2=idc2+1
              brat(idc2)=wd*1./3.  ! D0 -> p+  pi-

              brat(idc3)=wd*2./3.  ! D+ -> p+  pi0
              idc3=idc3+1
              brat(idc3)=wd*1./3.  ! D+ -> n0  pi+

            else if(idec.eq.3) then ! Delta pi

              brat(idc1)=wd*3./5.  ! D- -> D- + pi0
              idc1=idc1+1
              brat(idc1)=wd*2./5.  ! D- -> D0 + pi-
              brat(idc4)=wd*3./5.  ! D++ -> D++ + pi0
              idc4=idc4+1
              brat(idc4)=wd*2./5.  ! D++ -> D+ + pi+

              brat(idc2)=wd*1./15.  ! D0 -> D0 + pi0
              idc2=idc2+1
              brat(idc2)=wd*8./15.  ! D0 -> D+ + pi-
              idc2=idc2+1
              brat(idc2)=wd*2./5.   ! D0 -> D- + pi+

              brat(idc3)=wd*1./15.  ! D+ -> D+ + pi0
              idc3=idc3+1
              brat(idc3)=wd*8./15.  ! D+ -> D0 + pi+
              idc3=idc3+1
              brat(idc3)=wd*2./5.   ! D+ -> D++ + pi-


            else if(idec.eq.4) then ! N + rho
 
              brat(idc1)=wd  ! D- -> n0   rho-
              brat(idc4)=wd  ! D++ -> p+   rho+

c                mdme(idc1,1)=0
c                mdme(idc4,1)=0
c                mdme(idc2,1)=0
c                mdme(idc2+1,1)=0
c                mdme(idc3,1)=0
c                mdme(idc3+1,1)=0

              brat(idc2)=wd*2./3.  ! D0 -> n0   rho0
              idc2=idc2+1
              brat(idc2)=wd*1./3.  ! D0 -> p+   rho-

              brat(idc3)=wd*2./3.  ! D0 -> p+   rho0
              idc3=idc3+1
              brat(idc3)=wd*1./3.  ! D0 -> n0   rho+


            else if(idec.eq.5) then ! S + K

              brat(idc1)=wd  ! D- -> sig- + k0
              brat(idc4)=wd  ! D++ -> sig+ + k+

              brat(idc2)=wd*2./3.  ! D0 -> sig0   K0
              idc2=idc2+1
              brat(idc2)=wd*1./3.  ! D0 -> sig-   K+

              brat(idc3)=wd*2./3.  ! D+ -> sig0   K+
              idc3=idc3+1
              brat(idc3)=wd*1./3.  ! D+ -> sig+   K0

            endif
          endif

       enddo

        if(abs(brsum).gt.0.0005.and.abs(brsum-1.).gt.0.0005) then
           print *,'Sum of branching ratio not 1.00 N*',ip,brsum
           stop
        endif

      enddo

      end
