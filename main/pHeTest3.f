      include 'jam1.inc'
      include 'jam2.inc'
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      character frame*8,proj*8,targ*8,cwin*15
      real kinElist(50), ratioNNpiNDelta(50)
      character strMomList(50)*9
      character strKinElist(50)*6

c----Kinetic energies
      data kinElist/0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4,
     & 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 
     & 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0,
     & 16.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0,
     & 150.0, 200.0, 300.0, 500.0, 700.0, 1000.0/
      data strMomList/'8.827e-01','9.541e-01','1.090e+00','1.219e+00',
     & '1.343e+00','1.463e+00','1.581e+00','1.696e+00','1.921e+00',
     & '2.142e+00','2.358e+00','2.572e+00','2.784e+00','2.995e+00',
     & '3.203e+00','3.411e+00','3.618e+00','3.825e+00','4.030e+00',
     & '4.235e+00','4.440e+00','4.644e+00','4.848e+00','5.356e+00',
     & '5.863e+00','6.369e+00','6.874e+00','7.379e+00','7.882e+00',
     & '8.889e+00','9.894e+00','1.090e+01','1.290e+01','1.491e+01',
     & '1.691e+01','1.891e+01','2.092e+01','2.592e+01','3.092e+01',
     & '4.093e+01','5.093e+01','6.093e+01','8.093e+01','1.009e+02',
     & '1.509e+02','2.009e+02','3.009e+02','5.009e+02','7.009e+02',
     & '1.001e+03'/
      data strKinElist/'0.35','0.4','0.5','0.6','0.7','0.8','0.9','1.0',
     & '1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0',
     & '3.2','3.4','3.6','3.8','4.0','4.5','5.0','5.5','6.0','6.5',
     & '7.0','8.0','9.0','10.0','12.0','14.0','16.0','18.0','20.0',
     & '25.0','30.0','40.0','50.0','60.0','80.0','100.0','150.0',
     & '200.0','300.0','500.0','700.0','1000.0'/
      data ratioNNpiNDelta/0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.35, 0.3,
     & 0.25, 0.2, 0.15, 0.1, 0.05, 36*0.0/
c      data strECM/'2.177','2.311'/

      call jamnewparam

c...Initialize analysis part for overall run.
      call anal0

      nevent=10000
      mevent=nevent
      bmin=0.0D0
      bmax=-3.0D0
      dt=100.0D0
      nstep=1
      frame='lab'
      proj='p'
      targ='4He'
  
      mstc(8)=0 !This job mode 3 prints out all intermediate steps: was 0 in 2013
      mstc(21)=0   ! This is needed for multi run

      do iE=1,50  ! Loop for one iE

        write(6,*)'ie=',ie
c =0: no s-wave pion =1:s-wave pion
c =2:s-wave pion by reducing one-delta production
c     mstc(66)=2   ! 

c...the ratio of s-wave production cross section to delta cross section
c     parc(80)=ratioNNpiNDelta(iE) ! Add direct NNpi at 10% level
     
c     fname(2)='ppTest3Data/ppSummary'//strKinElist(iE)  !'GeV.dat' Summary of run
c     fname(3)='ppTest3Data/ppEvtList'//strKinElist(iE)  ! Detailed inf for all events
      cwin=strMomList(iE)//'gevc'
      call jaminit(mevent,bmin,bmax,dt,nstep,
     $                             frame,proj,targ,cwin)

c     print *,hint1(72),hint1(76)

c...Initialze event analysis
       call anal1

c...Loop over simulation
c==================================================
      do iev=1,nevent

c....Generate one event
        call jamevt(iev)

c       if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...Event analysis
         call anal2
c       call jamlist(4)
c....End simulation
      end do
c==================================================

c...Finsih event
      call jamfin

c...Output event analysis
       call anal3

      enddo  ! End for one iE

      call anal4

      end

c***********************************************************************

      subroutine jamnewparam

c...Initialize parameters for JAM

      implicit double precision(a-h, o-z)
      include 'jam2.inc'

c...transition range from JAM to HIJING in srt (GeV)
      parc(71)=8.0d0
      parc(72)=30.0d0

c...allow weak decay after simulation
c     mstc(42)=0
c     call jamdecsw(1)
c     mdcy(jamcomp(111),1)=0  ! no decay pi0
c     mdcy(jamcomp(3122),1)=1


c...N(1440)0  12112
      kc=jamcomp(12112)
      pmas(kc,1)=1.44  ! mass
      pmas(kc,2)=0.3   ! decay width
      pmas(kc,3)=pmas(kc,2)

c...N(1440)+  12212
      kc=jamcomp(12212)
      pmas(kc,1)=1.44  ! mass
      pmas(kc,2)=0.3   ! decay width
      pmas(kc,3)=pmas(kc,2)

      end
c***********************************************************************

      subroutine anal0

c...Initialize analysis for overall run.

      include 'jam1.inc'
      include 'jam2.inc'
      parameter(ncdet=19,npdet=3)
      character simfile(ncdet)*15
      character cfile(ncdet)*15
      dimension ncount(ncdet)
      dimension ntyp(npdet), mchan(npdet,ncdet), ktyp(npdet)

      save ncount,noreac

      data simfile/
     &             'pi1-.dat', 'pi10.dat','pi1+.dat'
     &            ,'pi-0.dat',   'pi-+.dat'
     &            ,'pi0+.dat',  'pi--.dat'
     &            ,'pi00.dat', 'pi++.dat'
     &            ,'pi---.dat', 'pi000.dat','pi+++.dat'
     &            ,'pi--0.dat', 'pi--+.dat','pi-0+.dat'
     &            ,'pi-++.dat', 'pi-00.dat','pi00+.dat'
     &            ,'pi0++.dat'/

      data cfile/
     &             'pi- + X',    'pi0 + X'
     &            ,'pi+ + X', 'pi- pi0 + X'
     &            ,'pi- pi+ + X',  'pi0 pi+ + X'
     &            ,'2pi- + X', '2pi0 + X'
     &            ,'2pi+ + X'
     &            ,'3pi- + X','3pi0 + X','3pi+ + X'
     &            ,'2pi-pi0 + X','2pi-pi+ + X','pi-pi0pi+ + X'
     &            ,'pi-2pi+ + X','pi-2pi0 + X','2pi0pi+ + X'
     &            ,'pi02pi+ + X'/

c........... pi- pi0 pi+
      data mchan/
     &        1, 0, 0,  !　pi-
     &        0, 1, 0,  !　pi0
     &        0, 0, 1,  !　pi+
     &        1, 1, 0,  !　pi- pi0
     &        1, 0, 1,  !　pi- pi+
     &        0, 1, 1,  !　pi0 pi+
     &        2, 0, 0,  !　2pi- 
     &        0, 2, 0,  !　2pi0
     &        0, 0, 2,  !　2pi+
     &        3, 0, 0,  !　3pi-
     &        0, 3, 0,  !　3pi0
     &        0, 0, 3,  !　3pi+
     &        2, 1, 0,  !　2pi- pi0
     &        2, 0, 1,  !　2pi- pi+
     &        1, 1, 1,  !　pi- pi0 pi+
     &        1, 0, 2,  !　pi-  2pi+
     &        1, 2, 0,  !　pi-  2pi0
     &        0, 2, 1,  !　2pi0 pi+
     &        0, 1, 2/  !　pi0 2pi+

c................. pi- pi0 pi+  K-   Kb0 K0 K+ eta
c     data ktyp / -211,111,211,-321,-311,311,321,221/
      data ktyp / -211,111,211/

      ifile=10
      do i=1,ncdet
        ifile=ifile+1
        open(ifile,file=simfile(i))
        write(ifile,'(a)')'# '//cfile(i)
        write(ifile,
     &   '(''# tlab(GeV) srt(GeV)  sigma(mb)  error mult.   count'')')
      enddo

      return

c***********************************************************************

      entry anal1

      do i=1,ncdet
        ncount(i)=0
      enddo

      noreac=0

      return

c***********************************************************************

      entry anal2

c...No inelasitc event
      if(nv.eq.5) then
         noreac = noreac+1
         return
      endif

      do i=1,npdet
        ntyp(i)=0
      enddo

c...Loop over all particles.
      notdef=0
      do i=1,nv

        kf=k(2,i)
        itag=0
        do ipdet=1,npdet
          if (kf.eq.ktyp(ipdet)) then
            itag=1
            ntyp(ipdet)=ntyp(ipdet)+1
          endif
        enddo
c       if(itag.eq.0) notdef=1

3000  end do  ! end loop over all particles

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
      if(notdef.eq.0.and.(ichanel.ge.1.and.ichanel.le.ncdet))
     $   ncount(ichanel)=ncount(ichanel)+1


c      write(6,*)'ichanel=',ichanel
c      io=mstc(38)
c      mstc(38)=6
c      call jamlist(1)
c      print *,ichanel,ntyp(1),ntyp(2),ntyp(3)
c      mstc(38)=io


      return

c***********************************************************************

      entry anal3

c...Output results.

c...Event weight
      wei=1.D0/dble(mstc(2))
      ratio = dble(mstc(2)-noreac)*wei
      print *,'noreac=',noreac,' r= ',ratio
      fac=parc(4)**2*paru(1)*10*wei

c...Inclusive data
      ifn=10
      do i=1,ncdet
        ifn=ifn+1
        e=0.0
        if(ncount(i).ge.1) e=1.0/sqrt(dble(ncount(i)))
        write(ifn,800)pard(14),pard(16),fac*ncount(i),
     &       fac*ncount(i)*e,
     &      ncount(i)*wei,ncount(i)
      enddo

800   format(2(f8.3,1x),f12.7,1x,f12.7,1x,f9.4,i5)
810   format(2(f8.3,1x),10(f10.7,1x))

      return

c***********************************************************************

      entry anal4

      ifn=10
      do i=1,ncdet
        ifn=ifn+1
        close(ifn)
      enddo

      end

