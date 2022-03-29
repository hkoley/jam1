c   p(10GeV/c)+p
c
      include 'jam1.inc'
      include 'jam2.inc'
      character frame*8,proj*8,targ*8,cwin*15
      real kinElist(12), ratioNNpiNDelta(12)
      character strMomList(12)*8
      character strKinElist(12)*8

c----Kinetic energies
      data kinElist/0.65, 0.97, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 
     &8.0, 9.0, 10.0/
      data strMomList/'1.2815','1.6617','1.6960','2.7845','3.8249',
     &'4.8483','5.8637','6.8746','7.8827','8.8889','9.8939','10.8980'/
      data strKinElist/'0.65','0.97','1.0','2.0','3.0','4.0','5.0',
     &'6.0','7.0','8.0','9.0','10.0'/
      data ratioNNpiNDelta/0.5, 0.4, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0,
     & 0.0, 0.0, 0.0, 0.0/
c      data strECM/'2.177','2.311'/

c....Initialize JAM.
c      fname(1)='ppTest1.cfg'  ! input file name. Do not use because 
                               ! of minor bugs in JAM which need to be 
                               ! fixed by Yasushi Nara

c     do iE=1,12  ! Loop for one iE
      do iE=1,1  ! Loop for one iE
  
      mstc(8)=3    ! Job mode
      mstc(17)=1   ! only inelastic collision
c     mstc(18)=1   ! do not call jamedit
      mstc(21)=0   ! This is needed for multi run

c =0: no s-wave pion =1:s-wave pion
c =2:s-wave pion by reducing one-delta production
      mstc(66)=2

c...the ratio of s-wave production cross section to delta cross section
c     parc(80)=ratioNNpiNDelta(iE) ! Add direct NNpi at 10% level
      parc(80)=1.0

      fname(2)='ppTest2Data/ppEvt'//strKinElist(iE)  !'GeV.dat' Summary of events
c      fname(3)='ppChk'//strKinElist(iE)   ! Detailed inf for 10evt
c     nevent=10000
      nevent=100
      mevent=nevent
      bmin=0.0D0
      bmax=-2.0D0
      dt=100.0D0
      nstep=1
c      win=10.D0
c      win=2.7845D0
      frame='lab'
      proj='p'
      targ='p'
      cwin=strMomlist(iE)//'gevc'
      cwin='2.3gevc'
      call jaminit(mevent,bmin,bmax,dt,nstep,
     $                             frame,proj,targ,cwin)

c...Initialze event analysis
c      call anal1

c...Loop over simulation
c==================================================
      do iev=1,nevent

c....Generate one event
        call jamevt(iev)

        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...Event analysis
c        call anal2
c        call jamlist(1)
c        call jamlist(2)
c        call jamlist(3)
        call jamlist(4)
c....End simulation
      end do
c==================================================

c...Finsih event
      call jamfin

c...Output event analysis
c      call anal3

      enddo  ! End for one iE

      end
