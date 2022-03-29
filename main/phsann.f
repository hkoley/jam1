************************************************************************
	program phsann
************************************************************************
*     Purpose:
*	Analyze jam phase space data, and calculate various observables.
*     Ver: 1.30
*
c for isw=1, following files are output.
*----------------------------------------------------------------------*
*	FileExt		X	Y
*----------------------------------------------------------------------*
*	'mTy-N' ,'dmT','dN/dmT^2dY'	#
*	'mTe-N' ,'dmT','dN/dmT^2dE'	#
*	'pTy-N' ,'pT' ,'dN/dpT^2dY'	#
*	'pTe-N' ,'pT' ,'dN/dpT^2dE'	#
*	'pTy-v2','pT' ,'v2(pT,ycut)'	#
*	'pTe-v2','pT' ,'v2(pT,etacut)'	#
*	'y-N'   ,'Y'  ,'dN/dY'		#
*	'y-pT'  ,'Y'  ,'<p_T>(Y)'	#
*	'y-px'  ,'Y'  ,'<p_x>(Y)'	#
*	'y-v2'  ,'Y'  ,'v2(Y)'		#
*	'y-v1'  ,'Y'  ,'v1(Y)'		#
*	'eta-N' ,'eta','dN/deta'	#
*	'eta-v2','eta','v2(eta)'	#
*	'eta-v1','eta','v1(eta)'	#
*	'v2c'   ,'pT' ,'v2c(pT)'	#
c
c for isw=2, following files are output.
*----------------------------------------------------------------------*
*	FileExt	  X	Y
*----------------------------------------------------------------------*
c	'y_pi'  ,'y' ,'Various': for chaged pi
c	'y_pro' ,'y' ,'Various': for proton
c	'pT_pi' ,'pT','Various': for chaged pi
c	'pT_pro','pT','Various': for proton
c  pT_pi, pT_pro, 
c	pT(dmT)	,'dN/dmt(Y)','dN/dpt(Y)','v2(pt,Y)','v1(pt,Y)'
c     &		,'dN/dmt(e)','dN/dpt(e)','v2(pt,e)','v1(pt,e)'
c  y_pi, y_pro, 
c	y(eta)	,'dN/dY','<P_x>','v2(Y)','v1(Y)'
c     &	  	,'dN/deta','v2(eta)','v1(eta)'
************************************************************************
	implicit double precision(a-h,o-z)
	parameter(mxy=100,mxk=6)
	parameter(mxptcl=100)
	parameter(mpT=100)
	common /par1/pi,srtnn,dy,dymid,ypcm,dm,bmin,bmax,isw
	common /par2/r0,r1,r2,r3,ipk
	common /par3/v2ev,nc,ncm
	common /dim0/ObsY(0:mxk)
	common /dim1/dmtY(mpT,mxptcl),dptY(mpT,mxptcl),elpY(mpT,mxptcl)
	common /dimY/dmtE(mpT,mxptcl),dptE(mpT,mxptcl),elpE(mpT,mxptcl)
	common /dimV/v1pY(mpT,mxptcl),v1pE(mpT,mxptcl)
	common /dim2/drp(0:mxk,-mxy:mxy,mxptcl)
	common /dim3/det(0:mxk,-mxy:mxy,mxptcl)
	common /v2corrEV/partEV,pairEV,v2cEV
c
	common /partinf/em,p1,p2,p3,wei,rap,eta,id,iQ
c...
	common /parM1/dy0,dymid0,dm0,bmin0,bmax0
	common /parM2/nevent0,nevent,nev,ifileEv
	common /parM3/isym,iyrat
	common/pydatr/mrpy(6),rrpy(100)
c...[ For IO
	common /IOch1/FnameEv(20)
	common /IOevp/bminEv(20),bmaxEv(20),WeiEv(20)
     &		   ,nevEv(20),ipkEv(20),nfileEv,iprec(mxptcl),npt
	character FnameEv*80
c...] For IO
c ---------------------------------------------------------------------*
c...Local Variables
	character ch1*1
c ---------------------------------------------------------------------*
	iseed=0
	isw=2
        call jamcpu(6,0,iseed,isec)
        mrpy(1)=iseed
	call annini(1)

c...[ New Incident Energy ---------------------------------------------*
 2000	continue
	call annini(2)
	call vcann(1)
	nevent0=0

c...[ Event File
	do 1000 ifileEv=1,nfileEv
	      call annini(3)
	      nev=0
	      write(*,*) FnameEv(ifileEv)(1:lengch(FnameEv(ifileEv)))
	 do iev=1,nevent
	      call vcann(2)
	      call headline(jev,ntot,nbar,nmes,bimp,ipk,icon)
	      if(icon.ne.0) goto 200
	      if(mod(iev,100).eq.0) write(*,*) 'N_event=',iev
	  if(bimp.ge.bmin.and.bimp.le.bmax) then
	      nev=nev+1
	      nc=0
	      ncm=0
	      v2ev=0.0d0
	   do i=1,ntot
	      call readptcl(id,em,p1,p2,p3,r1,r2,r3,r0,wei,ipk,iok)
	    if(icon.eq.0) then
	      call addbin(id,em,p1,p2,p3,wei,rap,eta,iQ)
c	      call vcann(3)
	    endif
	   enddo
	      call vcann(4)
	      v2ev1=0.0d0
	      if(ncm.gt.0) v2ev1=v2ev/ncm
	      if(isw.eq.1) 
     &	      write(42,811) bimp,v2ev1,v2cEV,ncm,nc,ntot,nbar,nmes
811	      format(3(1x,1pe10.3),5(1x,i5))
	  else
	   do i=1,ntot
	      read(15,'(a)') ch1
	   enddo
	  endif
	 enddo
 200	 continue
	    call annini(4)
 1000	continue
c...] Event File
	call annini(5)
c...[ Midrapidity -----------------------------------------------------*
c...[ Mt, Pt Dist.
c...dN/dY/Mt/dMt/(2*pi)
c...dN/dY/Pt/dPt/(2*pi)
	call transDMT(dmtY,1,100,mxptcl)
	call transDMT(dptY,1,100,mxptcl)
	call transDMT(elpY,1,100,mxptcl)
	call transDMT(v1pY,1,100,mxptcl)
	call transDMT(dmtE,1,100,mxptcl)
	call transDMT(dptE,1,100,mxptcl)
	call transDMT(elpE,1,100,mxptcl)
	call transDMT(v1pE,1,100,mxptcl)
	dymid1=dymid
	if(dymid.eq.0) dymid1=0.5d0
	fact=1.0d0/pi/2/dm/dymid1/2/nevent0
	do i=1,100
	  do jp=1,npt
	    ip=iprec(jp)
	    dmtY(i,ip)=dmtY(i,ip)*fact
	    dptY(i,ip)=dptY(i,ip)*fact
	    elpY(i,ip)=elpY(i,ip)*fact
	    v1pY(i,ip)=v1pY(i,ip)*fact
	    dmtE(i,ip)=dmtE(i,ip)*fact
	    dptE(i,ip)=dptE(i,ip)*fact
	    elpE(i,ip)=elpE(i,ip)*fact
	    v1pE(i,ip)=v1pE(i,ip)*fact
	    if(dptY(i,ip).ne.0) elpY(i,ip)=elpY(i,ip)/dptY(i,ip)
	    if(dptE(i,ip).ne.0) elpE(i,ip)=elpE(i,ip)/dptE(i,ip)
	    if(dptY(i,ip).ne.0) v1pY(i,ip)=v1pY(i,ip)/dptY(i,ip)
	    if(dptE(i,ip).ne.0) v1pE(i,ip)=v1pE(i,ip)/dptE(i,ip)
	  enddo
	enddo
c...] Mt, Pt Dist.
c...] Midrapidity -----------------------------------------------------*

c...[ Rapidity, PseudoRapdity Dist. -----------------------------------*
	call transDNDY(drp,mxk,mxy,mxptcl)
	call transDNDY(det,mxk,mxy,mxptcl)
c...[ Rapidity
	nxy=min(nint(ypcm*2/dy),mxy)
	do i=-nxy,nxy
	  do jp=1,npt
	    ip=iprec(jp)
	    do L=1,mxk
	      if(drp(0,i,ip).ne.0) drp(L,i,ip)=drp(L,i,ip)/drp(0,i,ip)
	      if(det(0,i,ip).ne.0) det(L,i,ip)=det(L,i,ip)/det(0,i,ip)
	    enddo
	    drp(0,i,ip)=drp(0,i,ip)/dy/nevent0
	    det(0,i,ip)=det(0,i,ip)/dy/nevent0
	  enddo
	enddo
c...] Rapidity
	if(isym.eq.1) then
	  do i=0,mxy
	    do jp=1,npt
	      ip=iprec(jp)
	      do L=0,6
		isgn=1
		if(L.eq.2.or.L.eq.4) isgn=-1
	        drp(L, i,ip)=(drp(L,i,ip)+isgn*drp(L,-i,ip))/2
	        drp(L,-i,ip)=isgn*drp(L,i,ip)
	        det(L, i,ip)=(det(L,i,ip)+isgn*det(L,-i,ip))/2
	        det(L,-i,ip)=isgn*det(L,i,ip)
	      enddo
	    enddo
	  enddo
	endif
c...] Rapidity, PseudoRapdity Dist. -----------------------------------*

c...[ Record
	if(isw.eq.1) then	! General Purpose
	  call recDMT(dmtY,1,'mTy-N' ,'dmT','dN/dmT^2dY'	)
	  call recDMT(dmtE,1,'mTe-N' ,'dmT','dN/dmT^2dE'	)
	  call recDMT(dptY,1,'pTy-N' ,'pT' ,'dN/dpT^2dY'	)
	  call recDMT(dptE,1,'pTe-N' ,'pT' ,'dN/dpT^2dE'	)
	  call recDMT(elpY,1,'pTy-v2','pT' ,'v2(pT,ycut)'	)
	  call recDMT(elpE,1,'pTe-v2','pT' ,'v2(pT,etacut)')
	  call recDNDY(drp,0,'y-N'   ,'Y'  ,'dN/dY'	)
	  call recDNDY(drp,1,'y-pT'  ,'Y'  ,'<p_T>(Y)'	)
	  call recDNDY(drp,2,'y-px'  ,'Y'  ,'<p_x>(Y)'	)
	  call recDNDY(drp,3,'y-v2'  ,'Y'  ,'v2(Y)'	)
	  call recDNDY(drp,4,'y-v1'  ,'Y'  ,'v1(Y)'	)
	  call recDNDY(det,0,'eta-N' ,'eta','dN/deta'	)
	  call recDNDY(det,3,'eta-v2','eta','v2(eta)'	)
	  call recDNDY(det,4,'eta-v1','eta','v1(eta)'	)
c...
	  call recini(31,nevent0,2,'v2c','pT' ,'v2c(pT)'	)
	  call vcann(5)
	  close(31)

c...[Proton and chaged pion observables
	elseif(isw.eq.2) then
	  call recini(31,nevent0,2,'pT_pro','pT' ,'Various')
	  write(31,802) 'pt or dmt'
     &	  	,'dN/dmt(Y)','dN/dpt(Y)','v2(pt,Y)','v1(pt,Y)'
     &	  	,'dN/dmt(e)','dN/dpt(e)','v2(pt,e)','v1(pt,e)'
	  ip=iprec(10)
	  nxpt=1
	  do i=1,mpT
	    if(dmtY(i,ip).ne.0.and.i.gt.nxpt) nxpt=i
	    if(dptY(i,ip).ne.0.and.i.gt.nxpt) nxpt=i
	    if(dmtE(i,ip).ne.0.and.i.gt.nxpt) nxpt=i
	    if(dptE(i,ip).ne.0.and.i.gt.nxpt) nxpt=i
	  enddo
	  do i=1,nxpt
	    demt=i*dm-dm/2
	    write(31,801) demt
     &	  	,dmtY(i,ip),dptY(i,ip),elpY(i,ip),v1pY(i,ip)
     &	  	,dmtE(i,ip),dptE(i,ip),elpE(i,ip),v1pE(i,ip)
	  enddo
	  close(31)
c
	  call recini(31,nevent0,2,'y_pro','Y' ,'Various')
	  write(31,802) 'Y or eta'
     &	  	,'dN/dY','<P_x>','v2(Y)','v1(Y)'
     &	  	,'dN/deta','v2(eta)','v1(eta)'
	  ip=iprec(10)
	  nxy=min(nint(ypcm*2/dy),mxy)
	  do i=-nxy,nxy
	    yrap=dy*i
	    if(iyrat.eq.1) yrap=dy/ypcm*i
	    write(31,801) yrap
     &	  	,drp(0,i,ip),drp(2,i,ip),drp(3,i,ip),drp(4,i,ip)
     &	  	,det(0,i,ip),det(3,i,ip),det(4,i,ip)
	  enddo
	  close(31)

c...pion   observables
	  call recini(31,nevent0,2,'pT_pi','pT' ,'Various')
	  write(31,802) 'pt or dmt'
     &	  	,'dN/dmt(Y)','dN/dpt(Y)','v2(pt,Y)','v1(pt,Y)'
     &	  	,'dN/dmt(e)','dN/dpt(e)','v2(pt,e)','v1(pt,e)'
	  ip=iprec(12)
	  nxpt=1
	  do i=1,mpT
	    if(dmtY(i,ip).ne.0.and.i.gt.nxpt) nxpt=i
	    if(dptY(i,ip).ne.0.and.i.gt.nxpt) nxpt=i
	    if(dmtE(i,ip).ne.0.and.i.gt.nxpt) nxpt=i
	    if(dptE(i,ip).ne.0.and.i.gt.nxpt) nxpt=i
	  enddo
	  do i=1,nxpt
	    demt=i*dm-dm/2
	    write(31,801) demt
     &	  	,dmtY(i,ip),dptY(i,ip),elpY(i,ip),v1pY(i,ip)
     &	  	,dmtE(i,ip),dptE(i,ip),elpE(i,ip),v1pE(i,ip)
	  enddo
	  close(31)
c
	  call recini(31,nevent0,2,'y_pi','Y' ,'Various')
	  write(31,802) 'Y or eta'
     &	  	,'dN/dY','<P_x>','v2(Y)','v1(Y)'
     &	  	,'dN/deta','v2(eta)','v1(eta)'
	  ip=iprec(12)
	  nxy=min(nint(ypcm*2/dy),mxy)
	  do i=-nxy,nxy
	    yrap=dy*i
	    if(iyrat.eq.1) yrap=dy/ypcm*i
	    write(31,801) yrap
     &	  	,drp(0,i,ip),drp(2,i,ip),drp(3,i,ip),drp(4,i,ip)
     &	  	,det(0,i,ip),det(3,i,ip),det(4,i,ip)
	  enddo
	  close(31)
	endif
c...]
c
801	format(21(1x,1pe10.3))
802	format("# ",21(a9,2x))
c...] Record
	goto 2000
c...] New Incident Energy ---------------------------------------------*
        call jamcpu(6,1,iseed,isec)
	end
************************************************************************
	subroutine annini(ic)
************************************************************************
	implicit double precision(a-h,o-z)
	parameter(mxy=100,mxk=6)
	parameter(mxptcl=100)
	parameter(mpT=100)
	common /par2/r0,r1,r2,r3,ipk
	common /par3/v2ev,nc,ncm
	common /dim0/ObsY(0:mxk)
	common /dim1/dmtY(mpT,mxptcl),dptY(mpT,mxptcl),elpY(mpT,mxptcl)
	common /dimY/dmtE(mpT,mxptcl),dptE(mpT,mxptcl),elpE(mpT,mxptcl)
	common /dimV/v1pY(mpT,mxptcl),v1pE(mpT,mxptcl)
	common /dim2/drp(0:mxk,-mxy:mxy,mxptcl)
	common /dim3/det(0:mxk,-mxy:mxy,mxptcl)
c...
	common /par1/pi,srtnn,dy,dymid,ypcm,dm,bmin,bmax,isw
	common /parM1/dy0,dymid0,dm0,bmin0,bmax0
	common /parM2/nevent0,nevent,nev,ifileEv
	common /parM3/isym,iyrat
c...[ For IO
	common /IOch1/FnameEv(20)
	common /IOevp/bminEv(20),bmaxEv(20),WeiEv(20)
     &		   ,nevEv(20),ipkEv(20),nfileEv,iprec(mxptcl),npt
	character FnameEv*80,DirName*15
	common /headname/head
	character head*26
	logical exex
	save /headname/
	save DirName
	data DirName/'phs'/
c...] For IO
c ---------------------------------------------------------------------*
c...Local Variables
	character fname*50,fnamei*65
c ---------------------------------------------------------------------*
	goto (1000,2000,3000,4000,5000) ic

c...[ Real Initial Call (ic=1) ----------------------------------------*
1000	continue
	write(41,'(a)') '#  Mid Rapidity Observables'
	write(41,'(a)') '#  0: dN/dY'
	write(41,'(a)') '#  1: <Pt>'
	write(41,'(a)') '#  2: v2'
	write(41,'(a)') '#  3: <Mt>'
	write(41,'(a)') '# 10: dN/deta'
	write(41,'(a)') '# 11: <Pt>'
	write(41,'(a)') '# 12: v2'
	write(41,'(a)') '# 13: <Mt>'
	pi=4*atan(1.0d0)

c...Common Info. Inputs
	write(*,*) 'dY, dM_t, dY(M_t), isym=?'
	read(*,*) dy0,dm0,dymid0,isym
	write(*,*) 'b_min, b_max=?'
	read(*,*) bmin0,bmax0
	write(*,*) 'Directory Name for Phase Space Data=?'
	read(*,*) DirName
c
	return
	
c...[ New Incident Energy (ic=2) --------------------------------------*
2000	continue
c...Incident Energy Inputs
	write(*,*) 'srt_nn =? (Negative: Inc. E/A)'
	read(*,*) srtnn
	if(srtnn.eq.0) then
	  write(*,*) 'Good Bye !'
	  stop
	endif
	write(*,*) 'Header of Output File Name=?'
	read(*,*) head
	if(head(1:1).eq.'0'.or.head(1:3).eq.'end') then
	  write(*,*) 'Good Bye !'
	  stop
	endif
c
	if(isw.eq.1) then
	  open(42,file=head(1:lengch(head))//'.nbev',status='unknown')
	endif
c
	dy=dy0
	dm=dm0
	dymid=dymid0
	bmin=bmin0
	bmax=bmax0
c
	do k=1,mxptcl
	  do i=1,100
	    dmtY(i,k)=0.0d0
	    dptY(i,k)=0.0d0
	    elpY(i,k)=0.0d0
	    v1pY(i,k)=0.0d0
	    dmtE(i,k)=0.0d0
	    dptE(i,k)=0.0d0
	    elpE(i,k)=0.0d0
	    v1pE(i,k)=0.0d0
	  enddo
	  do l=0,mxk
	    do i=-mxy,mxy
	      drp(l,i,k)=0.0d0
	      det(l,i,k)=0.0d0
	    enddo
	  enddo
	enddo
c...[
	srtnn0=srtnn
	emn=0.938d0
	   if(srtnn.lt.0) then
	einc=-srtnn
	srtnn=sqrt(2*emn*(einc+2*emn))
	   endif
        srtn1=srtnn/2
        pn1=sqrt(srtn1**2-emn**2)
        ypcm=log((srtn1+pn1)/(srtn1-pn1))/2
c
	if(dymid.lt.0) dymid=-dymid
	if(dy.lt.0) then
		iyrat=1
		dy=ypcm*abs(dy)
	endif
	if(iyrat.eq.1) then
		dymid=dymid*ypcm
	endif
c...]

	write(*,'(a,a)')     '# Output Header  =',head
	if(srtnn0.lt.0) then
	  write(*,'(a,f10.3)') '# Incident Energy(AGeV)=',-srtnn0
	else
	  write(*,'(a,f10.3)') '# sqrt(s_nn) (GeV) =',srtnn
	endif
c...[ Event File ------------------------------------------------------*
	do ifileEv=1,20
	  read(*,*,err=999,end=1100)
     &		fnamei,nevEv(ifileEv),ipkEv(ifileEv)
	  if(fnamei(1:1).eq.'0'
     &   .or.fnamei(1:3).eq.'end'
     &   .or.fnamei(1:4).eq.'next'
     &	) then
	  	goto 1100
	  endif
	  FnameEv(ifileEv)=DirName(1:lengch(DirName))//fnamei
	  nfileEv=ifileEv
	  nch=lengch(FnameEv(ifileEv))
c	write(*,*) nch
	  write(*,'(1x,i2,":",a)') ifileEv
     &		,DirName(1:lengch(DirName))//fnamei(1:lengch(fnamei))
c     &		,FnameEv(ifileEv)(l:nch)
	enddo
1100	continue
	return
c...] Event File

c...[ At the Beginning of Each Event File (ic=3) ----------------------*
3000	continue
	write(*,*)
	open(15,file=FnameEv(ifileEv),status='old',err=999)
	nevent=nevEv(ifileEv)
	ipk=ipkEv(ifileEv)
	write(*,'("# ",a,", ",i5," Events")')
     &	FnameEv(ifileEv)(1:lengch(FnameEv(ifileEv))),nevent
	if(ipk.eq.1.or.ipk.eq.2) read(15,'(a)') fname
	if(ipk.eq.2) read(15,'(a)') fname

	return

c...[ At the End       of Each Event File (ic=4) ----------------------*
4000	continue
	close(15)
	ifile=6
	write(ifile,'(a,a)')     '# Input Event File=',FnameEv(ifileEv)
	write(ifile,'(a,1x,i5)') '# Number of Events=',nev
	nevent0=nevent0+nev

	return

c...[ At the End of Each Incident Energy (ic=5) -----------------------*
5000	continue

	ifile=6
	write(ifile,'(a,1x,i5)') '# Total Events    =',nevent0
c...
	return

c...]]]]] ic=1,2,3 (goto 1000,2000,3000) ------------------------------*

999	continue
	write(*,*) 'Abnormal End'
	stop

	end
************************************************************************
	subroutine recDMT(DMT,n1,FileExt,Xname,Yname)
************************************************************************
	implicit double precision(a-h,o-z)
	parameter(mxptcl=100)
	parameter(mpT=100)
	common /par1/pi,srtnn,dy,dymid,ypcm,dm,bmin,bmax,isw
	common /parM2/nevent0,nevent,nev,ifileEv
	common /IOevp/bminEv(20),bmaxEv(20),WeiEv(20)
     &		   ,nevEv(20),ipkEv(20),nfileEv,iprec(mxptcl),npt
	dimension DMT(n1:mpT,mxptcl)
	character*(*) Xname,Yname,FileExt
c
	call recini(16,nevent0,1,FileExt,Xname,Yname)
	do i=1,mpT
	  demt=i*dm-dm/2
	  write(16,801) demt,(DMT(i,iprec(jp)),jp=1,npt)
	enddo
	close(16)
801	format(21(1x,1pe10.3))
	end
************************************************************************
	subroutine recDNDY(DNDY,IrecY,FileExt,Xname,Yname)
************************************************************************
	implicit double precision(a-h,o-z)
	parameter(mxptcl=100)
	parameter(mxy=100,mxk=6)
	common /par1/pi,srtnn,dy,dymid,ypcm,dm,bmin,bmax,isw
	common /parM2/nevent0,nevent,nev,ifileEv
	common /parM3/isym,iyrat
	common /IOevp/bminEv(20),bmaxEv(20),WeiEv(20)
     &		   ,nevEv(20),ipkEv(20),nfileEv,iprec(mxptcl),npt
	dimension DNDY(0:mxk,-mxy:mxy,mxptcl)
	character*(*) Xname,Yname,FileExt
c
	call recini(16,nevent0,1,FileExt,Xname,Yname)
	do i=-mxy,mxy
	  yrap=dy*i
	  if(iyrat.eq.1) yrap=dy/ypcm*i
	  write(16,801) yrap,(DNDY(IrecY,i,iprec(jp)),jp=1,npt)
	enddo
	close(16)
801	format(21(1x,1pe10.3))
	end
************************************************************************
	subroutine recini(ifile,nevent0,irectyp,FileExt,Xname,Yname)
************************************************************************
	implicit double precision(a-h,o-z)
	parameter(mxptcl=100)
	common /par1/pi,srtnn,dy,dymid,ypcm,dm,bmin,bmax,isw
	common /parM3/isym,iyrat
c...[ For IO
	common /IOch1/FnameEv(20)
	common /IOevp/bminEv(20),bmaxEv(20),WeiEv(20)
     &		   ,nevEv(20),ipkEv(20),nfileEv,iprec(mxptcl),npt
	character FnameEv*80
	common /headname/head
	character head*26
	save /headname/
c...] For IO
	character*(*) Xname,Yname,FileExt
c ---------------------------------------------------------------------*
c...Local Variables
	character fname*50,fnamei*80
	dimension chap(mxptcl)
	character chap*5
c ---------------------------------------------------------------------*
c       1     2     3      4     5     6      7      8      9       0
      data chap/
     & 'pi-','pi0','pi+' ,'K-' ,'K0b','K0'  ,'K+'  ,'0'   ,'0'    ,'0'
     1,'n'  ,'p'  ,'L'   ,'S-' ,'S0' ,'S+'  ,'Xi-' ,'Xi0' ,'Omg-' ,'0'
     2,'an' ,'ap' ,'aL'  ,'aS+','aS0','aS-' ,'aXi+','aXi0','aOmg+','0'
     3,'M-' ,'M0' ,'M+'  ,'B-' ,'B0' ,'B+'  ,'aB+' ,'aB0' ,'aB-'  ,'0'
     4,'l-' ,'lg0','l+'  ,'0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
     5,'Tot','h+-','h-'  ,'h0' ,'h+' ,'pi+-','K+-' ,'netp','p+ap' ,'N'
     6,'B'  ,'aB' ,'netB','0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
     7,'0'   ,'0' ,'0'   ,'0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
     8,'0'   ,'0' ,'0'   ,'0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
     9,'0'   ,'0' ,'0'   ,'0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
     &/
c
c	data (FileExt(k),Xname(k),Yname(k),k=nf1,nf2)/
c     &		 'mt' ,'dmT','dN/dmT^2dY'	! 16
c     &		,'pt' ,'pT' ,'dN/dpT^2dY'	! 17
c     &		,'elp','pT' ,'Elip(pT)'		! 18
c     &		,'mte','dmT','dN/dmT^2dE'	! 19
c     &		,'pte','pT' ,'dN/dpT^2dE'	! 20
c     &		,'epe','pT' ,'Elip(pT)'		! 21
c     &		,'rap','Y'  ,'dN/dY'		! 22
c     &		,'pav','Y'  ,'<p_T>(Y)'		! 23
c     &		,'dir','Y'  ,'<p_x>(Y)'		! 24
c     &		,'ely','Y'  ,'Elip(Y)'		! 25
c     &		,'bet','Y'  ,'beta_t(Y)'	! 26
c     &		,'eta','eta','dN/deta'		! 27
c     &		,'ele','eta','Elip(eta)'	! 28
c     &		,'btr','rt' ,'beta'		! 29
c     &		,'mid','Obs','Ave.'		! 30
c     &		,'v2c','pT' ,'v2c(pT)'		! 31
c     &	/
c ---------------------------------------------------------------------*
	fname=head(1:lengch(head))//'.'//FileExt(1:lengch(FileExt))
	open(ifile,file=fname,status='unknown')
	write(ifile,'(a,a,a,a)') '# ',Xname(1:lengch(Xname))
     &					,'-',Yname(1:lengch(Yname))
	write(ifile,'(a,1x,i5)') '# Total Events    =',nevent0
	write(ifile,'(a,f9.3)')  '# srtnn            =',srtnn
	write(ifile,'(a,f9.3)')  '# Y_cm(proj)       =',ypcm
	write(ifile,'(a,f9.3)')  '# dY(rapidity dist)=',dy
	write(ifile,'(a,f9.3)')  '# dM_t(mt     dist)=',dm
	write(ifile,'(a,f9.3)')  '# dY(mid-Y   cut)=+-',dymid
	write(ifile,'(a,f9.3)')  '# dY     /Yproj(cm)=',dy/ypcm
	write(ifile,'(a,f9.3)')  '# dY(mid)/Yproj(cm)=',dymid/ypcm
	write(ifile,'(a,f9.3)')  '# bmin             =',bmin
	write(ifile,'(a,f9.3)')  '# bmax             =',bmax
	if(iyrat.eq.1) then
	  write(ifile,'(a/a)')
     &		 '# rapidity/pseudorapidity is relative to ypcm'
	else
	  write(ifile,'(a/a)')
     &		 '# rapidity/pseudorapidity is absolute values'
	endif
	write(ifile,'(a/a)') '#','# Event Files'
	do jfile=1,nfileEv
	  fnamei=FnameEv(jfile)
	  write(ifile,'("# ",i2,": ",a)') jfile,fnamei(1:lengch(fnamei))
	enddo
	write(ifile,'(a)') '#'

	if(irectyp.eq.1) then
	  write(ifile,802) 1,Xname,(i+1,chap(iprec(i)),i=1,npt)
802	  format('# ',20(i2,':',a5,1x))
	elseif(irectyp.eq.0) then
	  write(ifile,'(a)') '#  0: dN/dY'
	  write(ifile,'(a)') '#  1: <Pt>'
	  write(ifile,'(a)') '#  2: v2'
	  write(ifile,'(a)') '#  3: <Mt>'
	  write(ifile,'(a)') '# 10: dN/deta'
	  write(ifile,'(a)') '# 11: <Pt>'
	  write(ifile,'(a)') '# 12: v2'
	  write(ifile,'(a)') '# 13: <Mt>'
	endif

	end
************************************************************************
	function id2q(id)
************************************************************************
*	id: PDG Monte Carlo code (input)
*	iq: charge
************************************************************************
	dimension IQq(6)
	save IQq
	data IQq / -1,2,-1,2,-1,2 /
c...
	id1=abs(id)
c...
	I1=mod(id1/10		,10)
	I2=mod(id1/100		,10)
	I3=mod(id1/1000		,10)
c...
	iq = 0		! Charge * 3

		if(I2.eq.0) then	! Leptons or gamma
	Ij=mod(id1		,10)
	if(I1.eq.1.and.mod(Ij,2).eq.1) iq=-1
		elseif(I3.eq.0) then	! Mesons
	Iodd=mod(I2,2)
	iq=(IQq(I2)-IQq(I1))/3
	if(Iodd.eq.1) iq=-iq
		else			! Baryons
	iq=(IQq(I1)+IQq(I2)+IQq(I3))/3
		endif
	if(id.lt.0) iq=-iq
	id2q=iq
	end
************************************************************************
      function lengch(char)
************************************************************************
      character*(*) char
c...Record
      nums=len(char)
      nchar=nums
         do 30 n=nums,1,-1
      nchar=n
      if(char(n:n).ne.' ') goto 31
 30      continue
      nchar=0
 31      continue
      lengch=nchar
      end
************************************************************************
	block data anndef
************************************************************************
	implicit double precision(a-h,o-z)
	parameter(mxy=100,mxk=6)
	parameter(mxptcl=100)
	parameter(mpT=100)
	common /par1/pi,srtnn,dy,dymid,ypcm,dm,bmin,bmax,isw
	common /par2/r0,r1,r2,r3,ipk
	common /par3/v2ev,nc,ncm
	common /dim0/ObsY(0:mxk)
	common /dim1/dmtY(mpT,mxptcl),dptY(mpT,mxptcl),elpY(mpT,mxptcl)
	common /dimY/dmtE(mpT,mxptcl),dptE(mpT,mxptcl),elpE(mpT,mxptcl)
	common /dimV/v1pY(mpT,mxptcl),v1pE(mpT,mxptcl)
	common /dim2/drp(0:mxk,-mxy:mxy,mxptcl)
	common /dim3/det(0:mxk,-mxy:mxy,mxptcl)
c...
	common /parM1/dy0,dymid0,dm0,bmin0,bmax0
	common /parM2/nevent0,nevent,nev,ifileEv
	common /parM3/isym,iyrat
c...[ For IO
	common /IOch1/FnameEv(20)
	common /IOevp/bminEv(20),bmaxEv(20),WeiEv(20)
     &		   ,nevEv(20),ipkEv(20),nfileEv,iprec(mxptcl),npt
	character FnameEv*80
c...] For IO
c ---------------------------------------------------------------------*
c...
	data dy0  ,dymid0,dm0  ,bmin0,bmax0 /
     &	     0.5d0,1.0d0 ,0.1d0,0.0d0,20.0d0/
	data npt,ipk,isym,iyrat/
     &	     18 ,3  ,0   ,0    /
c...
      data iprec/51,52,53,55,1,2,3,4,7,12,22,56,57,58,59,60,61,63,72*0/
c    &	 'Total','h+-','h-','h+','pi-','pi0','pi+','K-','K+','p','pbar'
c    &	,'pi+-','K+-','net p','p/pbar','N','B','net B'
c
	end
************************************************************************
	subroutine headline(jev,ntot,nbar,nmes,bimp,ipk,icon)
************************************************************************
	implicit double precision(a-h,o-z)
	character ch1*1,chline*150,ch2*2
	icon=0
c...[ Event Header line -----------------------------------------------*
c...[ PKS format
	if(ipk.eq.1) then
        read(15,*,end=200,err=200) jev,ntot,nbar,nmes,bimp
c...[ AO format
	else if(ipk.eq.0.or.ipk.eq.3) then
        read(15,'(a,a)',end=200,err=200) ch1,chline
	    if(ch1.eq.'#') then
        	read(chline,*) jev, ntot, nbar, nmes, bimp
	    else
		write(*,*) 'Something is funny.'
		stop
	    endif
c...[ NO format
        else if(ipk.eq.2)  then
        read(15,'(a,a)',end=200,err=200) ch2,chline
c...
	    if(ch2.eq.' #') then
        	read(chline,*) jev, ntot, nbar, nmes, bimp
	    else
		write(*,*) 'Something is funny.'
		stop
	    endif
c...
	endif
c...]]]
c...] Event Header line -----------------------------------------------*
	return

200	continue
	icon=999
	end
************************************************************************
	subroutine readptcl(id,em,p1,p2,p3,r1,r2,r3,r0,wei,ipk,icon)
************************************************************************
	implicit double precision(a-h,o-z)
	character chline*150,ch5*5
	icon=0
c...[ 
	if(ipk.eq.1) then
		read(15,*,end=999) id,ee,p1,p2,p3
		ems=max(0.0d0,ee**2-p1**2-p2**2-p3**2)
		em=sqrt(ems)
	elseif(ipk.eq.0) then
		read(15,*,end=999) id,em,p1,p2,p3
	elseif(ipk.eq.3) then
		read(15,*,end=999) id,em,p1,p2,p3,r1,r2,r3,r0
	elseif(ipk.eq.2) then
		read(15,'(a,a)',end=999) ch5,chline
		if(ch5(1:1).eq.'*') then
			icon=1
			return
		endif
		read(ch5,*) id
		read(chline,*) p1,p2,p3,ee, r1, r2, r3, r0
		ems=max(0.0d0,ee**2-p1**2-p2**2-p3**2)
		em=sqrt(ems)
	endif
c...]
	wei=1.0d0
c
	emLam=1.1157d0
c	emSp=1.1157d0
	emp=0.93830d0
	empic=0.13957d0
	empi0=0.13498d0
c
c Major baryons decaying to protons
c Lambda	3122 1.1157 0.639 (p pi-)
c Sigma+	3222 1.1894 0.5157 (p pi0)
c Sigma0	3212 1.1925 Lambda gamma 100 %
c Sigma-	3112 1.1974 0
c Xi0		3322 1.3149 Lambda pi0 100 %
c Xi-		3312 1.3213 Lambda pi- 100 %
c
	ida=abs(id)
	if(
     &	    ida.eq.3122
     &	.or.ida.eq.3222.or.ida.eq.3212
     &	.or.ida.eq.3322.or.ida.eq.3312
     &	) then
c...
	  if(ida.eq.3212.or.ida.eq.3322.or.ida.eq.3312) then
	    em1=emLam	! Lambda
	    em2=0.0d0			! Sigma0 --> Lambda gamma
	    if(ida.eq.3322) em2=empi0	! Xi0    --> Lambda pi0
	    if(ida.eq.3312) em2=empic	! Xi-    --> Lambda pi-
            call dec_iso(em,em1,em2,p1,p2,p3,px1,py1,pz1,icon)
	    em=em1
	    p1=px1
	    p2=py1
	    p3=pz1
	    id1=3122	! Lambda
	    if(id.lt.0) id1=-3122
	    id=id1
	    ida=3122
	  endif
c
	  em1=emp	! proton
	  em2=empic			! Lambda --> p pi-
	  if(ida.eq.3222) em2=empi0	! Sigma+ --> p pi0
          call dec_iso(em,em1,em2,p1,p2,p3,px1,py1,pz1,icon)
	  em=em1
	  p1=px1
	  p2=py1
	  p3=pz1
	  id1=2212
	  if(id.lt.0) id1=-2212
	  id=id1
	  wei=0.639d0			! B.R.(Lambda -> p pi-)
	  if(ida.eq.3222) wei=0.5157d0	! B.R.(Sigma+ --> p pi0)
	endif
	return

999	continue
	write(*,*) 'Read Error for particles'
	stop
	end
************************************************************************
	subroutine transDMT(DMT,n1,n2,mxptcl)
************************************************************************
	implicit double precision(a-h,o-z)
	dimension DMT(n1:n2,mxptcl)
c      data iprec/51,52,53,55,1,2,3,4,7,12,22,56,57,58,59,60,61,63,72*0/
cc    &	 'Total','h+-','h-','h+','pi-','pi0','pi+','K-','K+','p','pbar'
cc    &	,'pi+-','K+-','net p','p/pbar','N','B','net B'
cc
cc       1     2     3      4     5     6      7      8      9       0
c      data chap/
c     & 'pi-','pi0','pi+' ,'K-' ,'K0b','K0'  ,'K+'  ,'0'   ,'0'    ,'0'
c     1,'n'  ,'p'  ,'L'   ,'S-' ,'S0' ,'S+'  ,'Xi-' ,'Xi0' ,'Omg-' ,'0'
c     2,'an' ,'ap' ,'aL'  ,'aS+','aS0','aS-' ,'aXi+','aXi0','aOmg+','0'
c     3,'M-' ,'M0' ,'M+'  ,'B-' ,'B0' ,'B+'  ,'aB+' ,'aB0' ,'aB-'  ,'0'
c     4,'l-' ,'lg0','l+'  ,'0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
c     5,'Tot','h+-','h-'  ,'h0' ,'h+' ,'pi+-','K+-' ,'netp','p+ap' ,'N'
c     6,'B'  ,'aB' ,'netB','0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
c     7,'0'   ,'0' ,'0'   ,'0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
c     8,'0'   ,'0' ,'0'   ,'0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
c     9,'0'   ,'0' ,'0'   ,'0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
c     &/
c...[ Translation -----------------------------------------------------*
	do i=n1,n2
	DMT(i,31)=DMT(i,31)+DMT(i, 1)+DMT(i, 4)                    ! m-
	DMT(i,32)=DMT(i,32)+DMT(i, 2)+DMT(i, 5)+DMT(i, 6)          ! m0
	DMT(i,33)=DMT(i,33)+DMT(i, 3)+DMT(i, 7)                    ! m+
	DMT(i,34)=DMT(i,34)          +DMT(i,19)+DMT(i,14)+DMT(i,17)! b-
	DMT(i,35)=DMT(i,35)+DMT(i,11)+DMT(i,13)+DMT(i,15)+DMT(i,18)! b0
	DMT(i,36)=DMT(i,36)+DMT(i,12)          +DMT(i,16)          ! b+
	DMT(i,37)=DMT(i,37)+DMT(i,22)          +DMT(i,26)          ! a-
	DMT(i,38)=DMT(i,38)+DMT(i,21)+DMT(i,23)+DMT(i,25)+DMT(i,28)! a0
	DMT(i,39)=DMT(i,39)          +DMT(i,29)+DMT(i,24)+DMT(i,27)! a+
c
	DMT(i,53)=DMT(i,31)+DMT(i,34)+DMT(i,37)	! h-
	DMT(i,54)=DMT(i,32)+DMT(i,35)+DMT(i,38)	! h0
	DMT(i,55)=DMT(i,33)+DMT(i,36)+DMT(i,39)	! h+
	DMT(i,52)=DMT(i,53)+DMT(i,55)		! h+-
	DMT(i,51)=DMT(i,52)+DMT(i,54)		! hadrons
     &           +DMT(i,41)+DMT(i,42)+DMT(i,43)	! lep/gam
c
	DMT(i,56)=DMT(i,1)+DMT(i,3)		! pi+-
	DMT(i,57)=DMT(i,4)+DMT(i,7)		! K+-
	DMT(i,58)=DMT(i,12)-DMT(i,22)		! net p
	DMT(i,59)=DMT(i,12)+DMT(i,22)		! p/pbar
	DMT(i,60)=DMT(i,11)+DMT(i,12)		! N
c
	DMT(i,61)=DMT(i,34)+DMT(i,35)+DMT(i,36)	! B
	DMT(i,62)=DMT(i,37)+DMT(i,38)+DMT(i,39)	! aB
	DMT(i,63)=DMT(i,61)-DMT(i,62)		! net B
c
	enddo
	end
************************************************************************
	subroutine transDNDY(DNDY,mxk,mxy,mxptcl)
************************************************************************
	implicit double precision(a-h,o-z)
	dimension DNDY(0:mxk,-mxy:mxy,mxptcl)
c      data iprec/51,52,53,55,1,2,3,4,7,12,22,56,57,58,59,60,61,63,72*0/
cc    &	 'Total','h+-','h-','h+','pi-','pi0','pi+','K-','K+','p','pbar'
cc    &	,'pi+-','K+-','net p','p/pbar','N','B','net B'
cc
cc       1     2     3      4     5     6      7      8      9       0
c      data chap/
c     & 'pi-','pi0','pi+' ,'K-' ,'K0b','K0'  ,'K+'  ,'0'   ,'0'    ,'0'
c     1,'n'  ,'p'  ,'L'   ,'S-' ,'S0' ,'S+'  ,'Xi-' ,'Xi0' ,'Omg-' ,'0'
c     2,'an' ,'ap' ,'aL'  ,'aS+','aS0','aS-' ,'aXi+','aXi0','aOmg+','0'
c     3,'M-' ,'M0' ,'M+'  ,'B-' ,'B0' ,'B+'  ,'aB+' ,'aB0' ,'aB-'  ,'0'
c     4,'l-' ,'lg0','l+'  ,'0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
c     5,'Tot','h+-','h-'  ,'h0' ,'h+' ,'pi+-','K+-' ,'netp','p+ap' ,'N'
c     6,'B'  ,'aB' ,'netB','0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
c     7,'0'   ,'0' ,'0'   ,'0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
c     8,'0'   ,'0' ,'0'   ,'0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
c     9,'0'   ,'0' ,'0'   ,'0'  ,'0'  ,'0'   ,'0'   ,'0'   ,'0'    ,'0'
c     &/
c...[ Translation -----------------------------------------------------*
	do i=-mxy,mxy
	do L=0,mxk
c
	DNDY(L,i,31)=DNDY(L,i,31)+DNDY(L,i, 1)+DNDY(L,i, 4)            ! m-
	DNDY(L,i,32)=DNDY(L,i,32)+DNDY(L,i, 2)+DNDY(L,i, 5)+DNDY(L,i, 6)! m0
	DNDY(L,i,33)=DNDY(L,i,33)+DNDY(L,i, 3)+DNDY(L,i, 7)            ! m+
	DNDY(L,i,34)=DNDY(L,i,34)+DNDY(L,i,14)+DNDY(L,i,17)+DNDY(L,i,19)! b-
	DNDY(L,i,35)=DNDY(L,i,35)+DNDY(L,i,15)+DNDY(L,i,18)
     &             +DNDY(L,i,11)+DNDY(L,i,13)				! b0
	DNDY(L,i,36)=DNDY(L,i,36)+DNDY(L,i,16)
     &             +DNDY(L,i,12)					! b+
	DNDY(L,i,37)=DNDY(L,i,37)+DNDY(L,i,26)				! aB-
     &             +DNDY(L,i,22)
	DNDY(L,i,38)=DNDY(L,i,38)+DNDY(L,i,25)+DNDY(L,i,28)
     &             +DNDY(L,i,21)+DNDY(L,i,23)				! aB0
	DNDY(L,i,39)=DNDY(L,i,39)+DNDY(L,i,24)+DNDY(L,i,27)+DNDY(L,i,29)! aB+
c
	DNDY(L,i,53)=DNDY(L,i,31)+DNDY(L,i,34)+DNDY(L,i,37)	! h-
	DNDY(L,i,54)=DNDY(L,i,32)+DNDY(L,i,35)+DNDY(L,i,38)	! h0
	DNDY(L,i,55)=DNDY(L,i,33)+DNDY(L,i,36)+DNDY(L,i,39)	! h+
	DNDY(L,i,52)=DNDY(L,i,53)+DNDY(L,i,55)		! h+-
	DNDY(L,i,51)=DNDY(L,i,52)+DNDY(L,i,54)		! hadrons
     &           +DNDY(L,i,41)+DNDY(L,i,42)+DNDY(L,i,43)	! lep/gam
c
	DNDY(L,i,56)=DNDY(L,i,1)+DNDY(L,i,3)		! pi+-
	DNDY(L,i,57)=DNDY(L,i,4)+DNDY(L,i,7)		! K+-
	DNDY(L,i,58)=DNDY(L,i,12)-DNDY(L,i,22)		! net p
	DNDY(L,i,59)=DNDY(L,i,12)+DNDY(L,i,22)		! p/pbar
	DNDY(L,i,60)=DNDY(L,i,11)+DNDY(L,i,12)		! N
c
	DNDY(L,i,61)=DNDY(L,i,34)+DNDY(L,i,35)+DNDY(L,i,36)	! B
	DNDY(L,i,62)=DNDY(L,i,37)+DNDY(L,i,38)+DNDY(L,i,39)	! aB
	DNDY(L,i,63)=DNDY(L,i,61)-DNDY(L,i,62)		! net B
c
	enddo
	enddo
c...] Translation -----------------------------------------------------*
	end
************************************************************************
	subroutine addbin(id,em,p1,p2,p3,wei,rap,eta,iQ)
************************************************************************
	implicit double precision(a-h,o-z)
	parameter(mxy=100,mxk=6)
	parameter(mxptcl=100)
	parameter(mpT=100)
	common /par1/pi,srtnn,dy,dymid,ypcm,dm,bmin,bmax,isw
	common /par2/r0,r1,r2,r3,ipk
	common /par3/v2ev,nc,ncm
	common /dim0/ObsY(0:mxk)
	common /dim1/dmtY(mpT,mxptcl),dptY(mpT,mxptcl),elpY(mpT,mxptcl)
	common /dimY/dmtE(mpT,mxptcl),dptE(mpT,mxptcl),elpE(mpT,mxptcl)
	common /dimV/v1pY(mpT,mxptcl),v1pE(mpT,mxptcl)
	common /dim2/drp(0:mxk,-mxy:mxy,mxptcl)
	common /dim3/det(0:mxk,-mxy:mxy,mxptcl)
	common /dim6/vcsum(0:2,0:mxptcl),sumNCS(0:2,0:mxptcl)
c ---------------------------------------------------------------------*

c...[ Kinematics
	pts=p1**2+p2**2
	ps=pts+p3**2
	es=em**2+ps
	e=sqrt(es)
c...Rapidity
	rap=log((e+p3)/(e-p3))/2
	irap=nint(rap/dy)
c...Pseudo Rapidity
	pp=sqrt(ps)
	costht=p3/pp
	if(costht.lt.1.0d0.and.costht.gt.-1.0d0) then
		eta=log((1.0d0+costht)/(1.0d0-costht))/2
		ieta=nint((eta+dy*mxy*2)/dy)-2*mxy
	else 
		eta=10000.0d0
		if(costht.lt.0) eta=-eta
		ieta=mxy*2	! out-of-bins
	endif

	if(abs(irap).gt.mxy.and.abs(ieta).gt.mxy
     &	.and.dymid.ne.0
     &	.and.abs(eta).gt.dymid
     &	.and.abs(rap).gt.dymid
     &	) return 	! Out-of-Bins
c...[
	call partype(id,ip,iQ)
	pt=sqrt(pts)
	emt=sqrt(em**2+pts)
c	rr=sqrt(r1**2+r2**2)
	ObsY(0)=1.0d0			! Rapidity Dist.
	ObsY(1)=pt			! <Pt>
	ObsY(2)=p1			! <Px>: Directed Flow
	ObsY(3)=(p1**2-p2**2)/pts	! <(Px**2-Py**2)/Pt**2>: v2
	ObsY(4)=p1/pt			! <Px/Pt>: v1
	ObsY(5)=0.0d0
	ObsY(6)=0.0d0
c...[ Event Analysis
c	Iadd=1
	if(iQ.ne.0) nc=nc+1
c...] Event Analysis

c...[ Rapidity/Pseudo Rapidity Distribution ---------------------------*
	if(abs(irap).gt.mxy.and.abs(ieta).gt.mxy) goto 3000
	if(abs(irap).le.mxy) then
	   do L=0,mxk
	      drp(L, irap,ip)=drp(L, irap,ip)+wei*ObsY(L)
	   enddo
	endif
	if(abs(ieta).le.mxy) then
	   do L=0,mxk
	      det(L, ieta,ip)=det(L, ieta,ip)+wei*ObsY(L)
	   enddo
	endif
c...] Rapidity/Pseudo Rapidity Distribution ---------------------------*

3000	continue
c...[ MidRapidity Observables -----------------------------------------*
c...[ Pt or Mt Distribution -------------------------------------------*
	imt=int((emt-em)/dm)+1
	ipt=int(pt/dm)+1
	delp=(p1**2-p2**2)/(p1**2+p2**2)
	pdir=p1
	if(p3.lt.0) pdir=-p1
	if(dymid.eq.0.or.abs(rap).lt.dymid) then
	  if(imt.le.100) dmtY(imt,ip)=dmtY(imt,ip)+wei/emt
	  if(ipt.le.100) dptY(ipt,ip)=dptY(ipt,ip)+wei/pt
	  if(ipt.le.100) elpY(ipt,ip)=elpY(ipt,ip)+wei*delp/pt
	  if(ipt.le.100) v1pY(ipt,ip)=v1pY(ipt,ip)+wei*pdir/pt/pt
	endif
	if(dymid.eq.0.or.abs(eta).lt.dymid) then
	  if(imt.le.100) dmtE(imt,ip)=dmtE(imt,ip)+wei/emt
	  if(ipt.le.100) dptE(ipt,ip)=dptE(ipt,ip)+wei/pt
	  if(ipt.le.100) elpE(ipt,ip)=elpE(ipt,ip)+wei*delp/pt
	  if(ipt.le.100) v1pE(ipt,ip)=v1pE(ipt,ip)+wei*pdir/pt/pt
	  if(iQ.ne.0) then
	    ncm=ncm+1
	    v2ev=v2ev+ObsY(3)*wei
	    if(ipt.le.100) then	! AO:050427
	      phi=atan2(p2,p1)
	      sumNCS(0,ipt)=sumNCS(0,ipt)+wei
	      sumNCS(1,ipt)=sumNCS(1,ipt)+wei*cos(2*phi)
	      sumNCS(2,ipt)=sumNCS(2,ipt)+wei*sin(2*phi)
	    endif
	  endif
	endif

c...] Pt or Mt Distribution -------------------------------------------*
c...] MidRapidity Observables -----------------------------------------*
	end
************************************************************************
	subroutine partype(id,iptyp,iQ)
************************************************************************
	implicit double precision(a-h,o-z)
	save Ilep,Imes,Ibar,Iaba
	data Ilep,Imes,Ibar,Iaba/42,32,35,38/
c...
c        1    2     3    
c	 pi-  pi0   pi+  
c        4    5     6     7
c        K-   K0b   K0    K+
c        11   12    13    14    15    16    17    18
c        n    p     Lam   Sig-  Sig0  Sig+  Xi-   Xi0
c        21   22    23    24    25    26    27    28
c        an   ap    aLam  aSig- aSig0 aSig+ aXi-  aXi0
c        31   32    33
c        m-   m0    m+
c        34   35    36
c        b-   b0    b+
c        37   38    39
c        ab-  ab0   ab+
c        41   42    43
c        l-   lg0   l+
c
c	 1: pi-
c	 2: pi0
c	 3: pi+
c	 4: K-
c	 5: K0bar
c	 6: K0
c	 7: K+
c	11: n
c	12: p
c	13: Lam
c	14: Sigma-
c	15: Sigma0
c	16: Sigma+
c	17: Xi0
c	18: Xi-
c	19: Omega-
c	21: an
c	22: ap
c	23: aLam
c	24: aSigma+
c	25: aSigma0
c	26: aSigma-
c	27: aXi0
c	28: aXi+
c	29: aOmg+
c
c	31: Other Negative Mesons
c	32: Other Neutral  Mesons
c	33: Other Positive Mesons
c	34: Other Negative Baryons
c	35: Other Neutral  Baryons
c	36: Other Positive Baryons
c	37: Other Negative Anti-Baryons
c	38: Other Neutral  Anti-Baryons
c	39: Other Positive Anti-Baryons
c
c	41: Negative Leptons
c	42: Neutral  Leptons and Gammas
c	43: Positive Leptons
c...[
	Irec=0
c
	iQ=0
	Isft=Ilep
	if(abs(id).lt.100) goto 1000	! leptons and gammas
	imbx=abs(id)/1000
	imb=mod(imbx,10)+1		! 1: Meson, 2-7:Baryons
	if(imbx.ge.10) then		! Excited Hadrons, Should Not Exist
		Irec=1
		Isft=Imes
		if(imb.gt.1) then
		   Isft=Ibar
		   if(id.lt.0) Isft=Iaba
		endif
		goto 1000
	endif
	goto (200,300,300,400,500,500,500) imb ! Meson,N,N,Y,C-bar.,B-bar.
c...]

c...[ Particle Type ---------------------------------------------------*
c...[ Mesons
 200	continue
	iflv=mod(abs(id)/100,10)
	if(iflv.eq.1) then	! pi0
c        1    2     3    
c	 pi-  pi0   pi+  
	   if(id.eq.111) then
	        iptyp=2
	   else
		Isft=Imes
		Irec=1
		goto 1000
	   endif
	else if(iflv.eq.2) then	! pi+-
c        1    2     3    
c	 pi-  pi0   pi+  
	   if(id.eq.-211) then		! pi-
	        iptyp=1
		iQ=-1
	   elseif(id.eq.211) then	! pi+
	        iptyp=3
		iQ=1
	   elseif(id.eq.221) then	! eta, stable... 
		iptyp=Imes
	   else
		Isft=Imes
		Irec=1
		goto 1000
	   endif
c...Kaons
	else if(iflv.eq.3) then	! K
c        4    5     6     7
c        K-   K0b   K0    K+
	   if(id.eq.-321) then	! K-
		iptyp=4
		iQ=-1
	   elseif(id.eq.-311) then	! K0bar
		iptyp=5
	   elseif(id.eq.311) then	! K0
		iptyp=6
	   elseif(id.eq.321) then		! K+
		iptyp=7
		iQ=1
	   else
		Isft=Imes
		Irec=1
		goto 1000
	   endif

c...Heavy Meson with c,b,t
	else
c        31   32    33
c        m-   m0    m+
		Isft=Imes
		goto 1000
	endif
	return
c...] Mesons

c...[ Nucleon (Light-Flavored Baryon)
 300	continue
c        11   12    13    14    15    16    17    18
c        n    p     Lam   Sig-  Sig0  Sig+  Xi-   Xi0
c        21   22    23    24    25    26    27    28
c        an   ap    aLam  aSig- aSig0 aSig+ aXi-  aXi0
c        34   35    36
c        b-   b0    b+
c        37   38    39
c        ab-  ab0   ab+
	jd=abs(id)
	if(jd.eq.2112) then
		iptyp=11
	elseif(jd.eq.2212) then
		iptyp=12
		iQ=1
	else
		Isft=Ibar
		if(id.lt.0) Isft=Iaba
		Irec=1
		goto 1000
	endif
	if(id.lt.0) then
		iptyp=iptyp+10
		iQ=-iQ
	endif
	return
c...] Nucleon (Light-Flavored Baryon)

c...[ Hyperons
 400	continue
c        11   12    13    14    15    16    17    18
c        n    p     Lam   Sig-  Sig0  Sig+  Xi-   Xi0
c        21   22    23    24    25    26    27    28
c        an   ap    aLam  aSig- aSig0 aSig+ aXi-  aXi0
c        34   35    36
c        b-   b0    b+
c        37   38    39
c        ab-  ab0   ab+
	jd=abs(id)
	if(jd.eq.3122) then		! Lam
		iptyp=13
	elseif(jd.eq.3112) then		! Sigma-
		iptyp=14
		iQ=-1
	elseif(jd.eq.3212) then		! Sigma0
		iptyp=15
	elseif(jd.eq.3222) then		! Sigma+
		iptyp=16
		iQ=1
	elseif(jd.eq.3312) then		! Xi-
		iptyp=17
		iQ=-1
	elseif(jd.eq.3322) then		! Xi0
		iptyp=18
	elseif(jd.eq.3334) then		! Omega-
		iptyp=19
		iQ=-1
	else
		Isft=Ibar
		if(id.lt.0) Isft=Iaba
		Irec=1
		goto 1000
	endif
	if(id.lt.0) then
		iptyp=iptyp+10
		iQ=-iQ
	endif
	return
c...] Hyperons

c...[ Heavy-Flavor Baryons
 500	continue
c        11   12    13    14    15    16    17    18
c        n    p     Lam   Sig-  Sig0  Sig+  Xi-   Xi0
c        21   22    23    24    25    26    27    28
c        an   ap    aLam  aSig- aSig0 aSig+ aXi-  aXi0
c        34   35    36
c        b-   b0    b+
c        37   38    39
c        ab-  ab0   ab+
	Isft=Ibar
	if(id.lt.0) Isft=Iaba
	goto 1000
c...] Heavy-Flavor Baryons

c...[ Unspecified
 1000	continue
	iQ=id2q(id)
	iptyp=iQ+Isft
	if(Irec.eq.1) then
		write(*,*) 'Funny ID, id=',id
	endif

	end
************************************************************************
	subroutine vcann(ic)
************************************************************************
	implicit double precision(a-h,o-z)
	parameter(mxptcl=100)
	common /par1/pi,srtnn,dy,dymid,ypcm,dm,bmin,bmax,isw
	common /parM2/nevent0,nevent,nev,ifileEv
	common /dim6/vcsum(0:2,0:mxptcl),sumNCS(0:2,0:mxptcl)
	common /v2corrEV/partEV,pairEV,v2cEV
	character*1 ch1

c...ic=1: New Incident Energy
c...ic=2: At the Beginning of Each Event
c...ic=4: At the End       of Each Event
c...ic=5: At the End of Each Incident Energy
	partEV=0.0d0
	pairEV=0.0d0
	v2cEV=0.0d0
	goto (1000,2000,3000,4000,5000) ic
	goto 3000
c ---------------------------------------------------------------------*

 1000	continue
	do k=0,2
	  do ipt=0,mxptcl
	    vcsum(k,ipt)=0.0d0
	  enddo
	enddo
	return

 2000	continue
	do k=0,2
	  do ipt=0,mxptcl
	    sumNCS(k,ipt)=0.0d0
	  enddo
	enddo
	return

 4000	continue
	do k=0,2
	    sumNCS(k,0)=0.0d0
	  do ipt=1,mxptcl
	    sumNCS(k,0)=sumNCS(k,0)+sumNCS(k,ipt)
	  enddo
	enddo
	do ipt=0,mxptcl
	    pair=sumNCS(0,ipt)*(sumNCS(0,0)-1.0d0)
	  if(pair.gt.0) then
	    vcsum(0,ipt)=vcsum(0,ipt)+pair
	    dsum=sumNCS(1,ipt)*sumNCS(1,0)
     &          +sumNCS(2,ipt)*sumNCS(2,0)
     &          -sumNCS(0,ipt)
	    vcsum(1,ipt)=vcsum(1,ipt)+dsum
	    vcsum(2,ipt)=vcsum(2,ipt)+dsum**2
	  endif
	enddo
	dsumEV=sumNCS(1,0)**2+sumNCS(2,0)**2-sumNCS(0,0)
	pairEV=sumNCS(0,0)*(sumNCS(0,0)-1.0d0)
	partEV=sumNCS(0,0)
	v2cEV=0.0d0
	if(pairEV.gt.0.and.dsumEV.ge.0) then
	  v2cEV=sqrt(dsumEV/pairEV)
	endif
	return

 5000	continue
	pair0=vcsum(0,0)
	v20=0.0d0
	if(pair0.gt.0.and.vcsum(1,0).gt.0) then
	  v20=sqrt(vcsum(1,0)/pair0)
	endif
	write(31,'("#  pT         v2(corr)   D-v2(corr) <Pair>")')
	do ipt=0,mxptcl
	  pt=ipt*dm-dm/2
	  pair=vcsum(0,ipt)
	  v2=0.0d0
	  v2sgm=0.0d0
	  if(pair.gt.0.and.v20.gt.0) then
		v2=vcsum(1,ipt)/pair/v20
		v2sgm=vcsum(2,ipt)/pair-(vcsum(1,ipt)/pair)**2
		v2sgm=sqrt(v2sgm)/v20/sqrt(pair)
	  endif
	  ch1=' '
	  if(ipt.eq.0) ch1='#'
	  write(31,'(a1," ",10(1pe10.3,1x))')
     &		 ch1,pt,v2,v2sgm,vcsum(0,ipt)/nevent0
	enddo
	return

 3000	continue

	end
c***********************************************************************
      subroutine dec_iso(em0,em1,em2,px0,py0,pz0,px1,py1,pz1,icon)
c***********************************************************************
c...Purpose: to perform isotropic binary decay
      implicit double precision(a-h, o-z)
      save pi,ifirst
      data ifirst/1/
      if(ifirst.eq.1) then
         ifirst=0
         pi=4*atan(1.0d0)
      endif
      icon=999
      px1=0.0d0
      py1=0.0d0
      pz1=0.0d0
      if(em1+em2.gt.em0) return
      pcmq=(em0*em0-(em1+em2)**2)*(em0*em0-(em1-em2)**2)
      pcmq=sqrt(pcmq)/(2.d0*em0)
      costh=1.0d0-2*rn(0)
      sinth=sqrt(1.0d0-costh**2)
      phi=2*pi*rn(0)
      px1=pcmq*sinth*cos(phi)
      py1=pcmq*sinth*sin(phi)
      pz1=pcmq*costh
      ee1=sqrt(em1**2+px1**2+py1**2+pz1**2)
      call boost(px0,py0,pz0,em0,px1,py1,pz1,ee1)
      end
c***********************************************************************

      subroutine boost(pxboost,pyboost,pzboost,emboost,p1,p2,p3,p4)    
    
c...Purpose: to perform boost (p1,p2,p3,p4(=E))

      implicit double precision(a-h, o-z)
      dimension dp(4)
    
c...Boost, typically from rest to momentum/energy=beta. 
      psboost=pxboost**2+pyboost**2+pzboost**2
      eeboost=sqrt(emboost**2+psboost)
      bx=pxboost/eeboost
      by=pyboost/eeboost
      bz=pzboost/eeboost
      ga=eeboost/emboost	! =1d0/sqrt(1d0-beta**2) 
      dp(1)=p1
      dp(2)=p2
      dp(3)=p3
      dp(4)=p4
      bp=bx*dp(1)+by*dp(2)+bz*dp(3)   
      gabp=ga*(ga*bp/(1.0d0+ga)+dp(4)) 
      p1=dp(1)+gabp*bx  
      p2=dp(2)+gabp*by  
      p3=dp(3)+gabp*bz  
      p4=ga*(dp(4)+bp)  
      end   
C*********************************************************************
      block data randseed
      implicit double precision(a-h, o-z)
      common/pydatr/mrpy(6),rrpy(100)
      data mrpy/19780503,0,0,97,33,0/
      end
C*********************************************************************
C...PYRND
C...Generates random numbers uniformly distributed between
C...0 and 1, excluding the endpoints.
 
      function pyrnd(idummy)
 
C...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
C...Commonblocks.
      common/pydatr/mrpy(6),rrpy(100)
      save /pydatr/
C...Equivalence between commonblock and local variables.
      equivalence (mrpy1,mrpy(1)),(mrpy2,mrpy(2)),(mrpy3,mrpy(3)),
     &(mrpy4,mrpy(4)),(mrpy5,mrpy(5)),(mrpy6,mrpy(6)),
     &(rrpy98,rrpy(98)),(rrpy99,rrpy(99)),(rrpy00,rrpy(100))
 
C...Initialize generation from given seed.
      if(mrpy2.eq.0) then
        ij=mod(mrpy1/30082,31329)
        kl=mod(mrpy1,30082)
        i=mod(ij/177,177)+2
        j=mod(ij,177)+2
        k=mod(kl/169,178)+1
        l=mod(kl,169)
        do 110 ii=1,97
          s=0d0
          t=0.5d0
          do 100 jj=1,48
            m=mod(mod(i*j,179)*k,179)
            i=j
            j=k
            k=m
            l=mod(53*l+1,169)
            if(mod(l*m,64).ge.32) s=s+t
            t=0.5d0*t
  100     continue
          rrpy(ii)=s
  110   continue
        twom24=1d0
        do 120 i24=1,24
          twom24=0.5d0*twom24
  120   continue
        rrpy98=362436d0*twom24
        rrpy99=7654321d0*twom24
        rrpy00=16777213d0*twom24
        mrpy2=1
        mrpy3=0
        mrpy4=97
        mrpy5=33
      endif
 
C...Generate next random number.
  130 runi=rrpy(mrpy4)-rrpy(mrpy5)
      if(runi.lt.0d0) runi=runi+1d0
      rrpy(mrpy4)=runi
      mrpy4=mrpy4-1
      if(mrpy4.eq.0) mrpy4=97
      mrpy5=mrpy5-1
      if(mrpy5.eq.0) mrpy5=97
      rrpy98=rrpy98-rrpy99
      if(rrpy98.lt.0d0) rrpy98=rrpy98+rrpy00
      runi=runi-rrpy98
      if(runi.lt.0d0) runi=runi+1d0
      if(runi.le.0d0.or.runi.ge.1d0) goto 130
 
C...Update counters. Random number to output.
      mrpy3=mrpy3+1
      if(mrpy3.eq.1000000000) then
        mrpy2=mrpy2+1
        mrpy3=0
      endif
      pyrnd=runi
 
      return
      end
c***********************************************************************

      subroutine jamcpu(iunit,ic,iseed,isec)

c...Purpose:   count the cpu time  U77 library
c...Variables: IC     -  0->start, 1->end
c...For DTIME see  %man 3f dtime

      external dtime
      real tarray(2)
      character today*9,timestr*8
      data cptime1/0.0/
      data tarray/0.0e0,0.0e0/
      save  stime,cptime1

      if(ic.eq.0) then
        cptime1=dtime(tarray)
        call date_and_time(today)
        call time(timestr)
        write(iunit,'(''***********************************'')')
        write(iunit,'(''Starting time = '',A8,''  '',A9)')timestr,today
        stime=secnds(0.0)
c...Set a different random seed each time 
        if(iseed.eq.0) then
          read(timestr,100) iseed1,iseed2,iseed3
  100     format(i2,1x,i2,1x,i2)
          iseed = iseed1*10000 + iseed2*100 + iseed3
          if(iseed/2*2.eq.iseed) iseed = iseed + 1
        end if
      end if

C...Compute Elapse and CPU time
      if(ic.eq.1) then
        cptime2=dtime(tarray)
        call date_and_time(today)
        call time(timestr)
        write(iunit,'(''  ending time = '',a8,''  '',a9)')timestr,today
        stime=secnds(stime)
        isec=stime
        imin=isec/60
        ihrs=imin/60
        imin=imin-ihrs*60
        isec=isec-ihrs*3600-imin*60
        write(iunit,901)ihrs,imin,isec
901     format(' * Elapse time =',i3,' h ',i3,' m ',i3,' s')
        isec=cptime2-cptime1
        imin=isec/60
        ihrs=imin/60
        imin=imin-ihrs*60
        isec=isec-ihrs*3600-imin*60
        write(iunit,902)ihrs,imin,isec
902     format(' *    CUP time =',i3,' h ',i3,' m ',i3,' s')
        write(iunit,'(''***********************************'')')
        isec=cptime2-cptime1
      end if

      end

c***********************************************************************

      function rn(idumm)

c...Purpose:  link random number generator.
      double precision rn,pyrnd,rrpy
      common/pydatr/mrpy(6),rrpy(100)
      common/rseed/iseed
      save /rseed/
c...Program must be compiled with either +e or +E1 option
c...in order to use  function RAN on HP.
c...for HP
c     rn = ran(iseed)
      if(idumm.eq.-1) then
	iseed=0
        call jamcpu(6,0,iseed,isec)
        mrpy(1)=iseed
      endif
      rn = pyrnd(0)    ! random number from pythia
      end
