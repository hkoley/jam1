************************************************************************
	program phsann
************************************************************************
*     Purpose:
*	Analyze jam phase space data.
*     Ver: 1.22
*	Calculate 
*	Unit	Filename	X	Y
*     	16	$HEAD.mt	dmT	dN/dmT^2dY
*     	17	$HEAD.pt	pT	dN/dpT^2dY
*     	18	$HEAD.elp	pT	Elip(pT)
*     	19	$HEAD.mte	dmT	dN/dmT^2dE
*     	20	$HEAD.pte	pT	dN/dpT^2dE
*     	21	$HEAD.epe	pT	Elip(pT)
*     	22	$HEAD.rap	Y	dN/dY
*     	23	$HEAD.pav	Y	<p_T>(Y)
*     	24	$HEAD.dir	Y	<p_x>(Y)
*     	25	$HEAD.ely	Y	Elip(Y)
*     	26	$HEAD.bet	Y	beta_t(Y)
*     	27	$HEAD.eta	eta	dN/deta
*     	28	$HEAD.ele	eta	Elip(eta)
*     	29	$HEAD.btr	rt	beta
*     	30	$HEAD.mid	Obs	Ave.
*     	31	$HEAD.v2c	pT	v2c(pT)
*************************************************************************
	implicit real*8(a-h,o-z)
	parameter(mxy=100,mxk=6)
	parameter(mpt=100)
	common /par1/pi,srtnn,dy,dymid,ypcm,deta,dm,bmin,bmax,rtmax,dlr
	common /par2/r0,r1,r2,r3,ipk
	common /par3/v2ev,nc,ncm
	common /dim0/ObsY(0:mxk)
	common /dim1/dmtY(100,mpt),dptY(100,mpt),elpY(100,mpt)
	common /dimY/dmtE(100,mpt),dptE(100,mpt),elpE(100,mpt)
	common /dim2/drp(0:mxk,-mxy:mxy,mpt)
	common /dim3/det(0:mxk,-mxy:mxy,mpt)
	common /dim4/betr(6,100)
	common /dim5/drpm(0:mxk,mpt),detm(0:mxk,mpt)
c
	common /partinf/em,p1,p2,p3,rap,eta,id,iQ
c...
	common /parM1/dy0,dymid0,dm0,bmin0,bmax0,rtmax0
	common /parM2/nevent0,nevent,nev,ifileEv
	common /parM3/isym,inext,iyrat
c...[ For IO
	common /IOch1/FnameEv(20)
	common /IOevp/bminEv(20),bmaxEv(20),WeiEv(20)
     &		   ,nevEv(20),ipkEv(20),nfileEv,iprec(mpt),npt
	character FnameEv*80
c...] For IO
c ---------------------------------------------------------------------*
c...Local Variables
	character ch1*1
c ---------------------------------------------------------------------*
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
	      call readptcl(id,em,p1,p2,p3,r1,r2,r3,r0,ipk,iok)
	    if(icon.eq.0) then
	      call addbin(id,em,p1,p2,p3,rap,eta,iQ)
c	      call vcann(3)
	    endif
	   enddo
	      call vcann(4)
	      v2ev1=0.0d0
	      if(ncm.gt.0) v2ev1=v2ev/ncm
	      write(42,811) bimp,v2ev1,ncm,nc,ntot,nbar,nmes
811	      format(2(1x,1pe10.3),5(1x,i5))
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
	call vcann(5)
c...[ Midrapidity -----------------------------------------------------*

c...[ Mt, Pt Dist.
c...dN/dY/Mt/dMt/(2*pi)
c...dN/dY/Pt/dPt/(2*pi)
	call transDMT(dmtY,1,100,mpt)
	call transDMT(dptY,1,100,mpt)
	call transDMT(elpY,1,100,mpt)
	call transDMT(dmtE,1,100,mpt)
	call transDMT(dptE,1,100,mpt)
	call transDMT(elpE,1,100,mpt)
	dymid1=dymid
	if(dymid.eq.0) dymid1=0.5d0
		do i=1,100
	demt=i*dm-dm/2
c	emtpr=demt+emnuc
c	emtpi=demt+empi
		do jp=1,npt
	ip=iprec(jp)
	dmtY(i,ip)=dmtY(i,ip)/pi/2/dm/dymid1/2/nevent0
	dptY(i,ip)=dptY(i,ip)/pi/2/dm/dymid1/2/nevent0
	elpY(i,ip)=elpY(i,ip)/pi/2/dm/dymid1/2/nevent0
	dmtE(i,ip)=dmtE(i,ip)/pi/2/dm/dymid1/2/nevent0
	dptE(i,ip)=dptE(i,ip)/pi/2/dm/dymid1/2/nevent0
	elpE(i,ip)=elpE(i,ip)/pi/2/dm/dymid1/2/nevent0
	if(dptY(i,ip).ne.0) elpY(i,ip)=elpY(i,ip)/dptY(i,ip)
	if(dptE(i,ip).ne.0) elpE(i,ip)=elpE(i,ip)/dptE(i,ip)
		enddo
	write(16,801) demt,(dmtY(i,iprec(jp)),jp=1,npt)
	write(17,801) demt,(dptY(i,iprec(jp)),jp=1,npt)
	write(18,801) demt,(elpY(i,iprec(jp)),jp=1,npt)
	write(19,801) demt,(dmtE(i,iprec(jp)),jp=1,npt)
	write(20,801) demt,(dptE(i,iprec(jp)),jp=1,npt)
	write(21,801) demt,(elpE(i,iprec(jp)),jp=1,npt)
		enddo
c...] Mt, Pt Dist.

c...[ Mid-Rapidity Integrated Observables (dN/dY, <Pt>, v2, ...)
	call transDMT(drpm,0,mxk,mpt)
	call transDMT(detm,0,mxk,mpt)
	write(41,'(a,f9.3)') '# srtnn            =',srtnn
c
	do jp=1,npt
	ip=iprec(jp)
c...Normalization
		do L=0,mxk
	drpm(L,ip)=drpm(L,ip)/dymid1/2/nevent0
	detm(L,ip)=detm(L,ip)/dymid1/2/nevent0
		enddo
		do L=1,mxk
	if(drpm(0,ip).ne.0) drpm(L,ip)=drpm(L,ip)/drpm(0,ip)
	if(detm(0,ip).ne.0) detm(L,ip)=detm(L,ip)/detm(0,ip)
		enddo
	if(drpm(6,ip).gt.0) drpm(5,ip)=drpm(5,ip)/drpm(6,ip)
	if(detm(6,ip).gt.0) detm(5,ip)=detm(5,ip)/detm(6,ip)
		enddo
c...
	write(30,801)  0.0d0,(drpm(0,iprec(jp)),jp=1,npt) ! dN/dY
	write(30,801)  1.0d0,(drpm(1,iprec(jp)),jp=1,npt) ! <Pt>
	write(30,801)  2.0d0,(drpm(3,iprec(jp)),jp=1,npt) ! <Mt>
	write(30,801)  3.0d0,(drpm(4,iprec(jp)),jp=1,npt) ! v2
	write(30,801) 10.0d0,(detm(0,iprec(jp)),jp=1,npt) ! dN/deta
	write(30,801) 11.0d0,(detm(1,iprec(jp)),jp=1,npt) ! <Pt>
	write(30,801) 12.0d0,(detm(3,iprec(jp)),jp=1,npt) ! <Mt>
	write(30,801) 13.0d0,(detm(4,iprec(jp)),jp=1,npt) ! v2
c
	write(41,801)  srtnn,(drpm(0,iprec(jp)),jp=1,npt) ! dN/dY
	write(41,801)  srtnn,(drpm(1,iprec(jp)),jp=1,npt) ! <Pt>
	write(41,801)  srtnn,(drpm(3,iprec(jp)),jp=1,npt) ! <Mt>
	write(41,801)  srtnn,(drpm(4,iprec(jp)),jp=1,npt) ! v2
	write(41,801)  srtnn,(detm(0,iprec(jp)),jp=1,npt) ! dN/deta
	write(41,801)  srtnn,(detm(1,iprec(jp)),jp=1,npt) ! <Pt>
	write(41,801)  srtnn,(detm(3,iprec(jp)),jp=1,npt) ! <Mt>
	write(41,801)  srtnn,(detm(4,iprec(jp)),jp=1,npt) ! v2
c...] Mid-Rapidity Integrated Observables (dN/dY, <Pt>, v2, ...)
c...] Midrapidity -----------------------------------------------------*

c...[ Rapidity, PseudoRapdity Dist. -----------------------------------*
	call transDNDY(drp,mxk,mxy,mpt)
	call transDNDY(det,mxk,mxy,mpt)
c...[ Rapidity
	nxy=min(nint(ypcm*2/dy),mxy)
	do i=-nxy,nxy
c...[ particle
	  do jp=1,npt
	ip=iprec(jp)
c...Normalization
	    do L=0,mxk
	drp(L,i,ip)=drp(L,i,ip)/dy/nevent0
	det(L,i,ip)=det(L,i,ip)/dy/nevent0
	    enddo
	    do L=1,mxk
	if(drp(0,i,ip).ne.0) drp(L,i,ip)=drp(L,i,ip)/drp(0,i,ip)
	if(det(0,i,ip).ne.0) det(L,i,ip)=det(L,i,ip)/det(0,i,ip)
	    enddo
	if(drp(6,i,ip).gt.0) drp(5,i,ip)=drp(5,i,ip)/drp(6,i,ip)
	if(det(6,i,ip).gt.0) det(5,i,ip)=det(5,i,ip)/det(6,i,ip)
	  enddo
c...] particle
	enddo
c...] Rapidity
	if(isym.eq.1) then
	  do i=0,mxy
	    do jp=1,npt
	      ip=iprec(jp)
	      do L=0,6
		isgn=1
		if(L.eq.2) isgn=-1
	        drp(L, i,ip)=(drp(L,i,ip)+isgn*drp(L,-i,ip))/2
	        drp(L,-i,ip)=isgn*drp(L,i,ip)
	        det(L, i,ip)=(det(L,i,ip)+isgn*det(L,-i,ip))/2
	        det(L,-i,ip)=isgn*det(L,i,ip)
	      enddo
	    enddo
	  enddo
	endif

	do i=-mxy,mxy
	yrap=dy*i/ypcm
	  write(22,801) yrap,(drp(0,i,iprec(jp)),jp=1,npt) !dN/dY
	  write(23,801) yrap,(drp(1,i,iprec(jp)),jp=1,npt) !<pav>
	  write(24,801) yrap,(drp(2,i,iprec(jp)),jp=1,npt) !<px>
	  write(25,801) yrap,(drp(3,i,iprec(jp)),jp=1,npt) !v2(Y)
	  write(26,801) yrap,(drp(5,i,iprec(jp)),jp=1,npt) !<btr>
	  write(27,801) yrap,(det(0,i,iprec(jp)),jp=1,npt) !dN/deta
	  write(28,801) yrap,(det(3,i,iprec(jp)),jp=1,npt) !v2(eta)
	enddo
c...] Rapidity, PseudoRapdity Dist. -----------------------------------*

c...[ Radial Flow
	dymid1=dymid
	if(dymid.eq.0) dymid1=0.5d0
		do i=1,100
	if(betr(1,i).gt.0) betr(2,i)=betr(2,i)/betr(1,i)
	if(betr(3,i).gt.0) betr(4,i)=betr(4,i)/betr(3,i)
	if(betr(5,i).gt.0) betr(6,i)=betr(6,i)/betr(5,i)
	betr(1,i)=betr(1,i)/dymid1/2/nevent0
	betr(3,i)=betr(3,i)/dymid1/2/nevent0
	betr(5,i)=betr(5,i)/dymid1/2/nevent0
	write(29,801) (i-0.5d0)*dlr,(betr(k,i),k=1,6)
		enddo
c...] Radial Flow

801	format(21(1x,1pe10.3))
	if(inext.eq.1) goto 2000
c...] New Incident Energy ---------------------------------------------*
	end
************************************************************************
	subroutine annini(ic)
************************************************************************
	implicit real*8(a-h,o-z)
	parameter(mxy=100,mxk=6)
	parameter(mpt=100)
	common /par2/r0,r1,r2,r3,ipk
	common /par3/v2ev,nc,ncm
	common /dim0/ObsY(0:mxk)
	common /dim1/dmtY(100,mpt),dptY(100,mpt),elpY(100,mpt)
	common /dimY/dmtE(100,mpt),dptE(100,mpt),elpE(100,mpt)
	common /dim2/drp(0:mxk,-mxy:mxy,mpt)
	common /dim3/det(0:mxk,-mxy:mxy,mpt)
	common /dim4/betr(6,100)
	common /dim5/drpm(0:mxk,mpt),detm(0:mxk,mpt)
c...
	common /par1/pi,srtnn,dy,dymid,ypcm,deta,dm,bmin,bmax,rtmax,dlr
	common /parM1/dy0,dymid0,dm0,bmin0,bmax0,rtmax0
	common /parM2/nevent0,nevent,nev,ifileEv
	common /parM3/isym,inext,iyrat
	parameter(nf1=16,nf2=31)
c...[ For IO
	common /IOch1/FnameEv(20)
	common /IOevp/bminEv(20),bmaxEv(20),WeiEv(20)
     &		   ,nevEv(20),ipkEv(20),nfileEv,iprec(mpt),npt
	character FnameEv*80
c...] For IO
c ---------------------------------------------------------------------*
c...Local Variables
	character head*26
	character fname*50,fnamei*80
	save head
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
	write(*,*) 'rtmax =?'
	read(*,*) rtmax0
c
	return
	
c...[ New Incident Energy (ic=2) --------------------------------------*
2000	continue
	if(inext.eq.1) then
	  do ifile=nf1,nf2
	    close(ifile)
	  enddo
	endif
c...Incident Energy Inputs
	write(*,*) 'srt_nn =? (Negative: Inc. E/A)'
	read(*,*) srtnn
	write(*,*) 'Header of Output File Name=?'
	read(*,'(a)') head
c
	open(42,file=head(1:lengch(head))//'.nbev',status='unknown')
c
	dy=dy0
	dm=dm0
	dymid=dymid0
	bmin=bmin0
	bmax=bmax0
	rtmax=rtmax0
c
		do i=1,100
		do k=1,6
	betr(k,i)=0.0d0
		enddo
		enddo
c
		do k=1,mpt
c
		do i=1,100
	dmtY(i,k)=0.0d0
	dptY(i,k)=0.0d0
	elpY(i,k)=0.0d0
	dmtE(i,k)=0.0d0
	dptE(i,k)=0.0d0
	elpE(i,k)=0.0d0
		enddo
c
		do l=0,mxk
	drpm(l,k)=0.0d0
	detm(l,k)=0.0d0
		do i=-mxy,mxy
	drp(l,i,k)=0.0d0
	det(l,i,k)=0.0d0
		enddo
		enddo
c
		enddo
c...[
	emn=0.938d0
	   if(srtnn.lt.0) then
	einc=-srtnn
	srtnn=sqrt(2*emn*(einc+2*emn))
	   endif
        srtn1=srtnn/2
        pn1=sqrt(srtn1**2-emn**2)
        ypcm=log((srtn1+pn1)/(srtn1-pn1))/2
c
	deta=dymid
	if(dymid.lt.0) dymid=-dymid
	if(dy.lt.0) then
		iyrat=1
		dy=ypcm*abs(dy)
	endif
	if(iyrat.eq.1) then
		dymid=dymid*ypcm
	endif
c...]
c...
c...[ Event File ------------------------------------------------------*
	   do ifileEv=1,20
	write(*,*)
     &	'Event File, Number of Events, PKS format(0-3), bmin, bmax =?'
	read(*,*,err=999,end=1100)
     &	 fnamei,nevEv(ifileEv),ipkEv(ifileEv)
	if(fnamei(1:1).eq.'0'.or.fnamei(1:3).eq.'end') then
		inext=0
		goto 1100
	endif
	if(fnamei(1:4).eq.'next') then
		inext=1
		goto 1100
	endif
	FnameEv(ifileEv)=fnamei
	nfileEv=ifileEv
	   enddo
1100	return
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
	do ifile=nf1,nf2
	  irectyp=1
	  if(ifile.eq.30) irectyp=0
	  if(ifile.eq.31) irectyp=2
	  call recini(ifile,head,nevent0,irectyp)
	enddo
c...
	return

c...]]]]] ic=1,2,3 (goto 1000,2000,3000) ------------------------------*

999	continue
	write(*,*) 'Abnormal End'
	stop

	end
************************************************************************
	subroutine recini(ifile,head,nevent0,irectyp)
************************************************************************
	implicit real*8(a-h,o-z)
	parameter(mpt=100)
	common /par1/pi,srtnn,dy,dymid,ypcm,deta,dm,bmin,bmax,rtmax,dlr
c...[ For IO
	common /IOch1/FnameEv(20)
	common /IOevp/bminEv(20),bmaxEv(20),WeiEv(20)
     &		   ,nevEv(20),ipkEv(20),nfileEv,iprec(mpt),npt
	character FnameEv*80
c...] For IO
	character*(*) head
c ---------------------------------------------------------------------*
c...Local Variables
	character fname*50,fnamei*80
	parameter(nf1=16,nf2=31)
	dimension FileExt(nf1:nf2),Xname(nf1:nf2),Yname(nf1:nf2)
	character FileExt*3,Xname*3,Yname*10
	dimension chap(mpt)
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
	data (FileExt(k),Xname(k),Yname(k),k=nf1,nf2)/
     &		 'mt' ,'dmT','dN/dmT^2dY'	! 16
     &		,'pt' ,'pT' ,'dN/dpT^2dY'	! 17
     &		,'elp','pT' ,'Elip(pT)'		! 18
     &		,'mte','dmT','dN/dmT^2dE'	! 19
     &		,'pte','pT' ,'dN/dpT^2dE'	! 20
     &		,'epe','pT' ,'Elip(pT)'		! 21
     &		,'rap','Y'  ,'dN/dY'		! 22
     &		,'pav','Y'  ,'<p_T>(Y)'		! 23
     &		,'dir','Y'  ,'<p_x>(Y)'		! 24
     &		,'ely','Y'  ,'Elip(Y)'		! 25
     &		,'bet','Y'  ,'beta_t(Y)'	! 26
     &		,'eta','eta','dN/deta'		! 27
     &		,'ele','eta','Elip(eta)'	! 28
     &		,'btr','rt' ,'beta'		! 29
     &		,'mid','Obs','Ave.'		! 30
     &		,'v2c','pT' ,'v2c(pT)'		! 31
     &	/
c ---------------------------------------------------------------------*
	fname=head(1:lengch(head))//'.'//FileExt(ifile)
	open(ifile,file=fname,status='unknown')
	write(ifile,'(a,a,a,a)') '# ',Xname(ifile),'-',Yname(ifile)
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
	write(ifile,'(a,f9.3)')  '# rtmax            =',rtmax
	write(ifile,'(a,f9.3)')  '# dlr              =',dlr
	write(ifile,'(a/a)') '#','# Event Files'
	   do jfile=1,nfileEv
	fnamei=FnameEv(jfile)
	write(ifile,'("# ",i2,": ",a)') jfile,fnamei(1:lengch(fnamei))
	   enddo
	write(ifile,'(a)') '#'

	if(irectyp.eq.1) then
	  write(ifile,802) 1,Xname(ifile),(i+1,chap(iprec(i)),i=1,npt)
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
	implicit real*8(a-h,o-z)
	parameter(mxy=100,mxk=6)
	parameter(mpt=100)
	common /par1/pi,srtnn,dy,dymid,ypcm,deta,dm,bmin,bmax,rtmax,dlr
	common /par2/r0,r1,r2,r3,ipk
	common /par3/v2ev,nc,ncm
	common /dim0/ObsY(0:mxk)
	common /dim1/dmtY(100,mpt),dptY(100,mpt),elpY(100,mpt)
	common /dimY/dmtE(100,mpt),dptE(100,mpt),elpE(100,mpt)
	common /dim2/drp(0:mxk,-mxy:mxy,mpt)
	common /dim3/det(0:mxk,-mxy:mxy,mpt)
	common /dim4/betr(6,100)
	common /dim5/drpm(0:mxk,mpt),detm(0:mxk,mpt)
c...
	common /parM1/dy0,dymid0,dm0,bmin0,bmax0,rtmax0
	common /parM2/nevent0,nevent,nev,ifileEv
	common /parM3/isym,inext,iyrat
c...[ For IO
	common /IOch1/FnameEv(20)
	common /IOevp/bminEv(20),bmaxEv(20),WeiEv(20)
     &		   ,nevEv(20),ipkEv(20),nfileEv,iprec(mpt),npt
	character FnameEv*80
c...] For IO
c ---------------------------------------------------------------------*
c...
	data dlr,  dy0,  dymid0,dm0,  bmin0,bmax0, rtmax0  /
     &	     0.5d0,0.5d0,1.0d0, 0.1d0,0.0d0,20.0d0,-100.0d0/
	data npt,ipk,isym,inext,iyrat/
     &	     18 ,3  ,0,   0    ,0    /
c...
      data iprec/51,52,53,55,1,2,3,4,7,12,22,56,57,58,59,60,61,63,72*0/
c    &	 'Total','h+-','h-','h+','pi-','pi0','pi+','K-','K+','p','pbar'
c    &	,'pi+-','K+-','net p','p/pbar','N','B','net B'
c
	end
************************************************************************
	subroutine headline(jev,ntot,nbar,nmes,bimp,ipk,icon)
************************************************************************
	implicit real*8(a-h,o-z)
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
	subroutine readptcl(id,em,p1,p2,p3,r1,r2,r3,r0,ipk,icon)
************************************************************************
	implicit real*8(a-h,o-z)
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
	return

999	continue
	write(*,*) 'Read Error for particles'
	stop
	end
************************************************************************
	subroutine transDMT(dmt,n1,n2,mpt)
************************************************************************
	implicit real*8(a-h,o-z)
	dimension dmt(n1:n2,mpt)
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
	dmt(i,31)=dmt(i,31)+dmt(i, 1)+dmt(i, 4)                    ! m-
	dmt(i,32)=dmt(i,32)+dmt(i, 2)+dmt(i, 5)+dmt(i, 6)          ! m0
	dmt(i,33)=dmt(i,33)+dmt(i, 3)+dmt(i, 7)                    ! m+
	dmt(i,34)=dmt(i,34)          +dmt(i,19)+dmt(i,14)+dmt(i,17)! b-
	dmt(i,35)=dmt(i,35)+dmt(i,11)+dmt(i,13)+dmt(i,15)+dmt(i,18)! b0
	dmt(i,36)=dmt(i,36)+dmt(i,12)          +dmt(i,16)          ! b+
	dmt(i,37)=dmt(i,37)+dmt(i,22)          +dmt(i,26)          ! a-
	dmt(i,38)=dmt(i,38)+dmt(i,21)+dmt(i,23)+dmt(i,25)+dmt(i,28)! a0
	dmt(i,39)=dmt(i,39)          +dmt(i,29)+dmt(i,24)+dmt(i,27)! a+
c
	dmt(i,53)=dmt(i,31)+dmt(i,34)+dmt(i,37)	! h-
	dmt(i,54)=dmt(i,32)+dmt(i,35)+dmt(i,38)	! h0
	dmt(i,55)=dmt(i,33)+dmt(i,36)+dmt(i,39)	! h+
	dmt(i,52)=dmt(i,53)+dmt(i,55)		! h+-
	dmt(i,51)=dmt(i,52)+dmt(i,54)		! hadrons
     &           +dmt(i,41)+dmt(i,42)+dmt(i,43)	! lep/gam
c
	dmt(i,56)=dmt(i,1)+dmt(i,3)		! pi+-
	dmt(i,57)=dmt(i,4)+dmt(i,7)		! K+-
	dmt(i,58)=dmt(i,12)-dmt(i,22)		! net p
	dmt(i,59)=dmt(i,12)+dmt(i,22)		! p/pbar
	dmt(i,60)=dmt(i,11)+dmt(i,12)		! N
c
	dmt(i,61)=dmt(i,34)+dmt(i,35)+dmt(i,36)	! B
	dmt(i,62)=dmt(i,37)+dmt(i,38)+dmt(i,39)	! aB
	dmt(i,63)=dmt(i,61)-dmt(i,62)		! net B
c
	enddo
	end
************************************************************************
	subroutine transDNDY(drp,mxk,mxy,mpt)
************************************************************************
	implicit real*8(a-h,o-z)
	dimension drp(0:mxk,-mxy:mxy,mpt)
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
	drp(L,i,31)=drp(L,i,31)+drp(L,i, 1)+drp(L,i, 4)            ! m-
	drp(L,i,32)=drp(L,i,32)+drp(L,i, 2)+drp(L,i, 5)+drp(L,i, 6)! m0
	drp(L,i,33)=drp(L,i,33)+drp(L,i, 3)+drp(L,i, 7)            ! m+
	drp(L,i,34)=drp(L,i,34)+drp(L,i,14)+drp(L,i,17)+drp(L,i,19)	! b-
	drp(L,i,35)=drp(L,i,35)+drp(L,i,15)+drp(L,i,18)
     &             +drp(L,i,11)+drp(L,i,13)				! b0
	drp(L,i,36)=drp(L,i,36)+drp(L,i,16)
     &             +drp(L,i,12)						! b+
	drp(L,i,37)=drp(L,i,37)+drp(L,i,26)				! aB-
     &             +drp(L,i,22)
	drp(L,i,38)=drp(L,i,38)+drp(L,i,25)+drp(L,i,28)
     &             +drp(L,i,21)+drp(L,i,23)				! aB0
	drp(L,i,39)=drp(L,i,39)+drp(L,i,24)+drp(L,i,27)+drp(L,i,29)	! aB+
c
	drp(L,i,53)=drp(L,i,31)+drp(L,i,34)+drp(L,i,37)	! h-
	drp(L,i,54)=drp(L,i,32)+drp(L,i,35)+drp(L,i,38)	! h0
	drp(L,i,55)=drp(L,i,33)+drp(L,i,36)+drp(L,i,39)	! h+
	drp(L,i,52)=drp(L,i,53)+drp(L,i,55)		! h+-
	drp(L,i,51)=drp(L,i,52)+drp(L,i,54)		! hadrons
     &           +drp(L,i,41)+drp(L,i,42)+drp(L,i,43)	! lep/gam
c
	drp(L,i,56)=drp(L,i,1)+drp(L,i,3)		! pi+-
	drp(L,i,57)=drp(L,i,4)+drp(L,i,7)		! K+-
	drp(L,i,58)=drp(L,i,12)-drp(L,i,22)		! net p
	drp(L,i,59)=drp(L,i,12)+drp(L,i,22)		! p/pbar
	drp(L,i,60)=drp(L,i,11)+drp(L,i,12)		! N
c
	drp(L,i,61)=drp(L,i,34)+drp(L,i,35)+drp(L,i,36)	! B
	drp(L,i,62)=drp(L,i,37)+drp(L,i,38)+drp(L,i,39)	! aB
	drp(L,i,63)=drp(L,i,61)-drp(L,i,62)		! net B
c
	enddo
	enddo
c...] Translation -----------------------------------------------------*
	end
************************************************************************
	subroutine addbin(id,em,p1,p2,p3,rap,eta,iQ)
************************************************************************
	implicit real*8(a-h,o-z)
	parameter(mxy=100,mxk=6)
	parameter(mpt=100)
	common /par1/pi,srtnn,dy,dymid,ypcm,deta,dm,bmin,bmax,rtmax,dlr
	common /par2/r0,r1,r2,r3,ipk
	common /par3/v2ev,nc,ncm
	common /dim0/ObsY(0:mxk)
	common /dim1/dmtY(100,mpt),dptY(100,mpt),elpY(100,mpt)
	common /dimY/dmtE(100,mpt),dptE(100,mpt),elpE(100,mpt)
	common /dim2/drp(0:mxk,-mxy:mxy,mpt)
	common /dim3/det(0:mxk,-mxy:mxy,mpt)
	common /dim4/betr(6,100)
	common /dim5/drpm(0:mxk,mpt),detm(0:mxk,mpt)
	common /dim6/vcsum(0:2,0:mpt),sumNCS(0:2,0:mpt)
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
	rr=sqrt(r1**2+r2**2)
	ObsY(0)=1.0d0			! Rapidity Dist.
	ObsY(1)=pt			! <Pt>
	ObsY(2)=p1			! <Px>: Directed Flow
	ObsY(3)=(p1**2-p2**2)/pts	! <(Px**2-Py**2)/Pt**2>: v2
	ObsY(4)=emt
	ObsY(5)=0.0d0
	ObsY(6)=0.0d0
c...[ Event Analysis
c	Iadd=1
	if(iQ.ne.0) nc=nc+1
c...] Event Analysis
c	if(ipk.eq.2.and.rr.gt.0.and.rr.lt.10.0d0) then
	if((ipk.eq.2.or.ipk.eq.3).and.rr.gt.0) then
		btr=(p1*r1+p2*r2)/rr/emt !rad.flow
	   if(rtmax.lt.0.or.rr.lt.rtmax) then
		ObsY(5)=btr
		ObsY(6)=1.0d0
	   endif
	   if(abs(rap).le.dymid.or.dymid.eq.0) then
c     		write(98,800) id,em,rr,ObsY(5)
c800	format(1x,i9,8(1x,f10.3))
		ir=int(rr/dlr)+1
		   if(ir.lt.100) then
		betr(1,ir)=betr(1,ir)+1.0d0
		betr(2,ir)=betr(2,ir)+btr
			if(abs(id).eq.211.or.abs(id).eq.111) then !pion
		betr(3,ir)=betr(3,ir)+1.0d0
		betr(4,ir)=betr(4,ir)+btr
			elseif(abs(id).eq.2212.or.abs(id).eq.2112) then	! N
		betr(5,ir)=betr(5,ir)+1.0d0
		betr(6,ir)=betr(6,ir)+btr
			endif
	           endif
	   endif
	endif
c...]

c...[ Rapidity/Pseudo Rapidity Distribution ---------------------------*
	if(abs(irap).gt.mxy.and.abs(ieta).gt.mxy) goto 3000
	if(abs(irap).le.mxy) then
	   do L=0,mxk
	      drp(L, irap,ip)=drp(L, irap,ip)+ObsY(L)
	   enddo
	endif
	if(abs(ieta).le.mxy) then
	   do L=0,mxk
	      det(L, ieta,ip)=det(L, ieta,ip)+ObsY(L)
	   enddo
	endif
c...] Rapidity/Pseudo Rapidity Distribution ---------------------------*

3000	continue
c...[ MidRapidity Observables -----------------------------------------*

c...[ Pt or Mt Distribution -------------------------------------------*
	imt=int((emt-em)/dm)+1
	ipt=int(pt/dm)+1
	delp=(p1**2-p2**2)/(p1**2+p2**2)
	if(dymid.eq.0.or.abs(rap).lt.dymid) then
	  if(imt.le.100) dmtY(imt,ip)=dmtY(imt,ip)+1.0d0/emt
	  if(ipt.le.100) dptY(ipt,ip)=dptY(ipt,ip)+1.0d0/pt
	  if(ipt.le.100) elpY(ipt,ip)=elpY(ipt,ip)+delp/pt
	endif
	if(dymid.eq.0.or.abs(eta).lt.dymid) then
	  if(imt.le.100) dmtE(imt,ip)=dmtE(imt,ip)+1.0d0/emt
	  if(ipt.le.100) dptE(ipt,ip)=dptE(ipt,ip)+1.0d0/pt
	  if(ipt.le.100) elpE(ipt,ip)=elpE(ipt,ip)+delp/pt
	  if(ipt.le.100) then
	    phi=atan2(p2,p1)
	    sumNCS(0,ipt)=sumNCS(0,ipt)+1.0d0
	    sumNCS(1,ipt)=sumNCS(1,ipt)+cos(2*phi)
	    sumNCS(2,ipt)=sumNCS(2,ipt)+sin(2*phi)
	  endif
	endif
c...] Pt or Mt Distribution -------------------------------------------*

c...[ MidRapidity Selection
	if(dymid.ne.0) then
	if(deta.lt.0) then
	   if(abs(eta).gt.dymid) return
	else
	   if(abs(rap).gt.dymid) return
	endif
	endif
c...] MidRapidity Selection

c	Imid=1
c...[ Event V2 Analysis
	if(iQ.ne.0) then
	    ncm=ncm+1
	    v2ev=v2ev+ObsY(3)
	endif
c...] Event V2 Analysis
c...[ MidRapidity Observables
	do L=0,mxk
	   drpm(L, ip)=drpm(L, ip)+ObsY(L)
	   detm(L, ip)=detm(L, ip)+ObsY(L)
	enddo
c...] MidRapidity Observables

c...] MidRapidity Observables -----------------------------------------*

	end
************************************************************************
	subroutine partype(id,iptyp,iQ)
************************************************************************
	implicit real*8(a-h,o-z)
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
	implicit real*8(a-h,o-z)
	parameter(mpt=100)
	common /par1/pi,srtnn,dy,dymid,ypcm,deta,dm,bmin,bmax,rtmax,dlr
	common /parM2/nevent0,nevent,nev,ifileEv
	common /dim6/vcsum(0:2,0:mpt),sumNCS(0:2,0:mpt)

	goto (1000,2000,3000,4000,5000) ic
	goto 3000
c ---------------------------------------------------------------------*

 1000	continue
	do k=0,2
	  do ipt=0,mpt
	    vcsum(k,ipt)=0.0d0
	  enddo
	enddo
	return

 2000	continue
	do k=0,2
	  do ipt=0,mpt
	    sumNCS(k,ipt)=0.0d0
	  enddo
	enddo
	return

 4000	continue
	do k=0,2
	    sumNCS(k,0)=0.0d0
	  do ipt=1,mpt
	    sumNCS(k,0)=sumNCS(k,0)+sumNCS(k,ipt)
	  enddo
	enddo
	do ipt=0,mpt
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
	return

 5000	continue
	pair0=vcsum(0,0)
	v20=0.0d0
	if(pair0.gt.0.and.vcsum(1,0).gt.0) then
	  v20=sqrt(vcsum(1,0)/pair0)
	endif
	write(31,'("#  pT         v2(corr)   D-v2(corr) <Pair>")')
	do ipt=0,mpt
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
