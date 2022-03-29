c...use Cornelius by C. Nonaka  (modified by Y. Nara)
c*************************************************
      subroutine c_freezeout(i,it,htime) 
c     i =1 : initial, i=2 : calculation 
      implicit none
      include 'fluid.inc'
      integer i,it 
      double precision htime

c     write(*,*) 'now c_freezeout'

      if (i .eq. 1) then
      fe(0,:,:,:)  = el(:,:,:)
      fp(0,:,:,:)  = pl(:,:,:)
      fd(0,:,:,:)  = bl(:,:,:)
      fvx(0,:,:,:) = vx(:,:,:)
      fvy(0,:,:,:) = vy(:,:,:)
      fvz(0,:,:,:) = vz(:,:,:)
      nfreezeout=0
      ifreeze=0
c     write(*,*) 'i=',i,fvx(0,1,1,1),fvx(0,50,50,50)
     
      else if (i .eq. 2) then

        if(mod(it,Ftime).eq.0) then
          fe(1,:,:,:)  = el(:,:,:)
          fp(1,:,:,:)  = pl(:,:,:)
          fd(1,:,:,:)  = bl(:,:,:)
          fvx(1,:,:,:) = vx(:,:,:)
          fvy(1,:,:,:) = vy(:,:,:)
          fvz(1,:,:,:) = vz(:,:,:)

c         write(*,*) 'before freezeout'
c         write(*,*) 'i=',i,fvx(0,1,1,1),fvx(0,50,50,50)
c         write(*,*) 'i=',i,fvx(1,1,1,1),fvx(1,50,50,50)

          call freezeout(htime)

          ! for next freezeout
          fe(0,:,:,:)  =  fe(1,:,:,:)
          fp(0,:,:,:)  =  fp(1,:,:,:)
          fd(0,:,:,:)  =  fd(1,:,:,:)
          fvx(0,:,:,:) = fvx(1,:,:,:)
          fvy(0,:,:,:) = fvy(1,:,:,:)
          fvz(0,:,:,:) = fvz(1,:,:,:)
        end if
       end if  

       end 

c--------------------------------------------------------
      subroutine freezeout(htime) 
      implicit none 
      include 'fluid.inc'
      double precision htime,dgt,dgx,dgy,dgz,rt,rx,ry,rz
      double precision HyperCube(0:1,0:1,0:1,0:1)
      double precision psite(0:1,0:1,0:1,0:1),dsite(0:1,0:1,0:1,0:1) 
      double precision vxsite(0:1,0:1,0:1,0:1),vysite(0:1,0:1,0:1,0:1) 
      double precision vzsite(0:1,0:1,0:1,0:1)
      double precision dSigma(0:3,8),Vmid(0:3,8)
      double precision tmid,xmid,ymid,zmid,pmid,dmid,emid
      double precision vxmid,vymid,vzmid
      double precision Tempmid,mumid,smumid,smid
      double precision getpden,frmultlookup
      real*8 gam,vol,sig,vsig
      integer Nsurf,Nambi,Ndisc
      integer Ngp,n
      integer iz,iy,ix,l,k,j,i,m

c     write(*,*) 'now freezeout'
    
      dgt=dt*Ftime 
      dgx=dx*FX
      dgy=dy*FY
      dgz=dz*FZ

c     write(*,*) 'dgt',dgt

      n=2
      do iz=n,maxz-n-FZ
        do iy = n,maxy-n-FY
          do ix = n, maxx-n-FX

c          if(ifreeze(ix,iy,iz).eq.1) cycle

           Ngp = 0

           if ((mod(ix,FX) .eq. 0) .and. (mod(iy,FY) .eq. 0) .and. 
     $          (mod(iz,FZ) .eq. 0) ) then 
              do l = 0, 1 ! z
                 do k = 0, 1 ! y
                    do j = 0, 1 ! x
                      do i = 0, 1 ! t
                      if( fe(i,ix+j*FX,iy+k*FY,iz+l*FZ) .ge.efreeze) 
     $                Ngp=Ngp+1
                      end do
                    end do
                 end do
              end do


              if ((Ngp .gt. 0) .and. (Ngp .lt. 16)) then

                do l = 0, 1
                 do k = 0, 1
                  do j = 0, 1
                   do i = 0, 1
                     HyperCube(i,j,k,l) = fe(i,ix+j*FX,iy+k*FY,iz+l*FZ)
                     psite(i,j,k,l)  = fp(i,ix+j*FX,iy+k*FY,iz+l*FZ)
                     dsite(i,j,k,l)  = fd(i,ix+j*FX,iy+k*FY,iz+l*FZ)
                     vxsite(i,j,k,l) = fvx(i,ix+j*FX,iy+k*FY,iz+l*FZ)
                     vysite(i,j,k,l) = fvy(i,ix+j*FX,iy+k*FY,iz+l*FZ)
                     vzsite(i,j,k,l) = fvz(i,ix+j*FX,iy+k*FY,iz+l*FZ)
                    end do
                   end do
                  end do
                 end do

             ! for version 1.2 
c            call Cornelius(efreeze,HyperCube,Ngp,dSigma,Nsurf,Vmid,
c    $                      dgt,dgx,dgy,dgz,Nambi,Ndisc)

c            for version 1.4 
             call Cornelius(efreeze,HyperCube,dSigma,Nsurf,Vmid,       
     $                      dgt,dgx,dgy,dgz,Nambi,Ndisc)
             do k=1,Nsurf
               rt = Vmid(0,k)/dgt
               rx = Vmid(1,k)/dgx
               ry = Vmid(2,k)/dgy 
               rz = Vmid(3,k)/dgz

c              write(*,*) 'Nsurf=',Nsurf,k,rt,rx,ry,rz
               tmid = htime -Ftime*dt + rt*dgt
               xmid = dx*(ix-origx) 
               ymid = dy*(iy-origy) 
               zmid = dz*(iz-origz) 

c              call interpolation(HyperCube, emid, rt, rx, ry, rz)
c              call interpolation(psite, pmid, rt, rx, ry, rz)
c              call interpolation(dsite, dmid, rt, rx, ry, rz)
c              call interpolation(vxsite, vxmid, rt, rx, ry, rz)
c              call interpolation(vysite, vymid, rt, rx, ry, rz)
c              call interpolation(vzsite, vzmid, rt, rx, ry, rz)

               call interpolation4(HyperCube, emid, rt, rx, ry, rz)
               call interpolation4(psite, pmid, rt, rx, ry, rz)
               call interpolation4(dsite, dmid, rt, rx, ry, rz)
               call interpolation4(vxsite, vxmid, rt, rx, ry, rz)
               call interpolation4(vysite, vymid, rt, rx, ry, rz)
               call interpolation4(vzsite, vzmid, rt, rx, ry, rz)

c              gam = 1.0/sqrt(1.0-vxmid*vxmid-vymid*vymid-vzmid*vzmid)
c              vol=gam*(dSigma(0,k)+dSigma(1,k)*vxmid+dSigma(2,k)*vymid
c    &         + dSigma(3,k)*vzmid)
c              sig=sqrt(vol**2 -
c    & (dSigma(0,k)**2 - dSigma(1,k)**2-dSigma(2,k)**2-dSigma(3,k)**2))

c             if(vol.le.0.0d0.and.abs(vol/sig).gt.1.0d0) cycle
  
               call geteos2(emid,dmid,Tempmid,mumid,smumid,smid)
   
c              write(65,965) tmid,xmid,ymid,zmid,
c    $                       dSigma(0,k),dSigma(1,k),dSigma(2,k),
c    $                       dSigma(3,k),
c    $                       vxmid,vymid,vzmid,
c    $                       pmid,emid,Tempmid,mumid,smumid

c         if((emid.gt.efreeze*1.2) .or.(Tempmid.lt.tfreezecut)) then
c            print *,'Tf=',Tempmid*hbc,tfreezecut*hbc,'e=',emid*hbc,
c    &        'b=',dmid
c         endif
             
c               if(emid.gt.efreeze*1.2) cycle
c               if(emid.lt.efreeze*0.8) cycle
c               if(Tempmid.lt.tfreezecut) cycle
c               if(Tempmid.lt.tfreezecut) then
c                print *,'Tf=',Tempmid*hbc,tfreezecut*hbc,'e=',emid*hbc,
c    &        'b=',dmid
c                 Tempmid=tfreezecut
c               endif


                nfreezeout=nfreezeout+1
                frt(nfreezeout)=tmid
                frr(1,nfreezeout)=ix
                frr(2,nfreezeout)=iy
                frr(3,nfreezeout)=iz
                frv(1,nfreezeout)=vxmid
                frv(2,nfreezeout)=vymid
                frv(3,nfreezeout)=vzmid
                frtmp(nfreezeout)=max(Tempmid,tfreezecut)
                frmuq(nfreezeout)=mumid
                frmus(nfreezeout)=smumid
                frdsgm(0,nfreezeout)=dSigma(0,k)
                frdsgm(1,nfreezeout)=dSigma(1,k)
                frdsgm(2,nfreezeout)=dSigma(2,k)
                frdsgm(3,nfreezeout)=dSigma(3,k)

                frbdn(nfreezeout)=dmid
                ifreeze(ix,iy,iz)=1


c               if(abs(vol/sig).le.1.0) then
c               write(65,965)vol,sig
c               endif

c               write(65,965)tmid,emid*hbc,dmid,tempmid*hbc,
c    &                mumid*hbc,smumid

c               if(emid*hbc.ge.0.8) write(66,965)tmid,emid*hbc

                if(opt_table.eq.2) then 
                  frnum(nfreezeout)=getpden(emid,dmid)
                else
                  frnum(nfreezeout)=frmultlookup(emid*hbc,dmid)
                endif


             end do  ! Nsurf

          
               endif
           endif

          end do
        end do
      end do

965   format(16e15.6)

      
      end

c--------------------------------------------------------
      subroutine interpolation(qsite, qmid, rt, rx, ry, rz)
      implicit none
      double precision qsite(0:1,0:1,0:1,0:1)
      double precision qmid, rt, rx, ry, rz
      double precision st, sx, sy, sz

      st = 1.d0 - rt;  sx = 1.d0 - rx;  sy = 1.d0 - ry;  sz = 1.d0 - rz;

c    ! 0 = s, 1 = r

      qmid  = st*sx*sy*sz*qsite(0,0,0,0) 
     $  + st*sx*sy*rz*qsite(0,0,0,1) 
     $  + st*sx*ry*sz*qsite(0,0,1,0) 
     $  + st*sx*ry*rz*qsite(0,0,1,1) 
     $  + st*rx*sy*sz*qsite(0,1,0,0) 
     $  + st*rx*sy*rz*qsite(0,1,0,1) 
     $  + st*rx*ry*sz*qsite(0,1,1,0) 
     $  + st*rx*ry*rz*qsite(0,1,1,1) 
     $  + rt*sx*sy*sz*qsite(1,0,0,0) 
     $  + rt*sx*sy*rz*qsite(1,0,0,1) 
     $  + rt*sx*ry*sz*qsite(1,0,1,0) 
     $  + rt*sx*ry*rz*qsite(1,0,1,1) 
     $  + rt*rx*sy*sz*qsite(1,1,0,0) 
     $  + rt*rx*sy*rz*qsite(1,1,0,1) 
     $  + rt*rx*ry*sz*qsite(1,1,1,0) 
     $  + rt*rx*ry*rz*qsite(1,1,1,1)

      end 

c--------------------------------------------------------
      subroutine interpolation4(qsite, qmid, rt, rx, ry, rz)
      implicit none
      double precision qsite(0:1,0:1,0:1,0:1)
      double precision qmid, rt, rx, ry, rz
      integer j1,j2,j3,j4

      qmid=0d0
      do j4=0,1
      do j3=0,1
      do j2=0,1
      do j1=0,1
        qmid = qmid
     $      + (j1*rt+(1-j1)*(1-rt))*(j2*rx+(1-j2)*(1-rx))
     $       *(j3*ry+(1-j3)*(1-ry))*(j4*rz+(1-j4)*(1-rz))
     $       *qsite(j1,j2,j3,j4)
      end do
      end do
      end do
      end do

      end

c--------------------------------------------------------
      subroutine isochronous_freezeout(htime,opt) 
      implicit none 
      include 'fluid.inc'
      double precision htime,tcut
      double precision e,d,tf,muq,mus,sc
      double precision getpden,frmultlookup
      integer iz,iy,ix,n,opt
      real*8 mult,a,eden,pre,pden,bden,pnow,csnow

      n=2
      do iz=n,maxz-n
      do iy=n,maxy-n
      do ix=n,maxx-n
        if(u(0,ix,iy,iz).le.1d-7) cycle
c       if(tl(ix,iy,iz).le.tfreezecut) cycle ! Skip too low temperature.
        e=el(ix,iy,iz)
        d=bl(ix,iy,iz)

c       if(e.le.0.0001) cycle
        if(opt.eq.0.and.e.gt.efreeze) cycle

        call geteos2(e,d,tf,muq,mus,sc)
        if(tf.le.tfreezecut) then
c         print *,'isochronos_freezeout tf?',tf,tl(ix,iy,iz)
          tf=tfreezecut
        endif

        nfreezeout=nfreezeout+1
        frt(nfreezeout)=htime
        frr(1,nfreezeout)=ix
        frr(2,nfreezeout)=iy
        frr(3,nfreezeout)=iz
        frv(1,nfreezeout)=vx(ix,iy,iz)
        frv(2,nfreezeout)=vy(ix,iy,iz)
        frv(3,nfreezeout)=vz(ix,iy,iz)
        frtmp(nfreezeout)=tf
        frmuq(nfreezeout)=muq
        frmus(nfreezeout)=mus
        frdsgm(0,nfreezeout)=dx*dy*dz
        frdsgm(1,nfreezeout)=0d0
        frdsgm(2,nfreezeout)=0d0
        frdsgm(3,nfreezeout)=0d0

        frbdn(nfreezeout)=d

        if(opt_table.eq.2) then 
          frnum(nfreezeout)=getpden(e,d)
         else
          frnum(nfreezeout)=frmultlookup(e*hbc,d)
         endif

      end do
      end do
      end do

965   format(16e15.6)

      end

c**************************************************************************
      real*8 function frmultlookup(e0,b0)
      implicit none
      integer i,j
      real*8 e,b,x2,y2,x1,y1,e0,b0

      integer nemax,nbmax
      parameter(nemax=50,nbmax=25)
      real*8 frmtable,emin,bmin,de,db
      common/frmulttable/frmtable(0:nbmax,0:nemax),emin,bmin,de,db

      b=b0
      e=e0
      i=int((b-bmin)/db)
      j=int((e-emin)/de)
      if(i.ge.nbmax) then
        print *,'frmultlookup b=',b,i
        b=nbmax*db+bmin
      endif
      if(j.ge.nemax) then
        print *,'frmultlookup e=',e,j
        e=nemax*de+emin
      endif



      i=min(max(0,int((b-bmin)/db)),nbmax-1)
      j=min(max(0,int((e-emin)/de)),nemax-1)

      x2=(b-(bmin+i*db))/db
      y2=(e-(emin+j*de))/de
      x1=1d0-x2
      y1=1d0-y2

      frmultlookup=x1*y1*frmtable(i,j)
     &       +x2*y1*frmtable(i+1,j)
     &       +x1*y2*frmtable(i,j+1)
     &       +x2*y2*frmtable(i+1,j+1)

      if(frmultlookup.lt.0d0) then
        print *,'mult<0?',frmultlookup,i,j
        frmultlookup=0.0d0
      endif

      end

c**************************************************************************
      subroutine read_frezeout_mult_table2(fname)
      implicit none
      character  fname*(*),tmp*90
      real*8 K,B
      integer nemax,nbmax
      parameter(nemax=50,nbmax=25)
      real*8 frmtable,emin,bmin,de,db
      common/frmulttable/frmtable(0:nbmax,0:nemax),emin,bmin,de,db
      integer i,j
      real*8 bden,e,t,mub,bmu,smu,pdens

      open(82,file=fname,status='old')
      emin=0.0d0
      bmin=0.0d0

      de=10.0d0/499
      db=3.0d0/149

      read(82,'(a)')tmp
      read(82,'(a)')tmp
      read(82,'(a)')tmp
      read(82,*) K,B

      do i=0,nbmax-1
      do j=0,nemax-1
        read(82,*) bden,e,t,bmu,smu,pdens
        frmtable(i,j)=pdens
      end do
      end do

      close(82)

      end

c**************************************************************************
      subroutine edit_frezeout_surface
      implicit none
      include 'fluid.inc'
      integer i,i1,j,ix,iy,iz
c...Remove unwanted freeze-out surface.

      do i=1,nfreezeout
        if(frtmp(i).le.tfreezecut) then
          ix=frr(1,i)
          iy=frr(2,i)
          iz=frr(2,i)
          do j=1,nfreezeout
            if(ix.eq.frr(1,j).and.iy.eq.frr(2,j).and.iz.eq.frr(3,j))
     &      frtmp(j)=0.0d0             
          end do
        endif
      enddo

      i1=0
      do 110 i=1,nfreezeout
        if(frtmp(i).le.tfreezecut) goto 110
        i1=i1+1 
        frt(i1)=frt(i)
        frr(:,i1)=frr(:,i)
        frv(:,i1)=frv(:,i)
        frtmp(i1)=frtmp(i)
        frmuq(i1)=frmuq(i)
        frmus(i1)=frmus(i)
        frdsgm(:,i1)=frdsgm(:,i)
        frbdn(i1)=frbdn(i)
        frnum(i1)=frnum(i)
  110 continue 
      nfreezeout=i1

      end

