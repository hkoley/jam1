c*********************************************************************
      subroutine init_fluid

      implicit none
      include 'fluid.inc'
      include 'jam2.inc'

      h_mode=mstc(140)
c     if(h_mode.eq.2) dt=0.4*dx
c     if(h_mode.ge.3) dt=0.7*dx

      if(mstc(50).lt.20) then
        eos_mode=mstc(50)
      else if(mstc(50).ge.21) then
        eos_mode=mod(mstc(50),10)
      else
        print *,'invalid eos mode',eos_mode
        stop
      endif

      second_order=mstc(136)

      ntest=mstc(5)
      atomic=mstd(2)
      ycm=pard(18)
      t_pass=pard(140)
      efreeze=parc(140)/hbc
      dcut=parc(147)

      dx=0.3
      dy=0.3
      dz=0.3

c     dz=0.15

c     dx=0.2
c     dy=0.2
c     dz=0.2

      dt=parc(2)
      pard(143)=dx
      pard(144)=dy
      pard(145)=dz

      mstd(143)=maxx
      mstd(144)=maxy
      mstd(145)=maxz

      job_mode=mstc(8)

c...Flux is modified so that all Us are timelike by K.Murase's method.
      opt_timelike=1

c...Rescale momentum (=2) or energy (=1) to make U timelike.
      opt_rescale=0

c...Gaussian width for smearing particle.
      gsigma=parc(142)

c..Option for different form of Gaussian.
c     opt_gauss = mstc(146)

      opt_freezeout=mstc(144)
      opt_table=mstc(145)

      n_gaussp=mstc(138)

      if(opt_freezeout.ge.1) then
        if(opt_table.eq.0) then 
          call make_mult_table
        else
c         call read_frezeout_mult_table(fname(10))
          if(mstc(145).eq.1) call read_frezeout_mult_table2(fname(10))
        endif
      endif

      end

c***************************************************************************
      double precision function edminfluid(binv)
      implicit none
      include 'jam2.inc'
      real*8 binv,edminQGP

      if(parc(144).gt.0.0d0) then
        edminfluid=parc(144)
        return
      endif

      edminfluid = abs(parc(144)) + edminQGP(binv)

      end

c...This is called from jam.f when mstc(181)=1 or larger.
c***************************************************************************
      subroutine dens0output(kdt,dtp,msel)

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'fluid.inc'
      real*8 uv(4),dtp
      integer nev,msel,kdt
      data nev/1/
      save nev
      dimension cur(4),curb(4),tens(4,4)
      common/jamdattens/cur,curb,tens
      save/jamdattens/
c     parameter(optl=2)
      opt=mstc(181)

      if(kdt.eq.0) then
       write(30,800)nev
       write(31,800)nev
       nev=nev+1
      endif
 800  format('# nev = ',i7)

      if(msel.eq.1) then
      call jamdens(icon,0,cc,rho,rhob,einv,uv,0)

c...Local energy and baryon density.
c     if(optl.eq.1) then

      write(30,820)pard(1),
     &  bl(maxx/2,maxy/2,maxz/2),
     &  el(maxx/2,maxy/2,maxz/2)*hbc,
     &  rhob,einv,1

c     else

      write(31,820)pard(1),
     &  u(4,maxx/2,maxy/2,maxz/2),
     &  u(0,maxx/2,maxy/2,maxz/2)*hbc,
     &  curb(4),tens(4,4),1

c     endif


      else
       pard1=pard(1)
       n=pard1/dt-1
       do i=1,n
        pard(1)=i*dtp
        call jamdens(icon,0,cc,rho,rhob,einv,uv,0)
        write(30,820)pard(1),0.0, 0.0, rhob,einv,0
        write(31,820)pard(1),0.0, 0.0, curb(4),tens(4,4),0
       end do
       pard(1)=pard1
      endif

      flush(30)
      flush(31)
820   format(f6.3,4(1x,f10.4),1x,i1)

c     call outputjamzx
c     call outputjamxy

      end

c***************************************************************************
      subroutine outputjamxy
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'fluid.inc'
      integer ix,iy,iz,it
      double precision x,y,z
      real*8 uv(4)
      dimension cur(4),curb(4),tens(4,4)
      common/jamdattens/cur,curb,tens

      iz = origz
      z0=0.0
      t00=0.0d0
      t00f=0.0d0
      do ix = 0,maxx
      do iy = 0,maxy
        x0 = dx*(ix-origx)
        y0 = dy*(iy-origy)
        call jamdens2(icon,x0,y0,z0,cc,rho,rhob,einv,uv,0)
        t00=t00+tens(4,4)*dx*dy
        t00f=t00f+u(0,ix,iy,iz)*dx*dy*hbc
      end do
      end do

      total=t00+t00f
      rp=0.0
      rf=0.0
      if(total.gt.0d0) then
        rp=t00/total
        rf=t00f/total
      endif
      write(32,800)pard(1),t00,t00f,rp,rf
 800  format(5(f8.3,1x),3x,4(e15.7,6x))

      flush(32)

       end

c***************************************************************************
      subroutine outputjamzx
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'fluid.inc'
      integer ix,iy,iz,it
      double precision x,y,z
      real*8 uv(4)
      dimension cur(4),curb(4),tens(4,4)
      common/jamdattens/cur,curb,tens

      open(unit=11,file='densp.dat',status='unknown')


      write(11,801)pard(1)
 801  format('# time= ',f8.3)

      iy = origy         
      y0=0.0
      do ix = 0,maxx
      do iz = 0,maxz
        x0 = dx*(ix-origx)
        z0 = dz*(iz-origz)
        call jamdens2(icon,x0,y0,z0,cc,rho,rhob,einv,uv,0)
       write(11,800) z0,x0,cur(4),tens(4,4)
      end do
       write(11,*)' '
      end do
 800  format(2(f8.3,1x),3x,3(e15.7,6x))

      close(11)

       end

c...This is called from jam.f when mstc(182)=1 or larger.
c***************************************************************************
      subroutine hydrofraction_output(kdt,dtp,msel)

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'fluid.inc'
      real*8 uv(4),dtp
      integer nev,msel,kdt,opt
      data nev/1/
      save nev
c     parameter(optfra=1)
      parameter(nmaxt=300)
      real*8 efluid_tau(0:nmaxt),epart_tau(0:nmaxt)

      optfra=mstc(182)
      if(kdt.eq.0) then
       write(32,800)nev
       nev=nev+1
      endif
 800  format('# nev = ',i7)

      do i=0,nmaxt
        efluid_tau(i)=0.0
        epart_tau(i)=0.0
      end do

        epart1=0.0
        epart=0.0

        efluid1=0.0
        efluid=0.0

      do i=1,nv
          if(k(1,i).eq.0.or.k(1,i).ge.11) cycle
c         if(k(7,i).eq.-1500.or.k(7,i).eq.-1200) cycle
          if(abs(k(7,i)).eq.1) cycle   ! not yet interact
          if(p(5,i).eq.0d0) cycle

          pe=p(4,i)
          pz=p(3,i)
          if(k(2,i).eq.92) then ! String
            pe= (vq(4,i)+vq(9,i))
            pz= (vq(3,i)+vq(8,i))
          endif
           
c         call jamdens(icon,i,cc,rho,rhob,einv,u,1)
c         ep=einv
          ep=pe
          epart = epart + ep

c         if(pe.eq.pz) cycle
c         y=0.5*log(max(pe+pz,1d-8)/max(pe-pz,1d-8))
c         if(abs(y).le.1.0) epart1 = epart1 + ep

          rz=r(3,i)+(pard(1)-r(4,i))*pz/pe
          if(abs(rz).le.parc(182)) epart1 = epart1 + ep

          tau=r(4,i)**2-rz**2
          if(tau.ge.0d0) then
          tau=sqrt(tau)
          itau=nint(tau/dz)
          eta=0.5*log(max(r(4,i)+rz,1d-8)/max(r(4,i)-rz,1d-8))
          if(itau.ge.0.and.itau.le.nmaxt) then
c           epart_tau(itau)=epart_tau(itau)+pe*cosh(eta)-pz*sinh(eta)
            epart_tau(itau)=epart_tau(itau)+pe
          endif
          endif

      end do

        vol=dx*dy*dz
        volh=vol*hbc
        do iz = 0,maxz
          h=dz*(iz-origz)
          tau=pard(1)**2-h**2
          if(tau.ge.0d0) tau=sqrt(tau)
        do ix = 0,maxx
        do iy = 0,maxy

c         ef=el(ix,iy,iz)*volh
          ef=u(0,ix,iy,iz)*volh

c         if(opt_taueta.ne.0) then
c           h=dz*(iz-origz)
c           ef = volh*(cosh(h)*u(0,ix,iy,iz)+sinh(h)*u(3,ix,iy,iz))
c         endif

c        v3 = vz(ix,iy,iz)
c        rap=0.5*log( (1.0+v3)/(1.0-v3) )
         rap=h

          if(optfra.eq.1) then
            efluid  = efluid + ef
            if(abs(rap).le.parc(182)) efluid1  = efluid1 + ef
          else
            if(el(ix,iy,iz).ge.parc(140)) then
              efluid  = efluid + ef
              if(abs(rap).le.1.0) efluid1  = efluid1 + ef
            else
             epart=epart+ef
             if(abs(rap).le.parc(182)) epart1  = epart1 + ef
            endif
          endif

          if(tau.ge.0d0) then
            itau=nint(tau/dz)
            if(itau.ge.0.and.itau.le.nmaxt) then
c             ef=ef*cosh(h) - volh*sinh(h)*u(3,ix,iy,iz)
              efluid_tau(itau)=efluid_tau(itau)+ef
            endif
          endif

        end do
        end do
        end do


      et=epart+efluid
      et1=epart1+efluid1
      ep=0.0
      ef=0.0
      ep1=0.0
      ef1=0.0
      eft=0.0
      if(et.gt.0d0) then
        ep=epart/et
        ef=efluid/et
        eft=pe_tot/et
      endif
      if(et1.gt.0d0) then
        ep1=epart1/et1
        ef1=efluid1/et1
      endif
      write(32,820)pard(1),epart,efluid,ep,ef,ep1,ef1,eft
      flush(32)

c     i29=0
c     do i=0,nmaxt
c       ep=0.0; ef=0.0
c       et=epart_tau(i)+efluid_tau(i)
c       if(et.gt.0d0) then
c         ep=epart_tau(i)/et
c         ef=efluid_tau(i)/et
c         write(29,820)i*dz,ep,ef,pard(1)
c         i29=1
c       endif
c     end do
c     if(i29.eq.1) flush(29)

820   format(f6.3,8(1x,f10.4))

      end

c***************************************************************************
      subroutine jam2hydroall

c...Convert all jam particles into fluid elements.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      include 'fluid.inc'
      parameter(ngmax=25*25*25)
      real*8 uh(0:4),uv(4),gex(ngmax)

      if(nv.le.0) then
        print *,'(jam2hydroall:) nv=',nv
        stop
      endif

      iopt_core=mstc(147)
      isel=0
      if(iopt_core.eq.2) isel=1

c...Initial energy density at (0,0,0)
      if(mstc(8).ge.1) then
        call jamdens(icon,0,cc,rho,rhob0,einv0,uv,isel)
      endif

      tau0=pard(1)
      tau=1.0d0
      if(opt_taueta.eq.1) tau=tau0

      vol=dx*dy*dz
      volh=vol*hbc
      sigma=parc(142)
c     gnorm=(1.0/(2*paru(1))**(1.5))/sigma**3

      do i=1,nv  ! loop over all particles.

        if(abs(k(7,i)).eq.1) cycle   ! not yet interact
        if(k(1,i).eq.0.or.k(1,i).ge.11) cycle
        if(p(5,i).le.1d-5) cycle     ! photons

        if(k(2,i).eq.92) then ! String
          px=(vq(1,i)+vq(6,i))
          py=(vq(2,i)+vq(7,i))
          pz=(vq(3,i)+vq(8,i))
          pe=(vq(4,i)+vq(9,i))
        else
          px=p(1,i)
          py=p(2,i)
          pz=p(3,i)
          pe=p(4,i)
        endif

        dtp=pard(1)-r(4,i)
        r1=r(1,i)+dtp*px/pe
        r2=r(2,i)+dtp*py/pe
        r3=r(3,i)+dtp*pz/pe

        inside=isinside(r1,r2,r3,ix,iy,iz)
        if(inside.eq.0) cycle  ! outside the cell.

c       if(abs(k(1,i)).ge.11) cycle ! const.quark

        if(iopt_core.ge.1) then
          if(abs(k(7,i)).eq.1) cycle   ! Exclude spectator
          if(iopt_core.eq.2) then
            if(pard(1).lt.r(5,i)) cycle  !...pre-formed hadrons.
          else if(iopt_core.eq.3) then
            if(dtp.lt.0.0d0) cycle  ! not formed
            if(pard(1).lt.r(5,i)) cycle  !...pre-formed hadrons.
          endif
          call jamdens(icon,i,cc,rho,rhob,einv,u,isel)
          if(icon.ne.0) cycle  ! density is too low.

c         if(einv.ge.parc(144)) k(7,i)=-1500
          if(einv.ge.edminfluid(rhob)) k(7,i)=-1500

        else
          k(7,i)=-1500
        endif

        if(k(7,i).eq.-1500) then
          r(1,i)=r1
          r(2,i)=r2
          r(3,i)=r3
          r(4,i)=pard(1)
        endif

      end do


c...Cut off parameter for Gaussian smearing.
      ip=int(parc(147)/dx+0.5)
      iq=int(parc(147)/dy+0.5)
      ir=int(parc(147)/dz+0.5)

      call reset_fluid_val

      nvsave=0
      etotj=0.0
      pxj=0.0
      pyj=0.0
      pzj=0.0
      bjam=0.0
      sjam=0.0
      cjam=0.0
      do i=1,nv  ! loop over all particles.

        if(k(7,i).ne.-1500) goto 100

        xp=r(1,i)
        yp=r(2,i)
        zp=r(3,i)
        inside=isinside(xp,yp,zp,ix,iy,iz)
        if(inside.eq.0) then
          print *,'jam2hydro: particle is already outside',ix,iy,iz
          goto 100
        endif

        if(k(2,i).eq.92) then ! String
          px=(vq(1,i)+vq(6,i))
          py=(vq(2,i)+vq(7,i))
          pz=(vq(3,i)+vq(8,i))
          pe=(vq(4,i)+vq(9,i))
        else
          px=p(1,i)
          py=p(2,i)
          pz=p(3,i)
          pe=p(4,i)
        endif

        ux=0d0
        uy=0d0
        uz=0d0
        gg=1d0
        if(mstc(146).eq.2) then
          uz=pz/sqrt(pe*pe-pz*pz)
          gg=pe/sqrt(pe*pe-pz*pz)  ! gamma_z
        else if(mstc(146).eq.3) then
          em=sqrt(pe**2-px**2-py**2-pz**2)
          ux=px/em
          uy=py/em
          uz=pz/em
          gg=pe/em
        endif

c       call jamcupda(j,-1,0,0,0) ! remove collision
        k(1,i)=13                 ! Flag for dead particle
        bar=k(9,i)/3.0/mstc(5)        ! Baryon number
        ch  = jamk(1,i)/3.0/mstc(5)   ! Charge
        str = jamk(2,i)/mstc(5)       ! Strangeness
        px=px/mstc(5)
        py=py/mstc(5)
        pz=pz/mstc(5)
        pe=pe/mstc(5)

        mstd(30)=mstd(30)+1
        nvsave=nvsave+1
        etotj=etotj+pe
        pxj=pxj+px
        pyj=pyj+py
        pzj=pzj+pz
        bjam=bjam+bar
        cjam=cjam+ch
        sjam=sjam+str

        if(mstc(146).eq.0) then
          u(1,ix,iy,iz) = u(1,ix,iy,iz) + px/volh
          u(2,ix,iy,iz) = u(2,ix,iy,iz) + py/volh
          u(3,ix,iy,iz) = u(3,ix,iy,iz) + pz/volh
          u(0,ix,iy,iz) = u(0,ix,iy,iz) + pe/volh
          u(4,ix,iy,iz) = u(4,ix,iy,iz) + bar/vol
        else

c....Gauss smearing. First compute normalization so that total energy-momentum
c...of the fluid elements are the same as JAM particles.

        jx1=max(0,ix-ip); jx2=min(maxx,ix+ip)
        jy1=max(0,iy-iq); jy2=min(maxy,iy+iq)
        jz1=max(0,iz-ir); jz2=min(maxz,iz+ir)

        enorm=0d0
        ng=0
        do jx=jx1,jx2
        do jy=jy1,jy2
        do jz=jz1,jz2
          ng=ng+1

c         x1=dx*( jx - origx )
c         y1=dy*( jy - origy )
c         z1=dz*( jz - origz )
c         ex = (xp-x1)**2 + (yp-y1)**2 + (zp-z1)**2 
c    &      + ((xp-x1)*ux+(yp-y1)*uy+(zp-z1)*uz)**2
c         enorm = enorm + exp( -ex/(2*sigma**2) )

c         gex(ng)=gaussinter(jx,jy,jz,xp,yp,zp,ux,uy,uz)
c         enorm=enorm+gex(ng)
          enorm=enorm+gaussinter(jx,jy,jz,xp,yp,zp,ux,uy,uz)

        end do
        end do
        end do

       if(enorm.eq.0d0) then
          print *,'jam2hydroall::enorm?',enorm
          stop
        endif

        ng=0
        do jx=jx1,jx2
        do jy=jy1,jy2
        do jz=jz1,jz2
          ng=ng+1

c         x1=dx*( jx - origx )
c         y1=dy*( jy - origy )
c         z1=dz*( jz - origz )
c         ex = (xp-x1)**2 + (yp-y1)**2 + (zp-z1)**2 
c    &      + ((xp-x1)*ux+(yp-y1)*uy+(zp-z1)*uz)**2
c         ex = exp( -ex/(2*sigma**2) )/enorm/vol

c         ex=gex(ng)/enorm/vol
          ex=gaussinter(jx,jy,jz,xp,yp,zp,ux,uy,uz)/enorm/vol

          if(ex.ge.1d-12) then
          u(1,jx,jy,jz) = u(1,jx,jy,jz) + px/hbc * ex
          u(2,jx,jy,jz) = u(2,jx,jy,jz) + py/hbc * ex
          u(3,jx,jy,jz) = u(3,jx,jy,jz) + pz/hbc * ex
          u(0,jx,jy,jz) = u(0,jx,jy,jz) + pe/hbc * ex
          u(4,jx,jy,jz) = u(4,jx,jy,jz) + bar * ex
          endif

        end do
        end do
        end do

        endif
100   end do

      if(nvsave.eq.0) then
        if(mstc(8).ge.1) then
        print *,'no particle is converted into fluid',nvsave,etotj
        endif
        return
      endif

      call jamedit

c....Update collision list if some particles remain. 
c...(needed afterjamedit).
c     mentry=0
      mstd(29)=mstd(29)+1
      call jamclist

      if(mstc(8).ge.1) then
        write(6,*)'afte jamedit nv=',nv,' nbary=',nbary,
     &   ' nmeson=',nmeson, ' mentry=',mentry
      endif


      pe_tot=0.0
      px_tot=0.0
      py_tot=0.0
      pz_tot=0.0
      b_tot=0.0
      en_tot=0.0
      do ix = 0,maxx
      do iy = 0,maxy
      do iz = 0,maxz
        if(u(0,ix,iy,iz).lt.1d-12) then
         u(:,ix,iy,iz)=0d0
         cycle
        endif

        uh(:)=u(:,ix,iy,iz)/tau
        det = uh(0)*uh(0)-uh(1)*uh(1)-uh(2)*uh(2)-uh(3)*uh(3)
        if(det.le.0d0) then
          print *,'jam2hydroall det<0',det
          print *,'u4h=',uh(0),'u1h=',uh(1),uh(2),uh(3)
        endif

        call thermal(uh,ene,ba,pre,sva,vel)
        call geteos2(ene,ba,tplc,xmulc,smu,slc)

        el(ix,iy,iz) = ene
        bl(ix,iy,iz) = ba
        pl(ix,iy,iz) = pre
        tl(ix,iy,iz) = tplc
        ml(ix,iy,iz) = xmulc
        fh = uh(0)+tau*pre  ! gam^2*(e+p)
        gam = 1d0-(uh(1)*uh(1)+uh(2)*uh(2)+uh(3)*uh(3))/fh/fh
        if((fh .gt.0d0) .and. (gam .gt. 0d0)) then
          gam = sqrt(1d0/gam)
        else
          gam = 1d0
        endif
        bl(ix,iy,iz) = uh(4)/gam/tau
c       if(abs(uh(4)/gam-ba).gt.1d-6) then
c         print *,'ba=',ba,uh(4)/gam,abs(uh(4)/gam-ba),gam
c         stop
c       endif
        vx(ix,iy,iz)=uh(1)/fh
        vy(ix,iy,iz)=uh(2)/fh
        vz(ix,iy,iz)=uh(3)/fh

        px_tot=px_tot+u(1,ix,iy,iz)*volh
        py_tot=py_tot+u(2,ix,iy,iz)*volh
        pz_tot=pz_tot+u(3,ix,iy,iz)*volh
        pe_tot=pe_tot+u(0,ix,iy,iz)*volh
        b_tot =b_tot +u(4,ix,iy,iz)*vol
        en_tot = en_tot + getentropy(ene,bl(ix,iy,iz))*vol*gam
 

c       write(77,965)0d0,dble(ix),dble(iy),dble(iz),
c    &    dx*dy*dz,0d0,0d0,0d0,
c    &    vx(ix,iy,iz),vy(ix,iy,iz),vz(ix,iy,iz),
c    &    pl(ix,iy,iz),ene,tplc,xmulc,smu
c965    format(16e15.6)

200   end do
      end do
      end do

      s_tot=sjam
      c_tot=cjam

c....Flag for hydro
      mste(41)=1

      if(mstc(8).ge.1) then
      print *,'initial jam energy',etotj/mstd(2),' nvsave=',nvsave
      print *,'initial condition for hydro e=',pe_tot/mstd(2)
      print *,'initial hydro e(%) ',abs(etotj-pe_tot)/etotj*100,
     &  '# of jam particle remain=',nv
      print *,'B(jam)=',bjam,'B(hy)=',b_tot
      print *,'p(jam)=',pxj,pyj,pzj
      print *,'p(ini)=',px_tot,py_tot,pz_tot

      print *,'initial cc=',cc,rho
      print *,'initial rhob=',rhob0,
     & ' fluid=',bl(origx,origy,origz),u(4,origx,origy,origz)
      print *,'initial eden=',einv0,' fluid=',
     & el(origx,origy,origz)*hbc,u(0,origx,origy,origz)*hbc
      print *,'entropy = ',en_tot

      endif

c     if(mstc(181).ge.1) then
c     write(30,820)pard(1),
c    &  bl(maxx/2,maxy/2,maxz/2),
c    &  el(maxx/2,maxy/2,maxz/2)*hbc,
c    &  rhob0,einv0,2
c     flush(30)
c820   format(f6.3,4(1x,f10.4),1x,i1)
c      endif
 
      end

c***************************************************************************
      real*8 function gaussinter(jx,jy,jz,xp,yp,zp,ux,uy,uz)
      implicit none
      include 'fluid.inc'
      integer jx,jy,jz,i,j,k
      real*8 xp,yp,zp,ux,uy,uz
      real*8 vol,eg,x1,y1,z1,ex
      real*8 x(5),w(5)

      if(n_gaussp.eq.1) then
        x1=dx*( jx - origx )
        y1=dy*( jy - origy )
        z1=dz*( jz - origz )
        ex = (xp-x1)**2 + (yp-y1)**2 + (zp-z1)**2 
     &      + ((xp-x1)*ux+(yp-y1)*uy+(zp-z1)*uz)**2
        gaussinter= exp( -ex/(2*gsigma**2) )*dx*dy*dz
        return
      endif

      if(n_gaussp.eq.1) then
       x(1)=0d0
       w(1)=2.0d0
      else if(n_gaussp.eq.2) then
        x(1)= .57735026918962576451 ! 1/sqrt(3)
        x(2)=-x(1)
        w(1)=1.0d0
        w(2)=1.0d0
      else if(n_gaussp.eq.3) then
        x(1)=0
        x(2)=-sqrt(3.0d0/5d0)
        x(3)=sqrt(3.0d0/5d0)
        w(1)=8.0d0/9.0d0
        w(2)=5.0d0/9.0d0
        w(3)=5.0d0/9.0d0
      else if(n_gaussp.eq.4) then
        x(1)= 0.33998104358485626480
        x(2)= -x(1)
        x(3)=0.86113631159405257521
        x(4)= -x(3)
        w(1)= 0.65214515486254614262
        w(2)= w(1)
        w(3)= 0.34785484513745385737
        w(4)= w(3)
      else if(n_gaussp.eq.5) then
        x(1)=0
        x(2)=0.53846931010568309103
        x(3)=-x(2)
        x(4)=0.90617984593866399278
        x(5)=-x(4)
        w(1)=0.56888888888888888888
        w(2)=0.47862867049936646804
        w(3)=w(2)
        w(4)=0.23692688505618908751
        w(5)=w(4)
      else
       print *,'gaussinter wrong n_gaussp=',n_gaussp
       stop
      endif

      eg=0.0d0
      do i=1,n_gaussp
      do j=1,n_gaussp
      do k=1,n_gaussp
        x1=dx*( jx - origx ) + x(i)*dx/2
        y1=dy*( jy - origy ) + x(j)*dy/2
        z1=dz*( jz - origz ) + x(k)*dz/2
        ex = (xp-x1)**2 + (yp-y1)**2 + (zp-z1)**2 
     &      + ((xp-x1)*ux+(yp-y1)*uy+(zp-z1)*uz)**2
        eg = eg + w(i)*w(j)*w(k)*exp( -ex/(2*gsigma**2) )
      end do
      end do
      end do
      gaussinter=eg*dx*dy*dz/8


c     np=1
c     np2=2
c     eg=0.0d0
c     do i=-np,np
c     do j=-np,np
c     do k=-np,np
c       x1=dx*( jx - origx ) + i*dx/np2
c       y1=dy*( jy - origy ) + j*dy/np2
c       z1=dz*( jz - origz ) + k*dz/np2
c       ex = (xp-x1)**2 + (yp-y1)**2 + (zp-z1)**2 
c    &      + ((xp-x1)*ux+(yp-y1)*uy+(zp-z1)*uz)**2
c       eg = eg + exp( -ex/(2*gsigma**2) )
c     end do
c     end do
c     end do

c     np=3
c     np2=2
c     eg=0.0d0
c     do i=1,np
c     do j=1,np
c     do k=1,np
c       x1=dx*( jx - origx )+(i-1)*dx/np2
c       y1=dy*( jy - origy )+(j-1)*dy/np2
c       z1=dz*( jz - origz )+(k-1)*dz/np2
c       ex = (xp-x1)**2 + (yp-y1)**2 + (zp-z1)**2 
c    &      + ((xp-x1)*ux+(yp-y1)*uy+(zp-z1)*uz)**2
c       eg = eg + exp( -ex/(2*gsigma**2) )
c     end do
c     end do
c     end do

c     gaussinter=eg/np**3

      end

c***************************************************************************
      subroutine jam2fluid_aformation(ip)

c...Convert particle into fluid element after formation time.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      real*8 jamdtim,u(4)

c....Too late; All fluid was already converted to particle.
      if(mste(41).eq.-1) goto 1000

      if(k(2,ip).eq.92) then
       print *,'jam2fluid_aformation string?',k(1,ip),k(2,ip),k(7,ip)
          pe=vq(4,ip)+vq(9,ip)
          px=vq(1,ip)+vq(6,ip)
          py=vq(2,ip)+vq(7,ip)
          pz=vq(3,ip)+vq(8,ip)
          print *,'m=',sqrt(pe**2-px**2-py**2-pz**2)
      endif

c...Check density.
      inside=isinside(r(1,ip),r(2,ip),r(3,ip),ix,iy,iz)
      if(inside.eq.0) goto 1000

      efluid=eden_fluid(ix,iy,iz)
      bfluid=bden_fluid(ix,iy,iz)
      call jamdens(icon,ip,cc,rho,rhob,einv,u,1)
      bden=bfluid+rhob
c     if(icon.ne.0.or.einv+efluid.lt.parc(144)) goto 1000
      if(icon.ne.0.or.einv+efluid.lt.edminfluid(bden)) goto 1000

      k(1,ip)=13 ! flag for dead particle
      k(7,ip)=0
      bar = k(9,ip)/3.0/mstc(5)     ! baryon number
      kc=jamcomp(k(2,ip))
      ch  = kchg(kc,1)/3*isign(1,k(2,ip))/mstc(5)
      str = kchg(kc,7)*isign(1,k(2,ip))/mstc(5)
      px=p(1,ip)/mstc(5)
      py=p(2,ip)/mstc(5)
      pz=p(3,ip)/mstc(5)
      pe=p(4,ip)/mstc(5)
      if(pe**2-px**2-py**2-pz**2.le.0d0) then
        print *,'jam2fluid_aformation spacelike? kf=',k(2,ip),p(5,ip)
      endif
      call p2fluid(mstc(146),ix,iy,iz,px,py,pz,pe,bar,str,ch)
      mste(41)=mste(41)+1
      call jamcupda(j,-1,0,0,0) ! remove collision

c     call jam2fluid_adecay(ix,iy,iz)

c     print *,'particle conversion ',k(2,ip)

      return

1000  continue
c.....Set resonance decay time.
      k(7,ip)=2
      v(5,ip)=1.d+35
      if(k(1,ip).eq.2) then
        kc=jamcomp(k(2,ip))
        v(5,ip)=r(5,ip)+jamdtim(1,k(2,ip),kc,k(1,ip),p(5,ip),p(4,ip))
      endif
      call jamcupda(ip,-1,0,0,1) ! update collision


      end

c***************************************************************************
      integer function isconv(j,opt)
      implicit none
      include 'jam1.inc'
      integer opt,j

c...Only mesons are converted into fluid
        if(opt.eq.1) then
          isconv=0
          if(k(9,j).eq.0) isconv=1
c....Exclude reading baryons.
        else if(opt.eq.2) then
          isconv=1
          if(abs(k(1,j)).ge.11.and.k(9,j).eq.3) isconv=0 
c...Exclude reading hadrons.
        else if(opt.eq.3) then
          isconv=1
          if(abs(k(1,j)).ge.11) isconv=0
        else
          isconv=1
        endif
        if(p(5,j).le.1d-5) isconv=0     ! photons

        end

c***************************************************************************
      subroutine jam2fluid_decay(einv,efluid,indd,nadd,ind)

c...Convert hadrons from decay into fluid element.

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      integer indd(100),nadd,ind,isconv
      integer opt,i,j,iconv,inside,ix,iy,iz,kc,isinside,jamcomp
      double precision einv,efluid,em,ttime,ftime,dtp,r1,r2,r3
      double precision bar,ch,str,px,py,pz,pe

      opt=mod(mstc(142),10)

      do i=1,nadd  ! Loop over decaying particles.
        j=indd(i)
        iconv=isconv(j,opt)

c       if(mstc(142).gt.10.and.mste(7).eq.0) iconv=0

c...mstc(7)=0 means baryons from resonance decay.
c       if(mste(7).eq.0) then
c         if(k(9,j).eq.3.and.efluid.lt.parc(144)) iconv=0
c         if(k(9,j).eq.3.and.abs(k(7,j)).le.3) iconv=0
c         if(abs(k(7,j)).le.5) iconv=0
c         if(abs(k(7,j)).le.mstc(139)) iconv=0
c       endif

c....This particle will be converted into fluid after formation time.
        if(iconv.eq.1) then

c       ttime=0d0
c       if(parc(146).gt.0d0) then
c       em=sqrt(p(4,j)**2-p(1,j)**2-p(2,j)**2-p(3,j)**2)
c       ttime=parc(146)*p(4,j)/em
c       endif
c       ftime=r(4,j)
c       if(opt.eq.4.and.abs(k(1,j)).ge.11) ftime=r(5,j)

        ttime=0.0
        ftime=r(4,j)
c...Set thermalization time in case of string decay.
        if(parc(146).ge.0d0) then
          if(mste(7).eq.1) then
          em=sqrt(p(4,j)**2-p(1,j)**2-p(2,j)**2-p(3,j)**2)
          ttime=parc(146)*p(4,j)/em
          ftime=pard(1)
c         write(45,*) ttime
          endif

c...Particles will be converted into fluid after their formation time.
        else
c         write(44,*) ftime-pard(1)
c....opt=4: leading hadron will be converted into hadron after formation time.
          if(opt.eq.4.or.opt.eq.5) then
            if(abs(k(1,j)).ge.11) ftime=r(5,j)
          endif
c....opt=5: leading meson will be converted into hadron after formation time.
          if(opt.eq.6.and.k(9,j).eq.0.and.abs(k(1,j)).ge.11)ftime=r(5,j)
        endif

        dtp=pard(1)-ftime-ttime
        if(dtp.lt.0d0) then
          k(7,j)=-1200
          v(5,j)=ftime+ttime

c...opt=6:this particle can collide before fluization within formation time.
c         if(opt.eq.5.and.abs(k(1,j)).ge.11.and.k(9,j).eq.3) 
          if(opt.eq.5.and.abs(k(1,j)).ge.11) 
     &          call jamcupda(j,-1,ind,0,1) ! update collision
          cycle
        else
          r1=r(1,j)+dtp*p(1,j)/p(4,j)
          r2=r(2,j)+dtp*p(2,j)/p(4,j)
          r3=r(3,j)+dtp*p(3,j)/p(4,j)
          inside=isinside(r1,r2,r3,ix,iy,iz)
          if(inside.eq.0) iconv=0
        endif

        endif


        if(iconv.eq.1) then
          k(1,j)=13 ! flag for dead particle
          bar = k(9,j)/3.0/mstc(5)     ! baryon number
          kc=jamcomp(k(2,j))
          ch  = kchg(kc,1)/3*isign(1,k(2,j))/mstc(5)
          str = kchg(kc,7)*isign(1,k(2,j))/mstc(5)
          px=p(1,j)/mstc(5)
          py=p(2,j)/mstc(5)
          pz=p(3,j)/mstc(5)
          pe=p(4,j)/mstc(5)
          if(pe**2-px**2-py**2-pz**2.le.0d0) then
           print *,'jam2fluid_decay lightlike? kf=',k(2,j),p(5,j)
          endif
          call p2fluid(mstc(146),ix,iy,iz,px,py,pz,pe,bar,str,ch)
          mste(41)=mste(41)+1
          call jamcupda(j,-1,0,0,0) ! remove collision

c         call jam2fluid_adecay(ix,iy,iz)

c...We do not kill this particle.
        else
          call jamcupda(j,-1,ind,0,1) ! update collision
        endif
      end do

      end

c***************************************************************************
      subroutine jam2fluid_adecay(jx,jy,jz)

c...Convert particles which are located in the same fluid cell for
c...decaying particle.

      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      integer i,jamk,nadd,isconv,opt,inside,isinside,ix,iy,iz,kc,jamcomp
      integer jx,jy,jz
      real*8 ctime,r1,r2,r3,dt,px,py,pz,pe,bar,ch,str

      nadd=0
      ctime=pard(1)
      opt=mod(mstc(142),10)

      do i=1,nv  ! loop over all particles.

        if(abs(k(7,i)).eq.1) cycle   ! not yet interact
        if(k(1,i).gt.10.or.k(1,i).eq.0) cycle       ! dead particle
        if(p(5,i).le.1d-5) cycle     ! photons
        if(abs(k(1,i)).ge.11) cycle ! const.quark
        if(r(5,i).gt.ctime) cycle !.....pre-formed hadrons.
        if(isconv(i,opt).eq.0) cycle

        dt=ctime-r(4,i)
        if(dt.lt.0.0d0) cycle  ! not formed
        r1=r(1,i)+dt*p(1,i)/p(4,i)
        r2=r(2,i)+dt*p(2,i)/p(4,i)
        r3=r(3,i)+dt*p(3,i)/p(4,i)
        inside=isinside(r1,r2,r3,ix,iy,iz)
        if(inside.eq.0) cycle  ! outside the cell.
        if(ix.ne.jx) cycle
        if(iy.ne.jy) cycle
        if(iz.ne.jz) cycle

c...Convert particle to fluid element.
          px=p(1,i)/mstc(5)
          py=p(2,i)/mstc(5)
          pz=p(3,i)/mstc(5)
          pe=p(4,i)/mstc(5)
          bar = k(9,i)/3.0/mstc(5)     ! baryon number
          kc=jamcomp(k(2,i))
          ch  = kchg(kc,1)/3*isign(1,k(2,i))/mstc(5)
          str = kchg(kc,7)*isign(1,k(2,i))/mstc(5)
          call p2fluid(mstc(146),ix,iy,iz,px,py,pz,pe,bar,str,ch)
          call jamcupda(i,-1,0,0,0) ! remove this from collision list
          k(1,i)=13 ! flag for dead particle
          k(7,i)=0
c         mstd(30)=mstd(30)+1          ! update number of dead particles.
          nadd=nadd+1
      end do

c...Flag for fluid evolution. =1 meas that there are fluid elements.
      mste(41)=mste(41)+nadd

      end

c***************************************************************************
c...Under construction.
      subroutine jam2fluid_dynamical2

c..Compute source of fluid element.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      real*8 u(4),pini(5)
      integer opt,jamk

      opt=2

c     it=mstd(140)
c      print *,'it=',it,nadd

      dtp=parc(2)
      nadd=0
      nc=0
      ctime=pard(1)

      pxtot=0d0
      pytot=0d0
      pztot=0d0
      petot=0d0
      btoto=0d0
      do i=1,nv  ! loop over all particles.

        if(abs(k(7,i)).eq.1) cycle   ! not yet interact
        if(k(1,i).gt.10) cycle   ! dead particle
        if(p(5,i).le.1d-5) cycle

c.....pre-formed hadrons.
        if(mstc(89).eq.1.and.r(5,i).gt.ctime) cycle
        if(ctime-r(4,i).lt.0.0d0) cycle

        inside=isinside(r(1,i),r(2,i),r(3,i),ix,iy,iz)
        if(inside.eq.0) cycle  ! outside the cell.

        call jamdens(icon,i,cc,rho,rhob,einv,u,1)
        if(icon.ne.0) cycle  ! density is too low.
c       if(rho.le.0.2) cycle
c       if(einv.lt.0.5) cycle

c       call qbfac(i,facq,facb)

        if(opt.eq.1) then
          be1=u(1)/u(4)
          be2=u(2)/u(4)
          be3=u(3)/u(4)
          gam=u(4)
          px0=p(1,i)
          py0=p(2,i)
          pz0=p(3,i)
          pe0=p(4,i)
          pini(:)=p(:,i)
c         y=0.5*log((pe0+pz0)/(pe0-pz0))
          call jamrobo(0d0,0d0,-be1,-be2,-be3,gam,px0,py0,pz0,pe0)
          psq=px0**2+py0**2+pz0**2
          emsq=pe0**2-psq
          em=sqrt(emsq)
          pp=sqrt(psq)
          de=parc(143)*rho*dtp
          pe = pe0-de
          iconv=1
          if(pe.le.p(5,i)+0.0001) iconv=0
        else
          iconv=0
           efluid=eden_fluid(ix,iy,iz)
          if(einv+efluid.lt.0.5) cycle
        endif

c...Convert particle to fluid element.
        if(iconv.eq.0) then
          px=p(1,i)/mstc(5)
          py=p(2,i)/mstc(5)
          pz=p(3,i)/mstc(5)
          pe=p(4,i)/mstc(5)
          bar = k(9,i)/3.0/mstc(5)     ! baryon number
c         ch  = jamk(1,i)/3.0/mstc(5)  ! charge
c         str  = jamk(2,i)/mstc(5)      ! strangeness
          kc=jamcomp(k(2,i))
          ch  = kchg(kc,1)/3*isign(1,k(2,i))/mstc(5)
          str = kchg(kc,7)*isign(1,k(2,i))/mstc(5)
          call p2fluid(mstc(146),ix,iy,iz,px,py,pz,pe,bar,str,ch)
          call jamcupda(i,-1,0,0,0) ! remove this from collision list
          k(1,i)=13 ! flag for dead particle
c         mstd(30)=mstd(30)+1          ! update number of dead particles.
          nadd=nadd+1
          nc=nc+1
          pxtot=pxtot+px
          pytot=pytot+py
          pztot=pztot+pz
          petot=petot+pe
          btot =btot+bar

c         print *,'p->f',k(2,i),einv,efluid,pe

c....Part of energy is put into fluid.
        else

          bar=0d0
          str=0d0
          ch=0d0
          dp = sqrt(pe**2-p(5,i)**2)
          px = px0*dp/pp
          py = py0*dp/pp
          pz = pz0*dp/pp
          dpx=px0-px
          dpy=py0-py
          dpz=pz0-pz
          dpe=pe0-pe
          print *,'fluid m=',dpe**2-dpx**2-dpy**2-dpz**2
          call jamrobo(0d0,0d0,be1,be2,be3,gam,px,py,pz,pe)
          p(1,i)=px
          p(2,i)=py
          p(3,i)=pz
          p(4,i)=pe

c         print *,'ei=',pini(4),sqrt(pini(1)**2+pini(2)**2+pini(3)**2)

c         ppf=px**2+py**2+pz**2
c         print *,'ep=',sqrt(emsq+ppf),pe,sqrt(ppf),
c    &  sqrt(pe**2-ppf),p(5,i),k(2,i)
c         y2=0.5*log((pe+pz)/(pe-pz))
c         print *,'y=',y,' y2=',y2

          call jamcupda(i,-1,i,0,1) ! update collision
          nadd=nadd+1
          dpx=(pini(1)-px)/mstc(5)
          dpy=(pini(2)-py)/mstc(5)
          dpz=(pini(3)-pz)/mstc(5)
          dpe=(pini(4)-pe)/mstc(5)

          print *,'ef=',dpe,sqrt(dpx**2+dpy**2+dpz**2)
          det=dpe**2-dpx**2-dpy**2-dpz**2
          if(det.lt.0d0) then
            print *,'fluid element is spacelike',det
            stop
          endif

          call p2fluid(mstc(146),ix,iy,iz,dpx,dpy,dpz,dpe,bar,str,ch)

          pxtot=pxtot+ dpx
          pytot=pytot+ dpy
          pztot=pztot+ dpz
          petot=petot+ dpe

c         print *,'energy loss',k(1,i),r(1,i),r(2,i),r(3,i),de
        endif

      end do

c...Flag for fluid evolution. =1 meas that there are fluid elements.
      if(nadd.ge.1) mste(41)=1

c     if(nadd.ge.1) then
c       call hydroTotalEnergy(it,1d0,eh,pxh,pyh,pzh,bh)
c       print *,'jam2fluid_dynamical',nc,nadd
c       print *,pxtot,pytot,pztot,petot
c       print *,pxh,pyh,pzh,eh
c       call jamdisp(6,2)
c     endif


      end

c***************************************************************************
      subroutine jam2fluid_dynamical1

c...Compute source of fluid element at which energy density exceed
c...a critical value.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      integer jamk
      real*8 u(4)

      nadd=0
      ctime=pard(1)

      pxtot=0d0
      pytot=0d0
      pztot=0d0
      petot=0d0
      btoto=0d0

c...First check which particle will be converted into fluid.
      do i=1,nv  ! loop over all particles.

        if(abs(k(7,i)).eq.1) cycle   ! not yet interact
        if(k(1,i).gt.10.or.k(1,i).eq.0) cycle       ! dead particle
        if(p(5,i).le.1d-5) cycle     ! photons
        if(abs(k(1,i)).ge.11) cycle ! const.quark

c.....pre-formed hadrons.
        if(r(5,i).gt.ctime) cycle
        dt=ctime-r(4,i)
        if(dt.lt.0.0d0) cycle  ! not formed
        r1=r(1,i)+dt*p(1,i)/p(4,i)
        r2=r(2,i)+dt*p(2,i)/p(4,i)
        r3=r(3,i)+dt*p(3,i)/p(4,i)

        inside=isinside(r1,r2,r3,ix,iy,iz)
        if(inside.eq.0) cycle  ! outside the cell.

        call jamdens(icon,i,cc,rho,rhob,einv,u,1)
        if(icon.ne.0) cycle  ! density is too low.

c...Fluid conversion above critical energy density.
         efluid=eden_fluid(ix,iy,iz)
         bfluid=bden_fluid(ix,iy,iz)

c        if(einv+efluid.ge.parc(144)) k(7,i)=-1000
         if(einv+efluid.ge.edminfluid(bfluid+rhob)) k(7,i)=-1000

      end do

      do i=1,nv  ! loop over all particles.
c...Convert particle to fluid element.
      if(k(7,i).eq.-1000) then
          px=p(1,i)/mstc(5)
          py=p(2,i)/mstc(5)
          pz=p(3,i)/mstc(5)
          pe=p(4,i)/mstc(5)
          bar = k(9,i)/3.0/mstc(5)     ! baryon number
          kc=jamcomp(k(2,i))
          ch  = kchg(kc,1)/3*isign(1,k(2,i))/mstc(5)
          str = kchg(kc,7)*isign(1,k(2,i))/mstc(5)

          dt=ctime-r(4,i)
          r1=r(1,i)+dt*p(1,i)/p(4,i)
          r2=r(2,i)+dt*p(2,i)/p(4,i)
          r3=r(3,i)+dt*p(3,i)/p(4,i)
          inside=isinside(r1,r2,r3,ix,iy,iz)

          call p2fluid(mstc(146),ix,iy,iz,px,py,pz,pe,bar,str,ch)
          call jamcupda(i,-1,0,0,0) ! remove this from collision list
          k(1,i)=13 ! flag for dead particle
          k(7,i)=0
c         mstd(30)=mstd(30)+1          ! update number of dead particles.
          nadd=nadd+1
          pxtot=pxtot+px
          pytot=pytot+py
          pztot=pztot+pz
          petot=petot+pe
          btot =btot+bar
      endif
      end do

c...Flag for fluid evolution. =1 meas that there are fluid elements.
      if(nadd.ge.1) mste(41)=1

      end

c***************************************************************************
      subroutine jam2fluid_dynamical3

c...Absorb particle into fluid when a particle pass through fluid element.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      integer jamk,jamchge
      real*8 u(4)

      nadd=0
      ctime=pard(1)
      iopt=mod(mstc(148),10)

      pxtot=0d0
      pytot=0d0
      pztot=0d0
      petot=0d0
      btoto=0d0
      do i=1,nv  ! loop over all particles.

        if(abs(k(7,i)).eq.1) cycle   ! not yet interact
        if(k(1,i).gt.10.or.k(1,i).eq.0) cycle       ! dead particle
        if(p(5,i).le.1d-5) cycle     ! photons
c       if(abs(k(1,i)).ge.11) cycle ! const.quark
        if(r(5,i).gt.ctime) cycle !.....pre-formed hadrons.
        if(isconv(i,iopt).eq.0) cycle

c.....pre-formed hadrons.
c       if(mstc(89).eq.1.and.r(5,i).gt.ctime) cycle
        dt=ctime-r(4,i)
        if(dt.lt.0.0d0) cycle  ! not formed
        if(ctime.lt.r(5,i)) cycle

c       if(mstc(143).ge.2.and.abs(k(1,i)).ge.11) cycle ! const.quark

c       if(mstc(143).eq.2) then
c         kc=jamcomp(k(2,i))
c         if(k(1,i).eq.2.or.mdcy(kc,1).eq.1)  cycle
c         if(k(1,i).ne.1)  cycle
c       endif

        rx=r(1,i)
        ry=r(2,i)
        rz=r(3,i)
        if(k(1,i).eq.3) then
          pe=vq(4,i)+vq(9,i)
          px=vq(1,i)+vq(6,i)
          py=vq(2,i)+vq(7,i)
          pz=vq(3,i)+vq(8,i)
        else if(k(1,i).eq.4) then
          rx=0.0d0
          ry=0.0d0
          rz=0.0d0
          px=0.0d0
          py=0.0d0
          pz=0.0d0
          pe=0.0d0
          mpa=k(11,i)-k(10,i)+1
          do j=k(10,i),k(11,i)
            rx=rx+r(1,j)/mpa
            ry=ry+r(2,j)/mpa
            rz=rz+r(3,j)/mpa
            px=px+p(1,j)
            py=py+p(2,j)
            pz=pz+p(3,j)
            pe=pe+p(4,j)
          end do
        else
         pe=p(4,i)
         px=p(1,i)
         py=p(2,i)
         pz=p(3,i)
        endif

        r1=rx+dt*px/pe
        r2=ry+dt*py/pe
        r3=rz+dt*pz/pe

        inside=isinside(r1,r2,r3,ix,iy,iz)
        if(inside.eq.0) cycle  ! outside the cell.
c       if(pe.gt.parc(142)) cycle

         efluid=eden_fluid(ix,iy,iz)
         bfluid=bden_fluid(ix,iy,iz)

         if(efluid.lt.edminfluid(bfluid)) cycle

c....Particle is absorbed by the fluid.
         if(mstc(143).le.2) then

          if(mstc(143).eq.2) then
             yp=0.5*log((pe+pz)/(pe-pz))
             yf=rapidity_fluid(ix,iy,iz)
             if(abs(yp-yf).gt.1.0) cycle
          endif

          if(k(1,i).eq.3) then
            ch=jamchge(kq(1,i))+jamchge(kq(2,i))/3
            str=kfprop(kq(1,i),3)+kfprop(kq(2,i),3)
          else if(k(1,i).eq.4) then
            ch=jamchge(k(2,i))/3
            str=kfprop(k(2,i),3)
          else
            kc=jamcomp(k(2,i))
            ch  = kchg(kc,1)/3*isign(1,k(2,i))
            str = kchg(kc,7)*isign(1,k(2,i))
          endif

          px=px/mstc(5)
          py=py/mstc(5)
          pz=pz/mstc(5)
          pe=pe/mstc(5)
          bar = k(9,i)/3.0/mstc(5)     ! baryon number
          ch=ch/mstc(5)
          str=str/mstc(5)

          call p2fluid(mstc(146),ix,iy,iz,px,py,pz,pe,bar,str,ch)
          call jamcupda(i,-1,0,0,0) ! remove this from collision list
          if(k(1,i).eq.4) then
            do j=k(10,i),k(11,i)
              k(1,j)=13
            end do
          endif
          k(1,i)=13 ! flag for dead particle

c         mstd(30)=mstd(30)+1          ! update number of dead particles.
          nadd=nadd+1
          pxtot=pxtot+px
          pytot=pytot+py
          pztot=pztot+pz
          petot=petot+pe
          btot =btot+bar


          else

          call fluid_tensor(ix,iy,iz,pe,px,py,pz)
          a=0.2
          bar=0.0
          str=0.0
          ch=0.0
          px=-a*px
          py=-a*py
          pz=-a*pz
          pe=-a*pe
          px1=p(1,i)-px
          py1=p(2,i)-py
          pz1=p(3,i)-pz
          pe1=p(4,i)-pe
          pp2=px1**2+py1**2+pz1**2
          if(pp2.ge.pe1*pe1) then
c           print *,'dynamical3 pp2> e',pe1,sqrt(pp2)
            cycle
          endif
          em1=sqrt(pe1*pe1-pp2)
          kf1=k(2,i)
          call jamidres(2,kf1,em1,icon)
c         print *,'dynamical3 m=',p(5,i),k(2,i)
          if(icon.ne.0) then
c           print *,'dynamical3 after jamidres icon=',icon,
c    &            'm=',em1,kf1
            cycle
          endif

c....Update fluid element.
          call p2fluid(0,ix,iy,iz,px,py,pz,pe,bar,str,ch)

c.....Update particle.
          p(1,i)=px1
          p(2,i)=py1
          p(3,i)=pz1
          p(4,i)=pe1
          p(5,i)=em1
          k(2,i)=kf1
          call jamcupda(i,-1,i,0,1) ! update collision
          nadd=nadd+1

          endif

c         print *,'p->f',k(2,i),einv,efluid,pe

c       endif

      end do  ! end loop over particles.

c...Flag for fluid evolution. =1 meas that there are fluid elements.
      if(nadd.ge.1) mste(41)=1

c     if(nadd.ge.1) then
c       call hydroTotalEnergy(it,1d0,eh,pxh,pyh,pzh,bh)
c       print *,'jam2fluid_dynamical',nadd,nv
c       print *,pxtot,pytot,pztot,petot
c       print *,pxh,pyh,pzh,eh
c       call jamdisp(6,2)
c     endif


      end

c***************************************************************************
      integer function isinside(x,y,z,ix,iy,iz)
      implicit none
      include 'fluid.inc'
      integer ix,iy,iz
      double precision x,y,z

      isinside=0
      ix = int(( x + dx/2d0)/dx+origx)
      iy = int(( y + dy/2d0)/dy+origy)        
      iz = int(( z + dz/2d0)/dz+origz)
      if(ix.lt.0.or.ix.gt.maxx) return
      if(iy.lt.0.or.iy.gt.maxy) return
      if(iz.lt.0.or.iz.gt.maxz) return
      isinside=1
      end

c***************************************************************************
      double precision function eden_fluid(ix,iy,iz)
      implicit none
      include 'fluid.inc'
      integer ix,iy,iz

      eden_fluid=el(ix,iy,iz)*hbc

      end

c***************************************************************************
      double precision function bden_fluid(ix,iy,iz)
      implicit none
      include 'fluid.inc'
      integer ix,iy,iz

      bden_fluid=bl(ix,iy,iz)

      end

c***************************************************************************
      double precision function rapidity_fluid(ix,iy,iz)
      implicit none
      include 'fluid.inc'
      integer ix,iy,iz

      rapidity_fluid=0.5*log( (1.0+vz(ix,iy,iz))/(1.0-vz(ix,iy,iz)) )

      end

c***************************************************************************
      subroutine fluid_tensor(ix,iy,iz,pe,px,py,pz)
      implicit none
      include 'fluid.inc'
      integer ix,iy,iz
      real*8 pe,px,py,pz,vol

      vol=dx*dy*dz
      pe=u(0,ix,iy,iz)*hbc*vol
      px=u(1,ix,iy,iz)*hbc*vol
      py=u(2,ix,iy,iz)*hbc*vol
      pz=u(3,ix,iy,iz)*hbc*vol

      end

c***************************************************************************
      subroutine p2fluid(opt_gauss,ix,iy,iz,px,py,pz,pe,bar,str,ch)
      implicit none
      include 'fluid.inc'
      integer inside,opt_gauss
      integer ix,iy,iz,ip,iq,ir,jx,jy,jz,jx1,jx2,jy1,jy2,jz1,jz2
      double precision bar,str,ch,fg,va,det,vel,ux,uy,uz,gg
      double precision vol,px,py,pz,pe,tau
      double precision u1h,u2h,u3h,u4h,u5h,fh,em
      double precision ene,pre,nba,sva,tplc,mulc,smu,slc
      double precision enorm,ex,gg2,x1,y1,z1,xp,yp,zp,m
      double precision gaussinter

      tau=1.0
c     if(opt_taueta.eq.1) tau=tau0 + it*dtt/3.0
      vol=dx*dy*dz

c...Update total fluid energy-momentum and baryon number.
      px_tot=px_tot+px
      py_tot=py_tot+py
      pz_tot=pz_tot+pz
      pe_tot=pe_tot+pe
      b_tot =b_tot +bar
      s_tot =s_tot +str
      c_tot =c_tot +ch

      if(opt_gauss.eq.0) then

      u(1,ix,iy,iz) = u(1,ix,iy,iz) + px/vol/hbc
      u(2,ix,iy,iz) = u(2,ix,iy,iz) + py/vol/hbc
      u(3,ix,iy,iz) = u(3,ix,iy,iz) + pz/vol/hbc
      u(0,ix,iy,iz) = u(0,ix,iy,iz) + pe/vol/hbc
      u(4,ix,iy,iz) = u(4,ix,iy,iz) + bar/vol

      call thermal(u(:,ix,iy,iz),ene,nba,pre,sva,vel)
      call geteos2(ene,nba,tplc,mulc,smu,slc)
      fg = u(0,ix,iy,iz) + tau*pre
      vx(ix,iy,iz) = u(1,ix,iy,iz)/fg

c     vy(iy,iy,iz) = u(2,ix,iy,iz)/fg
c     vz(iz,iy,iz) = u(3,ix,iy,iz)/fg
      vy(ix,iy,iz) = u(2,ix,iy,iz)/fg
      vz(ix,iy,iz) = u(3,ix,iy,iz)/fg

      cs(ix,iy,iz) = sva
      el(ix,iy,iz) = ene
      bl(ix,iy,iz) = nba
      pl(ix,iy,iz) = pre
      tl(ix,iy,iz) = tplc
      ml(ix,iy,iz) = mulc

      else

c...Cut off parameter for Gaussian smearing.
      ip=int(dcut/dx+0.5)
      iq=int(dcut/dy+0.5)
      ir=int(dcut/dz+0.5)

      xp=dx*( ix - origx )
      yp=dy*( iy - origy )
      zp=dz*( iz - origz )

c....Gauss smearing.
        ux=0d0
        uy=0d0
        uz=0d0
        gg=1.0d0
        if(opt_gauss.eq.2) then
          uz=pz/sqrt(pe*pe -pz*pz)
          gg=pe/sqrt(pe*pe - pz*pz)
        else if(opt_gauss.eq.3) then
          em=sqrt(pe**2-px**2-py**2-pz**2)
          ux=px/em
          uy=py/em
          uz=pz/em
          gg=pe/em
        endif

        jx1=max(0,ix-ip); jx2=min(maxx,ix+ip)
        jy1=max(0,iy-iq); jy2=min(maxy,iy+iq)
        jz1=max(0,iz-ir); jz2=min(maxz,iz+ir)
        enorm=0d0
        do jx=jx1,jx2
        do jy=jy1,jy2
        do jz=jz1,jz2
c         x1=dx*( jx - origx )
c         y1=dy*( jy - origy )
c         z1=dz*( jz - origz )
c         ex = (xp-x1)**2 + (yp-y1)**2 + (zp-z1)**2 
c    &      + ( ((xp-x1)*ux+(yp-y1)*uy+(zp-z1)*uz) )**2
c         ex = exp( -ex/(2*gsigma**2) )
c         enorm = enorm + ex*vol
          enorm=enorm+gaussinter(jx,jy,jz,xp,yp,zp,ux,uy,uz)
        end do
        end do
        end do

       if(enorm.eq.0d0) then
          print *,'p2hydro::enorm?',enorm
          stop
        endif

        do jx=jx1,jx2
        do jy=jy1,jy2
        do jz=jz1,jz2
c         x1=dx*( jx - origx )
c         y1=dy*( jy - origy )
c         z1=dz*( jz - origz )
c         ex = (xp-x1)**2 + (yp-y1)**2 + (zp-z1)**2 
c    &      + ((xp-x1)*ux+(yp-y1)*uy+(zp-z1)*uz)**2
c         ex = exp( -ex/(2*gsigma**2) )/enorm
          ex=gaussinter(jx,jy,jz,xp,yp,zp,ux,uy,uz)/enorm/vol

          if(ex.lt.1d-12) goto 106

          u(1,jx,jy,jz) = u(1,jx,jy,jz) + px/hbc * ex
          u(2,jx,jy,jz) = u(2,jx,jy,jz) + py/hbc * ex
          u(3,jx,jy,jz) = u(3,jx,jy,jz) + pz/hbc * ex
          u(0,jx,jy,jz) = u(0,jx,jy,jz) + pe/hbc * ex
          u(4,jx,jy,jz) = u(4,jx,jy,jz) + bar * ex

          det=u(0,jx,jy,jz)**2-u(1,jx,jy,jz)**2-u(2,jx,jy,jz)**2
     &          -u(3,jx,jy,jz)**2
          if(det.le.0d0) then
            print *,'p2fluid lightlike? =',bar,str,ch
            print *,'det',det,pe**2-px**2-py**2-pz**2
            print *,'u0=',u(0,jx,jy,jz)
            print *,'u=',
     &    sqrt(u(1,jx,jy,jz)**2+u(2,jx,jy,jz)**2+u(3,jx,jy,jz)**2)
           stop
          endif

          call thermal(u(:,jx,jy,jz),ene,nba,pre,sva,vel)
          call geteos2(ene,nba,tplc,mulc,smu,slc)

          fg = u(0,jx,jy,jz) + tau*pre
          vx(jx,jy,jz) = u(1,jx,jy,jz)/fg
          vy(jx,jy,jz) = u(2,jx,jy,jz)/fg
          vz(jx,jy,jz) = u(3,jx,jy,jz)/fg

          va=sqrt(vx(jx,jy,jz)**2 + vy(jx,jy,jz)**2 + vz(jx,jy,jz)**2)
          if(va.ge.1.0d0) then
            print *,'p2fulid v>=1?',va,fg,vel
            print *,'vx=',vx(jx,jy,jz),vy(jx,jy,jz),vz(jx,jy,jz)
          endif

          cs(jx,jy,jz)=sva
          el(jx,jy,jz) = ene
          bl(jx,jy,jz) = nba
          pl(jx,jy,jz) = pre
          tl(jx,jy,jz) = tplc
          ml(jx,jy,jz) = mulc

106     end do
105     end do
104     end do

       endif

      end 

c***********************************************************************
      subroutine jamdens2(icon,x0,y0,z0,cc,rho,rhob,einv,u,iopt)

c...Purpose: compute energy-momentum tensor and hydrodynamics velocity.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      dimension pv(4)
      dimension g(4,4),cur(4),curb(4),tens(4,4),u(4)
      dimension gg(4,4),dl(4,4),z(4),uu(0:4)
      data g/ -1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,1/
      data gg/ 1,1,1,-1,  1,1,1,-1, 1,1,1,-1, -1,-1,-1,1/
      data dl/ 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1/
      real*8 xg3(3),wg3(3)
      data xg3/0.0d0, -.77459666924148337703,.77459666924148337703/
      data wg3/.8888888888888888,.55555555555555,.555555555555555/
      common/jamdattens/cur,curb,tens
      save/jamdattens/

c     xg3(1)=0
c     xg3(2)=-sqrt(3.0d0/5d0)
c     xg3(3)=sqrt(3.0d0/5d0)
c     wg3(1)=8.0d0/9.0d0
c     wg3(2)=5.0d0/9.0d0
c     wg3(3)=5.0d0/9.0d0

      widg=2*parc(84)**2 ! Gaussian width
c     widcof=(1.d0/(3.14159d0*widg))**1.5d0)
      widcof=pard(90)

      icon=1
      einv=0.0d0
      rho=0d0
      rhob=0d0
      pfree=0.0d0
      pxx=0.0d0
      pyy=0.0d0
      pzz=0.0d0
      pperp=0d0
      plong=0d0
      ctime=pard(1)

      dcut=parc(147)
      dx=pard(143)
      dy=pard(144)
      dz=pard(145)
      ip=int(dcut/dx+0.5)
      iq=int(dcut/dy+0.5)
      ir=int(dcut/dz+0.5)


      do i=1,4
        cur(i)=0.0d0
        curb(i)=0.0d0
        do j=1,4
          tens(i,j)=0.0d0
      end do
      end do

      isfluid=isinside(x0,y0,z0,ix0,iy0,iz0)
      vol=pard(143)*pard(144)*pard(145)

      ncount=0


c....Loop over all particles
      do i=1,nv
 
        if(k(1,i).gt.10.or.k(1,i).eq.0) goto 100   ! dead particle
        if(p(5,i).le.1d-5) goto 100

c.....pre-formed hadrons.
        if(mstc(89).eq.1.and.r(5,i).gt.ctime) goto 100

        if(abs(k(7,i)).eq.1) goto 100   ! not yet interact

        dt=ctime-r(4,i)
        if(iopt.eq.1.and.dt.lt.0.0d0) goto 100

        if(k(2,i).eq.92) then
          do j=1,4
            pv(j)=vq(j,i)+vq(j+5,i)
          end do
        else
          pv(1)=p(1,i)
          pv(2)=p(2,i)
          pv(3)=p(3,i)
          pv(4)=p(4,i)
        endif

c       y1=0.5*log((pv(4)+pv(3))/(pv(4)-pv(3)))
c       if(abs(y1-ycm).gt.parc(30)) goto 100

        bar=k(9,i)/3.0d0
        facq=1.0d0
        facb=1.0d0
        if(iopt.eq.1.and.r(5,i).gt.time) then
          call qbfac(i,facq,facb)
        endif

        x1=r(1,i) + dt*pv(1)/pv(4)
        y1=r(2,i) + dt*pv(2)/pv(4)
        z1=r(3,i) + dt*pv(3)/pv(4)

        if(mstc(146).eq.0) then

          isfluid=isinside(x1,y1,z1,ix,iy,iz)
          if(ix.eq.ix0.and.iy.eq.iy0.and.iz.eq.iz0) then
            den = 1.0/vol
          else
           goto 100
          endif

        else


        ux=0d0
        uy=0d0
        uz=0d0
        gam=1.0d0
        if(mstc(146).eq.2) then
          uz=pv(3)/sqrt(pv(4)**2-pv(3)**2)
          gam=pv(4)/sqrt(pv(4)**2-pv(3)**2)
        else if(mstc(146).eq.3) then
          emd=pv(4)**2-pv(1)**2-pv(2)**2-pv(3)**2
          if(emd.le.1d-5) goto 100
          emd=sqrt(emd)
          ux=pv(1)/emd
          uy=pv(2)/emd
          uz=pv(3)/emd
          gam=pv(4)/emd
       endif
       enorm=widcof*gam

        ioptg1=0
       if(ioptg1.eq.1) then
          isfluid=isinside(x1,y1,z1,ix,iy,iz)
          if(isfluid.eq.0) cycle
           jx1=max(0,ix-ip); jx2=min(mstd(143),ix+ip)
           jy1=max(0,iy-iq); jy2=min(mstd(144),iy+iq)
           jz1=max(0,iz-ir); jz2=min(mstd(145),iz+ir)
           enorm=0d0
           do jx=jx1,jx2
           do jy=jy1,jy2
           do jz=jz1,jz2
c           x3=dx*( jx - mstd(143)/2 )
c           y3=dy*( jy - mstd(144)/2 )
c           z3=dz*( jz - mstd(145)/2 )
c           ex = (x1-x3)**2 + (y1-y3)**2 + (z1-z3)**2 
c    &      + ((x1-x3)*ux+(y1-y3)*uy+(z1-z3)*uz)**2
c           enorm = enorm + exp( -ex/widg )*vol
            enorm=enorm+gaussinter(jx,jy,jz,x1,y1,z1,ux,uy,uz)
           end do
           end do
           end do

           if(enorm.eq.0d0) then
             cycle
             print *,'jamdens enorm=0?',enorm,jx1,jx2,jy1,jy2,jz1,jz2
             print *,'x1=',x1,y1,z1
           endif

           enorm=1.0d0/enorm

           endif

        ioptg=2
        if(ioptg.eq.1) then
          den= gaussinter(ix0,iy0,iz0,x1,y1,z1,ux,uy,uz)/vol

        else if(ioptg.eq.2) then
        den=0.0d0
        do ig=1,3
        do jg=1,3
        do kg=1,3
          x2=x0 + xg3(ig)*dx/2
          y2=y0 + xg3(jg)*dy/2
          z2=z0 + xg3(kg)*dz/2
          x3= x1 - x2
          y3= y1 - y2
          z3= z1 - z2
          xtra=x3**2 + y3**2 + z3**2+(x3*ux + y3*uy + z3*uz)**2
         if(xtra/widg.gt.30.d0) cycle
         den=den+exp(-xtra/widg)*wg3(ig)*wg3(jg)*wg3(kg)/8
        end do
        end do
        end do

        else 
         xtra=(x1-x0)**2+(y1-y0)**2+(z1-z0)**2
     &       +((x1-x0)*ux + (y1-y0)*uy + (z1-z0)*uz)**2
         if(xtra/widg.gt.30.d0) cycle
         den=exp(-xtra/widg)
        endif

        endif

        bden=enorm*den*bar
        den=enorm*den*facq


c...Compute current and energy-momentum tensor.
        do im=1,4
          cur(im)  = cur(im)  + pv(im)/pv(4)*den
          curb(im) = curb(im) + pv(im)/pv(4)*bden
        do ik=1,4
          tnsmn = pv(im)*pv(ik)/pv(4)
          tens(im,ik) = tens(im,ik)  + tnsmn*den
        end do
        end do
        ncount=ncount+1

100   end do   ! End loop over particles.

      if(ncount.eq.0) return
      if(cur(4).le.1d-7) return

      cc=cur(4)**2-(cur(1)**2+cur(2)**2+cur(3)**2)
      if(cc.le.1d-7) return
      cc=sqrt(cc)

c...Compute local energy density and baryon density from EoS.
      if(mstc(137).eq.1) then
        uu(0)=tens(4,4)
        uu(1)=tens(1,4)
        uu(2)=tens(2,4)
        uu(3)=tens(3,4)
        uu(4)=curb(4)
        call thermal(uu,einv,rhob,pre,sva,vel)
        tau=1.0d0
        fg = uu(0) + tau*pre
        u(1) = uu(1)/fg
        u(2) = uu(2)/fg
        u(3) = uu(3)/fg
        u(4)=sqrt(1d0-u(1)**2-u(2)**2-u(3)**2)
        icon=0
        return
      endif


      u(1)=cur(1)/cc
      u(2)=cur(2)/cc
      u(3)=cur(3)/cc
      u(4)=cur(4)/cc

c     if(mstc(47).eq.2) call landauframe(u,tens)

c...Lorentz invariant Scalar Number Density
      rho=cur(4)*u(4)-cur(1)*u(1)-cur(2)*u(2)-cur(3)*u(3)

c...Lorentz invariant Baryon density
      rhob=curb(4)*u(4)-curb(1)*u(1)-curb(2)*u(2)-curb(3)*u(3)

c...Lorentz invariant pressure and energy density
      gam=u(4)
      vz=u(3)/u(4)
      gamz=1.0d0/sqrt(1.0d0-vz**2)
      z(4)=gamz*vz
      z(3)=gamz
      z(2)=0d0
      z(1)=0d0
      do i=1,4
      do j=1,4
       tmp=gg(i,j)*u(i)*u(j)
       einv=einv + tens(i,j)*tmp                    ! energy density
       pfree=pfree - 1.d0/3.d0*tens(i,j)*(g(i,j)-tmp) ! pressure

       pl=gg(i,j)*z(i)*z(j)
       plong=plong+tens(i,j)*pl
       pperp=pperp-0.5d0*tens(i,j)*(g(i,j)-tmp+pl)

c...Txx,Tyy,Tzz at the local rest frame.
       ci=-1.0d0
       cj=-1.0d0
       if(i.le.3) ci=u(i)/(1+gam)
       if(j.le.3) cj=u(j)/(1+gam)
       pxx=pxx+tens(i,j)*(dl(1,i)+u(1)*ci)*(dl(1,j)+u(1)*cj)
       pyy=pyy+tens(i,j)*(dl(2,i)+u(2)*ci)*(dl(2,j)+u(2)*cj)
       pzz=pzz+tens(i,j)*(dl(3,i)+u(3)*ci)*(dl(3,j)+u(3)*cj)

      end do
      end do

      icon=0

      end

c***********************************************************************
      subroutine jamdens(icon,i1,cc,rho,rhob,einv,u,iopt)

c...Purpose: compute energy-momentum tensor and hydrodynamics velocity.
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      dimension pv(4)
      dimension g(4,4),cur(4),curb(4),tens(4,4),u(4)
      dimension gg(4,4),dl(4,4),z(4),uu(0:4)
      data g/ -1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,1/
      data gg/ 1,1,1,-1,  1,1,1,-1, 1,1,1,-1, -1,-1,-1,1/
      data dl/ 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1/
      real*8 xg3(3),wg3(3)
      data xg3/0.0d0, -.77459666924148337703,.77459666924148337703/
      data wg3/.8888888888888888,.55555555555555,.555555555555555/
      common/jamdattens/cur,curb,tens
      save/jamdattens/

c     xg3(1)=0
c     xg3(2)=-sqrt(3.0d0/5d0)
c     xg3(3)=sqrt(3.0d0/5d0)
c     wg3(1)=8.0d0/9.0d0
c     wg3(2)=5.0d0/9.0d0
c     wg3(3)=5.0d0/9.0d0

      widg=2*parc(84)**2 ! Gaussian width
c     widcof=(1.d0/(3.14159d0*widg))**1.5d0)
      widcof=pard(90)

      icon=1
      einv=0.0d0
      rho=0d0
      rhob=0d0
      pfree=0.0d0
      pxx=0.0d0
      pyy=0.0d0
      pzz=0.0d0
      pperp=0d0
      plong=0d0
      ctime=pard(1)

      dcut=parc(147)
      dx=pard(143)
      dy=pard(144)
      dz=pard(145)
      ip=int(dcut/dx+0.5)
      iq=int(dcut/dy+0.5)
      ir=int(dcut/dz+0.5)

      if(i1.ge.1) then
        dt=ctime-r(4,i1)
        rx=r(1,i1)
        ry=r(2,i1)
        rz=r(3,i1)

        if(k(1,i1).eq.3) then

        e0=vq(4,i1)+vq(9,i1)
        vx=(vq(1,i1)+vq(6,i1))/e0
        vy=(vq(2,i1)+vq(7,i1))/e0
        vz=(vq(3,i1)+vq(8,i1))/e0

        else if(k(1,i1).eq.4) then
c        print *,',jamdens parton ?',k(1,i1),k(10,i1),k(11,i1)
c        print *,' kf=',k(2,i1)
          rx=0.0d0
          ry=0.0d0
          rz=0.0d0
          vx=0.0d0
          vy=0.0d0
          vz=0.0d0
          e0=0.0d0
          mpa=k(11,i1)-k(10,i1)+1
          do i=k(10,i1),k(11,i1)
            rx=rx+r(1,i)/mpa
            ry=ry+r(2,i)/mpa
            rz=rz+r(3,i)/mpa
            vx=vx+p(1,i)
            vy=vy+p(2,i)
            vz=vz+p(3,i)
            e0=e0+p(4,i)
          end do
          vx=vx/e0
          vy=vy/e0
          vz=vz/e0

        else
         vx=p(1,i1)/p(4,i1)
         vy=p(2,i1)/p(4,i1)
         vz=p(3,i1)/p(4,i1)
        endif

        x0=rx+dt*vx
        y0=ry+dt*vy
        z0=rz+dt*vz

      else
       x0=0d0; y0=0d0; z0=0d0
      endif

      do i=1,4
        cur(i)=0.0d0
        curb(i)=0.0d0
        do j=1,4
          tens(i,j)=0.0d0
      end do
      end do

      isfluid=isinside(x0,y0,z0,ix0,iy0,iz0)
      vol=pard(143)*pard(144)*pard(145)

      ncount=0


c....Loop over all particles
      do i=1,nv
 
        if(i.eq.i1) goto 100
        if(k(1,i).gt.10.or.k(1,i).eq.0) goto 100   ! dead particle
        if(p(5,i).le.1d-5) goto 100

c.....pre-formed hadrons.
        if(mstc(89).eq.1.and.r(5,i).gt.ctime) goto 100

        if(abs(k(7,i)).eq.1) goto 100   ! not yet interact

        dt=ctime-r(4,i)
        if(iopt.eq.1.and.dt.lt.0.0d0) goto 100

        if(k(2,i).eq.92) then
          do j=1,4
            pv(j)=vq(j,i)+vq(j+5,i)
          end do
        else
          pv(1)=p(1,i)
          pv(2)=p(2,i)
          pv(3)=p(3,i)
          pv(4)=p(4,i)
        endif

c       y1=0.5*log((pv(4)+pv(3))/(pv(4)-pv(3)))
c       if(abs(y1-ycm).gt.parc(30)) goto 100

        bar=k(9,i)/3.0d0
        facq=1.0d0
        facb=1.0d0
        if(iopt.eq.1.and.r(5,i).gt.ctime) then
          call qbfac(i,facq,facb)
        endif

        x1=r(1,i) + dt*pv(1)/pv(4)
        y1=r(2,i) + dt*pv(2)/pv(4)
        z1=r(3,i) + dt*pv(3)/pv(4)

        if(mstc(146).eq.0) then

          isfluid=isinside(x1,y1,z1,ix,iy,iz)
          if(ix.eq.ix0.and.iy.eq.iy0.and.iz.eq.iz0) then
            den = 1.0/vol
          else
           goto 100
          endif

        else


        ux=0d0
        uy=0d0
        uz=0d0
        gam=1.0d0
        if(mstc(146).eq.2) then
          uz=pv(3)/sqrt(pv(4)**2-pv(3)**2)
          gam=pv(4)/sqrt(pv(4)**2-pv(3)**2)
        else if(mstc(146).eq.3) then
          emd=pv(4)**2-pv(1)**2-pv(2)**2-pv(3)**2
          if(emd.le.1d-5) goto 100
          emd=sqrt(emd)
          ux=pv(1)/emd
          uy=pv(2)/emd
          uz=pv(3)/emd
          gam=pv(4)/emd
       endif
       enorm=widcof*gam

        ioptg1=0
       if(ioptg1.eq.1) then
          isfluid=isinside(x1,y1,z1,ix,iy,iz)
          if(isfluid.eq.0) cycle
           jx1=max(0,ix-ip); jx2=min(mstd(143),ix+ip)
           jy1=max(0,iy-iq); jy2=min(mstd(144),iy+iq)
           jz1=max(0,iz-ir); jz2=min(mstd(145),iz+ir)
           enorm=0d0
           do jx=jx1,jx2
           do jy=jy1,jy2
           do jz=jz1,jz2
c           x3=dx*( jx - mstd(143)/2 )
c           y3=dy*( jy - mstd(144)/2 )
c           z3=dz*( jz - mstd(145)/2 )
c           ex = (x1-x3)**2 + (y1-y3)**2 + (z1-z3)**2 
c    &      + ((x1-x3)*ux+(y1-y3)*uy+(z1-z3)*uz)**2
c           enorm = enorm + exp( -ex/widg )*vol
            enorm=enorm+gaussinter(jx,jy,jz,x1,y1,z1,ux,uy,uz)
           end do
           end do
           end do

           if(enorm.eq.0d0) then
             cycle
             print *,'jamdens enorm=0?',enorm,jx1,jx2,jy1,jy2,jz1,jz2
             print *,'x1=',x1,y1,z1
           endif

           enorm=1.0d0/enorm

           endif

        ioptg=2
        if(ioptg.eq.1) then
          den= gaussinter(ix0,iy0,iz0,x1,y1,z1,ux,uy,uz)/vol

        else if(ioptg.eq.2) then
        den=0.0d0
        do ig=1,3
        do jg=1,3
        do kg=1,3
          x2=x0 + xg3(ig)*dx/2
          y2=y0 + xg3(jg)*dy/2
          z2=z0 + xg3(kg)*dz/2
          x3= x1 - x2
          y3= y1 - y2
          z3= z1 - z2
          xtra=x3**2 + y3**2 + z3**2+(x3*ux + y3*uy + z3*uz)**2
         if(xtra/widg.gt.30.d0) cycle
         den=den+exp(-xtra/widg)*wg3(ig)*wg3(jg)*wg3(kg)/8
        end do
        end do
        end do

        else 
         xtra=(x1-x0)**2+(y1-y0)**2+(z1-z0)**2
     &       +((x1-x0)*ux + (y1-y0)*uy + (z1-z0)*uz)**2
         if(xtra/widg.gt.30.d0) cycle
         den=exp(-xtra/widg)
        endif

        endif

        bden=enorm*den*bar
        den=enorm*den*facq


c...Compute current and energy-momentum tensor.
        do im=1,4
          cur(im)  = cur(im)  + pv(im)/pv(4)*den
          curb(im) = curb(im) + pv(im)/pv(4)*bden
        do ik=1,4
          tnsmn = pv(im)*pv(ik)/pv(4)
          tens(im,ik) = tens(im,ik)  + tnsmn*den
        end do
        end do
        ncount=ncount+1

100   end do   ! End loop over particles.

      if(ncount.eq.0) return
      if(cur(4).le.1d-7) return

      cc=cur(4)**2-(cur(1)**2+cur(2)**2+cur(3)**2)
      if(cc.le.1d-7) return
      cc=sqrt(cc)

c...Compute local energy density and baryon density from EoS.
      if(mstc(137).eq.1) then
        uu(0)=tens(4,4)
        uu(1)=tens(1,4)
        uu(2)=tens(2,4)
        uu(3)=tens(3,4)
        uu(4)=curb(4)
        call thermal(uu,einv,rhob,pre,sva,vel)
        tau=1.0d0
        fg = uu(0) + tau*pre
        u(1) = uu(1)/fg
        u(2) = uu(2)/fg
        u(3) = uu(3)/fg
        u(4)=sqrt(1d0-u(1)**2-u(2)**2-u(3)**2)
        icon=0
        return
      endif


      u(1)=cur(1)/cc
      u(2)=cur(2)/cc
      u(3)=cur(3)/cc
      u(4)=cur(4)/cc

c     if(mstc(47).eq.2) call landauframe(u,tens)

c...Lorentz invariant Scalar Number Density
      rho=cur(4)*u(4)-cur(1)*u(1)-cur(2)*u(2)-cur(3)*u(3)

c...Lorentz invariant Baryon density
      rhob=curb(4)*u(4)-curb(1)*u(1)-curb(2)*u(2)-curb(3)*u(3)

c...Lorentz invariant pressure and energy density
      gam=u(4)
      vz=u(3)/u(4)
      gamz=1.0d0/sqrt(1.0d0-vz**2)
      z(4)=gamz*vz
      z(3)=gamz
      z(2)=0d0
      z(1)=0d0
      do i=1,4
      do j=1,4
       tmp=gg(i,j)*u(i)*u(j)
       einv=einv + tens(i,j)*tmp                    ! energy density
       pfree=pfree - 1.d0/3.d0*tens(i,j)*(g(i,j)-tmp) ! pressure

       pl=gg(i,j)*z(i)*z(j)
       plong=plong+tens(i,j)*pl
       pperp=pperp-0.5d0*tens(i,j)*(g(i,j)-tmp+pl)

c...Txx,Tyy,Tzz at the local rest frame.
       ci=-1.0d0
       cj=-1.0d0
       if(i.le.3) ci=u(i)/(1+gam)
       if(j.le.3) cj=u(j)/(1+gam)
       pxx=pxx+tens(i,j)*(dl(1,i)+u(1)*ci)*(dl(1,j)+u(1)*cj)
       pyy=pyy+tens(i,j)*(dl(2,i)+u(2)*ci)*(dl(2,j)+u(2)*cj)
       pzz=pzz+tens(i,j)*(dl(3,i)+u(3)*ci)*(dl(3,j)+u(3)*cj)

      end do
      end do

      icon=0

      end


c***************************************************************************
      subroutine hydro2jamall(it,icheck)

c....Convert fluid element to jam particle.
c...icheck=0: final call, i.e. all fluid elements are below
c...particlization energy density e_p.
c...icheck>0: some of fluid elements are below e_p.

      implicit none
      include 'fluid.inc'
      include 'jam2.inc'
      integer it,icon,npos,ntot,nbar,nch,nstr,icheck
      real*8 ptot(5),etot,dtot,pet,pxt,pyt,pzt
      real*8 cc,rho,rhob,einv,uv(4)
      integer opt_mode,opt_econ,n,ix,iy,iz

      mste(43)=icheck

c...Option for sampling. =1: only total momentum 
c...=2: energy and momentum conservation are recovered
      opt_econ=mstc(149)

c...=0: Conventional sampling (no conservation of quantum number)
c...=1: Mode sampling
c...=2: Single Particle Rejection with Exponential Weights (SPREW) sampling
      opt_mode=mstc(148)

      if(mstc(8).ge.1) print *,'hydro2jamall',nfreezeout

c...Compute total energy-momentum that will be converted into particles.
      fconv=0.0d0

c     if(icheck.ne.0) then
c     do n=1,nfreezeout
c       ix=frr(1,n)
c       iy=frr(2,n)
c       iz=frr(3,n)
c       fconv =fconv+u(:,ix,iy,iz)*dx*dy*dz*hbc
c     end do
c     endif

      if(nfreezeout.le.0) then
        if(icheck.ne.0) return
c       print *,'no freezeout surface',nfreezeout
c       print *,'pe=',pe_tot,'p=',px_tot,py_tot,pz_tot
c       print *,'b_tot=',b_tot,c_tot,s_tot
        nfreezeout=0
        call isochronous_freezeout(pard(1),1)
        if(nfreezeout.eq.0) then
          print *,'still nfreezeout=0?'
          return
        endif
c       print *,'nfreeze=',nfreezeout
        call sample1p(npos,nbar,nch,nstr,ptot,etot,dtot,ntot)
        goto 1200
      endif


      if(opt_mode.le.9) then
        call sampling(npos,nbar,nch,nstr,ptot,etot,dtot,ntot)
      else if(opt_mode.eq.10) then
        call modesampling(npos,nbar,nch,nstr,ptot,etot,dtot,ntot)
      else if(opt_mode.lt.20) then
        call sampling2(npos,nbar,nch,nstr,ptot,etot,dtot,ntot)

c...Sample particle according to the positive part of Cooper-Frye
c...using Murase's integral formula.
      else if(opt_mode.lt.30) then
        call sampling2pos(npos,nbar,nch,nstr,ptot,etot,dtot,ntot)

      else if(opt_mode.lt.40) then
        call sampling3(npos,nbar,nch,nstr,ptot,etot,dtot,ntot)

c...Using Scott Pratt method.
      else if(opt_mode.lt.50) then ! SP1
        call sampling4(npos,nbar,nch,nstr,ptot,etot,dtot,ntot)
      else if(opt_mode.lt.60) then ! SP2
        call sampling4b(npos,nbar,nch,nstr,ptot,etot,dtot,ntot)
      else if(opt_mode.lt.70) then  ! SP3
        call  sampling4c(npos,nbar,nch,nstr,ptot,etot,dtot,ntot)
      else if(opt_mode.lt.80) then  ! SP4 modesampling
        call modesampling_pos(npos,nbar,nch,nstr,ptot,etot,dtot,ntot)
      else
        write(6,*)'wrong option mstc(148)=',mstc(148)
        stop
      endif


1200  continue
      if(ntot.eq.0) then
        if(dtot.ge.1.0) then
          print *,'hydro2jam no particle ntot',ntot,dtot
          return
        endif
        call sample1p(npos,nbar,nch,nstr,ptot,etot,dtot,ntot)
        if(ntot.ge.1) goto 1500
        if(ntot.eq.0.and.icheck.ne.0) return

        nfreezeout=0
        call isochronous_freezeout(pard(1),1)
        if(nfreezeout.eq.0) then
          print *,'still nfreezeout=0 at the end of hydro?',
     &    pe_tot,px_tot,py_tot,pz_tot
          return
        endif
        print *,'nfreeze=',nfreezeout
        call sample1p(npos,nbar,nch,nstr,ptot,etot,dtot,ntot)
        if(ntot.eq.0) then
           print *,'no particle after hydro?',ntot
           print *,'pe_tot=',pe_tot,px_tot,py_tot,pz_tot
           return
        endif

      endif

1500  continue
c     if(mstc(144).le.2) then
      if(icheck.eq.0) then
        pet=pe_tot
        pxt=px_tot
        pyt=py_tot
        pzt=pz_tot
      else if(mstc(144).eq.11) then
        fconv=0.0
        do n=1,nfreezeout
          fconv=fconv+u(:,frr(1,n),frr(2,n),frr(3,n))*dx*dy*dz*hbc
          call reset_fluid_element(frr(1,n),frr(2,n),frr(3,n))
        end do
        pet=fconv(0)
        pxt=fconv(1)
        pyt=fconv(2)
        pzt=fconv(3)
      else
        pet=fconv(0)
        pxt=fconv(1)
        pyt=fconv(2)
        pzt=fconv(3)
        nbar = fconv(4)
c       pet=ptot(4)
c       pxt=ptot(1)
c       pyt=ptot(2)
c       pzt=ptot(3)
c       call edit_frezeout_surface
      endif

      if(opt_econ.ge.1) then
        call recover_mom(opt_econ,ntot,ptot, pet,pxt,pyt,pzt)
      endif

c...Dynamical freeze-out.
c     if(mstc(144).ge.10) then
        px_tot = px_tot - pxt
        py_tot = py_tot - pyt
        pz_tot = pz_tot - pzt
        pe_tot = pe_tot - pet
        b_tot= b_tot - nbar
        c_tot= c_tot - nch
        s_tot= s_tot - nstr
c     endif

c...Reset freeze-out vectors.
c     if(mstc(144).ne.12) then

      nfreezeout=0
      if(mstc(144).ge.1) call c_freezeout(1,0,0d0)

c     endif

c------------------------------------------------------------------------


      if(mstc(8).ge.1) then

c     pftot=0.0
c     do ix = 0,maxx
c     do iy = 0,maxy
c     do iz = 0,maxz
c       pftot = pftot + u(:,ix,iy,iz)*dx*dy*dz*hbc
c     end do
c     end do
c     end do
c     print *,'pftot0-pftot',pftot0-pftot

      print *,'number of hyper surface',nfreezeout, 'npos=',npos
      print *,'number of sampled particle',ntot,dtot

      print *,'total energy=',etot,ptot(4),pet
c    & ,(etot-pet)/pe_tot*100,'%',
c    & (ptot(4)-pet)/pet*100,'%'

      print *,'baryon number=',nbar,b_tot,
     &  abs(nbar-b_tot)/max(1d0,b_tot)*100,'%'

      print *,'momentum f=',pet,pxt,pyt,pzt
      print *,'momentum p=',ptot(4),ptot(1),ptot(2),ptot(3)
      if(abs(pxt-ptot(1)).gt.1d-5.or.
     &   abs(pyt-ptot(2)).gt.1d-5.or.
     &   abs(pzt-ptot(3)).gt.1d-5) then
        print *,'momentum strange'
        stop
      endif

      print *,'total strangeness=',nstr,s_tot,
     & abs(nstr-s_tot)/max(1d0,s_tot)*100,'%'
      print *,'total charge=',nch,c_tot,
     & abs(nch-c_tot)/max(1d0,c_tot)*100,'%'

c...energy density at (0,0,0)
      call jamdens(icon,0,cc,rho,rhob,einv,uv,0)
      print *,'icon=',icon
      print *,'baryon density=',rhob,' fluid=',bl(origx,origy,origz)
     &,u(4,origx,origy,origz)
      print *,'enegy density=',einv,' fluid=',el(origx,origy,origz)*hbc
     &,u(0,origx,origy,origz)*hbc


c     call hydroTotalEnergy(it,0.0d0,pet,pxt,pyt,pzt,nbar)

      endif


c------------------------------------------------------------------------
c...Output energy-momentum, baryon number, charge, strangeness.
      if(mstc(180).eq.1) then
      write(80,800)ptot(4)-pet,
     & ptot(1)-pxt,ptot(2)-pyt,ptot(3)-pzt,
     &  nbar-nint(b_tot),nch-nint(c_tot),nstr-nint(s_tot)
      flush(80)
 800  format(4(e14.7,1x),1x,i7,1x,i6,1x,i4)
      endif
c------------------------------------------------------------------------


      end

c***************************************************************************
      subroutine sampling(npos,nbar,nch,nstr,ptot,etot,dtot,ntot)

c...Particles are sampled by checking all hyper-surface.
      implicit none
      integer i,n,npos,npar,nbar,nch,nstr,kc,kf,ibary,str,cha,ntot
      integer mode,mode1,mode2,nmode(0:7),istat
      real*8 rn,v1,v2,v3,dst,dsx,dsy,dsz,gam,vol,tdns,pdns
      real*8 tf,bmu,smu,cmu,ene,etot,dtot,rc(5),pc(5),ptot(5)
      real*8 multlookup,mult,pjmass
      integer generate_poisson,generate_poisson1
      real*8 sp_potential
      include 'fluid.inc'
      include 'jam2.inc'

      npos=0
      nbar=0
      nch=0
      nstr=0
      ptot=0d0
      etot=0d0
      dtot=0d0
      ntot=0

      if(nfreezeout.le.0) return

      do 1000 n=1,nfreezeout     

      v1=frv(1,n)
      v2=frv(2,n)
      v3=frv(3,n)
      dst=frdsgm(0,n)
      dsx=frdsgm(1,n)
      dsy=frdsgm(2,n)
      dsz=frdsgm(3,n)
      gam=1.0d0/sqrt(1.0d0-v1*v1-v2*v2-v3*v3)
      vol=gam*(dst + dsx*v1 + dsy*v2 + dsz*v3)
      if(vol.le.0d0) goto 1000

      tf=frtmp(n)*hbc
      bmu=frmuq(n)*hbc
      smu=frmus(n)*hbc
      pare(51)=tf
      pare(52)=bmu
      pare(53)=smu
c...2018/10/19 bmu in the new EoS table already contains the
c...potential,so we do not need to add potential anymore.
c     bmu = bmu - pard(142)*frbdn(n)
c     bmu = bmu - sp_potential(frbdn(n))

      if(opt_table.ge.1) then
        pdns=frnum(n)
      else
        pdns=mult(tf,bmu,smu)
      endif
      if(pdns.le.0d0) goto 1000
 
      dtot = dtot + pdns*vol*ntest
      npos=npos+1

c...Determine number of particle to be generated.
      npar=generate_poisson(pdns*vol)*ntest
      if(npar.eq.0) goto 1000

      if(opt_table.ge.1) pdns=mult(tf,bmu,smu)

      do i=1,npar
        call selectparticle(tf,bmu,smu,pdns,nstr,nbar,nch,kc,kf)
        ibary=kchg(kc,6)*isign(1,kf)
        str=kchg(kc,7)*isign(1,kf)
        cha=kchg(kc,1)*isign(1,kf)/3
        istat=1
        if(ibary.eq.0) istat=-1
        cmu = str*smu + ibary*bmu/3
        pc(5)=pjmass(kf)
        call cooperfrye(dst,dsx,dsy,dsz,v1,v2,v3,tf,cmu,istat,pc)

c...Set coordinate of particle.
        rc(1)=dx*( frr(1,n) - origx - 0.5 + rn(0) )
        rc(2)=dy*( frr(2,n) - origy - 0.5 + rn(0) )
        rc(3)=dz*( frr(3,n) - origz - 0.5 + rn(0) )

c...Set time and formation time of particle.
        rc(4)=frt(n)
        rc(5)=rc(4)
        call p2jam(kc,kf,ibary,rc,pc)
c       pe_tot=pe_tot-pc(4)
c       px_tot=px_tot-pc(1)
c       py_tot=py_tot-pc(2)
c       pz_tot=pz_tot-pc(3)
c       b_tot=b_tot-ibary
        ptot=ptot+pc

c....Reset fluid elements.
        call reset_fluid_element(frr(1,n),frr(2,n),frr(3,n))

        ntot = ntot + 1
        nbar = nbar + ibary/3/ntest
        nstr = nstr + str/ntest
        nch = nch + cha/ntest
        etot = etot + pc(4)/ntest

      end do

 1000 continue

      end

c***************************************************************************
      subroutine sampling4(npos,nbar,nch,nstr,ptot,etot,nftot,ntot)

c...MC sampling of particles from freeze-out surface
c...by the Scott Pratt method I.
      implicit none
      integer i,n,npos,npar,nbar,nch,nstr,kc,kf,ibary,str,cha,ntot
      integer mode,mode1,mode2,nmode(0:7),istat
      real*8 rn,v1,v2,v3,dst,dsx,dsy,dsz,gam,vol,tdns,pdns,nftot
      real*8 tf,bmu,smu,cmu,ene,etot,rc(5),ptot(5)
      real*8 cost,sint,phi,er,pc(5),thermaldist3,pds,pp,prob
      real*8 pjmass,vopt,spin,fac,pdensity,pdensity2,tdns2
      real*8 probb,probs,probc,db,ds,dc
      real*8 sp_potential
      integer j,k,opt,sgn
      integer generate_poisson
      include 'fluid.inc'
      include 'jam2.inc'

      fac=1.0d0/(2*paru(1)**2*paru(3)**3)
      npos=0
      nbar=0
      nch=0
      nstr=0
      ptot=0d0
      etot=0d0
      ntot=0
      opt=mod(mstc(148),10)

      do n=1,nfreezeout

      v1=frv(1,n)
      v2=frv(2,n)
      v3=frv(3,n)
      dst=frdsgm(0,n)
      dsx=frdsgm(1,n)
      dsy=frdsgm(2,n)
      dsz=frdsgm(3,n)
      gam=1.0d0/sqrt(1.0d0-v1*v1-v2*v2-v3*v3)
      vol=gam*(dst + dsx*v1 + dsy*v2 + dsz*v3)
      vopt=vol+sqrt(vol**2 - (dst**2 - dsx**2-dsy**2-dsz**2))
      tf=frtmp(n)*hbc
      bmu=frmuq(n)*hbc
      smu=frmus(n)*hbc
      pare(51)=tf
      pare(52)=bmu
      pare(53)=smu

c...Include potential.
c     bmu = bmu - pard(142)*frbdn(n)
c     bmu = bmu - sp_potential(frbdn(n))

c...Loop over particle species.
      do kc=132,366
        if(kc.ge.195.and.kc.le.251) cycle ! exclude charm
        if(kc.eq.171.or.kc.eq.172) cycle ! K_L0 or K_S0
      do k=1,2  ! particle-anti-particle loop
        if(k.eq.1) sgn=1.0d0
        if(k.eq.2) then
          if(kchg(kc,3).eq.0) cycle
          sgn=-1.0d0
        endif
        ibary=kchg(kc,6)*sgn
        str=kchg(kc,7)*sgn
        cha=kchg(kc,1)*sgn/3
        if(opt.eq.2) then
          db=nbar-b_tot
          ds=nstr-s_tot
          dc=nch-c_tot
          probb=exp(-min(50d0,abs(db)))
          probs=exp(-min(50d0,abs(ds)))
          probc=exp(-min(50d0,abs(dc)))
          if(ibary*db.gt.0d0.and.rn(0).gt.probb) cycle
          if(str*ds.gt.0d0.and.rn(0).gt.probs) cycle
          if(cha*dc.gt.0d0.and.rn(0).gt.probc) cycle
        endif
        kf=kchg(kc,4)
        spin=max(1,mod(kf,10))
        kf=kf*sgn
        istat=1
        if(ibary.eq.0) istat=-1
        cmu = str*smu + ibary*bmu/3
        pc(5)=pmas(kc,1)
        tdns=vopt*pdensity(pc(5),tf,cmu,ibary,istat)*spin*fac
c       tdns2=vopt*pdensity2(pc(5),tf,cmu,istat)*spin*fac

        npar=generate_poisson(tdns)*ntest
        if(npar.le.0) cycle

c...Sample momentum.
        do i=1,npar
        pc(5)=pjmass(kf)
100     pp=thermaldist3(pc(5),tf,cmu,istat)
        cost=2d0*rn(0)-1d0
        sint=sqrt(1d0-cost**2)
        phi=paru(2)*rn(0)
        pc(1)=pp*sint*cos(phi)
        pc(2)=pp*sint*sin(phi)
        pc(3)=pp*cost
        pc(4)=sqrt(pp**2+pc(5)**2)
        er=pc(4)
        call jamrobo(0d0,0d0,v1,v2,v3,gam,pc(1),pc(2),pc(3),pc(4))
        pds = pc(4)*dst+pc(1)*dsx+pc(2)*dsy+pc(3)*dsz
        if(pds.le.0d0) cycle
        prob=pds/(er*vopt)
        if(prob.lt.0.0.or.prob.gt.1.0) then
          print *,'prob?',prob,vopt
          vopt=vopt*1.2d0
          goto 100
        endif
        if(rn(0).gt.prob) cycle

c...Set coordinate of particle.
        rc(1)=dx*( frr(1,n) - origx - 0.5 + rn(0) )
        rc(2)=dy*( frr(2,n) - origy - 0.5 + rn(0) )
        rc(3)=dz*( frr(3,n) - origz - 0.5 + rn(0) )

c...Set time and formation time of particle.
        rc(4)=frt(n)
        rc(5)=rc(4)
        call p2jam(kc,kf,ibary,rc,pc)
c       pe_tot=pe_tot-pc(4)
c       px_tot=px_tot-pc(1)
c       py_tot=py_tot-pc(2)
c       pz_tot=pz_tot-pc(3)
c       b_tot=b_tot-ibary
        ptot=ptot+pc

c....Reset fluid elements.
c       call reset_fluid_element(frr(1,n),frr(2,n),frr(3,n))

        npos=npos+1
        ntot = ntot + 1
        nftot=ntot
        nbar = nbar + ibary/3/ntest
        nstr = nstr + str/ntest
        nch = nch + cha/ntest
        etot = etot + pc(4)/ntest
        end do

      end do
      end do

      end do

      end


c***************************************************************************
      subroutine setlocalq(ix,iy,iz)
      implicit none
      include 'fluid.inc'
      integer ix,iy,iz
      double precision tau,ene,nba,pre,sav,vel
      double precision tplc,mulc,smu,slc,fg,sva,va

      call thermal(u(:,ix,iy,iz),ene,nba,pre,sva,vel)
      call geteos2(ene,nba,tplc,mulc,smu,slc)

      tau=1.0
      fg = u(0,ix,iy,iz) + tau*pre
      vx(ix,iy,iz) = u(1,ix,iy,iz)/fg
      vy(ix,iy,iz) = u(2,ix,iy,iz)/fg
      vz(ix,iy,iz) = u(3,ix,iy,iz)/fg

      va=sqrt(vx(ix,iy,iz)**2 + vy(ix,iy,iz)**2 + vz(ix,iy,iz)**2)
      if(va.ge.1.0d0) then
            print *,'setlocalq v>=1?',va,fg,vel
            print *,'vx=',vx(ix,iy,iz),vy(ix,iy,iz),vz(ix,iy,iz)
      endif

      cs(ix,iy,iz)=sva
      el(ix,iy,iz) = ene
      bl(ix,iy,iz) = nba
      pl(ix,iy,iz) = pre
      tl(ix,iy,iz) = tplc
      ml(ix,iy,iz) = mulc

      end

c***************************************************************************
      subroutine sample1p(npos,nbar,nch,nstr,ptot,etot,nftot,ntot)

c...Sample one particle from freeze-out surface.
      implicit none
      integer i,n,npos,npar,nbar,nch,nstr,kc,kf,ibary,str,cha,ntot
      integer mode,mode1,mode2,nmode(0:7),istat
      real*8 rn,v1,v2,v3,dst,dsx,dsy,dsz,gam,vol,tdns,pdns,nftot
      real*8 tf,bmu,smu,cmu,ene,etot,rc(5),ptot(5)
      real*8 cost,sint,phi,er,pc(5),thermaldist3,pds,pp,prob
      real*8 pjmass,vopt,spin,pdensity
      real*8 mult,probb,probs,probc,db,ds,dc
      real*8 sp_potential,utmp(0:4),det
      integer ix,iy,iz,it
      integer j,k,opt,ntry,itry
      integer generate_poisson
      include 'fluid.inc'
      include 'jam2.inc'

      ntry=0
1000  continue
      npos=0
      nbar=0
      nch=0
      nstr=0
      ptot=0d0
      etot=0d0
      ntot=0
      nftot=0.0
      opt=mod(mstc(148),10)
c...Chose freeze-out hyper-surface.
200   n=int(1+(nfreezeout-1)*rn(0))
      ntry=ntry+1
      if(ntry.ge.100) then
c       print *,'sample1p no positive freezeout surpface?',nfreezeout
        return
      endif

      if(frtmp(n).lt.tfreezecut) goto 200
      v1=frv(1,n)
      v2=frv(2,n)
      v3=frv(3,n)
      dst=frdsgm(0,n)
      dsx=frdsgm(1,n)
      dsy=frdsgm(2,n)
      dsz=frdsgm(3,n)
      gam=1.0d0/sqrt(1.0d0-v1*v1-v2*v2-v3*v3)
      vol=gam*(dst + dsx*v1 + dsy*v2 + dsz*v3)
      vopt=vol+sqrt(vol**2 - (dst**2 - dsx**2-dsy**2-dsz**2))
      if(vopt.le.0.0) goto 200
      tf=frtmp(n)*hbc
      bmu=frmuq(n)*hbc
      smu=frmus(n)*hbc
      pare(51)=tf
      pare(52)=bmu
      pare(53)=smu
c     print *,'tf=',tf,'mu=',bmu
      nftot = nftot + frnum(n)*vol*ntest
      npar=1
      pdns=mult(tf,bmu,smu)


        do i=1,npar
        call selectparticle(tf,bmu,smu,pdns,nstr,nbar,nch,kc,kf)
        ibary=kchg(kc,6)*isign(1,kf)
        str=kchg(kc,7)*isign(1,kf)
        cha=kchg(kc,1)*isign(1,kf)/3
        istat=1
        if(ibary.eq.0) istat=-1
        cmu = str*smu + ibary*bmu/3
        pc(5)=pjmass(kf)
        itry=0
100     pp=thermaldist3(pc(5),tf,cmu,istat)
        itry=itry+1
        if(itry.ge.100) then
          print *,'sample1p no phase space for partile?'
          return
        endif
        cost=2d0*rn(0)-1d0
        sint=sqrt(1d0-cost**2)
        phi=paru(2)*rn(0)
        pc(1)=pp*sint*cos(phi)
        pc(2)=pp*sint*sin(phi)
        pc(3)=pp*cost
        pc(4)=sqrt(pp**2+pc(5)**2)
        er=pc(4)
        call jamrobo(0d0,0d0,v1,v2,v3,gam,pc(1),pc(2),pc(3),pc(4))
        pds = pc(4)*dst+pc(1)*dsx+pc(2)*dsy+pc(3)*dsz
c       print *,'pds=',pds
        if(pds.le.0d0) goto 100

        prob=pds/(er*vopt)
c       print *,'prob=',prob
        if(rn(0).gt.prob) goto 100

c...Set coordinate of particle.
        rc(1)=dx*( frr(1,n) - origx - 0.5 + rn(0) )
        rc(2)=dy*( frr(2,n) - origy - 0.5 + rn(0) )
        rc(3)=dz*( frr(3,n) - origz - 0.5 + rn(0) )

c...Set time and formation time of particle.
        rc(4)=frt(n)
        rc(5)=rc(4)

c...Real-time freeze-out by the Cooper-Frey formula
        if(mstc(144).eq.12.and.mste(43).ne.0) then
          vol=dx*dy*dz*hbc
          ix=frr(1,n)
          iy=frr(2,n)
          iz=frr(3,n)
          utmp(0)=u(0,ix,iy,iz)-pc(4)/vol
          utmp(1)=u(1,ix,iy,iz)-pc(1)/vol
          utmp(2)=u(2,ix,iy,iz)-pc(2)/vol
          utmp(3)=u(3,ix,iy,iz)-pc(3)/vol
          det=utmp(0)**2-utmp(1)**2-utmp(2)**2-utmp(3)**2
          if(utmp(0).le.0.0d0.or.det.le.0d0) then

c           cycle  ! Skip this particle because of too much energy.

            fconv(:) = fconv(:) + u(:,ix,iy,iz)*vol
c           pc(1)=u(1,ix,iy,iz)*vol
c           pc(2)=u(2,ix,iy,iz)*vol
c           pc(3)=u(3,ix,iy,iz)*vol
c           pc(4)=sqrt(pc(5)**2+pc(1)**2+pc(2)**2+pc(3)**2)
            call reset_fluid_element(frr(1,n),frr(2,n),frr(3,n))
            frtmp(n)=0.0d0
          else
            fconv(:) = fconv(:) + u(:,ix,iy,iz)-utmp(:)
            u(:,ix,iy,iz)=utmp(:)
            call setlocalq(ix,iy,iz)
          endif
          fconv(4) = fconv(4) + ibary
        endif

        call p2jam(kc,kf,ibary,rc,pc)

c....Reset fluid elements.
c       call reset_fluid_element(frr(1,n),frr(2,n),frr(3,n))

c       print *,'Sample 1 particle!',n
c       print *,'u=',(u(k,frr(1,n),frr(2,n),frr(3,n)),k=1,3)
c       print *,'v=',v1,v2,v3
c       print *,'p=',(pc(k),k=1,3)
c       print *,'e=',u(0,frr(1,n),frr(2,n),frr(3,n)),pc(4)

        ptot = ptot + pc
        npos=npos+1
        ntot = ntot + 1
        nbar = nbar + ibary/3/ntest
        nstr = nstr + str/ntest
        nch = nch + cha/ntest
        etot = etot + pc(4)/ntest
        end do

      end

c***************************************************************************
      subroutine sampling4c(npos,nbar,nch,nstr,ptot,etot,nftot,ntot)

c...MC sampling of particles from freeze-out surface by the method of
c...Scott Pratt III.
      implicit none
      integer i,n,npos,npar,nbar,nch,nstr,kc,kf,ibary,str,cha,ntot
      integer mode,mode1,mode2,nmode(0:7),istat
      real*8 rn,v1,v2,v3,dst,dsx,dsy,dsz,gam,vol,volf,tdns,pdns,nftot
      real*8 tf,bmu,smu,cmu,ene,etot,rc(5),ptot(5)
      real*8 cost,sint,phi,er,pc(5),thermaldist3,pds,pp,prob
      real*8 pjmass,vopt,spin,pdensity
      real*8 mult,probb,probs,probc,db,ds,dc
      real*8 sp_potential,utmp(0:4),det
      integer j,k,opt,sgn,ntry
      integer ix,iy,iz
      integer generate_poisson
      include 'fluid.inc'
      include 'jam2.inc'

      volf=dx*dy*dz*hbc

      ntry=0
1000  continue
c     fac=1.0d0/(2*paru(1)**2*paru(3)**3)
      npos=0
      nbar=0
      nch=0
      nstr=0
      ptot=0d0
      etot=0d0
      ntot=0
      nftot=0.0
      opt=mod(mstc(148),10)
c...Loop over freeze-out hyper-surface.
      do n=1,nfreezeout

      if(frtmp(n).lt.tfreezecut) cycle
      v1=frv(1,n)
      v2=frv(2,n)
      v3=frv(3,n)
      dst=frdsgm(0,n)
      dsx=frdsgm(1,n)
      dsy=frdsgm(2,n)
      dsz=frdsgm(3,n)
      gam=1.0d0/sqrt(1.0d0-v1*v1-v2*v2-v3*v3)
      vol=gam*(dst + dsx*v1 + dsy*v2 + dsz*v3)
      vopt=vol+sqrt(vol**2 - (dst**2 - dsx**2-dsy**2-dsz**2))
      tf=frtmp(n)*hbc
      bmu=frmuq(n)*hbc
      smu=frmus(n)*hbc
      pare(51)=tf
      pare(52)=bmu
      pare(53)=smu
c...Include potential.
c     bmu = bmu - pard(142)*frbdn(n)
c     bmu = bmu - sp_potential(frbdn(n))

      nftot = nftot + frnum(n)*vol*ntest
c...Determine number of particle to be generated.
      if(mstc(145).ge.1) then
        pdns=frnum(n)
        if(pdns.le.0d0) cycle
        npar=generate_poisson(pdns*vopt)*ntest
        if(npar.eq.0) cycle
        pdns=mult(tf,bmu,smu)
      else
        pdns=mult(tf,bmu,smu)
        if(pdns.le.0d0) cycle
        npar=generate_poisson(pdns*vopt)*ntest
        if(npar.eq.0) cycle
      endif

        do i=1,npar
        call selectparticle(tf,bmu,smu,pdns,nstr,nbar,nch,kc,kf)
        ibary=kchg(kc,6)*isign(1,kf)
        str=kchg(kc,7)*isign(1,kf)
        cha=kchg(kc,1)*isign(1,kf)/3
        istat=1
        if(ibary.eq.0) istat=-1
        cmu = str*smu + ibary*bmu/3
        pc(5)=pjmass(kf)
100     pp=thermaldist3(pc(5),tf,cmu,istat)
        cost=2d0*rn(0)-1d0
        sint=sqrt(1d0-cost**2)
        phi=paru(2)*rn(0)
        pc(1)=pp*sint*cos(phi)
        pc(2)=pp*sint*sin(phi)
        pc(3)=pp*cost
        pc(4)=sqrt(pp**2+pc(5)**2)
        er=pc(4)
        call jamrobo(0d0,0d0,v1,v2,v3,gam,pc(1),pc(2),pc(3),pc(4))
        pds = pc(4)*dst+pc(1)*dsx+pc(2)*dsy+pc(3)*dsz
        if(pds.le.0d0) cycle

        prob=pds/(er*vopt)
        if(prob.lt.0.0.or.prob.gt.1.0) then
          print *,'prob?',prob,vopt
          vopt=vopt*1.2d0
          goto 100
        endif
        if(rn(0).gt.prob) cycle

c...Set coordinate of particle.
        rc(1)=dx*( frr(1,n) - origx - 0.5 + rn(0) )
        rc(2)=dy*( frr(2,n) - origy - 0.5 + rn(0) )
        rc(3)=dz*( frr(3,n) - origz - 0.5 + rn(0) )

c...Set time and formation time of particle.
        rc(4)=frt(n)
        rc(5)=rc(4)

c...Real-time freeze-out by the Cooper-Frye formula:remove
c...energy-momentum of particle from fluid element.
        if(mstc(144).eq.12.and.mste(43).ne.0) then
          ix=frr(1,n)
          iy=frr(2,n)
          iz=frr(3,n)
          utmp(0)=u(0,ix,iy,iz)-pc(4)/volf
          utmp(1)=u(1,ix,iy,iz)-pc(1)/volf
          utmp(2)=u(2,ix,iy,iz)-pc(2)/volf
          utmp(3)=u(3,ix,iy,iz)-pc(3)/volf
          det=utmp(0)**2-utmp(1)**2-utmp(2)**2-utmp(3)**2
          if(utmp(0).le.0.0d0.or.det.le.0d0) then
            cycle  ! Skip this particle because of too much energy.

c...Nevertheless let's include this particle setting fluid element to be zero
            fconv(:) = fconv(:) + u(:,ix,iy,iz)*volf
            call reset_fluid_element(frr(1,n),frr(2,n),frr(3,n))
            frtmp(n)=0.0d0

c...No problem. Subtract particle energy-momentum from the fluid element. 
          else
            fconv(:) = fconv(:) + u(:,ix,iy,iz)-utmp(:)
            u(:,ix,iy,iz)=utmp
            call setlocalq(ix,iy,iz)
          endif
          fconv(4) = fconv(4) + ibary
        endif

        call p2jam(kc,kf,ibary,rc,pc)

        ptot=ptot+pc

c....Reset fluid elements.
c       call reset_fluid_element(frr(1,n),frr(2,n),frr(3,n))

        npos=npos+1
        ntot = ntot + 1
c       nftot=ntot
        nbar = nbar + ibary/3/ntest
        nstr = nstr + str/ntest
        nch = nch + cha/ntest
        etot = etot + pc(4)/ntest
        end do

3000  end do


      ntry=ntry+1
      if(ntot.eq.0) then
        print *,'sampling4c ntot=0 nftot=',nftot,ntot,nfreezeout
        if(nftot.lt.1.0) return

        if(ntry.le.20) goto 1000
c...Create one particle.
c       if(nftot.lt.1.0) then
c         call sample1p(npos,nbar,nch,nstr,ptot,etot,nftot,ntot)
c         return
c       endif

        print *,'no particle sampling4c ? ',ntry,nftot,ntot
      endif

      end

c***************************************************************************
      subroutine sampling4b(npos,nbar,nch,nstr,ptot,etot,nftot,ntot)

c...MC sampling of particles from freeze-out surface
c...by the Scott Pratt method II.
      implicit none
      integer i,n,npos,npar,nbar,nch,nstr,kc,kf,ibary,str,cha,ntot
      integer mode,mode1,mode2,nmode(0:7),istat
      real*8 rn,v1,v2,v3,dst,dsx,dsy,dsz,gam,vol,tdns,pdns,nftot
      real*8 tf,bmu,smu,cmu,ene,etot,rc(5),ptot(5)
      real*8 cost,sint,phi,er,pc(5),thermaldist3,pds,pp,prob
      real*8 pjmass,vopt,spin,fac,pdensity,rfac
      real*8 probb,probs,probc,db,ds,dc,a
      integer j,k,opt,sgn,i0,ka
      real*8 sp_potential
      integer generate_poisson
      include 'fluid.inc'
      include 'jam2.inc'
      real*8 scsum, mult
      integer kcsum,maxkc
      common/bolz4/maxkc,kcsum(0:500),scsum(0:500)

      fac=1.0d0/(2*paru(1)**2*paru(3)**3)
      npos=0
      nbar=0
      nch=0
      nstr=0
      ptot=0d0
      etot=0d0
      ntot=0
      opt=mod(mstc(148),10)

      do n=1,nfreezeout

c...Determine number of particle to be generated.
      v1=frv(1,n)
      v2=frv(2,n)
      v3=frv(3,n)
      dst=frdsgm(0,n)
      dsx=frdsgm(1,n)
      dsy=frdsgm(2,n)
      dsz=frdsgm(3,n)
      gam=1.0d0/sqrt(1.0d0-v1*v1-v2*v2-v3*v3)
      vol=gam*(dst + dsx*v1 + dsy*v2 + dsz*v3)
      vopt=vol+sqrt(vol**2 - (dst**2 - dsx**2-dsy**2-dsz**2))
      tf=frtmp(n)*hbc
      bmu=frmuq(n)*hbc
      smu=frmus(n)*hbc
      pare(51)=tf
      pare(52)=bmu
      pare(53)=smu

c...Include potential.
c     bmu = bmu - pard(142)*frbdn(n)
c     bmu = bmu - sp_potential(frbdn(n))

      i0=1
      a=0.0d0
5000  continue
      a = a -log(1d0-rn(0))/vopt
c     if(a.ge.tdns) cycle
      if(a.ge.frnum(n)) cycle
      if(i0.eq.1) then
        tdns=mult(tf,bmu,smu)
        rfac=frnum(n)/tdns
        a=a/rfac
      endif

      do i=i0,maxkc
        if(a.ge.scsum(i-1).and. a.lt.scsum(i)) then
         ka=kcsum(i)
         i0=i
         goto 1000
        endif
      end do
      print *,'no particle?',a,scsum(maxkc),tdns,frnum(n),i0
      cycle

1000  kc=abs(ka) 
      sgn=isign(1,ka)
      ibary=kchg(kc,6)*sgn
      str=kchg(kc,7)*sgn
      cha=kchg(kc,1)*sgn/3
      if(opt.eq.2) then
          db=nbar-b_tot
          ds=nstr-s_tot
          dc=nch-c_tot
          probb=exp(-min(50d0,abs(db)))
          probs=exp(-min(50d0,abs(ds)))
          probc=exp(-min(50d0,abs(dc)))
          if(ibary*db.gt.0d0.and.rn(0).gt.probb) goto 5000 
          if(str*ds.gt.0d0.and.rn(0).gt.probs)  goto 5000
          if(cha*dc.gt.0d0.and.rn(0).gt.probc)  goto 5000
      endif
      kf=kchg(kc,4)
      kf=kf*sgn
      istat=1
      if(ibary.eq.0) istat=-1
      cmu = str*smu + ibary*bmu/3
      pc(5)=pjmass(kf)
100   pp=thermaldist3(pc(5),tf,cmu,istat)
      cost=2d0*rn(0)-1d0
      sint=sqrt(1d0-cost**2)
      phi=paru(2)*rn(0)
      pc(1)=pp*sint*cos(phi)
      pc(2)=pp*sint*sin(phi)
      pc(3)=pp*cost
      pc(4)=sqrt(pp**2+pc(5)**2)
      er=pc(4)
      call jamrobo(0d0,0d0,v1,v2,v3,gam,pc(1),pc(2),pc(3),pc(4))
      pds = pc(4)*dst+pc(1)*dsx+pc(2)*dsy+pc(3)*dsz
      if(pds.le.0d0) goto 5000
      prob=pds/(er*vopt)
      if(prob.lt.0.0.or.prob.gt.1.0) then
          print *,'prob?',prob,vopt
          vopt=vopt*1.2d0
          goto 100
      endif
      if(rn(0).gt.prob) goto 5000

c     if(etot+pc(4)/ntest.gt.pe_tot) then
c       if(abs(etot+pc(4)/ntest-pe_tot).gt.abs(etot-pe_tot)) return
c     endif

c...Set coordinate of particle.
      rc(1)=dx*( frr(1,n) - origx - 0.5 + rn(0) )
      rc(2)=dy*( frr(2,n) - origy - 0.5 + rn(0) )
      rc(3)=dz*( frr(3,n) - origz - 0.5 + rn(0) )

c...Set time and formation time of particle.
        rc(4)=frt(n)
        rc(5)=rc(4)
        call p2jam(kc,kf,ibary,rc,pc)
c       pe_tot=pe_tot-pc(4)
c       px_tot=px_tot-pc(1)
c       py_tot=py_tot-pc(2)
c       pz_tot=pz_tot-pc(3)
c       b_tot=b_tot-ibary
        ptot=ptot+pc

c....Reset fluid elements.
c       call reset_fluid_element(frr(1,n),frr(2,n),frr(3,n))

        npos=npos+1
        ntot = ntot + 1
        nftot=ntot
        nbar = nbar + ibary/3/ntest
        nstr = nstr + str/ntest
        nch = nch + cha/ntest
        etot = etot + pc(4)/ntest

        end do


      end

c***************************************************************************
      subroutine sampling2pos(npos,nbar,nch,nstr,ptot,etot,nftot,ntot)

c...MS sampling of particles from freeze-out surface by KM integral.
      implicit none
      integer i,n,npos,npar,nbar,nch,nstr,kc,kf,ibary,str,cha,ntot
      integer mode,mode1,mode2,nmode(0:7),istat
      real*8 rn,v1,v2,v3,dst,dsx,dsy,dsz,gam,vol,tdns,nftot
      real*8 tf,bmu,smu,cmu,ene,etot,dtot,rc(5),pc(5),ptot(5)
      real*8 mult_pos,pjmass,xran,sigma
      real*8 sp_potential
      integer generate_poisson,generate_poisson1
      include 'fluid.inc'
      include 'jam2.inc'

      npos=0
      nbar=0
      nch=0
      nstr=0
      ptot=0d0
      etot=0d0
      dtot=0d0
      ntot=0

      do n=1,nfreezeout

      v1=frv(1,n)
      v2=frv(2,n)
      v3=frv(3,n)
      dst=frdsgm(0,n)
      dsx=frdsgm(1,n)
      dsy=frdsgm(2,n)
      dsz=frdsgm(3,n)
      gam=1.0d0/sqrt(1.0d0-v1*v1-v2*v2-v3*v3)
      vol=gam*(dst + dsx*v1 + dsy*v2 + dsz*v3)
      sigma=sqrt(vol**2-dst**2+dsx**2+dsy**2+dsz**2)
      tf=frtmp(n)*hbc
      bmu=frmuq(n)*hbc
      smu=frmus(n)*hbc
c     bmu = bmu - pard(142)*frbdn(n)
c     bmu = bmu - sp_potential(frbdn(n))
      pare(51)=tf
      pare(52)=bmu
      pare(53)=smu

      tdns=mult_pos(tf,bmu,smu,vol,sigma)
      npar=generate_poisson(tdns)*ntest
      if(npar.eq.0) cycle

      dtot = dtot + tdns*ntest
      npos=npos+1

      do i=1,npar

        call selectparticle(tf,bmu,smu,tdns,nstr,nbar,nch,kc,kf)
        ibary=kchg(kc,6)*isign(1,kf)
        str=kchg(kc,7)*isign(1,kf)
        cha=kchg(kc,1)*isign(1,kf)/3
        istat=1
        if(ibary.eq.0) istat=-1
        cmu = str*smu + ibary*bmu/3
        pc(5)=pjmass(kf)
        call cooperfrye(dst,dsx,dsy,dsz,v1,v2,v3,tf,cmu,istat,pc)

c...Set coordinate of particle.
        rc(1)=dx*( frr(1,n) - origx - 0.5 + rn(0) )
        rc(2)=dy*( frr(2,n) - origy - 0.5 + rn(0) )
        rc(3)=dz*( frr(3,n) - origz - 0.5 + rn(0) )

c...Set time and formation time of particle.
        rc(4)=frt(n)
        rc(5)=rc(4)
        call p2jam(kc,kf,ibary,rc,pc)
c       pe_tot=pe_tot-pc(4)
c       px_tot=px_tot-pc(1)
c       py_tot=py_tot-pc(2)
c       pz_tot=pz_tot-pc(3)
c       b_tot=b_tot-ibary
        ptot=ptot+pc

c....Reset fluid elements.
c       call reset_fluid_element(frr(1,n),frr(2,n),frr(3,n))

        ntot = ntot + 1
        nbar = nbar + ibary/3/ntest
        nstr = nstr + str/ntest
        nch = nch + cha/ntest
        etot = etot + pc(4)/ntest

      end do

      end do


      end

c***************************************************************************
      subroutine sampling3(npos,nbar,nch,nstr,ptot,etot,nftot,ntot)

c...MC sampling of particles from freeze-out surface by the KM formula.

      implicit none
      integer i,n,npos,npar,nbar,nch,nstr,kc,kf,ibary,str,cha,ntot
      integer mode,mode1,mode2,nmode(0:7),istat
      real*8 rn,v1,v2,v3,dst,dsx,dsy,dsz,gam,vol,tdns,pdns,nftot
      real*8 tf,bmu,smu,cmu,ene,etot,rc(5),ptot(5),tdns2
      real*8 cost,sint,phi,er,pc(5),thermaldist3,pds,pp,prob
      real*8 pjmass,vopt,spin,fac,pdensity,pdensity_pos,sigma
      real*8 pdensity_pos2,sigma2,sigma0
      real*8 dsxr,dsyr,dszr,dstr
      real*8 probb,probs,probc,db,ds,dc
      real*8 mult_pos,tdtot,mult
      real*8 sp_potential
      integer j,k,opt,sgn,opt_positive,np
      integer generate_poisson
      include 'fluid.inc'
      include 'jam2.inc'

      fac=1.0d0/(2*paru(1)**2*paru(3)**3)
      npos=0
      nbar=0
      nch=0
      nstr=0
      ptot=0d0
      etot=0d0
      ntot=0
      opt=mod(mstc(148),10)
      opt_positive=1
      if(mstc(148).eq.33) opt_positive=0

      do 3000 n=1,nfreezeout

      v1=frv(1,n)
      v2=frv(2,n)
      v3=frv(3,n)
      dst=frdsgm(0,n)
      dsx=frdsgm(1,n)
      dsy=frdsgm(2,n)
      dsz=frdsgm(3,n)
      gam=1.0d0/sqrt(1.0d0-v1*v1-v2*v2-v3*v3)
      vol=gam*(dst + dsx*v1 + dsy*v2 + dsz*v3)
      if(opt_positive.eq.0.and.vol.le.0d0) goto 3000
      sigma=sqrt(vol**2-dst**2+dsx**2+dsy**2+dsz**2)
      tf=frtmp(n)*hbc
      bmu=frmuq(n)*hbc
      smu=frmus(n)*hbc
      pare(51)=tf
      pare(52)=bmu
      pare(53)=smu

c...Include potential.
c     bmu = bmu - pard(142)*frbdn(n)
c     bmu = bmu - sp_potential(frbdn(n))

      tdtot=0.0d0
      np=0
c...Loop over particle species.
      do kc=132,366
        if(kc.ge.195.and.kc.le.251) cycle ! exclude charm
        if(kc.eq.171.or.kc.eq.172) cycle ! K_L0 or K_S0
      do k=1,2
        if(k.eq.1) then
          i=kc
          sgn=1.0d0
        else if(k.eq.2) then
          if(kchg(kc,3).eq.0) cycle
          i=-kc
          sgn=-1.0d0
        endif
        ibary=kchg(kc,6)*sgn
        str=kchg(kc,7)*sgn
        cha=kchg(kc,1)*sgn/3
        if(opt.eq.2) then
          db=nbar-b_tot
          ds=nstr-s_tot
          dc=nch-c_tot
          probb=exp(-min(50d0,abs(db)))
          probs=exp(-min(50d0,abs(ds)))
          probc=exp(-min(50d0,abs(dc)))
          if(ibary*db.gt.0d0.and.rn(0).gt.probb) cycle
          if(str*ds.gt.0d0.and.rn(0).gt.probs) cycle
          if(cha*dc.gt.0d0.and.rn(0).gt.probc) cycle
        endif
        kf=kchg(kc,4)
        spin=max(1,mod(kf,10))
        kf=kf*sgn
        istat=1
        if(ibary.eq.0) istat=-1
        cmu = str*smu + ibary*bmu/3
        pc(5)=pmas(kc,1)

        if(opt_positive.eq.1) then
          tdns=pdensity_pos(pc(5),tf,cmu,ibary,istat,vol,sigma)*spin*fac
c         tdns2=pdensity_pos2(pc(5),tf,cmu,istat,vol,sigma)*spin*fac
        else
          tdns=vol*pdensity(pc(5),tf,cmu,ibary,istat)*spin*fac
        endif
        npar=generate_poisson(tdns)*ntest
        tdtot=tdtot+tdns
        np=np+1
        if(npar.le.0) cycle

        do i=1,npar
        pc(5)=pjmass(kf)
        call cooperfrye(dst,dsx,dsy,dsz,v1,v2,v3,tf,cmu,istat,pc)

c...Set coordinate of particle.
        rc(1)=dx*( frr(1,n) - origx - 0.5 + rn(0) )
        rc(2)=dy*( frr(2,n) - origy - 0.5 + rn(0) )
        rc(3)=dz*( frr(3,n) - origz - 0.5 + rn(0) )

c...Set time and formation time of particle.
        rc(4)=frt(n)
        rc(5)=rc(4)
        call p2jam(kc,kf,ibary,rc,pc)
c       pe_tot=pe_tot-pc(4)
c       px_tot=px_tot-pc(1)
c       py_tot=py_tot-pc(2)
c       pz_tot=pz_tot-pc(3)
c       b_tot=b_tot-ibary
        ptot=ptot+pc

c....Reset fluid elements.
c       call reset_fluid_element(frr(1,n),frr(2,n),frr(3,n))

        npos=npos+1
        ntot = ntot + 1
        nftot=ntot
        nbar = nbar + ibary/3/ntest
        nstr = nstr + str/ntest
        nch = nch + cha/ntest
        etot = etot + pc(4)/ntest

        end do

      end do
      end do

c     print *,tdtot,mult_pos(tf,bmu,smu,vol,sigma)
c     print *,np,tdtot,mult(tf,bmu,smu)*vol

3000  continue

      end

c***************************************************************************
      subroutine sampling2(npos,nbar,nch,nstr,ptot,etot,nftot,ntot)

c...Particles are randomly chosen from freeze-out surface until total energy
c...is filled.
      implicit none
      integer i,n,npos,npar,nbar,nch,nstr,kc,kf,ibary,str,cha,ntot
      integer mode,mode1,mode2,nmode(0:7),istat
      real*8 rn,v1,v2,v3,dst,dsx,dsy,dsz,gam,vol,tdns,pdns,nftot
      real*8 tf,bmu,smu,cmu,ene,etot,dtot,rc(5),pc(5),ptot(5)
      real*8 mult,pjmass,xran
      real*8 sp_potential
      integer generate_poisson,generate_poisson1
      integer opt_select
      include 'fluid.inc'
      include 'jam2.inc'

      npos=0
      nbar=0
      nch=0
      nstr=0
      ptot=0d0
      etot=0d0
      dtot=0d0
      ntot=0

      opt_select=mstc(150)

      if(opt_select.eq.1) then
        call totalmult(nftot)
        if(nftot.le.0d0) return
      endif

c....Loop over freeze-out surface.
 1000 continue

      if(opt_select.eq.1) then
 2100   n=0
        xran=nftot*rn(0)
 2000   n=n+1
        if(n.gt.nfreezeout) goto 2100
        xran=xran-frnum(n)
        if(xran.gt.0d0) goto 2000
c...At least one particle is generated.
        tdns=frnum(n)
c       npar=generate_poisson1(tdns)*ntest
        npar=1

      else

 3000   n=int(1+(nfreezeout-1)*rn(0))
c...Determine number of particle to be generated.
        tdns=frnum(n)
        if(tdns.le.0d0) goto 3000
        npar=generate_poisson(tdns)*ntest
        if(npar.eq.0) goto 3000
      endif


      v1=frv(1,n)
      v2=frv(2,n)
      v3=frv(3,n)
      dst=frdsgm(0,n)
      dsx=frdsgm(1,n)
      dsy=frdsgm(2,n)
      dsz=frdsgm(3,n)
      gam=1.0d0/sqrt(1.0d0-v1*v1-v2*v2-v3*v3)
      vol=gam*(dst + dsx*v1 + dsy*v2 + dsz*v3)
      tf=frtmp(n)*hbc
      bmu=frmuq(n)*hbc
      smu=frmus(n)*hbc
      pare(51)=tf
      pare(52)=bmu
      pare(53)=smu

      pdns=tdns/vol
      dtot = dtot + tdns*ntest

      npos=npos+1

c...Include potential.
c     if(opt_table.eq.2) then
c     bmu = bmu - pard(142)*frbdn(n)
c     bmu = bmu - sp_potential(frbdn(n))
c     endif
      pdns=mult(tf,bmu,smu)

      do i=1,npar

        call selectparticle(tf,bmu,smu,pdns,nstr,nbar,nch,kc,kf)
        ibary=kchg(kc,6)*isign(1,kf)
        str=kchg(kc,7)*isign(1,kf)
        cha=kchg(kc,1)*isign(1,kf)/3

        istat=1
        if(ibary.eq.0) istat=-1
        cmu = str*smu + ibary*bmu/3

        pc(5)=pjmass(kf)
        call cooperfrye(dst,dsx,dsy,dsz,v1,v2,v3,tf,cmu,istat,pc)

        if(etot+pc(4)/ntest.gt.pe_tot) then
          if(abs(etot+pc(4)/ntest-pe_tot).gt.abs(etot-pe_tot)) return
        endif

c...Set coordinate of particle.
        rc(1)=dx*( frr(1,n) - origx - 0.5 + rn(0) )
        rc(2)=dy*( frr(2,n) - origy - 0.5 + rn(0) )
        rc(3)=dz*( frr(3,n) - origz - 0.5 + rn(0) )

c...Set time and formation time of particle.
        rc(4)=frt(n)
        rc(5)=rc(4)
        call p2jam(kc,kf,ibary,rc,pc)
c       pe_tot=pe_tot-pc(4)
c       px_tot=px_tot-pc(1)
c       py_tot=py_tot-pc(2)
c       pz_tot=pz_tot-pc(3)
c       b_tot=b_tot-ibary
        ptot=ptot+pc

c....Reset fluid elements.
c       call reset_fluid_element(frr(1,n),frr(2,n),frr(3,n))

        ntot = ntot + 1
        nbar = nbar + ibary/3/ntest
        nstr = nstr + str/ntest
        nch = nch + cha/ntest
        etot = etot + pc(4)/ntest

      end do

      if(etot.lt.pe_tot) goto 1000


      end

c***************************************************************************
      subroutine totalmult(nftot)
c...Compute multiplicity  for each freeze-out cell and sum up them.
      implicit none
      include 'fluid.inc'
      integer n,opt_surface
      real*8 tmax,mubmax,musmax,bmu,smu,tf
      real*8 multlookup,mult
      real*8 nftot,vol,gam,rz,dsgm2
      parameter(opt_surface=0)


      tmax=0d0
      mubmax=0d0
      musmax=0d0
      nftot=0d0
      do n=1,nfreezeout
        tf=frtmp(n)*hbc
        bmu=frmuq(n)*hbc
        smu=frmus(n)*hbc
        tmax=max(tf,tmax)
        mubmax=max(mubmax,bmu)
        musmax=max(musmax,smu)
        gam=1.0d0/sqrt(1.0d0-frv(1,n)**2-frv(2,n)**2-frv(3,n)**2)
        vol=gam*( frdsgm(0,n) + frdsgm(1,n)*frv(1,n) + 
     &        frdsgm(2,n)*frv(2,n) + frdsgm(3,n)*frv(3,n) )
        if(vol.le.0d0) then
          frnum(n)=0d0
          cycle
        endif

        rz=dx*(frr(3,n) - origz)
        if(frt(n).lt.abs(rz) ) then
          if(abs(frt(n)-abs(rz)).lt.1d-8) then
            rz=frt(n)
          else
c           print *,'acausal? ',n,frt(n),rz
            frnum(n)=0d0
            cycle
          endif
        endif

        if(opt_surface.ge.1) then
c...Exclude timelike surface.
          dsgm2=frdsgm(0,n)**2 -frdsgm(1,n)**2 - frdsgm(2,n)**2
     &           - frdsgm(3,n)**2
          if(dsgm2.le.0d0) then
            frnum(n)=0d0
            cycle
          endif
c...Exclude negative dsigma_t
          if(frdsgm(0,n).le.0d0) then
            frnum(n)=0d0
            cycle
          endif
        endif
         
        if(opt_table.ge.1) then
c         frnum(n)=multlookup(tf,bmu,smu)*vol
c       else if(opt_table.eq.2) then
          frnum(n)=frnum(n)*vol
        else
          frnum(n)=mult(tf,bmu,smu)*vol
        endif
        nftot=nftot+frnum(n)
        if(frnum(n).lt.0d0.or.frnum(n).gt.50d0) then
           print *,'frnum?',frnum(n),tf,bmu,smu,vol
           stop
        endif
      end do

      if(nftot.le.0d0) then
        print *,'totalmult nftot=0?',nftot
      endif

      if(job_mode.ge.1) then
 
      print *,'total n=',nftot,' Tmax=',tmax,' muB_max=',mubmax,
     &  'muS_max=',musmax

      endif

      end

c***************************************************************************
      subroutine modesampling(npos,nbar,nch,nstr,ptot,etot,nftot,ntot)

c...Particles are sampled according to a method "mode sampling".
      implicit none
      integer i,n,npos,npar,nbar,nch,nstr,kc,kf,ibary,str,cha,ntot
      integer mode,mode1,mode2,nmode(0:7),istat
      real*8 nftot,dsgm2
      real*8 xran,rn,v1,v2,v3,dst,dsx,dsy,dsz,gam,vol,tdns,pdns
      real*8 tf,bmu,smu,cmu,ene,etot,rc(5),pc(5),ptot(5)
      integer generate_poisson,generate_poisson1
      integer istot,ictot,ibtot
      real*8 multlookup,mult,pjmass
      real*8 nftot2
      real*8 sp_potential
      integer opt_select
      include 'fluid.inc'
      include 'jam2.inc'

      npos=0
      nbar=0
      nch=0
      nstr=0
      ptot=0d0
      etot=0d0
      ntot=0

      istot=nint(s_tot)
      ictot=nint(c_tot)
      ibtot=nint(b_tot)

      opt_select=mstc(150)

c....Compute total multiplicity by sum up all freeze-out surface.
      call totalmult(nftot)

      if(nftot.le.0d0) return

      if(nftot.lt.b_tot) then
       print *,'modesampling::nftot<b_tot?',nftot,b_tot
      endif


c...mode=1,2 strangeness conservation.
c...mode=3,4 baryon number conservation.
c...mode=5,6 charge conservation.
c...mode=7   energy conservation.
      nmode=0
      mode1=1
      mode2=7

c...Loop over mode sampling.
      do mode=mode1,mode2

        ene=0d0

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        if(mode.eq.3) then
          if(nstr.ne.istot) then
            print *,'after mode 2 strangeness does not conserve',
     &    nstr,istot
           print *,'npar=',npar,kf,str
          endif
        endif
        if(mode.eq.7) then
          if(nch.ne.ictot) then
            print *,'after mode 6 charge does not conserve',
     &    nch,ictot
           print *,'npar=',npar,kf,cha
          endif
        endif
        if(mode.eq.6) then
          if(nch.gt.ictot) then
            print *,'after mode 5 charge>total ch',
     &    nch,ictot
           print *,'npar=',npar,kf,cha
          endif
        endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

c....Loop over freeze-out surface.
 1000 continue

      if(opt_select.eq.1) then

 2100   n=0
        xran=nftot*rn(0)
 2000   n=n+1
        if(n.gt.nfreezeout) goto 2100
        xran=xran-frnum(n)
        if(xran.gt.0d0) goto 2000
        tdns=frnum(n)
c...At least one particle is generated.
c       npar=generate_poisson1(tdns)*ntest
        npar=1


      else

2200    n=int(1+(nfreezeout-1)*rn(0))
c...Determine number of particle to be generated.
        tdns=frnum(n)

        if(tdns.lt.0.01) then
          if(rn(0).ge.tdns) goto 2200
          npar=1
        else
          npar=generate_poisson(tdns)*ntest
          if(npar.eq.0) goto 2200
        endif

      endif

      v1=frv(1,n)
      v2=frv(2,n)
      v3=frv(3,n)
      dst=frdsgm(0,n)
      dsx=frdsgm(1,n)
      dsy=frdsgm(2,n)
      dsz=frdsgm(3,n)
      gam=1.0d0/sqrt(1.0d0-v1*v1-v2*v2-v3*v3)
      vol=gam*(dst + dsx*v1 + dsy*v2 + dsz*v3)
      if(vol.le.0d0) goto 1000

c...Exclude timelike surface.
c     dsgm2=dst*dst -dsx*dsx - dsy*dsy - dsz*dsz
c     if(dsgm2.le.0d0) goto 1000
c...Exclude negative dst.
c     if(dst.le.0d0) goto 1000

      npos=npos+1

c     pdns=tdns/vol
c     frnum(n)=0d0  ! Remove this cell not to select next.

      tf=frtmp(n)*hbc
      bmu=frmuq(n)*hbc
      smu=frmus(n)*hbc
      pare(51)=tf
      pare(52)=bmu
      pare(53)=smu
c...Include potential.
c     if(opt_table.eq.2) then
c       bmu = bmu - pard(142)*frbdn(n)
c     bmu = bmu - sp_potential(frbdn(n))
c     endif
      pdns=mult(tf,bmu,smu)

      do i=1,npar
        call selectparticle(tf,bmu,smu,pdns,nstr,nbar,nch,kc,kf)
        ibary=kchg(kc,6)*isign(1,kf)
        str=kchg(kc,7)*isign(1,kf)
        cha=kchg(kc,1)*isign(1,kf)/3

c...2) S<0 particle
        if(mode.eq.2) then
          if(str.ge.0) cycle
          if(nstr+str.lt.istot) cycle ! too many S
        endif

c...4) Non-strange baryon
        if(mode.eq.4) then
          if(str.ne.0.or.ibary.ne.3) cycle
          if(nbar+ibary/3.gt.ibtot) cycle
        endif

c...6) positively charged non-strange meson.
        if(mode.eq.6) then
          if(ibary.ne.0.or.str.ne.0.or.cha.le.0) cycle
          if(nch+cha.gt.ictot) cycle
        endif

c...7) non-strange neutral meson.
        if(mode.eq.7.and.(ibary.ne.0.or.str.ne.0.or.cha.ne.0)) cycle

        istat=1
        if(ibary.eq.0) istat=-1
        cmu = str*smu + ibary*bmu/3
        pc(5)=pjmass(kf)
        call cooperfrye(dst,dsx,dsy,dsz,v1,v2,v3,tf,cmu,istat,pc)

        if(mode.eq.1.or.mode.eq.3.or.mode.eq.5) ene=ene+pc(4)/ntest

c...1) S>0 particle.
        if(mode.eq.1.and.str.le.0) cycle

c...3) None-strane anti-baryon.
        if(mode.eq.3.and.(ibary.ne.-3.or.str.ne.0)) cycle

c...5) Non-strange negative mesons.
        if(mode.eq.5) then
          if(ibary.ne.0) cycle
          if(str.ne.0) cycle
          if(cha.ge.0) cycle
        endif

        if(mode.eq.7.and.etot+pc(4)/ntest.gt.pe_tot) then
          if(abs(etot+pc(4)/ntest-pe_tot).gt.abs(etot-pe_tot)) goto 200
        endif

c...Set coordinate of particle.
        rc(1)=dx*( frr(1,n) - origx - 0.5 + rn(0) )
        rc(2)=dy*( frr(2,n) - origy - 0.5 + rn(0) )
        rc(3)=dz*( frr(3,n) - origz - 0.5 + rn(0) )

c...Set time and formation time of particle.
        rc(4)=frt(n)
        rc(5)=rc(4)
        call p2jam(kc,kf,ibary,rc,pc)
c       pe_tot=pe_tot-pc(4)
c       px_tot=px_tot-pc(1)
c       py_tot=py_tot-pc(2)
c       pz_tot=pz_tot-pc(3)
c       b_tot=b_tot-ibary

c....Reset fluid elements.
c       call reset_fluid_element(frr(1,n),frr(2,n),frr(3,n))

c       dtot = dtot + tdns/ntest
        ntot = ntot + 1
        nbar = nbar + ibary/3/ntest
        nstr = nstr + str/ntest
        nch = nch + cha/ntest
        etot = etot + pc(4)/ntest
        ptot=ptot+pc/ntest

        nmode(mode)=nmode(mode)+1

      end do


      if(mode.eq.1) then
        if(ene.lt.pe_tot) goto 1000
        if(nstr.le.s_tot) goto 1000
      endif
      if(mode.eq.3) then
        if(ene.lt.pe_tot) goto 1000
        if(nbar.ge.ibtot) goto 1000
      endif
      if(mode.eq.5) then
        if(ene.lt.pe_tot) goto 1000
        ! mode 5: current charge is larger than total charge.
        ! so we generate more negative charge.
        if(nch.ge.ictot) goto 1000
      endif

      if(mode.eq.7.and.etot.lt.pe_tot) goto 1000

      if(mode.eq.2.and.nstr.gt.istot) goto 1000
      if(mode.eq.4.and.nbar.lt.ibtot) goto 1000
      if(mode.eq.6.and.nch.lt.ictot) goto 1000

      end do ! mode loop
 200  continue
c...Sampling done.

      if(mstc(8).ge.1) then

      do i=0,7
       print *,'mode=',i,nmode(i)
      end do
    
      if(abs(etot-pe_tot)/pe_tot.ge.0.1) then
        print *,'energy does not conserve after hydro2jam',etot,pe_tot
        stop
      endif
      if(nbar.ne.ibtot) then
        print *,'baryon number does not conserve',nbar,ibtot
        stop
      endif
      if(nstr.ne.istot) then
        print *,'strangeness does not conserve',nstr,istot,s_tot
        stop
      endif
      if(nch.ne.ictot) then
        print *,'charge does not conserve',nch,ictot
        stop
      endif

      endif

      end

c***************************************************************************
      subroutine modesampling_pos(npos,nbar,nch,nstr,ptot,etot,
     &      nftot,ntot)

c...Particles are sampled according to a method "mode sampling".
      implicit none
      integer i,n,npos,npar,nbar,nch,nstr,kc,kf,ibary,str,cha,ntot
      integer mode,mode1,mode2,nmode(0:7),istat
      real*8 nftot,dsgm2,vopt
      real*8 xran,rn,v1,v2,v3,dst,dsx,dsy,dsz,gam,vol,tdns,pdns
      real*8 tf,bmu,smu,cmu,ene,etot,rc(5),pc(5),ptot(5)
      integer generate_poisson,generate_poisson1
      integer istot,ictot,ibtot
      real*8 multlookup,mult,pjmass
      real*8 nftot2,pds,er,prob,thermaldist3,sint,cost,phi,pp
      real*8 sp_potential
      integer opt_select
      include 'fluid.inc'
      include 'jam2.inc'

      npos=0
      nbar=0
      nch=0
      nstr=0
      ptot=0d0
      etot=0d0
      ntot=0

      istot=nint(s_tot)
      ictot=nint(c_tot)
      ibtot=nint(b_tot)

      opt_select=mstc(150)

c....Compute total multiplicity by sum up all freeze-out surface.
      if(opt_select.eq.1) then

      nftot=0d0
      do n=1,nfreezeout
        gam=1.0d0/sqrt(1.0d0-frv(1,n)**2-frv(2,n)**2-frv(3,n)**2)
        vol=gam*( frdsgm(0,n) + frdsgm(1,n)*frv(1,n) + 
     &        frdsgm(2,n)*frv(2,n) + frdsgm(3,n)*frv(3,n) )
        dsgm2=frdsgm(0,n)**2 -frdsgm(1,n)**2 - frdsgm(2,n)**2
     &           - frdsgm(3,n)**2
        vopt=vol+sqrt(vol**2 - dsgm2)
        if(opt_table.ge.1) then
          frnum(n)=frnum(n)*vopt
        else
          tf=frtmp(n)*hbc
          bmu=frmuq(n)*hbc
          smu=frmus(n)*hbc
          frnum(n)=mult(tf,bmu,smu)*vopt
        endif
        nftot=nftot+frnum(n)
      end do

      if(job_mode.ge.1) then
      print *,'total n=',nftot
      endif

      if(nftot.le.0d0) return

      if(nftot.lt.b_tot) then
       print *,'modesampling::nftot<b_tot?',nftot,b_tot
      endif

      endif


c...mode=1,2 strangeness conservation.
c...mode=3,4 baryon number conservation.
c...mode=5,6 charge conservation.
c...mode=7   energy conservation.
      nmode=0
      mode1=1
      mode2=7

c...Loop over mode sampling.
      do mode=mode1,mode2

        ene=0d0

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        if(mode.eq.3) then
          if(nstr.ne.istot) then
            print *,'after mode 2 strangeness does not conserve',
     &    nstr,istot
           print *,'npar=',npar,kf,str
          endif
        endif
        if(mode.eq.7) then
          if(nch.ne.ictot) then
            print *,'after mode 6 charge does not conserve',
     &    nch,ictot
           print *,'npar=',npar,kf,cha
          endif
        endif
        if(mode.eq.6) then
          if(nch.gt.ictot) then
            print *,'after mode 5 charge>total ch',
     &    nch,ictot
           print *,'npar=',npar,kf,cha
          endif
        endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

c....Loop over freeze-out surface.
 1000 continue

      if(opt_select.eq.1) then

 2100   n=0
        xran=nftot*rn(0)
 2000   n=n+1
        if(n.gt.nfreezeout) goto 2100
        xran=xran-frnum(n)
        if(xran.gt.0d0) goto 2000
        tdns=frnum(n)
c...At least one particle is generated.
        npar=generate_poisson1(tdns)*ntest

      else

2200    n=int(1+(nfreezeout-1)*rn(0))
c...Determine number of particle to be generated.
        tdns=frnum(n)

c       if(tdns.lt.0.01) then
c         if(rn(0).ge.tdns) goto 2200
c         npar=1
c       else
          npar=generate_poisson(tdns)*ntest
          if(npar.eq.0) goto 2200
c       endif

      endif

      v1=frv(1,n)
      v2=frv(2,n)
      v3=frv(3,n)
      dst=frdsgm(0,n)
      dsx=frdsgm(1,n)
      dsy=frdsgm(2,n)
      dsz=frdsgm(3,n)
      gam=1.0d0/sqrt(1.0d0-v1*v1-v2*v2-v3*v3)
      vol=gam*(dst + dsx*v1 + dsy*v2 + dsz*v3)
      vopt=vol+sqrt(vol**2 - (dst**2 - dsx**2-dsy**2-dsz**2))

      npos=npos+1

      tf=frtmp(n)*hbc
      bmu=frmuq(n)*hbc
      smu=frmus(n)*hbc
      pare(51)=tf
      pare(52)=bmu
      pare(53)=smu
c...Include potential.
c     bmu = bmu - pard(142)*frbdn(n)
c     bmu = bmu - sp_potential(frbdn(n))
      pdns=mult(tf,bmu,smu)

      do i=1,npar
        call selectparticle(tf,bmu,smu,pdns,nstr,nbar,nch,kc,kf)
        ibary=kchg(kc,6)*isign(1,kf)
        str=kchg(kc,7)*isign(1,kf)
        cha=kchg(kc,1)*isign(1,kf)/3

c...2) S<0 particle
        if(mode.eq.2) then
          if(str.ge.0) cycle
          if(nstr+str.lt.istot) cycle ! too many S
        endif

c...4) Non-strange baryon
        if(mode.eq.4) then
          if(str.ne.0.or.ibary.ne.3) cycle
          if(nbar+ibary/3.gt.ibtot) cycle
        endif

c...6) positively charged non-strange meson.
        if(mode.eq.6) then
          if(ibary.ne.0.or.str.ne.0.or.cha.le.0) cycle
          if(nch+cha.gt.ictot) cycle
        endif

c...7) non-strange neutral meson.
        if(mode.eq.7.and.(ibary.ne.0.or.str.ne.0.or.cha.ne.0)) cycle

        istat=1
        if(ibary.eq.0) istat=-1
        cmu = str*smu + ibary*bmu/3
        pc(5)=pjmass(kf)
c       call cooperfrye(dst,dsx,dsy,dsz,v1,v2,v3,tf,cmu,istat,pc)
100     pp=thermaldist3(pc(5),tf,cmu,istat)
        cost=2d0*rn(0)-1d0
        sint=sqrt(1d0-cost**2)
        phi=paru(2)*rn(0)
        pc(1)=pp*sint*cos(phi)
        pc(2)=pp*sint*sin(phi)
        pc(3)=pp*cost
        pc(4)=sqrt(pp**2+pc(5)**2)
        er=pc(4)
        call jamrobo(0d0,0d0,v1,v2,v3,gam,pc(1),pc(2),pc(3),pc(4))
        pds = pc(4)*dst+pc(1)*dsx+pc(2)*dsy+pc(3)*dsz
        if(pds.le.0d0) cycle
        prob=pds/(er*vopt)
        if(rn(0).gt.prob) cycle

        if(mode.eq.1.or.mode.eq.3.or.mode.eq.5) ene=ene+pc(4)/ntest

c...1) S>0 particle.
        if(mode.eq.1.and.str.le.0) cycle

c...3) None-strane anti-baryon.
        if(mode.eq.3.and.(ibary.ne.-3.or.str.ne.0)) cycle

c...5) Non-strange negative mesons.
        if(mode.eq.5) then
          if(ibary.ne.0) cycle
          if(str.ne.0) cycle
          if(cha.ge.0) cycle
        endif

        if(mode.eq.7.and.etot+pc(4)/ntest.gt.pe_tot) then
          if(abs(etot+pc(4)/ntest-pe_tot).gt.abs(etot-pe_tot)) goto 200
        endif

c...Set coordinate of particle.
        rc(1)=dx*( frr(1,n) - origx - 0.5 + rn(0) )
        rc(2)=dy*( frr(2,n) - origy - 0.5 + rn(0) )
        rc(3)=dz*( frr(3,n) - origz - 0.5 + rn(0) )

c...Set time and formation time of particle.
        rc(4)=frt(n)
        rc(5)=rc(4)
        call p2jam(kc,kf,ibary,rc,pc)
c       pe_tot=pe_tot-pc(4)
c       px_tot=px_tot-pc(1)
c       py_tot=py_tot-pc(2)
c       pz_tot=pz_tot-pc(3)
c       b_tot=b_tot-ibary

c....Reset fluid elements.
c       call reset_fluid_element(frr(1,n),frr(2,n),frr(3,n))

c       dtot = dtot + tdns/ntest
        ntot = ntot + 1
        nbar = nbar + ibary/3/ntest
        nstr = nstr + str/ntest
        nch = nch + cha/ntest
        etot = etot + pc(4)/ntest
        ptot=ptot+pc/ntest

        nmode(mode)=nmode(mode)+1

      end do


      if(mode.eq.1) then
        if(ene.lt.pe_tot) goto 1000
        if(nstr.le.s_tot) goto 1000
      endif
      if(mode.eq.3) then
        if(ene.lt.pe_tot) goto 1000
        if(nbar.ge.ibtot) goto 1000
      endif
      if(mode.eq.5) then
        if(ene.lt.pe_tot) goto 1000
        ! mode 5: current charge is larger than total charge.
        ! so we generate more negative charge.
        if(nch.ge.ictot) goto 1000
      endif

      if(mode.eq.7.and.etot.lt.pe_tot) goto 1000

      if(mode.eq.2.and.nstr.gt.istot) goto 1000
      if(mode.eq.4.and.nbar.lt.ibtot) goto 1000
      if(mode.eq.6.and.nch.lt.ictot) goto 1000

      end do ! mode loop
 200  continue
c...Sampling done.

      if(mstc(8).ge.1) then

      do i=0,7
       print *,'mode=',i,nmode(i)
      end do
    
      if(abs(etot-pe_tot)/pe_tot.ge.0.1) then
        print *,'energy does not conserve after hydro2jam',etot,pe_tot
        stop
      endif
      if(nbar.ne.ibtot) then
        print *,'baryon number does not conserve',nbar,ibtot
        stop
      endif
      if(nstr.ne.istot) then
        print *,'strangeness does not conserve',nstr,istot,s_tot
        stop
      endif
      if(nch.ne.ictot) then
        print *,'charge does not conserve',nch,ictot
        stop
      endif

      endif

      end


c***********************************************************************
      subroutine selectparticle(tf,bmu,smu,tdns,ns,nb,nc,kc,kf)

c...Before calling this, mult() should be called to compute pmult().     
      implicit none
      include 'fluid.inc'
      include 'jam2.inc'
      integer i,j,k,kc,kf,ns,nb,nc,opt,sgn,iba,str,cha
      real*8 pc(5),cmu,xpart,tf,bmu,smu,vol,rn,tdns2
      real*8 pmult,mult,smult,probb,probs,probc,tdns,db,ds,dc
      common/bolz3/pmult(-500:500)

      opt=mod(mstc(148),10)
      if(opt.eq.2) then
        db=nb-b_tot
        ds=ns-s_tot
        dc=nc-c_tot
        probb=exp(-min(50d0,abs(db)))
        probs=exp(-min(50d0,abs(ds)))
        probc=exp(-min(50d0,abs(dc)))
      endif

c....Use the value from the table.
c200  xpart=pdns*rn(0)
c...Compute multiplicities for each particle.
c     tdns=mult(tf,bmu,smu)

 200  xpart=tdns*rn(0)
      tdns2=0d0
c     do i=-366,366
      do j=132,366
        if(j.ge.195.and.j.le.251) cycle
        if(j.eq.171.or.j.eq.172) cycle ! K_L0 or K_S0
c         pmult(i)= smult(i,tf,bmu,smu)
      do k=1,2
          if(k.eq.1) i=j
          if(k.eq.2) then
            if(kchg(j,3).eq.0) cycle
            i=-j
          endif
          xpart=xpart-pmult(i)
          tdns2=tdns2+pmult(i)
          if(xpart.le.0d0) then
            kc=abs(i)
            if(opt.eq.2) then
              sgn=isign(1,i)
              iba=kchg(kc,6)*sgn
              str=kchg(kc,7)*sgn
              cha=kchg(kc,1)*sgn
              if(iba*db.gt.0d0.and.rn(0).gt.probb) goto 200
              if(str*ds.gt.0d0.and.rn(0).gt.probs) goto 200
              if(cha*dc.gt.0d0.and.rn(0).gt.probc) goto 200
            endif
            kf=kchg(kc,4)*isign(1,i)
            goto 100
          endif
      end do
      end do
      if(mstc(148).eq.0) then
        print *,'tdns?',tdns,tdns2
      endif
      go to 200
100   continue


      end


c***********************************************************************
      subroutine p2jam(kc,kf,ibary,rc,pc)

c....Put particle into JAM.

      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      integer mdel,indd(100),idel(100)
      real*8 jamdtim
      real*8 pc(5),rc(5)

      mdel=0
      npar=0

c...Find a slot of this particle.
        if(ibary.eq.0) then
          nmeson=nmeson+1
          nv=nv+1
          ip=nv
          if(nv.gt.mxv) then
            print *,'jam2hydro::sampleparticle nv > mxv',nv
            stop
          endif
        else
           call jamindb(ip,indd,npar,mdel,idel,2)
        endif

c...Zero the vector.
        call jamzero(ip)

c....Resonance or not?
        if(pmas(kc,2).le.1d-7.or.mdcy(kc,1).eq.0
     $              .or.mdcy(kc,2).eq.0.or.mdcy(kc,3).eq.0)then
          k(1,ip)=1
        else
          k(1,ip)=2
        endif

        k(2,ip)=kf
        k(3,ip)=0
        k(4,ip)=0
        k(5,ip)=-1
        k(6,ip)=0
c       k(7,ip)=2
        k(7,ip)=-2000
        k(8,ip)=mstd(8)
        k(9,ip)=ibary
        k(10,ip)=0
        k(11,ip)=0

        p(:,ip)=pc(:)
        r(:,ip)=rc(:)

c...Vertex.
        v(1,ip)=r(1,ip)
        v(2,ip)=r(2,ip)
        v(3,ip)=r(3,ip)
        v(4,ip)=r(4,ip)

c.....Set resonance decay time.
        v(5,ip)=1.d+35
        if(k(1,ip).eq.2)
     $  v(5,ip)=r(5,ip)+jamdtim(1,kf,kc,k(1,ip),p(5,ip),p(4,ip))

c         write(65,*)sqrt(r(1,ip)**2+r(2,ip)**2+r(3,ip)**2)
c       print *,'particle is created!',nv,k(2,ip),p(5,ip),
c    & sqrt(r(1,ip)**2+r(2,ip)**2+r(3,ip)**2)

c...Update collision list.
c     call jamcupda(ip,-1,ip,0,1)

c     if(mdel.ge.1) then
c       do i=1,mdel
c         call jamcupda(idel(i),-1,ip,0,1)
c        end do
c     endif

c     i=mstc(38)
c     mstc(38)=6
c     call jamlist(1)
c     print *,'kf=',kf, ibary,pc(5)
c     print *,'nv=',nv,nbary,nmeson
c     mstc(38)=i

      end

c***************************************************************************

      subroutine recover_mom(iopt,ntot,ptot,pe_tot,px_tot,py_tot,pz_tot)
      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      integer iopt
      real*8 ptot(5),pf(5)

c..ntot: total number of particles produced from fluid.
c..ptot: total momentum of particles produced from fluid.
c..px_tot, py_tot, pz_tot: momentum

c...Shift momenta to match the required total momentum.
      pf=0d0
      nn=0
      do i=1,nv
       if(k(7,i).eq.-2000) then
         nn=nn+1
         p(1,i)=p(1,i)-1.0d0/ntot*(ptot(1)-px_tot)
         p(2,i)=p(2,i)-1.0d0/ntot*(ptot(2)-py_tot)
         p(3,i)=p(3,i)-1.0d0/ntot*(ptot(3)-pz_tot)
         p(4,i)=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
         pf(1)=pf(1)+p(1,i)
         pf(2)=pf(2)+p(2,i)
         pf(3)=pf(3)+p(3,i)
         pf(4)=pf(4)+p(4,i)
       endif
      end do

c     print *,'nv ntot=',nv,ntot,nn
c     print *,'ini    tot p=',px_tot,py_tot,pz_tot,pe_tot
c     print *,'after  tot p=',ptot(1),ptot(2),ptot(3),ptot(4)
c     print *,'after1 tot p=',pf(1),pf(2),pf(3),pf(4)
c     print *,'after1  p=',abs(px_tot-pf(1)),
c    &  abs(py_tot-pf(2)),
c    &  abs(pz_tot-pf(3)),
c    &  abs(pe_tot-pf(4))

c...Momentum conservation is recovered now.
      ptot=pf
      if(iopt.eq.1) goto 3000

c...Go to the C.M. frame of all particles.
      vx=px_tot/pe_tot
      vy=py_tot/pe_tot
      vz=pz_tot/pe_tot
      v2=vx**2+vy**2+yz**2
      if(v2.ge.1.0d0-1d-15) then
        print *,'recover_mom v>1?',v2
        return
      endif
      gam=1.0d0/sqrt(1.0d0-v2)
      do i=1,nv
      if(k(7,i).eq.-2000) then
       call jamrobo(0d0,0d0,-vx,-vy,-vz,gam,p(1,i),p(2,i),p(3,i),p(4,i))
      endif
      end do

c...Find a factor to correct total energy.
      pe=sqrt(pe_tot**2-px_tot**2-py_tot**2-pz_tot**2)
      a=1.0d0
      do iter=1,50
        f=-pe
        df=0.0d0
        do i=1,nv
        if(k(7,i).eq.-2000) then
          e=sqrt(p(5,i)**2+a*a*(p(1,i)**2+p(2,i)**2+p(3,i)**2))
          f = f + e
          df = df + a/e
        endif
        end do
        if(abs(f).lt.1d-8) goto 1000
        a = a - f/df
      end do

      print *,'recover_mom iter',iter,f,f/df,a

c...Go back to the computational frame after scaling the momenta.
 1000 continue
      pf=0.0d0
      do i=1,nv
      if(k(7,i).eq.-2000) then
        p(1,i)=a*p(1,i)
        p(2,i)=a*p(2,i)
        p(3,i)=a*p(3,i)
        p(4,i)=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
        call jamrobo(0d0,0d0,vx,vy,vz,gam,p(1,i),p(2,i),p(3,i),p(4,i))
        pf(1)=pf(1)+p(1,i)
        pf(2)=pf(2)+p(2,i)
        pf(3)=pf(3)+p(3,i)
        pf(4)=pf(4)+p(4,i)
      endif
      end do

c     print *,'a=',a,iter
c     print *,'nv ntot=',nv,ntot,nn
c     print *,'after2 p=',abs(px_tot-pf(1)),
c    &  abs(py_tot-pf(2)),
c    &  abs(pz_tot-pf(3)),
c    &  abs(pe_tot-pf(4))

c     if(abs(py_tot-pf(2)) .ge.1e-5
c    & .or.abs(pz_tot-pf(3)).ge.1e-5
c    & .or.abs(pe_tot-pf(4)).ge.1e-5) then
c     print *,'a=',a,iter
c     print *,'nv ntot=',nv,ntot,nn
c     print *,'after2 p=',abs(px_tot-pf(1)),
c    &  abs(py_tot-pf(2)),
c    &  abs(pz_tot-pf(3)),
c    &  abs(pe_tot-pf(4))
c     endif


c...This is not really needed.
      ptot=0d0
      do i=1,nv
       if(k(7,i).eq.-2000) then
         p(1,i)=p(1,i)-1.0d0/ntot*(pf(1)-px_tot)
         p(2,i)=p(2,i)-1.0d0/ntot*(pf(2)-py_tot)
         p(3,i)=p(3,i)-1.0d0/ntot*(pf(3)-pz_tot)
         p(4,i)=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
         ptot(1)=ptot(1)+p(1,i)
         ptot(2)=ptot(2)+p(2,i)
         ptot(3)=ptot(3)+p(3,i)
         ptot(4)=ptot(4)+p(4,i)
       endif
      end do

c     print *,'after3 % p=',abs(px_tot-pf(1)),
c    &  abs(py_tot-pf(2)),
c    &  abs(pz_tot-pf(3)),
c    &  abs(pe_tot-pf(4))

c     do i=1,nv
c       write(77,*)r(1,i),r(2,i),r(3,i),r(4,i)
c     end do

3000  continue
      do i=1,nv
        if(k(7,i).eq.-2000) k(7,i)=-2001
      enddo

      end


c**************************************************************************

c...Poisson distribution for n>=1
      integer function generate_poisson1(a)
      implicit none
      integer n
      real*8 r,rn,a

        n=0
        r=-log(exp(-a)+rn(0)*(1.0-exp(-a)))
100     n=n+1
c       if(n.gt.100) then
c          print *,'too much jet? a r=',a,r
c          generate_poisson1=n
c          return
c      endif
       r=r-log(max(1d-10,rn(0)))
       if(r.lt.a) go to 100
       generate_poisson1=n

      end

c**************************************************************************
c...Poisson distribution for n>=0
      integer function generate_poisson(a)
      implicit none
      real*8 a, rn,expmean,pir,jamgasdev,rr
      integer N,nmax

      if(a.lt.100)then
        expmean=exp(-a)   !mean of Poisson
        pir = 1
        N = -1
        do while(.true.)
          N=N+1
          pir = pir*rn(0)
          if(pir.le.expmean)goto 9
        enddo
 9            nmax=N
      else
        rr=jamgasdev()
        nmax=nint(a)+sqrt(a)*rr
      endif

      generate_poisson=nmax

      end

c**************************************************************************
      real*8 function jamgasdev()
      implicit none
      INTEGER iset
      real*8 v1,v2,rsq,fac,gset,rn
      SAVE iset,gset
      DATA iset/0/

      if (iset.eq.0) then
1       v1=2.*rn(0)-1.
        v2=2.*rn(0)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.) goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        jamgasdev=v2*fac
        iset=1
      else
        jamgasdev=gset
        iset=0
      endif

      end

c**************************************************************************
      double precision function pdensity(pm,tf,mu,bar,istat)
      implicit none
      integer istat,i,j,nb,opt_integ,bar
      double precision pm,tf,mu,stat,pden,x,dis,p,e,dw
      double precision besk0,besk1,besk2
      parameter(opt_integ=0)

      integer NGP
      parameter(NGP=38)
      real*8 xg(NGP),wg(NGP)
      common/gaussgrid/xg,wg
      real*8 xg1(NGP),wg1(NGP)
      common/gaussgrid1/xg1,wg1

      real*8 spot
      common/mdepot/spot(NGP)

c...Bosons.
      if(istat.eq.-1) then


      if(opt_integ.eq.0) then
      stat=-istat
      pden=0d0
      do j=1,10
        x=j*pm/tf
        besk2 = besk0(x) + 2/x*besk1(x)
        pden=pden+(stat**(j+1)/j)*besk2*exp(j*mu/tf)
      end do
      pdensity=pm**2*tf*pden

      else

c...Numerical integration.
      pdensity=0.0d0
      do i=1, NGP
        p=xg1(i)
        dw=wg1(i)
        e=sqrt(p*p+pm**2)
        dis=exp(-(e-mu)/tf)
        dis = dis/(1.0+istat*dis)
        pdensity = pdensity + p*p*dw*dis
      end do

      endif

      if(pdensity.lt.0d0) then
       print *,'bosons pdensity<0?',pdensity,pm,tf,mu
      endif

c...Fermions.
      else if(istat.eq.1) then

c...Numerical integration.
      pdensity=0.0d0
      do i=1, NGP
        p=xg1(i)
        dw=wg1(i)
        e=sqrt(p*p+pm**2)+bar*spot(i)
        dis=exp(-(e-mu)/tf)
        dis = dis/(1.0+istat*dis)
        pdensity = pdensity + p*p*dw*dis
      end do

c     stat=-istat
c     nb=3
c     if(pm.le.1.0d0.and.mu.gt.0d0) nb=10
c     pden=0d0
c     do j=1,nb
c       x=j*pm/tf
c       besk2 = besk0(x) + 2/x*besk1(x)
c       pden=pden+(stat**(j+1)/j)*besk2*exp(j*mu/tf)
c     end do
c     pdensity=pm**2*tf*pden

      if(pdensity.lt.0d0) then
       print *,'fermion pdensity<0?',pdensity,pm,tf,mu
      endif
c     print *,'m=',pm,mu,tf,mu/tf
c     print *,'d=',pdensity,pden,abs(pdensity-pden)/pden*100

      else
        x=pm/tf
        besk2 = besk0(x) + 2/x * besk1(x)
        pdensity=pm**2*tf*besk2*exp(mu/tf)
      endif

      end

c************************************************************************
      real*8 function fBulk(t)
      implicit none
      real*8 t,tantt,x,jacob,bp2,exp_,f0,pE
      real*8 lam,bmass,xsig,stat,vsig
      common/pdata/bmass,xsig,lam,stat,vsig

      tantt = tan(t * t)
      x = tantt + xsig
      if (x -lam >= 100.0) then
          fBulk=0.0; 
          return;
      endif
      jacob = 2 * t * (tantt * tantt + 1);
      bp2 = x * x - bmass * bmass;
      exp_ = exp(x-lam);
      f0 = 1.0 / (exp_ + stat);
      pE = x * sqrt(bp2);
      fBulk   = jacob * pE * f0;

      end

c************************************************************************
      real*8  function fSurface(t)
      implicit none
      real*8 t, x, bp2, exp_, f0,f, bp, prod, jacob
      real*8 lam,bmass,xsig,stat,vsig,tantt
      common/pdata/bmass,xsig,lam,stat,vsig

      fSurface=0.0;
      tantt = tan(t * t)
      x = tantt + xsig
      if (x -lam >= 100.0) return;

      jacob = 2 * t * (tantt * tantt + 1);
      bp2 = x * x - bmass * bmass;
      exp_ = exp(x-lam);
      f0 = 1.0 / (exp_ + stat);

      ! (x*vsig-bp)^2
      bp = sqrt(bp2);
      prod = x * vsig - bp;
      fSurface = jacob * f0 * prod * prod;

      end

c************************************************************************
      real*8 function pdensity2(mass,tf, cmu,istat)
      implicit none
      integer NGP2
      parameter(NGP2=50)
      real*8 xg2(NGP2),wg2(NGP2)
      common/gaussgrid2/xg2,wg2
      real*8 mass,tf,cmu
      integer i,istat
      real*8 bmass, xsig, lam,stat,vsig, fBulk
      common/pdata/bmass,xsig,lam,stat,vsig
      real*8 lower,upper,center,dxds,i1,i2
      real*8 pi
      parameter(pi=3.141592653589793d0)

      bmass=mass/tf
      lam=cmu/tf
      stat = istat
      xsig= bmass

      lower  = 0.0;
      upper  = sqrt(pi/2)
      center = 0.5 * (upper + lower);
      dxds   = 0.5 * (upper - lower);
      pdensity2 = 0.0

      do i=1,NGP2
        i1=fBulk(center + dxds * xg2(i))
        i2=fBulk(center - dxds * xg2(i))
        pdensity2 =  pdensity2 + wg2(i) * (i1 + i2);
      end do
      pdensity2  =  pdensity2 * dxds*tf**3;

      end

c**************************************************************************
c...Compute positive part of integral in the Cooper-Frye formula
      double precision function pdensity_pos2(pm,tf,mu,istat,vol,sigma)
      implicit none
      integer istat,i,j,nb,opt_integ
      double precision pm,tf,mu,pden,x,dis,p,e,dw,vol,sigma,vsig
      double precision pmin,pmax,pdensity2,tantt,tsig
      double precision pdensity,a,pd0,pd1,pd,pb,pd3,t,jacob

      real*8 bmass, xsig, lam,stat,fBulk,fSurface
      common/pdata/bmass,xsig,lam,stat,vsig
      real*8 lower,upper,center,dxds,i1,i2

      real*8 pi
      parameter(pi=3.141592653589793d0)

      integer NGP2
      parameter(NGP2=50)
      real*8 xg2(NGP2),wg2(NGP2)
      common/gaussgrid2/xg2,wg2

      integer NGP
      parameter(NGP=38)
      real*8 xg(NGP),wg(NGP)
      common/gaussgrid/xg,wg
      real*8 xg1(NGP),wg1(NGP)
      common/gaussgrid1/xg1,wg1

      integer NGP3
      parameter(NGP3=100)
      real*8 xg3(NGP3),wg3(NGP3)
      common/gaussgrid3/xg3,wg3

      pdensity_pos2=0.0d0
      bmass=pm/tf
      lam=mu/tf
      stat = istat
      xsig= bmass

      lower  = 0.0;
      upper  = sqrt(pi/2)
      center = 0.5 * (upper + lower);
      dxds   = 0.5 * (upper - lower);

      if(vol.gt.0d0) then

      pdensity2 = 0.0
      do i=1,NGP2
        i1=fBulk(center + dxds * xg2(i))
        i2=fBulk(center - dxds * xg2(i))
        pdensity2 =  pdensity2 + wg2(i) * (i1 + i2);
      end do
      pdensity2  =  pdensity2 * dxds*tf**3;
      pdensity_pos2 = vol*pdensity2

      pd0=pdensity2

      endif

      vsig=abs(vol/sigma)
      if(vsig.ge.1.0d0) return

      xsig=bmass/sqrt(1.0d0-vsig**2)
c     tsig=sqrt(atan(xsig-bmass))

c     if(tsig .le. upper*(2.3d0-1d-16)) then
c     tantt=tan(tsig*tsig)
c     xsig=tantt + bmass;
c     vsig = sqrt(1.0d0 - (bmass/xsig)**2)

      pdensity2 = 0.0
      do i=1,NGP2
        i1=fSurface(center + dxds * xg2(i))
        i2=fSurface(center - dxds * xg2(i))
        pdensity2 =  pdensity2 + wg2(i) * (i1 + i2);
      end do
      pdensity2  =  pdensity2 * dxds*tf**3;
      pdensity_pos2 = pdensity_pos2 + sigma/4*pdensity2

c     endif

c...Numerical integration.
c     pdensity=0.0d0
      pmin=pm*vsig/sqrt(1.0d0-vsig*vsig)

c     pmax=100*tf
c     pd=0d0
c     do i=1, NGP
c       p=xg1(i)
c       dw=wg1(i)
c       e=sqrt(p*p+pm**2)
c       dis=exp(-(e-mu)/tf)
c       dis = dis/(1.0+istat*dis)
c       pdensity = pdensity + (e*vsig - p)**2*p/e*dw*dis

c       p=xg(i)*pmin
c       dw=wg(i)*pmin
c       e=sqrt(p*p+pm**2)
c       dis=exp(-(e-mu)/tf)
c       dis = dis/(1.0+istat*dis)
c       pd = pd + (e*vsig - p)**2*p/e*dw*dis
c     end do

c     a=1.0d0
c     do i=1, NGP3
c       p=(pmin+a*xg3(i))/(1-xg3(i))
c       dw=(pmin+a)/(1-xg3(i))**2*wg3(i)
c       e=sqrt(p*p+pm**2)
c       dis=exp(-(e-mu)/tf)
c       dis = dis/(1.0+istat*dis)
c       pdensity = pdensity + (e*vsig - p)**2*p/e*dw*dis
c     end do

c     xsig=bmass/sqrt(1.0d0-vsig*vsig)
c     pd=0.0d0
c     do i=1, NGP3
c       x=(xsig+a*xg3(i))/(1-xg3(i))
c       dw=(xsig+a)/(1-xg3(i))**2*wg3(i)
c       dis=exp(-x+lam)
c       dis = dis/(1.0d0+dis*istat)
c       pb = x*vsig - sqrt(x*x - bmass*bmass)
c       pd = pd + pb**2*dw*dis
c     end do
c     pd=pd*tf**3

c     pd3=0.0
c     do i=1, NGP3
c       t=xg3(i)*upper
c       dw=upper*wg3(i)
c       tantt = tan(t * t)
c       x = tantt + xsig
c       if (x -lam >= 100.0)cycle;
c       jacob = 2 * t * (tantt * tantt + 1);
c       dis = 1.0 / (exp(x-lam) + istat);
c       pb=x*x -bmass*bmass
c       if(pb.lt.0d0) cycle
c       pb = x*vsig - sqrt(x*x - bmass*bmass)
c       pd3 = pd3 + jacob * dis * pb**2 * dw
c     end do
c     pd3=pd3*tf**3





c     pd1=0.0d0
c     do i=1, NGP
c       p=xg1(i)
c       dw=wg1(i)
c       e=sqrt(p*p+pm**2)
c       dis=exp(-(e-mu)/tf)
c       dis = dis/(1.0+istat*dis)
c       pd1= pd1 + p*p*dw*dis
c     end do

c     print *,'p1=',pd0,pd1
c     print *,'p2=',pdensity2,pdensity
c     print *,'p3=',pd,pd3

c     pd=0.0d0
c     do i=1, NGP
c       x=(xsig+a*xg(i))/(1-xg(i))
c       dw=(xsig+a)/(1-xg(i))**2*wg(i)
c       dis=exp(-x+lam)
c       dis = dis/(1.0d0+dis*istat)
c       pb = x*vsig - sqrt(x*x - bmass*bmass)
c       pd = pd + pb**2*dw*dis
c     end do
c     pd=pd*tf**3

c     pd3=0.0
c     do i=1, NGP
c       t=xg(i)*upper
c       dw=upper*wg(i)
c       tantt = tan(t * t)
c       x = tantt + xsig
c       if (x -lam >= 100.0)cycle;
c       jacob = 2 * t * (tantt * tantt + 1);
c       dis = 1.0 / (exp(x-lam) + istat);
c       pb=x*x -bmass*bmass
c       if(pb.lt.0d0) cycle
c       pb = x*vsig - sqrt(x*x - bmass*bmass)
c       pd3 = pd3 + jacob * dis * pb**2 * dw
c     end do
c     pd3=pd3*tf**3

c     print *,'p3=',pd,pd3

      end

c**************************************************************************
c...Compute positive part of integral in the Cooper-Frye formula
      double precision function pdensity_pos(pm,tf,mu,bar,istat,
     & vol,sigma)
      implicit none
      integer istat,i,j,nb,opt_integ,bar
      double precision pm,tf,mu,stat,pden,x,dis,p,e,dw,vol,sigma,vsig
      double precision besk0,besk1,besk2,pdensity
      double precision pmin,pmax,pd,pb
      double precision bmass,lam,xsig
      parameter(opt_integ=0)

      integer NGP
      parameter(NGP=38)
      real*8 xg(NGP),wg(NGP)
      common/gaussgrid/xg,wg
      real*8 xg1(NGP),wg1(NGP)
      common/gaussgrid1/xg1,wg1

      real*8 spot
      common/mdepot/spot(NGP)

      pdensity_pos=0.0d0

      if(vol.gt.0d0) then

c...Bosons.
      if(istat.eq.-1) then


      if(opt_integ.eq.0) then
      stat=-istat
      pden=0d0
      do j=1,10
        x=j*pm/tf
        besk2 = besk0(x) + 2/x*besk1(x)
        pden=pden+(stat**(j+1)/j)*besk2*exp(j*mu/tf)
      end do
      pdensity=pm**2*tf*pden

      else

c...Numerical integration.
      pdensity=0.0d0
      do i=1, NGP
        p=xg1(i)
        dw=wg1(i)
        e=sqrt(p*p+pm**2)
        dis=exp(-(e-mu)/tf)
        dis = dis/(1.0+istat*dis)
        pdensity = pdensity + p*p*dw*dis
      end do

      endif

      if(pdensity.lt.0d0) then
       print *,'bosons pdensity<0?',pdensity,pm,tf,mu
      endif

c...Fermions.
      else if(istat.eq.1) then

c...Numerical integration.
      pdensity=0.0d0
      do i=1, NGP
        p=xg1(i)
        dw=wg1(i)
        e=sqrt(p*p+pm**2) + bar*spot(i)
        dis=exp(-(e-mu)/tf)
        dis = dis/(1.0+istat*dis)
        pdensity = pdensity + p*p*dw*dis
      end do

      if(pdensity.lt.0d0) then
       print *,'fermion pdensity<0?',pdensity,pm,tf,mu
      endif
c     print *,'m=',pm,mu,tf,mu/tf
c     print *,'d=',pdensity,pden,abs(pdensity-pden)/pden*100

      else
        x=pm/tf
        besk2 = besk0(x) + 2/x * besk1(x)
        pdensity=pm**2*tf*besk2*exp(mu/tf)
      endif

      endif

      pdensity_pos = vol*pdensity
      vsig=abs(vol/sigma)
      if(vsig.ge.1.0d0) return

c...Numerical integration.

      bmass=pm/tf
      lam=mu/tf
      xsig=bmass/sqrt(1.0d0-vsig**2)

      pdensity=0.0d0
      do i=1, NGP
        x=(xsig+1.0d0*xg(i))/(1-xg(i))
        dw=(xsig+1.0d0)/(1-xg(i))**2*wg(i)
        dis=exp(-x+lam)
        dis = dis/(1.0d0+dis*istat)
        pb = x*vsig - sqrt(x*x - bmass*bmass)
        pdensity = pdensity + pb**2*dw*dis
      end do
      pdensity=pdensity*tf**3

      pdensity_pos = pdensity_pos + sigma*pdensity/4

c      print *,'p=',pdensity_pos,sigma*pdensity/4

      end

c**************************************************************************
      real*8 function multlookup(t,mub,mus)
      implicit none
      integer ntmp,nmub,nmus
      parameter(ntmp=170,nmub=60,nmus=40)
      real*8 ptable,tmin,mubmin,musmin,dtmp,dmub,dmus,mult
      common/multtable/ptable(0:ntmp,0:nmub,0:nmus),tmin,mubmin,musmin,
     & dtmp,dmub,dmus
      real*8 t,mub,mus,x1,y1,z1,x2,y2,z2
      integer i,j,k

      i=int((t-tmin)/dtmp)
      j=int((mub-mubmin)/dmub)
      k=int((mus-musmin)/dmus)

      if(i.lt.0.or.i.ge.ntmp) then
       write(6,*)'T out of range',t,i
       multlookup=mult(t,mub,mus)
       return
      endif
      if(j.lt.0.or.j.ge.nmub) then
       write(6,*)'mu_B out of range',mub,j
       multlookup=mult(t,mub,mus)
       return
      endif
      if(k.lt.0.or.k.ge.nmus) then
       write(6,*)'mu_S out of range',mus,k
       multlookup=mult(t,mub,mus)
       return
      endif

      x2=(t-(tmin+i*dtmp))/dtmp
      y2=(mub-(mubmin+j*dmub))/dmub
      z2=(mus-(musmin+k*dmus))/dmus
      x1=1d0-x2
      y1=1d0-y2
      z1=1d0-z2

      multlookup=x1*y1*z1*ptable(i,j,k)
     &       +x2*y1*z1*ptable(i+1,j,k)
     &       +x1*y2*z1*ptable(i,j+1,k)
     &       +x1*y1*z2*ptable(i,j,k+1)
     &       +x2*y1*z2*ptable(i+1,j,k+1)
     &       +x1*y2*z2*ptable(i,j+1,k+1)
     &       +x2*y2*z1*ptable(i+1,j+1,k)
     &       +x2*y2*z2*ptable(i+1,j+1,k+1)

      if(multlookup.lt.0d0) then
        multlookup=0d0
        print *,'multlookup<0?',multlookup,i,j,k
        print *,x1,y1,z1
        print *,x2,y2,z2
        print *,t,mub,mus
        print *,ptable(i,j,k)
        print *,ptable(i+1,j,k)
        print *,ptable(i,j+1,k)
        print *,ptable(i,j,k+1)
        print *,ptable(i+1,j+1,k)
        print *,ptable(i+1,j,k+1)
        print *,ptable(i,j+1,k+1)
        print *,ptable(i+1,j+1,k+1)
      endif

      end

c**************************************************************************
      subroutine read_frezeout_mult_table(fname)
      implicit none
      character  fname*(*)
      integer ntmp,nmub,nmus
      integer ntmp1,nmub1,nmus1
      parameter(ntmp=170,nmub=60,nmus=40)
      real*8 ptable,tmin,mubmin,musmin,dtmp,dmub,dmus
      common/multtable/ptable(0:ntmp,0:nmub,0:nmus),tmin,mubmin,musmin,
     & dtmp,dmub,dmus
      integer i,j,k
      real*8 tf,mub,mus,a

      open(80,file=fname,status='old')
      tmin=0.003d0
      mubmin=0.0d0
      musmin=0.0d0

      dtmp=0.001d0
      dmub=0.02d0
      dmus=0.02d0
      read(80,*)ntmp1,nmub1,nmus1
      if(ntmp1.ne.ntmp.or.nmub1.ne.nmub.or.nmus1.ne.nmus) then
        print *,'ntmp size is wrong',ntmp1,ntmp
        print *,'nmub',nmub1,nmub
        print *,'nmus',nmus1,nmus
        stop
      endif
      read(80,*)tmin,dtmp
      read(80,*)mubmin,dmub
      read(80,*)musmin,dmus


      do j=0,nmub
      do k=0,nmus
      do i=0,ntmp
        read(80,*)tf,mub,mus,a
        ptable(i,j,k)=a
        if(a.lt.0d0) then
         print *,'freezeout table mult<0?',a,tf,mub,mus
         print *,'i j k',i,j,k
         stop
        endif
      end do
      end do
      end do

      close(80)

      end

c**************************************************************************
      subroutine make_mult_table_file

      implicit none
      integer i,j,k
      real*8 mult,a,b,tf,bmu,smu
      integer ntmp,nmub,nmus
      parameter(ntmp=170,nmub=60,nmus=40)
      real*8 ptable,tmin,mubmin,musmin,dtmp,dmub,dmus
      common/multtable/ptable(0:ntmp,0:nmub,0:nmus),tmin,mubmin,musmin,
     & dtmp,dmub,dmus

      tmin=0.003d0
      mubmin=0.0d0
      musmin=0.0d0

      dtmp=0.001d0
      dmub=0.02d0
      dmus=0.02d0

      write(10,800)ntmp,nmub,nmus
      write(10,810)tmin,dtmp
      write(10,810)mubmin,dmub
      write(10,810)musmin,dmus

       do j=0,nmub
         bmu=mubmin+j*dmub
       do k=0,nmus
         smu=musmin+k*dmus
      do i=0,ntmp
       tf=tmin+i*dtmp
         if((bmu+smu)/tf.lt.50.0d0) then
           a=mult(tf,bmu,smu)
         else
           a=0d0
         endif
c     write(6,820)tf,bmu,smu,a
      write(10,820)tf,bmu,smu,a
      end do
      end do
      end do
 800  format(3(i4,1x))
 810  format(3(e15.8,1x))
 820  format(4(e15.8,1x))

      end

c**************************************************************************
      subroutine make_mult_table
      implicit none
      integer ntmp,nmub,nmus
      parameter(ntmp=170,nmub=60,nmus=40)
      real*8 ptable,tmin,mubmin,musmin,dtmp,dmub,dmus
      common/multtable/ptable(0:ntmp,0:nmub,0:nmus),tmin,mubmin,musmin,
     & dtmp,dmub,dmus
      integer i,j,k
      real*8 t,mub,mus,tf,mult

      tmin=0.003d0
      mubmin=0.0d0
      musmin=0.0d0

      dtmp=0.001d0
      dmub=0.02d0
      dmus=0.02d0
c       tf=0.166 - 0.14*mub**2-0.53*mub**4

      do i=0,ntmp
         t=tmin+i*dtmp
      do j=0,nmub
        mub=mubmin+j*dmub
      do k=0,nmus
        mus=musmin+k*dmus
        ptable(i,j,k)=0d0
        if((mub+mus)/t.lt.40.0d0) then
          ptable(i,j,k)=mult(t,mub,mus)
        endif
      end do
      end do
      end do

      end

c**************************************************************************
      real*8 function smult(i,tch,bmu,smu)
      implicit none
      include 'jam2.inc'
      integer i,kc,kf,ns,ibary,istat,spin
      real*8 cmu,smu,bmu,tch,pm,pdensity
      integer kcmin,kcmax
      parameter(kcmin=132,kcmax=366)

      smult=0d0

      if(i.eq.0) return
      kc=abs(i)
      if(kchg(kc,3).eq.0.and.i.lt.0) return
      if(kc.lt.kcmin.or.kc.gt.kcmax) return
      if(kc.ge.195.and.kc.le.251) return
      kf=kchg(kc,4)

      if(kf.eq.130.or.kf.eq.310) return
c....Exclude heavy flavour hadrons.
      if(mod(kf/1000,10).ge.4) return
      if(mod(kf/100,10).ge.4) return
      if(mod(kf/10,10).ge.4) return

      spin=max(1,mod(kf,10))
      ns=kchg(kc,7)
      ibary=kchg(kc,6)/3
      cmu = (ns*smu + ibary*bmu)*isign(1,i)
      pm=pmas(kc,1)

      istat=1
      if(ibary.eq.0) istat=-1

c     if(mstc(146).eq.1) then
c...Boltzmann approximation except pions.
c       istat=0
c       if(abs(kf).eq.211.or.kf.eq.111) istat=-1
c     endif

c....Get hadron density.
      smult=pdensity(pm,tch,cmu,ibary,istat)
     &  *spin/(2*paru(1)**2*paru(3)**3)

      end 


c**************************************************************************
      real*8 function mult(tch,bmu,smu)
      implicit none
      include 'jam2.inc'
      integer i,kcmin,kcmax,kc,kf,ispin,ns,ibary,istat,isp
      real*8 pmult,cmu,smu,bmu,pm,dns,pdensity,tch,fac
      common/bolz3/pmult(-500:500)
      real*8 scsum
      integer kcsum,maxkc
      common/bolz4/maxkc,kcsum(0:500),scsum(0:500)

      integer NGP
      parameter(NGP=38)
      real*8 spot
      common/mdepot/spot(NGP)

c...Compute momentum dependent potential.
      spot=0.0
      if(mstd(101).eq.1) call momdeppot(tch,bmu,smu)

      fac=1.0d0/(2*paru(1)**2*paru(3)**3)

c...Initialize particle multiplicities.
      pmult=0    ! hadron density.
      scsum=0    ! cumulative sum for each species
      kcsum=0
      mult=0d0
      maxkc=0

      kcmin=132
      kcmax=366
      isp=0
c...Loop over particles.
      do kc=kcmin,kcmax

c....Exclude heavy flavour hadrons.
        if(kc.ge.195.and.kc.le.251) cycle

        kf=kchg(kc,4)
        if(kf.eq.130.or.kf.eq.310) cycle ! K_L0 or K_S0

        ispin=max(1,mod(kf,10)) ! spin 2J+1
        ns=kchg(kc,7)*isign(1,kf)
        ibary=kchg(kc,6)/3               ! baryon number times 3.
        istat=1
        if(ibary.eq.0) istat=-1
        cmu = ns*smu + ibary*bmu
        pm=pmas(kc,1)                 ! particle mass
c...Neglect mass integral
c       width=pmas(kc,2)                 ! particle width

c...Boltzmann approximation except pions.
c       if(mstc(146).eq.1) then
c         istat=0
c         if(kf.eq.211.or.kf.eq.111) istat=-1
c       endif

c....Get hadron density.
        dns=pdensity(pm,tch,cmu,ibary,istat)*ispin*fac
        pmult(kc) = dns
        mult=mult+dns
        isp=isp+1
        scsum(isp)=scsum(isp-1)+dns
        kcsum(isp)=kc
        maxkc=isp

c......Anti-particle if it exists.
        if(kchg(kc,3).eq.1) then
          dns=pdensity(pm,tch,-cmu,-ibary,istat)*ispin*fac
          pmult(-kc) = dns
          mult=mult+dns
          isp=isp+1
          scsum(isp)=scsum(isp-1)+dns
          kcsum(isp)=-kc
          maxkc=isp
        endif

      end do

c     if(mult.ne.scsum(maxkc)) then
c     print *,'scsum(maxkc)=',scsum(maxkc),mult
c     stop
c     endif

c     print *,'isp=',isp

      end 

c**************************************************************************
      real*8 function mult_pos(tch,bmu,smu,vol,sigma)
      implicit none
      include 'jam2.inc'
      integer i,kcmin,kcmax,kc,kf,ispin,ns,ibary,istat
      real*8 pmult,cmu,smu,bmu,pm,dns,pdensity_pos,tch,fac,vol,sigma
      common/bolz3/pmult(-500:500)

      integer NGP
      parameter(NGP=38)
      real*8 spot
      common/mdepot/spot(NGP)

c...Compute momentum dependent potential.
      spot=0.0
      if(mstd(101).eq.1) call momdeppot(tch,bmu,smu)

      fac=1.0d0/(2*paru(1)**2*paru(3)**3)

c...Initialize particle multiplicities.
      do i=-500,500
        pmult(i)=0    ! hadron density.
      end do
      mult_pos=0d0

      kcmin=132
      kcmax=366
c...Loop over particles.
      do kc=kcmin,kcmax

c....Exclude heavy flavour hadrons.
        if(kc.ge.195.and.kc.le.251) cycle
        kf=kchg(kc,4)
        if(kf.eq.130.or.kf.eq.310) cycle ! K_L0 or K_S0

        ispin=max(1,mod(kf,10)) ! spin 2J+1
        ns=kchg(kc,7)*isign(1,kf)
        ibary=kchg(kc,6)/3               ! baryon number times 3.
        istat=1
        if(ibary.eq.0) istat=-1
        cmu = ns*smu + ibary*bmu
        pm=pmas(kc,1)                 ! particle mass
c...Neglect mass integral
c       width=pmas(kc,2)                 ! particle width

c...Boltzmann approximation except pions.
c       if(mstc(146).eq.1) then
c         istat=0
c         if(kf.eq.211.or.kf.eq.111) istat=-1
c       endif

c....Get hadron density.
        dns=pdensity_pos(pm,tch,cmu,ibary,istat,vol,sigma)*ispin*fac
        pmult(kc) = dns
        mult_pos=mult_pos+dns

c......Anti-particle if it exists.
        if(kchg(kc,3).eq.1) then
          dns=pdensity_pos(pm,tch,-cmu,ibary,istat,vol,sigma)*ispin*fac
          pmult(-kc) = dns
          mult_pos=mult_pos+dns
        endif

      end do

      end 

c**************************************************************************
      real*8 function mult2(tch,bmu,smu)
      implicit none
      include 'jam2.inc'
      integer i,kcmin,kcmax,ibary,istat,kc,kf,ispin,ns,jamflav
      real*8 tch,bmu,smu,fac,cmu,pm,dns,pdensity
      real*8 pmult,m
      common/bolz3/pmult(-500:500)

      fac=1.0d0/(2*paru(1)**2*paru(3)**3)

c...Initialize particle multiplicities.
      do i=-500,500
        pmult(i)=0    ! hadron density.
      end do
      mult2=0d0

      ibary=0
      istat=-1
      kcmin=132 ! pi0
      kcmax=194 ! K_2(1820)+
c     kcmax=366
c...Loop over mesons.
      do kc=kcmin,kcmax
        kf=kchg(kc,4)
        if(kf.eq.130.or.kf.eq.310) cycle ! K_L0 or K_S0
        if(kc.ge.195.and.kc.le.251) cycle
        ispin=mod(kf,10) ! spin 2J+1
        if(kf.eq.10220) ispin=1  ! sigma meson
c       if(ispin.eq.0) goto 1000
        ibary=kchg(kc,6)/3               ! baryon number times 3.
c       ns=jamflav(kf,3)                 ! net s
        ns=kchg(kc,7)
        cmu = ns*smu + ibary*bmu
        pm=pmas(kc,1)                 ! particle mass

c...Boltzmann approximation except pions.
c       if(mstc(146).eq.1) then
c         istat=0
c         if(kf.eq.211.or.kf.eq.111) istat=-1
c       endif

c....Get hadron density.
        dns=pdensity(pm,tch,cmu,ibary,istat)*ispin*fac
        pmult(kc) = dns
        mult2=mult2+dns

c......Anti-particle if it exists.
        if(kchg(kc,3).eq.1) then
          dns=pdensity(pm,tch,-cmu,ibary,istat)*ispin*fac
          pmult(-kc) = dns
          mult2=mult2+dns
        endif

 1000 end do


c...N/N*
      call baryon_mult(m,tch,smu,bmu,fac,252,252,1,0) ! n0
      mult2=mult2+m
      call baryon_mult(m,tch,smu,bmu,fac,253,253,1,0) ! n0
      mult2=mult2+m

      kcmin=mstc(22)
      kcmax=mstc(23)
c     kcmin=252
c     kcmax=272
      call baryon_mult(m,tch,smu,bmu,fac,kcmin,kcmax,2,0)
      mult2=mult2+m

c...Delta/Delta*
      call baryon_mult(m,tch,smu,bmu,fac,274,274,1,0)
      mult2=mult2+m
      call baryon_mult(m,tch,smu,bmu,fac,275,275,1,0)
      mult2=mult2+m
      call baryon_mult(m,tch,smu,bmu,fac,276,276,1,0)
      mult2=mult2+m
      call baryon_mult(m,tch,smu,bmu,fac,277,277,1,0)
      mult2=mult2+m

      kcmin=mstc(24)
      kcmax=mstc(25)
c     kcmin=274
c     kcmax=310
      callbaryon_mult(m,tch,smu,bmu,fac,kcmin,kcmax,4,0)
      mult2=mult2+m

c...L/L*
      kcmin=mstc(26)-1
      kcmax=mstc(27)
      kcmin=314
      kcmax=326
      call baryon_mult(m,tch,smu,bmu,fac,kcmin,kcmax,1,-1)
      mult2=mult2+m

c...Sigma
      call baryon_mult(m,tch,smu,bmu,fac,327,327,1,-1)
      mult2=mult2+m
      call baryon_mult(m,tch,smu,bmu,fac,328,328,1,-1)
      mult2=mult2+m
      call baryon_mult(m,tch,smu,bmu,fac,329,329,1,-1)
      mult2=mult2+m
      call baryon_mult(m,tch,smu,bmu,fac,330,330,1,-1)
      mult2=mult2+m
      call baryon_mult(m,tch,smu,bmu,fac,331,331,1,-1)
      mult2=mult2+m
      call baryon_mult(m,tch,smu,bmu,fac,332,332,1,-1)
      mult2=mult2+m

c...Sigma*
c     kcmin=327  !mstc(28)
c     kcmax=mstc(29)
      kcmin=333
      kcmax=351
      call baryon_mult(m,tch,smu,bmu,fac,kcmin,kcmax,3,-1)
      mult2=mult2+m
  
c...Xi
      call baryon_mult(m,tch,smu,bmu,fac,354,354,1,-2)
      mult2=mult2+m
      call baryon_mult(m,tch,smu,bmu,fac,355,355,1,-2)
      mult2=mult2+m
      call baryon_mult(m,tch,smu,bmu,fac,356,356,1,-2)
      mult2=mult2+m
      call baryon_mult(m,tch,smu,bmu,fac,357,357,1,-2)
      mult2=mult2+m

c...Xi*
      kcmin=mstc(30)
      kcmax=mstc(31)
      kcmin=358
      call baryon_mult(m,tch,smu,bmu,fac,kcmin,kcmax,2,-2)
      mult2=mult2+m

c...Omega
      call baryon_mult(m,tch,smu,bmu,fac,366,366,1,-3)
      mult2=mult2+m


      end 

c*******************************************************************
      subroutine baryon_mult(bm,tch,smu,bmu,fac,kcmin,kcmax,isp,ns)
      implicit none
      integer istat,ibary,kc,j,isp,ns,kcmin,kcmax,kf,ispin
      real*8 pm,dns,pdensity,pmult,bmu,smu,cmu,tch,fac,bm
      include 'jam2.inc'
      common/bolz3/pmult(-500:500)

      istat=1  ! Fermi distribution.
c...Boltzmann approximation except pions.
c     if(mstc(146).eq.1) istat=0

      bm=0d0
      ibary=1
      do kc=kcmin,kcmax,isp
        kf=kchg(kc,4)
        ispin=mod(kf,10) ! spin 2J+1
        cmu = ns*smu + ibary*bmu
        pm=pmas(kc,1)
        dns=pdensity(pm,tch,cmu,ibary,istat)*ispin*fac
        do j=0,isp-1
          pmult(kc+j) = dns
        end do
        bm=bm+dns*isp

c......Anti-particle.
        dns=pdensity(pm,tch,-cmu,ibary,istat)*ispin*fac
        do j=0,isp-1
          pmult(-kc-j)=dns
        end do
        bm=bm+dns*isp
      end do


      end
c*******************************************************************
      subroutine cooperfrye(dst,dsx,dsy,dsz,vx,vy,vz,tf,mu,istat,p)
      implicit none
      integer istat,ntry
      real*8 dst,dsx,dsy,dsz,vx,vy,vz,tf,mu,p(5),prob
      real*8 rfac,gam,vol,pm,dstr,dsxr,dsyr,dszr,pp,cost,sint,er,pds
      real*8 thermaldist3,phi,pi,rn
      parameter(pi=3.141592653589793d0)
      real*8 thermaldist0

      gam = 1.0/sqrt(1.0-vx*vx-vy*vy-vz*vz)
      vol=gam*(dst + dsx*vx + dsy*vy + dsz*vz)

c     rfac=7.0d0
c     if(dsx*dsy*dsz.eq.0d0) rfac=2.0d0
c     rfac=rfac*vol
      rfac=vol+sqrt(vol**2 - (dst**2 - dsx**2-dsy**2-dsz**2))

      ntry=0
 100  continue
      ntry=ntry+1
      if(ntry.ge.100) then
        print *,'cooperfrye::infinit loop?',ntry
        return
      endif

      pp=thermaldist3(p(5),tf,mu,istat)
c     pp=thermaldist0(p(5),tf,mu,istat)

      cost=2d0*rn(0)-1d0
      sint=sqrt(1d0-cost**2)
      phi=2*pi*rn(0)
      p(1)=pp*sint*cos(phi)
      p(2)=pp*sint*sin(phi)
      p(3)=pp*cost
      p(4)=sqrt(pp**2+p(5)**2)
      er=p(4)
      call jamrobo(0d0,0d0,vx,vy,vz,gam,p(1),p(2),p(3),p(4))

      pds = p(4)*dst+p(1)*dsx+p(2)*dsy+p(3)*dsz
c     pds=abs(pds)
      if(pds.le.0d0) goto 100
      prob=pds/(er*rfac)

      if(prob.lt.0.0.or.prob.gt.1.0) then
        print *,'prob?',prob,rfac
        rfac=rfac*1.2d0
        goto 100
      endif

      if(rn(0).gt.prob) goto 100

      end

c*******************************************************************
      subroutine cooperfrye2(dst,dsx,dsy,dsz,vx,vy,vz,tf,mu,istat,p)
      implicit none
      integer istat,ntry
      real*8 dst,dsx,dsy,dsz,vx,vy,vz,tf,mu,p(5),prob
      real*8 rfac,gam,vol,dstr,dsxr,dsyr,dszr,pp,cost,sint,er,pds
      real*8 thermaldist3,phi,pi,rn
      parameter(pi=3.141592653589793d0)

      rfac=3.5d0
      gam = 1.0/sqrt(1.0-vx*vx-vy*vy-vz*vz)
      vol=gam*(dst + dsx*vx + dsy*vy + dsz*vz)
      dstr=dst
      dsxr=dsx
      dsyr=dsy
      dszr=dsz
      call jamrobo(0d0,0d0,vx,vy,vz,gam,dsxr,dsyr,dszr,dstr)

      ntry=0
100   continue
      ntry=ntry+1
      if(ntry.ge.200) then
        print *,'cooperfrye::infinit loop?',ntry
        goto  200
      endif
      pp=thermaldist3(p(5),tf,mu,istat)
      cost=2d0*rn(0)-1d0
      sint=sqrt(1d0-cost**2)
      phi=2*pi*rn(0)
      p(1)=pp*sint*cos(phi)
      p(2)=pp*sint*sin(phi)
      p(3)=pp*cost
      p(4)=sqrt(pp**2+p(5)**2)
      er=p(4)
      pds = p(4)*dstr+p(1)*dsxr+p(2)*dsyr+p(3)*dszr
c     call jamrobo(0d0,0d0,vx,vy,vz,gam,p(1),p(2),p(3),p(4))
c     pds = p(4)*dst+p(1)*dsx+p(2)*dsy+p(3)*dsz
c     pds=abs(pds)
      if(pds.le.0d0) goto 100
      prob=pds/(er*vol*rfac)
      if(prob.lt.0.0.or.prob.gt.1.0) then
        print *,'prob?',prob,rfac
        rfac=rfac*1.2
        goto 100
      endif

      if(rn(0).gt.prob) goto 100
200   call jamrobo(0d0,0d0,vx,vy,vz,gam,p(1),p(2),p(3),p(4))

      end

c******************************************************************
      real*8 function thermaldist3(pm,tf,mu,a)

c....Generate the momentum from the thermal distribution.
      implicit none
      integer a
      real*8 rn,p,tf,pm,mu,w0,w,fmax,x,e,pmax,emax
      real*8 m,xm,xm1,pmaxsq,beta,g,gfac,gmax,ex
      real*8 spotmom

      if(pm.eq.0d0) then
        print *,'thermaldist3 m=',pm
        stop
      endif

      m=pm/tf
      xm=mu/tf
      w0=1.0d0
      if(a.eq.-1) w0 = 1.0 - exp(-m+xm)
      pmaxsq = 2 + 2*sqrt(1.0+m*m)      ! Peak mass squared
      emax=sqrt(pmaxsq+m*m)
      pmax=sqrt(pmaxsq)

c....Log-Logistic with alpha=2: g(x)=2*(x/b)/(b*(1+(x/b)**2)**2)
      beta=pmax*sqrt(3.)          ! Peak is chosen to be the same as MB
      gmax=9./(beta*8*sqrt(3.0))
      gfac=1.0/gmax*2/beta**2     ! Normalization factor to one

 1000 continue

c...Uniformly distributed.
c2000  p=30*rn(0)
c      g=1.0

 2000  p=rn(0)
       p=beta*(p/(1.0-p))**0.5d0
       g=gfac*p/(1.0+(p/beta)**2)**2

      e=sqrt(p*p+m*m)
      ex=2*log(p/pmax)-e+emax
      if(g*rn(0).gt.exp(max(-50d0,min(50d0,ex)))) goto 2000
c     if(g*rn(0).gt.exp(2*log(p/pmax)-e+emax)) goto 2000

      xm1=xm
c2018/10/27
c....Add momentum dependent potential for baryons
      if(a.eq.1)  xm1 = xm - spotmom(p*tf)/tf

      w=1.0 + a*exp(-e+xm1)
      if(rn(0).gt.w0/w) goto 1000

      thermaldist3=p*tf

      end

c******************************************************************
      real*8 function thermaldist(pm,tf,mu,a)
      implicit none
      integer a
      real*8 rn,p,tf,pm,mu,w0,w

      w0=0.0d0
      if(a.eq.-1) w0 = exp(-(pm-mu)/tf)

c1000 e=2.0*rn(0)**(1.0/3.0)
c     if(rn(0).gt.1.0/(exp(sqrt(e**2+pm**2)/tf)+1)/fmax) goto 1000

c...First Generate Boltzmann distribution.
c...B. Zhang, M. Gyulassy and Y. Pang, Phys. Rev. C58,1175 (1998).
 1000 p=rn(0)
      p=p*rn(0)
      p=p*rn(0)
      p=-tf*log(p)
      if(rn(0).gt.exp((p-sqrt(p**2+pm**2))/tf)) goto 1000

      w=a*exp(-(sqrt(p**2+pm**2)-mu)/tf)
      if(rn(0).gt.(1.0-w0)/(1.0+w)) goto 1000
      thermaldist=p

      end

c*************************************************************************
c...Compute single particle potential
      real*8 function sp_potential(rhob)
      implicit none
      include 'jam2.inc'
      real*8 rhob

      sp_potential=0.0d0

      if(pard(142).eq.0d0) return

c...EoS-Q
      if(pard(142).gt.0.d0) then
        sp_potential=pard(142)*rhob
        return
      endif

c...QMD potential.
      sp_potential= pard(101)*(rhob/parc(21))
     &             +pard(102)*(rhob/parc(21))**pard(103)

c     bmu = bmu - pard(142)*frbdn(n)

      end
c*************************************************************************
      block data gausspoint

      implicit none
      integer NGP
      parameter(NGP=38)
      real*8 xg(NGP),wg(NGP)
      common/gaussgrid/xg,wg
      real*8 xg1(NGP),wg1(NGP)
      common/gaussgrid1/xg1,wg1

      integer NGP2
      parameter(NGP2=50)
      real*8 xg2(NGP2),wg2(NGP2)
      common/gaussgrid2/xg2,wg2

      integer NGP3
      parameter(NGP3=100)
      real*8 xg3(NGP3),wg3(NGP3)
      common/gaussgrid3/xg3,wg3

c...Useage: x  = xmin + xg[i]*(xmax-xmin)
c...        dx =        wg[i]*(xmax-xmin)
      data xg/9.750347321562d-04,5.130272866807d-03,1.257683570492d-02,
     & 2.326683453324d-02, 3.712933397571d-02, 5.407213049768d-02,
     & 7.398248903382d-02, 9.672791619734d-02, 1.221570481230d-01,
     & 1.501006598104d-01, 1.803727920852d-01, 2.127719894761d-01,
     & 2.470826410360d-01, 2.830764152838d-01, 3.205137797603d-01,
     & 3.591455951049d-01, 3.987147730539d-01, 4.389579873311d-01,
     & 4.796074260477d-01, 5.203925739523d-01, 5.610420126689d-01,
     & 6.012852269461d-01, 6.408544048951d-01, 6.794862202397d-01,
     & 7.169235847162d-01, 7.529173589640d-01, 7.872280105239d-01,
     & 8.196272079148d-01, 8.498993401896d-01, 8.778429518770d-01,
     & 9.032720838027d-01, 9.260175109662d-01, 9.459278695023d-01,
     & 9.628706660243d-01, 9.767331654668d-01, 9.874231642951d-01,
     & 9.948697271332d-01, 9.990249652678d-01/

      data wg/2.501440374819d-03,5.806722358234d-03,9.078288854807d-03,
     & 1.228986986912d-02, 1.541975027258d-02, 1.844704079701d-02,
     & 2.135157925234d-02, 2.411403093038d-02, 2.671600995517d-02,
     & 2.914019957350d-02, 3.137046669607d-02, 3.339196898957d-02,
     & 3.519125353345d-02, 3.675634629237d-02, 3.807683177422d-02,
     & 3.914392232911d-02, 3.995051662176d-02, 4.049124688530d-02,
     & 4.076251464019d-02, 4.076251464019d-02, 4.049124688530d-02,
     & 3.995051662176d-02, 3.914392232911d-02, 3.807683177422d-02,
     & 3.675634629237d-02, 3.519125353345d-02, 3.339196898957d-02,
     & 3.137046669607d-02, 2.914019957350d-02, 2.671600995517d-02,
     & 2.411403093038d-02, 2.135157925234d-02, 1.844704079701d-02,
     & 1.541975027258d-02, 1.228986986912d-02, 9.078288854807d-03,
     & 5.806722358234d-03, 2.501440374819d-03/

c...Integration from 0 to infinity
      data xg1/ 9.7598635274825029E-004, 5.1567282900350925E-003,
     &  1.2737027203429915E-002, 2.3821075556614637E-002,
     &  3.8561081239514093E-002, 5.7163058876923696E-002,
     &  7.9893185774238318E-002, 0.10708613487768781,
     &  0.13915592517069370, 0.17660992627306188,
     &  0.22006686740430886, 0.27028000354623533,
     &  0.32816701341037852, 0.39484879744313967,
     &  0.47170019084008463, 0.56041683159486888,
     &  0.66310422273140357, 0.78239771250444978,
     &  0.92162619155991743, 1.0850386080146326,
     &   1.2781223462412827, 1.5080585611125856,
     &   1.7843860919632593, 2.1199906623294522,
     &   2.5326150325784020, 3.0472288777832852,
     &   3.6998667562505618, 4.5440734073012035,
     &   5.6621958974937234, 7.1861834037850976,
     &   9.3382770901404282, 12.516712036315512,
     &   17.493815405383994, 25.932882788963024,
     &   41.979632599852103, 78.511255729336369,
     $   193.92140592949386, 1024.6044908149897/
      data wg1/ 2.5063255009009689E-003 , 5.8667641486190823E-003 ,
     &   9.3110224667196529E-003 , 1.2882359515084497E-002 ,
     &   1.6631883263887846E-002 , 2.0616297188131877E-002 ,
     &   2.4899556074333709E-002 , 2.9555113877367468E-002 ,
     &   3.4668730861567069E-002 , 4.0342010459859773E-002 ,
     &   4.6696921065189613E-002 , 5.3881658070395952E-002 ,
     &   6.2078343057469607E-002 , 7.1513263780566799E-002 ,
     &   8.2470665216640246E-002 , 9.5311563422862275E-002 ,
     &  0.11049975937612265      , 0.12863832689014129      ,
     &  0.15052158636994833      , 0.17721038534567848      ,
     &  0.21014315041270704      , 0.25130304168722140      ,
     &  0.30347523233916146      , 0.37065289260111578      ,
     &  0.45869600731274773      , 0.57643490027775091      ,
     &  0.73758677242364035      , 0.96422619051202441      ,
     &   1.2933835087462366      , 1.7903359705458926      ,
     &   2.5773069794025916      , 3.9009656479992461      ,
     &   6.3092791820017702      , 11.185181156179484      ,
     &   22.702447594239192      , 57.393303343920699      ,
     &   220.62266770103969      , 2631.1765081920767      /

      data wg2/
     & 3.125542345386416E-002, 3.122488425484931E-002,
     & 3.116383569621074E-002, 3.107233742756551E-002,
     & 3.095047885049043E-002, 3.079837903115330E-002,
     & 3.061618658397945E-002, 3.040407952645461E-002,
     & 3.016226510516754E-002, 2.989097959333334E-002,
     & 2.959048805991241E-002, 2.926108411063884E-002,
     & 2.890308960112587E-002, 2.851685432239523E-002,
     & 2.810275565910164E-002, 2.766119822079209E-002,
     & 2.719261344657662E-002, 2.669745918357103E-002,
     & 2.617621923954466E-002, 2.562940291020853E-002,
     & 2.505754448157973E-002, 2.446120270795577E-002,
     & 2.384096026596899E-002, 2.319742318525421E-002,
     & 2.253122025633597E-002,2.184300241624767E-002,
     & 2.113344211252747E-002,2.040323264620976E-002,
     & 1.965308749443533E-002,1.888373961337431E-002,
     & 1.809594072212793E-002,1.729046056832494E-002,
     & 1.646808617614444E-002,1.562962107754551E-002,
     & 1.477588452744186E-002,1.390771070371846E-002,
     & 1.302594789297084E-002,1.213145766298006E-002,
     & 1.122511402318562E-002,1.030780257486862E-002,
     & 9.380419653694766E-003,8.443871469668482E-003,
     & 7.499073255464270E-003,6.546948450845874E-003,
     & 5.588428003865582E-003,4.624450063422539E-003,
     & 3.655961201326225E-003,2.683925371553286E-003,
     & 1.709392653518899E-003,7.346344905058351E-004/

      data xg2/
     &1.562898442154312E-002, 4.687168242159179E-002,
     &7.806858281343659E-002, 1.091892035800607E-001,
     &1.402031372361136E-001, 1.710800805386032E-001,
     &2.017898640957357E-001, 2.323024818449742E-001,
     &2.625881203715034E-001, 2.926171880384720E-001,
     &3.223603439005293E-001, 3.517885263724220E-001,
     &3.808729816246301E-001, 4.095852916783023E-001,
     &4.378974021720319E-001, 4.657816497733580E-001,
     &4.932107892081916E-001, 5.201580198817631E-001,
     &5.465970120650943E-001, 5.725019326213812E-001,
     &5.978474702471789E-001, 6.226088602037082E-001,
     &6.467619085141294E-001, 6.702830156031415E-001,
     &6.931491993558019E-001, 7.153381175730568E-001,
     &7.368280898020211E-001, 7.575981185197069E-001,
     &7.776279096494952E-001, 7.968978923903148E-001,
     &8.153892383391762E-001, 8.330838798884006E-001,
     &8.499645278795909E-001, 8.660146884971646E-001,
     &8.812186793850185E-001, 8.955616449707267E-001,
     &9.090295709825298E-001, 9.216092981453337E-001,
     &9.332885350430794E-001, 9.440558701362557E-001,
     &9.539007829254919E-001, 9.628136542558153E-001,
     &9.707857757637066E-001, 9.778093584869182E-001,
     &9.838775407060567E-001, 9.889843952429919E-001,
     &9.931249370374432E-001, 9.962951347331254E-001,
     &9.984919506395956E-001, 9.997137267734412E-001/


      data xg3/
     &    1.4313661327935989E-004 ,   7.5402468020208113E-004 ,
     &    1.8524326334374286E-003 ,   3.4375314812782887E-003 ,
     &    5.5078023785041230E-003 ,   8.0612296469714795E-003 ,
     &    1.1095320756540850E-002 ,   1.4607112118146859E-002 ,
     &    1.8593172872092223E-002 ,   2.3049608537254129E-002 ,
     &    2.7972064931871987E-002 ,   3.3355732478460243E-002 ,
     &    3.9195350927333006E-002 ,   4.5485214508735161E-002 ,
     &    5.2219177514636506E-002 ,   5.9390660307490795E-002 ,
     &    6.6992655751417662E-002 ,   7.5017736060204343E-002 ,
     &    8.3458060055799588E-002 ,   9.2305380830411898E-002 ,
     &   0.10155105380484275      ,  0.11118604517525227      ,
     &   0.12120094074014642      ,  0.13158595509898963      ,
     &   0.14233094121347178      ,  0.15342540032209900      ,
     &   0.16485849219842952      ,  0.17661904574293535      ,
     &   0.18869556989814612      ,  0.20107626487641062      ,
     &   0.21374903368930942      ,  0.22670149396745293      ,
     &   0.23992099005911849      ,  0.25339460539590453      ,
     &   0.26710917511332100      ,  0.28105129891398428      ,
     &   0.29520735416084920      ,  0.30956350918768505      ,
     &   0.32410573681378912      ,  0.33881982804973543      ,
     &   0.35369140598076398      ,  0.36870593981424826      ,
     &   0.38384875907751304      ,  0.39910506795213202      ,
     &   0.41445995973069838      ,  0.42989843138194300      ,
     &   0.44540539820996944      ,  0.46096570859328168      ,
     &   0.47656415878920416      ,  0.49218550778922848      ,
     &   0.50781449221077157      ,  0.52343584121079578      ,
     &   0.53903429140671832      ,  0.55459460179003051      ,
     &   0.57010156861805694      ,  0.58554004026930162      ,
     &   0.60089493204786804      ,  0.61615124092248696      ,
     &   0.63129406018575174      ,  0.64630859401923602      ,
     &   0.66118017195026457      ,  0.67589426318621082      ,
     &   0.69043649081231495      ,  0.70479264583915080      ,
     &   0.71894870108601572      ,  0.73289082488667900      ,
     &   0.74660539460409547      ,  0.76007900994088151      ,
     &   0.77329850603254702      ,  0.78625096631069058      ,
     &   0.79892373512358938      ,  0.81130443010185394      ,
     &   0.82338095425706470      ,  0.83514150780157048      ,
     &   0.84657459967790105      ,  0.85766905878652822      ,
     &   0.86841404490101037      ,  0.87879905925985358      ,
     &   0.88881395482474779      ,  0.89844894619515725      ,
     &   0.90769461916958805      ,  0.91654193994420041      ,
     &   0.92498226393979566      ,  0.93300734424858234      ,
     &   0.94060933969250926      ,  0.94778082248536344      ,
     &   0.95451478549126478      ,  0.96080464907266694      ,
     &   0.96664426752153976      ,  0.97202793506812801      ,
     &   0.97695039146274587      ,  0.98140682712790772      ,
     &   0.98539288788185320      ,  0.98890467924345915      ,
     &   0.99193877035302847      ,  0.99449219762149588      ,
     &   0.99656246851872177      ,  0.99814756736656252      ,
     &   0.99924597531979797      ,  0.99985686338672064      /
      data wg3/
     &    3.6731724525279299E-004 ,   8.5469632675902094E-004 ,
     &    1.3419626857767407E-003 ,   1.8279806006590005E-003 ,
     &    2.3122250317102171E-003 ,   2.7942140019325745E-003 ,
     &    3.2734742254225979E-003 ,   3.7495366277323469E-003 ,
     &    4.2219357348344897E-003 ,   4.6902098268472130E-003 ,
     &    5.1539012874344866E-003 ,   5.6125570115929824E-003 ,
     &    6.0657288314897168E-003 ,   6.5129739464857726E-003 ,
     &    6.9538553518593933E-003 ,   7.3879422637206649E-003 ,
     &    7.8148105387730350E-003 ,   8.2340430880725977E-003 ,
     &    8.6452302841617968E-003 ,   9.0479703610640786E-003 ,
     &    9.4418698066874446E-003 ,   9.8265437472176628E-003 ,
     &    1.0201616323104732E-002 ,   1.0566721056263847E-002 ,
     &    1.0921501208123663E-002 ,   1.1265610128168154E-002 ,
     &    1.1598711592627035E-002 ,   1.1920480132984091E-002 ,
     &    1.2230601353978524E-002 ,   1.2528772240789786E-002 ,
     &    1.2814701455104074E-002 ,   1.3088109619772845E-002 ,
     &    1.3348729591785475E-002 ,   1.3596306723288451E-002 ,
     &    1.3830599110396229E-002 ,   1.4051377829550590E-002 ,
     &    1.4258427161197572E-002 ,   1.4451544800562611E-002 ,
     &    1.4630542055319155E-002 ,   1.4795244029956318E-002 ,
     &    1.4945489796666437E-002 ,   1.5081132552584564E-002 ,
     &    1.5202039763227438E-002 ,   1.5308093291990243E-002 ,
     &    1.5399189515576296E-002 ,   1.5475239425245528E-002 ,
     &    1.5536168713783277E-002 ,   1.5581917848104971E-002 ,
     &    1.5612442127424587E-002 ,   1.5627711726931670E-002 ,
     &    1.5627711726931670E-002 ,   1.5612442127424587E-002 ,
     &    1.5581917848104971E-002 ,   1.5536168713783277E-002 ,
     &    1.5475239425245528E-002 ,   1.5399189515576296E-002 ,
     &    1.5308093291990243E-002 ,   1.5202039763227438E-002 ,
     &    1.5081132552584564E-002 ,   1.4945489796666437E-002 ,
     &    1.4795244029956318E-002 ,   1.4630542055319155E-002 ,
     &    1.4451544800562611E-002 ,   1.4258427161197572E-002 ,
     &    1.4051377829550590E-002 ,   1.3830599110396229E-002 ,
     &    1.3596306723288451E-002 ,   1.3348729591785475E-002 ,
     &    1.3088109619772845E-002 ,   1.2814701455104074E-002 ,
     &    1.2528772240789786E-002 ,   1.2230601353978524E-002 ,
     &    1.1920480132984091E-002 ,   1.1598711592627035E-002 ,
     &    1.1265610128168154E-002 ,   1.0921501208123663E-002 ,
     &    1.0566721056263847E-002 ,   1.0201616323104732E-002 ,
     &    9.8265437472176628E-003 ,   9.4418698066874446E-003 ,
     &    9.0479703610640786E-003 ,   8.6452302841617968E-003 ,
     &    8.2340430880725977E-003 ,   7.8148105387730350E-003 ,
     &    7.3879422637206649E-003 ,   6.9538553518593933E-003 ,
     &    6.5129739464857726E-003 ,   6.0657288314897168E-003 ,
     &    5.6125570115929824E-003 ,   5.1539012874344866E-003 ,
     &    4.6902098268472130E-003 ,   4.2219357348344897E-003 ,
     &    3.7495366277323469E-003 ,   3.2734742254225979E-003 ,
     &    2.7942140019325745E-003 ,   2.3122250317102171E-003 ,
     &    1.8279806006590005E-003 ,   1.3419626857767407E-003 ,
     &    8.5469632675902094E-004 ,   3.6731724525279299E-004 /

      end

*************************************************************************
*
* $Id: besk0.F,v 1.1.1.1 1996/02/15 17:49:09 mclareni Exp $
*
* $Log: besk0.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:09  mclareni
* Kernlib
*
*
*#include "kernnum/pilot.h"
      REAL*8 FUNCTION BESK0(X)
      LOGICAL LEX
      parameter(LEX=.FALSE.)
      DOUBLE PRECISION X,Y,R,A,A0,A1,A2,B,B0,B1,B2,T(10)
      DOUBLE PRECISION U0,U1,U2,U3,U4,U5,U6,U7,U8,U9
      DOUBLE PRECISION F,F1,F2,F3,C,C0,PI1,CE,EPS,H,ALFA,D
      DOUBLE PRECISION ZERO,ONE,TWO,FOUR,FIVE,SIX,SEVEN,EIGHT,NINE,HALF
      DOUBLE PRECISION C1(0:14),C2(0:15),C3(0:12)
 
      DATA ZERO /0.0D0/, ONE /1.0D0/, TWO /2.0D0/
      DATA FOUR /4.0D0/, FIVE /5.0D0/, SIX /6.0D0/, SEVEN /7.0D0/
      DATA EIGHT /8.0D0/, NINE /9.0D0/, HALF /0.5D0/
 
      DATA T /16.0D0,368.0D0,43.0D0,75.0D0,400.0D0,40.0D0,
     1        48.0D0,12.0D0,20.0D0,28.0D0/
 
      DATA PI1 /1.25331 41373 155D0/, CE /0.57721 56649 0153D0/
      DATA EPS /1.0D-14/
 
      DATA C1( 0) /0.12773 34398 1218D3/
      DATA C1( 1) /0.19049 43201 7274D3/
      DATA C1( 2) /0.82489 03274 4024D2/
      DATA C1( 3) /0.22274 81924 2462D2/
      DATA C1( 4) /0.40116 73760 1793D1/
      DATA C1( 5) /0.50949 33654 3998D0/
      DATA C1( 6) /0.04771 87487 9817D0/
      DATA C1( 7) /0.00341 63317 6601D0/
      DATA C1( 8) /0.00019 24693 5969D0/
      DATA C1( 9) /0.00000 87383 1550D0/
      DATA C1(10) /0.00000 03260 9105D0/
      DATA C1(11) /0.00000 00101 6973D0/
      DATA C1(12) /0.00000 00002 6883D0/
      DATA C1(13) /0.00000 00000 0610D0/
      DATA C1(14) /0.00000 00000 0012D0/
 
      DATA C2( 0) /0.24027 70596 4072D3/
      DATA C2( 1) /0.36947 40739 7287D3/
      DATA C2( 2) /0.16997 34116 9840D3/
      DATA C2( 3) /0.49020 46377 7263D2/
      DATA C2( 4) /0.93884 97325 2684D1/
      DATA C2( 5) /0.12594 79763 6677D1/
      DATA C2( 6) /0.12377 69641 1492D0/
      DATA C2( 7) /0.00924 43098 6287D0/
      DATA C2( 8) /0.00054 06238 9649D0/
      DATA C2( 9) /0.00002 53737 9603D0/
      DATA C2(10) /0.00000 09754 7830D0/
      DATA C2(11) /0.00000 00312 4957D0/
      DATA C2(12) /0.00000 00008 4643D0/
      DATA C2(13) /0.00000 00000 1963D0/
      DATA C2(14) /0.00000 00000 0039D0/
      DATA C2(15) /0.00000 00000 0001D0/
 
      DATA C3( 0) /+0.98840 81742 3083D0/
      DATA C3( 1) /-0.01131 05046 4693D0/
      DATA C3( 2) /+0.00026 95326 1276D0/
      DATA C3( 3) /-0.00001 11066 8520D0/
      DATA C3( 4) /+0.00000 06325 7511D0/
      DATA C3( 5) /-0.00000 00450 4734D0/
      DATA C3( 6) /+0.00000 00037 9300D0/
      DATA C3( 7) /-0.00000 00003 6455D0/
      DATA C3( 8) /+0.00000 00000 3904D0/
      DATA C3( 9) /-0.00000 00000 0458D0/
      DATA C3(10) /+0.00000 00000 0058D0/
      DATA C3(11) /-0.00000 00000 0008D0/
      DATA C3(12) /+0.00000 00000 0001D0/

 
    9 IF(X .LE. ZERO) THEN
       print *,'C313.1 besk0',x
       BESK0=ZERO
       RETURN
      ENDIF

      IF(X .LT. HALF) THEN
       Y=X/EIGHT
       H=TWO*Y**2-ONE
       ALFA=-TWO*H
       B1=ZERO
       B2=ZERO
       DO 1 I = 14,0,-1
       B0=C1(I)-ALFA*B1-B2
       B2=B1
    1  B1=B0
       R=B0-H*B2
       B1=ZERO
       B2=ZERO
       DO 2 I = 15,0,-1
       B0=C2(I)-ALFA*B1-B2
       B2=B1
    2  B1=B0
       B1=-(CE+LOG(HALF*X))*R+B0-H*B2
       IF(LEX) B1=EXP(X)*B1
      ELSE IF(X .GT. FIVE) THEN
       R=ONE/X
       Y=FIVE*R
       H=TWO*Y-ONE
       ALFA=-TWO*H
       B1=ZERO
       B2=ZERO
       DO 3 I = 12,0,-1
       B0=C3(I)-ALFA*B1-B2
       B2=B1
    3  B1=B0
       B1=PI1*SQRT(R)*(B0-H*B2)
       IF(.NOT.LEX) B1=EXP(-X)*B1
      ELSE
       Y=(T(1)*X)**2
       A0=ONE
       A1=(T(1)*X+SEVEN)/NINE
       A2=(Y+T(2)*X+T(3))/T(4)
       B0=ONE
       B1=(T(1)*X+NINE)/NINE
       B2=(Y+T(5)*X+T(4))/T(4)
       U1=ONE
       U4=T(6)
       U5=T(7)
       C=ZERO
       F=TWO
    4  C0=C
       F=F+ONE
       U0=T(8)*F**2-ONE
       U1=U1+TWO
       U2=U1+TWO
       U3=U1+FOUR
       U4=U4+T(9)
       U5=U5+T(10)
       U6=ONE/U3**2
       U7=U2*U6
       U8=-U7/U1
       U9=T(1)*U7*X
       F1=U9-(U0-U4)*U8
       F2=U9-(U0-U5)*U6
       F3=-U8*(U3-SIX)**2
       A=F1*A2+F2*A1+F3*A0
       B=F1*B2+F2*B1+F3*B0
       C=A/B
       IF(ABS((C0-C)/C) .GE. EPS) THEN
        A0=A1
        A1=A2
        A2=A
        B0=B1
        B1=B2
        B2=B
        GO TO 4
       ENDIF
       B1=PI1*C/SQRT(X)
       IF(.NOT.LEX) B1=EXP(-X)*B1
      ENDIF

      BESK0=B1

c     IF(LEX)  THEN
c        IF(ENAME .EQ. 'EBESK0')  THEN
c           EBESK0=ROUND(B1)
c        ELSE
c           DEBSK0=B1
c        ENDIF
c     ELSE
c        IF(ENAME .EQ. ' BESK0')  THEN
c           BESK0=ROUND(B1)
c        ELSE
c           DBESK0=B1
c        ENDIF
c     ENDIF
      RETURN
 
  100 FORMAT(7X,A6,' ... NON-POSITIVE ARGUMENT X = ',E16.6)
      END

****************************************************************************
*
* $Id: besk1.F,v 1.1.1.1 1996/02/15 17:49:09 mclareni Exp $
*
* $Log: besk1.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:09  mclareni
* Kernlib
*
*
*#include "kernnum/pilot.h"
      REAL*8 FUNCTION BESK1(X)
      LOGICAL LEX
      parameter(LEX=.FALSE.)
      DOUBLE PRECISION X,Y,R,A,A0,A1,A2,B,B0,B1,B2,T(12)
      DOUBLE PRECISION U0,U1,U2,U3,U4,U5,U6,U7,U8,U9
      DOUBLE PRECISION F,F1,F2,F3,C,C0,PI1,CE,EPS,H,ALFA,D
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,EIGHT,HALF
      DOUBLE PRECISION C1(0:14),C2(0:14),C3(0:11)
 
      DATA ZERO /0.0D0/, ONE /1.0D0/, TWO /2.0D0/, THREE /3.0D0/
      DATA FOUR /4.0D0/, FIVE /5.0D0/, SIX /6.0D0/, EIGHT /8.0D0/
      DATA HALF /0.5D0/
 
      DATA T /16.0D0,3.2D0,2.2D0,432.0D0,131.0D0,35.0D0,336.0D0,
     1        40.0D0,48.0D0,12.0D0,20.0D0,28.0D0/
 
      DATA PI1 /1.25331 41373 155D0/, CE /0.57721 56649 0153D0/
      DATA EPS /1.0D-14/
 
      DATA C1( 0) /0.22060 14269 2352D3/
      DATA C1( 1) /0.12535 42668 3715D3/
      DATA C1( 2) /0.42865 23409 3128D2/
      DATA C1( 3) /0.94530 05229 4349D1/
      DATA C1( 4) /0.14296 57709 0762D1/
      DATA C1( 5) /0.15592 42954 7626D0/
      DATA C1( 6) /0.01276 80490 8173D0/
      DATA C1( 7) /0.00081 08879 0069D0/
      DATA C1( 8) /0.00004 10104 6194D0/
      DATA C1( 9) /0.00000 16880 4220D0/
      DATA C1(10) /0.00000 00575 8695D0/
      DATA C1(11) /0.00000 00016 5345D0/
      DATA C1(12) /0.00000 00000 4048D0/
      DATA C1(13) /0.00000 00000 0085D0/
      DATA C1(14) /0.00000 00000 0002D0/
 
      DATA C2( 0) /0.41888 94461 6640D3/
      DATA C2( 1) /0.24989 55490 4287D3/
      DATA C2( 2) /0.91180 31933 8742D2/
      DATA C2( 3) /0.21444 99505 3962D2/
      DATA C2( 4) /0.34384 15392 8805D1/
      DATA C2( 5) /0.39484 60929 4094D0/
      DATA C2( 6) /0.03382 87455 2688D0/
      DATA C2( 7) /0.00223 57203 3417D0/
      DATA C2( 8) /0.00011 71310 2246D0/
      DATA C2( 9) /0.00000 49754 2712D0/
      DATA C2(10) /0.00000 01746 0493D0/
      DATA C2(11) /0.00000 00051 4329D0/
      DATA C2(12) /0.00000 00001 2890D0/
      DATA C2(13) /0.00000 00000 0278D0/
      DATA C2(14) /0.00000 00000 0005D0/
 
      DATA C3( 0) /+1.03595 08587 724D0/
      DATA C3( 1) /+0.03546 52912 433D0/
      DATA C3( 2) /-0.00046 84750 282D0/
      DATA C3( 3) /+0.00001 61850 638D0/
      DATA C3( 4) /-0.00000 08451 720D0/
      DATA C3( 5) /+0.00000 00571 322D0/
      DATA C3( 6) /-0.00000 00046 456D0/
      DATA C3( 7) /+0.00000 00004 354D0/
      DATA C3( 8) /-0.00000 00000 458D0/
      DATA C3( 9) /+0.00000 00000 053D0/
      DATA C3(10) /-0.00000 00000 007D0/
      DATA C3(11) /+0.00000 00000 001D0/

 
    9 IF(X .LE. ZERO) THEN
       BESK1=ZERO
       RETURN
      ENDIF

      IF(X .LT. HALF) THEN
       Y=X/EIGHT
       H=TWO*Y**2-ONE
       ALFA=-TWO*H
       B1=ZERO
       B2=ZERO
       DO 1 I = 14,0,-1
       B0=C1(I)-ALFA*B1-B2
       B2=B1
    1  B1=B0
       R=Y*(B0-B2)
       B1=ZERO
       B2=ZERO
       DO 2 I = 14,0,-1
       B0=C2(I)-ALFA*B1-B2
       B2=B1
    2  B1=B0
       B1=(CE+LOG(HALF*X))*R+ONE/X-Y*(B0-B2)
       IF(LEX) B1=EXP(X)*B1
      ELSE IF(X .GT. FIVE) THEN
       R=ONE/X
       Y=FIVE*R
       H=TWO*Y-ONE
       ALFA=-TWO*H
       B1=ZERO
       B2=ZERO
       DO 3 I = 11,0,-1
       B0=C3(I)-ALFA*B1-B2
       B2=B1
    3  B1=B0
       B1=PI1*SQRT(R)*(B0-H*B2)
       IF(.NOT.LEX) B1=EXP(-X)*B1
      ELSE
       Y=(T(1)*X)**2
       A0=ONE
       A1=T(2)*X+T(3)
       A2=(Y+T(4)*X+T(5))/T(6)
       B0=ONE
       B1=T(2)*X+ONE
       B2=(Y+T(7)*X+T(6))/T(6)
       U1=ONE
       U4=T(8)
       U5=T(9)
       C=ZERO
       F=TWO
    4  C0=C
       F=F+ONE
       U0=T(10)*F**2+THREE
       U1=U1+TWO
       U2=U1+TWO
       U3=U1+FOUR
       U4=U4+T(11)
       U5=U5+T(12)
       U6=ONE/(U3**2-FOUR)
       U7=U2*U6
       U8=-U7/U1
       U9=T(1)*U7*X
       F1=U9-(U0-U4)*U8
       F2=U9-(U0-U5)*U6
       F3=U8*(FOUR-(U3-SIX)**2)
       A=F1*A2+F2*A1+F3*A0
       B=F1*B2+F2*B1+F3*B0
       C=A/B
       IF(ABS((C0-C)/C) .GE. EPS) THEN
        A0=A1
        A1=A2
        A2=A
        B0=B1
        B1=B2
        B2=B
        GO TO 4
       ENDIF
       B1=PI1*C/SQRT(X)
       IF(.NOT.LEX) B1=EXP(-X)*B1
      ENDIF

      BESK1=B1

c     IF(LEX)  THEN
c        IF(ENAME .EQ. 'EBESK1')  THEN
c           EBESK1=ROUND(B1)
c        ELSE
c           DEBSK1=B1
c        ENDIF
c     ELSE
c        IF(ENAME .EQ. ' BESK1')  THEN
c           BESK1=ROUND(B1)
c        ELSE
c           DBESK1=B1
c        ENDIF
c     ENDIF
      RETURN
 
  100 FORMAT(7X,A6,' ... NON-POSITIVE ARGUMENT X = ',E16.6)
      END


c...momentum dependent potential part
c**************************************************************************
c divided by 4pi
      real*8 function potMomdepT0(p,pf,lam)
c...single particle potential at T=0
      implicit none
      real*8 p,pf,lam,x,y

      x=p/lam
      y=pf/lam
      potMomdepT0=1.0/4.0*lam**3*( 
     & 2*y + 0.5/x*(y**2 - x**2 +1.0)*log(((x+y)**2+1.0)/((x-y)**2+1.0))
     & - 2*(atan(x+y) - atan(x-y)) )

      end

c**************************************************************************
      subroutine spot_momt0(bmu,spot)
c...Compute momentum dependent single particle potential at T=0.
      implicit none
      include 'jam2.inc'
      integer i,NGP
      real*8 bmu,g,fac,pf,potmomdept0,c1,c2
      parameter(NGP=38)
      real*8 xg1(NGP),wg1(NGP)
      common/gaussgrid1/xg1,wg1
      real*8 spot(NGP)

      c1=pard(107)/parc(21)
      c2=pard(108)/parc(21)

      g=4.0
      fac=g/(2*paru(1)**2*paru(3)**3)
      pf=sqrt(max(0.0,bmu**2 - 0.938**2))
c     pf=(6*paru(1)*2*rhob/g)**(1.0/3.0)*paru(3)
      do i=1,NGP
        spot(i)=fac*(c1*potMomdepT0(xg1(i),pf,pard(105))
     &              +c2*potMomdepT0(xg1(i),pf,pard(106)))
      end do

      end

c**************************************************************************
      real*8 function spotmom2(pin,tf,mu,pm,ibar)
c...Compute momentum dependent single particle potential for given pin.
      implicit none
      include 'jam2.inc'
      real*8 tf,mu,mu1,mu2,c1,c2,pm,pin,d,dis,pot1,pot2,p,e
      integer j,istat,ibar
      integer NGP
      parameter(NGP=38)
      real*8 xg1(NGP),wg1(NGP)
      common/gaussgrid1/xg1,wg1
      real*8 spot(NGP)
      common/mdepot/spot

      mu1=pard(105)**2
      mu2=pard(106)**2
      c1=pard(107)/parc(21)
      c2=pard(108)/parc(21)

      istat=1
c...Numerical integration.
      d=0.0d0
      do j=1, NGP
        p=xg1(j)
        e=( sqrt(p**2+pm**2) + ibar*spot(j) - mu )/tf
        dis=exp(-e/tf)
        dis = dis/(1.0+istat*dis)
c       dis = 1.0/(exp(min(50.0d0,e))+istat)
        pot1=c1*mu1*log( (mu1+(pin+p)**2)/(mu1+(pin-p)**2) )
        pot2=c2*mu2*log( (mu2+(pin+p)**2)/(mu2+(pin-p)**2) )
        d = d + p*wg1(j)*dis*(pot1+pot2)
      end do
      spotmom2=d/(8.0*paru(1)**2*paru(3)**3)/pin

      end

c**************************************************************************
      real*8 function spotmom(pin)
c...Compute momentum dependent single particle potential for given pin.
      implicit none
      include 'jam2.inc'
      integer i,kc,kf,ispin,ns,ibary,isp
      real*8 spotmom2,cmu,smu,bmu,pm,tch,pin,slope
      integer kcmin,kcmax
      parameter(kcmin=132, kcmax=366)
      integer NGP
      parameter(NGP=38)
      real*8 xg1(NGP),wg1(NGP)
      common/gaussgrid1/xg1,wg1
      real*8 spot(NGP)
      common/mdepot/spot

      spotmom=0.0
      if(mstd(101).eq.0) return

c...Linear interpolation.
      if(mstc(135).eq.1) then
        if(pin.lt.xg1(1)) then
           slope=(spot(2)-spot(1))/(xg1(2)-xg1(1))
           spotmom=spot(1)+(pin-xg1(1))*slope
           return
        endif
        do i=1,NGP-1
          if(pin.ge.xg1(i).and.pin.le.xg1(i+1)) then
             slope=(spot(i+1)-spot(i))/(xg1(i+1)-xg1(i))
             spotmom=spot(i)+(pin-xg1(i))*slope
             return
          endif
        end do
        slope=(spot(NGP)-spot(NGP-1))/(xg1(NGP)-xg1(NGP-1))
        spotmom=spot(NGP-1)+(pin-xg1(NGP-1))*slope
        return
      endif

      tch=pare(51)
      bmu=pare(52)
      smu=pare(53)
c...Loop over particles.
      do kc=kcmin,kcmax

c....Exclude heavy flavour hadrons.
        if(kc.ge.195.and.kc.le.251) cycle
        ibary=kchg(kc,6)/3*isign(1,kf)   ! baryon number times 3.
        if(ibary.eq.0) cycle             ! exclude mesons

        kf=kchg(kc,4)
        ispin=max(1,mod(kf,10)) ! spin 2J+1
        ns=kchg(kc,7)*isign(1,kf)
        cmu = ns*smu + ibary*bmu
        pm=pmas(kc,1)                 ! particle mass

c....Compute potential.
        spotmom = spotmom + ispin*spotmom2(pin,tch,cmu,pm,ibary)

c......Anti-particle if it exists.
        if(kchg(kc,3).eq.1) then
          spotmom = spotmom + ispin*spotmom2(pin,tch,-cmu,pm,ibary)
        endif

      end do

      end

c**************************************************************************
      subroutine spotentialm(tf,mu,pm,ibar,potm)
c...Compute momentum dependent single particle potential.
      implicit none
      include 'jam2.inc'
      real*8 tf,mu,mu1,mu2,c1,c2,pm,pin,d,dis,pot1,pot2,p,dw,e
      integer i,j,istat,ibar
      integer NGP
      parameter(NGP=38)
      real*8 potm(NGP)
      real*8 xg1(NGP),wg1(NGP)
      common/gaussgrid1/xg1,wg1
      real*8 spot(NGP)
      common/mdepot/spot

      mu1=pard(105)**2
      mu2=pard(106)**2
      c1=pard(107)/parc(21)
      c2=pard(108)/parc(21)

      istat=1
c...Numerical integration.
      do i=1, NGP
      d=0.0d0
      pin=xg1(i)
      do j=1, NGP
        p=xg1(j)
        dw=wg1(j)
        e=( sqrt(p**2+pm**2) + ibar*spot(j) - mu )/tf
        dis=exp(-e/tf)
        dis = dis/(1.0+istat*dis)
c       dis = 1.0/(exp(min(50.0d0,e))+istat)
        pot1=c1*mu1*log( (mu1+(pin+p)**2)/(mu1+(pin-p)**2) )
        pot2=c2*mu2*log( (mu2+(pin+p)**2)/(mu2+(pin-p)**2) )
        d = d + p*dw*dis*(pot1+pot2)
      end do
      potm(i)=d/(8.0*paru(1)**2*paru(3)**3)/pin
      end do


      end

c**************************************************************************
      subroutine spot_mom(tch,bmu,smu)
c...Compute momentum dependent single particle potential.
      implicit none
      include 'jam2.inc'
      integer i,kc,kf,ispin,ns,ibary,isp
      real*8 pmult,cmu,smu,bmu,pm,dns,pdensity,tch,fac
      integer NGP
      parameter(NGP=38)
      real*8 potm(NGP),potm2(NGP)
      integer kcmin,kcmax
      parameter(kcmin=132, kcmax=366)
      real*8 spot(NGP)
      common/mdepot/spot

      potm=0.0

c...Loop over particles.
      do kc=kcmin,kcmax

c....Exclude heavy flavour hadrons.
        if(kc.ge.195.and.kc.le.251) cycle
        ibary=kchg(kc,6)/3               ! baryon number times 3.
        if(ibary.eq.0) cycle             ! exclude mesons

        kf=kchg(kc,4)
        ispin=max(1,mod(kf,10)) ! spin 2J+1
        ns=kchg(kc,7)*isign(1,kf)
        cmu = ns*smu + ibary*bmu
        pm=pmas(kc,1)                 ! particle mass

c....Compute potential.
        call spotentialm(tch,cmu,pm,ibary,potm2)
        potm = potm + potm2*ispin

c......Anti-particle if it exists.
        if(kchg(kc,3).eq.1) then
          call spotentialm(tch,-cmu,pm,ibary,potm2)
          potm = potm + potm2*ispin
        endif

      end do

      spot=potm

      end

c*************************************************************************
      subroutine momdeppot(tch,bmu,smu)
c...Compute momentum dependent single particle potential.
      implicit none
      include 'jam2.inc'
      real*8 tch,bmu,smu,err
      integer i,NGP
      parameter(NGP=38)
      real*8 spot2(NGP)
      real*8 spot(NGP)
      common/mdepot/spot

c...get initial condition for potential.
      call spot_momt0(bmu,spot)

      do i=1,100
      spot2=spot
c     print *,'spot2=',spot2
      call spot_mom(tch,bmu,smu)
      err = sum(spot-spot2)
c      print *,i,err
c      print *,'spot=',spot
c      print *,'error=',err
      if(abs(err).lt.1d-5) exit
      end do

      if(abs(err).gt.1d-5) then
        print *,'warning momdeppot not converge error=',err
      endif

      end 


