c======================================================================
c...last modified Y.Nara Sep.8  2017
c...last modified Y.Nara Oct.8  2017
c...last modified Y.Nara Oct.30 2017
c...last modified Y.Nara Nov.27 2017 version 0.005
c...last modified Y.Nara Dec. 7 2017 version 0.008a
c======================================================================
c  U(0) = g^2*(e+p)-p
c  U(1) = g^2*(e+p)*v_x
c  U(2) = g^2*(e+p)*v_y
c  U(3) = g^2*(e+p)*v_z 
c  U(4) = g*n_B
c*********************************************************************
      block data hydrodata
      implicit none
      include 'fluid.inc'

      data h_mode/1/  ! PPM by Hirano
      data op_splitting/.true./
      data second_order/0/
      data dx,dy,dz,dt/0.3d0,0.3d0,0.3d0,0.3d0/
      data tau0/0.6d0/
      data gsigma/1.0d0/
      data opt_taueta/0/
      data opt_timelike,opt_rescale/0,1/
      data opt_freezeout/0/
      data job_mode/0/
      data efreeze/2.02709208572572430534d0/ ! 0.4 GeV/fm^3
      data tfreezecut/0.015203190642942932290/   ! 0.003 GeV

      end

c*********************************************************************
      subroutine fluid_initialtime(t)
      implicit none
      include 'fluid.inc'
      real*8 t
      tau0=t
      end

c*********************************************************************
      subroutine fluid_evolution1(it,check)
      implicit none
      include 'fluid.inc'
      integer ntime,it,check
      double precision htime,pxdir(2)
      double precision etot,btot,pxtot,pytot,pztot
      real*8 emax,tmax,pmax,bmax,t,aveb,r
      common/fluid1/emax,tmax,pmax,bmax,aveb

c     htime = it* dt + tau0
      htime = it* dt

c...3D evolution by the Operator splitting.
c       if(op_splitting) then
          call  fluid_evolution3dos(it,htime,dt,opt_taueta)
c       else
c         call  fluid_evolution3d(h_mode)
c       endif

        if(mod(opt_freezeout,10).eq.2) then
          call c_freezeout(2,it,htime)
        else if(opt_freezeout.eq.11) then
          call isochronous_freezeout(htime,0)
        endif

        call checkout(htime,check)

        if(job_mode.ge.1) then
        if(mod(it,5).eq.0) then
c         call directedFlow(1,it,pxdir)
c         call directedFlow(2,it,pxdir)
c         call directedFlow(3,it,pxdir)
c         write(30,820)htime,aveb,bmax,emax,pxdir(1),pxdir(2)
          call outputhzx(it)
        endif
        endif

c         write(30,820)htime+dt,bl(maxx/2,maxy/2,maxz/2),
c    &          el(maxx/2,maxy/2,maxz/2)*hbc
c    &          ,u(4,maxx/2,maxy/2,maxz/2)
c    &          ,u(0,maxx/2,maxy/2,maxz/2)*hbc
c    &          ,vx(maxx/2,maxy/2,maxz/2)
c    &          ,vy(maxx/2,maxy/2,maxz/2)
c    &          ,vz(maxx/2,maxy/2,maxz/2)
c         flush(30)

820   format(f6.3,8(1x,f10.4))
      end

c************************************************************
      subroutine reset_fluid_val
      implicit none
      include 'fluid.inc'

      u=0; tl=0; el=0; bl=0; pl=0; ml=0
      cs=0; vx=0; vy=0; vz=0
      uzero=0; uout=0
      ifreeze=0

      px_tot=0.0d0
      py_tot=0.0d0
      pz_tot=0.0d0
      pe_tot=0.0d0
      b_tot=0.0d0
      c_tot=0.0d0
      s_tot=0.0d0

      end

c************************************************************
      subroutine reset_fluid_element(x,y,z)
      implicit none
      include 'fluid.inc'
      integer x,y,z

      u(:,x,y,z)=0; tl(x,y,z)=0; el(x,y,z)=0; bl(x,y,z)=0;
      pl(x,y,z)=0; ml(x,y,z)=0
      cs(x,y,z)=0; vx(x,y,z)=0; vy(x,y,z)=0; vz(x,y,z)=0
c     uzero=0; uout=0

      end

c*********************************************************************
      subroutine fluid_evolution(ntime)
      implicit none
      include 'fluid.inc'
      integer ntime,it,check
      double precision htime,pxdir(2)
      double precision etto,pxtot,pytot,pztot,btot
      real*8 emax,tmax,pmax,bmax,t,aveb,r
      common/fluid1/emax,tmax,pmax,bmax,aveb

      htime=tau0

      if(opt_freezeout.eq.1) then
        write(*,*) 'before c_freezeout, htime=',htime
        call c_freezeout(1,0,htime)              ! for freezeout
      endif

      call outputhzx(0)
c     call outputhxy(0)

c...Loop over time step.
      do it=0,ntime

        call hydroTotalEnergy(it,htime,etto,pxtot,pytot,pztot,btot)

c...3D evolution by the Operator splitting.
c       if(op_splitting) then
          if(second_order.le.1) then
            call  fluid_evolution3dos(it,htime,dt,opt_taueta)
          else 
            call  fluid_evolution3d2(it,htime)
          endif
c       else
c         call  fluid_evolution3d(h_mode)
c       endif

        htime = htime + dt
c       call hydroTotalEnergy(it,htime,etto,pxtot,pytot,pztot,btot)
        if(opt_freezeout.eq.1) call c_freezeout(2,it,htime)
        call checkout(htime,check)
        if(mod(it,5).eq.0) then
          call directedFlow(1,it,pxdir)
          call directedFlow(2,it,pxdir)
          call directedFlow(3,it,pxdir)
          call outputhzx(it)
c         call outputhxy(it)
          write(30,820)htime,aveb,bmax,emax,pxdir(1),pxdir(2)
          flush(30)
        endif
        if(check.eq.0) exit
      end do

      if(mod(opt_freezeout,10).eq.2) call isochronous_freezeout(htime,0)

820   format(f6.3,5(1x,f10.4))
      end

c************************************************************
c....hydrodynamical evolution in full 3-dimension by the operator
c...splitting method.
      subroutine fluid_evolution3d2(it,tau)

      implicit none
      include 'fluid.inc'
      integer it,is,iz
      real*8 tau,dt0

      is = mod(it,6)
      iz=0
      if(opt_taueta.ge.1) iz=1

      dt0=dt
      if(is.eq.0) then
        dt=0.5*dt0;                 call update_z(tau,0)
        dt=0.5*dt0; tau=tau+dt/3d0; call update_x(tau,0)
        dt=dt0;   ; tau=tau+dt/3d0; call update_y(tau,0)
        dt=0.5*dt0; tau=tau+dt/3d0; call update_x(tau,0)
        dt=0.5*dt0; tau=tau+dt/3d0; call update_z(tau,iz)
      else if(is.eq.1) then
        dt=0.5*dt0;               ; call update_x(tau,0)
        dt=0.5*dt0; tau=tau+dt/3d0; call update_z(tau,0)
        dt=    dt0; tau=tau+dt/3d0; call update_y(tau,0)
        dt=0.5*dt0; tau=tau+dt/3d0; call update_z(tau,0)
        dt=0.5*dt0; tau=tau+dt/3d0; call update_x(tau,iz)
      else if(is.eq.2) then
        dt=0.5*dt0;                   call update_y(tau,0)
        dt=0.5*dt0; tau=tau+dt/3.0d0; call update_x(tau,0)
        dt=dt0;     tau=tau+dt/3.0d0; call update_z(tau,0)
        dt=0.5*dt0; tau=tau+dt/3.0d0; call update_x(tau,0)
        dt=0.5*dt0; tau=tau+dt/3.0d0; call update_y(tau,iz)
      else if(is.eq.3) then
        dt=0.5*dt0;                  call update_z(tau,0)
        dt=0.5*dt0; tau=tau+dt/3.0d0;call update_y(tau,0)
        dt=dt0;     tau=tau+dt/3.0d0;call update_x(tau,0)
        dt=0.5*dt0; tau=tau+dt/3.0d0;call update_y(tau,0)
        dt=0.5*dt0; tau=tau+dt/3.0d0;call update_z(tau,iz)
      else if(is.eq.4) then
        dt=0.5*dt0;                  call update_x(tau,0)
        dt=0.5*dt0; tau=tau+dt/3.0d0;call update_y(tau,0)
        dt=dt0;     tau=tau+dt/3.0d0;call update_z(tau,0)
        dt=0.5*dt0; tau=tau+dt/3.0d0;call update_y(tau,0)
        dt=0.5*dt0; tau=tau+dt/3.0d0;call update_x(tau,iz)
      else if(is.eq.5) then
        dt=0.5*dt0;                   call update_y(tau,0)
        dt=0.5*dt0; tau=tau+dt/3.0d0; call update_z(tau,0)
        dt=dt0;     tau=tau+dt/3.0d0; call update_x(tau,0)
        dt=0.5*dt0; tau=tau+dt/3.0d0; call update_z(tau,0)
        dt=0.5*dt0; tau=tau+dt/3.0d0; call update_y(tau,iz)
      endif
      dt=dt0

      end

c************************************************************
c....hydrodynamical evolution in full 3-dimension by the operator
c...splitting method.
      subroutine fluid_evolution3dos(it,tau,dt,opt_taueta)

      implicit none
      integer it,is,opt_taueta,s1,s2,s3
      real*8 tau,dt,tau1,tau2,tau3

      s1=0; s2=0; s3=0
      tau1=1.0d0
      tau2=1.0d0
      tau3=1.0d0
      if(opt_taueta.ge.1) then
c       s1=3; s2=3; s3=1  ! Include source term for all direction
        s1=0; s2=0; s3=1  ! Include source term only for z-direction
        tau1=tau
        tau2=tau+dt/3d0
        tau3=tau+2d0*dt/3d0
      endif

      is = mod(it,6)

      if(is.eq.0) then
        call update_z(tau1,s3)
        call update_x(tau2,s1)
        call update_y(tau3,s2)
      else if(is.eq.1) then
        call update_x(tau1,s1)
        call update_z(tau2,s3)
        call update_y(tau3,s2)
      else if(is.eq.2) then
        call update_y(tau1,s2)
        call update_x(tau2,s1)
        call update_z(tau3,s3)
      else if(is.eq.3) then
        call update_z(tau1,s3)
        call update_y(tau2,s2)
        call update_x(tau3,s1)
      else if(is.eq.4) then
        call update_x(tau1,s1)
        call update_y(tau2,s2)
        call update_z(tau3,s3)
      else if(is.eq.5) then
        call update_y(tau1,s2)
        call update_z(tau2,s3)
        call update_x(tau3,s1)
      endif

      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fluid_evolution3dos2(it0,tau,dt,opt_taueta)

c....Hydrodynamical evolution in full 3-dimension by the operator
c...splitting method.

      implicit none
      integer it0,i,it,is,iz,opt_taueta
      double precision tau,dt

      it=3*it0

c...Operator splitting.
      do i=1,3

      is = mod(it,6)
      iz=0

c...x-direction.
      if((is .eq. 0) .or. (is .eq. 4)) then
        if(opt_taueta.ge.1) iz=1
        call update_z(tau,iz)
        tau=tau+dt/3d0
c...y-direction.
      else if((is .eq. 1) .or. (is .eq. 3)) then
        call update_x(tau,iz)
        tau=tau+dt/3d0
c...z-direction.
      else if((is .eq. 2) .or. (is .eq. 5)) then
        call update_y(tau,iz)
        tau=tau+dt/3d0
      endif

      it = it +1
      end do

      end

c******************************************************************
      subroutine update_x(tau,is)
      implicit none
      real*8 tau,lam,lhm,tau1,tau2,tau3
      include 'fluid.inc'
      integer iy,iz,n,is
      real*8 uu(0:4,0:maxx),uh(0:4,0:maxx),g(0:4,0:maxx)
      real*8 uv(0:4),uo(0:4)
      real*8 vv(3,0:maxx),sv(0:maxx)
      real*8 er(0:maxx),pr(0:maxx),nr(0:maxx)
      real*8 tr(0:maxx),mr(0:maxx)

      tau1=1.0d0
      tau2=1.0d0
      tau3=1.0d0
      if(opt_taueta.ge.1) then
        tau1=tau
        tau2=tau+dt/6.0d0
        tau3=tau+dt/3.0d0
      endif
      lam=dt/dx
      lhm=0.5*lam
      n=2
      do iy = n,maxy-n
      do iz = n,maxz-n
        uu(1,:)=u(1,:,iy,iz)
        uu(2,:)=u(2,:,iy,iz)
        uu(3,:)=u(3,:,iy,iz)
        uu(0,:)=u(0,:,iy,iz)
        uu(4,:)=u(4,:,iy,iz)
        vv(1,:)=vx(:,iy,iz)
        vv(2,:)=vy(:,iy,iz)
        vv(3,:)=vz(:,iy,iz)
        sv=cs(:,iy,iz)
        pr=pl(:,iy,iz)

        if(h_mode.eq.1) then
          uh=uu
          if(second_order.eq.1) then
            call ppm1d(tau1,maxx,dx,0.5*dt,uh,vv(1,:),sv,pr,
     &            g,opt_timelike,0)
            uh=uh+g
            call update_local_val(tau2,maxx,uh,vv,sv,er,pr,nr,tr,mr,uv,
     &                             uo,opt_rescale)
          endif
          call ppm1d(tau1,maxx,dx,dt,uh,vv(1,:),sv,pr,g,opt_timelike,0)
          call add_source(tau1,dt,dx,maxx,uu,g,pr,vv(1,:),vv(3,:),
     &    is,opt_taueta,opt_timelike)

        else if(h_mode.eq.2) then
          nr=bl(:,iy,iz); er=el(:,iy,iz)
          if(second_order.eq.1) then
            call RHS(maxx,nr,vv(1,:),vv(2,:),vv(3,:),er,sv,lhm,g,
     &             opt_timelike)
            call update_local_val(1d0,maxx,uu,vv,sv,er,pr,nr,tr,mr,uv,
     &                             uo,opt_rescale)
          endif

          call RHS(maxx,nr,vv(1,:),vv(2,:),vv(3,:),er,sv,lam,g,
     &             opt_timelike)
          g=g*tau1
          call add_source(tau1,dt,dx,maxx,uu,g,pr,vv(1,:),vv(3,:),
     &    is,opt_taueta,opt_timelike)


        else
          print *,'update_x wrong h_mode',h_mode
          stop
        endif


        call update_local_val(tau3,maxx,uu,vv,sv,er,pr,nr,tr,mr,uv,
     &                             uo,opt_rescale)
      
        uzero = uzero+uv
        uout  = uout + uo

        u(1,:,iy,iz) = uu(1,:)
        u(2,:,iy,iz) = uu(2,:)
        u(3,:,iy,iz) = uu(3,:)
        u(0,:,iy,iz) = uu(0,:)
        u(4,:,iy,iz) = uu(4,:)
        bl(:,iy,iz) = nr
        el(:,iy,iz) = er
        pl(:,iy,iz) = pr
        tl(:,iy,iz) = tr
        ml(:,iy,iz) = mr
        vx(:,iy,iz) = vv(1,:)
        vy(:,iy,iz) = vv(2,:)
        vz(:,iy,iz) = vv(3,:)
        cs(:,iy,iz) = sv
      end do
      end do

      call checkboundary_x

      end

c******************************************************************
      subroutine update_y(tau,is)
      implicit none
      real*8 tau,lam,lhm,tau1,tau2,tau3
      include 'fluid.inc'
      integer ix,iz,n,is
      real*8 uu(0:4,0:maxy),uh(0:4,0:maxy),g(0:4,0:maxy)
      real*8 uv(0:4),uo(0:4)
      real*8 vv(3,0:maxy),sv(0:maxy)
      real*8 er(0:maxy),pr(0:maxy),nr(0:maxy)
      real*8 tr(0:maxy),mr(0:maxy)

      tau1=1.0d0
      tau2=1.0d0
      tau3=1.0d0
      if(opt_taueta.ge.1) then
        tau1=tau
        tau2=tau+0.5*dt/3.0d0
        tau3=tau+dt/3.0d0
      endif

      lam=dt/dy
      lhm=0.5*lam
      n=2
      do ix = n,maxx-n
      do iz = n,maxz-n
        uu(1,:)=u(2,ix,:,iz)
        uu(2,:)=u(1,ix,:,iz)
        uu(3,:)=u(3,ix,:,iz)
        uu(0,:)=u(0,ix,:,iz)
        uu(4,:)=u(4,ix,:,iz)
        vv(1,:)=vy(ix,:,iz)
        vv(2,:)=vx(ix,:,iz)
        vv(3,:)=vz(ix,:,iz)
        sv=cs(ix,:,iz)
        pr=pl(ix,:,iz)
        if(h_mode.eq.1) then
          uh=uu
          if(second_order.eq.1) then
            call ppm1d(tau1,maxy,dy,0.5*dt,uh,vv(1,:),sv,pr,
     &            g,opt_timelike,0)
            uh=uh+g
            call update_local_val(tau2,maxy,uh,vv,sv,er,pr,nr,tr,mr,uv,
     &                             uo,opt_rescale)
          endif
          call ppm1d(tau1,maxy,dy,dt,uu,vv(1,:),sv,pr,g,opt_timelike,0)
          call add_source(tau1,dt,dy,maxy,uu,g,pr,vv(1,:),vv(3,:),
     &    is,opt_taueta,opt_timelike)

        else if(h_mode.eq.2) then
          nr=bl(ix,:,iz); er=el(ix,:,iz)
          if(second_order.eq.1) then
            call RHS(maxy,nr,vv(1,:),vv(2,:),vv(3,:),er,sv,lhm,g,
     &             opt_timelike)
            uh = uh + g
            call update_local_val(1d0,maxy,uu,vv,sv,er,pr,nr,tr,mr,uv,
     &                             uo,opt_rescale)
          endif
          call RHS(maxy,nr,vv(1,:),vv(2,:),vv(3,:),er,sv,lam,g,
     &             opt_timelike)
          g=g*tau1
          call add_source(tau1,dt,dy,maxy,uu,g,pr,vv(1,:),vv(3,:),
     &    is,opt_taueta,opt_timelike)

        else
          print *,'update_y wrong h_mode',h_mode
          stop
        endif


        call update_local_val(tau3,maxy,uu,vv,sv,er,pr,nr,tr,mr,uv,
     &                             uo,opt_rescale)
        uzero= uzero+uv
        uout  = uout + uo
        u(2,ix,:,iz) = uu(1,:)
        u(1,ix,:,iz) = uu(2,:)
        u(3,ix,:,iz) = uu(3,:)
        u(0,ix,:,iz) = uu(0,:)
        u(4,ix,:,iz) = uu(4,:)
        bl(ix,:,iz) = nr
        el(ix,:,iz) = er
        pl(ix,:,iz) = pr
        tl(ix,:,iz) = tr
        ml(ix,:,iz) = mr
        vy(ix,:,iz) = vv(1,:)
        vx(ix,:,iz) = vv(2,:)
        vz(ix,:,iz) = vv(3,:)
        cs(ix,:,iz) = sv
      end do
      end do

      call checkboundary_y

      end

c******************************************************************
      subroutine update_z(tau,is)
      implicit none
      real*8 tau,lam,lhm,tau1,tau2,tau3
      include 'fluid.inc'
      integer ix,iy,n,is
      real*8 uu(0:4,0:maxz),uh(0:4,0:maxz),g(0:4,0:maxz)
      real*8 uv(0:4),uo(0:4)
      real*8 vv(3,0:maxz),sv(0:maxz)
      real*8 er(0:maxz),pr(0:maxz),nr(0:maxz)
      real*8 tr(0:maxz),mr(0:maxz)

      tau1=1.0d0
      tau2=1.0d0
      tau3=1.0d0
      if(opt_taueta.ge.1) then
        tau1=tau
        tau2=tau+0.5*dt/3.0d0
        tau3=tau+dt/3.0d0
      endif

      lam=dt/(dz*tau1)
      lhm=0.5*lam
      n=2

c     is=0
c     if(opt_taueta.eq.1) is=1

      do ix = n,maxx-n
      do iy = n,maxy-n
        uu(1,:)=u(3,ix,iy,:)
        uu(2,:)=u(2,ix,iy,:)
        uu(3,:)=u(1,ix,iy,:)
        uu(0,:)=u(0,ix,iy,:)
        uu(4,:)=u(4,ix,iy,:)
        vv(1,:)=vz(ix,iy,:)
        vv(2,:)=vy(ix,iy,:)
        vv(3,:)=vx(ix,iy,:)
        sv=cs(ix,iy,:)
        pr=pl(ix,iy,:)
        if(h_mode.eq.1) then
          uh=uu
          if(second_order.eq.1) then
            call ppm1d(tau1,maxz,dz*tau1,0.5*dt,uh,vv(1,:),sv,pr,g,
     &            opt_timelike,0)
            uh=uh+g
            call update_local_val(tau2,maxz,uh,vv,sv,er,pr,nr,tr,mr,uv,
     &                             uo,opt_rescale)
          endif
          call ppm1d(tau1,maxz,dz*tau1,dt,uh,vv(1,:),sv,pr,g,
     &       opt_timelike,0)
          call add_source(tau1,dt,dz,maxz,uu,g,pr,vv(1,:),vv(1,:),
     &    is,opt_taueta,opt_timelike)

        else if(h_mode.eq.2) then
          nr=bl(ix,iy,:); er=el(ix,iy,:)
          uh=uu
          if(second_order.eq.1) then
            call RHS(maxz,nr,vv(1,:),vv(2,:),vv(3,:),er,sv,lhm,g,
     &        opt_timelike)
            uh = uh + g
            call update_local_val(1d0,maxz,uh,vv,sv,er,pr,nr,tr,mr,uv,
     &                             uo,opt_rescale)
          endif
          call RHS(maxz,nr,vv(1,:),vv(2,:),vv(3,:),er,sv,lam,g,
     &        opt_timelike)
          g=g*tau1
          call add_source(tau1,dt,dz,maxz,uu,g,pr,vv(1,:),vv(1,:),
     &    is,opt_taueta,opt_timelike)

        else
          print *,'update_z wrong h_mode',h_mode
          stop
        endif


        call update_local_val(tau3,maxz,uu,vv,sv,er,pr,nr,tr,mr,uv,
     &                             uo,opt_rescale)
        uzero = uzero + uv
        uout  = uout  + uo
        u(3,ix,iy,:) = uu(1,:)
        u(2,ix,iy,:) = uu(2,:)
        u(1,ix,iy,:) = uu(3,:)
        u(0,ix,iy,:) = uu(0,:)
        u(4,ix,iy,:) = uu(4,:)
        bl(ix,iy,:) = nr
        el(ix,iy,:) = er
        pl(ix,iy,:) = pr
        tl(ix,iy,:) = tr
        ml(ix,iy,:) = mr
        vz(ix,iy,:) = vv(1,:)
        vy(ix,iy,:) = vv(2,:)
        vx(ix,iy,:) = vv(3,:)
        cs(ix,iy,:) = sv
      end do
      end do

      call checkboundary_z

      end

c*************************************************************************
      subroutine rescaleu(maxx,uu)

      integer i
      real*8 uu(0:4,0:maxx),u2
      real*8 eps
      parameter(eps=1d-12)

      do i = 0,maxx
         u2=sqrt(uu(1,i)**2+uu(2,i)**2+uu(3,i)**2)
         if(uu(0,i).le.eps) then
           uu(:,i)=0d0
         else if (uu(0,i).lt.u2) then
            if (u2.gt.1d-16) then
c              uu(0,i) = u2+eps
               uu(1,i) = uu(0,i)*(1d0-eps)*uu(1,i)/u2
               uu(2,i) = uu(0,i)*(1d0-eps)*uu(2,i)/u2
               uu(3,i) = uu(0,i)*(1d0-eps)*uu(3,i)/u2
            else
               uu(:,i) = 0d0
            end if
         end if
      end do

      end

c***********************************************************************
      subroutine add_source(tau,dt,dx,mxx,u,f,pr,vx,vz,iz,opt1,opt2)

c...Update conserved quantities u(:) by using flux f(:), and optionally
c...add source term in case of tau-eta coordinate.
      implicit none
      integer mxx,iz,i,opt1,opt2
      real*8 u(0:4,0:mxx),f(0:4,0:mxx),pr(0:mxx)
      real*8 vz(0:mxx),vx(0:mxx)
      real*8 dx,tau,dt,ch,sh

      real*8 f01,fz1,f0,fz,ch1,ch2,sh1,sh2,etot,etot2,h,ptot,ptot2,u0,uz

      if(opt1.eq.0.or.iz.eq.0) then

        if(opt2.eq.1) call limiter_timelike(mxx,f,u,vx,dt/(tau*dx))
        u(:,0) = u(:,0) -f(:,0)
        do i=1,mxx
          u(:,i) = u(:,i) + f(:,i-1) - f(:,i)
        end do

      else if(opt1.eq.1) then

c...Add source term in the Milne coordinate.
      do i = 0,mxx
        u(0,i)=u(0,i)-u(iz,i)*vz(i)/tau*dt - pr(i)*dt
c       u(0,i)=u(0,i)-u(0,i)*vz(i)**2/tau*dt - pr(i)*(1d0+vz(i)**2)*dt
        u(iz,i)=u(iz,i) -u(iz,i)/tau*dt
      end do

      if(opt2.eq.1) call limiter_timelike(mxx,f,u,vx,dt/(tau*dx))

      u(:,0) = u(:,0) -f(:,0)
      do i=1,mxx
        u(:,i) = u(:,i) + f(:,i-1) - f(:,i)
      end do

      else if(opt1.eq.2) then

c....We assume this is called only for z-direction update,and
c...assume that 1 = z-direction, 2=y-direction, 3=x-direction

        if(opt2.eq.1)
     &        call limiter_timelike_taueta(mxx,f,u,vx,dt/(tau*dx),dx)

c...K. Murase's method. See his D-thesis for a general formulation.
        ch=cosh(dx/2)
        sh=sinh(dx/2)
        u(2:,0) = u(2:,0) - f(2:,0)
        u(3:,0) = u(3:,0) - f(3:,0)
        u(4:,0) = u(4:,0) - f(4:,0)
        u(0:,0) = u(0,0) - ch*f(0,0) - sh*f(1,0)
        u(1:,0) = u(1,0) - sh*f(0,0) - ch*f(1,0)

        do i=1,mxx
          u(2,i) = u(2,i) + f(2,i-1) - f(2,i)
          u(3,i) = u(3,i) + f(3,i-1) - f(3,i)
          u(4,i) = u(4,i) + f(4,i-1) - f(4,i)
          u(0,i) = u(0,i) + ( ch*f(0,i-1) - sh*f(1,i-1))
     &                    - ( ch*f(0,i)   + sh*f(1,i  ))
          u(1,i) = u(1,i) + (-sh*f(0,i-1) + ch*f(1,i-1)) 
     &                    - ( sh*f(0,i  ) + ch*f(1,i  )) 
        end do

c       ch=cosh(dx)
c       sh=sinh(dx)
c       do i=1,mxx
c         u(3,i) = u(3,i) + f(3,i-1) - f(3,i)
c         u(2,i) = u(2,i) + f(2,i-1) - f(2,i)
c         u(0,i) = u(0,i) + ( ch*f(0,i-1) - sh*f(1,i-1)) - f(0,i)
c         u(1,i) = u(1,i) + (-sh*f(0,i-1) + ch*f(1,i-1)) - f(1,i)
c         u(4,i) = u(4,i) + f(4,i-1) - f(4,i)
c       end do

      endif


      end

c*******************************************************************
      subroutine limiter_timelike_taueta(mxx,f,uu,vx,lam,dx)

c...Flux limiter to preserve timelikeness of energy-momentum tensor
c...by K.Murase method.
      implicit none
      integer mxx,i
      double precision f(0:4,0:mxx),uu(0:4,0:mxx),vx(0:mxx),lam,v1,v2
      double precision ur(0:4),ul(0:4),ur1(0:4),ul1(0:4),ur0,ul0,dx
      double precision fr(0:4),fl(0:4),ch,sh
      double precision ur4,ul4,alpha,eps,eps1,eps3,scalef,fac
c     parameter(alpha=0.5d0,eps=0.02d0,eps3=1d-10)
      parameter(alpha=0.5d0,eps=0.00d0,eps3=1d-10)
      parameter(eps1=0d0)
      real*8 alp,alp1,alp2

c...Check.
      real*8 facsave(0:mxx)

      ch=cosh(dx/2)
      sh=sinh(dx/2)

      do i=0,mxx-1

c....Murase original.
c       ur(:)=(alpha-eps)*uu(:,i)
c       ul(:)=(1.0-alpha-eps)*uu(:,i+1)

c...Murase new.
        ur(:)=min(1.0d0,max(0.0d0,0.5+vx(i)*lam))*uu(:,i)
        ul(:)=min(1.0d0,max(0.0d0,0.5-vx(i+1)*lam))*uu(:,i+1)

c...Murase new II use the velocity of the total energy density.
c       v1=0.0; v2=0.0;
c       if(u(0,i).ne.0d0) v1=uu(1,i)/uu(0,i)
c       ur(:)=min(1.0d0,max(0.0d0,0.5+v1*lam))*uu(:,i)
c       if(u(0,i+1).ne.0d0) v2=uu(1,i+1)/uu(0,i+1)
c       ul(:)=min(1.0d0,max(0.0d0,0.5-v2*lam))*uu(:,i+1)

c       ur(:)=alp1*uu(:,i)
c       ul(:)=alp2*uu(:,i+1)


        fr(0) = ch*f(0,i) + sh*f(1,i)
        fl(0) = ch*f(0,i) - sh*f(1,i)

        ur0=ur(0)-fr(0)
        ul0=ul(0)+fl(0)
        if(ur0.lt.eps1) then
          f(:,i) = ur(:)
        else if(ul0.lt.eps1) then
          f(:,i) = -ul(:)
        endif

        fr(1)=  sh*f(0,i) + ch*f(1,i)
        fl(1)= -sh*f(0,i) + ch*f(1,i)
        fr(2)=f(2,i)
        fr(3)=f(3,i)
        fl(2)=f(2,i)
        fl(3)=f(3,i)

         ur1(0)=ur(0) - fr(0)
         ur1(1)=ur(1) - fr(1)
         ur1(2)=ur(2) - fr(2)
         ur1(3)=ur(3) - fr(3)

         ul1(0)=ul(0) + fl(0)
         ul1(1)=ul(1) + fl(1)
         ul1(2)=ul(2) + fl(2)
         ul1(3)=ul(3) + fl(3)

c       ur1(:)=ur(:)-f(:,i)
c       ul1(:)=ul(:)+f(:,i)

        ur4=ur1(0)**2-ur1(1)**2-ur1(2)**2-ur1(3)**2
        ul4=ul1(0)**2-ul1(1)**2-ul1(2)**2-ul1(3)**2
        if(ur4.ge.0d0.and.ul4.ge.0d0) cycle

        fac=max(0d0,min(scalef(fr,-ur),scalef(fl,ul))-eps3)
c       f(:,i) = fac*f(:,i)
        f(0,i) = fac*f(0,i)
        f(1,i) = fac*f(1,i)
        f(2,i) = fac*f(2,i)
        f(3,i) = fac*f(3,i)

c....Check.
        ur1(:)=ur(:)-fac*fr(:)
        ul1(:)=ul(:)+fac*fl(:)
        ur4=ur1(0)**2-ur1(1)**2-ur1(2)**2-ur1(3)**2
        ul4=ul1(0)**2-ul1(1)**2-ul1(2)**2-ul1(3)**2
        if(ur4.lt.0d0.or.ul4.lt.0d0) then
         print *,'limiter tau-eta ur<0 or ul<0?',ur4,ul4
        endif

c....Check.
c       facsave(i)=fac
c       if(i.gt.0) then
c       ur1(:) = uu(:,i) + f(:,i-1)-f(:,i)
c       ur4=ur1(0)**2-ur1(1)**2-ur1(2)**2-ur1(3)**2
c       if(ur4.lt.0d0) then
c        print *,'u<0 ?',i,ur4,ur1(0),ur1(1),ur1(2),ur1(3)
c        print *,'fac=',facsave(i-1),fac
c        print *,'alp1 alp2=',alp1,alp2
c       endif
c       endif

      end do

      return

c...check
      do i=1,mxx-1
        ur1(:) = uu(:,i) + f(:,i-1)-f(:,i)
        ur4=ur1(0)**2-ur1(1)**2-ur1(2)**2-ur1(3)**2
        if(ur1(0).eq.0d0) cycle
c       if(ur1(0).le.1d-15.or.ur4.lt.0d0) then
        if(ur4.lt.0d0) then
         print *,'timelike: u<0 ?',i,ur4,ur1(0),ur1(1),ur1(2),ur1(3)
         print *,'u=',uu(0,i),uu(1,i),uu(2,i),uu(3,i),uu(4,i)
         print *,'f=',f(0,i-1),f(1,i-1),f(2,i-1),f(3,i-1),f(4,i-1)
         print *,'f2=',f(0,i),f(1,i),f(2,i),f(3,i),f(4,i)
        endif
      end do

      end
c***********************************************************************
c...Update thermo dynamical quantities from uu
      subroutine update_local_val(tau,mxx,uu,vv,sv,er,pr,nr,tr,mr,uv,uo,
     &                             opt_rescale)
      implicit none
      integer mxx,ix,i,j,opt_rescale
      real*8 tau,vel
      real*8 uu(0:4,0:mxx),uv(0:4),uo(0:4)
      real*8 vv(3,0:mxx),sv(0:mxx)
      real*8 er(0:mxx),pr(0:mxx),nr(0:mxx),tr(0:mxx),mr(0:mxx)
      real*8 uh(0:4),det,va,gam,fg,u2,ut2,ut
      real*8 mim,eps
      parameter(mim = 0.00001d0)
      real*8 ene,nba,pre,sva,tplc,mulc,mus,slc

      uv=0.0d0; uo=0.0d0
c     uo(:) = uo(:) + uu(:,1)
c     uo(:) = uo(:) + uu(:,mxx-1)
c     uo(:) = uo(:) + uu(:,0)
c     uo(:) = uo(:) + uu(:,mxx)
c     uu(:,0)=0d0; uu(:,1)=0d0; uu(:,mxx)=0d0; uu(:,mxx-1)=0d0
c-----------------------------------------------------------------------      
      if(opt_rescale.ge.1) then
      eps=1d-15
c     eps=1d-12
      do i = 0,mxx
         u2=sqrt(uu(1,i)**2+uu(2,i)**2+uu(3,i)**2)
         if (uu(0,i).lt.u2) then
            if (u2.gt.1d-16) then
               uv(:)=uv(:)+uu(:,i)
c...Scale energy.
               if(opt_rescale.eq.1) then
c                uu(0,i)=u2*(1d0+eps)
                 uu(0,i)=u2+1d-12
               else if(opt_rescale.eq.2) then
c...Scale momentum.
                 uu(1,i)=uu(0,i)*(1d0-eps)*uu(1,i)/u2
                 uu(2,i)=uu(0,i)*(1d0-eps)*uu(2,i)/u2
                 uu(3,i)=uu(0,i)*(1d0-eps)*uu(3,i)/u2
c...Scale transverse momentum.
               else
                 ut2=uu(0,i)**2-uu(3,i)**2
                 if(ut2.ge.0d0) then
                   ut=sqrt(uu(1,i)**2+uu(2,i)**2)
                   uu(1,i)=sqrt(ut2)*(1d0-eps)*uu(1,i)/ut
                   uu(2,i)=sqrt(ut2)*(1d0-eps)*uu(2,i)/ut
                 else
                   print *,'u0<uz?',uu(0,i)**2-uu(3,i)**2
                   uu(0,i)=u2+1d-12
                 endif
               endif

               uv(:)=uv(:)-uu(:,i)
            else
               uv(:)=uv(:)+uu(:,i)
               uu(:,i) = 0d0
            end if
         end if
         if(uu(0,i).gt.1d+5) then
          print *,'energy big?',i,(uu(j,i),j=0,4)
         endif
      end do

      endif
c-----------------------------------------------------------------------      

      do ix = 0,mxx

        uh(:) = uu(:,ix)/tau

c       if (uh(0).lt.0d0.or.sqrt(uh(0)*uh(0)) .lt. mim) uh=0d0
        if(uh(0).le.1d-15) uh=0d0

        if(uh(4).lt.0d0) then
c         print *,'baryon density<0?',uu(4,ix),uu(0,ix)
          uu(4,ix)=0d0; uh(4)=0d0
        endif
  
        va=uh(1)*uh(1)+uh(2)*uh(2)+uh(3)*uh(3)
        det = uh(0)*uh(0)-va
        if (det .le. 0d0) then
c         if(det.ne.0d0) print *,'det<0?',ix,det,uh(0),uu(0,ix)
          uv(:)=uv(:)+uu(:,ix)
          er(ix) = 0d0
          pr(ix) = 0d0
          nr(ix) = 0d0
          tr(ix) = 0d0
          mr(ix) = 0d0
          uu(:,ix) = 0d0
          vv(:,ix) = 0d0
          sv(ix) = 0d0
 
        else
 
        call thermal(uh,ene,nba,pre,sva,vel)
        call geteos2(ene,nba,tplc,mulc,mus,slc)
        er(ix) = ene
        nr(ix) = nba
        pr(ix) = pre
        tr(ix) = tplc
        mr(ix) = mulc
c       fg = uh(0) + tau*pre
        fg = uh(0) + pre

        gam = 1d0-va/fg/fg
        if((fg .gt.0d0) .and. (gam .gt. 0d0)) then
          gam = sqrt(1d0/gam)
        else
          gam = 1d0
        endif
c       nr(ix) = uh(4)/gam/tau
        nr(ix) = uh(4)/gam

        if(fg .ne. 0d0) then
          vv(1,ix) = uh(1)/fg
          vv(2,ix) = uh(2)/fg
          vv(3,ix) = uh(3)/fg
          sv(ix)=sva
        else
          vv(1,ix) = 0d0
          vv(2,ix) = 0d0
          vv(3,ix) = 0d0
          sv(ix) = 0d0
        end if

c       gam = 1d0/sqrt(1d0-va/fg/fg)
c       nr(ix) = uh(4)/gam/tau
        if(abs(nr(ix)-nba).gt.1d-5) then
          u2=sqrt(uu(1,ix)**2+uu(2,ix)**2+uu(3,ix)**2)
          print *,ix,'nr?',nr(ix),nba,' gam=',gam,
     &    ' v=',sqrt(va/fg/fg),
     &    'det=',uu(0,ix)**2-u2**2,'u0=',uu(0,ix),'u2=',u2,
     &    'pre=',pre,'ene=',ene,'fg=',fg 
        endif


        endif

      end do

      end

c***********************************************************************
c...Update thermo dynamical quantities from uu
      subroutine update_all(tau)
      implicit none
      include 'fluid.inc'
      integer ix,iy,iz
      real*8 tau
      real*8 uh(0:4),det,va,gam,fg,vel
      real*8 mim
      parameter(mim = 0.00001d0)
      real*8 ene,nba,pre,sva,tplc,mulc,mus,slc

      do ix = 0,maxx
      do iy = 0,maxy
      do iz = 0,maxz

        uh(:) = u(:,ix,iy,iz)/tau
        if (uh(0).lt.0d0.or.sqrt(uh(0)*uh(0)) .lt. mim) then
          uh=0d0
        endif
        va=uh(1)*uh(1)+uh(2)*uh(2)+uh(3)*uh(3)
        det = uh(0)*uh(0)-va
        if (det .le. 0d0) then
          el(ix,iy,iz) = 0d0
          pl(ix,iy,iz) = 0d0
          bl(ix,iy,iz) = 0d0
          tl(ix,iy,iz) = 0d0
          ml(ix,iy,iz) = 0d0
          u(:,ix,iy,iz) = 0d0
          vx(ix,iy,iz) = 0d0
          vy(ix,iy,iz) = 0d0
          vz(ix,iy,iz) = 0d0
          cs(ix,iy,iz) = 0d0
          fg=0d0
 
        else
 
c       taufutr=tau
c       if(opt_taueta.eq.1) taufutr=tau+dt/3d0
        call thermal(uh,ene,nba,pre,sva,vel)
        call geteos2(ene,nba,tplc,mulc,mus,slc)
        el(ix,iy,iz) = ene
        bl(ix,iy,iz) = nba
        pl(ix,iy,iz) = pre
        tl(ix,iy,iz) = tplc
        ml(ix,iy,iz) = mulc
c       fg = uh(0) + tau*pre
        fg = uh(0) + pre
        gam = 1d0-va/fg/fg
        if((fg .gt.0d0) .and. (gam .gt. 0d0)) then
          gam = sqrt(1d0/gam)
        else
          gam = 1d0
        endif
c       bl(ix,iy,iz) = uh(4)/gam/tau
        bl(ix,iy,iz) = uh(4)/gam

        endif

        if(fg .ne. 0d0) then
          vx(ix,iy,iz) = uh(1)/fg
          vy(ix,iy,iz) = uh(2)/fg
          vz(ix,iy,iz) = uh(3)/fg
          cs(ix,iy,iz)=sva
        else
          vx(ix,iy,iz) = 0d0
          vy(ix,iy,iz) = 0d0
          vz(ix,iy,iz) = 0d0
          cs(ix,iy,iz) = 0d0
        end if

      end do
      end do
      end do

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine thermal(uu,elc,nblc,plc,cslc,v)

c...Compute thermodynamical values from uu()

      implicit none
      double precision uu(0:4)
      double precision vv1,vv2,vv3,vv4,vv5
      double precision elc,nblc,plc,cslc,m,v,vp,vm,fv
      double precision tplc,mulc
      double precision v0,eps,vf
      integer iv0,iv,mxit,icon
      parameter(mxit=30,eps=1d-8)
      real*8 g1,g2,fe

      icon=0
      vv1 = uu(1)
      vv2 = uu(2)
      vv3 = uu(3)
      vv4 = uu(0)
      vv5 = uu(4)
      m = sqrt(vv1*vv1+vv2*vv2+vv3*vv3)

      if(m.eq.0d0) then
c     if(m.lt.1d;-8) then
        elc=vv4
        nblc=vv5
        call geteos3(vv4,nblc,plc,cslc)
        return
      endif
c     if(vv4.lt.1d-8) then
      if(vv4.lt.1d-12) then
        elc=0d0
        nblc=0d0
        plc=0d0
        cslc=0d0
c       print *,'thermal::vv4 too small?',vv4
        return
      endif

c------------------------------------------------------
      call iteratev(vv4,vv5,m,elc,nblc,plc,cslc,v,icon)

c     if(plc.eq.0d0) then
c        print *,'after iterative p=0',plc,elc,nblc
c        print *,'vv4=',vv4,vv5,m,vf
c        print *,'icon=',icon,vv4-m
c        elc=0d0; nblc=0d0; plc=0d0; cslc=0d0
c        return
c     endif

c     plc=m/vf-uu(0) ! just in case
        fe=vv4+plc
        g1=1d0-m*m/fe/fe
        if(g1.le.0d0) then
         print *,'g1=0?',g1,m,fe,' p=',plc,' vv0=',uu(0),
     6   uu(1),uu(2),uu(3),
     &   'v=',v,
     &   ' v=',  m / (vv4  + plc ),
     &   ' v2=',  m / fe
        
         return
        endif

        if(v.ge.1d0-1d-20) then
          print *,'after iteratev v=',v,icon
        endif

c       g1=1.0d0/sqrt(max(0d0,1d0-m*m/fe/fe))
c       g2=1.0d0/sqrt(1d0-v*v)
c       if(abs(g1-g2)/g1.gt.1d-3) then
c       print *,'iter ene=',elc,plc,v
c       print *,'vv4=',vv4,m,vv4-m
c        print *,'g1=',g1,' g2=',g2,icon
c        print *,'v=',m/fe,v
c       endif

      if(icon.eq.0) return

c     if((vv4+plc)*v-m.gt.0d0) then
c       vp=v
c     else
c       vm=v
c     endif
c------------------------------------------------------
200   continue
      vp = 1d0
      vm = 0d0
      v0=-1.0
      do 100 iv = 1,mxit
      v = (vp+vm)/2.
c     if(abs(v-v0).lt.eps) return
      if(abs(v-v0).lt.eps) then
      endif

      v0=v
      elc = vv4-m*v
      nblc = vv5*sqrt(1d0-v*v)
      if(elc.lt.0d0) then
        print *,'thermal: elc<0',elc,'vv4=',vv4,'m=',m,'iv=',iv
        stop
      endif
      call geteos3(elc,nblc,plc,cslc)
      if(vv4+plc .eq. 0d0) then
        print *,'thermal vv4+plc is zero',vv4+plc
      stop
      endif
      fv = v-m/(vv4+plc)
      if(fv .gt. 0d0)then
      vp = v      
      else
      vm = v
      endif

      if(abs(vp-vm).lt.eps) then
        if(abs(v-m/(uu(0)+plc)).gt.1d-2) then
          print *,'plc?',plc,m/v-uu(0),abs(m/v-uu(0)-plc)
        endif
c       fe=uu(0)+plc
c       g1=1.0d0/sqrt(1d0-m*m/fe/fe)
c       g2=1.0d0/sqrt(1d0-v*v)
c       if(abs(g1-g2).gt.1d-2) then
c       print *,'bisec ene=',elc,plc,v,m
c        print *,'g1=',g1,' g2=',g2
c       endif
        return
      endif

 100  continue

      print *,'not converge? v=',iv0,v,vp,vm

      end

c***********************************************************************
      subroutine iteratev(u0,u4,m,elc,nblc,plc,cslc,vf,icon)
      implicit none
      real*8 u0,u4,eps,vf,elc,m,v,nblc,plc,cslc
      integer iv,icon,mxit
c     parameter(eps=1d-7,mxit=50)
      parameter(eps=1d-7,mxit=20)

      icon=0
      iv = 0
      vf = 0d0
 10   v = vf
      elc = u0-m*v
      if(elc.lt.0d0) then
        print *,'iteratev::elc<0',elc,'u0=',u0,'m=',m,'iv=',iv
        icon=2
        return
      endif
      nblc = u4*sqrt(max(0d0, 1d0 - v*v) )
      call geteos3(elc,nblc,plc,cslc)
      vf = m / (u0  + plc )
      iv = iv + 1
      if (abs((v-vf)/vf).lt.eps) return
      if(iv.gt.mxit) then
        icon=1
        return
      endif
      goto 10

      end

c***********************************************************************
      subroutine hydroTotalEnergy(it,tau,etot,pxtot,pytot,pztot,btot)
      implicit none
      include 'fluid.inc'
      double precision etot,btot,tau
      integer ix,iy,iz,mnx,mny,mxy,it
      double precision yv,hh,fh
      double precision etot0,px0,py0,pz0,btot0,entro0
      double precision pxtot,pytot,pztot,entro,vol,volh
      double precision esave,pxsave,pysave,pzsave,bsave
      double precision pxerr,pyerr,pzerr,peerr,berr,serr
      double precision y,h,ch
      double precision v2,gam,getentropy
c     integer first
c     data first/0/
c     save first
c     save etot0,px0,py0,pz0,btot0


      etot0 = pe_tot
      px0   = px_tot
      py0   = py_tot
      pz0   = pz_tot
      btot0 = b_tot
      entro0 = en_tot

      etot = 0d0
      pxtot = 0d0
      pytot = 0d0
      pztot = 0d0
      btot  = 0d0

      entro=0d0
      vol=dx*dy*dz
      volh=vol*hbc

      do ix = 0,maxx
      do iy = 0,maxy
      do iz = 0,maxz
        pxtot = pxtot + u(1,ix,iy,iz)*volh
        pytot = pytot + u(2,ix,iy,iz)*volh
        btot  = btot  + u(4,ix,iy,iz)*vol

        if(opt_taueta.eq.0) then
          etot  = etot  + u(0,ix,iy,iz)*volh
          pztot = pztot + u(3,ix,iy,iz)*volh
        else 

          h=dz*(iz-origz)

c         y=0.5*log((1d0+vz(ix,iy,iz))/(1d0-vz(ix,iy,iz))) + h
c         ch=cosh(y)/cosh(y-h)
c         etot = etot + u(0,ix,iy,iz)*ch*volh
c    &                + tau*pl(ix,iy,iz)*sinh(y)*tanh(y-h)*volh
c... cosh(y)/cosh(y-h) - cosh(h) = sinh(y)*tanh(y-h)
c         pztot = pztot + u(0,ix,iy,iz)*ch*tanh(y-h)*volh
c    &                + tau*pl(ix,iy,iz)*(ch*tanh(y-h)+sinh(y))*volh

          etot  = etot + volh*(
     &      cosh(h)*u(0,ix,iy,iz)+sinh(h)*u(3,ix,iy,iz) )

          pztot = pztot + volh*(
     &      sinh(h)*u(0,ix,iy,iz)+cosh(h)*u(3,ix,iy,iz) )

c         etot  = etot  + u(0,ix,iy,iz)*volh
c         pztot = pztot + u(3,ix,iy,iz)*volh

c           if(u(0,ix,iy,iz).ne.0d0) then
c         print *,'e',u(0,ix,iy,iz)*ch
c    &                + tau*pl(ix,iy,iz)*sinh(y)*tanh(y-h)
c    &     ,cosh(h)*u(0,ix,iy,iz)+sinh(h)*u(3,ix,iy,iz),
c    &  u(0,ix,iy,iz)*ch+ tau*pl(ix,iy,iz)*sinh(y)*tanh(y-h)
c    &  -cosh(h)*u(0,ix,iy,iz)-sinh(h)*u(3,ix,iy,iz)

c         print *,'p',u(0,ix,iy,iz)*ch*tanh(y-h)
c    &                + tau*pl(ix,iy,iz)*(ch*tanh(y-h)+sinh(y))
c    &     -sinh(h)*u(0,ix,iy,iz)-cosh(h)*u(3,ix,iy,iz)

c         print *,u(0,ix,iy,iz)*ch*tanh(y-h)
c    &                + tau*pl(ix,iy,iz)*(ch*tanh(y-h)+sinh(y))
c    &     ,sinh(h)*u(0,ix,iy,iz)+tau*cosh(h)*u(3,ix,iy,iz)

c            endif

        endif

        if(u(0,ix,iy,iz).gt.1d-12) then
          v2=vx(ix,iy,iz)**2+vy(ix,iy,iz)**2+vz(ix,iy,iz)**2
          if(v2.lt.1d0) then
          gam=1.0/sqrt(1.0d0-v2)
          entro = entro + getentropy(el(ix,iy,iz),bl(ix,iy,iz))*vol*gam
          else
            print *,'hydroTotalenergy v>1?',v2
          endif
        endif

      end do
      end do
      end do


c     first=first+1
c     if(first.eq.1) then
c       etot0=etot
c       px0=pxtot
c       py0=pytot
c       pz0=pztot
c       btot0=btot
c       entro0=entro
c     endif

        peerr = abs((etot-etot0)/etot0)*100
        pzerr = abs((pztot-pz0)/max(1d0,pz0))*100
        berr  = abs((btot-btot0)/max(1d0,btot0))*100
        serr  = abs((entro-entro0)/max(1d0,entro0))*100

      write(6,*)'======================================================'
        write(6,9805)it,tau,etot,btot,berr,serr
9805  format('*h',i4,1x,f6.3,'(fm/c)',' E=',e9.4,' B=',e9.4,
     & ' B (%)',e9.4,' S (%)',e12.6)

c       write(6,9803)entro0,entro,serr
c9803 format('entropy0= ',e10.5,' entropy=',e10.5, ' S (%)',f10.6)

        write(6,804)uzero(0)*volh,uzero(1)*volh,uzero(2)*volh,
     & uzero(3)*volh
804   format('pzero = ',4(e12.5,1x))

        write(6,805)uout(0)*volh,uout(1)*volh,uout(2)*volh,uout(3)*volh
805   format('p_out = ',4(e12.5,1x))

        write(6,801)pxtot,pytot,pztot,peerr
801   format('p_total = ',3(e12.5,1x),' econ(%)',e12.6)

      peerr = abs((etot+uzero(0)*volh-etot0)/etot0)*100
      write(6,803)px0-pxtot+uzero(1)*volh,
     &            py0-pytot+uzero(2)*volh,
     &            pz0-pztot+uzero(3)*volh,peerr
803   format('p_tot+pzero       = ',3(e12.5,1x),'econ(%)',e12.6)

        peerr = abs((etot+uout(0)*volh-etot0)/etot0)*100
      write(6,806)px0-pxtot+uout(1)*volh,
     &            py0-pytot+uout(2)*volh,
     &            pz0-pztot+uout(3)*volh,peerr
806   format('p_tot+p_out       = ',3(e12.5,1x),'econ(%)',e12.6)

       etot=etot+(uout(0)+uzero(0))*volh
       pxtot=pxtot+(uout(1)+uzero(1))*volh
       pytot=pytot+(uout(2)+uzero(2))*volh
       pztot=pztot+(uout(3)+uzero(3))*volh

        peerr = abs((etot-etot0)/etot0)*100
      write(6,807)px0-pxtot,
     &            py0-pytot,
     &            pz0-pztot,peerr
807   format('p_tot+p_out+pzero = ',3(e12.5,1x),'econ(%)',e12.6)

      end

c************************************************************************

      subroutine checkout(t,check)
      implicit none
      include 'fluid.inc'
      integer check
      integer ix,iy,iz,n
      real*8 emax,tmax,pmax,bmax,t,aveb,r,edens0,avet,avemu
      common/fluid1/emax,tmax,pmax,bmax,aveb

      check=0

      emax=-100d0
      tmax=-100d0
      pmax=-100d0
      bmax=-100d0
      aveb=0d0
      avet=0d0
      avemu=0d0
      r=0d0
      n=2
      do ix = 2,maxx-n
      do iy = 2,maxy-n
      do iz = 2,maxz-n
         pmax=max(pmax,pl(ix,iy,iz))
         tmax=max(tmax,tl(ix,iy,iz))
         emax=max(emax,el(ix,iy,iz))
         bmax=max(bmax,bl(ix,iy,iz))
         aveb=aveb+bl(ix,iy,iz)*bl(ix,iy,iz)
         avet=avet+tl(ix,iy,iz)*bl(ix,iy,iz)
         avemu=avemu+ml(ix,iy,iz)*bl(ix,iy,iz)
         r=r+bl(ix,iy,iz)
        if(el(ix,iy,iz) .ge. efreeze) then
          check = check +1
        endif
      end do
      end do
      end do

      if(t.le.t_pass) check=1

      if(r.gt.0d0) then
        aveb=aveb/r
        avet=avet/r*hbc
        avemu=avemu/r*hbc
      endif
      emax=emax*hbc
      pmax=pmax*hbc
      tmax=tmax*hbc
      edens0=el(maxx/2,maxy/2,maxz/2)*hbc

      if(job_mode.ge.1) then

      write(6,800)emax,bmax/0.16,pmax,tmax*1000,check
800   format('e_max=',f9.4,' GeV/fm^3, b_max/b0=',f8.4,', p_max=',f8.4,
     & ' Tmax=',f9.4,' MeV ',i5)

      write(6,810) edens0,aveb/0.16,avet,avemu,t_pass
810   format('e(0,0,0)=',f9.4,' ave. B dens/n0 = ',f6.3,
     &    ' <T>=',f8.3,' <mu>=',f8.3,
     &   ' passing time=',f8.3)

       endif

c     print *,'check=',check,efreeze*hbc,t_pass,t

      end


c************************************************************************

      subroutine checkboundary_x
      implicit none
      include 'fluid.inc'
      integer i,ix,iy,iz
      double precision h

c     do i=0,1
c       ix=(maxx-2)*i+1
c     do iy = 1,maxy-1
c     do iz = 1,maxz-1
c       uout(:) = uout(:) + u(:,ix,iy,iz)
c       u(:,ix,iy,iz)=0.0d0
c     end do
c     end do
c     end do

      do i=0,1
       ix=maxx*i
      do iy = 0,maxy
      do iz = 0,maxz
        if(opt_taueta.eq.0) then
        uout(:) = uout(:) + u(:,ix,iy,iz)
        else 
          h=dz*(iz-origz)
        uout(1) = uout(1) + u(1,ix,iy,iz)
        uout(2) = uout(2) + u(2,ix,iy,iz)
        uout(4) = uout(4) + u(4,ix,iy,iz)
        uout(0) = uout(0) + cosh(h)*u(0,ix,iy,iz)+sinh(h)*u(3,ix,iy,iz)
        uout(3) = uout(3) + sinh(h)*u(0,ix,iy,iz)+cosh(h)*u(3,ix,iy,iz)
        endif
        u(:,ix,iy,iz)=0.0d0
      end do
      end do
      end do

      end

c************************************************************************

      subroutine checkboundary_y
      implicit none
      include 'fluid.inc'
      integer i,ix,iy,iz
      double precision h

c     do i=0,1
c       iy=(maxy-2)*i+1
c     do ix = 1,maxx-1
c     do iz = 1,maxz-1
c       uout(:) = uout(:) + u(:,ix,iy,iz)
c       u(:,ix,iy,iz)=0.0d0
c     end do
c     end do
c     end do

      do i=0,1
       iy=maxy*i
      do ix = 0,maxx
      do iz = 0,maxz
        if(opt_taueta.eq.0) then
        uout(:) = uout(:) + u(:,ix,iy,iz)
        else 
          h=dz*(iz-origz)
        uout(1) = uout(1) + u(1,ix,iy,iz)
        uout(2) = uout(2) + u(2,ix,iy,iz)
        uout(4) = uout(4) + u(4,ix,iy,iz)
        uout(0) = uout(0) + cosh(h)*u(0,ix,iy,iz)+sinh(h)*u(3,ix,iy,iz)
        uout(3) = uout(3) + sinh(h)*u(0,ix,iy,iz)+cosh(h)*u(3,ix,iy,iz)
        endif
        u(:,ix,iy,iz)=0.0d0
      end do
      end do
      end do

      end

c************************************************************************

      subroutine checkboundary_z
      implicit none
      include 'fluid.inc'
      integer i,ix,iy,iz
      double precision h

c     do i=0,1
c       iz=(maxz-4)*i+2
c     do ix = 2,maxx-2
c     do iy = 2,maxy-2
c       uout(:) = uout(:) + u(:,ix,iy,iz)
c       u(:,ix,iy,iz)=0.0d0
c     end do
c     end do
c     end do

c     do i=0,1
c       iz=(maxz-2)*i+1
c     do ix = 1,maxx-1
c     do iy = 1,maxy-1
c       uout(:) = uout(:) + u(:,ix,iy,iz)
c       u(:,ix,iy,iz)=0.0d0
c     end do
c     end do
c     end do

      do i=0,1
       iz=maxz*i
          h=dz*(iz-origz)
      do ix = 0,maxx
      do iy = 0,maxy
        if(opt_taueta.eq.0) then
        uout(:) = uout(:) + u(:,ix,iy,iz)
        else 
        uout(1) = uout(1) + u(1,ix,iy,iz)
        uout(2) = uout(2) + u(2,ix,iy,iz)
        uout(4) = uout(4) + u(4,ix,iy,iz)
        uout(0) = uout(0) + cosh(h)*u(0,ix,iy,iz)+sinh(h)*u(3,ix,iy,iz)
        uout(3) = uout(3) + sinh(h)*u(0,ix,iy,iz)+cosh(h)*u(3,ix,iy,iz)
        endif
        u(:,ix,iy,iz)=0.0d0
      end do
      end do
      end do

      end

c************************************************************************

      subroutine checkboundary
      implicit none
      include 'fluid.inc'
      integer i,ix,iy,iz,j
      real*8 x,y,z

      do i=0,1
        ix=(maxx-2)*i+1
      do iy = 1,maxy-1
      do iz = 1,maxz-1
        uout(:) = uout(:) + u(:,ix,iy,iz)
        u(:,ix,iy,iz)=0.0d0
      end do
      end do
      end do

      do i=0,1
        iy=(maxy-2)*i+1
      do ix = 1,maxx-1
      do iz = 1,maxz-1

c       if(u(0,ix,iy,iz).gt.0d0) then
c         x = dx*(ix-origx)
c         y = dy*(iy-origy)
c         z = dz*(iz-origz)
c         print *,'out y',x,y,z
c         print *,ix,iy,iz
c         print *,(u(j,ix,iy,iz),j=0,3)
c         print *,u(0,ix,iy,iz)+u(3,ix,iy,iz)
c         pause
c       endif

        uout(:) = uout(:) + u(:,ix,iy,iz)
        u(:,ix,iy,iz)=0.0d0
      end do
      end do
      end do

      do i=0,1
        iz=(maxz-2)*i+1
      do ix = 1,maxx-1
      do iy = 1,maxy-1

c       if(u(0,ix,iy,iz).gt.0d0) then
c         x = dx*(ix-origx)
c         y = dy*(iy-origy)
c         z = dz*(iz-origz)
c         print *,'out z',x,y,z,u(0,ix,iy,iz),u(3,ix,iy,iz)
c       endif

        uout(:) = uout(:) + u(:,ix,iy,iz)
        u(:,ix,iy,iz)=0.0d0
      end do
      end do
      end do


      do i=0,1
       ix=maxx*i
      do iy = 0,maxy
      do iz = 0,maxz
        u(:,ix,iy,iz)=0.0d0
      end do
      end do
      end do

      do i=0,1
       iy=maxy*i
      do ix = 0,maxx
      do iz = 0,maxz
        u(:,ix,iy,iz)=0.0d0
      end do
      end do
      end do

      do i=0,1
       iz=maxz*i
      do ix = 0,maxx
      do iy = 0,maxy
        u(:,ix,iy,iz)=0.0d0
      end do
      end do
      end do

      end

c************************************************************************

      subroutine directedFlow(msel,it,pxdirtot)
      implicit none
      include 'fluid.inc'
      real*8 v1,v2,v3,fh,gam,rap,nbar
      real*8 ymin,ymax,emnuc,px,v,drap
      real*8 pxdirtot(2),ntot(2)
      integer ny,it,iy,ix,iz,ir,msel
      logical first
      parameter(ny=30)
      parameter(emnuc=0.938)
      real*8 pxdir(0:ny),r(0:ny)
      save pxdir,r

      if(msel.eq.1) then
        ymax=ycm*1.3
        ymin=-ymax
        drap=(ymax-ymin)/ny
        do ir=1,ny
          pxdir(ir)=0.0
          r(ir)=0.0
        end do

      else if(msel.eq.2) then

      do 113 ix = 0,maxx
      do 213 iy = 0,maxy
      do 313 iz = 0,maxz
         nbar = u(4,ix,iy,iz)
         if(nbar .lt.1d-8) goto 313
         v1 = vx(ix,iy,iz)
         v2 = vy(ix,iy,iz)
         v3 = vz(ix,iy,iz)
         v=v1*v1+v2*v2+v3*v3
c        if(v.ge.1d0) goto 313
         gam=1.0/sqrt(1.0d0-v)
         rap=0.5*log( (1.0+v3)/(1.0-v3) )
         ir=nint((rap-ymin)/drap)
         if(ir.ge.0.and.ir.le.ny) then
           r(ir)=r(ir)+nbar
           pxdir(ir)=pxdir(ir)+nbar*emnuc*gam*v1
         endif

313   end do
213   end do
113   end do

      else

      open(unit=12,file='pdir.dat',status='unknown')
      write(12,*)'# it= ',it,' time= ',it*dt
      write(12,*)'# rapidity    <p_x/N>'
      write(6,*)'# it= ',it,' time= ',it*dt
      write(6,*)'# rapidity    <p_x/N>'
      pxdirtot(1)=0d0
      pxdirtot(2)=0d0
      ntot(1)=0d0
      ntot(2)=0d0
      do iy=0,ny
        rap=iy*drap + ymin
        px=0.0
        if(r(iy).gt.0d0) px=pxdir(iy)/r(iy)
        write(12,800) rap,px
        write(6,800) rap,px
        if(abs(rap).le.ycm) then
          ntot(1)=ntot(1)+r(iy) 
          pxdirtot(1)=pxdirtot(1)+r(iy)*px*sign(1d0,rap)
        endif
        if(abs(rap).le.0.25) then
          ntot(2)=ntot(2)+r(iy) 
          pxdirtot(2)=pxdirtot(2)+r(iy)*px*sign(1d0,rap)
        endif
      end do
      if(ntot(1).gt.0d0) pxdirtot(1)=pxdirtot(1)/ntot(1)
      if(ntot(2).gt.0d0) pxdirtot(2)=pxdirtot(2)/ntot(2)
      write(12,810)pxdirtot(1),pxdirtot(2)
      write(6,810)pxdirtot(1),pxdirtot(2)

      close(12)

      endif
800   format(f8.3,2x,e13.5)
810   format('# <px/N>(y<|y_cm|) = ',f8.4,' GeV',
     &   ' <px/N>(y<|0.25|)= ',f8.4,' GeV')

      end

c***************************************************************************
      subroutine outputhxy(it)
      implicit none
      include 'fluid.inc'
      integer ix,iy,iz,it
      double precision x,y,z

      open(unit=10,file='dens.dat',status='unknown')

      write(10,801)it*dt
 801  format('# time= ',f8.3)

      iz = origz
      do ix = 0,maxx
      do iy = 0,maxy
        x = dx*(ix-origx)
        y = dy*(iy-origy)
       write(10,800) x,y,u(4,ix,iy,iz),u(0,ix,iy,iz),bl(ix,iy,iz)
      end do
       write(10,*)' '
      end do
 800  format(2(f8.3,1x),3x,3(e15.7,6x))

      close(10)

       end

c***************************************************************************
      subroutine outputhzx(it)
      implicit none
      include 'fluid.inc'
      integer ix,iy,iz,it
      double precision x,y,z

      open(unit=10,file='dens.dat',status='unknown')

      write(10,801)it*dt
 801  format('# time= ',f8.3)

      iy = origy         
      do ix = 0,maxx
      do iz = 0,maxz
        x = dx*(ix-origx)
        z = dz*(iz-origz)
       write(10,800) z,x,u(4,ix,iy,iz),u(0,ix,iy,iz),bl(ix,iy,iz)
      end do
       write(10,*)' '
      end do
 800  format(2(f8.3,1x),3x,3(e15.7,6x))

      close(10)

       end

