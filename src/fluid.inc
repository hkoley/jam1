
      integer maxx,maxy,maxz,origx,origy,origz
c     parameter (maxx=80,maxy=80,maxz=80)
      parameter (maxx=100,maxy=100,maxz=100)
c     parameter (maxx=100,maxy=100,maxz=135)
c     parameter (maxx=100,maxy=100,maxz=200)

c     parameter (maxx=130,maxy=130,maxz=130)
c     parameter (maxx=150,maxy=150,maxz=150)
c     parameter (maxx=200,maxy=200,maxz=200)
      parameter (origx=maxx/2,origy=maxy/2,origz=maxz/2)

      double precision u,uzero,uout
      integer ifreeze
      common/hydro3/u(0:4,0:maxx,0:maxy,0:maxz),
     &             uzero(0:4),uout(0:4),
     &             ifreeze(0:maxx,0:maxy,0:maxz)

c....Local thermodaynamic variables
      double precision el,pl,bl,vx,vy,vz,cs,ml,tl
      common/hydro4/bl(0:maxx,0:maxy,0:maxz), ! baryon density
     $              el(0:maxx,0:maxy,0:maxz), ! energy density
     $              pl(0:maxx,0:maxy,0:maxz), ! pressure
     $              tl(0:maxx,0:maxy,0:maxz), ! temperature
     $              ml(0:maxx,0:maxy,0:maxz), ! chemical potential
     $              vx(0:maxx,0:maxy,0:maxz), ! v_x
     $              vy(0:maxx,0:maxy,0:maxz), ! v_y
     $              vz(0:maxx,0:maxy,0:maxz), ! v_z
     $              cs(0:maxx,0:maxy,0:maxz)  ! sound velocity

c...Total energy and baryon number in the fluid.
      real*8 pe_tot,px_tot,py_tot,pz_tot,b_tot,s_tot,c_tot,en_tot
      common/hydro5/pe_tot,px_tot,py_tot,pz_tot,b_tot,s_tot,c_tot,en_tot

c...Total energy and baryon number which are converted into particles.
      real*8 fconv
      common/hydro6/fconv(0:4)

c....Constant values.
      double precision hbc
      parameter(hbc=0.197327d0)

c...Parameters for hydro
      double precision dx,dy,dz,dt,tau0,efreeze,t_pass,tfreezecut
      common/param_hydro1/dx,dy,dz,dt,tau0,efreeze,t_pass,tfreezecut

      double precision ycm,atomic,gsigma,dcut
      common/param_hydro2/ycm,atomic,gsigma,dcut

      integer ntest,h_mode,eos_mode,second_order,
     & opt_timelike,opt_rescale,
     & opt_freezeout,opt_taueta,opt_table,
     & job_mode,n_gaussp

      logical op_splitting
      common/param_hydro2/ntest,h_mode,eos_mode,second_order,
     & op_splitting,
     & opt_timelike,opt_rescale,opt_freezeout,opt_taueta,opt_table,
     & job_mode,n_gaussp

c... for freezeout
      integer Ftime,FX,FY,FZ
c     parameter (Ftime=2,FX=2,FY=2,FZ=2)
c     parameter (Ftime=1,FX=2,FY=2,FZ=2)
      parameter (Ftime=1,FX=1,FY=1,FZ=1)
      double precision fe,fp,fd,fvx,fvy,fvz
      common/param_freezeout/fe(0:1,0:maxx,0:maxy,0:maxz),  ! energy densiy
     $                       fp(0:1,0:maxx,0:maxy,0:maxz),  ! pressure
     $                       fd(0:1,0:maxx,0:maxy,0:maxz),  ! baryon density
     $                       fvx(0:1,0:maxx,0:maxy,0:maxz), ! v_x
     $                       fvy(0:1,0:maxx,0:maxy,0:maxz), ! v_y
     $                       fvz(0:1,0:maxx,0:maxy,0:maxz)  ! v_z


      integer nfreezeout,frr,nfmax
      parameter(nfmax=2*maxx*maxy*maxz)

      double precision frt,frdsgm,frtmp,frv,frmuq,frmus,frnum,frbdn
      common/param_freezeout2/nfreezeout,frr(3,nfmax)

      common/param_freezeout3/frt(nfmax),         ! freezeout time
     $                        frv(3,nfmax),       ! flow velocity
     $                        frdsgm(0:3,nfmax),  ! surface
     $                        frtmp(nfmax),       ! Temperature
     $                        frmuq(nfmax),       ! mu_B
     $                        frmus(nfmax),       ! mu_S
     $                        frnum(nfmax),       ! particle density
     $                        frbdn(nfmax)        ! baryon density



