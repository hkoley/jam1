c...parameters and common block for RQMD/S
c     integer mx
c     real*8 rho0,pi,hc
c     parameter (mx = 1500)
c     parameter (rho0 = 0.168d0)
c     parameter (pi=3.141592653589793d0)
c     parameter (hc=0.19732705d0)

      real*8 t1,t2,t3,t3f,t4,t4f,pmu1,pmu2,vex1,vex2
      common /potpar/t1,t2,t3,t3f,t4,t4f,pmu1,pmu2,vex1,vex2
      save /potpar/

