      integer maxe,maxn
      real*8 emin,emax,bmax,de,dn
      real*8 eostable,eostabt,eostabmu,eostabsmu,eostablam,eostabs
      real*8 eostabp

c     parameter(emin=0.0003,de=0.05,dn=0.05,maxe=200, maxn=40)
c     parameter(emin=0.00,emax=10d0,maxe=500,bmax=3d0,maxn=150)
c2018/10/18
      parameter(emin=0.00,emax=20d0,maxe=1000,bmax=3d0,maxn=150)

      common/jameostable/eostable(0:maxn,0:maxe),
     & eostabt(0:maxn,0:maxe),eostabmu(0:maxn,0:maxe),
     & eostabsmu(0:maxn,0:maxe),eostablam(0:maxn,0:maxe),
     & eostabs(0:maxn,0:maxe),
     & eostabp(0:maxn,0:maxe),
     &      de,dn
      save /jameostable/

      real*8 Bagconst
      common/jameostable2/Bagconst
      save/jameostable2/

c     real*8 pnow,tnow,bmnow,smnow,fnow,snow,csnow
c     common/jameos1/pnow,tnow,bmnow,smnow,fnow,snow


