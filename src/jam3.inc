c...common block for buu 

c...Number of spatial grid.
      integer maxx,maxz
c     parameter(maxx=20,maxz=24)
      parameter(maxx=1,maxz=1)

c...Force acting on hadrons.
      real*8 force,forcer,pots,potv,dtfree,tforce
      common/jamevf/force(3,mxv),forcer(3,mxv),
     $ pots(mxv),potv(0:3,mxv),dtfree(mxv),
     $ tforce

      integer mxw
      parameter (mxw = 1500)
      real*8 rho,rhoi,rhom,vmoms,vmom,rhoj,rhog
      common/jamevf2/rho(mxv),rhoi(mxv),rhog(mxv),rhom(mxw,mxw),
     $      vmoms(mxv),vmom(0:3,mxv),rhoj(0:3,mxv)

      ! for Coulomb potential
      real*8 rhoc,rhc
      common/jamevf3/rhoc(mxv),rhc(mxw,mxw)

      ! for symmetry energy
      real*8 rhos,rhs
      common/jamevf4/rhos(mxv),rhs(mxw,mxw)

      ! for Yukawa potential
      real*8 rhoy,rhy
      common/jamevf5/rhoy(mxv),rhy(mxw,mxw)

c...Baryon density.
      real*8 rhob
      common/jambuu1/rhob(-maxx:maxx,-maxx:maxx,-maxz:maxz)

c...Charge density.
      real*8 rhoz
      common/jambuu2/rhoz(-maxx:maxx,-maxx:maxx,-maxz:maxz)

c....Coulomb potential.
      real*8 cop
      common/jambuu3/cop(-maxx:maxx,-maxx:maxx,-maxz:maxz)

c...4-momentum of the grid.
      real*8 avp
      common/jambuu4/avp(4,-maxx:maxx,-maxx:maxx,-maxz:maxz)

c....end of jam3.inc

