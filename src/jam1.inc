c***********************************************************************
c...Include file of JAM
c***********************************************************************

C...Double precision.
c     implicit double precision(a-h, o-z)

      integer mxv
c...Maximum number of all particles.
c     parameter(mxv=800)   ! AGS Au+Au
c     parameter(mxv=2000)
c     parameter(mxv=30000)
c     parameter(mxv=100000)
c     parameter(mxv=3000)
c     parameter(mxv=5000)
c     parameter(mxv=10000)
      parameter(mxv=30000)

c...The event record.
      real*8 r,p,v
      integer k,nv,nbary,nmeson
      common/jamevnt1/r(5,mxv),p(5,mxv),v(5,mxv),k(11,mxv)
      common/jamevnt2/nv,nbary,nmeson
      save /jamevnt1/,/jamevnt2/

c...MF status/Can Collide or Not
      integer MF_on,IcanColl
      common/jamRQMF/MF_on(mxv),IcanColl(mxv)
      save /jamRQMF/

c...Parton properties or decay orientation.
      real*8 vq
      integer kq
      common/jamjet1/vq(10,mxv),kq(2,mxv)
      save /jamjet1/

c...Deuteron coalescence.
      real*8 weid
      common/jamevnt3/weid(mxv)
      save /jamevnt3/

c...Save two-body collision information.
      real*8 rcp,pcp,vcp,vqcp,fcp,frcp,spcp,vpcp,dtfcp
      integer kcp,kqcp,MFcp
      common/jamsave1/rcp(5,6),pcp(5,6),vcp(5,6),kcp(11,6)
      common/jamsave2/vqcp(10,6),kqcp(2,6)
      common/jamsave3/fcp(3,6),frcp(3,6),spcp(6),vpcp(0:3,6),
     &       dtfcp(6),MFcp(6)
      save /jamsave1/,/jamsave2/

c...Max. collision number in each time step.
      integer mxcoll
c     parameter (mxcoll=50000)
      parameter (mxcoll=500000)
c     parameter (mxcoll=1000)

c...Collision matrix.
      real*8 coll
      integer icoll,mentry
      common/jamcol/coll(6,mxcoll),icoll(3,mxcoll),mentry
      save /jamcol/

c...Box by maru  cell.
      real*8 rcell
      common/jamcell/rcell(3,27)
      save /jamcell/

c=======================================================================
