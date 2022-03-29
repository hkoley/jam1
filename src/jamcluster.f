! YN  last modified 2022.02.13 nuclear cluster formation added
! YN  last modified 2022.03.23 nuclear cluster formation revised
c***********************************************************************
      subroutine jamdistance1(i1,i2,rr,pp)
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      integer i1,i2
      real*8 r1(5),r2(5),dx(5),dp(5),pcm(5),s,rp,p1sq,p2sq,rr,pp,tc

c...Set global time for clusterlization.
      tc=max(pare(1),r(4,i1),r(4,i2))
      r1(:)=r(:,i1)+(tc-r(4,i1))*p(:,i1)/p(4,i1)
      r2(:)=r(:,i2)+(tc-r(4,i2))*p(:,i2)/p(4,i2)

c...clustering time in the computational frame.
      r(5,i1)=tc
      r(5,i2)=tc

c...compute relative distance and momentum in their two-body c.m.
      dx(:) = r1(:) - r2(:)
      pcm(:)= p(:,i1) + p(:,i2)
      s=pcm(4)**2 - pcm(1)**2 - pcm(2)**2 - pcm(3)**2
      rp= -dx(4)*pcm(4) + dx(1)*pcm(1) + dx(2)*pcm(2) + dx(3)*pcm(3)
      rr = -dx(4)**2 + dx(1)**2 + dx(2)**2 + dx(3)**2 + rp**2/s

      dp(:) = p(:,i1) - p(:,i2)
      p1sq=p(4,i1)**2 - p(1,i1)**2 - p(2,i1)**2 - p(3,i1)**2
      p2sq=p(4,i2)**2 - p(1,i2)**2 - p(2,i2)**2 - p(3,i2)**2
      pp = -dp(4)**2 + dp(1)**2 + dp(2)**2 + dp(3)**2 + (p1sq-p2sq)**2/s

      end

c***********************************************************************
      subroutine jamdistance2(i1,i2,rr,pp)
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      integer i1,i2,j
      real*8 pcm1,pcm2,pcm3,pcm4,srt,bex,bey,bez,gam
      real*8 p1(5),p2(5),r1(5),r2(5),t,rr,pp

      pcm1=p(1,i1)+p(1,i2)
      pcm2=p(2,i1)+p(2,i2)
      pcm3=p(3,i1)+p(3,i2)
      pcm4=p(4,i1)+p(4,i2)
      srt=sqrt(pcm4**2-pcm1**2-pcm2**2-pcm3**2)
      bex=pcm1/pcm4
      bey=pcm2/pcm4
      bez=pcm3/pcm4
      gam=pcm4/srt
      p1(:)=p(:,i1)
      p2(:)=p(:,i2)
      r1(:)=r(:,i1)
      r2(:)=r(:,i2)
      call jamrobo(0d0,0d0,-bex,-bey,-bez,gam,p1(1),p1(2),p1(3),p1(4))
      call jamrobo(0d0,0d0,-bex,-bey,-bez,gam,p2(1),p2(2),p2(3),p2(4))
      call jamrobo(0d0,0d0,-bex,-bey,-bez,gam,r1(1),r1(2),r1(3),r1(4))
      call jamrobo(0d0,0d0,-bex,-bey,-bez,gam,r2(1),r2(2),r2(3),r2(4))
      t=max(r1(4),r2(4))
      r1(4)=t
      r2(4)=t
      do j=1,3
       r1(j)=r1(j)+(t-r1(4))*p1(j)/p1(4)
       r2(j)=r2(j)+(t-r2(4))*p2(j)/p2(4)
      end do
      rr=(r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2
      pp=(p2(1)-p1(1))**2+(p2(2)-p1(2))**2+(p2(3)-p1(3))**2

c...boost back to the original frame.
      call jamrobo(0d0,0d0,bex,bey,bez,gam,r1(1),r1(2),r1(3),r1(4))
      call jamrobo(0d0,0d0,bex,bey,bez,gam,r2(1),r2(2),r2(3),r2(4))

c...clustering time in the computational frame.
      r(5,i1)=r1(4)
      r(5,i2)=r2(4)

      end

c***********************************************************************
      subroutine jamdistance3(i1,i2,bsq,pp)
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      integer i1,i2
      real*8 bsq,pp,dt,dx,dy,dz,rsqare,dx12,dxp1,dxp2,dp12,dn,dt1,dt2,
     &  pc1(5),pc2(5),pcm(5),dr(5),s,rp,em1sq,em2sq,t01,t02,b12,
     &  tc,tcol1,tcol2

      dx=r(1,i2)-r(1,i1)
      dy=r(2,i2)-r(2,i1)
      dz=r(3,i2)-r(3,i1)
      rsqare=dx**2+dy**2+dz**2

      pc1(:)=p(:,i1)
      pc2(:)=p(:,i2)
      em1sq=pc1(4)**2-pc1(1)**2-pc1(2)**2-pc1(3)**2
      em2sq=pc2(4)**2-pc2(1)**2-pc2(2)**2-pc2(3)**2

      pcm(:)= pc1(:) + pc2(:)
      s=pcm(4)**2 - pcm(1)**2 - pcm(2)**2 - pcm(3)**2

      t01=r(4,i1)
      t02=r(4,i2)
      dt=t02-t01
      dx12=dt**2-rsqare
      dxp1=dt*pc1(4)-dx*pc1(1)-dy*pc1(2)-dz*pc1(3)
      dxp2=dt*pc2(4)-dx*pc2(1)-dy*pc2(2)-dz*pc2(3)
      dp12=pc1(4)*pc2(4)-pc1(1)*pc2(1)-pc1(2)*pc2(2)-pc1(3)*pc2(3)
      dn=dp12*dp12-em1sq*em2sq
      if(dn.ge.1d-9) then
        dt1=-pc1(4)*(dxp1*em2sq-dxp2*dp12)/dn
        dt2= pc2(4)*(dxp2*em1sq-dxp1*dp12)/dn
        if(dt1>0.0.and.dt2>0.0) then
          b12=dxp1**2*em2sq+dxp2**2*em1sq-2.d0*dxp1*dxp2*dp12
          bsq=-dx12-b12/dn
          tcol1=t01+dt1
          tcol2=t02+dt2
        else
          bsq=-dx12+(dxp1+dxp2)**2/s
          tcol1=t01
          tcol2=t02
        endif
      else
        tc=max(t01,t02)
        tcol1=tc
        tcol2=tc
        dr(:) = r(:,i1)+(tc-t01)*p(:,i1)/p(4,i1)
     &        - r(:,i2)+(tc-t02)*p(:,i2)/p(4,i2)
        rp= -dr(4)*pcm(4) + dr(1)*pcm(1) + dr(2)*pcm(2) + dr(3)*pcm(3)
        bsq = -dr(4)**2 + dr(1)**2 + dr(2)**2 + dr(3)**2 + rp**2/s
      endif

c...relative momentum in the two-body c.m.
      pp = -(pc1(4)-pc2(4))**2+(pc1(1)-pc2(1))**2+(pc1(2)-pc2(2))**2
     &  +(pc1(3)-pc2(3))**2+(em1sq-em2sq)**2/s

c...clustering time in the computational frame.
      r(5,i1)=tcol1
      r(5,i2)=tcol2

      end

************************************************************************
      logical function clust(i1,i2)
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      integer i1,i2
      real*8 rr,pp

      clust=.false.
      if(k(9,i1)*k(9,i2)<0) return ! B-antiB

      if(mstc(131).eq.1) then
        call jamdistance1(i1,i2,rr,pp)
      else if(mstc(131).eq.2) then
        call jamdistance2(i1,i2,rr,pp)
      else
        call jamdistance3(i1,i2,rr,pp)
      endif

      if(rr.le.parc(131)**2.and.pp.le.parc(132)**2) then 
        clust=.true.
      endif

      end
************************************************************************
      subroutine jamclusterform
      implicit none
      include 'jam1.inc'
      include 'jam2.inc'
      integer,allocatable :: mscl(:),num(:)
      integer nclust,icheck,i1,i2,i,j,lp,iz,inn,is,il,ix,io,ii,nt
      integer jamcomp,ic,id,kf,kc,jj,j0,istr
      logical unknown,clust
      real*8 rc(5),pc(5),lc(3)
      real*8 tclust
      real*8 ptot(5),ptot1(5)

      ptot(:)=0.0
      do i=1,nv
       ptot(:)=ptot(:)+p(:,i)
      end do

c...better to use nbary instead of nv?
      allocate( mscl(nv) )
      allocate( num(nv) )
      mscl(:)=1
      do i=1,nv
        num(i)=i
      end do
      nclust=1
      icheck=1


c...Find nuclear clusters
      do i=1,nv-1
        i1=num(i)
        if(k(1,i1).gt.10) cycle     ! dead particle
        if(k(9,i1).eq.0) cycle      ! meson
        j0=icheck+1
        do j=j0,nv
          i2=num(j)
          if(k(1,i2).gt.10) cycle    ! dead particle
          if(k(9,i2).eq.0) cycle     ! meson
          if(clust(i1,i2)) then
            lp=num(icheck+1)
            num(icheck+1)=i2
            num(j)=lp
            icheck=icheck+1
            mscl(nclust)=mscl(nclust)+1
          endif
        end do
        if(icheck == i) then
          nclust = nclust + 1
          icheck = icheck + 1
        endif
      end do

c...Loop over clusters.
      ii=1
      do ic=1,nclust
        nt=mscl(ic)
        jj=ii
        if(nt.eq.1) cycle
        if(nt.eq.2) then
          i1=num(ii)
          i2=num(ii+1)
          if(k(2,i1)==2112 .and. k(2,i2)==2112) cycle ! (n,n)
          if(k(2,i1)==2212 .and. k(2,i2)==2212) cycle ! (p,p)
        endif
        iz=0
        inn=0
        il=0
        is=0
        ix=0
        io=0
        istr=0
        unknown=.false.
c....loop over nucleons in the cluster.
        do i=1,nt
          j=num(ii)
          ii=ii+1
          kf=k(2,j)
          kc=jamcomp(kf)
          id=kchg(kc,5)
          istr=istr + abs(kchg(kc,7))
          if(abs(kf)==2112) then
            inn=inn+1
          else if(abs(kf)==2212) then
            iz=iz+1
          else if(id==id_lamb .or.id==id_lambs) then
            il=il+1
          else if(id==id_sigm .or.id==id_sigms) then
            is=is+1
          else if(id==id_xi .or.id==id_xis) then
            ix=ix+1
          else if(id==id_omega) then
            io=io+1
c....In case you want to investigate exotic cluster,
c....this part should be modified.
          else 
            unknown=.true.
c           print *,'jamcluster unkown id=',kf
          endif
        end do
        if(unknown) cycle
        if(nt>=3.and.iz==0) cycle
        if(nt>=3.and.inn==0) cycle

        rc(:)=0.0
        pc(:)=0.0
        lc(:)=0.0
        nv=nv+1
        kf=k(2,num(jj))

c...Find global time for clusterlization.
        tclust=-1.0
        do i=jj,nt+jj-1
          j=num(i)
          tclust=max(tclust,r(5,j))
        end do

        do i=1,nt
          j=num(jj)
          jj=jj+1
          k(1,j)=30
          rc(:)=rc(:)+r(:,j)+(tclust-r(4,j))*p(:,j)/p(4,j)
          pc(:)=pc(:)+p(:,j)
          lc(1)=lc(1)+r(2,j)*p(3,j) - r(3,j)*p(2,j)
          lc(2)=lc(2)+r(3,j)*p(1,j) - r(1,j)*p(3,j)
          lc(3)=lc(3)+r(1,j)*p(2,j) - r(2,j)*p(1,j)
        end do
       mstd(30)=mstd(30)+1
       k(:,nv)=0
       k(1,nv)=1
       kf=1 ! for now 
       k(2,nv)=(1000000000 + nt*10 + iz*10000
     &  + istr*10000000 )*isign(1,kf)
       k(3,nv)=il
       k(4,nv)=is
       k(5,nv)=ix
       k(6,nv)=io
       k(7,nv)=nt
       r(:,nv)=rc(:)/nt
       p(:,nv)=pc(:)
       p(5,nv)=sqrt(p(4,nv)**2-p(1,nv)**2-p(2,nv)**2-p(3,nv)**2)
       v(:,nv)=p(:,nv)
       kq(:,nv)=0
       vq(1,nv)=lc(1)
       vq(2,nv)=lc(2)
       vq(3,nv)=lc(3)
       vq(:,nv)=0.0
      end do ! end loop over clusters

      deallocate(mscl)
      deallocate(num)


      if(mstc(132).eq.1) call jamedit

      ptot1(:)=0.0
      do i=1,nv
       if(k(1,i).gt.10) cycle
       ptot1(:)=ptot1(:)+p(:,i)
      end do
      if(abs((ptot(1)-ptot1(1))**2+(ptot(2)-ptot1(2))**2
     & +(ptot(3)-ptot1(3))**2).gt.1e-8) then
        print *,'jamcluster momtnum does not converve'
        print *,'ptot0=',ptot(:)
        print *,'ptot1=',ptot1(:)
        stop
      endif

      end

