
c...written by Yukinao Akamatsu 2017
c...modified by Y. Nara 2017-2018
      subroutine RHS(mx,d,vx,vy,vz,e,c,lam,fi,opt_timelike)
      implicit none
      integer mx,i,opt_timelike
      double precision d(0:mx),vx(0:mx),vy(0:mx),vz(0:mx),e(0:mx)
c     double precision u(0:4,0:mx),g(0:4,0:mx)
      double precision ql(0:4,0:mx),qr(0:4,0:mx),c(0:mx)
      double precision grdQ(0:4,0:mx),fi(0:4,0:mx)
      double precision depend, lam, facdep
      parameter(facdep=0.5d0)
      real*8 u1(0:4),u2

      call Limiting(mx, e,grdQ(0,:))
      call Limiting(mx,vx,grdQ(1,:))
      call Limiting(mx,vy,grdQ(2,:))
      call Limiting(mx,vz,grdQ(3,:))
      call Limiting(mx, d,grdQ(4,:))

      call QLIMIT(mx,d,vx,vy,vz,e,grdQ)

      do I = 1, mx-1
        depend = (1.d0 - c(i)*lam) * facdep   ! Domain of dependence
        QL(0,I) =  e(I) - depend * grdQ(0,I)
        QL(1,I) = vx(I) - depend * grdQ(1,I)
        QL(2,I) = vy(I) - depend * grdQ(2,I)
        QL(3,I) = vz(I) - depend * grdQ(3,I)
        QL(4,I) =  d(I) - depend * grdQ(4,I)

        QR(0,I) =  e(I) + depend * grdQ(0,I)
        QR(1,I) = vx(I) + depend * grdQ(1,I)
        QR(2,I) = vy(I) + depend * grdQ(2,I)
        QR(3,I) = vz(I) + depend * grdQ(3,I)
        QR(4,I) =  d(I) + depend * grdQ(4,I)
      end do

      do I = 1, mx-2
        call hlle(QL(:,I+1), QR(:,I), FI(:,I),lam)
      end do

c...K. Murase method to preserve timelike nature of energy-momentum
c...tensor.
c     if(opt_timelike.ge.1) call limiter_timelike(mx,fi,u,vx,lam)

c     do I = 2, mx-2
c       g(:,I) = lam * (FI(:,I-1) - FI(:,I))
c       g(:,I) = FI(:,I-1) - FI(:,I)
c     end do

c...Boundary condition I.
c     g(:,0)=g(:,2)
c     g(:,1)=g(:,2)
c     g(:,mx-1)=g(:,mx-2)
c     g(:,mx)=g(:,mx-2)

c...Boundary condition II.
      fi(:,0)=fi(:,1)
      fi(:,mx)=0d0
      fi(:,mx-1)=fi(:,mx-2)

c     g(:,1)=0d0
c     g(:,0)= -fi(:,1)
c     g(:,mx) = fi(:,mx-2)
c     g(:,mx-1)=0d0

c...Boundary condition III.
c     g(:,0)=0d0
c     g(:,1)= -fi(:,1)
c     g(:,mx-1) = fi(:,mx-2)
c     g(:,mx)=0d0


c     u1(:)=u(:,1)+g(:,1)
c     u2=sqrt(u1(1)**2+u1(2)**2+u1(3)**2)
c     if(u1(0).lt.u2) g(0,1)=u2-u(0,1)+1d-12 
c     if(u1(0).lt.u2) g(:,1)=0d0

c     u1(:)=u(:,mx-1)+g(:,mx-1)
c     u2=sqrt(u1(1)**2+u1(2)**2+u1(3)**2)
c     if(u1(0).lt.u2) g(0,mx-1)=u2-u(0,mx-1)+1d-12 
c     if(u1(0).lt.u2) g(0,mx-1)=0d0

c     do I = 2, mx-2
c       U(:,I) = U(:,I) + lam * (FI(:,I-1) - FI(:,I))
c       if(u(0,i)**2-u(1,i)**2-u(2,i)**2-u(3,i)**2 < 0d0) u(:,i)=0d0
c     end Do
c     u(:,0)=u(:,2)
c     u(:,1)=u(:,2)
c     u(:,mx-1)=u(:,mx-2)
c     u(:,mx)=u(:,mx-2)

  
      end


      subroutine hlle(QR,QL,FI,lam)
      implicit none
      integer  I
! left & right
! Q = e, vx, vy, vz, d; F = numerical flux for E, Mx, My, Mz, D
      double precision QL(0:4), QR(0:4), FI(0:4)
! thermodynamic variables
      double precision pL, c2L, pR, c2R,lam
      double precision WL, WR                     ! gamma factors
      double precision UL(0:4), UR(0:4), FL(0:4), FR(0:4)
      double precision bL, bR                     ! signal velocities
      double precision eLrt, eRrt                 ! sqrt of eL and eR
      double precision vxI, WI, c2I,q

      FI(:)=0d0
      UL(:)=0d0
      FL(:)=0d0
      UR(:)=0d0
      FR(:)=0d0
      c2L=0d0
      c2R=0d0

      q = QL(1)**2 + QL(2)**2 + QL(3)**2
      if(q.lt.1d0) then
        WL = 1.d0/sqrt(1.d0 - q)
        call geteos3(QL(0),QL(4),pL,c2L)
        UL(0) = (pL + QL(0))*WL**2 - pL     ! E
        UL(1) = (pL + QL(0))*WL**2*QL(1)    ! Mx
        UL(2) = (pL + QL(0))*WL**2*QL(2)    ! My
        UL(3) = (pL + QL(0))*WL**2*QL(3)    ! Mz
        UL(4) = WL*QL(4)                    ! D

        FL(0) = UL(1)                ! flux for E
        FL(1) = UL(1)*QL(1) + pL     ! flux for Mx       
        FL(2) = UL(2)*QL(1)          ! flux for My
        FL(3) = UL(3)*QL(1)          ! flux for Mz
        FL(4) = UL(4)*QL(1)          ! flux for D
        q = ul(0)**2-ul(1)**2-ul(2)**2-ul(3)**2
        if(ul(0).le.0d0 .or. q .le.0d0) then
          ul(:)=0d0; fl(:)=0d0; c2L=0d0
        endif
      endif

      q = QR(1)**2 + QR(2)**2 + QR(3)**2
      if(q.lt.1d0) then
        WR = 1.d0/sqrt(1.d0 - q)
        call geteos3(QR(0),QR(4),pR,c2R)
        UR(0) = (pR + QR(0))*WR**2 - pR     ! E
        UR(1) = (pR + QR(0))*WR**2 *QR(1)   ! Mx
        UR(2) = (pR + QR(0))*WR**2 *QR(2)   ! My
        UR(3) = (pR + QR(0))*WR**2 *QR(3)   ! Mz
        UR(4) = WR*QR(4)                    ! D

        FR(0) = UR(1)                ! flux for E
        FR(1) = UR(1)*QR(1) + pR     ! flux for Mx       
        FR(2) = UR(2)*QR(1)          ! flux for My
        FR(3) = UR(3)*QR(1)          ! flux for Mz
        FR(4) = UR(4)*QR(1)          ! flux for D
        q = ur(0)**2-ur(1)**2-ur(2)**2-ur(3)**2
        if(ur(0).le.0d0 .or.q.le.0d0) then
          ur(:)=0d0; fr(:)=0d0; c2R=0d0
        endif
      endif

      eLrt = sqrt(UL(0));   eRrt = sqrt(UR(0))
      if(elrt+errt.eq.0d0) return

      !! 0.5d0 is the suggested value
      c2I = (eLrt*c2L*c2L + eRrt*c2R*c2R)/(eLrt + eRrt)
     & + 0.5d0*(eLrt*eRrt)*(QR(1)-QL(1))**2/(eLrt + eRrt)**2
      c2I=sqrt(c2I)

      vxI = (eLrt*QL(1) + eRrt*QR(1))/(eLrt + eRrt)

      bL = min(0.d0, (vxI - c2I)/(1.d0 - vxI*c2I), 
     &     (QL(1) - c2L)/(1.d0 - QL(1)*c2L))


      bR = max(0.d0, (vxI + c2I)/(1.d0 + vxI*c2I), 
     &     (QR(1) + c2R)/(1.d0 + QR(1)*c2R))

      if(bR-bL.ne.0.0) then
       FI(:) = lam*(bR*FL(:) - bL*FR(:) + bL*bR*(UR(:)-UL(:)))/(bR - bL)
      endif

      end
!======================================================================*
!                      Van Leer's Limiter
!======================================================================*
      function VLlimiter(a,b,c)
      implicit none
      double precision VLlimiter,a,b,c
      double precision, parameter :: alpha = 2.d0

      VLlimiter = DSIGN(1.D0, c)*(DSIGN(0.5D0,a*b) + 0.5D0)
     &   *dmin1(alpha*DABS(a), alpha*DABS(b), DABS(c))

      end

!======================================================================*
!                      Compute grad
!======================================================================*
      subroutine Limiting(mxx,U,grdU)
      implicit none
      integer i,mxx
      double precision dul,dur,dum
      double precision U(0:mxx),grdu(0:mxx)
      double precision, parameter :: alpha = 2.d0
      double precision VLlimiter

      do i=1,mxx-1
        dul=u(i)   - u(i-1)
        dur=u(i+1) - u(i)
        dum=0.5*(u(i+1)-u(i-1))
c       grdU(i)=VLlimiter(dul,dur,dum)
        grdU(i)=sign(1d0,dum)*(sign(0.5d0,dul*dur) + 0.5d0)
     &     *min(alpha*abs(dul), alpha*abs(dur), abs(dum)) ! Van Leer
      end do

      end

!======================================================================*
!           L and R state of primitive variables
!======================================================================*
      subroutine QLIMIT(m,d,vx,vy,vz,e,grdQ)
      implicit none
      integer i,k,m
      double precision  d(0:m),vx(0:m),vy(0:m),vz(0:m),e(0:m)
      double precision  grdQ(0:4,0:m) 
      double precision  QL(0:4,0:m), QR(0:4,0:m)
      double precision :: vL, vR, vmax, vmin
 
      QL(0,1:m-1) =  e(1:m-1) - 0.5d0 * grdQ(0,1:m-1)
      QL(1,1:m-1) = vx(1:m-1) - 0.5d0 * grdQ(1,1:m-1)
      QL(2,1:m-1) = vy(1:m-1) - 0.5d0 * grdQ(2,1:m-1)
      QL(3,1:m-1) = vz(1:m-1) - 0.5d0 * grdQ(3,1:m-1)
      QL(4,1:m-1) =  d(1:m-1) - 0.5d0 * grdQ(4,1:m-1)

      QR(0,1:m-1) =  e(1:m-1) + 0.5d0 * grdQ(0,1:m-1)
      QR(1,1:m-1) = vx(1:m-1) + 0.5d0 * grdQ(1,1:m-1)
      QR(2,1:m-1) = vy(1:m-1) + 0.5d0 * grdQ(2,1:m-1)
      QR(3,1:m-1) = vz(1:m-1) + 0.5d0 * grdQ(3,1:m-1)
      QR(4,1:m-1) =  d(1:m-1) + 0.5d0 * grdQ(4,1:m-1)

      do i = 1, m-2

      if(QL(4,I+1) < 0.d0) then
        grdQ(4,I+1) = 0.d0
        grdQ(4,I  ) = 0.d0
      else if(QR(4,I) < 0.d0) then
        grdQ(4,I  ) = 0.d0
        grdQ(4,I+1) = 0.d0
      end if
      if(QL(0,I+1) < 0.d0) then
        grdQ(0,I+1) = 0.d0
        grdQ(0,I  ) = 0.d0
      else if(QR(0,I) < 0.d0) then
        grdQ(0,I  ) = 0.d0
        grdQ(0,I+1) = 0.d0
      end if
     
      vL = sqrt(QL(1,I+1)**2 + QL(2,I+1)**2 + QL(3,I+1)**2)
      vR = sqrt(QR(1,I  )**2 + QR(2,I  )**2 + QR(3,I  )**2)
     
      vmax = max(abs(vx(I)),abs(vx(I+1)), 
     &       abs(vy(I)),abs(vy(I+1)), abs(vz(I)),abs(vz(I+1)) )
      vmin = min(abs(vx(I)),abs(vx(I+1)),
     &       abs(vy(I)),abs(vy(I+1)), abs(vz(I)),abs(vz(I+1)) )
     
      if(vL > vmax) then
        grdQ(1:3,I+1) = 0.d0
        grdQ(1:3,I  ) = 0.d0
      else if(vR > vmax) then
        grdQ(1:3,I ) = 0.d0;
        grdQ(1:3,I+1) = 0.d0
      end if
     
      if(vL < vmin) then
        grdQ(1:3,I+1) = 0.d0
        grdQ(1:3,I  ) = 0.d0
      else if(vR < vmin) then
        grdQ(1:3,I  ) = 0.d0
        grdQ(1:3,I+1) = 0.d0
      end if

      end do
  
      end

