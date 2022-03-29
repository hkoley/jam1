      implicit none
      include 'jam2.inc'
      integer i,j,k
      real*8 mult,a,b,getmult2
      integer ntmp,nmub,nmus
      parameter(ntmp=170,nmub=60,nmus=40)
      real*8 ptable,tmin,mubmin,musmin,dtmp,dmub,dmus
      common/multtable/ptable(0:ntmp,0:nmub,0:nmus),tmin,mubmin,musmin,
     & dtmp,dmub,dmus
      real*8 mub,mus,tf,tmax,mubmax,musmax
      real*8 a3,a4,a1,a2,tx,x

c     tf=0.02
c     mub=0.96
c     mus=0.0
c     print *,mult(tf,mub,mus)
c     print *,mult(tf,1.0,mus)
c     read(5,*)

c...Boltzmann approximation except pions. 
      mstc(146)=1

      call make_mult_table2

      tmax=tmin+dtmp*ntmp
      mubmax=dmub*nmub 
      musmax=dmus*nmus

      do i=1,100
        tf=tmin+(tmax-tmin)*ran(0)
        mub=(mubmax-mubmin)*ran(0)
        mus=(musmax-musmin)*ran(0)
        if((mub+mus)/tf.le.50d0) then
          a=getmult2(tf,mub,mus)
          b=mult(tf,mub,mus)
          write(6,800)tf,mub,mus
          write(6,810)a,b,abs(a-b)/b*100
        endif
      end do

        mub=0.0d0
        mus=0.0d0
        tf=tmin+dtmp*10
        a1=mult(tf,mub,mus)
        tf=tmin+dtmp*11
        a2=mult(tf,mub,mus)
        tf=tmin+dtmp*(10+ran(0))
        a3=mult(tf,mub,mus)
        x=(tf-(tmin+10*dtmp))/dtmp
        a4=x*a2+(1.0-x)*a1
        print *,a3,a4,abs(a3-a4)/a3*100

 800  format('T=',f12.7,' muB=',f12.7,' muS=',f12.7)
 810  format(e16.8,1x,e16.8,2x,f12.7)

      end

c**************************************************************************
      subroutine make_mult_table2
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

      write(10,800)ntmp,nmub,nmus
      write(10,810)tmin,dtmp
      write(10,810)mubmin,dmub
      write(10,810)musmin,dmus

      do j=0,nmub
        mub=mubmin+j*dmub
        print *,j,mub
      do k=0,nmus
        mus=musmin+k*dmus
      do i=0,ntmp
         tf=tmin+i*dtmp
        ptable(i,j,k)=0d0
        if((mub+mus)/tf.lt.50.0d0) then
          ptable(i,j,k)=mult(tf,mub,mus)
        endif
        if(ptable(i,j,k).lt.0d0) then
          print *,'ptable<0?',ptable(i,j,k),tf,mub,mus
          stop
        endif
        write(10,820)tf,mub,mus,ptable(i,j,k)
      end do
      end do
      end do

 800  format(3(i4,1x))
 810  format(3(e15.8,1x))
 820  format(4(e15.8,1x))

      end

c**************************************************************************
      real*8 function getmult2(t,mub,mus)
      implicit none
      integer ntmp,nmub,nmus
      parameter(ntmp=170,nmub=60,nmus=40)
      real*8 ptable,tmin,mubmin,musmin,dtmp,dmub,dmus
      common/multtable/ptable(0:ntmp,0:nmub,0:nmus),tmin,mubmin,musmin,
     & dtmp,dmub,dmus
      real*8 t,mub,mus,x1,y1,z1,x2,y2,z2
      integer i,j,k

      i=int((t-tmin)/dtmp)
      j=int((mub-mubmin)/dmub)
      k=int((mus-musmin)/dmus)

      if(i.lt.0.or.i.ge.ntmp) then
       print *,'T out of range',t,i
       stop
      endif
      if(j.lt.0.or.j.ge.nmub) then
       print *,'mu_B out of range',mub,j
       stop
      endif
      if(k.lt.0.or.k.ge.nmus) then
       print *,'mu_S out of range',mus,k
       stop
      endif

      x2=(t-(tmin+i*dtmp))/dtmp
      y2=(mub-(mubmin+j*dmub))/dmub
      z2=(mus-(musmin+k*dmus))/dmus
      x1=1d0-x2
      y1=1d0-y2
      z1=1d0-z2

      getmult2=x1*y1*z1*ptable(i,j,k)
     &       +x2*y1*z1*ptable(i+1,j,k)
     &       +x1*y2*z1*ptable(i,j+1,k)
     &       +x1*y1*z2*ptable(i,j,k+1)
     &       +x2*y1*z2*ptable(i+1,j,k+1)
     &       +x1*y2*z2*ptable(i,j+1,k+1)
     &       +x2*y2*z1*ptable(i+1,j+1,k)
     &       +x2*y2*z2*ptable(i+1,j+1,k+1)

      if(getmult2.eq.0d0) then
        print *,i,j,k
      endif
      end

