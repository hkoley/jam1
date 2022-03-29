c     include 'jam1.inc'
c     include 'jam2.inc'
      implicit double precision(a-h, o-z)

      character chau*16
      character frame*8,proj*8,targ*8,cwin*15
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2d0*a)

      common/jyjets/n,npad,k(1000,5),p(1000,5),v(1000,5)
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      dimension npa(3)

c....Initialize JAM.
      nevent=10000

      b1=3.0
      b2=3.0
      frame='collider'        ! comp. frame
      dt=0.3d0          ! collision time(fm/c)
      nstp=100            ! time step (i.e. no time step)
      cwin='4.93gev'      ! initial c.m. energy per nucl, in this case,
      proj='p'       ! projectile
      targ='p'       ! target
      call jaminit(nev,b1,b2,dt,nstp,frame,proj,targ,cwin)

c     nevent=mstc(2)
      mstu(11)=6
      mstj(21)=1
        mstu(10)=2

      npa(1)=0
      npa(2)=0
      npa(3)=0
c     ecm=4.93
      ecm=2.0
      kf1=2
      kf2=2101
c==================================================
      do iev=1,nevent

        call pj2ent(0,kf1,kf2,ecm)
c       call pjlist(1)
c       read(5,*)

         do i=1,n
         if(k(i,1).ge. 11 .or.k(i,1).eq.0) cycle
         if(k(i,2).eq. 2212) npa(1)=npa(1)+1
         if(k(i,2).eq. -211) npa(2)=npa(2)+1
         if(k(i,2).eq.  211) npa(3)=npa(3)+1
        end do


c....End simulation
      end do

      wei=1.0/nevent
      print *,'proton= ',npa(1)*wei
      print *,'pi- ',npa(2)*wei
      print *,'pi+= ',npa(3)*wei

      call jamfin

      end
