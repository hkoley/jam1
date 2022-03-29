
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      include 'jameos.inc'
      include 'fluid.inc'
c     common/jameostable2/ttable(0:maxn,0:maxe),
c    & stable(0:maxn,0:maxe),
c    & bmutable(0:maxn,0:maxe)

       character*80 dir

      iprint=1

      eos_mode=11
c     fname(9)='/export/ynara/lib/eos_MF_fullB220.dat'
c     fname(9)='/export/ynara/lib/eos_MF_HG.dat'
c     fname(9)='/export/ynara/lib/eos_MF_HGK045.dat'
c     fname(9)='/export/ynara/lib/eos_MF_freeHG.dat'
      fname(9)='/export/ynara/lib/eos_MF_HGK045JAM.dat'
      call readEOStable(fname(9))
c     call readEOStable(fname(9))
      dir='/export/ynara/lib/'
      iopt=2 ! HG
c     iopt=3 ! 1st order
c     iopt=4 ! cross over
c     iopt=5
      call readeos(iopt,dir)

c     rho0 = 0.15891d0
      rho0 = 0.168
      b=0*rho0
c     b=1*rho0

      srho=10.0
c     write(6,*)'# input s/rho'
c     read(5,*)srho

      eemax=2.5 
      dee=0.05
      nne=eemax/dee

      db=0.001
      nnb=bmax/db

c....factor needed for urqmd eos.
      fac = 146.51751415742*1d-3
      b0 = 0.15891d0 

      if(iprint.eq.1) then

      do i=1,nne
        e=i*dee+0.03
        eos_mode=11
        call geteos3(e/hbc,b,p1,cs)
        eos_mode=1
        call geteos3(e/hbc,b,p2,cs)
        print *,e,b,p1*hbc,p2*hbc
c       print *,e,p1,press(e/fac,b/b0)*fac
      end do

      return
      endif


      do i=1,nne
        e=i*dee+0.1

        do j=1,nnb
          b=j*db
          s=b*srho
c         s0=getEoStab(b,e,stable)
          call geteos2(e,b,t,bmu,smu,s0)
          if(abs(s-s0).le.0.05) then
            s1=s0
            b1=b
            goto 100
          endif
c           print *,j,' e=',e,'b=',b,'s=',s,'s1=',s1
        end do
        print *,'no solution?',j,b
        stop
100     continue

        do j=1,nnb
          b=j*db
          s=b*srho
          s1=sress(e/fac,b/b0)*b0
c         s1=entro(e/fac,b/b0)*b0
          if(abs(s-s1).le.0.05) then
c           print *,'e=',e,'b=',b,'s=',s,'s1=',s1
c           read(5,*)
            s2=s1
            b2=b
            goto 200
          endif
        end do
        print *,'(2) no solution?',b
        stop
200     continue


        p1=getP(b1,e)
        p2=press(e/fac,b2/b0)*fac
        write(6,800)e,p1,p2,s1,s2,b1,b2
      end do
800   format(8(1x,e12.5))

      end


c----------------------------------------------------------------------
c
      function sress(e,n)
c     INCLUDE 'defs.f'
      integer ngr
       PARAMETER (ngr=200)

c
c  This function-subprogram determines the pressure of the
c  underlying EoS.
c  It is used in subroutines sinit, tinit, untang, velo.
c
c  Type declarations for variables in common-blocks.
c
      real*8 e0,n0,B,hc3,pi2
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)


      integer eos,stabil,anti
c
c  Type declarations for variables used in function press.
c
      real*8 e,n,sress,et0
      real*8 de,dn,p1,p2,p3,p4,p13,p24
c      real*8 p12,p34
c
c  Common-blocks.
c
      common /grstate/ e0,n0,B,hc3,pi2
      common /eqofst/ eos,stabil,anti
      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3


c
c  Fermi energy density in the QGP for given baryon density
c  (in units of e0).
c
      et0 = 1d0/54d0/pi2*dabs(40.5d0*pi2*n*n0*hc3)**(4d0/3d0)+B
      et0 = et0/e0/hc3

      anti = 0


      if (n.lt.0.0)then
         n=-n
         anti = 1
      end if

c
c     Calculate pressure (in units of e0).
c     If flag eos=[else], use tabellized equation of state.
c     
c     
c     arithmetic mean interpolation (yields better values for
c     interpolation at small e,n, but is worse for all other
c     thermodynamic quantities, also at larger e,n)
c     
c     If e and n are both smaller than e0, n0, the nuclei are in the
c     ground state. The pressure is set to zero by hand, if flag
c     stabil = 1.
c     
      if (eos.eq.3) then
         if (e.le.20d0) then
            de = 0.1d0
            dn = 0.05d0
            p1 = stab(idint(e/de),idint(n/dn))
            p2 = stab(idint(e/de)+1,idint(n/dn))
            p3 = stab(idint(e/de),idint(n/dn)+1)
            p4 = stab(idint(e/de)+1,idint(n/dn)+1)
            p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
            p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
            p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
            p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
      
            sress = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     sress = 0.25d0*(p13+p24+p12+p34)
         else if (e.ge.20d0) then
            sress = (e-4d0*B/e0/hc3)/3d0
         end if
         
      else if ((eos.eq.4).or.(eos.eq.2).or.(eos.eq.5)) then
         if (e.le.1000.0d0) then

            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = stab3(idint(e/de),idint(n/dn))
               p2 = stab3(idint(e/de)+1,idint(n/dn))
               p3 = stab3(idint(e/de),idint(n/dn)+1)
               p4 = stab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2 
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               sress = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     write(99,*) p1, stab3(0,0) 
            end if
            
            if (((e.lt.10.0d0).and.(n.lt.2.0d0)).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = stab2(idint(e/de),idint(n/dn))
               p2 = stab2(idint(e/de)+1,idint(n/dn))
               p3 = stab2(idint(e/de),idint(n/dn)+1)
               p4 = stab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2 
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               sress = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then

               de = 0.5d0
               dn = 0.1d0
 
               p1 = stab(idint(e/de),idint(n/dn))
               p2 = stab(idint(e/de)+1,idint(n/dn))
               p3 = stab(idint(e/de),idint(n/dn)+1)
               p4 = stab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2 
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               sress = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13  
            end if
         else if (e.ge.1000.0d0) then
c     JS muss bestimmt werden
            sress = e/3.0d0
         end if
      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1).and.(n.ge.-1d0))
     $     then
         sress = 0d0
      end if
      
      if (sress.lt.0d0) sress=0d0
      
      
      
      
      
      if (anti.eq.1)then
         n=-n
      end if
      return
      end

