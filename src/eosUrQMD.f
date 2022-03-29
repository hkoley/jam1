c     real*8 p,e,n,press,de
c     integer i,nn,eos
c     real*8 e0,n0,Bag,hc3,pi2
c     real*8 fac
c     common /grstate/ e0,n0,Bag,hc3,pi2

c     write(6,*)'# input eos'
c     read(5,*)eos
c     print *,'# eos=',eos
c     call readeos(eos)

c     fac = 146.51751415742*1d-3
c     n=5.0d0
c     de=0.05
c     nnn=100
c     do i=1,nnn
c        e=i*de+0.03
c        p=press(e/fac,n)*fac
c        print *,e,p,n*n0
c     end do

c     end

      subroutine readeos(optEoS,dir)

      character*80 dir
      integer eos,stabil,optEoS
      common /eqofst/ eos,stabil,anti
      real*8 e0,n0,Bag,hc3,pi2
      common /grstate/ e0,n0,Bag,hc3,pi2
      pi2 = dacos(-1d0)**2
      hc3 = 197.327053d0**3d0 
      e0 = 146.51751415742d0 
      n0 = 0.15891d0 
      Bag = 235d0**4 

      EoS=optEoS
      stabil = 0
c     EoS=CTOption(47)     ! Equation of state
      if(EoS.eq.1)EoS=5
c EoS: 1=5, 2= HG, 3=Bag model w/PT, 4=chiral w/PT, 5=chiral+HG
c
c     Read in equations of state. eos=2: hadronic EoS, 
c     eos=3: EoS with p.t. to QGP (bag modell mu = mu_b) eos=4:
c     chiral EoS with critical endpoint (mu = mu_q)
c     eos=5: chiral EoS but with hadronic T and mu_q.
c    

      if (eos.eq.2) then
         call readeos1(dir)
      else if (eos.eq.3) then
         call readeos2(dir)
      else if (eos.eq.4) then
         call readeos3(dir)
      else if (eos.eq.5) then
         call readeos3(dir)
      end if

      end

c------------------------------------------------------------------------
c     
c     subroutine-subprograms which read in EoS matrices.
c     
c     readeos1 reads in pure hadronic EoS, readeos2 reads in the 
c     EoS with phase transition to the QGP.
c     Both are used in the main program.
c
      subroutine readeos1(dir)
      character*80 dir
      integer leng
      real*8 t,mu,e,p,n,s,mstar,mus,lam
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400),msttab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)
      integer in, ie, j
c     
c     Common-Blocks.
c     
      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3
      
 

      leng=index(dir,' ')-1

      open(unit=53,
     $     file=dir(:leng)//'eosfiles/hadgas_eos.dat')
      open(unit=54,
     $     file=dir(:leng)//'eosfiles/hg_eos_small.dat')     
      open(unit=55,
     $     file=dir(:leng)//'eosfiles/hg_eos_mini.dat') 

c  Read in EoS.
c
      do 1169 in = 0,400,1
         j = 53
         do 1116 ie = 0,2000,1
            read(j,7787) t,mu,e,p,n,s,mus,lam
            ptab(ie,in) = p
            ttab(ie,in) = t
            mutab(ie,in) = mu
            stab(ie,in) = s
            mustab(ie,in) = mus
            cstab(ie,in) = lam

 1116    continue
 1169 continue
      close(53)

      do 1269 in = 0,200,1
         j = 54
         do 1216 ie = 0,200,1
            read(j,7787) t,mu,e,p,n,s,mus,lam
            ptab2(ie,in) = p
            ttab2(ie,in) = t
            mutab2(ie,in) = mu
            stab2(ie,in) = s
            mustab2(ie,in) = mus
            cstab2(ie,in) = lam

 1216    continue
 1269 continue
      close(54)      

      do 1569 in = 0,200,1
         j = 55
         do 1516 ie = 0,200,1
            read(j,7787) t,mu,e,p,n,s,mus,lam
            ptab3(ie,in) = p
            ttab3(ie,in) = t
            mutab3(ie,in) = mu
            stab3(ie,in) = s
            mustab3(ie,in) = mus
            cstab3(ie,in) = lam
 1516    continue
 1569 continue
      close(55)


      return
 7787 format(2(1x,f8.3),6(1x,e15.7))

      end
c
c-------------------------------------------------------------------
c     
      subroutine readeos2(dir)
      character*80 dir
      integer leng
      real*8 t,mu,e,p,n,s,mstar,lam,mus
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400),msttab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200)
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)
      integer in, ie, j
c
c     Common-blocks.
c     

      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3      
      

c     
c     Open files of EoS table.
c     
      leng=index(dir,' ')-1

      open(unit=47,
     $     file=dir(:leng)//'eosfiles/hadgas_eos.dat')
      open(unit=48,
     $     file=dir(:leng)//'eosfiles/hg_eos_small.dat')     
      open(unit=49,
     $     file=dir(:leng)//'eosfiles/hg_eos_mini.dat') 

      open(unit=51,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn00.dat')
      open(unit=52,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn05.dat')
      open(unit=53,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn10.dat')
      open(unit=54,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn15.dat')
      open(unit=55,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn20.dat')
      open(unit=56,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn25.dat')
      open(unit=57,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn30.dat')
      open(unit=58,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn35.dat')
      open(unit=59,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn40.dat')
      open(unit=60,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn45.dat')
      open(unit=61,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn50.dat')
      open(unit=62,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn55.dat')
      open(unit=63,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn60.dat')
      open(unit=64,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn65.dat')
      open(unit=65,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn70.dat')
      open(unit=66,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn75.dat')
      open(unit=67,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn80.dat')
      open(unit=68,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn85.dat')
      open(unit=69,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn90.dat')
      open(unit=70,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn95.dat')
      open(unit=71,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn100.dat')
      open(unit=72,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn105.dat')
      open(unit=73,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn110.dat')
      open(unit=74,
     $     file=dir(:leng)//'eosfiles/qgpeos/qgpn115.dat')
c
c     Read in EoS.
c
      do 1009 in = 0,239,1
         j = 51 + in/10
         do 1010 ie = 0,200,1
            read(j,7777) t,mu,e,p,n,s,mstar,lam
            ptab(ie,in) = p
            ttab(ie,in) = t
            mutab(ie,in) = mu
            stab(ie,in) = s
            msttab(ie,in) = mstar
            lamtab(ie,in) = lam
 1010   continue
 1009 continue

         do 4669 in = 0,400,1
            j = 47
            do 4616 ie = 0,2000,1
               read(j,7778) t,mu,e,p,n,s,mus,lam
               ttab(ie,in) = t
               mutab(ie,in) = mu
               mustab(ie,in) = mus
 4616       continue
 4669    continue
         
         do 4769 in = 0,200,1
            j = 48
            do 4716 ie = 0,200,1
               read(j,7778) t,mu,e,p,n,s,mus,lam
               ttab2(ie,in) = t
               mutab2(ie,in) = mu
               mustab2(ie,in) = mus
 4716       continue
 4769    continue
 
           do 4869 in = 0,200,1
            j = 49
            do 4816 ie = 0,200,1
               read(j,7778) t,mu,e,p,n,s,mus,lam
               ttab3(ie,in) = t
               mutab3(ie,in) = mu
               mustab3(ie,in) = mus
 4816       continue
 4869    continue
  

      close(47)
      close(48) 
      close(49)

      close(51)
      close(52) 
      close(53)
      close(54)
      close(55)
      close(56)
      close(57)
      close(58)
      close(59)
      close(60)
      close(61)
      close(62)
      close(63)
      close(64)
      close(65)
      close(66)
      close(67)
      close(68)
      close(69)
      close(70)
      close(71)
      close(72)
      close(73)
      close(74)
      return
 7777 format(2(1x,f8.3),4(1x,e15.7),2(1x,f6.3))
 7778 format(2(1x,f8.3),6(1x,e15.7))
      end
c     
c------------------------------------------------------------------------
c

      subroutine readeos3(dir)
      character*80 dir
      integer leng
      real*8 t,mu,e,p,n,s,mstar,lam,mus
      real*8 ptab(0:2000,0:400),ttab(0:2000,0:400),lamtab(0:200,0:239)
      real*8 mutab(0:2000,0:400),stab(0:2000,0:400),msttab(0:2000,0:400)
      real*8 ptab2(0:200,0:200),ttab2(0:200,0:200)
      real*8 mutab2(0:200,0:200),stab2(0:200,0:200),msttab2(0:200,0:200)
      real*8 mustab(0:2000,0:400),mustab2(0:200,0:200) 
      real*8 cstab(0:2000,0:400),cstab2(0:200,0:200)
      real*8 ptab3(0:200,0:200),ttab3(0:200,0:200)
      real*8 mutab3(0:200,0:200),stab3(0:200,0:200),msttab3(0:200,0:200)
      real*8 mustab3(0:200,0:200),cstab3(0:200,0:200)
      integer in, ie, j, eos

c     
c     Common-blocks.
c
      common /eqofst/ eos,stabil,anti
      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3
c     
c     Open files of EoS table.
c     
      leng=index(dir,' ')-1
      open(unit=51,
     $     file=dir(:leng)//'eosfiles/chiraleos.dat')
      open(unit=52,
     $     file=dir(:leng)//'eosfiles/chiralsmall.dat')     
      open(unit=56,
     $     file=dir(:leng)//'eosfiles/chiralmini.dat')  



c...Temperature [MeV],
c...1/3 baryo chemical potential [MeV],
c...energy density [e_0=146.517 MeV/fm^3],
c...pressure [e_0],
c...baryon density [r_0=0.15891/fm^3],
c...entropy density [r_0],
c...strange chemical potential [MeV], 
c...QGP fraction

c...The tables are ordered in energy density and baryon density. Step size
c...is equal de= 0.5 and drho=0.1 with e_max=1000 and rho_max=40.




c  Read in EoS.
c
         do 1109 in = 0,400,1
            j = 51
            do 1110 ie = 0,2000,1
            read(j,7777) t,mu,e,p,n,s,mstar,lam
            ptab(ie,in) = p
            ttab(ie,in) = t
            mutab(ie,in) = mu
            stab(ie,in) = s
            msttab(ie,in) = mstar
            cstab(ie,in) = lam
 1110    continue
 1109 continue
      close(51)

      do 1209 in = 0,200,1
         j = 52
         do 1210 ie = 0,200,1
            read(j,7777) t,mu,e,p,n,s,mstar,lam
            ptab2(ie,in) = p
            ttab2(ie,in) = t
            mutab2(ie,in) = mu
            stab2(ie,in) = s
            msttab2(ie,in) = mstar
            cstab2(ie,in) = lam
 1210    continue
 1209 continue
      close(52)

      do 1509 in = 0,200,1
         j = 56
         do 1510 ie = 0,200,1
            read(j,7777) t,mu,e,p,n,s,mstar,lam
            ptab3(ie,in) = p
            ttab3(ie,in) = t
            mutab3(ie,in) = mu
            stab3(ie,in) = s
            msttab3(ie,in) = mstar
            cstab3(ie,in) = lam
 1510    continue
 1509 continue
      close(56)
      
      if(eos.eq.5) then
         open(unit=53,
     $        file=dir(:leng)//'eosfiles/hadgas_eos.dat')
         open(unit=54,
     $        file=dir(:leng)//'eosfiles/hg_eos_small.dat')     
         open(unit=55,
     $        file=dir(:leng)//'eosfiles/hg_eos_mini.dat')     
                
         
         do 1669 in = 0,400,1
            j = 53
            do 1616 ie = 0,2000,1
               read(j,7777) t,mu,e,p,n,s,mus,lam
               ttab(ie,in) = t
               mutab(ie,in) = mu
               mustab(ie,in) = mus
 1616       continue
 1669    continue
         close(53)
         
         do 1769 in = 0,200,1
            j = 54
            do 1716 ie = 0,200,1
               read(j,7777) t,mu,e,p,n,s,mus,lam
               ttab2(ie,in) = t
               mutab2(ie,in) = mu
               mustab2(ie,in) = mus
 1716       continue
 1769    continue
         close(54)
 
           do 1869 in = 0,200,1
            j = 55
            do 1816 ie = 0,200,1
               read(j,7777) t,mu,e,p,n,s,mus,lam
               ttab3(ie,in) = t
               mutab3(ie,in) = mu
               mustab3(ie,in) = mus
 1816       continue
 1869    continue
         close(55)   
      endif
         
      return

 7777 format(2(1x,f8.3),6(1x,e15.7))
      
      end
c
c----------------------------------------------------------------------
c
      function press(e,n)
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
      real*8 e,n,press,et0
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
            p1 = ptab(idint(e/de),idint(n/dn))
            p2 = ptab(idint(e/de)+1,idint(n/dn))
            p3 = ptab(idint(e/de),idint(n/dn)+1)
            p4 = ptab(idint(e/de)+1,idint(n/dn)+1)
            p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
            p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
            p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
            p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
      
            press = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     press = 0.25d0*(p13+p24+p12+p34)
         else if (e.ge.20d0) then
            press = (e-4d0*B/e0/hc3)/3d0
         end if
         
      else if ((eos.eq.4).or.(eos.eq.2).or.(eos.eq.5)) then
         if (e.le.1000.0d0) then

            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = ptab3(idint(e/de),idint(n/dn))
               p2 = ptab3(idint(e/de)+1,idint(n/dn))
               p3 = ptab3(idint(e/de),idint(n/dn)+1)
               p4 = ptab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2 
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               press = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     write(99,*) p1, ptab3(0,0) 
            end if
            
            if (((e.lt.10.0d0).and.(n.lt.2.0d0)).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = ptab2(idint(e/de),idint(n/dn))
               p2 = ptab2(idint(e/de)+1,idint(n/dn))
               p3 = ptab2(idint(e/de),idint(n/dn)+1)
               p4 = ptab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2 
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               press = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then

               de = 0.5d0
               dn = 0.1d0
 
               p1 = ptab(idint(e/de),idint(n/dn))
               p2 = ptab(idint(e/de)+1,idint(n/dn))
               p3 = ptab(idint(e/de),idint(n/dn)+1)
               p4 = ptab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2 
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               press = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13  
            end if
         else if (e.ge.1000.0d0) then
c     JS muss bestimmt werden
            press = e/3.0d0
         end if
      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1).and.(n.ge.-1d0))
     $     then
         press = 0d0
      end if
      
      if (press.lt.0d0) press=0d0
      
      
      
      
      
      if (anti.eq.1)then
         n=-n
      end if
      return
      end

c     
c----------------------------------------------------------------------
c
      function entro(e,n)
c     INCLUDE 'defs.f'
      integer ngr
      PARAMETER (ngr=200)
c
c  This function-subprogram determines the entropy density of the
c  underlying EoS.
c  It is used in subroutines fileo, prop3d, and the main program.
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
      integer eos,stabil,antie
      real*8 chem
c
c  Type declarations for variables used in function entro.
c
      real*8 e,n,entro,et0,temp
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
      
      
      antie = 0

      if (n.lt.0)then
         n=-n
         antie = 1
      end if
      
c
c  Calculate entropy density (in units of n0).
c  If flag eos=1, use ideal gas equation of state.   
c  If flag eos=[else], use tabellized equation of state.
c     
c     arithmetic mean interpolation (yields better values for
c     interpolation at small e,n, but is worse for all other
c     thermodynamic quantities, also at larger e,n)
c     
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
            
            entro = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     entro = 0.25d0*(p13+p24+p12+p34)    
            
          else if (e.gt.20d0) then
             entro = 4d0*(e-B/e0/hc3)*e0/3d0/temp(e,n)/n0-
     $            n*chem(e,n)/temp(e,n)
          else
             entro = 0d0
          end if
          
       else if ((eos.eq.4).or.(eos.eq.5).or.(eos.eq.2)) then
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

               entro = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if          
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $             ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
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

               entro = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
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
               
               entro = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if            
            
         else if (e.ge.1000.0d0) then
c     JS muss bestimmt werden
            entro = 100d0
         end if
         
      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1).and.(n.ge.-1d0)) 
     $     then
         entro = 0d0
      end if
      
c     if (entro.lt.0d0) entro=0d0
      
      
      
      if (antie.eq.1)then
         n=-n
      end if
      
      return
      end
c
c----------------------------------------------------------------------
c
      function temp(e,n)
c     INCLUDE 'defs.f'
      integer ngr
      PARAMETER (ngr=200)
c     
c     This function-subprogram determines the temperature of the
c     underlying EoS.
c     It is used in subroutines fileo and entro.
c
c     Type declarations for variables in common-blocks.
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
      integer eos,stabil,antit
c
c     Type declarations for variables used in function temp.
c     
      real*8 e,n,temp,et0,mu
      real*8 de,dn,p1,p2,p3,p4,p13,p24
c     real*8 p12,p34
c     
c     Common-blocks.
c     
      common /grstate/ e0,n0,B,hc3,pi2
      common /eqofst/ eos,stabil,anti
      common /eos/ ptab,ttab,mutab,stab,lamtab,ptab2,ttab2,mutab2,stab2,
     $     mustab,mustab2,cstab,cstab2,cstab3,ptab3,ttab3,mutab3
     $     ,stab3,mustab3


c     
c     Fermi energy density in the QGP for given baryon density
c     (in units of e0).
c     
      et0 = 1d0/54d0/pi2*dabs(40.5d0*pi2*n*n0*hc3)**(4d0/3d0)+B
      et0 = et0/e0/hc3
      
      
      antit = 0
      if (n.lt.0)then
         n=-n
         antit = 1
      end if
c
c     Calculate temperature (in MeV).   
c     If flag eos=[else], use tabellized equation of state.
c     
c     
c     If e and n are both smaller than e0, n0, the nuclei are in the
c     ground state. The pressure is set to zero by hand, if flag
c     stabil = 1.
c     
      if (eos.eq.3) then
         if (e.le.20d0) then

            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = ttab3(idint(e/de),idint(n/dn))
               p2 = ttab3(idint(e/de)+1,idint(n/dn))
               p3 = ttab3(idint(e/de),idint(n/dn)+1)
               p4 = ttab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = ttab2(idint(e/de),idint(n/dn))
               p2 = ttab2(idint(e/de)+1,idint(n/dn))
               p3 = ttab2(idint(e/de),idint(n/dn)+1)
               p4 = ttab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then
               de = 0.5d0
               dn = 0.1d0
               p1 = ttab(idint(e/de),idint(n/dn))
               p2 = ttab(idint(e/de)+1,idint(n/dn))
               p3 = ttab(idint(e/de),idint(n/dn)+1)
               p4 = ttab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
         else if (e.gt.20d0) then
            call findqgp(e,n,temp,mu)
         else
            temp = 0d0
         end if
      else if ((eos.eq.4).or.(eos.eq.5).or.(eos.eq.2)) then      
         if (e.le.1000.0d0) then

            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = ttab3(idint(e/de),idint(n/dn))
               p2 = ttab3(idint(e/de)+1,idint(n/dn))
               p3 = ttab3(idint(e/de),idint(n/dn)+1)
               p4 = ttab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $              ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = ttab2(idint(e/de),idint(n/dn))
               p2 = ttab2(idint(e/de)+1,idint(n/dn))
               p3 = ttab2(idint(e/de),idint(n/dn)+1)
               p4 = ttab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then
               de = 0.5d0
               dn = 0.1d0

               p1 = ttab(idint(e/de),idint(n/dn))
               p2 = ttab(idint(e/de)+1,idint(n/dn))
               p3 = ttab(idint(e/de),idint(n/dn)+1)
               p4 = ttab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               temp = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if 
            
         else if (e.ge.1000.0d0) then
c     JS muss bestimmt werden
            Temp = 400d0
c     else
            
c     temp = 0d0
            
         end if
         
      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1).and.(n.ge.-1d0))
     $     then
         temp = 0d0
      end if
      
      
      if (antit.eq.1)then
         n=-n
      end if    
      return
      end     
      
      
      
      
        

c
c----------------------------------------------------------------------
c
      function chem(e,n)
c     INCLUDE 'defs.f'
      integer ngr
      PARAMETER (ngr=200)
c
c  This function-subprogram determines the chemical potential of the
c  underlying EoS.
c  It is used in subroutine fileo.
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
      integer eos,stabil,antic
c
c  Type declarations for variables used in function temp.
c
      real*8 e,n,chem,et0,temp
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
      

      antic = 0
      if (n.lt.0)then
         n=-n
         antic = 1
      end if
c
c  Calculate chemical potential (in MeV).

c     
c  If flag eos=[else], use tabellized equation of state.
c     

c     If e and n are both smaller than e0, n0, the nuclei are in the
c     ground state. The pressure is set to zero by hand, if flag
c     stabil = 1.
c     
      if (eos.eq.3) then
         if (e.le.20d0) then
            
            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = mutab3(idint(e/de),idint(n/dn))
               p2 = mutab3(idint(e/de)+1,idint(n/dn))
               p3 = mutab3(idint(e/de),idint(n/dn)+1)
               p4 = mutab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = mutab2(idint(e/de),idint(n/dn))
               p2 = mutab2(idint(e/de)+1,idint(n/dn))
               p3 = mutab2(idint(e/de),idint(n/dn)+1)
               p4 = mutab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then       
               de = 0.5d0
               dn = 0.1d0
               p1 = mutab(idint(e/de),idint(n/dn))
               p2 = mutab(idint(e/de)+1,idint(n/dn))
               p3 = mutab(idint(e/de),idint(n/dn)+1)
               p4 = mutab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13

            end if
            
         else if (e.ge.20d0) then
            call findqgp(e,n,temp,chem)
         else
            chem = 0d0
         end if

      else if ((eos.eq.4).or.(eos.eq.5).or.(eos.eq.2)) then        
         if (e.lt.1000.0d0) then


            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = mutab3(idint(e/de),idint(n/dn))
               p2 = mutab3(idint(e/de)+1,idint(n/dn))
               p3 = mutab3(idint(e/de),idint(n/dn)+1)
               p4 = mutab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = mutab2(idint(e/de),idint(n/dn))
               p2 = mutab2(idint(e/de)+1,idint(n/dn))
               p3 = mutab2(idint(e/de),idint(n/dn)+1)
               p4 = mutab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then       
               de = 0.5d0
               dn = 0.1d0

               p1 = mutab(idint(e/de),idint(n/dn))
               p2 = mutab(idint(e/de)+1,idint(n/dn))
               p3 = mutab(idint(e/de),idint(n/dn)+1)
               p4 = mutab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               chem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
               
            end if
         else if (e.gt.1000.0d0) then
c     JS muss bestimmt werden
            chem = 0.0d0
c     else
c     chem = 0d0
         end if
         
      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1)) then
         chem = e0/n0/3
      end if
      
      
      
      
      if (antic.eq.1)then
         n=-n
         chem = -chem
      end if
      return
      end
c     
c----------------------------------------------------------------------
c
      function schem(e,n)
c     INCLUDE 'defs.f'
      integer ngr
      PARAMETER (ngr=200)
c
c  This function-subprogram determines the chemical potential of the
c  underlying EoS.
c  It is used in subroutine fileo.
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
      integer eos,stabil,antic
c
c  Type declarations for variables used in function schem.
c
      real*8 e,n,chem,et0,temp,schem
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
      

      antic = 0
      if (n.lt.0)then
         n=-n
         antic = 1
      end if
c
c  Calculate chemical potential (in MeV).     
c  If flag eos=[else], use tabellized equation of state.
c     

c     If e and n are both smaller than e0, n0, the nuclei are in the
c     ground state. The pressure is set to zero by hand, if flag
c     stabil = 1.
c     
      if (eos.eq.3) then
         if (e.le.20d0) then
            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = mustab3(idint(e/de),idint(n/dn))
               p2 = mustab3(idint(e/de)+1,idint(n/dn))
               p3 = mustab3(idint(e/de),idint(n/dn)+1)
               p4 = mustab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               schem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = mustab2(idint(e/de),idint(n/dn))
               p2 = mustab2(idint(e/de)+1,idint(n/dn))
               p3 = mustab2(idint(e/de),idint(n/dn)+1)
               p4 = mustab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               schem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then        
               
               de = 0.5d0
               dn = 0.1d0
               p1 = mustab(idint(e/de),idint(n/dn))
               p2 = mustab(idint(e/de)+1,idint(n/dn))
               p3 = mustab(idint(e/de),idint(n/dn)+1)
               p4 = mustab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               schem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13         
              
            end if
         else if (e.ge.20d0) then
            schem = 0d0
         end if
         
      else if ((eos.eq.4).or.(eos.eq.5).or.(eos.eq.2)) then        
         if (e.lt.1000.0d0) then
            
            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = mustab3(idint(e/de),idint(n/dn))
               p2 = mustab3(idint(e/de)+1,idint(n/dn))
               p3 = mustab3(idint(e/de),idint(n/dn)+1)
               p4 = mustab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
                  
               schem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = mustab2(idint(e/de),idint(n/dn))
               p2 = mustab2(idint(e/de)+1,idint(n/dn))
               p3 = mustab2(idint(e/de),idint(n/dn)+1)
               p4 = mustab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               

               schem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then      
               de = 0.5d0
               dn = 0.1d0

               p1 = mustab(idint(e/de),idint(n/dn))
               p2 = mustab(idint(e/de)+1,idint(n/dn))
               p3 = mustab(idint(e/de),idint(n/dn)+1)
               p4 = mustab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
               
               schem = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     chem = 0.25d0*(p13+p24+p12+p34)    
            end if
         else if (e.gt.1000.0d0) then
c     JS muss bestimmt werden
            schem = 0d0
c     else
c     chem = 0d0
         end if
         
      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1)) then
         schem = 0
      end if
      
      
      
      
      if (antic.eq.1)then
         n=-n
         schem = schem
      end if
      return
      end
c 
c----------------------------------------------------------------------
c
      function cs(e,n)
c     INCLUDE 'defs.f'
      integer ngr
      PARAMETER (ngr=200)
c
c  This function-subprogram determines the chemical potential of the
c  underlying EoS.
c  It is used in subroutine fileo.
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
      integer eos,stabil,antic
c
c  Type declarations for variables used in function schem.
c
      real*8 e,n,chem,et0,temp,schem,cs
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
      

      antic = 0
      if (n.lt.0)then
         n=-n
         antic = 1
      end if
c
c  Calculate chemical potential (in MeV).
c     
c  If flag eos=[else], use tabellized equation of state.
c     

c     If e and n are both smaller than e0, n0, the nuclei are in the
c     ground state. The pressure is set to zero by hand, if flag
c     stabil = 1.
c     
      if (eos.eq.3) then
         if (e.le.20d0) then
            de = 0.2d0
            dn = 0.1d0
            p1 = cstab(idint(e/de),idint(n/dn))
            p2 = cstab(idint(e/de)+1,idint(n/dn))
            p3 = cstab(idint(e/de),idint(n/dn)+1)
            p4 = cstab(idint(e/de)+1,idint(n/dn)+1)
            p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
            p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
            p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
            p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3

            cs = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c            chem = 0.25d0*(p13+p24+p12+p34)         
            
         else if (e.ge.20d0) then
            call findqgp(e,n,temp,chem)
         else
            cs = 0d0
         end if

      else if ((eos.eq.4).or.(eos.eq.5).or.(eos.eq.2)) then        
         if (e.lt.1000.0d0) then

            if((e.lt.0.1d0).and.(n.lt.0.02d0)) then
               de = 0.0005d0
               dn = 0.0001d0
               p1 = cstab3(idint(e/de),idint(n/dn))
               p2 = cstab3(idint(e/de)+1,idint(n/dn))
               p3 = cstab3(idint(e/de),idint(n/dn)+1)
               p4 = cstab3(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3


               cs = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if((e.lt.10.0d0).and.(n.lt.2.0d0).and.
     $           ((e.ge.0.1d0).or.(n.ge.0.02d0))) then
               de = 0.05d0
               dn = 0.01d0
               p1 = cstab2(idint(e/de),idint(n/dn))
               p2 = cstab2(idint(e/de)+1,idint(n/dn))
               p3 = cstab2(idint(e/de),idint(n/dn)+1)
               p4 = cstab2(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            


               cs = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
            end if
            if ((e.ge.10.0d0).or.(n.ge.2.0d0)) then 
               de = 0.5d0
               dn = 0.1d0
               
               p1 = cstab(idint(e/de),idint(n/dn))
               p2 = cstab(idint(e/de)+1,idint(n/dn))
               p3 = cstab(idint(e/de),idint(n/dn)+1)
               p4 = cstab(idint(e/de)+1,idint(n/dn)+1)
               p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
               p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
               p12 = (p2-p1)/de*(e-dfloat(idint(e/de))*de) + p1
               p34 = (p4-p3)/de*(e-dfloat(idint(e/de))*de) + p3
            
               cs = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     chem = 0.25d0*(p13+p24+p12+p34)    
            end if
         else if (e.gt.1000.0d0) then
c     JS muss bestimmt werden
            cs = 1/sqrt(3d0)
c     else
c     chem = 0d0
         end if
         
      end if
      
      if ((e.le.1d0).and.(n.le.1d0).and.(stabil.eq.1)) then
         cs = 0
      end if
     
 
      
      
      
      if (antic.eq.1)then
         n=-n
         cs = cs
      end if
      return
      end
c     
c----------------------------------------------------------------------    
c----------------------------------------------------------------------
c
      function lambda(e,n)
c     INCLUDE 'defs.f'
      integer ngr
      PARAMETER (ngr=200)
c
c  This function-subprogram determines the fraction of QGP of the
c  underlying EoS.
c  It is used in subroutine prop3d and the main program.
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
      integer eos,stabil,antil
c
c  Type declarations for variables used in function lambda.
c
      real*8 e,n,et0,lambda
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
      

      antil = 0

      if (n.lt.0)then
         n=-n
         antil = 1
      end if
      
c     
c     Calculate fraction of QGP.
c     If eos = 1 or eos = 2, lambda is zero.
c     
      if ((eos.eq.2).or.(eos.eq.4)) then
         lambda = 0d0
c     
c     If eos=[else], use tabellized equation of state.
c     
      else
         if (e.le.20d0) then
            de = 0.1d0
            dn = 0.05d0
            p1 = lamtab(idint(e/de),idint(n/dn))
            p2 = lamtab(idint(e/de)+1,idint(n/dn))
            p3 = lamtab(idint(e/de),idint(n/dn)+1)
            p4 = lamtab(idint(e/de)+1,idint(n/dn)+1)
            p13 = (p3-p1)/dn*(n-dfloat(idint(n/dn))*dn) + p1
            p24 = (p4-p2)/dn*(n-dfloat(idint(n/dn))*dn) + p2
c     

c     
c     Linear interpolation.
c     
            lambda = (p24-p13)/de*(e-dfloat(idint(e/de))*de) + p13
c     
c     arithmetic mean interpolation (yields better values for
c     interpolation at small e,n, but is worse for all other
c     thermodynamic quantities, also at larger e,n)
c     
c     lambda = 0.25d0*(p13+p24+p12+p34) 
c     
         else if (e.gt.20d0) then
            lambda = 1d0
         else
            lambda = 0d0
         end if
      end if
      if (lambda.gt.1d0) lambda=1d0
      
      if (antil.eq.1)then
         n=-n
      end if
      return
      end
c     
c------------------------------------------------------------------------
c
      subroutine findqgp(energ,dens,t,mu)
c     INCLUDE 'defs.f'
      integer ngr
      PARAMETER (ngr=200)
c
c  This subroutine-subprogram finds T, mu in the QGP for given e,n.
c  It is used in subroutines temp and chem.
c
      real*8 e0,n0,B,hc3,pi2
      real*8 energ,dens,t,mu
      real*8 pi
      real*8 a1,a2,a3,a4
      real*8 et0
      integer l
c     
c     Common-blocks.
c
      common /grstate/ e0,n0,B,hc3,pi2
c     
c  Define constants.
c     
      pi = dsqrt(pi2)
      l = 0
c     
c     Convert energy density and baryon density, 
c     so far given in units of e0 and n0, into MeV**4 
c     resp. MeV**3, so that T and mu emerge in units of MeV.
c     
      energ = energ*e0*hc3
      dens = dens*n0*hc3
c     
c     If e < e(T=0) (the Fermi energy density in the QGP for given
c     baryon density), set everything to zero.
c     
      et0 = 1d0/54d0/pi2*dabs(40.5d0*pi2*dens)**(4d0/3d0)+B
      if (energ.lt.et0) then
         t = 0d0
        mu = 0d0
      else
c     
c     If density is larger than zero, complicated sixth
c     order polynomial in mu has to be solved...
c     
         if (dens.gt.0d0) then
c     
c     Calculate coefficients of sixth order polynomial.
c     
          a1 = 4d0/1215d0/pi2
          a2 = - 4d0/15d0*dens
          a3 = energ - B
          a4 = - 999d0*pi2*dens*dens/40d0
          call froot(mu,a1,a2,a3,a4,1d-7,l)
c
c     After having found mu, calculate T.
c     
          t = dsqrt(4.5d0*dens/mu - mu*mu/9d0/pi2)
c     
c     ...if density is zero, mu = 0 and t is fourth root
c     of e-B (up to a constant). 
c     
        else
           mu = 0d0
           if(energ.ge.B) then
              t = ((energ-B)*30d0/37d0/pi2)**0.25d0
           else
            t=0d0
         end if
      end if
      end if 
c     
c     Finally, convert e, n back into dimensionless
c     units, so that the original quantities will
c     be returned to the main program.
c     
      energ = energ/e0/hc3
      dens = dens/n0/hc3
      return
      end 
c     
c     
c------------------------------------------------------------------------
c
      subroutine froot(mu,a1,a2,a3,a4,acc,l)
c     INCLUDE 'defs.f'
      integer ngr
      PARAMETER (ngr=200)
c     
c  This subroutine-subprogram finds the root of the sixth order 
c  polynomial in mu.
c  It is used in subprogram findqgp.
c
      real*8 e0,n0,B,hc3,pi2
      real*8 mu,a1,a2,a3,a4,acc
      real*8 mua,mub,fa,fb,fc
      integer i,l
c
c  Common-blocks.
c
      common /grstate/ e0,n0,B,hc3,pi2
c
      i = 0
c
c  Lower boundary for interval search in mu is mua=0.
c
      mua = 0d0
      fa = a4
c
c  Upper boundary for interval search in mu corresponds
c  to the maximum possible mu, ie, that mu which
c  corresponds to T=0 and baryon density 12*n0.
c
      mub = (40.5d0*pi2*12d0*n0*hc3)**(1d0/3d0)
      fb = a4 + a3*mub*mub + a2*mub**3 + a1*mub**6
 22   mu = (mua + mub)*0.5d0
      fc = a4 + a3*mu*mu + a2*mu**3 + a1*mu**6
      i = i+1
      if (i.gt.1000) then
        write(6,*) ' froot: more than 1000 iterations to find mu'
        write(6,*) '  fc/e0= ',fc/e0/hc3
        goto 23
      end if
      if (dabs(fa)/e0/hc3.lt.acc) then
        mu = mua
        goto 23
      end if
      if (dabs(fb)/e0/hc3.lt.acc) then
        mu = mub
        goto 23
      end if  
      if (fb*fc.lt.0d0) then
        mua = mu
        fa = fc
      else if (fa*fc.lt.0d0) then
        mub = mu
        fb = fc
      else
c        write(6,*) ' froot: 2 or no zeros in search interval'   
        l = 1
        return
      end if
      if (dabs(fc)/e0/hc3.gt.acc) goto 22
 23   continue
      l = 0
      return
      end 
c
