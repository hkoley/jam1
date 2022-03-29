      implicit double precision(a-h, o-z)
      include 'jam1.inc'
      include 'jam2.inc'
      character chau*16
      character frame*8,proj*8,targ*8,cwin*15
      logical dump 
c     data dump/.false./
      data dump/.true./

      ioptbinary=2

      if(dump) then
        if(ioptbinary.eq.0) then
          open(33,file='phase.dat',status='new')
        else if(ioptbinary.eq.1) then
          open(33,file='phase.dat',status='new',
     &           form='unformatted',access='stream')
        else
         call jamccopen()
        endif
      endif

c....Initialize JAM.
      fname(1)='jam.cfg'
      call jaminit(nev,b1,b2,dt,nstp,frame,proj,targ,cwin)
      nevent=mstc(2)
      if(dump) then
        if(ioptbinary.eq.0) then
          write(33,800)nevent,pard(17),pard(5),pard(6),mstc(4)
        else if(ioptbinary.eq.1) then
          write(33)nevent,pard(17),pard(5),pard(6),mstc(4)
        else
          call jamccout1(nevent)
        endif
      endif

      if(mstc(50).ge.21) then
        iopt=mod(mstc(50),10)
        print *,'iopt=',iopt
        call readeos(iopt,fname(9))
      else if(mstc(50).ge.11) then
        call readEOStable(fname(9))
      endif

c==================================================
      do iev=1,nevent

        call jamevt(iev)
        if(mod(iev,100).eq.0) write(6,*)'event=',iev

c...Dump phase space data.
        if(dump) then
         ncol=mstd(41)+mstd(42)
         ncolbb=mstd(44)
         nonpart=0
         do i=1,nv
          if(abs(k(7,i)).eq.1) nonpart=nonpart+1
         end do
         npart=mstd(11)-nonpart
          if(ioptbinary.eq.0) then
            write(33,811)iev,nv,nbary,nmeson,pard(2),npart,ncol,ncolbb
            do i=1,nv
              write(33,813) k(1,i),k(2,i),k(7,i),p(5,i),(p(j,i),j=1,4),
     &         (r(j,i),j=1,5)
            end do

          else if(ioptbinary.eq.1) then

            write(33)iev,nv,nbary,nmeson,pard(2),npart,ncol,ncolbb
            do i=1,nv
              write(33) k(1,i),k(2,i),k(7,i),p(5,i),(p(j,i),j=1,4),
     &         (r(j,i),j=1,5)
            end do

          else
            call jamccout(iev,npart,ncol,ncolbb)
          endif

        endif

c....End simulation
      end do
 800  format('#',1x,i8,3(1x,f11.5),1x,i4)
 811  format('#',4(1x,i8),f6.2,1x,3(i6,1x))
 813  format(1x,i8,1x,i12,1x,i5,10(1x,1pe16.8))

      call jamfin
      print *,'mstd(199)=',mstd(199)


      if(dump) then
        if(ioptbinary.le.2) then
          close(33)
          if(ioptbinary.eq.1) call system('gzip -f '//'phase.dat')
        else
         call jamccclose()
        endif
      endif

      end
