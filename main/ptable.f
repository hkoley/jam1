      implicit none
      include 'jam2.inc'
      integer kcmin,kcmax,kc,kf,ispin,ns,ibary,ch,nc,nbot,id
      integer getiso,iso
      real*8 pm,width
      character*16 chap

      kcmin=132
      kcmax=366
      do kc=kcmin,kcmax
        kf=kchg(kc,4)
        if(kf.eq.130.or.kf.eq.310) cycle
        if(kc.ge.195.and.kc.le.251) cycle
        ispin=max(1,mod(kf,10)) ! spin 2J+1
        nc=0
        nbot=0
        ch=kchg(kc,1)/3
        ns=kchg(kc,7)
        ibary=kchg(kc,6)/3               ! baryon number times 3.
        pm=pmas(kc,1)                 ! particle mass
        width=pmas(kc,2)                 ! particle width
        id=kchg(kc,5)
        iso=getiso(id)

        call pjname(kf,chap)
        write(6,800)kf,chap,pm,width,ispin,ibary,ns,nc,nbot,iso,ch,0

        if(kchg(kc,3).eq.1.and.ibary.eq.0) then
          kf=-kf
          ns=-ns
          ch=-ch
        call pjname(kf,chap)
        write(6,800)kf,chap,pm,width,ispin,ibary,ns,nc,nbot,iso,ch,0
        endif

 800  format(i9,1x,a16,4x,2(f10.5,1x),8(i4,1x))

      end do

      write(6,*)7

      end

      integer function getiso(id)
      implicit none
      include 'jam2.inc'
      integer id
      getiso=0
      if(id.eq.id_pi) getiso=3
      if(id.eq.id_light0) getiso=1
      if(id.eq.id_light1) getiso=3
      if(id.eq.id_str) getiso=2
      if(id.eq.id_nucl) getiso=2
      if(id.eq.id_nucls) getiso=2
      if(id.eq.id_delt) getiso=4
      if(id.eq.id_delts) getiso=4
      if(id.eq.id_sigm) getiso=3
      if(id.eq.id_sigms) getiso=3
      if(id.eq.id_lamb) getiso=1
      if(id.eq.id_lambs) getiso=1
      if(id.eq.id_xi) getiso=2
      if(id.eq.id_xis) getiso=2
      if(id.eq.id_omega) getiso=1

      if(getiso.eq.0) then
        print *,'unknown id? id=',id,getiso
        stop
      endif

      end
