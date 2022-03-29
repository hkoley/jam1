
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension kf0(43)
      data kf0/2112,2212,3112,3122,3212,3222,3312,3322,3334,
     &  4122,4132,4232,
     &  111,211,-211,221,333,130,310,321,-321,311,
     &  411,-411,421,431,-431,443,513,523,533,
     &  4222,4212,4112,4312,4322,
     &  425,435,511,521,531,541,551/

      do i=1,43
       kf=kf0(i)
       kc=jamcomp(kf)
       ist=kchg(kc,7)*isign(1,kf)

c     print *,kf,jamchge(kf),kfprop(kf,1),ist,kfprop(kf,3)
c     print *,'    case: ',kf,'z = ',jamchge(kf)
      write(6,800)kf,jamchge(kf),ist
      end do

800   format(5x,'case',i4,': z= opt=1 ? ',i2,':',i2,';')

      end
