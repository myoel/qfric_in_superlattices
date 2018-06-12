! This small program will give the reflection coefficients of a Drude and Lorentz bulk,
! a superlattice with Drude-Lorentz and Lorentz-Drude layers and within the EMA.
! The data is stored separately for real and imaginary part for r(w,k=const.)
! and is put together for r(w=const.,k).


! In this module all shared parameters are defined:
module comref
 implicit none
 ! General parameters:
 real(8) :: pi ! Pi

 ! Defining the superlattice:
 real(8), public :: dd  ! Thickness of the superlattice
 real(8), public :: rho ! filling factor rho= dm/(dm+di),
                        ! with d1/2 is the thickness of the 
                        ! components
                       
 ! Defining the metal component (Drude model):
  real(8), public :: einf1! Background constant
 real(8), public :: wp1 ! Plasma frequency
  real(8), public :: w01 ! Central frequency
 real(8), public :: wsp1! Surface Plasmon freqeuncy
 real(8), public :: g1 ! Damping constant

 ! epD(w) = einf - wp^2/(w^2 + i*g*w)
 
 ! Defining the insulator component (multiple Lorentz model):
 real(8), public :: einf2! Background constant
 real(8), public :: wp2 ! Plasma frequency
 real(8), public :: w02 ! Central frequency
 real(8), public :: g2  ! Damping constant
 real(8), public :: wp3 ! Plasma frequency
 real(8), public :: w03 ! Central frequency
 real(8), public :: g3  ! Damping constant
 real(8), public :: wp4 ! Plasma frequency
 real(8), public :: w04 ! Central frequency
 real(8), public :: g4  ! Damping constant
 ! epL(w) = einf - wp^2/(w^2 - w0^2 + i*g*w)
 
 ! Doping which yields in Drude-Lorentz model:
 logical, public :: doped ! If TRUE then Drude term is added
 real(8), public :: wpd ! Plasma freqeuncy
 real(8), public :: gd  ! Damping constant
 
 ! epDL(w) = epL(w) - wp^2/(w^2 + i*g*w)

 ! Defining the order of the layers:
 logical :: toplog
 logical :: bulk ! if TRUE then only the bulk of the top-most
                 ! material is calculated.
  logical :: EMA  ! if TRUE then only EMA is performed
  logical :: rimag! Chooses either imaginary or real
                 ! part of the reflection coefficient
  logical:: near ! If TRUE then only near-field approximation 
                 ! is considered.
  logical::slab  ! If TRUE then only a slab with thickness d*rho is calculated
 ! File names:
 character(1), public :: para
 
  real(8), public :: dx2, dy2, dz2 ! Cartesian form of dipole
  
  real(8) ::dtheta, dphi
end module comref

program refcalc
 use comref
 implicit none
 Integer :: i, maxi
 real(8) :: za ! distance to the surface
 real(8) :: v,k,w, wcons
 real(8) :: sta, sto, spac ! Start and stop and spacing of the calculated points
 complex*16 :: resDLs, resDs, resLs, resLDs, resEMAp, resEMAs
 complex*16 :: resDLp, resDp, resLp, resLDp, resDslabp, resDslabs
 real(8), external :: alpha2, metricSInat,metricNatSI
 real :: t1,t2 ! variables to measure the CPU-time
 character(1) :: vari,answer
 character(29):: filenws, filenwp, filenks, filenkp
 complex*16 ::   epM
 
 ! Create File Name:

 print*,'Which input file ? (1,2,3,...)'
 read(*,*) para
 write(filenws,"(A16,A1,A4)") '../output/refws_', para, '.dat'
 print*,filenws
  write(filenwp,"(A16,A1,A4)") '../output/refwp_', para, '.dat'
 print*,filenwp
 write(filenks,"(A16,A1,A4)") '../output/refks_', para, '.dat'
 print*,filenks
  write(filenkp,"(A16,A1,A4)") '../output/refkp_', para, '.dat'
 print*,filenkp
 call input(v,za)
 print*,'Which variable shall be varied (w or k)?'
 read(*,*) answer
 if (answer == 'w') then
  print*,'At which k point? (in kpDrude)'
  read(*,*) k
  k = k*wp1
  print*,'Start value of w? (in wpDrude)'
  read(*,*) sta
  print*,'Stop  value of w? (in wpDrude)'
  read(*,*) sto
  print*,'How many points? (integer)'
  read(*,*) maxi
  spac = (sto/sta)**(1.D0/maxi)
 
 open(unit=11,file=filenws,status='new',action='write',&
     form='formatted')
     
  open(unit=12,file=filenwp,status='new',action='write',&
      form='formatted')
 ! At this k-point we expect a significant contribution    
!  k = (2*za)**(-1)

  do i=0,maxi
   w = wp1*sta*spac**i

   slab = .FALSE.
   bulk  =.FALSE.
   toplog=.TRUE.
   call refl(w,k,resDLs,resDLp)
   toplog=.FALSE.
   call refl(w,k,resLDs,resLDp)
   EMA=.TRUE.
   call refl(w,k,resEMAs,resEMAp)
   EMA=.FALSE.
   bulk = .TRUE.
   toplog=.TRUE.
   call refl(w,k,resDs,resDp)
   slab = .TRUE.
   call refl(w,k,resDslabs,resDslabp)
   slab = .FALSE.
   toplog=.FALSE.
   call refl(w,k,resLs,resLp)


   write(11,'(13E40.20)') w/wp1, resDLs, resLDs,resEMAs,resDs,resDslabs, resLs
   write(12,'(13E40.20)') w/wp1, resDLp, resLDp,resEMAp,resDp,resDslabp, resLp   

  enddo
 close(11)
  close(12)
 elseif (answer == 'k') then
  
  print*,'At which w point? (in wpDrude)'
  read(*,*) w
  w = w*wp1
  print*,'Start value of k? (in kpDrude)'
  read(*,*) sta
  print*,'Stop  value of k? (in kpDrude)'
  read(*,*) sto
  print*,'How many points? (integer)'
  read(*,*) maxi
  spac = (sto/sta)**(1.D0/maxi)
 
 open(unit=13,file=filenks,status='new',action='write',&
     form='formatted')
  open(unit=14,file=filenkp,status='new',action='write',&
     form='formatted')
     
     
    epM = DCMPLX(1D0) - wp1**2/DCMPLX(w**2,w*g1)    
  print*,'ksp=',DREAL(w*SQRT(epm/(epm+DCMPLX(1D0))))/wp1
  do i=0,maxi
  k = wp1*sta*spac**i
   
   slab = .FALSE.
   bulk  =.FALSE.
   toplog=.TRUE.
   call refl(w,k,resDLs,resDLp)
   toplog=.FALSE.
   call refl(w,k,resLDs,resLDp)
   EMA=.TRUE.
   call refl(w,k,resEMAs,resEMAp)
   EMA=.FALSE.
   bulk = .TRUE.
   toplog=.TRUE.
   call refl(w,k,resDs,resDp)
   slab = .TRUE.
   call refl(w,k,resDslabs,resDslabp)
   slab = .FALSE.
   toplog=.FALSE.
   call refl(w,k,resLs,resLp)


   write(13,'(13E40.20)') k/wp1, resDLs, resLDs,resEMAs,resDs, resDslabs, resLs
   write(14,'(13E40.20)') k/wp1, resDLp, resLDp,resEMAp,resDp, resDslabp, resLp 
  enddo
 
  close(13)
  close(14)
  else
  print*,'No such variable...only w or k'
  endif
end program
!=================================================================================
! SUBROUTINES:
! ============
subroutine input(v,za)
  use comref
  implicit none
  ! Input related variables
  character(len=100) :: buffer, label
  integer :: pos
  integer, parameter :: fh = 15
  integer :: ios = 0
  integer :: line = 0
  real(8) :: v, za
  real(8), external :: metricNatSI, metricSInat
  character(23) :: infile
  write(infile,"(A18,A1,A4)") '../input/parameter',para,'.txt'
  print*,infile
  open(fh, file=infile)
  

  
  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.

 do while (ios == 0)
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
        pos = scan(buffer, '    ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        select case (label)
        case ('v')
           read(buffer, *, iostat=ios) v
           print *, 'Read v  : ', v
        case ('za')
           read(buffer, *, iostat=ios) za
           print *, 'Read za : ', za      
        case ('near')
           read(buffer, *, iostat=ios) near
           print *, 'Read near : ', near
        case ('slab')
           read(buffer, *, iostat=ios) slab
           print *, 'Read slab : ', slab
        case ('einf1')
           read(buffer, *, iostat=ios) einf1
           print *, 'Read einf1 : ', einf1
        case ('wp1')
           read(buffer, *, iostat=ios) wp1
           print *, 'Read wp1: ', wp1
        case ('w01')
           read(buffer, *, iostat=ios) w01
           print *, 'Read w01: ', w01
        case ('g1')
           read(buffer, *, iostat=ios) g1
           print *, 'Read g1 : ', g1
        case ('einf2')
           read(buffer, *, iostat=ios) einf2
           print *, 'Read einf2 : ', einf2
        case ('wp2')
           read(buffer, *, iostat=ios) wp2
           print *, 'Read wp2: ', wp2
        case ('w02')
           read(buffer, *, iostat=ios) w02
           print *, 'Read w02: ', w02
        case ('g2')
           read(buffer, *, iostat=ios) g2
           print *, 'Read g2 : ', g2
        case ('wp3')
           read(buffer, *, iostat=ios) wp3
           print *, 'Read wp3: ', wp3
        case ('w03')
           read(buffer, *, iostat=ios) w03
           print *, 'Read w03: ', w03
        case ('g3')
           read(buffer, *, iostat=ios) g3
           print *, 'Read g3 : ', g3
        case ('wp4')
           read(buffer, *, iostat=ios) wp4
           print *, 'Read wp4: ', wp4
        case ('w04')
           read(buffer, *, iostat=ios) w04
           print *, 'Read w04: ', w04
        case ('g4')
           read(buffer, *, iostat=ios) g4
           print *, 'Read g4 : ', g4
        case ('doped')
           read(buffer, *, iostat=ios) doped
           print *, 'Read doped: ', doped
        case ('wpd')
           read(buffer, *, iostat=ios) wpd
           print *, 'Read wpd: ', wpd
        case ('gd')
           read(buffer, *, iostat=ios) gd
           print *, 'Read gd: ', gd
        case ('dd')
           read(buffer, *, iostat=ios) dd
           print *, 'Read dd : ', dd
        case ('rho')
           read(buffer, *, iostat=ios) rho
           print *, 'Read rho: ', rho
        case ('bulk')
           read(buffer, *, iostat=ios) bulk
           print *, 'Read bulk: ', bulk
        case ('EMA')
           read(buffer, *, iostat=ios) EMA
           print *, 'Read EMA: ', EMA
        case ('dphi')
           read(buffer, *, iostat=ios) dphi
           print *, 'Read dphi : ', dphi
        case ('dtheta')
           read(buffer, *, iostat=ios) dtheta
           print *, 'Read dtheta : ', dtheta
        case default
           print *, 'Skipping invalid label at line', line
        end select
     end if
  end do
 
  wsp1  = wp1/(DSQRT(2.D0))
  pi    = 4.D0*DATAN(1.D0)
  za    = metricSInat(za,'l')
  dd    = metricSInat(dd,'l')
  dphi  = dphi/180D0  * pi
  dtheta= dtheta/180D0  * pi
  dx2    = (DCOS(dphi)*DSIN(dtheta))**2
  dy2    = (DSIN(dphi)*DSIN(dtheta))**2
  dz2    = DCOS(dtheta)**2
end subroutine input

subroutine refl(w,k,rs,rp)
use comref
implicit none
 real(8) :: w,k, kap
 complex*16 :: epI, epM, epMw, epsx, epsz
 complex*16 :: ZZ, trhmM22, ZIZM, ZIZ0,ZMZ0, ZMZI
 complex*16 :: ZpZ0, ZsZ0
 complex*16 :: kM, kI, tI, tM, tt, cc, keffp,keffs
 complex*16 :: sq, kapc
 complex*16 :: trh, exp1, exp2, expr
 complex*16 :: cI,cM, M22, M21Z0
 logical :: dec(3)
 complex*16 :: epDL, epD
 complex*16 :: one, i
 real(8) :: abe
 complex*16 :: rp, rs, ref, phase

  one = DCMPLX(1D0)
  i = DCMPLX(0D0,1D0)
  kap = DSIGN(DSQRT(DABS(k**2-w**2)),k**2-w**2)
 
  if (kap < 0D0) then
   kapc = DCMPLX(0D0,kap)

  else
   kapc = DCMPLX(kap, 0D0)
  endif

  epI = DCMPLX(einf2) - wp2**2/DCMPLX(w**2 - w02**2,g2*w) &
      - wp3**2/DCMPLX(w**2 - w03**2,g3*w) &
      - wp4**2/DCMPLX(w**2 - w04**2,g4*w)

  if (doped) then ! For doping we add an additional
                  ! Drude background :
   epI = epI + DCMPLX(-(wpd)/(w**2+gd**2),&
                   wpd**2*gd/(w*(w**2+gd**2))) 
  endif
 
 
  epM = DCMPLX(einf1) - wp1**2/DCMPLX(w**2-w01**2,w*g1) 
  epMw = DCMPLX(einf1)*w - wp1**2/DCMPLX(w-(w01**2)/w,g1) 
   
   ! Near-field approximation:
   if (near) then
    kI = i*k
    kM = i*k
    kapc = DCMPLX(k)
   
   ! Including retardation effects:
   else
    kI = SQRT(epI*w**2 - DCMPLX((k)**2))
    ! Choose the correct sign of the wavevector
    if ((DIMAG(kI)<0D0)) then
    kI = -kI
    endif
   
   
    kM = SQRT(epMw*w - DCMPLX((k)**2))
    ! Choose the correct sign of the wavevector
    if ((DIMAG(kM)<0D0)) then
    kM = -kM
    endif

endif
 
!  
    cM = COS(kM*rho*dd)
    cI = COS(kI*(1.D0-rho)*dd)
!    
    tM = TAN(kM*rho*dd)
    tI = TAN(kI*(1.D0-rho)*dd)

!    
 ! Now we add the part of the structure of the Green's function
 ! where we already apply the trace in order to reduce the
 ! tensorial structure to a scalar.
   if (bulk) then
    if (toplog) then
     ZpZ0 = kM/(epM*i*kapc)
     ZsZ0 = i*kapc/kM
    else
     ZpZ0 = kI/(epI*i*kapc)
     ZsZ0 = i*kapc/kI
    endif
   else
    if (EMA) then
     epsx = epM*rho + epI*(1D0-rho)
     epsz = 1D0/(rho/epM + (1D0-rho)/epI)
     
     ! near field
     if (near) then
     keffp = i*SQRT(epsx/epsz)*k
     keffs = i*k
     
     else
     keffp = SQRT(epsx*w**2 - (epsx/epsz)*k**2) 
     keffs = SQRT(epsx*w**2 - DCMPLX(k**2)) 
     
     if ((DIMAG(keffp)<0D0)) then
      keffp = -keffp
     endif
     

     if ((DIMAG(keffs)<0D0)) then
      keffs = -keffs
     endif
     endif
    else
   ZpZ0 = -refZ(.TRUE.)!refZ(.TRUE.)
   ZsZ0 = refZ(.FALSE.)!-refZ(.FALSE.)
!    print*,'ZpZ0sup=',ZpZ0
!    print*,'ZsZ0sup=',ZsZ0
!    print*,'ZpZ0bul=',kM/(epM*i*kapc)
!    print*,'ZsZ0bul=',i*kapc/kM
    endif
   endif
   if (EMA) then
   rp = (i*kapc*epsx - keffp)/((i*kapc*epsx + keffp))
   rs = (i*kapc - keffs)/((i*kapc + keffs))
   else
! Now we calculate the p-polarized reflection coefficient...
   rp =  (one  -  ZpZ0 ) / ( one + ZpZ0 )
! ... and the s-polarized reflection coefficietn.
   rs = -(one  - ZsZ0 ) / ( one + ZsZ0 )
   
   ! Now we are able to calculate thin slabs of the bulk
   ! material. Here we implemented the M-material.

   endif
   if (slab) then
    rp = rp*(one - Exp(2*i*kM*dd*rho))&
           /(one - rp**2*Exp(2*i*kM*dd*rho))
    rs = rs*(one - Exp(2*i*kM*dd*rho))&
           /(one - rs**2*Exp(2*i*kM*dd*rho))

   endif
   
   
 contains
  function refZ(sig) 
   implicit none
   complex*16 :: refZ
   logical :: sig ! If TRUE, the p-polarization is calculated,
                ! if FALSE, thes s-polarization is calculated
   ! The s- and p-polarization differ in some definitions:
   if (sig) then
    ZIZM = (kI*epM)/(kM*epI)
    ZMZI = 1D0/ZIZM
!     ZMZI = (kM*epI)/(kI*epM)
    ZIZ0 = kI/(i*epI*kapc)
    ZMZ0 = kM/(i*epM*kapc)
    M21Z0 = i*(tM/ZMZ0  + tI/ZIZ0) 
   elseif (.not.sig) then
    ZIZM = (kM)/(kI)
    ZMZI = 1D0/ZIZM
!     ZMZI = (kI)/(kM)
    ZIZ0 = (i*kapc)/kI

    ZMZ0 = (i*kapc)/kM
    
    M21Z0 = -i*(tM/ZMZ0  + tI/ZIZ0) 
   endif

    ZZ = (ZIZM + ZMZI)
    tt = tI*tM
    
    trh = one-0.5D0*ZZ*tt

    ! The only term which differs when changing 
    ! the layering order is M22:
    if (toplog) then
    trhmM22 = 0.5D0*(-ZMZI+ZIZM)*tt
    else
    trhmM22 = 0.5D0*(ZMZI-ZIZM)*tt
    endif
 
   cc = cI*cM
   
   sq =i*SQRT(ZZ*tt - (0.5D0*ZZ*tt)**2 &
            + (tI**2 + tM**2 + (tt)**2))
   
   exp1 = trh + sq
   exp2 = trh - sq

   dec(1) = (dlog(abs(exp1))<=-(log(abs(cI))+log(abs(cM))))
   dec(2) = (dlog(abs(exp1))<=-(log(abs(cI))+log(abs(cM))))
   dec(3) = abs(exp1)==abs(exp2)
   
   if (dec(1).and..not.dec(2)) then
    sq = sq
   elseif (dec(2).and..not.dec(1)) then
    sq = -sq
   elseif (dec(3)) then
    if     (DIMAG(exp1*cc) > 0D0) then
     sq = sq
    elseif (DIMAG(exp2*cc) > 0D0) then
     sq = -sq
    else
     sq = sq
    endif
    else
     if (abs(exp1)<abs(exp2)) then
      sq = sq
     elseif (abs(exp2)<abs(exp1)) then
      sq = -sq
!      else
!      print*,'Niiiichts funktioniert!'
!      stop
     endif
   endif
    refZ = (trhmM22 - sq)/M21Z0
!     print*,'sqrt=',sq
   end function refZ
   
   
end subroutine refl



function metricSInat(SI, lf)
implicit none
real*8 :: metricSInat ! output in natural units
real*8 :: SI     ! SI-unit
character(1) :: lf  ! Either 'l' for length
                    ! or 'f' for frequency
if ( lf == 'l' ) then
 metricSInat = SI/(1.97412D-7)
elseif ( lf == 'f' ) then
 metricSInat = SI*6.58D-16
else
 print*,'Either l for length or f for frequency'
 stop
endif
return
end function

function metricNatSI(nat, lf)
implicit none
real*8 :: metricNatSI ! output in SI units
real*8 :: nat       ! natural unit
character(1) :: lf  ! Either 'l' for length
                    ! or 'f' for frequency
if ( lf == 'l' ) then
 metricNatSI = nat*(1.97412D-7)
elseif ( lf == 'f' ) then
 metricNatSI = nat/6.58D-16
else
 print*,'Either l for length or f for frequency'
 stop
endif
return
end function