! This module contains general parameters,
! the input routine and a converter from 
! natural units to SI units and vice versa.

module gen
implicit none
 ! General parameters:
 real(8) :: pi ! Pi
 real(8) :: F0 ! Scaling constant for QF
 real(8) :: v, za ! Velocity and distance to the surface
                  ! of the atom
                  
 real(8), public :: resis1 ! The low-frequency resistivity
                           ! of the metal (times the dielectric
                           !  constant)
 real(8), public :: resis2 ! The low-frequency resistivity
                           ! of the dielectric (times the
                           !  dielectric constant)  
                           
 ! Defining the superlattice:
 real(8), public :: dd  ! Thickness of the superlattice
 real(8), public :: rho ! filling factor rho= dm/(dm+di),
                        ! with d1/2 is the thickness of the 
                        ! components
                        
 ! Defining the metal component (Drude model):

 ! epD(w) = 1 - wp^2/(w^2 + i*g*w)
 real(8), public :: einf1! Background constant
 real(8):: wp1 ! Plasma frequency
 real(8):: wsp1! Surface Plasmon freqeuncy
 real(8):: g1 ! Damping constant

  real(8):: wsp! Surface Plasmon freqeuncy
 ! Defining the insulator component (Lorentz model):
 
 ! epL(w) = 1 - wp^2/(w^2 - w0^2 + i*g*w)
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
 
 ! Doping which yields in Drude-Lorentz model:
 
 ! epL(w) = 1 - wp^2/(w^2 - w0^2 + i*g*w) + (epDd(w)-1)
 
 logical:: doped ! If TRUE then Drude term is added
 real(8):: wpd ! Plasma freqeuncy
 real(8):: gd  ! Damping constant
 
  ! Defining the order of the layers:
 character(1) :: top ! if TRUE metal on top if false the
                     ! other material is on top
 logical :: toplog ! Converts 'top' to logical
 logical:: bulk ! if TRUE then only the bulk of the top-most
                 ! material is calculated.
 logical :: EMA  ! if TRUE then only EMA is performed
 
  ! Flags to perform different integral with the same
 ! integration scheme :
 logical:: heavi ! Add Heaviside function if TRUE
 logical:: doppsign ! sign of the Dopplershift 
 logical:: addkx ! multiplying with kx if TRUE
 logical :: rimag! Chooses either imaginary or real
                 ! part of the reflection coefficient
 
 
 ! Local thermal equilibrium:
 logical :: LTE  ! if TRUE, then only the
                 ! LTE calculation is perfomed
 
 ! Nano-particle or harmonic oscillator:
 logical :: nano ! if TRUE, then the calculation is done
                 ! for a nano-particle, meaning that
                 ! gamma=const. and delta=const.
 
 ! Defining the atomic dipole parameters:
 real(8), public :: a0 ! static polarizability
 real(8), public :: dphi, dtheta ! direction of the dipole
 real(8), public :: dx2, dy2, dz2 ! Cartesian form of dipole
 logical, public :: average ! If TRUE the average over all possible
                            ! angles of the dipole is taken
 real(8), public :: wa ! dipole frequency
 real(8), public :: ga ! internal damping
 real(8), public :: da ! frequency shift

 ! Cut-off for the integrations
 real(8), public :: kcut ! Cut-off due to Exp(-kcut) ~ 0
  
 ! File names:
 character(1), public :: para
 
   logical:: near ! If TRUE then only near-field approximation 
                 ! is considered.
  logical::slab  ! If TRUE then only a slab with thickness d*rho is calculated
 
 contains
!====================================================
! Following subroutines and function are contained:
! - input (reads input file)
! - convunits (converts between SI and natural units)
! - refl (calculates reflection coefficient)
! - Fasym (calculates asymptotics)
!====================================================

!----------------------------------------------------
! Input routine 
!----------------------------------------------------
! This subroutine imports all relevant parameters
! to the main program
 subroutine input
  implicit none
  ! Input related variables
  character(len=100) :: buffer, label
  integer :: pos
  integer, parameter :: fh = 15
  integer :: ios = 0
  integer :: line = 0
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
        case ('LTE')
           read(buffer, *, iostat=ios) LTE
           print *, 'Read LTE  : ', LTE
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
        case ('top')
           read(buffer, *, iostat=ios) top
           print *, 'Read top: ', top
        case ('bulk')
           read(buffer, *, iostat=ios) bulk
           print *, 'Read bulk: ', bulk
        case ('EMA')
           read(buffer, *, iostat=ios) EMA
           print *, 'Read EMA: ', EMA
        case ('nano')
           read(buffer, *, iostat=ios) nano
           print *, 'Read nano : ', nano
        case ('wa')
           read(buffer, *, iostat=ios) wa
           print *, 'Read wa : ', wa
        case ('ga')
           read(buffer, *, iostat=ios) ga
           print *, 'Read ga : ', ga
        case ('da')
           read(buffer, *, iostat=ios) da
           print *, 'Read da : ', da
        case ('a0')
           read(buffer, *, iostat=ios) a0
           print *, 'Read a0 : ', a0
        case ('average')
           read(buffer, *, iostat=ios) average
           print *, 'Read average: ', average
        case ('dphi')
           read(buffer, *, iostat=ios) dphi
           print *, 'Read dphi : ', dphi
        case ('dtheta')
           read(buffer, *, iostat=ios) dtheta
           print *, 'Read dtheta : ', dtheta
        case ('kcut')
           read(buffer, *, iostat=ios) kcut
           print *, 'Read kcut: ', kcut
        case default
           print *, 'Skipping invalid label at line', line
        end select
     end if
  end do
  wsp  = wp1/(DSQRT(einf1+einf2))
  wsp1  = wp1/(DSQRT(2D0))
  pi    = 4.D0*DATAN(1.D0)
  F0    = 3*wsp1**5*a0*2
  call convunits(za,'l',.FALSE.)
  call convunits(dd,'l',.FALSE.)
  dphi  = dphi/180D0  * pi
  dtheta= dtheta/180D0  * pi
  dx2    = (DCOS(dphi)*DSIN(dtheta))**2
  dy2    = (DSIN(dphi)*DSIN(dtheta))**2
  dz2    = DCOS(dtheta)**2
  resis1 =  g1/(wp1**2)
  resis2 = g2*(wp2/((einf2+1D0)*w02**2 + wp2**2))**2
  if (doped) then
  resis2 = gd/(wpd**2)
  endif
   if (top == 'M') then
    toplog =.TRUE.
   elseif (top == 'I') then
    toplog =.FALSE.
   else
    print*,'No correct ordering!'
    stop
   endif
end subroutine input
!====================================================


!----------------------------------------------------
! Unit Converter (nat -> SI or SI -> nat)
! ---------------------------------------------------
subroutine convunits(conv,lf,natSI)
 implicit none
 real(8) :: conv        ! variable which will be converted
 logical:: natSI 
 character(1) :: lf  ! Either 'l' for length
                    ! or 'f' for frequency
 if (natSI) then
 ! Converting natural units to SI system

  if ( lf == 'l' ) then
   conv = conv*(1.9732705D-7)
  elseif ( lf == 'f' ) then
   conv = conv/6.5821220D-16
  else
   print*,'Either l for length or f for frequency'
   stop
  endif
 
 else
 ! Converting from SI system to natural units

  if ( lf == 'l' ) then
   conv = conv/(1.9732705D-7)
  elseif ( lf == 'f' ) then
   conv = conv*6.5821220D-16
  else
   print*,'Either l for length or f for frequency'
   stop
  endif
 endif

return
end subroutine convunits
!====================================================

!----------------------------------------------------
! Calculation of the reflection coefficient
! ---------------------------------------------------
function refl(w,kap,k,cosp)
implicit none
 real(8) :: w,k, cosp, kap
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
 real(8) :: refl
  one = DCMPLX(1D0)
  i = DCMPLX(0D0,1D0)
 
 
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
 
 
  epM = DCMPLX(einf1) - wp1**2/DCMPLX(w**2,w*g1) 
  epMw = DCMPLX(einf1)*w - wp1**2/DCMPLX(w,g1) 
   
   
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
   
! All together we construct the Green's tensor
! traced with the static polarizability tensor
   ref = (rp*&
         ( (dz2*k**2) + dsign(kap**2,kap)*(dx2*cosp**2 + dy2*(1D0-cosp**2)))&
        + (w**2)*rs*&
        ( dx2*cosp**2 + dy2*(1D0-cosp**2) ))*1.5D0*EXP(-2*kapc*za)
        
     if (abs(ref)/=abs(ref)) then
     print*,'exp1=',exp1,abs(exp1)
    print*,'exp2=',exp2,abs(exp2)
       print*,'w=',w
    print*,'k=',k
    print*,'trh=',trh
    print*,'cI=',cI
    print*,'cM=',cM
    print*,'kI=',kI
    print*,'kM=',kM
    print*,'ZIZ0=',ZIZ0
    print*,'ZMZ0=',ZMZ0
    print*,'sq=',sq
    print*,'ref=',ref
    print*,'rp=',rp
    print*,'rs=',rs
    print*,'einf'
     stop
    endif

 
 if (rimag) then
  refl = DIMAG(ref)
 elseif (.not.rimag) then
  refl = DREAL(ref)
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
   
end function refl
!====================================================


!----------------------------------------------------
! Calculation of asymptotics
! ---------------------------------------------------
! This subroutine calculates all needed asymptotics
! Dipole dependent asymptotics:
! (low velocity: bulk1, bulk2, sup)
! (resonant    : bulk1, bulk2, sup)
subroutine Fasym(lowb1, lowb2, lowsup, resb1, resb2, ressup)
implicit none
real(8) :: lowb1, lowb2, lowsup
real(8) :: resb1, resb2, ressup
real(8) :: FLTE, FJ
real(8) :: Fbulk, Fres
real(8) :: sum1, sum2
real(8) :: x
real(8) :: Aa0, Aa2, ALTE, AJ
integer :: n

! Calculating the low-velocity asymptotics for the bulks:
    Fbulk  = (27D0/64D0*( &
             27D0 + (11*DCOS(2*dphi) -27D0)*DSIN(dtheta)**2 &
            +(13D0-14*DCOS(dphi)**2+3*DCOS(dphi)**4)*DSIN(dtheta)**4)) &
            *32/(pi**3)*a0**2*(4*pi)**2*v**3/(2*za)**10

    lowb1 = Fbulk*resis1**2
    lowb2 = Fbulk*resis2**2
    
! Calculating asymptotic for resonant atom for bulk materials:   
   Fres = 0.5D0*(3D0-3*DSin(dtheta)**2*DSin(dphi)**2) &
          *a0/(2*DSQRT(pi))*dEXP(-2*za*wa/v) 
          
   resb1 = Fres*(1D0+15*v/(4*za*wa))*DSQRT(wa**7/(za**5*v**3))*resis1
   resb2 = Fres*(1D0+15*v/(4*za*wa))*DSQRT(wa**7/(za**5*v**3))*resis2

! Calculating asymptotic for resonant atom for superlattice:   
   ressup = Fres*DSQRT(wa**5/(za**5*v))*(1D0+5*v/(2*za*wa))/(rho*dd)

! Calculating the low-velocity asymptotic for superlattice:
   x = dd*rho/za
   sum1 = 1D0
   sum2 = 1D0
   FJ   = 1D0
   Aa0  = 3D0/2D0 * &
       (1D0 + (3*DCOS(dphi)**2 -2D0 )*DSIN(dtheta)**2)
   Aa2  = 3D0/2D0 * (1D0-DCOS(dphi)**2*DSIN(dtheta)**2)
   ALTE = (Aa0 + 3*Aa2)*(5*Aa0 + 7*Aa2)/48D0
   AJ   = ((3*Aa0+5*Aa2)/8D0)**2
   
   do n=1,1000000
    sum1 = sum1 + 2.D0/(1D0 + n*x)**7
    sum2 = sum2 + 2.D0/(1D0 + n*x)**3
    FJ   = FJ   + 2.D0/(1D0 + n*x)**5
   enddo
   
   FLTE = 90D0*ALTE*sum1*sum2
   FJ   = 72D0*AJ*FJ**2

lowsup = FLTE + FJ
lowsup = lowsup/(pi**3)*a0**2*(g1*4*pi/(wp1**2))**2*v**3/(2*za)**10

   ! Normalization to the constant F0:
   lowb1 =lowb1 /F0
   lowb2 =lowb2 /F0
   lowsup=lowsup/F0
   resb1 =resb1 /F0
   resb2 =resb2 /F0
   ressup=ressup/F0

end subroutine Fasym
!====================================================

end module gen
