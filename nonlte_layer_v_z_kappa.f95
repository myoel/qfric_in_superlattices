!_______________________________________________________________________________!
!                                                                               !
! This program calculates the quantum friction without the LTE in a system      !
! where a harmonic oscillator flies over a semi-infinite multi-layer material   !
! composed out of alternating layers with \epsilon_I(\omega,k) and              !
! \epsilon_M(\omega,k) (see SUBROUTINES).                    !
! The scattering at the interfaces is described by Bloch waves and does not     !
! use the effective medium approximation (EMA).                                 !
! We work in natural units, i.e.:                                               !
! c = 1                                                                         !
! 4*pi*epsilon_0 = 1                                                            !
! hbar = 1                                                                      
!
!
! For a better visualization F_fric(z) is given as F_fric(z)*z^(10)
!
!_______________________________________________________________________________!

! In this module all shared parameters are defined:
module com
 implicit none
 ! Local thermal equilibrium:
 logical :: LTE  ! if TRUE, then only the
                 ! LTE calculation is perfomed
 
 ! Nano-particle or harmonic oscillator:
 logical :: nano ! if TRUE, then the calculation is done
                 ! for a nano-particle, meaning that
                 ! gamma=const. and delta=const.
 
 real(8), public :: Fcons ! collected constants drawn out of the Integration
 
end module com

program nonltebloch
 use gen
 use com
 implicit none
 Integer :: i, maxi
 real(8) :: sta, sto, spac ! Start and stop and spacing of the calculated points
 real(8) :: Ffric ! calculated quantum friction
 real(8) :: mures
 real(8) ::lowb1, lowb2, lowsup, resb1, resb2, ressup
 real :: t1,t2 ! variables to measure the CPU-time
 character(1) :: vari
 character(29):: filen
 
 ! Create File Name:
 print*,'Ffric(v;za=const.) or Ffric(v=const.,za) ? ( v / z )'
 read(*,*) vari
 print*,'Which input file ? (1,2,3,...)'
 read(*,*) para
 write(filen,"(A23,A1,A1,A4)") '../output/NonLTE_layer_', vari, para, '.dat'
 print*,filen
 ! Generating the new file:
  open(unit=11,file=filen,status='new',&
  action='write', form='formatted')
  close(11)
 ! Choose start and stop value and spacing
 print*,'Start value? (in units of c for v or meters for za)'
 read(*,*) sta
 print*,'Stop  value? (in units of c for v or meters for za)'
 read(*,*) sto
 print*,'How many points? (integer)'
 read(*,*) maxi
 spac = (sto/sta)**(1.D0/maxi)
  
 ! Load input data:
 call input
 
 ! Variation of velocity:
 ! ======================
 if ( vari == 'v' ) then 

  ! CALCULATIONS:
  ! =============
  do i=0,maxi
  
   ! Timestamp 1
   call CPU_TIME(t1)
   
   ! Open file for writing:
   open(unit=11,file=filen,status='old',action='write',&
     form='formatted',position='append') 
   
   ! For equal spacing in a logarithmic plot we space
   ! the velocity with an exponent:
    print*,i*100.D0/maxi,'%'
    v = sta*spac**i
    print*,'v= ',v
   
   ! Collecting drawn out constants:
    Fcons = 2*v/pi * (2*a0/(pi*(2*za)**4))**2
      
   ! Call the full calculation:
   call Ffricf(Ffric)
   call Fasym(lowb1, lowb2, lowsup, resb1, resb2, ressup)
   ! Adding drawn out constants and normalization:
    Ffric = Ffric * Fcons /F0
   
   ! OUTPUT:
   ! =======
   
   ! Display on screen:
   print*,'Ffric=',Ffric

   ! Timestamp 2
   call CPU_TIME(t2)
   
   print*,'Time:',t2-t1,'s'
   
   ! Write in file:
   write(11,'(8E40.20)') v, Ffric, lowb1, lowb2, lowsup,&
                           resb1, resb2, ressup
   close(11)
  enddo
  
 ! Variation of distance:
 ! ====================== 
  elseif (vari == 'z') then
   ! Converting to natural units
   call convunits(sta,'l',.FALSE.)
   call convunits(sto,'l',.FALSE.)
 
   do i=0,maxi
    ! Timestamp 1
    call CPU_TIME(t1)
    
    open(unit=11,file=filen,status='old',action='write',&
       form='formatted',position='append') 
  
    print*,i*100.D0/maxi,'%'
    
    za = sta*spac**i
    
   ! Collecting drawn out constants:
    Fcons =2* v/pi * (2*a0/(pi*(2*za)**4))**2
      
   ! Call the full calculation:
    call Ffricf(Ffric)
   call Fasym(lowb1, lowb2, lowsup, resb1, resb2, ressup)
   ! Adding drawn out constants and normalization:
    Ffric = Ffric * Fcons /F0
   
   ! OUTPUT:
   ! =======
   
   ! Display on screen:
   !call convunits(za,'l',.TRUE.)
   
   print*,'za',za/(2*pi/wp1),'in lambda_p(Drude)'
   print*,'Ffric=',Ffric

   ! Timestamp 2
   call CPU_TIME(t2)
   print*,'Time:',t2-t1,'s'

   write(11,'(8E40.20)') za/dd, Ffric, lowb1, lowb2, lowsup, resb1, resb2, ressup
   close(11)

   enddo
else
 print*,'Wrong input.'
endif

end program
!=================================================================================
! SUBROUTINES:
! ============
! This subroutine calculates all integrals needed for QF
! and yields by combining them the final result

subroutine Ffricf(res)
use gen
use intdef
implicit none
 real(8) :: res
 real(8) :: xa,xp
 xa = wa*(2*za)/v
 
if ( xa < kcut ) then
 res = integ3(integrand,0D0,xa) &
      +integ3(integrand,xa,kcut)
else 
 res = integ3(integrand,0D0,kcut)
endif


 contains
! Here we use wrapper to hand over a integration
! as function to another integration. Before the
! last integration above, we collect all needed
! subintegrals and arrange them in the correct
! form within the function integrand.

function integrand(xx)
 implicit none
 real(8), intent(in) :: xx
 real(8) :: integrand
 real(8) :: wgam, delta, alpha2
 real(8) :: muk, muq, mukk, muqq, mu
 complex*16 :: alpha
  real(8) :: phi,kap,xi
  COMMON /phikapxi/ phi,kap,xi
 
 xi = xx
 ! Firstly, we need the four integrals
 ! int dk^2 (..) * int dq^2 (..) *(qx-kx),
 ! where (..) is the common integrand.
 ! Since the q and k integration have 
 ! different limits we have to calculated
 ! each single integral:
  
 rimag   =.TRUE. 
 
 ! int dk^2 (..) =
 addkx   =.FALSE.
 heavi   =.TRUE.
 doppsign=.TRUE.
 muk = wrapxi(xi)
 
  ! int dk^2 (..)*kx
 addkx   =.TRUE.
 heavi   =.TRUE.
 doppsign=.TRUE.
 mukk = wrapxi(xi) 

 ! int dq^2 (..) =
 addkx   =.FALSE.
 heavi   =.TRUE.
 doppsign=.FALSE.
 muq = wrapxi(xi)

 ! int dq^2 (..)*qx
 addkx   =.TRUE.
 heavi   =.TRUE.
 doppsign=.FALSE.
 muqq = wrapxi(xi) 
 
 ! Secondly, we like to calculate the
 ! correction of the polarizabillity
 ! w*gamma as imaginary part and
 ! delta as real part:
 rimag =.FALSE.
 heavi =.FAlSE.
 addkx =.FALSE.
  delta = wrapxi(xi)*a0*2D0/( pi*(2*za)**3 )
 ! The w*gamma can be obtained by 
 ! adding up already calculated
 ! contributions:
  wgam = (muk+muq)*a0*2D0/( pi*(2*za)**3 )
 
 ! With this we can calculate the
 ! polarizabillity:
 
 alpha2 = 1D0/((1.D0 - delta - (xi/xa)**2)**2 + wgam**2)
 
 ! And the contributions from the 
 ! truncated Green's functions:
  mu = muk*muqq - mukk*muq

 integrand = mu *alpha2
 
end function integrand


! Second integration over k :
 function wrapxi(xx)
  implicit none
  real(8), intent(in) :: xx
  real(8) :: wrapxi
  real(8) :: limk1, limk2
  real(8) :: phi,kap,xi
  real(8) :: int1, int2
  COMMON /phikapxi/ phi,kap,xi
  xi = xx

  limk1 = 0D0
  limk2 = kcut

  ! In the analytical calcualtions
  ! kcut = inf, but in order to avoid
  ! inf we introduce a cut-off kcut>>1.
  ! The cut-off can be justified since
  ! the exp(-k) is strongly decaying.
  
  ! For the particular case of the
  ! int_{-inf}^{-xi}dkx integration
  ! the lower limit of the k integral
  ! changes to xi
  if (heavi.and..not.doppsign) then
   limk1 = xi
  endif
  
  
   int1 = integ2(wrapkap,limk1,limk2)
  if (heavi.and..not.doppsign) then
   wrapxi = int1
 else
  int2 = integ2(wrapkap,-xi*v,0D0)
  wrapxi = int1 + int2
 endif
  
  end function
  
! First integration over cos(phi) : 
 function wrapkap(kk)
  implicit none
  real(8), intent(in) :: kk
  real(8) :: wrapkap
  real(8) :: limc1, limc2
  real(8) :: phi,kap,xi
  COMMON /phikapxi/ phi,kap,xi
  kap = kk  
  ! For the full integration we use the
  ! limits :
  limc1 = -1D0
  limc2 =  1D0
  ! If the integral is truncated via a
  ! Heaviside function we obtain:
  if (heavi) then
   if (doppsign) then
    ! This case is for the int_{-xi}^inf dkx
    if (kap >= 0D0) then
    limc1 = maxval((/-1D0,-xi/kap/))
    else
    limc1 = -1D0
    endif
   else
    ! This case is for the int_{-inf}^{-xi} dkx
    limc2 = -xi/kap
   endif
  endif
  ! Performing the integration over cos(phi)
  wrapkap = integ1(wrapcos,DACOS(limc2),DACOS(limc1))

  end function
  
! Input function :
 function wrapcos(pp)
  implicit none
  real(8), intent(in) :: pp
  real(8) :: cosp
  real(8) :: vsin
  real(8) :: wrapcos
  real(8) :: phi,k,xi,kap
  COMMON /phikapxi/ phi,kap,xi
  phi = pp
  cosp = DCOS(phi)
  vsin = 1D0-(cosp*v)**2
  k = (cosp*v**2*xi &
   + DSQRT(dsign(kap**2,kap)*vsin&
    +(v*xi)**2 ))/vsin
     
  wrapcos =refl((k*cosp + xi)*v/(2*za),kap/(2*za),k/(2*za),cosp)&
           *(2*za)**2*( 1D0 + cosp*xi*v**2&
           /DSQRT(vsin*dsign(kap**2,kap)+(xi*v)**2))/vsin

  ! For some integrals we need to add a kx.
  ! Therefore we add here the cos(phi) and
  ! in the k integration a k.
  
  if (addkx) then
   wrapcos = wrapcos*cosp*k
  endif
 end function 
 
end subroutine Ffricf