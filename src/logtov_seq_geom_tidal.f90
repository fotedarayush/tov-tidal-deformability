program tov

!Programm zum Loesen der TOv gleichung fuer verschiedene Zentraldichten

!gfortran -o logtov_seq_geom_tidal.out logtov_seq_geom_tidal.f90

!  in Schwarzschild coordinates!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   things to be set before compiling:
!
!   * which eos to use
!
!   * redo a energyshift that is possibily made in the eos
!
!   * converge towards a certain mass 
!
!   * histgram output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none


!Funktionen
!real:: rhsp
!real:: rhsm
!real:: eps


!Variablen:
real:: r,p,nue,m,rho,h,rbar,pexit,rhoexit
real:: l1,l2,l3,l4,k1,k2,k3,k4,m1,m2,m3,m4,n1,n2,n3,n4,o1, o2, o3, o4,p1,p2,p3,p4

real::  a0,y,comp,caph,beta,tidalk2,nenner,lambda,kappa2  !see Hindere et al. 2010

real:: drbar,alp

!Konstanten
real:: r0, rstar, p0,rho0,nue0,rhostep,mrest

integer:: i, n,q,k,bn,l,zae,length
integer:: strange  !entscheidet ob mit strange eos oder mit Tabelle gerechnet wird

integer:: ipoly  !entscheidet ob mit poltrope gerechnet wird

real :: nf, ef, R1,R2,LL1,LL2, lambda0
integer:: kk

!eos-Parameter
real::a,eps0,c,G,rest0,dum
!eos-arrays
  real, dimension(30000):: harray, rhoarray,parray,earray,yparray,narray,e0array
real::factorr  !Faktor je nachdem ob strange oder shen muss die ruhemassen dichte in 1/fm^3 oder cgs gerechnet werden, dadurch andere umrechnungsfakor

!real:: r12,r13,r14,m12,m13,m14 ! save the radius and the grav. mass at a density of 10^12,10^13,.. g/cm^3
!real:: rhofac=6.176e17 !factore for transformation cgs-geom units

real::   pini, restini, oldmmax, rhoini, deltarho

integer:: nmodels

!fuer Polytrope:
real::kpoly,gamma


  !declarations for the histogram
     integer, parameter::     nhist=100
     real,dimension(nhist)::  yehist,mhist,pmhist,rhohist
     real::                   dye,summ,drho
     integer::                ll

     
! read infile

open(unit=1744,file='infile',status='unknown',form='formatted')
read(1744,*) rhoini, deltarho, nmodels
close(1744)

oldmmax=0.0

kpoly=123.6 !120.!25000.!85.!110.!87.
gamma=2. !2.05 !2.85  !2.


strange=0   !entscheidet ob mit MIT-Bagmodell(analytisch) oder mit Tabellen eos gerechnet wird

ipoly=0 !das laeuft noch nicht, da ich eos calls umschreiben muss
!Achtung: wird mit polytope gerechnet entsprcht das rho einer Ruhemassendichte

if (ipoly/=1) then
if (strange/=1) call eosinit2(rhoarray,parray,harray,earray,yparray,length,narray,e0array)
endif

!Parameter fuer MITBag-Modell
! a=0.324   !dimlos
! eps0=3.0563*1e14   !/6.167e17  !geom units
! rest0=0.19611

a=0.324   !dimlos
eps0=3.0563*1e14   !/6.167e17  !geom units
rest0=0.19611

!initial constant for beta and H to initialize values

!Naturkonstanet
c=1. !2.9979e10  !in cgs
G=1. !6.674e-8

!es wird in cgs Einheiten gerechnet nur unten wird beim rausschreiben in geom umgerechnet und 
! bei den aufrufen der tabellen eos

open(unit=1745,file='logtov_seq_geom_tidal.dat',status='replace', form='formatted',recl=1024)

!Anfangsdichte
rhoini=rhoini/6.176e17
deltarho=deltarho/6.176e17
rho0=rhoini         !1.e14/6.167e17
rhostep=deltarho    !0.4e14/6.176e17
rho0=rho0-rhostep



do q=1,nmodels  !Schleife ueber verschieden Zentraldichten sollte einsbleiebne, sonst
              !macht das profile rausschreiben keinen Sinn

write(*,*)"step",q

  dye=0.0005
  drho=0.02

  do l=1,nhist
     mhist(l)=0.0
    pmhist(l)=0.0
    yehist(l)=dye*real(l)
   rhohist(l)=13.+drho*real(l)
!   rhohist(l)=1.e13+1.e13*real(l-1)
  enddo
!this was just for the histogram


!Startradius, Feinheit, Masse: Starting radius, fineness, mass:
r0=10.*1e-5/1.746  !falscherumrechungsfaktor macht hier nix ! Startradius in geom
h=10.*1e-5/1.746    !feinheit in geom
m=0.     !grav masse
mrest=0.  !ruhemasse -> rest mass
n=2000000!floor(35*1e-5/h)  !Maximal 35 km Radius - Anzahl maximaler schleifenducrhlaeufe fur RK
r=r0   
nue=0.1 ! lapse bis auf faktor, der erst nach Integration feststeht -> lapse except for a factor that is only established after integration
a0=1.
caph=a0*r*r
beta=2.*a0*r

!Startdichte -> Starting density
rho0=rho0+rhostep

!################################
!p=a*c*c*(rho0-eps0)   !ein pout(rho0,p,energ0rho0,rhoarray,parray,harray,earray, yparray)
if (ipoly==1) then
    p=kpoly*rho0**gamma !Achtung: rho ist jetst Ruhemassendichte
else
    p=pout(rho0)  !rho0 Energiedichte

    pini = p
    restini = rest(pini)

endif
!##############################



m=4.*3.1415927654*eps(p)*r0**3/3.
write(*,*)"Start",rho0,p

!p=p/6.167e17

!#->
open(unit=1744,file='profile.dat',status='replace',form='formatted',recl=1024)

write(*,*)"drin",n
do i=1,n
!write(*,*)"drin"
   !Runge-Kutta-Schritte
   !1.RK
   k1=rhsp(r,p,m)
   l1=rhsm(r,p,m)
   m1=rhsn(r,p,m)
   n1=rhsr(r,p,m)
   o1=rhsh(r,p,m,beta)
   p1=rhsb(r,p,m,beta,caph)

   !2.RK
   k2=rhsp(r+h/2.,p+h*k1/2.,m+h*l1/2.)
   l2=rhsm(r+h/2.,p+h*k1/2.,m+h*l1/2.)
   m2=rhsn(r+h/2.,p+h*k1/2.,m+h*l1/2.)
   n2=rhsr(r+h/2.,p+h*k1/2.,m+h*l1/2.)
   o2=rhsh(r+h/2.,p+h*k1/2.,m+h*l1/2.,beta+h*p1/2.)
   p2=rhsb(r+h/2.,p+h*k1/2.,m+h*l1/2.,beta+h*p1/2.,caph+h*o1/2.)
   
   !3.RK
   k3=rhsp(r+h/2.,p+h*k2/2.,m+h*l2/2.)
   l3=rhsm(r+h/2.,p+h*k2/2.,m+h*l2/2.)
   m3=rhsn(r+h/2.,p+h*k2/2.,m+h*l2/2.)
   n3=rhsr(r+h/2.,p+h*k2/2.,m+h*l2/2.)
   o3=rhsh(r+h/2.,p+h*k2/2.,m+h*l2/2.,beta+h*p2/2.)
   p3=rhsb(r+h/2.,p+h*k2/2.,m+h*l2/2.,beta+h*p2/2.,caph+h*o2/2.)
   
   !4.RK
   k4=rhsp(r+h,p+h*k3,m+h*l3)
   l4=rhsm(r+h,p+h*k3,m+h*l3)
   m4=rhsn(r+h,p+h*k3,m+h*l3)
   n4=rhsr(r+h,p+h*k3,m+h*l3)
   o4=rhsh(r+h,p+h*k3,m+h*l3,beta+h*p3)
   p4=rhsb(r+h,p+h*k3,m+h*l3,beta+h*p3,caph+h*o3)
   
   !Abbruchbedingung -> Termination condition
!#->
   if (.not.(p+(h/6.)*(k1+2.*k2+2.*k3+k4))>6.e-17) exit!>6.e-17; e-13 for newls220 !org!!!
   !if (.not.(p+(h/6.)*(k1+2.*k2+2.*k3+k4))>6.e-17) exit
!    pexit=p+(h/6.)*(k1+2.*k2+2.*k3+k4)
!    rhoexit= 4e14 /6.176e17
!    if (rest(pexit)<rhoexit) exit

   !Gesamtschritt -> Overall step
   p=p+(h/6.)*(k1+2.*k2+2.*k3+k4)
   m=m+(h/6.)*(l1+2.*l2+2.*l3+l4)
   nue=nue+(h/6.)*(m1+2.*m2+2.*m3+m4)
   mrest=mrest+(h/6.)*(n1+2.*n2+2.*n3+n4)
   caph=caph+(h/6.)*(o1+2.*o2+2.*o3+o4)
   beta=beta+(h/6.)*(p1+2.*p2+2.*p3+p4)
   r=r+h


!if (mod(i,100)==1) then
       !fill histgram
       !drho=0.005 oben gesetzt
    !ll=floor(ye(i)/dye)
    ll=floor((log10(rest(p)*6.176e17)-13.)/drho)+1
!!    ll=nint(((rest(p)*6.176e17)-1.e13)/1.e13)

    if (ll>0.and.ll<=nhist) pmhist(ll)=pmhist(ll)+(h/6.)*(n1+2.*n2+2.*n3+n4)
!endif

!write(*,*)i,r,p
!read(*,*)dum
!#->
if (mod(i,100)==1) write(1744,*)r,p,eps(p),m,mrest,rest(p),(h/6.)*(n1+2.*n2+2.*n3+n4)

enddo

!#->
close(1744)
!#->    to search a configuration with a certain mass but then !!! comment stop at mmax
! if (m>=1.35) then
!   write(*,*)"last mass close to requested (->stop)",m,1.5
!   !stop
!   rho0=rho0-rhostep
!   rhostep=rhostep*0.5
! endif

  open(unit=1111,file='asshisto_tov.dat',status='replace',form='formatted')
  do ll=1,nhist
    write(1111,*)yehist(ll),mhist(ll),pmhist(ll),rhohist(ll),sum(pmhist(1:ll))
  enddo
  close(unit=1111)
  
  
  !note that there are small discontinuities in the histogram. These are a consequce of the EoS table and its interpolation. The jump occur at the points of the eos. Maybe it would be better to interpolate in logspace!!!


if (i/=n) then

! raus if i search a special mass
   ! stop when maximum m is reached
!    if (oldmmax>m.and.m>1.2) then
!       stop
!    else
!       oldmmax=m
!    endif


   nue0=log(sqrt(1.-2.*G*m/(r*c**2)))-nue  !so , dass lapse gg 1 geht im unendlichen
   nue=nue0+nue
  ! write(*,*)rho0,r,m,mrest,exp(nue),nue0
   
   ! compute tidal parameters and perform rescaling of a0
   
   y=r*beta/caph - 4.*3.141592654*r*r*r*eps(p)/m
   
   comp=m/r!*0.99
   
   !nenner is the denominator
   nenner = 2.*comp*(6.-3.*y+3.*comp*(5.*y-8.)) &
            + 4.*comp*comp*comp*(13.-11.*y+comp*(3.*y-2)+2.*comp*comp*(1.+y)) &
            + 3.*(1.-2.*comp)*(1.-2.*comp)*(2.-y+2.*comp*(y-1.))*log(1.-2.*comp)
   
   tidalk2 = 8.*(comp**5)/5.*(1.-2.*comp)*(1.-2.*comp) * (2.+2.*comp*(y-1.)-y) / nenner  ! love number
   
   
   kappa2 = 2.*(1./2.**5 * tidalk2/comp**5 + 1./2.**5 *tidalk2/comp**5) ! tidal coupling constant assuming q=1 = symmetric mergers
   
   lambda = 2./3.*tidalk2*(1./comp)**5  ! dimensionless tidal deformability
   
   ef=rho0*6.176e17*5.6175e-13
	
	 do kk=2,4000
    if (narray(kk) .gt. 0) then
      if ((e0array(kk) .gt. ef) .and. (e0array(kk-1) .le. ef)) then
        R2 = e0array(kk) 
        R1 = e0array(kk-1)
        LL2 = narray(kk)
        LL1 = narray(kk-1)
      endif
     ! write(*,*) e0array(k),' ',e0array(k-1),' ',narray(k),' ',narray(k-1)
     ! write(*,*) e0array(kk),' ',narray(kk)
      endif
    enddo
    
!write(*,*) R2,' ',R1,' ',L2,' ',L1
!The linear interpolation:
    lambda0 = (ef-R1)/(R2-R1)
    nf   = lambda0*LL2  + (1.0-lambda0)*LL1
    
    !write(*,*) ' ef= ', ef, ' nf= ', nf
   
   ! here the output is written
   
   
   write(1745,*)rho0,1.476*r,m,mrest,exp(nue),nue0,(m/r),restini,pini,tidalk2,kappa2,lambda,y,eps(p)*6.176e17,caph,beta,ef,nf
   write(*,*)rho0,1.476*r,m,mrest,exp(nue),nue0,(m/r),restini,pini,tidalk2,kappa2,lambda,y,eps(p)*6.176e17,caph,beta,ef,nf

else
   write(*,*)"Ende der Schleife ohne Konvergenz"
endif 


if (strange/=1) then
   factorr=1   !Umrechnung fuer ruhenmassendichte
else
   factorr=(1.658e15)  !ist evtl nicht ganz korrekt!!!!
endif

 enddo !ende der Schleife ueber verschiedene Zentraldichten
! write(1746,*)r,m,p,rbar

close(1745)


contains

integer function find_upper_index(array, x, nused)

    implicit none

    integer, intent(in) :: nused
    real, dimension(:), intent(in) :: array
    real, intent(in) :: x

    integer :: lo, hi, mid

    ! Same behaviour as the old search for values
    ! below the first table entry.
    if (x < array(1)) then
        find_upper_index = 1
        return
    endif

    ! Signal that x lies outside the upper table boundary.
    if (x >= array(nused)) then
        find_upper_index = nused + 1
        return
    endif

    lo = 1
    hi = nused

    ! Find the first index satisfying array(index) > x.
    do while (hi - lo > 1)

        mid = (lo + hi) / 2

        if (array(mid) > x) then
            hi = mid
        else
            lo = mid
        endif

    enddo

    find_upper_index = hi

end function find_upper_index

real function rhsp(rr,pp,mm)
real,intent(in)::rr,pp,mm
real::pi
real::ee

pi=3.1415927654

if (m==0) then
   rhsp=0.
else

   ee = eps(pp)

   rhsp=-(G*ee*mm/rr**2) &
        *(1.+pp/(ee*c**2)) &
        *(1.+4.*pi*rr**3*pp/(mm*c**2)) &
        /(1.-2.*G*mm/(rr*c**2))

endif
return

end function rhsp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function rhsm(rr,pp,mm)
real,intent(in)::rr,pp,mm
real::pi
pi=3.1415927654

rhsm=4.*pi*rr**2*eps(pp)
return
end function rhsm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function rhsn(rr,pp,mm)
real,intent(in)::rr,pp,mm
real::pi
pi=3.1415927654

rhsn=(G*m/(c**2*r**2))*(1.+4.*pi*r**3*pp/(m*c**2))/(1.-2.*G*m/(c**2*r))
return
end function rhsn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function rhsr(rr,pp,mm)
real,intent(in)::rr,pp,mm
real::pi
pi=3.1415927654

if (strange/=1) then
  rhsr=4.*pi*rr**2*rest(pp)/sqrt(1.-2.*G*mm/(rr*c**2)) 
else
   rhsr=4.*pi*rr**2*rest(pp)*(1.658e15)/sqrt(1.-2.*G*mm/(rr*c**2)) !incl umrechnung 1/fm^3 nach g/cm^3
endif
!Umrechnungsfaktor experimentell bestimmt!!!

return
end function rhsr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function rhsh(rr,pp,mm,bb)
real,intent(in)::rr,pp,mm,bb

  rhsh=bb

return
end function rhsh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function rhsb(rr,pp,mm,bb,hh)
real,intent(in)::rr,pp,mm,bb,hh
real::pi
real :: ee
real::dum !dummy
real::f
real::erg1,erg2
pi=3.1415927654

  dum = 1./(1.-2.*mm/rr)
  ee = eps(pp)
  f = (eps(pp+0.01*pp)-eps(pp-0.01*pp))/(0.02*pp)
  
  !f = (eps(pp+0.1*pp)-eps(pp-0.1*pp))/(0.2*pp)
    
  erg1 = 2.*dum*hh* &
         (-2.*pi*(5*ee+9.*pp+f*(ee+pp))+3./rr/rr+2.*dum &
           *(mm/rr/rr+4.*pi*rr*pp)*(mm/rr/rr+4.*pi*rr*pp))
         
  erg2 = 2.*bb/rr*dum*(-1.+mm/rr+2.*pi*rr*rr*(ee-pp))

  rhsb = erg1 + erg2

return
end function rhsb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function eps(pp)   !Energiedichte
real,intent(in)::pp
real::pi,rrho,energ,yp
!real::aa,eeps0
pi=3.1415927654

if (strange/=1) then

   if (ipoly/=1) then
         call eosshen2d(rrho,pp,energ,rhoarray,parray,harray,earray, yparray)
   else
         call EoSpoly(rrho,energ,pp,kpoly,gamma)
   endif

  eps=energ
else
  eps=pp/(a*c*c)+eps0
endif

return
end function eps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function pout(ee)   !Druck
real,intent(in)::ee
real::pi,rrho,energ,yp,ppp
!real::aa,eeps0
pi=3.1415927654

if (strange/=1) then
  call poutshen(rrho,ppp,ee,rhoarray,parray,harray,earray, yparray)
  pout=ppp
else
  pout=a*c*c*(rho0-eps0)
endif

return
end function pout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function rest(pp)   !Ruhemassendichte
real,intent(in)::pp
real::pi,rrho,energ,yp
!real::aa,eeps0
pi=3.1415927654

if (strange/=1) then

if (ipoly/=1) then
  call eosshen2d(rrho,pp,energ,rhoarray,parray,harray,earray, yparray)
else
  call EoSpoly(rrho,energ,pp,kpoly,gamma)
endif

  rest=rrho
else
  rest=rest0*(1.+((1.+a)/a)*(pp/(eps0*c**2)))**(1./(1.+a))
endif

return
end function rest

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eosinit2(rhoarray,parray,harray,earray,yparray,length,narray,e0array)

  real, dimension(:), intent(out):: rhoarray,parray,harray,earray,yparray,narray,e0array
  integer, intent(out):: length
  integer::i,idum
  real::   rdum1,rdum2,mb
  
  
  character*400::  eosfile,adum
  
open(unit=1746,file='infile',status='unknown',form='formatted')
read(1746,*) rdum1,rdum2, idum
read(1746,'(A)') eosfile
close(1746)
  
  
open(unit=10,file=trim(eosfile),form='formatted')
!open(unit=10,file='/home/bausweas/work/eos/tables/eos_akmalpr.166',form='formatted')

  i=0
  read(10,'(A)')adum
  write(*,*)'Reading Eos',adum
  do
     i=i+1

  !fuer neue eosse: read(10,*,end=99)harray(i),rhoarray(i),parray(i),earray(i), yparray(i)

!    !!!!!! my format  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!    read(10,fmt='(5e13.5)',end=99)harray(i),rhoarray(i),parray(i),earray(i), yparray(i)
! 
   !!!!!! David's format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    read(10,*,end=99)earray(i),parray(i),rhoarray(i),harray(i)
   
   narray(i)=rhoarray(i)
   e0array(i)=earray(i)
   ! convert to geom units
   earray(i)=earray(i)*1.783e12    !1.783e-27*1.e39=1.783e12
   earray(i)=earray(i)/6.176e17
   
   parray(i)=parray(i)*1.783e12    !1.783e-27*1.e39=1.783e12
   parray(i)=parray(i)/6.176e17   
   
   mb = 933.  ! fiducial baryon mass
   
   rhoarray(i)=rhoarray(i)*mb      ! mev/fm^3
   rhoarray(i)=rhoarray(i)*1.783e12    !1.783e-27*1.e39=1.783e12
   rhoarray(i)=rhoarray(i)/6.176e17 
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
!evtl offset rueckgaengig machen -> reverse offset
!#->
    !  earray(i)=earray(i)-9.8E-3*rhoarray(i)

    
    rhoarray(i) = log10(rhoarray(i))
    parray(i)   = log10(parray(i))
    earray(i)   = log10(earray(i))
    
  enddo
99 close(unit=10)
  write(*,*)i-1,rhoarray(i-1)*6.176e17,earray(i-1)*6.176e17,parray(i-1)*6.676e17*2.9979e10*2.9979e10
 write(*,*)i-1,rhoarray(i-1),earray(i-1),parray(i-1)
  length=i-1

end subroutine eosinit2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EosShen2D(rrho,pp,e,rhoarray,parray,harray,earray,yparray)

   real, intent(in):: pp
   real, dimension(:), intent(in):: rhoarray,parray,harray,earray,yparray
    real, intent(out):: rrho,e
real::dump,u  !dummy druck
    integer:: i,it,nn

!p von cgs auf geom transformieren:
dump=log10(pp)  !/(6.676e17*2.9979e10)
!dump=dump/(2.9979e10)  !irgendwiw kann ersdas nicht auf einmla rechnen
!write(*,*)"dump",pp,dump
!alles nachfolgende in geom einheiten

  nn=length
! Find the position in the table
  it = find_upper_index(parray, dump, nn)
  if (it>nn) then
     write(*,*)'In EOS-lookup: Bounds exceeded!!!',rho,rhoarray(1),rhoarray(n)
! Do the interpolation
  else if (it==1) then
     rrho=rhoarray(2)+(rhoarray(1)-rhoarray(2))*(dump-parray(2))/ &
          (parray(1)-parray(2))
     e=earray(2)+(earray(1)-earray(2))*(dump-parray(2))/ &
          (parray(1)-parray(2))

  else
     rrho=rhoarray(it-1)+(rhoarray(it)-rhoarray(it-1))*(dump-parray(it-1))/ &
          (parray(it)-parray(it-1))
     e=earray(it-1)+(earray(it)-earray(it-1))*(dump-parray(it-1))/ &
          (parray(it)-parray(it-1))
  endif

  rrho = 10.d0**rrho
  e    = 10.d0**e
  
!noch auf cgs transformieren:
!rrho=rrho*6.176e17  !von geom auf g/cm^3
!e=e*6.176e17  !von geom auf g/cm^3
!rest ist egal, weil nicht gebraucht

end subroutine EosShen2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine poutshen(rrho,pp,e,rhoarray,parray,harray,earray,yparray)

 real, intent(in):: e
   real, dimension(:), intent(in):: rhoarray,parray,harray,earray,yparray
    real, intent(out):: pp,rrho
real::dump,u  !dummy druck
    integer:: i,it,nn

!p von cgs auf geom transformieren:
dump=log10(e)  !/(6.676e17)

!alles nachfolgende in geom einheiten

  nn=length
! Find the position in the table
  it = find_upper_index(earray, dump, nn)
  write(*,*)i,nn,it
  if (it>=nn) then
     write(*,*)'In EOS-lookup: Bounds exceeded!!!',rho,dump!rhoarray(1),rhoarray(n),dump
     write(*,*)maxval(earray),dump
     stop
! Do the interpolation
  else if (it==1) then
     rrho=rhoarray(2)+(rhoarray(1)-rhoarray(2))*(dump-earray(2))/ &
          (earray(1)-earray(2))
     pp=parray(2)+(parray(1)-parray(2))*(dump-earray(2))/ &
          (earray(1)-earray(2))

  else
     rrho=rhoarray(it-1)+(rhoarray(it)-rhoarray(it-1))*(dump-earray(it-1))/ &
          (earray(it)-earray(it-1))
     pp=parray(it-1)+(parray(it)-parray(it-1))*(dump-earray(it-1))/ &
          (earray(it)-earray(it-1))
  endif

  rrho = 10.d0**rrho
  pp   = 10.d0**pp
  
!noch auf cgs transformieren:
!rrho=rrho*6.176e17  !von geom auf g/cm^3
!pp=pp*(6.676e17)
!pp=pp*2.9979e10
!pp=pp*2.9979e10
!rest ist egal, weil nicht gebraucht

end subroutine poutshen

!   subroutine poutEoSpoly(rrho,pp,e,k,gamma)
! 
!     real, intent(in):: e,k,gamma
!     real, intent(out):: pp,rrho
! 
!  
! 
!   end subroutine poutEoSpoly

  subroutine EoSpoly(rrho,e,pp,k,gama)

    real, intent(in):: pp,k,gama
    real, intent(out):: rrho,e

 rrho=(pp/k)**(1./gama)
 e=rrho+pp/(gama-1.) 


  end subroutine EoSpoly

end program






