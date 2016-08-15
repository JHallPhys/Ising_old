program test 

use constants
use rand 
use spinflip
use energy
use nn
use spin_prob


implicit none 

integer,dimension(:,:),allocatable :: lattice
real(kind=dp),dimension(:,:),allocatable :: p_s 
real(kind=dp),dimension(:,:),allocatable :: sum_en
real(kind=dp),dimension(:,:),allocatable :: flip_en
real(kind=dp),dimension(:,:),allocatable :: curie
 

real(kind=dp) :: t,den_s,p_ij,sum_s,m_t0,m_av,mct_0,en_ij
real(kind=dp) :: en_av,en_t0,ensqr_av,msqr_av,ensqr_t0,msqr_t0
real(kind=dp) :: chi,chi_prime,c_v,h,p_rev,m_sum,en_sum,norm,abs_m,abs_mav
integer :: n,i,j,step_max,s_new,k,p,s_sign,s_ij
integer :: time,tick,tock,status,t_0,t_1,counter

open(file='mt.dat',unit=24,status='replace')
open(file='en.dat',unit=25,status='replace')
open(file='mtsqr.dat',unit=26,status='replace')
open(file='ensqr.dat',unit=27,status='replace')
open(file='msus.dat',unit=28,status='replace')
open(file='cv.dat',unit=29,status='replace')
open(file='cvp.dat',unit=30,status='replace')
open(file='mabs.dat',unit=31,status='replace')
open(file='slow.dat',unit=32,status='replace')
! arrays params(nxn) and mc steps 
n =  4
t_0 = 0
t_1 = 1!*t_0

! physical variables 
read*, t
!t = 0.0_dp
h = 0.0_dp

! Array allocation 
allocate(lattice(1:n,1:n),stat=status)
if (status.ne.0) stop 'issue allocating lattice'
allocate(curie(1:500,1:2),stat=status)
if (status.ne.0) stop 'issue allocating curie'

lattice(:,:) = 1
step_max = n**2
counter = 0
norm = 1.0_dp/((t_1-t_0)*real(step_max))
!----------------------------------
!       Main program
!----------------------------------
!do tock = 1,50
! number of times temperature will be incremented
do tick = 1,2000
!----------------------------------
! Ensure Average Quantites Reset 
!----------------------------------
allocate(sum_en(1:n,1:n),stat=status)
if (status.ne.0) stop 'issue allocating energy array'
sum_en(:,:) = 0
m_av = 0
msqr_av = 0
abs_mav = 0
en_av = 0
ensqr_av = 0
m_t0 = 0
msqr_t0 = 0
abs_m = 0
en_t0 = 0
ensqr_t0 = 0
!----------------------------------
!    Monte Carlo loop
!----------------------------------
!print*, (t/5.0_dp)*100.0_dp,'% done'
do time = 1,t_1
! Record observables upon reaching equlibrium 
!------------------------------------------------
if (time.eq.t_0) then
  m_t0 = real(sum(lattice))
  msqr_t0 = m_t0**2
  abs_m = sqrt(msqr_t0)

  en_t0 = 0.5_dp*sum(sum_en)
  ensqr_t0 = (en_t0)**2
end if
!call prob(lattice,t,n,p_s,sum_en,flip_en,h)

do i = 1,n
  do j = 1,n
    sum_s = 0
    do k = 1,4
      p = 0
      call b_c(i,j,k,n,p)
      if (p==0) then
        print*, 'periodic boundary not working'
      end if
      if (k.le.2) then
        sum_s = sum_s + lattice(p,j)
      else 
        sum_s = sum_s + lattice(i,p)     
      end if
    end do 
    s_sign = lattice(i,j)
    call en(t,sum_s,s_sign,p_rev,en_ij,den_s,h)
    !p_s = p_rev  
    sum_en(i,j) = en_ij
    if (den_s.gt.0) then
      if (ran().lt.p_rev) then
        lattice(i,j) = -lattice(i,j)
      end if
    else 
      lattice(i,j)= -lattice(i,j)
    end if  
  end do
end do


!----------------------------
!  Obtain Observables 
!----------------------------
if (time.gt.t_0) then
  m_t0 = real(sum(lattice))
  msqr_t0 = m_t0**2
  abs_m = sqrt(msqr_t0)

  en_t0 = 0.5_dp*sum(sum_en)
  ensqr_t0 = (en_t0)**2
end if
! update running average
m_av = m_av + m_t0
msqr_av = msqr_av + msqr_t0
abs_mav = abs_mav + abs_m

en_av = en_av + en_t0
ensqr_av = ensqr_av + ensqr_t0

end do

!-------------------------
!     mc time av
!-------------------------

! Calculate averages 

m_av = (m_av)*norm
msqr_av = (msqr_av)*norm
abs_mav = (abs_mav)*norm

en_av = (en_av)*norm
ensqr_av = (ensqr_av)*norm

!--------------------------------------
! specific heat and mag susceptibility
!--------------------------------------
if (t.gt.0.1_dp) then
  counter = counter + 1
  c_v = (1.0_dp/t**2) * (ensqr_av - (real(step_max)*(en_av)**2)) 
  !curie(counter,1) = c_v
  !curie(counter,2) = t
  chi = (1.0_dp/t)*(msqr_av - (real(step_max)*(m_av)**2))
  chi_prime = (1.0_dp/t)*(msqr_av - (real(step_max)*(abs_mav)**2))
end if

!-------------------------------
! write results to files
!-------------------------------

!write(24,*) tick,sum(lattice)
write(24,*) t,m_av
write(25,*) t,en_av 
write(26,*) t,msqr_av
write(27,*) t,ensqr_av
write(28,*) t,chi
write(29,*) t,c_v
write(30,*) t,chi_prime
write(31,*) t,abs_m
write(32,*) tick,1.0_dp/m_av


!-------------------------------------
!incremental temperature increase and 
!-------------------------------------

!t = t + 0.1_dp
deallocate(sum_en)
end do

!----------------------------------
! find curie temperature
!----------------------------------
do i = 1,size(curie)
  if (curie(i,1).eq.maxval(curie(:,1))) then
    print*, 'curie temp : ',curie(i,2)
    exit
  end if
end do

! close everything
deallocate(lattice)
deallocate(curie)


close(24)
close(25)
close(26)
close(27)
close(28)
close(29)
close(30)

end program 
