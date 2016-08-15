module energy
use constants 
implicit none

contains 

subroutine en(t,sum_s,s_sign,p_rev,e_ij,delta_e,h)
 
real(kind=dp),intent(in) :: t,sum_s,h
real(kind=dp),intent(inout) :: p_rev,e_ij,delta_e
integer,intent(in) :: s_sign


! energy at given site (J=1)
e_ij = 1.0_dp*real(s_sign)*(sum_s+h)

! change in energy
delta_e = 2.0_dp*e_ij

! probabilityofreversal
p_rev = exp(-delta_e/t)


end subroutine

end module
