module spinflip
use constants
use rand

implicit none

contains 

subroutine flip(p_rev,e_ij,s_ij,n,s_new)

real(kind=dp),intent(in) :: e_ij,p_rev
integer, intent(in) :: s_ij,n
integer, intent(out) :: s_new 

real(kind=dp) :: r

if (e_ij.lt.0) then
  s_new = -s_ij
else if (e_ij.gt.0) then
  r = ran()
  if (r.lt.p_rev) then
    s_new = -s_ij
  else
    s_new = s_ij
  end if
else if (e_ij.eq.0) then
  s_new = s_ij
end if    

end subroutine

end module  
