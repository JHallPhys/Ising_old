module rand
use constants
implicit none

contains 

function ran()

implicit none 

real(kind=dp) :: ran
integer :: i

call system_clock(i)
call random_seed(i)
call random_number(ran)

end function 

end module 
