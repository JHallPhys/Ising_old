module spin_prob

use constants
use energy
use nn

implicit none 

contains 

!--------------------------------
!-----------------------------
! prob matrix for given temp
!-----------------------------
!--------------------------------


subroutine prob(A,t,n,p_s,sum_en,flip_en,h)
! add in a case select to reduce calc for all up or all down 
integer,intent(in) :: n
integer,dimension(1:n,1:n),intent(in) :: A
real(kind=dp),intent(in) :: t,h
real(kind=dp),dimension(1:n,1:n),intent(out) :: p_s,sum_en,flip_en


real(kind=dp) :: sum_s,en_ij,delta_e,p_rev
integer :: i,j,k,l,p,s_sign

! move sequentially across the grid 
do i = 1,n
  do j = 1,n
   sum_s = 0
   do k = 1,4
      p = 0
      call b_c(i,j,k,n,p)
      if (p==0) then
        print*, 'im breaking your probability subroutine'
        return
      end if
      if (k.le.2) then
        sum_s = sum_s + A(p,j)
      else 
        sum_s = sum_s + A(i,p)     
      end if
    end do 
    s_sign = A(i,j)
    call en(t,sum_s,s_sign,p_rev,en_ij,delta_e,h)
    p_s(i,j) = p_rev  
    sum_en(i,j) = en_ij
    flip_en(i,j) = delta_e
  end do
end do  

end subroutine


end module 
