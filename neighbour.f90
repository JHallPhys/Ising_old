module nn
use constants 
implicit none 

contains 


subroutine b_c(i,j,k,n,p)

integer,intent(in) :: i,j,k,n
integer,intent(out) :: p  

if (k==1) then
  p = i + 1  
else if (k==2) then
  p = i - 1  
else if (k==3) then 
  p = j + 1  
else
  p = j - 1  
end if 

if (p==0) then 
  p = n
else if (p==n) then
  p = n
else
  p = modulo(p,n)
end if 

end subroutine 


end module 
