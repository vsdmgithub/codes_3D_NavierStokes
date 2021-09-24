module tools
implicit none
contains
subroutine print_sec(m,dta,c_test)
  implicit none
  integer(kind=4),intent(in)::m
  double precision,dimension(0:m-1,0:m-1),intent(in)::dta
  integer(kind=4)::k,j
  character(len=*),intent(in)::c_test
  do j=0,m-1
  do k=0,m-1
    print*,dta(j,k),'|'
  end do
    print*,' '
end do
print*,trim(adjustl(c_test))
end
end module tools
program main
  USE tools
implicit none
integer(kind=4)::i,j,k,N
double precision,dimension(0:2,0:2,0:2)::val

N=3
do i=0,N-1
do j=0,N-1
do k=0,N-1
  val(i,j,k)=(i+1)*100+(j+1)*10+k+1
end do
end do
end do

call print_sec(N,val(2,:,:),'char_test')
end
