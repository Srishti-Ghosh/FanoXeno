program main
implicit none
double precision Pnx, x, dx
integer*4 n,i,imax
n = 0
read (5,*) n
x = -1.0d0
dx = 1.0d0/1024
imax = 2049
do i = 1, imax
write (6,*) i, x, Pnx(n,x)
x = x + dx       
enddo
end program
