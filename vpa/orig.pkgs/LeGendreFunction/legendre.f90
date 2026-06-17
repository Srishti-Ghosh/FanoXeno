function Pnx(n,x)
! Calculate Legendre polynomials
! using Bonnet's recursion formula.
implicit none
double precision Pnx,x,tpnx,tpnx0,tpnx1
integer*4 n,id
if (n == 0) then
Pnx = 1
elseif (n == 1) then
Pnx = x
else
tpnx0 = 1
tpnx1 = x
do id=2,n
tpnx = ((2.d0*dble(id)-1.d0)*x*tpnx1 - dble(id-1)*tpnx0)/dble(id)
tpnx0 = tpnx1
tpnx1 = tpnx
enddo
Pnx = tpnx
endif
end function
