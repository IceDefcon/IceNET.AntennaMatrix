module Definitions

	integer*2 i
	integer tc,ixyz

	real,dimension(3) :: xi_3p,eta_3p,w_3p
	real,dimension(7) :: xi,eta,w
	real,dimension(20) :: rfpx_,rfpy_,rfpz_
	real Sx1,Sx2,Sx3,Sy1,Sy2,Sy3,Sz1,Sz2,Sz3
	real k,cx,cy,ct,ctx,cty,aln,alm,pi,mu0,eps0,omega
	real zeta,zeta_3p,rDashX,rDashY,rDashZ,R
	real rfpx,rfpy,rfpz,comx
	real VPconst
	real execute

	complex SPconst
	complex j,Z,G

	complex current,current_xi,current_eta,current_x,current_y
	complex SP,VPx,VPy
	complex SP11,SP12,SP21,SP22
	complex VP11x,VP12x,VP21x,VP22x
	complex VP11y,VP12y,VP21y,VP22y

end module


