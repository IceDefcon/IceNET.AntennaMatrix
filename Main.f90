program main

    use Definitions
    use func

    complex SPconst
    complex SP13,SP23,SP14,SP24
    complex VP13C1x,VP23C1x,VP14C3x,VP24C3x
    complex VP13C1y,VP23C1y,VP14C3y,VP24C3y
    complex Z13
    complex j

    j=cmplx(0.0,1.0)

    pi=4.*atan2(1.,1.)
    k=2.*pi
    umu0=4.*pi*1.e-7
    eps0=8.854*1.e-12
    omega=2.*pi*300.*1.e6

    VPconst=umu0/(4.*pi)
    SPconst=-1./(4.*pi*j*omega*eps0)

! Assign weights and values for quadrature summation : 7 point

    Xi(1)=0.0
    Xi(2)=1./2.
    Xi(3)=1.
    Xi(4)=1./2.
    Xi(5)=0.
    Xi(6)=0.
    Xi(7)=1./3.

    Eta(1)=0.0
    Eta(2)=0.0
    Eta(3)=0.0
    Eta(4)=1./2.
    Eta(5)=1.
    Eta(6)=1./2.
    Eta(7)=1./3.

    w(1)=1./40.
    w(2)=1./15.
    w(3)=1./40.
    w(4)=1./15.
    w(5)=1./40.
    w(6)=1./15.
    w(7)=9./40.

!   Define source triangle vertices for triangle S3 in half wave dipole
!   Use convention of choosing vertices in clockwise direction starting at bottom (left)

    S3x1=0.05;S3x2=0.05;S3x3=0.1
    S3y1=0.0;S3y2=0.05;S3y3=0.0
    S3z1=0.0;S3z2=0.0;S3z3=0.0

!   Define source triangle vertices for triangle S4 in half wave dipole
!   Use convention of choosing vertices in clockwise direction starting at bottom (left)

    S4x1=0.1;S4x2=0.05;S4x3=0.1
    S4y1=0.0;S4y2=0.05;S4y3=0.05
    S4z1=0.0;S4z2=0.0;S4z3=0.0

!   Define test triangle vertices for triangle T1 in same way

    T1x1=0.0;T1x2=0.0;T1x3=0.05
    T1y1=0.0;T1y2=0.05;T1y3=0.0;
    T1z1=0.0;T1z2=0.0;T1z3=0.0

!   Define test triangle vertices for triangle T2 in same way

    T2x1=0.05;T2x2=0.0;T2x3=0.05
    T2y1=0.0;T2y2=0.05;T2y3=0.05;
    T2z1=0.0;T2z2=0.0;T2z3=0.0

!   Define test centroid for triangle T1 corner 1
    rfpx1=(T1x1+T1x2+T1x3)/3.; rfpy1=(T1y1+T1y2+T1y3)/3.; rfpz1=(T1z1+T1z2+T1z3)/3.

!   Define test centroid for triangle T2 corner 3
    rfpx2=(T2x1+T2x2+T2x3)/3.; rfpy2=(T2y1+T2y2+T2y3)/3.; rfpz2=(T2z1+T2z2+T2z3)/3.

!   perform 7 point integration to T1 from S3

    Sx1=S3x1;Sx2=S3x2;Sx3=S3x3
    Sy1=S3y1;Sy2=S3y2;Sy3=S3y3

!   source/test triangle area

    An=0.5*(Sx3-Sx1)*(Sy2-Sy1)

!   common diagonal

    aln=sqrt((Sx3-Sx1)**2+(Sy2-Sy1)**2)
    alm=aln

    rfpx=rfpx1;rfpy=rfpy1

!   first evaluate scalar potential from source triangle to test triangle

    ir=1 ! real SP
    realSP13=scalarPotential();

    ir=2 ! imag SP
    aimagSP13=scalarPotential()

    SP13=SPconst*2.*An*(aln/An)*cmplx(realSP13,aimagSP13) ! eq(33)

!   now vector potential - this is more complex since we need to choose one of three corners of source triangle
!   and for each corner there are 3 Cartesian components

    ir=1 ! real vp
    realVPeta=vectorPotentialeta(); realVPneta=vectorPotentialneta()

    ir=2 ! imag vp
    aimagVPeta=vectorPotentialeta(); aimagVPneta=vectorPotentialneta()

    realVPx=sx1*realVPeta+sx2*realVPneta+sx3*(realSP13-realVPeta-realVPneta)
    realVPx13corner1=realVPx-sx1*realSP13
    realVPx13corner2=realVPx-sx2*realSP13
    realVPx13corner3=realVPx-sx3*realSP13

    realVPy=sy1*realVPeta+sy2*realVPneta+sy3*(realSP13-realVPeta-realVPneta)
    realVPy13corner1=realVPy-sy1*realSP13
    realVPy13corner2=realVPy-sy2*realSP13
    realVPy13corner3=realVPy-sy3*realSP13

    ir=2 ! imag part
    aimagVPx=sx1*aimagVPeta+sx2*aimagVPneta+sx3*(aimagSP13-aimagVPeta-aimagVPneta)
    aimagVPx13corner1=aimagVPx-sx1*aimagSP13
    aimagVPx13corner2=aimagVPx-sx2*aimagSP13
    aimagVPx13corner3=aimagVPx-sx3*aimagSP13

    aimagVPy=sy1*aimagVPeta+sy2*aimagVPneta+sy3*(aimagSP13-aimagVPeta-aimagVPneta)
    aimagVPy13corner1=aimagVPy-sy1*aimagSP13
    aimagVPy13corner2=aimagVPy-sy2*aimagSP13
    aimagVPy13corner3=aimagVPy-sy3*aimagSP13

    VP13C1x=VPconst*2.*An*(aln/(2.*An))*(cmplx(realVPx13corner1,aimagVPx13corner1)) ! eq(32)
    VP13C1y=VPconst*2.*An*(aln/(2.*An))*(cmplx(realVPy13corner1,aimagVPy13corner1))

!   perform 7 point integration to T2 from S3

    Sx1=S3x1;Sx2=S3x2;Sx3=S3x3
    Sy1=S3y1;Sy2=S3y2;Sy3=S3y3

    rfpx=rfpx2;rfpy=rfpy2

!   first evaluate scalar potential from source triangle to test triangle

    ir=1 ! real SP
    realSP23=scalarPotential();

    ir=2 ! imag SP
    aimagSP23=scalarPotential()

    SP23=SPconst*2.*An*(aln/An)*cmplx(realSP23,aimagSP23) ! eq(33)

!   now vector potential - this is more complex since we need to choose one of three corners of source triangle
!   and for each corner there are 3 Cartesian components

    ir=1 ! real vp
    realVPeta=vectorPotentialeta(); realVPneta=vectorPotentialneta()

    ir=2 ! imag vp
    aimagVPeta=vectorPotentialeta(); aimagVPneta=vectorPotentialneta()

    realVPx=sx1*realVPeta+sx2*realVPneta+sx3*(realSP23-realVPeta-realVPneta)
    realVPx23corner1=realVPx-sx1*realSP23
    realVPx23corner2=realVPx-sx2*realSP23
    realVPx23corner3=realVPx-sx3*realSP23

    realVPy=sy1*realVPeta+sy2*realVPneta+sy3*(realSP23-realVPeta-realVPneta)
    realVPy23corner1=realVPy-sy1*realSP23
    realVPy23corner2=realVPy-sy2*realSP23
    realVPy23corner3=realVPy-sy3*realSP23

    ir=2 ! imag part
    aimagVPx=sx1*aimagVPeta+sx2*aimagVPneta+sx3*(aimagSP23-aimagVPeta-aimagVPneta)
    aimagVPx23corner1=aimagVPx-sx1*aimagSP23
    aimagVPx23corner2=aimagVPx-sx2*aimagSP23
    aimagVPx23corner3=aimagVPx-sx3*aimagSP23

    aimagVPy=sy1*aimagVPeta+sy2*aimagVPneta+sy3*(aimagSP23-aimagVPeta-aimagVPneta)
    aimagVPy23corner1=aimagVPy-sy1*aimagSP23
    aimagVPy23corner2=aimagVPy-sy2*aimagSP23
    aimagVPy23corner3=aimagVPy-sy3*aimagSP23

    VP23C1x=VPconst*2.*An*(aln/(2.*An))*(cmplx(realVPx23corner1,aimagVPx23corner1)) ! eq(32)
    VP23C1y=VPconst*2.*An*(aln/(2.*An))*(cmplx(realVPy23corner1,aimagVPy23corner1))


!   perform 7 point integration to T1 from S4

    Sx1=S4x1;Sx2=S4x2;Sx3=S4x3
    Sy1=S4y1;Sy2=S4y2;Sy3=S4y3

!   source/test triangle area

    An=0.5*(Sx3-Sx2)*(Sy2-Sy1)

!   common diagonal

    aln=sqrt((Sx3-Sx2)**2+(Sy2-Sy1)**2)

    rfpx=rfpx1;rfpy=rfpy1

!   first evaluate scalar potential from source triangle to test triangle

    ir=1 ! real SP
    realSP14=scalarPotential();

    ir=2 ! imag SP
    aimagSP14=scalarPotential()

    SP14=-SPconst*2.*An*(aln/An)*cmplx(realSP14,aimagSP14) ! eq(33)

!   now vector potential - this is more complex since we need to choose one of three corners of source triangle
!   and for each corner there are 3 Cartesian components

    ir=1 ! real vp
    realVPeta=vectorPotentialeta(); realVPneta=vectorPotentialneta()

    ir=2 ! imag vp
    aimagVPeta=vectorPotentialeta(); aimagVPneta=vectorPotentialneta()

    realVPx=sx1*realVPeta+sx2*realVPneta+sx3*(realSP14-realVPeta-realVPneta)
    realVPx14corner1=realVPx-sx1*realSP14
    realVPx14corner2=realVPx-sx2*realSP14
    realVPx14corner3=realVPx-sx3*realSP14

    realVPy=sy1*realVPeta+sy2*realVPneta+sy3*(realSP14-realVPeta-realVPneta)
    realVPy14corner1=realVPy-sy1*realSP14
    realVPy14corner2=realVPy-sy2*realSP14
    realVPy14corner3=realVPy-sy3*realSP14

    ir=2 ! imag part
    aimagVPx=sx1*aimagVPeta+sx2*aimagVPneta+sx3*(aimagSP14-aimagVPeta-aimagVPneta)
    aimagVPx14corner1=aimagVPx-sx1*aimagSP14
    aimagVPx14corner2=aimagVPx-sx2*aimagSP14
    aimagVPx14corner3=aimagVPx-sx3*aimagSP14

    aimagVPy=sy1*aimagVPeta+sy2*aimagVPneta+sy3*(aimagSP14-aimagVPeta-aimagVPneta)
    aimagVPy14corner1=aimagVPy-sy1*aimagSP14
    aimagVPy14corner2=aimagVPy-sy2*aimagSP14
    aimagVPy14corner3=aimagVPy-sy3*aimagSP14

    VP14C3x=-VPconst*2.*An*(aln/(2.*An))*(cmplx(realVPx14corner3,aimagVPx14corner3)) ! eq(32)
    VP14C3y=-VPconst*2.*An*(aln/(2.*An))*(cmplx(realVPy14corner3,aimagVPy14corner3))


!   perform 7 point integration to T2 from S4

    Sx1=S4x1;Sx2=S4x2;Sx3=S4x3
    Sy1=S4y1;Sy2=S4y2;Sy3=S4y3

    rfpx=rfpx2;rfpy=rfpy2

!   first evaluate scalar potential from source triangle to test triangle

    ir=1 ! real SP
    realSP24=scalarPotential();

    ir=2 ! imag SP
    aimagSP24=scalarPotential()

    SP24=-SPconst*2.*An*(aln/An)*cmplx(realSP24,aimagSP24) ! eq(33)


!   now vector potential - this is more complex since we need to choose one of three corners of source triangle
!   and for each corner there are 3 Cartesian components

    ir=1 ! real vp
    realVPeta=vectorPotentialeta(); realVPneta=vectorPotentialneta()

    ir=2 ! imag vp
    aimagVPeta=vectorPotentialeta(); aimagVPneta=vectorPotentialneta()

    realVPx=sx1*realVPeta+sx2*realVPneta+sx3*(realSP24-realVPeta-realVPneta)
    realVPx24corner1=realVPx-sx1*realSP24
    realVPx24corner2=realVPx-sx2*realSP24
    realVPx24corner3=realVPx-sx3*realSP24

    realVPy=sy1*realVPeta+sy2*realVPneta+sy3*(realSP24-realVPeta-realVPneta)
    realVPy24corner1=realVPy-sy1*realSP24
    realVPy24corner2=realVPy-sy2*realSP24
    realVPy24corner3=realVPy-sy3*realSP24

    ir=2 ! imag part
    aimagVPx=sx1*aimagVPeta+sx2*aimagVPneta+sx3*(aimagSP24-aimagVPeta-aimagVPneta)
    aimagVPx24corner1=aimagVPx-sx1*aimagSP24
    aimagVPx24corner2=aimagVPx-sx2*aimagSP24
    aimagVPx24corner3=aimagVPx-sx3*aimagSP24

    aimagVPy=sy1*aimagVPeta+sy2*aimagVPneta+sy3*(aimagSP24-aimagVPeta-aimagVPneta)
    aimagVPy24corner1=aimagVPy-sy1*aimagSP24
    aimagVPy24corner2=aimagVPy-sy2*aimagSP24
    aimagVPy24corner3=aimagVPy-sy3*aimagSP24

    VP24C3x=-VPconst*2.*An*(aln/(2.*An))*(cmplx(realVPx24corner3,aimagVPx24corner3)) ! eq(32)
    VP24C3y=-VPconst*2.*An*(aln/(2.*An))*(cmplx(realVPy24corner3,aimagVPy24corner3))


!   Assemble SP and VP components into impedance elements

    Z13=alm*(j*omega*0.5*(VP13C1x*(rfpx1-T1x1)+VP13C1y*(rfpy1-T1y1)+VP14C3x*(rfpx1-T1x1)+VP14C3y*(rfpy1-T1y1)+& ! eq(17)
                          VP23C1x*(T2x3-rfpx2)+VP23C1y*(T2y3-rfpy2)+VP24C3x*(T2x3-rfpx2)+VP24C3y*(T2y3-rfpy2))+&
                          SP23+SP24-SP13-SP14)
    write(*,107)real(Z13),aimag(Z13),sqrt((real(Z13))**2+(aimag(Z13))**2)

107 format(/,'Z13 = ',e12.4,'+j',e12.4,' ---> magnitude =',e12.4)

    stop
end

