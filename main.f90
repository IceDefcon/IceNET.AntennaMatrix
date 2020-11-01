program main

    use definitions
    use functions

    pi=4.*atan2(1.,1.)
    k=2.*pi

    ! Assign weights and values for quadrature summation : 7 point

    eta1(1)=0.0
    eta1(2)=1./2.
    eta1(3)=1.
    eta1(4)=1./2.
    eta1(5)=0.
    eta1(6)=0.
    eta1(7)=1./3.

    eta2(1)=0.0
    eta2(2)=0.0
    eta2(3)=0.0
    eta2(4)=1./2.
    eta2(5)=1.
    eta2(6)=1./2.
    eta2(7)=1./3.

    w(1)=1./40.
    w(2)=1./15.
    w(3)=1./40.
    w(4)=1./15.
    w(5)=1./40.
    w(6)=1./15.
    w(7)=9./40.


!   Define source triangle vertices
    Sx1=0.0;Sy1=-0.0;Sz1=0.0; Sx2=0.05;Sy2=-0.0;Sz2=0.0; Sx3=0.0;Sy3=0.05;Sz3=0.0

!   Define test triangle vertices
    Tx1=0.1;Ty1=0.0;Tz1=0.0; Tx2=0.15;Ty2=0.0;Tz2=0.0; Tx3=0.1;Ty3=0.05;Tz3=0.0

!   Define test centroid
    rfpx=(Tx1+Tx2+Tx3)/3.; rfpy=(Ty1+Ty2+Ty3)/3.; rfpz=(Tz1+Tz2+Tz3)/3.

!   Triangle area

    Ar=0.5*(Sx2-Sx1)*(Sy3-Sy1)

!   perform 7 point integration

!   Real result
    ir=1
    realRes=func1();

!   Imag result

    aImagRes=func2()

    print *,''
    print *,realRes,aImagRes
    write(*,103) realRes*2.*Ar,aImagRes*2.*Ar

103 format(/,'scalar potential = ',e12.4,' +j',e12.4,/)

    stop
end

