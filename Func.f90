module func

    use definitions

contains

!---------------------------------------------------------------------------------------------------
    complex function ipq()

        g=cmplx(0.0,0.0)

        do i=1,7

            zeta=1.-xi(i)-eta(i)

            rdashx=sx1*xi(i)+sx2*eta(i)+sx3*zeta
            rdashy=sy1*xi(i)+sy2*eta(i)+sy3*zeta
            rdashz=sz1*xi(i)+sz2*eta(i)+sz3*zeta

            r=sqrt((rfpx-rdashx)*(rfpx-rdashx)+&
                (rfpy-rdashy)*(rfpy-rdashy)+&
                (rfpz-rdashz)*(rfpz-rdashz))

            g=g+w(i)*cexp(cmplx(0.0,-k*r))/r

        enddo

    ipq=g

    end function
!---------------------------------------------------------------------------------------------------
    complex function ipq_xi()

        g=cmplx(0.0,0.0)

        do i=1,7

            zeta=1.-xi(i)-eta(i)

            rdashx=sx1*xi(i)+sx2*eta(i)+sx3*zeta
            rdashy=sy1*xi(i)+sy2*eta(i)+sy3*zeta
            rdashz=sz1*xi(i)+sz2*eta(i)+sz3*zeta

            r=sqrt((rfpx-rdashx)*(rfpx-rdashx)+&
                (rfpy-rdashy)*(rfpy-rdashy)+&
                (rfpz-rdashz)*(rfpz-rdashz))

            g=g+xi(i)*w(i)*cexp(cmplx(0.0,-k*r))/r

        enddo

    ipq_xi=g

    end function
!---------------------------------------------------------------------------------------------------
    complex function ipq_eta()

        g=cmplx(0.0,0.0)

        do i=1,7

            zeta=1.-xi(i)-eta(i)

            rdashx=sx1*xi(i)+sx2*eta(i)+sx3*zeta
            rdashy=sy1*xi(i)+sy2*eta(i)+sy3*zeta
            rdashz=sz1*xi(i)+sz2*eta(i)+sz3*zeta

            r=sqrt((rfpx-rdashx)*(rfpx-rdashx)+&
                (rfpy-rdashy)*(rfpy-rdashy)+&
                (rfpz-rdashz)*(rfpz-rdashz))

            g=g+eta(i)*w(i)*cexp(cmplx(0.0,-k*r))/r

        enddo

        ipq_eta=g

    end function
!---------------------------------------------------------------------------------------------------
    complex function ipq_3p()

        g=cmplx(0.0,0.0)

        do i=1,3

            zeta_3p=1.-xi_3p(i)-eta_3p(i)

            rdashx=sx1*xi_3p(i)+sx2*eta_3p(i)+sx3*zeta_3p
            rdashy=sy1*xi_3p(i)+sy2*eta_3p(i)+sy3*zeta_3p
            rdashz=sz1*xi_3p(i)+sz2*eta_3p(i)+sz3*zeta_3p

            r=sqrt((rfpx-rdashx)*(rfpx-rdashx)+&
                (rfpy-rdashy)*(rfpy-rdashy)+&
                (rfpz-rdashz)*(rfpz-rdashz))

            g=g+w_3p(i)*cexp(cmplx(0.0,-k*r))/r

        enddo

    ipq_3p=g

    end function
!---------------------------------------------------------------------------------------------------
    complex function ipq_xi_3p()

        g=cmplx(0.0,0.0)

        do i=1,3

            zeta_3p=1.-xi_3p(i)-eta_3p(i)

            rdashx=sx1*xi_3p(i)+sx2*eta_3p(i)+sx3*zeta_3p
            rdashy=sy1*xi_3p(i)+sy2*eta_3p(i)+sy3*zeta_3p
            rdashz=sz1*xi_3p(i)+sz2*eta_3p(i)+sz3*zeta_3p

            r=sqrt((rfpx-rdashx)*(rfpx-rdashx)+&
                (rfpy-rdashy)*(rfpy-rdashy)+&
                (rfpz-rdashz)*(rfpz-rdashz))

            g=g+xi_3p(i)*w_3p(i)*cexp(cmplx(0.0,-k*r))/r

        enddo

    ipq_xi_3p=g

    end function
!---------------------------------------------------------------------------------------------------
    complex function ipq_eta_3p()

        g=cmplx(0.0,0.0)

        do i=1,3

            zeta_3p=1.-xi_3p(i)-eta_3p(i)

            rdashx=sx1*xi_3p(i)+sx2*eta_3p(i)+sx3*zeta_3p
            rdashy=sy1*xi_3p(i)+sy2*eta_3p(i)+sy3*zeta_3p
            rdashz=sz1*xi_3p(i)+sz2*eta_3p(i)+sz3*zeta_3p

            r=sqrt((rfpx-rdashx)*(rfpx-rdashx)+&
                (rfpy-rdashy)*(rfpy-rdashy)+&
                (rfpz-rdashz)*(rfpz-rdashz))

            g=g+eta_3p(i)*w_3p(i)*cexp(cmplx(0.0,-k*r))/r

        enddo

    ipq_eta_3p=g

    end function
!---------------------------------------------------------------------------------------------------
    complex function scalar()
 
    comx=(sx1+sx2+sx3)/3.

    if (comx==rfpx) then

    current=ipq_3p();

    else

    current=ipq();

    end if

        sp=spconst*aln*current
        scalar=sp
 
    end function
!---------------------------------------------------------------------------------------------------
    complex function areavector_x()

    comx=(sx1+sx2+sx3)/3.
    if (comx==rfpx) then

    current=ipq_3p()
    current_xi=ipq_xi_3p()
    current_eta=ipq_eta_3p();

        current_x=sx1*current_xi+sx2*current_eta+sx3*(current-current_xi-current_eta)-cx*current

    else

    current=ipq()
    current_xi=ipq_xi()
    current_eta=ipq_eta()

        current_x=sx1*current_xi+sx2*current_eta+sx3*(current-current_xi-current_eta)-cx*current

    end if

        vpx=vpconst*aln*current_x
        areavector_x=vpx

    end function
!---------------------------------------------------------------------------------------------------
    complex function areavector_y()

    comx=(sx1+sx2+sx3)/3.
    if (comx==rfpx) then

    current=ipq_3p()
    current_xi=ipq_xi_3p()
    current_eta=ipq_eta_3p()

        current_y=sy1*current_xi+sy2*current_eta+sy3*(current-current_xi-current_eta)-cy*current

    else

    current=ipq()
    current_xi=ipq_xi()
    current_eta=ipq_eta()

        current_y=sy1*current_xi+sy2*current_eta+sy3*(current-current_xi-current_eta)-cy*current

    end if

        vpy=vpconst*aln*current_y
        areavector_y=vpy

    end function
!---------------------------------------------------------------------------------------------------

    real function source()


            j=cmplx(0.0,1.0)
            pi=4.*atan2(1.,1.)
            k=2.*pi
            mu0=4.*pi*1.e-7
            eps0=8.854*1.e-12
            omega=2.*pi*300.*1.e6
            vpconst=mu0/(4.*pi)
            spconst=-1./(2.*pi*j*omega*eps0)

            s1x1=0.0;  s1x2=0.0;  s1x3=0.05;
            s1y1=0.0;  s1y2=0.05; s1y3=0.0;
            s1z1=0.0;  s1z2=0.0;  s1z3=0.0

            s2x1=0.05; s2x2=0.0;  s2x3=0.05;
            s2y1=0.0;  s2y2=0.05; s2y3=0.05;
            s2z1=0.0;  s2z2=0.0;  s2z3=0.0;

            s3x1=0.05; s3x2=0.05; s3x3=0.1;
            s3y1=0.0;  s3y2=0.05; s3y3=0.0;
            s3z1=0.0;  s3z2=0.0;  s3z3=0.0;

            s4x1=0.1;  s4x2=0.05; s4x3=0.1;
            s4y1=0.0;  s4y2=0.05; s4y3=0.05;
            s4z1=0.0;  s4z2=0.0;  s4z3=0.0;

            s5x1=0.1;  s5x2=0.1;  s5x3=0.15;
            s5y1=0.0;  s5y2=0.05; s5y3=0.0;
            s5z1=0.0;  s5z2=0.0;  s5z3=0.0;

            s6x1=0.15; s6x2=0.1;  s6x3=0.15;
            s6y1=0.0;  s6y2=0.05; s6y3=0.05;
            s6z1=0.0;  s6z2=0.0;  s6z3=0.0;

            s7x1=0.15; s7x2=0.15; s7x3=0.2;
            s7y1=0.0;  s7y2=0.05; s7y3=0.0;
            s7z1=0.0;  s7z2=0.0;  s7z3=0.0;

            s8x1=0.2;  s8x2=0.15; s8x3=0.2;
            s8y1=0.0;  s8y2=0.05; s8y3=0.05;
            s8z1=0.0;  s8z2=0.0;  s8z3=0.0;

            s9x1=0.2;  s9x2=0.2;  s9x3=0.25;
            s9y1=0.0;  s9y2=0.05; s9y3=0.0;
            s9z1=0.0;  s9z2=0.0;  s9z3=0.0;

            s10x1=0.25; s10x2=0.2;  s10x3=0.25;
            s10y1=0.0;  s10y2=0.05; s10y3=0.05;
            s10z1=0.0;  s10z2=0.0;  s10z3=0.0;

            s11x1=0.25; s11x2=0.25; s11x3=0.3;
            s11y1=0.0;  s11y2=0.05; s11y3=0.0;
            s11z1=0.0;  s11z2=0.0;  s11z3=0.0;

            s12x1=0.3;  s12x2=0.25; s12x3=0.3;
            s12y1=0.0;  s12y2=0.05; s12y3=0.05;
            s12z1=0.0;  s12z2=0.0;  s12z3=0.0;

            s13x1=0.3;  s13x2=0.3;  s13x3=0.35;
            s13y1=0.0;  s13y2=0.05; s13y3=0.0;
            s13z1=0.0;  s13z2=0.0;  s13z3=0.0;

            s14x1=0.35; s14x2=0.3;  s14x3=0.35;
            s14y1=0.0;  s14y2=0.05; s14y3=0.05;
            s14z1=0.0;  s14z2=0.0;  s14z3=0.0;

            s15x1=0.35; s15x2=0.35; s15x3=0.4;
            s15y1=0.0;  s15y2=0.05; s15y3=0.0;
            s15z1=0.0;  s15z2=0.0;  s15z3=0.0;

            s16x1=0.4 ; s16x2=0.35; s16x3=0.4;
            s16y1=0.0;  s16y2=0.05; s16y3=0.05;
            s16z1=0.0;  s16z2=0.0;  s16z3=0.0;

            s17x1=0.4;  s17x2=0.4;  s17x3=0.45;
            s17y1=0.0;  s17y2=0.05; s17y3=0.0;
            s17z1=0.0;  s17z2=0.0;  s17z3=0.0;

            s18x1=0.45; s18x2=0.4;  s18x3=0.45;
            s18y1=0.0;  s18y2=0.05; s18y3=0.05;
            s18z1=0.0;  s18z2=0.0;  s18z3=0.0;

            s19x1=0.45; s19x2=0.45; s19x3=0.5;
            s19y1=0.0;  s19y2=0.05; s19y3=0.0;
            s19z1=0.0;  s19z2=0.0;  s19z3=0.0;

            s20x1=0.5;  s20x2=0.45; s20x3=0.5;
            s20y1=0.0;  s20y2=0.05; s20y3=0.05;
            s20z1=0.0;  s20z2=0.0;  s20z3=0.0;

        !*************************************************************************************************

            rfpx_(1)=(s1x1+s1x2+s1x3)/3.;     rfpy_(1)=(s1y1+s1y2+s1y3)/3.;     rfpz_(1)=(s1z1+s1z2+s1z3)/3.
            rfpx_(2)=(s2x1+s2x2+s2x3)/3.;     rfpy_(2)=(s2y1+s2y2+s2y3)/3.;     rfpz_(2)=(s2z1+s2z2+s2z3)/3.
            rfpx_(3)=(s3x1+s3x2+s3x3)/3.;     rfpy_(3)=(s3y1+s3y2+s3y3)/3.;     rfpz_(3)=(s3z1+s3z2+s3z3)/3.
            rfpx_(4)=(s4x1+s4x2+s4x3)/3.;     rfpy_(4)=(s4y1+s4y2+s4y3)/3.;     rfpz_(4)=(s4z1+s4z2+s4z3)/3.
            rfpx_(5)=(s5x1+s5x2+s5x3)/3.;     rfpy_(5)=(s5y1+s5y2+s5y3)/3.;     rfpz_(5)=(s5z1+s5z2+s5z3)/3.
            rfpx_(6)=(s6x1+s6x2+s6x3)/3.;     rfpy_(6)=(s6y1+s6y2+s6y3)/3.;     rfpz_(6)=(s6z1+s6z2+s6z3)/3.
            rfpx_(7)=(s7x1+s7x2+s7x3)/3.;     rfpy_(7)=(s7y1+s7y2+s7y3)/3.;     rfpz_(7)=(s7z1+s7z2+s7z3)/3.
            rfpx_(8)=(s8x1+s8x2+s8x3)/3.;     rfpy_(8)=(s8y1+s8y2+s8y3)/3.;     rfpz_(8)=(s8z1+s8z2+s8z3)/3.
            rfpx_(9)=(s9x1+s9x2+s9x3)/3.;     rfpy_(9)=(s9y1+s9y2+s9y3)/3.;     rfpz_(9)=(s9z1+s9z2+s9z3)/3.
            rfpx_(10)=(s10x1+s10x2+s10x3)/3.; rfpy_(10)=(s10y1+s10y2+s10y3)/3.; rfpz_(10)=(s10z1+s10z2+s10z3)/3.
            rfpx_(11)=(s11x1+s11x2+s11x3)/3.; rfpy_(11)=(s11y1+s11y2+s11y3)/3.; rfpz_(11)=(s11z1+s11z2+s11z3)/3.
            rfpx_(12)=(s12x1+s12x2+s12x3)/3.; rfpy_(12)=(s12y1+s12y2+s12y3)/3.; rfpz_(12)=(s12z1+s12z2+s12z3)/3.
            rfpx_(13)=(s13x1+s13x2+s13x3)/3.; rfpy_(13)=(s13y1+s13y2+s13y3)/3.; rfpz_(13)=(s13z1+s13z2+s13z3)/3.
            rfpx_(14)=(s14x1+s14x2+s14x3)/3.; rfpy_(14)=(s14y1+s14y2+s14y3)/3.; rfpz_(14)=(s14z1+s14z2+s14z3)/3.
            rfpx_(15)=(s15x1+s15x2+s15x3)/3.; rfpy_(15)=(s15y1+s15y2+s15y3)/3.; rfpz_(15)=(s15z1+s15z2+s15z3)/3.
            rfpx_(16)=(s16x1+s16x2+s16x3)/3.; rfpy_(16)=(s16y1+s16y2+s16y3)/3.; rfpz_(16)=(s16z1+s16z2+s16z3)/3.
            rfpx_(17)=(s17x1+s17x2+s17x3)/3.; rfpy_(17)=(s17y1+s17y2+s17y3)/3.; rfpz_(17)=(s17z1+s17z2+s17z3)/3.
            rfpx_(18)=(s18x1+s18x2+s18x3)/3.; rfpy_(18)=(s18y1+s18y2+s18y3)/3.; rfpz_(18)=(s18z1+s18z2+s18z3)/3.
            rfpx_(19)=(s19x1+s19x2+s19x3)/3.; rfpy_(19)=(s19y1+s19y2+s19y3)/3.; rfpz_(19)=(s19z1+s19z2+s19z3)/3.
            rfpx_(20)=(s20x1+s20x2+s20x3)/3.; rfpy_(20)=(s20y1+s20y2+s20y3)/3.; rfpz_(20)=(s20z1+s20z2+s20z3)/3.


        !**********************************************z_x1**********************************************

            aln=sqrt((0.05)**2+(0.05)**2)

        !   t1 from s1      t2 from s1

            sx1=s1x1;sx2=s1x2;sx3=s1x3
            sy1=s1y1;sy2=s1y2;sy3=s1y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx1; cy=sy1;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1);

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s2      t2 from s2

            sx1=s2x1;sx2=s2x2;sx3=s2x3
            sy1=s2y1;sy2=s2y2;sy3=s2y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1);

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,1)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            1 format('zx01 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x2**********************************************

            aln=0.05

        !   t1 from s2      t2 from s2

            sx1=s2x1;sx2=s2x2;sx3=s2x3
            sy1=s2y1;sy2=s2y2;sy3=s2y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx2; cy=sy2;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1);

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s3      t2 from s3

            sx1=s3x1;sx2=s3x2;sx3=s3x3
            sy1=s3y1;sy2=s3y2;sy3=s3y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1);

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,2)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            2 format('zx02 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x3**********************************************

            aln=sqrt((0.05)**2+(0.05)**2)

        !   t1 from s3      t2 from s3

            sx1=s3x1;sx2=s3x2;sx3=s3x3
            sy1=s3y1;sy2=s3y2;sy3=s3y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx1; cy=sy1;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1);

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s4      t2 from s4

            sx1=s4x1;sx2=s4x2;sx3=s4x3
            sy1=s4y1;sy2=s4y2;sy3=s4y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1);

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,3)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            3 format('zx03 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x4**********************************************

            aln=0.05

        !   t1 from s4      t2 from s4

            sx1=s4x1;sx2=s4x2;sx3=s4x3
            sy1=s4y1;sy2=s4y2;sy3=s4y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx2; cy=sy2;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx2; cy=sy2;

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s5      t2 from s5

            sx1=s5x1;sx2=s5x2;sx3=s5x3
            sy1=s5y1;sy2=s5y2;sy3=s5y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,4)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            4 format('zx04 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x5**********************************************

            aln=sqrt((0.05)**2+(0.05)**2)

        !   t1 from s5      t2 from s5

            sx1=s5x1;sx2=s5x2;sx3=s5x3
            sy1=s5y1;sy2=s5y2;sy3=s5y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx1; cy=sy1;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx1; cy=sy1

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s6      t2 from s6

            sx1=s6x1;sx2=s6x2;sx3=s6x3
            sy1=s6y1;sy2=s6y2;sy3=s6y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,5)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            5 format('zx05 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x6**********************************************

            aln=0.05

        !   t1 from s6      t2 from s6

            sx1=s6x1;sx2=s6x2;sx3=s6x3
            sy1=s6y1;sy2=s6y2;sy3=s6y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx2; cy=sy2;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx2; cy=sy2

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s7      t2 from s7

            sx1=s7x1;sx2=s7x2;sx3=s7x3
            sy1=s7y1;sy2=s7y2;sy3=s7y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,6)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            6 format('zx06 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x7**********************************************

            aln=sqrt((0.05)**2+(0.05)**2)

        !   t1 from s7      t2 from s7

            sx1=s7x1;sx2=s7x2;sx3=s7x3
            sy1=s7y1;sy2=s7y2;sy3=s7y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx1; cy=sy1;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx1; cy=sy1;

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s8      t2 from s8

            sx1=s8x1;sx2=s8x2;sx3=s8x3
            sy1=s8y1;sy2=s8y2;sy3=s8y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,7)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            7 format('zx07 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x8**********************************************

            aln=0.05

        !   t1 from s8      t2 from s8

            sx1=s8x1;sx2=s8x2;sx3=s8x3
            sy1=s8y1;sy2=s8y2;sy3=s8y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx2; cy=sy2;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx2; cy=sy2

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s9      t2 from s9

            sx1=s9x1;sx2=s9x2;sx3=s9x3
            sy1=s9y1;sy2=s9y2;sy3=s9y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,8)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            8 format('zx08 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x9**********************************************

            aln=sqrt((0.05)**2+(0.05)**2)

        !   t1 from s9      t2 from s9

            sx1=s9x1;sx2=s9x2;sx3=s9x3
            sy1=s9y1;sy2=s9y2;sy3=s9y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx1; cy=sy1;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx1; cy=sy1

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s10     t2 from s10

            sx1=s10x1;sx2=s10x2;sx3=s10x3
            sy1=s10y1;sy2=s10y2;sy3=s10y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,9)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            9 format('zx09 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x10**********************************************

            aln=0.05

        !   t1 from s10     t2 from s10

            sx1=s10x1;sx2=s10x2;sx3=s10x3
            sy1=s10y1;sy2=s10y2;sy3=s10y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx2; cy=sy2;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx2; cy=sy2

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s11     t2 from s11

            sx1=s11x1;sx2=s11x2;sx3=s11x3
            sy1=s11y1;sy2=s11y2;sy3=s11y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,10)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            10 format('zx10 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x11**********************************************

            aln=sqrt((0.05)**2+(0.05)**2)

        !   t1 from s11     t2 from s11

            sx1=s11x1;sx2=s11x2;sx3=s11x3
            sy1=s11y1;sy2=s11y2;sy3=s11y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx1; cy=sy1;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx1; cy=sy1

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s12     t2 from s12

            sx1=s12x1;sx2=s12x2;sx3=s12x3
            sy1=s12y1;sy2=s12y2;sy3=s12y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,11)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            11 format('zx11 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x12**********************************************

            aln=0.05

        !   t1 from s12     t2 from s12

            sx1=s12x1;sx2=s12x2;sx3=s12x3
            sy1=s12y1;sy2=s12y2;sy3=s12y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx2; cy=sy2;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx2; cy=sy2

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s13     t2 from s13

            sx1=s13x1;sx2=s13x2;sx3=s13x3
            sy1=s13y1;sy2=s13y2;sy3=s13y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,12)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            12 format('zx12 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x13**********************************************

            aln=sqrt((0.05)**2+(0.05)**2)

        !   t1 from s13     t2 from s13

            sx1=s13x1;sx2=s13x2;sx3=s13x3
            sy1=s13y1;sy2=s13y2;sy3=s13y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx1; cy=sy1;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx1; cy=sy1

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s14     t2 from s14

            sx1=s14x1;sx2=s14x2;sx3=s14x3
            sy1=s14y1;sy2=s14y2;sy3=s14y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,13)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            13 format('zx13 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x14**********************************************

            aln=0.05

        !   t1 from s14     t2 from s14

            sx1=s14x1;sx2=s14x2;sx3=s14x3
            sy1=s14y1;sy2=s14y2;sy3=s14y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx2; cy=sy2;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx2; cy=sy2

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s15     t2 from s15

            sx1=s15x1;sx2=s15x2;sx3=s15x3
            sy1=s15y1;sy2=s15y2;sy3=s15y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,14)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            14 format('zx14 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x15**********************************************

            aln=sqrt((0.05)**2+(0.05)**2)

        !   t1 from s15     t2 from s15

            sx1=s15x1;sx2=s15x2;sx3=s15x3
            sy1=s15y1;sy2=s15y2;sy3=s15y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx1; cy=sy1;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx1; cy=sy1

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s16     t2 from s16

            sx1=s16x1;sx2=s16x2;sx3=s16x3
            sy1=s16y1;sy2=s16y2;sy3=s16y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,15)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            15 format('zx15 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x16**********************************************

            aln=0.05

        !   t1 from s16     t2 from s16

            sx1=s16x1;sx2=s16x2;sx3=s16x3
            sy1=s16y1;sy2=s16y2;sy3=s16y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx2; cy=sy2;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx2; cy=sy2

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s17     t2 from s17

            sx1=s17x1;sx2=s17x2;sx3=s17x3
            sy1=s17y1;sy2=s17y2;sy3=s17y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,16)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            16 format('zx16 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x17**********************************************

            aln=sqrt((0.05)**2+(0.05)**2)

        !   t1 from s17     t2 from s17

            sx1=s17x1;sx2=s17x2;sx3=s17x3
            sy1=s17y1;sy2=s17y2;sy3=s17y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx1; cy=sy1;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx1; cy=sy1

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s18     t2 from s18

            sx1=s18x1;sx2=s18x2;sx3=s18x3
            sy1=s18y1;sy2=s18y2;sy3=s18y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,17)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            17 format('zx17 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x18**********************************************

            aln=0.05

        !   t1 from s18     t2 from s18

            sx1=s18x1;sx2=s18x2;sx3=s18x3
            sy1=s18y1;sy2=s18y2;sy3=s18y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx2; cy=sy2;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx2; cy=sy2

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s19     t2 from s19

            sx1=s19x1;sx2=s19x2;sx3=s19x3
            sy1=s19y1;sy2=s19y2;sy3=s19y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,18)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            18 format('zx18 = ',e12.4,e12.4,'j       magnitude =',e12.4)

        !**********************************************z_x19**********************************************

            aln=sqrt((0.05)**2+(0.05)**2)

        !   t1 from s19     t2 from s19

            sx1=s19x1;sx2=s19x2;sx3=s19x3
            sy1=s19y1;sy2=s19y2;sy3=s19y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx1; cy=sy1;

            sp11=scalar()
            vp11x=areavector_x()
            vp11y=areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx1; cy=sy1

            sp21=scalar()
            vp21x=areavector_x()
            vp21y=areavector_y()

        !   t1 from s20     t2 from s20

            sx1=s20x1;sx2=s20x2;sx3=s20x3
            sy1=s20y1;sy2=s20y2;sy3=s20y3

            rfpx=rfpx_(tc); rfpy=rfpy_(tc); cx=sx3; cy=sy3

            sp12=-scalar()
            vp12x=-areavector_x()
            vp12y=-areavector_y()

            rfpx=rfpx_(tc+1); rfpy=rfpy_(tc+1); cx=sx3; cy=sy3

            sp22=-scalar()
            vp22x=-areavector_x()
            vp22y=-areavector_y()

            z=alm*(j*omega*0.5*(ctx*(vp11x+vp12x+vp21x+vp22x)+cty*(vp11y+vp12y+vp21y+vp22y))+sp21+sp22-sp11-sp12)

            write(*,19)real(z),aimag(z),sqrt((real(z))**2+(aimag(z))**2)

            19 format('zx19 = ',e12.4,e12.4,'j       magnitude =',e12.4,/)

    end function

end module





