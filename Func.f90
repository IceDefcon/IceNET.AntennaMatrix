module Func

    use Definitions

contains

    real function scalarPotential()
	
        complex G

        real Zeta,rDashX,rDashY,rDashZ,R

        integer*2 i

        G=cmplx(0.0,0.0)

        do i=1,7

            Zeta=1.-Xi(i)-Eta(i)

            rDashX=Sx1*Xi(i)+Sx2*Eta(i)+Sx3*Zeta
            rDashY=Sy1*Xi(i)+Sy2*Eta(i)+Sy3*Zeta
            rDashZ=Sz1*Xi(i)+Sz2*Eta(i)+Sz3*Zeta
		
            R=sqrt((rfpx-rDashX)*(rfpx-rDashX)+&
                (rfpy-rDashY)*(rfpy-rDashY)+&
                (rfpz-rDashZ)*(rfpz-rDashZ))
		
            G=G+w(i)*cexp(cmplx(0.0,-k*R))/R

            ! Calculate the coefficients for the vector potential. Acorner defines the originating corner of basis ftn and rCoeff the coefficient of r.

        enddo

        if(ir==1) scalarPotential=real(G)
        if(ir==2) scalarPotential=Aimag(G)

    end function

    real function vectorPotentialeta()

        complex G

        real Zeta,rDashX,rDashY,rDashZ,R

        integer*2 i

        G=cmplx(0.0,0.0)

        do i=1,7

            Zeta=1.-Xi(i)-Eta(i)

            rDashX=Sx1*Xi(i)+Sx2*Eta(i)+Sx3*Zeta
            rDashY=Sy1*Xi(i)+Sy2*Eta(i)+Sy3*Zeta
            rDashZ=Sz1*Xi(i)+Sz2*Eta(i)+Sz3*Zeta

            R=sqrt((rfpx-rDashX)*(rfpx-rDashX)+&
                (rfpy-rDashY)*(rfpy-rDashY)+&
                (rfpz-rDashZ)*(rfpz-rDashZ))

            G=G+Xi(i)*w(i)*cexp(cmplx(0.0,-k*R))/R

        enddo

        if(ir==1) vectorPotentialeta=real(G)
        if(ir==2) vectorPotentialeta=Aimag(G)

    end function

    real function vectorPotentialneta()

        complex G

        real Zeta,rDashX,rDashY,rDashZ,R

        integer*2 i

        G=cmplx(0.0,0.0)

        do i=1,7

            Zeta=1.-Xi(i)-Eta(i)

            rDashX=Sx1*Xi(i)+Sx2*Eta(i)+Sx3*Zeta
            rDashY=Sy1*Xi(i)+Sy2*Eta(i)+Sy3*Zeta
            rDashZ=Sz1*Xi(i)+Sz2*Eta(i)+Sz3*Zeta

            R=sqrt((rfpx-rDashX)*(rfpx-rDashX)+&
                (rfpy-rDashY)*(rfpy-rDashY)+&
                (rfpz-rDashZ)*(rfpz-rDashZ))

            G=G+Eta(i)*w(i)*cexp(cmplx(0.0,-k*R))/R

        enddo

        if(ir==1) vectorPotentialneta=real(G)
        if(ir==2) vectorPotentialneta=Aimag(G)

    end function


end module
