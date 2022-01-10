module Func

    use Definitions

contains

    real function scalarPotential()
	
        complex G

        real psi,rDashX,rDashY,rDashZ,R

        integer*2 i

        G=cmplx(0.0,0.0)

        do i=1,7

            psi=1.-eta(i)-neta(i)

            rDashX=Sx1*eta(i)+Sx2*neta(i)+Sx3*psi
            rDashY=Sy1*eta(i)+Sy2*neta(i)+Sy3*psi
            rDashZ=Sz1*eta(i)+Sz2*neta(i)+Sz3*psi
		
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

        real psi,rDashX,rDashY,rDashZ,R

        integer*2 i

        G=cmplx(0.0,0.0)

        do i=1,7

            psi=1.-eta(i)-neta(i)

            rDashX=Sx1*eta(i)+Sx2*neta(i)+Sx3*psi
            rDashY=Sy1*eta(i)+Sy2*neta(i)+Sy3*psi
            rDashZ=Sz1*eta(i)+Sz2*neta(i)+Sz3*psi

            R=sqrt((rfpx-rDashX)*(rfpx-rDashX)+&
                (rfpy-rDashY)*(rfpy-rDashY)+&
                (rfpz-rDashZ)*(rfpz-rDashZ))

            G=G+eta(i)*w(i)*cexp(cmplx(0.0,-k*R))/R

        enddo

        if(ir==1) vectorPotentialeta=real(G)
        if(ir==2) vectorPotentialeta=Aimag(G)

    end function

    real function vectorPotentialneta()

        complex G

        real psi,rDashX,rDashY,rDashZ,R

        integer*2 i

        G=cmplx(0.0,0.0)

        do i=1,7

            psi=1.-eta(i)-neta(i)

            rDashX=Sx1*eta(i)+Sx2*neta(i)+Sx3*psi
            rDashY=Sy1*eta(i)+Sy2*neta(i)+Sy3*psi
            rDashZ=Sz1*eta(i)+Sz2*neta(i)+Sz3*psi

            R=sqrt((rfpx-rDashX)*(rfpx-rDashX)+&
                (rfpy-rDashY)*(rfpy-rDashY)+&
                (rfpz-rDashZ)*(rfpz-rDashZ))

            G=G+neta(i)*w(i)*cexp(cmplx(0.0,-k*R))/R

        enddo

        if(ir==1) vectorPotentialneta=real(G)
        if(ir==2) vectorPotentialneta=Aimag(G)

    end function


end module
