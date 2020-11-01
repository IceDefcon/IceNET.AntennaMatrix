module functions

    use definitions

contains

    real function func1()

        if(ir.eq.1) func1=real(A())
	
	return

    end function

    real function func2()

	if(ir.eq.1)func2=aimag(A())

        return

    end function


    complex function A()
	
        complex G, A1

        real eta3,rDashX,rDashY,rDashZ,R

        integer*2 i

        A1=cmplx(0.0,0.0)

        do i=1,7

            eta3=1.-eta1(i)-eta2(i)

            rDashX=Sx1*eta1(i)+Sx2*eta2(i)+Sx3*eta3
            rDashY=Sy1*eta1(i)+Sy2*eta2(i)+Sy3*eta3
            rDashZ=Sz1*eta1(i)+Sz2*eta2(i)+Sz3*eta3
		
            R=sqrt((rfpx-rDashX)*(rfpx-rDashX)+&
                (rfpy-rDashY)*(rfpy-rDashY)+&
                (rfpz-rDashZ)*(rfpz-rDashZ))	
            G=cmplx(0.0,-k*R)
            G=cexp(G)/R

            ! Calculate the coefficients for the vector potential. Acorner defines the originating corner of basis ftn and rCoeff the coefficient of r.
			
            A1=A1+G*eta1(i)*w(i)
            print *,'Ray:',R,'Function:',G,'Vector potential:',A1
        enddo

        A=A1

    end function

end module
