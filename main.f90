program prak1

real :: LEFT_BORDER, RIGHT_BORDER, h, sigma, eps=0.000000001
integer :: N, n0, m, i, k
real, allocatable, dimension (:) :: X, Y, polinom, gradient, V, tmp
logical :: first

    write (*, '(A)') 'Hello!'
    write (*, '(A)') 'Please, enter values:'
    write (*, '(A)') ''
    write (*, '(A)', ADVANCE='NO') 'c = '
    read *, LEFT_BORDER
    write (*, '(A)', ADVANCE='NO') 'd = '
    read *, RIGHT_BORDER
    write (*, '(A)', ADVANCE='NO') 'N = '
    read *, N
    write (*, '(A)', ADVANCE='NO') 'n0 = '
    read *, n0
    write (*, '(A)', ADVANCE='NO') 'm = '
    read *, m
    write (*, '(A)') ''

    h=(RIGHT_BORDER-LEFT_BORDER)/N
    allocate (X(N+1))

    do i=1, N+1
        X(i)=LEFT_BORDER+(i-1)*h
    end do

    allocate (Y(N+1))
    write (*, '(A)') 'Operate with function: y = cos(2x) - ch(x/5) '
    write (*, '(A)') ''

    call get_target_F_values (X, Y, N+1)

    allocate (polinom(n0+m))
    allocate (gradient(n0+m))
    allocate (V(n0+m))
    allocate (tmp(n0+m))

    do i=1, size(polinom)
        polinom(i)=0
        gradient(i)=0
        V(i)=0
        tmp(i)=0
    end do

   do k=n0, n0+m-1
     i=0
     first=.TRUE.
     sigma=1000

     write (*, '(A8, I1)') 'for n = ', k
     write (*, '(A)' ) ''

        do while (sigma>=eps)
            tmp=polinom
            i=i+1

            call get_V()

            polinom=polinom+lambda()*V
            if (first.eqv..true.) then
                first=.false.
            else
                sigma=get_sigma()
            end if

            write(*,'(A7, I3, A14, F10.5)')'iter = ', i, ': |grad(F)| = ', sqrt(dot_product(gradient,gradient))

        end do

        write (*, '(A)') ''

        write (*, '(A4, A9, A11, A11, A11)') 'iter', 'X(i)   ', 'Y(i)    ','P(i)    ','d    '

        do i=1, N+1
            write (*, '(I2, F10.4, F10.4, F10.4, F10.4)') i-1, X(i), Y(i), P(X(i), polinom), abs(Y(i)-P(X(i), polinom))
        end do

        write (*, '(A)') ''
        write (*, '(A)') 'polynomial coefficients:'
        do i=0, k
            write (*, '(A1, I1, A3, F9.6)') 'a', i, ' = ', polinom(i+1)
        end do
        write (*, '(A)') ''

   end do

   deallocate (X)
   deallocate (Y)
   deallocate (polinom)
   deallocate (V)
   deallocate (gradient)
   deallocate (tmp)

   write (*, '(A)') ''

contains

    real function P (point, array)

        real :: point
        integer :: i
        real, dimension (:) :: array

        P=0

        do i=0, k
            P= array(k+1-i)+P*point
        end do

    end function P

    subroutine get_gradient ()
        integer :: i, j

        do j=1, k+1
            gradient(j)=0
            do i=1, N+1
                gradient(j)=gradient(j)+2*(X(i)**(j-1))*(P(X(i), polinom)-Y(i))
            end do
        end do

    end subroutine get_gradient

    real function G (array)
        real, dimension (:) :: array
        integer :: i

        G=0

        do i=1, N+1
            G=G+(Y(i)-P(X(i), array))**2
        end do

    end function G

    subroutine get_V ()
        real :: tmp

        if (first.eqv..TRUE.) then
                call get_gradient ()
                V=(-1)*gradient
            else
                tmp=dot_product(gradient, gradient)
                call get_gradient ()
                V=(-1)*gradient+(dot_product(gradient, gradient)/tmp)*V
        end if

    end subroutine get_V

    real function lambda ()
        real :: m1, m2, m3

        m1=G(polinom-V)
        m2=G(polinom)
        m3=G(polinom+V)

        lambda=(m1-m3)/(2*(m1-2*m2+m3))

    end function lambda

    real function get_sigma ()
        integer :: i
        real :: q

        tmp=tmp-polinom
        get_sigma=1000

        do i=1, k+1
            if (polinom(i)/=0) then
                q=abs(tmp(i)/polinom(i))
                if (q<get_sigma) then
                    get_sigma=q
                end if
            end if
        end do

    end function get_sigma

end program prak1

subroutine my_print (X, N)

    integer :: N, i
    real, dimension (N) :: X

    do i=1, N
        print *, X(i)
    end do

end subroutine my_print

subroutine get_target_F_values (X, Y, N)

    integer :: N, i
    real, dimension (N) :: X, Y

    do i=1, N
        Y(i)=cos(2*X(i))-((exp(X(i)/5)+exp(-X(i)/5))/2)
    end do

end subroutine get_target_F_values
