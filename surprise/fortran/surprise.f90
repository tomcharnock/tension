module tension_mod
implicit none
integer, parameter :: DP = kind(1.0D0)
real(DP), parameter :: pi = 4.0d0*atan(1.0d0)
contains

        subroutine read_mean_cov(filename,means,covs)
        implicit none
        character (len=50) :: filename,param
        real(DP) :: mean,means(5),covs(5,5)
        real(DP), allocatable :: cov(:,:)
        integer :: i,j,num,covnum(5),io
        
        open(unit=501,file='data/'//trim(filename)//'.margestats',action='read')

        do i=1,3
                read(501,*)
        end do        
        num = 0
        do 
                read(501,*,iostat=io) param,mean
                num = num + 1
                if (io<0) exit
                if (io>0) stop 'Problem Reading'
                if (param.eq.'omegabh2') then 
                        means(1) = mean
                        covnum(1) = num
                end if
                if (param.eq.'omegach2') then
                        means(2) = mean
                        covnum(2) = num
                end if
                if (param.eq.'theta') then 
                        means(3) = mean
                        covnum(3) = num
                end if
                if (param.eq.'logA') then 
                        means(4) = mean
                        covnum(4) = num
                end if
                if (param.eq.'ns') then 
                        means(5) = mean
                        covnum(5) = num
                end if
        end do
        close(unit=501)

        allocate(cov(covnum(5),covnum(5)))

        open(unit=501,file='data/'//trim(filename)//'.covmat',action='read')
        read(501,*)
        do i=1,covnum(5)
                read(501,*) cov(:,i)
        end do
        close(unit=501)
        
        do i=1,5
                do j=1,5
                        covs(i,j) = cov(covnum(i),covnum(j))
                end do
        end do
        return
       
        end subroutine read_mean_cov
        
        subroutine inv(A, AINV)
        implicit none
        real(DP), intent(in) :: A(5,5)
        REAL(DP), intent(out) :: AINV(5,5)
        logical :: OK_FLAG
        real(DP), parameter :: EPS = 1.0D-40
        real(DP) :: DET, A11, A12, A13, A14, A15, A21, A22, A23, A24, &
           A25, A31, A32, A33, A34, A35, A41, A42, A43, A44, A45,   &
           A51, A52, A53, A54, A55
        real(DP) :: COFACTOR(5,5)

        A11=A(1,1); A12=A(1,2); A13=A(1,3); A14=A(1,4); A15=A(1,5)
        A21=A(2,1); A22=A(2,2); A23=A(2,3); A24=A(2,4); A25=A(2,5)
        A31=A(3,1); A32=A(3,2); A33=A(3,3); A34=A(3,4); A35=A(3,5)
        A41=A(4,1); A42=A(4,2); A43=A(4,3); A44=A(4,4); A45=A(4,5)
        A51=A(5,1); A52=A(5,2); A53=A(5,3); A54=A(5,4); A55=A(5,5)

        DET = A15*A24*A33*A42*A51-A14*A25*A33*A42*A51-A15*A23*A34*A42*A51+    &
           A13*A25*A34*A42*A51+A14*A23*A35*A42*A51-A13*A24*A35*A42*A51-       &
           A15*A24*A32*A43*A51+A14*A25*A32*A43*A51+A15*A22*A34*A43*A51-       &
           A12*A25*A34*A43*A51-A14*A22*A35*A43*A51+A12*A24*A35*A43*A51+       &
           A15*A23*A32*A44*A51-A13*A25*A32*A44*A51-A15*A22*A33*A44*A51+       &
           A12*A25*A33*A44*A51+A13*A22*A35*A44*A51-A12*A23*A35*A44*A51-       &
           A14*A23*A32*A45*A51+A13*A24*A32*A45*A51+A14*A22*A33*A45*A51-       &
           A12*A24*A33*A45*A51-A13*A22*A34*A45*A51+A12*A23*A34*A45*A51-       &
           A15*A24*A33*A41*A52+A14*A25*A33*A41*A52+A15*A23*A34*A41*A52-       &
           A13*A25*A34*A41*A52-A14*A23*A35*A41*A52+A13*A24*A35*A41*A52+       &
           A15*A24*A31*A43*A52-A14*A25*A31*A43*A52-A15*A21*A34*A43*A52+       &
           A11*A25*A34*A43*A52+A14*A21*A35*A43*A52-A11*A24*A35*A43*A52-       &
           A15*A23*A31*A44*A52+A13*A25*A31*A44*A52+A15*A21*A33*A44*A52-       &
           A11*A25*A33*A44*A52-A13*A21*A35*A44*A52+A11*A23*A35*A44*A52+       &
           A14*A23*A31*A45*A52-A13*A24*A31*A45*A52-A14*A21*A33*A45*A52+       &
           A11*A24*A33*A45*A52+A13*A21*A34*A45*A52-A11*A23*A34*A45*A52+       &
           A15*A24*A32*A41*A53-A14*A25*A32*A41*A53-A15*A22*A34*A41*A53+       &
           A12*A25*A34*A41*A53+A14*A22*A35*A41*A53-A12*A24*A35*A41*A53-       &
           A15*A24*A31*A42*A53+A14*A25*A31*A42*A53+A15*A21*A34*A42*A53-       &
           A11*A25*A34*A42*A53-A14*A21*A35*A42*A53+A11*A24*A35*A42*A53+       &
           A15*A22*A31*A44*A53-A12*A25*A31*A44*A53-A15*A21*A32*A44*A53+       &
           A11*A25*A32*A44*A53+A12*A21*A35*A44*A53-A11*A22*A35*A44*A53-       &
           A14*A22*A31*A45*A53+A12*A24*A31*A45*A53+A14*A21*A32*A45*A53-       &
           A11*A24*A32*A45*A53-A12*A21*A34*A45*A53+A11*A22*A34*A45*A53-       &
           A15*A23*A32*A41*A54+A13*A25*A32*A41*A54+A15*A22*A33*A41*A54-       &
           A12*A25*A33*A41*A54-A13*A22*A35*A41*A54+A12*A23*A35*A41*A54+       &
           A15*A23*A31*A42*A54-A13*A25*A31*A42*A54-A15*A21*A33*A42*A54+       &
           A11*A25*A33*A42*A54+A13*A21*A35*A42*A54-A11*A23*A35*A42*A54-       &
           A15*A22*A31*A43*A54+A12*A25*A31*A43*A54+A15*A21*A32*A43*A54-       &
           A11*A25*A32*A43*A54-A12*A21*A35*A43*A54+A11*A22*A35*A43*A54+       &
           A13*A22*A31*A45*A54-A12*A23*A31*A45*A54-A13*A21*A32*A45*A54+       &
           A11*A23*A32*A45*A54+A12*A21*A33*A45*A54-A11*A22*A33*A45*A54+       &
           A14*A23*A32*A41*A55-A13*A24*A32*A41*A55-A14*A22*A33*A41*A55+       &
           A12*A24*A33*A41*A55+A13*A22*A34*A41*A55-A12*A23*A34*A41*A55-       &
           A14*A23*A31*A42*A55+A13*A24*A31*A42*A55+A14*A21*A33*A42*A55-       &
           A11*A24*A33*A42*A55-A13*A21*A34*A42*A55+A11*A23*A34*A42*A55+       &
           A14*A22*A31*A43*A55-A12*A24*A31*A43*A55-A14*A21*A32*A43*A55+       &
           A11*A24*A32*A43*A55+A12*A21*A34*A43*A55-A11*A22*A34*A43*A55-       &
           A13*A22*A31*A44*A55+A12*A23*A31*A44*A55+A13*A21*A32*A44*A55-       &
           A11*A23*A32*A44*A55-A12*A21*A33*A44*A55+A11*A22*A33*A44*A55

        IF (ABS(DET) .LE. EPS) THEN
           AINV = 0.0D0
           OK_FLAG = .FALSE.
           RETURN
        END IF

        COFACTOR(1,1) = A25*A34*A43*A52-A24*A35*A43*A52-A25*A33*A44*A52+      &
           A23*A35*A44*A52+A24*A33*A45*A52-A23*A34*A45*A52-A25*A34*A42*A53+   &
           A24*A35*A42*A53+A25*A32*A44*A53-A22*A35*A44*A53-A24*A32*A45*A53+   &
           A22*A34*A45*A53+A25*A33*A42*A54-A23*A35*A42*A54-A25*A32*A43*A54+   &
           A22*A35*A43*A54+A23*A32*A45*A54-A22*A33*A45*A54-A24*A33*A42*A55+   &
           A23*A34*A42*A55+A24*A32*A43*A55-A22*A34*A43*A55-A23*A32*A44*A55+   &
           A22*A33*A44*A55

        COFACTOR(2,1) = -A15*A34*A43*A52+A14*A35*A43*A52+A15*A33*A44*A52-     &
           A13*A35*A44*A52-A14*A33*A45*A52+A13*A34*A45*A52+A15*A34*A42*A53-   &
           A14*A35*A42*A53-A15*A32*A44*A53+A12*A35*A44*A53+A14*A32*A45*A53-   &
           A12*A34*A45*A53-A15*A33*A42*A54+A13*A35*A42*A54+A15*A32*A43*A54-   &
           A12*A35*A43*A54-A13*A32*A45*A54+A12*A33*A45*A54+A14*A33*A42*A55-   &
           A13*A34*A42*A55-A14*A32*A43*A55+A12*A34*A43*A55+A13*A32*A44*A55-   &
           A12*A33*A44*A55

        COFACTOR(3,1) = A15*A24*A43*A52-A14*A25*A43*A52-A15*A23*A44*A52+      &
           A13*A25*A44*A52+A14*A23*A45*A52-A13*A24*A45*A52-A15*A24*A42*A53+   &
           A14*A25*A42*A53+A15*A22*A44*A53-A12*A25*A44*A53-A14*A22*A45*A53+   &
           A12*A24*A45*A53+A15*A23*A42*A54-A13*A25*A42*A54-A15*A22*A43*A54+   &
           A12*A25*A43*A54+A13*A22*A45*A54-A12*A23*A45*A54-A14*A23*A42*A55+   &
           A13*A24*A42*A55+A14*A22*A43*A55-A12*A24*A43*A55-A13*A22*A44*A55+   &
           A12*A23*A44*A55

        COFACTOR(4,1) = -A15*A24*A33*A52+A14*A25*A33*A52+A15*A23*A34*A52-     &
           A13*A25*A34*A52-A14*A23*A35*A52+A13*A24*A35*A52+A15*A24*A32*A53-   &
           A14*A25*A32*A53-A15*A22*A34*A53+A12*A25*A34*A53+A14*A22*A35*A53-   &
           A12*A24*A35*A53-A15*A23*A32*A54+A13*A25*A32*A54+A15*A22*A33*A54-   &
           A12*A25*A33*A54-A13*A22*A35*A54+A12*A23*A35*A54+A14*A23*A32*A55-   &
           A13*A24*A32*A55-A14*A22*A33*A55+A12*A24*A33*A55+A13*A22*A34*A55-   &
           A12*A23*A34*A55

        COFACTOR(5,1) = A15*A24*A33*A42-A14*A25*A33*A42-A15*A23*A34*A42+      &
           A13*A25*A34*A42+A14*A23*A35*A42-A13*A24*A35*A42-A15*A24*A32*A43+   &
           A14*A25*A32*A43+A15*A22*A34*A43-A12*A25*A34*A43-A14*A22*A35*A43+   &
           A12*A24*A35*A43+A15*A23*A32*A44-A13*A25*A32*A44-A15*A22*A33*A44+   &
           A12*A25*A33*A44+A13*A22*A35*A44-A12*A23*A35*A44-A14*A23*A32*A45+   &
           A13*A24*A32*A45+A14*A22*A33*A45-A12*A24*A33*A45-A13*A22*A34*A45+   &
           A12*A23*A34*A45

        COFACTOR(1,2) = -A25*A34*A43*A51+A24*A35*A43*A51+A25*A33*A44*A51-     &
           A23*A35*A44*A51-A24*A33*A45*A51+A23*A34*A45*A51+A25*A34*A41*A53-   &
           A24*A35*A41*A53-A25*A31*A44*A53+A21*A35*A44*A53+A24*A31*A45*A53-   &
           A21*A34*A45*A53-A25*A33*A41*A54+A23*A35*A41*A54+A25*A31*A43*A54-   &
           A21*A35*A43*A54-A23*A31*A45*A54+A21*A33*A45*A54+A24*A33*A41*A55-   &
           A23*A34*A41*A55-A24*A31*A43*A55+A21*A34*A43*A55+A23*A31*A44*A55-   &
           A21*A33*A44*A55

        COFACTOR(2,2) = A15*A34*A43*A51-A14*A35*A43*A51-A15*A33*A44*A51+      &
           A13*A35*A44*A51+A14*A33*A45*A51-A13*A34*A45*A51-A15*A34*A41*A53+   &
           A14*A35*A41*A53+A15*A31*A44*A53-A11*A35*A44*A53-A14*A31*A45*A53+   &
           A11*A34*A45*A53+A15*A33*A41*A54-A13*A35*A41*A54-A15*A31*A43*A54+   &
           A11*A35*A43*A54+A13*A31*A45*A54-A11*A33*A45*A54-A14*A33*A41*A55+   &
           A13*A34*A41*A55+A14*A31*A43*A55-A11*A34*A43*A55-A13*A31*A44*A55+   &
           A11*A33*A44*A55

        COFACTOR(3,2) = -A15*A24*A43*A51+A14*A25*A43*A51+A15*A23*A44*A51-     &
           A13*A25*A44*A51-A14*A23*A45*A51+A13*A24*A45*A51+A15*A24*A41*A53-   &
           A14*A25*A41*A53-A15*A21*A44*A53+A11*A25*A44*A53+A14*A21*A45*A53-   &
           A11*A24*A45*A53-A15*A23*A41*A54+A13*A25*A41*A54+A15*A21*A43*A54-   &
           A11*A25*A43*A54-A13*A21*A45*A54+A11*A23*A45*A54+A14*A23*A41*A55-   &
           A13*A24*A41*A55-A14*A21*A43*A55+A11*A24*A43*A55+A13*A21*A44*A55-   &
           A11*A23*A44*A55

        COFACTOR(4,2) = A15*A24*A33*A51-A14*A25*A33*A51-A15*A23*A34*A51+      &
           A13*A25*A34*A51+A14*A23*A35*A51-A13*A24*A35*A51-A15*A24*A31*A53+   &
           A14*A25*A31*A53+A15*A21*A34*A53-A11*A25*A34*A53-A14*A21*A35*A53+   &
           A11*A24*A35*A53+A15*A23*A31*A54-A13*A25*A31*A54-A15*A21*A33*A54+   &
           A11*A25*A33*A54+A13*A21*A35*A54-A11*A23*A35*A54-A14*A23*A31*A55+   &
           A13*A24*A31*A55+A14*A21*A33*A55-A11*A24*A33*A55-A13*A21*A34*A55+   &
           A11*A23*A34*A55

        COFACTOR(5,2) = -A15*A24*A33*A41+A14*A25*A33*A41+A15*A23*A34*A41-     &
           A13*A25*A34*A41-A14*A23*A35*A41+A13*A24*A35*A41+A15*A24*A31*A43-   &
           A14*A25*A31*A43-A15*A21*A34*A43+A11*A25*A34*A43+A14*A21*A35*A43-   &
           A11*A24*A35*A43-A15*A23*A31*A44+A13*A25*A31*A44+A15*A21*A33*A44-   &
           A11*A25*A33*A44-A13*A21*A35*A44+A11*A23*A35*A44+A14*A23*A31*A45-   &
           A13*A24*A31*A45-A14*A21*A33*A45+A11*A24*A33*A45+A13*A21*A34*A45-   &
           A11*A23*A34*A45

        COFACTOR(1,3) = A25*A34*A42*A51-A24*A35*A42*A51-A25*A32*A44*A51+      &
           A22*A35*A44*A51+A24*A32*A45*A51-A22*A34*A45*A51-A25*A34*A41*A52+   &
           A24*A35*A41*A52+A25*A31*A44*A52-A21*A35*A44*A52-A24*A31*A45*A52+   &
           A21*A34*A45*A52+A25*A32*A41*A54-A22*A35*A41*A54-A25*A31*A42*A54+   &
           A21*A35*A42*A54+A22*A31*A45*A54-A21*A32*A45*A54-A24*A32*A41*A55+   &
           A22*A34*A41*A55+A24*A31*A42*A55-A21*A34*A42*A55-A22*A31*A44*A55+   &
           A21*A32*A44*A55

        COFACTOR(2,3) = -A15*A34*A42*A51+A14*A35*A42*A51+A15*A32*A44*A51-     &
           A12*A35*A44*A51-A14*A32*A45*A51+A12*A34*A45*A51+A15*A34*A41*A52-   &
           A14*A35*A41*A52-A15*A31*A44*A52+A11*A35*A44*A52+A14*A31*A45*A52-   &
           A11*A34*A45*A52-A15*A32*A41*A54+A12*A35*A41*A54+A15*A31*A42*A54-   &
           A11*A35*A42*A54-A12*A31*A45*A54+A11*A32*A45*A54+A14*A32*A41*A55-   &
           A12*A34*A41*A55-A14*A31*A42*A55+A11*A34*A42*A55+A12*A31*A44*A55-   &
           A11*A32*A44*A55

        COFACTOR(3,3) = A15*A24*A42*A51-A14*A25*A42*A51-A15*A22*A44*A51+      &
           A12*A25*A44*A51+A14*A22*A45*A51-A12*A24*A45*A51-A15*A24*A41*A52+   &
           A14*A25*A41*A52+A15*A21*A44*A52-A11*A25*A44*A52-A14*A21*A45*A52+   &
           A11*A24*A45*A52+A15*A22*A41*A54-A12*A25*A41*A54-A15*A21*A42*A54+   &
           A11*A25*A42*A54+A12*A21*A45*A54-A11*A22*A45*A54-A14*A22*A41*A55+   &
           A12*A24*A41*A55+A14*A21*A42*A55-A11*A24*A42*A55-A12*A21*A44*A55+   &
           A11*A22*A44*A55

        COFACTOR(4,3) = -A15*A24*A32*A51+A14*A25*A32*A51+A15*A22*A34*A51-     &
           A12*A25*A34*A51-A14*A22*A35*A51+A12*A24*A35*A51+A15*A24*A31*A52-   &
           A14*A25*A31*A52-A15*A21*A34*A52+A11*A25*A34*A52+A14*A21*A35*A52-   &
           A11*A24*A35*A52-A15*A22*A31*A54+A12*A25*A31*A54+A15*A21*A32*A54-   &
           A11*A25*A32*A54-A12*A21*A35*A54+A11*A22*A35*A54+A14*A22*A31*A55-   &
           A12*A24*A31*A55-A14*A21*A32*A55+A11*A24*A32*A55+A12*A21*A34*A55-   &
           A11*A22*A34*A55

        COFACTOR(5,3) = A15*A24*A32*A41-A14*A25*A32*A41-A15*A22*A34*A41+      &
           A12*A25*A34*A41+A14*A22*A35*A41-A12*A24*A35*A41-A15*A24*A31*A42+   &
           A14*A25*A31*A42+A15*A21*A34*A42-A11*A25*A34*A42-A14*A21*A35*A42+   &
           A11*A24*A35*A42+A15*A22*A31*A44-A12*A25*A31*A44-A15*A21*A32*A44+   &
           A11*A25*A32*A44+A12*A21*A35*A44-A11*A22*A35*A44-A14*A22*A31*A45+   &
           A12*A24*A31*A45+A14*A21*A32*A45-A11*A24*A32*A45-A12*A21*A34*A45+   &
           A11*A22*A34*A45

        COFACTOR(1,4) = -A25*A33*A42*A51+A23*A35*A42*A51+A25*A32*A43*A51-     &
           A22*A35*A43*A51-A23*A32*A45*A51+A22*A33*A45*A51+A25*A33*A41*A52-   &
           A23*A35*A41*A52-A25*A31*A43*A52+A21*A35*A43*A52+A23*A31*A45*A52-   &
           A21*A33*A45*A52-A25*A32*A41*A53+A22*A35*A41*A53+A25*A31*A42*A53-   &
           A21*A35*A42*A53-A22*A31*A45*A53+A21*A32*A45*A53+A23*A32*A41*A55-   &
           A22*A33*A41*A55-A23*A31*A42*A55+A21*A33*A42*A55+A22*A31*A43*A55-   &
           A21*A32*A43*A55

        COFACTOR(2,4) = A15*A33*A42*A51-A13*A35*A42*A51-A15*A32*A43*A51+      &
           A12*A35*A43*A51+A13*A32*A45*A51-A12*A33*A45*A51-A15*A33*A41*A52+   &
           A13*A35*A41*A52+A15*A31*A43*A52-A11*A35*A43*A52-A13*A31*A45*A52+   &
           A11*A33*A45*A52+A15*A32*A41*A53-A12*A35*A41*A53-A15*A31*A42*A53+   &
           A11*A35*A42*A53+A12*A31*A45*A53-A11*A32*A45*A53-A13*A32*A41*A55+   &
           A12*A33*A41*A55+A13*A31*A42*A55-A11*A33*A42*A55-A12*A31*A43*A55+   &
           A11*A32*A43*A55

        COFACTOR(3,4) = -A15*A23*A42*A51+A13*A25*A42*A51+A15*A22*A43*A51-     &
           A12*A25*A43*A51-A13*A22*A45*A51+A12*A23*A45*A51+A15*A23*A41*A52-   &
           A13*A25*A41*A52-A15*A21*A43*A52+A11*A25*A43*A52+A13*A21*A45*A52-   &
           A11*A23*A45*A52-A15*A22*A41*A53+A12*A25*A41*A53+A15*A21*A42*A53-   &
           A11*A25*A42*A53-A12*A21*A45*A53+A11*A22*A45*A53+A13*A22*A41*A55-   &
           A12*A23*A41*A55-A13*A21*A42*A55+A11*A23*A42*A55+A12*A21*A43*A55-   &
           A11*A22*A43*A55

        COFACTOR(4,4) = A15*A23*A32*A51-A13*A25*A32*A51-A15*A22*A33*A51+      &
           A12*A25*A33*A51+A13*A22*A35*A51-A12*A23*A35*A51-A15*A23*A31*A52+   &
           A13*A25*A31*A52+A15*A21*A33*A52-A11*A25*A33*A52-A13*A21*A35*A52+   &
           A11*A23*A35*A52+A15*A22*A31*A53-A12*A25*A31*A53-A15*A21*A32*A53+   &
           A11*A25*A32*A53+A12*A21*A35*A53-A11*A22*A35*A53-A13*A22*A31*A55+   &
           A12*A23*A31*A55+A13*A21*A32*A55-A11*A23*A32*A55-A12*A21*A33*A55+   &
           A11*A22*A33*A55

        COFACTOR(5,4) = -A15*A23*A32*A41+A13*A25*A32*A41+A15*A22*A33*A41-     &
           A12*A25*A33*A41-A13*A22*A35*A41+A12*A23*A35*A41+A15*A23*A31*A42-   &
           A13*A25*A31*A42-A15*A21*A33*A42+A11*A25*A33*A42+A13*A21*A35*A42-   &
           A11*A23*A35*A42-A15*A22*A31*A43+A12*A25*A31*A43+A15*A21*A32*A43-   &
           A11*A25*A32*A43-A12*A21*A35*A43+A11*A22*A35*A43+A13*A22*A31*A45-   &
           A12*A23*A31*A45-A13*A21*A32*A45+A11*A23*A32*A45+A12*A21*A33*A45-   &
           A11*A22*A33*A45

        COFACTOR(1,5) = A24*A33*A42*A51-A23*A34*A42*A51-A24*A32*A43*A51+      &
           A22*A34*A43*A51+A23*A32*A44*A51-A22*A33*A44*A51-A24*A33*A41*A52+   &
           A23*A34*A41*A52+A24*A31*A43*A52-A21*A34*A43*A52-A23*A31*A44*A52+   &
           A21*A33*A44*A52+A24*A32*A41*A53-A22*A34*A41*A53-A24*A31*A42*A53+   &
           A21*A34*A42*A53+A22*A31*A44*A53-A21*A32*A44*A53-A23*A32*A41*A54+   &
           A22*A33*A41*A54+A23*A31*A42*A54-A21*A33*A42*A54-A22*A31*A43*A54+   &
           A21*A32*A43*A54

        COFACTOR(2,5) = -A14*A33*A42*A51+A13*A34*A42*A51+A14*A32*A43*A51-     &
           A12*A34*A43*A51-A13*A32*A44*A51+A12*A33*A44*A51+A14*A33*A41*A52-   &
           A13*A34*A41*A52-A14*A31*A43*A52+A11*A34*A43*A52+A13*A31*A44*A52-   &
           A11*A33*A44*A52-A14*A32*A41*A53+A12*A34*A41*A53+A14*A31*A42*A53-   &
           A11*A34*A42*A53-A12*A31*A44*A53+A11*A32*A44*A53+A13*A32*A41*A54-   &
           A12*A33*A41*A54-A13*A31*A42*A54+A11*A33*A42*A54+A12*A31*A43*A54-   &
           A11*A32*A43*A54

        COFACTOR(3,5) = A14*A23*A42*A51-A13*A24*A42*A51-A14*A22*A43*A51+      &
           A12*A24*A43*A51+A13*A22*A44*A51-A12*A23*A44*A51-A14*A23*A41*A52+   &
           A13*A24*A41*A52+A14*A21*A43*A52-A11*A24*A43*A52-A13*A21*A44*A52+   &
           A11*A23*A44*A52+A14*A22*A41*A53-A12*A24*A41*A53-A14*A21*A42*A53+   &
           A11*A24*A42*A53+A12*A21*A44*A53-A11*A22*A44*A53-A13*A22*A41*A54+   &
           A12*A23*A41*A54+A13*A21*A42*A54-A11*A23*A42*A54-A12*A21*A43*A54+   &
           A11*A22*A43*A54

        COFACTOR(4,5) = -A14*A23*A32*A51+A13*A24*A32*A51+A14*A22*A33*A51-     &
           A12*A24*A33*A51-A13*A22*A34*A51+A12*A23*A34*A51+A14*A23*A31*A52-   &
           A13*A24*A31*A52-A14*A21*A33*A52+A11*A24*A33*A52+A13*A21*A34*A52-   &
           A11*A23*A34*A52-A14*A22*A31*A53+A12*A24*A31*A53+A14*A21*A32*A53-   &
           A11*A24*A32*A53-A12*A21*A34*A53+A11*A22*A34*A53+A13*A22*A31*A54-   &
           A12*A23*A31*A54-A13*A21*A32*A54+A11*A23*A32*A54+A12*A21*A33*A54-   &
           A11*A22*A33*A54

        COFACTOR(5,5) = A14*A23*A32*A41-A13*A24*A32*A41-A14*A22*A33*A41+      &
           A12*A24*A33*A41+A13*A22*A34*A41-A12*A23*A34*A41-A14*A23*A31*A42+   &
           A13*A24*A31*A42+A14*A21*A33*A42-A11*A24*A33*A42-A13*A21*A34*A42+   &
           A11*A23*A34*A42+A14*A22*A31*A43-A12*A24*A31*A43-A14*A21*A32*A43+   &
           A11*A24*A32*A43+A12*A21*A34*A43-A11*A22*A34*A43-A13*A22*A31*A44+   &
           A12*A23*A31*A44+A13*A21*A32*A44-A11*A23*A32*A44-A12*A21*A33*A44+   &
           A11*A22*A33*A44

        AINV = TRANSPOSE(COFACTOR) / DET

        OK_FLAG = .true.

        return
        end subroutine inv

        function multivariate_normal(x,mean,norm,invcov)
        implicit none
        real(DP) :: multivariate_normal
        real(DP) :: x(5),mean(5),norm,invcov(5,5)

        multivariate_normal = norm*dexp(-0.5d0*dot_product((x-mean),matmul(invcov,(x-mean))))
        return
        end function multivariate_normal 

        function det(A)
        implicit none
        real(DP), intent(in) :: A(5,5)
        real(DP) :: det, A11, A12, A13, A14, A15, A21, A22, A23, A24, &
                    A25, A31, A32, A33, A34, A35, A41, A42, A43, A44, A45,   &
                    A51, A52, A53, A54, A55


        A11=A(1,1); A12=A(1,2); A13=A(1,3); A14=A(1,4); A15=A(1,5)
        A21=A(2,1); A22=A(2,2); A23=A(2,3); A24=A(2,4); A25=A(2,5)
        A31=A(3,1); A32=A(3,2); A33=A(3,3); A34=A(3,4); A35=A(3,5)
        A41=A(4,1); A42=A(4,2); A43=A(4,3); A44=A(4,4); A45=A(4,5)
        A51=A(5,1); A52=A(5,2); A53=A(5,3); A54=A(5,4); A55=A(5,5)

        det = A15*A24*A33*A42*A51-A14*A25*A33*A42*A51-A15*A23*A34*A42*A51+ &
              A13*A25*A34*A42*A51+A14*A23*A35*A42*A51-A13*A24*A35*A42*A51- &
              A15*A24*A32*A43*A51+A14*A25*A32*A43*A51+A15*A22*A34*A43*A51- &
              A12*A25*A34*A43*A51-A14*A22*A35*A43*A51+A12*A24*A35*A43*A51+ &
              A15*A23*A32*A44*A51-A13*A25*A32*A44*A51-A15*A22*A33*A44*A51+ &
              A12*A25*A33*A44*A51+A13*A22*A35*A44*A51-A12*A23*A35*A44*A51- &
              A14*A23*A32*A45*A51+A13*A24*A32*A45*A51+A14*A22*A33*A45*A51- &
              A12*A24*A33*A45*A51-A13*A22*A34*A45*A51+A12*A23*A34*A45*A51- &
              A15*A24*A33*A41*A52+A14*A25*A33*A41*A52+A15*A23*A34*A41*A52- &
              A13*A25*A34*A41*A52-A14*A23*A35*A41*A52+A13*A24*A35*A41*A52+ &
              A15*A24*A31*A43*A52-A14*A25*A31*A43*A52-A15*A21*A34*A43*A52+ &
              A11*A25*A34*A43*A52+A14*A21*A35*A43*A52-A11*A24*A35*A43*A52- &
              A15*A23*A31*A44*A52+A13*A25*A31*A44*A52+A15*A21*A33*A44*A52- &
              A11*A25*A33*A44*A52-A13*A21*A35*A44*A52+A11*A23*A35*A44*A52+ &
              A14*A23*A31*A45*A52-A13*A24*A31*A45*A52-A14*A21*A33*A45*A52+ &
              A11*A24*A33*A45*A52+A13*A21*A34*A45*A52-A11*A23*A34*A45*A52+ &
              A15*A24*A32*A41*A53-A14*A25*A32*A41*A53-A15*A22*A34*A41*A53+ &
              A12*A25*A34*A41*A53+A14*A22*A35*A41*A53-A12*A24*A35*A41*A53- &
              A15*A24*A31*A42*A53+A14*A25*A31*A42*A53+A15*A21*A34*A42*A53- &
              A11*A25*A34*A42*A53-A14*A21*A35*A42*A53+A11*A24*A35*A42*A53+ &
              A15*A22*A31*A44*A53-A12*A25*A31*A44*A53-A15*A21*A32*A44*A53+ &
              A11*A25*A32*A44*A53+A12*A21*A35*A44*A53-A11*A22*A35*A44*A53- &
              A14*A22*A31*A45*A53+A12*A24*A31*A45*A53+A14*A21*A32*A45*A53- &
              A11*A24*A32*A45*A53-A12*A21*A34*A45*A53+A11*A22*A34*A45*A53- &
              A15*A23*A32*A41*A54+A13*A25*A32*A41*A54+A15*A22*A33*A41*A54- &
              A12*A25*A33*A41*A54-A13*A22*A35*A41*A54+A12*A23*A35*A41*A54+ &
              A15*A23*A31*A42*A54-A13*A25*A31*A42*A54-A15*A21*A33*A42*A54+ &
              A11*A25*A33*A42*A54+A13*A21*A35*A42*A54-A11*A23*A35*A42*A54- &
              A15*A22*A31*A43*A54+A12*A25*A31*A43*A54+A15*A21*A32*A43*A54- &
              A11*A25*A32*A43*A54-A12*A21*A35*A43*A54+A11*A22*A35*A43*A54+ &
              A13*A22*A31*A45*A54-A12*A23*A31*A45*A54-A13*A21*A32*A45*A54+ &
              A11*A23*A32*A45*A54+A12*A21*A33*A45*A54-A11*A22*A33*A45*A54+ &
              A14*A23*A32*A41*A55-A13*A24*A32*A41*A55-A14*A22*A33*A41*A55+ &
              A12*A24*A33*A41*A55+A13*A22*A34*A41*A55-A12*A23*A34*A41*A55- &
              A14*A23*A31*A42*A55+A13*A24*A31*A42*A55+A14*A21*A33*A42*A55- &
              A11*A24*A33*A42*A55-A13*A21*A34*A42*A55+A11*A23*A34*A42*A55+ &
              A14*A22*A31*A43*A55-A12*A24*A31*A43*A55-A14*A21*A32*A43*A55+ &
              A11*A24*A32*A43*A55+A12*A21*A34*A43*A55-A11*A22*A34*A43*A55- &
              A13*A22*A31*A44*A55+A12*A23*A31*A44*A55+A13*A21*A32*A44*A55- &
              A11*A23*A32*A44*A55-A12*A21*A33*A44*A55+A11*A22*A33*A44*A55

        return

        end function det

        function shr3 ( jsr )
        implicit none
        integer :: jsr
        integer :: jsr_input
        integer :: shr3
        
        jsr_input = jsr
        
        jsr = ieor ( jsr, ishft ( jsr,   13 ) )
        jsr = ieor ( jsr, ishft ( jsr, - 17 ) )
        jsr = ieor ( jsr, ishft ( jsr,    5 ) )
        
        shr3 = jsr_input + jsr
        
        return
        end function shr3

        function r4_uni ( jsr )
        implicit none
        integer :: jsr
        integer :: jsr_input
        real(DP) :: r4_uni

        jsr_input = jsr

        jsr = ieor ( jsr, ishft ( jsr,   13 ) )
        jsr = ieor ( jsr, ishft ( jsr, - 17 ) )
        jsr = ieor ( jsr, ishft ( jsr,    5 ) )

        r4_uni = 0.5d0 + real(jsr_input+jsr)/65536.0d0/65536.0d0

        return
        end function r4_uni

        function trace(M)
        implicit none
        real(DP) :: M(5,5)
        real(DP) :: trace
        integer :: i

        trace = 0.0d0
        do i=1,5
                trace = trace + M(i,i)
        end do
        return
        end function

        function identity()
        implicit none
        real(DP) :: identity(5,5)
        integer :: i

        identity = 0.0d0
        do i = 1,5
                identity(i,i) = 1.0d0
        end do
        return
        end function identity
       
        function surprise_update(P1_mean,P2_mean,P1_invcov,P2_cov)
        implicit none
        real(DP) :: surprise_update
        real(DP) :: P1_mean(5),P2_mean(5),P1_invcov(5,5),P2_cov(5,5)

        surprise_update = 0.5d0*(dot_product((P1_mean-P2_mean),matmul(P1_invcov,(P1_mean-P2_mean)))-trace(identity()-matmul(P2_cov,P1_invcov)))
        return
        end function surprise_update
        
        function surprise_independent(P1_mean,P2_mean,P1_invcov,P2_cov)
        implicit none
        real(DP) :: surprise_independent
        real(DP) :: P1_mean(5),P2_mean(5),P1_invcov(5,5),P2_cov(5,5)

        surprise_independent = 0.5d0*(dot_product((P1_mean-P2_mean),matmul(P1_invcov,(P1_mean-P2_mean)))-trace(identity()+matmul(P2_cov,P1_invcov)))
        return
        end function surprise_independent

        function rand(jsr,down,up)
        implicit none
        integer :: i,jsr
        real(DP) :: x(5),rand(5),up(5),down(5)

        do i = 1,5
                x(i) = r4_uni(jsr)
        end do
        rand = x*(up-down)+down
        return
        end function rand

        function D_out(P1,P2)
        implicit none
        real(DP) :: D_out,P1,P2,log_coeff

        if (P1.eq.0.0d0) then
                D_out = 0.0d0
        else
                log_coeff = P2/P1
                if (log_coeff.eq.0.0d0) then
                        D_out = 0.0d0
                else
                        D_out = P2*dlog(log_coeff)
                end if
        end if
        return
        end function D_out

        subroutine get_P(filename,mean,cov,invcov,norm)
        implicit none
        character (len=50) :: filename
        real(DP) :: mean(5),cov(5,5),invcov(5,5),norm
 
        call read_mean_cov(filename,mean,cov)
        call inv(cov,invcov)
        norm = ((2.0d0*pi)**-2.5d0)*(det(cov)**-0.5d0)
        end subroutine get_P
        
end module tension_mod

program tension
        use tension_mod
        use omp_lib
        implicit none
        character (len=50) :: P1_file,P2_file,filename,type_calc
        integer :: i
        integer(kind=8) :: nmc
        integer :: threads,jsr,id
        integer, allocatable :: seed(:)
        logical :: ms
        real(DP) :: P1_mean(5), P1_cov(5,5),P1_invcov(5,5),P1_norm
        real(DP) :: P2_mean(5), P2_cov(5,5),P2_invcov(5,5),P2_norm
        real(DP) :: up(5),down(5),domainsize
        real(DP) :: S,D1D2,D1D2_sq,D1D2_res,D2D1,D2D1_sq,D2D1_res,P1_temp,P2_temp,error_D1D2,error_D2D1
        real(DP) :: means(5)
        real(DP), allocatable :: D1D2_temp(:),D1D2_sq_temp(:),D2D1_temp(:),D2D1_sq_temp(:)

        call getarg(1,P1_file)
        call getarg(2,P2_file)
        call getarg(3,filename)
        call getarg(4,type_calc)

        call get_P(P1_file,P1_mean,P1_cov,P1_invcov,P1_norm)
        call get_P(P2_file,P2_mean,P2_cov,P2_invcov,P2_norm)


        domainsize = 1
        do i=1,5
                up(i) = max(P1_mean(i)+5.0d0*dsqrt(P1_cov(i,i)),P2_mean(i)+5.0d0*dsqrt(P2_cov(i,i)))
                down(i) = min(P1_mean(i)-5.0d0*dsqrt(P1_cov(i,i)),P2_mean(i)-5.0d0*dsqrt(P2_cov(i,i)))
                domainsize = domainsize*(up(i)-down(i))
        end do

        open(unit=501,file="output/"//trim(filename)//".txt")

        write(501,*) 'Surprise D1 updates D2 = ',surprise_update(P1_mean,P2_mean,P1_invcov,P2_cov)
        write(501,*) 'Surprise D2 updates D1 = ',surprise_update(P2_mean,P1_mean,P2_invcov,P1_cov)
        write(501,*) 'Surprise D1 compared to D2 = ',surprise_independent(P1_mean,P2_mean,P1_invcov,P2_cov)
        write(501,*) 'Surprise D2 compared to D1 = ',surprise_independent(P2_mean,P1_mean,P2_invcov,P1_cov)

        !$omp parallel
                threads = omp_get_num_threads()
        !$omp end parallel
        call omp_set_num_threads(threads)
        jsr = 123456789
        allocate(seed(0:threads-1))
        allocate(D1D2_temp(0:threads-1),D1D2_sq_temp(0:threads-1))
        allocate(D2D1_temp(0:threads-1),D2D1_sq_temp(0:threads-1))

        do i = 0,threads-1
                seed(i) = shr3(jsr)
        end do

        write(501,*) "Number of threads = ",threads
                write(501,*) 
        nmc = 10
        ms = .true.
        

        do while (ms)        
                D1D2_temp(:) = 0.0d0
                D2D1_temp(:) = 0.0d0
                D1D2_sq_temp(:) = 0.0d0
                D2D1_sq_temp(:) = 0.0d0
                !$omp parallel shared(seed,D1D2_temp,D1D2_sq_temp,D2D1_temp,D2D1_sq_temp),private(means,jsr,id,i,D1D2_res,D2D1_res,P1_temp,P2_temp)
                !$omp do schedule(guided)
                do i = 1,nmc
                        id = omp_get_thread_num()
                        jsr = seed(id)
                        means = rand(jsr,down,up)
                        P1_temp = multivariate_normal(means,P1_mean,P1_norm,P1_invcov)
                        P2_temp = multivariate_normal(means,P2_mean,P2_norm,P2_invcov)
                        D1D2_res = D_out(P1_temp,P2_temp)
                        D2D1_res = D_out(P2_temp,P1_temp)
                        D1D2_temp(id) = D1D2_temp(id) + D1D2_res
                        D2D1_temp(id) = D2D1_temp(id) + D2D1_res
                        D1D2_sq_temp(id) = D1D2_sq_temp(id) + D1D2_res**2.0d0
                        D2D1_sq_temp(id) = D2D1_sq_temp(id) + D2D1_res**2.0d0
                        seed(id) = jsr
                end do
                !$omp end do
                !$omp end parallel
                D1D2 = domainsize*sum(D1D2_temp)/nmc
                D2D1 = domainsize*sum(D2D1_temp)/nmc
                D1D2_sq = sum(D1D2_sq_temp)
                D2D1_sq = sum(D2D1_sq_temp)
                error_D1D2 = 100.0d0*(dsqrt((D1D2_sq-((D1D2/domainsize)**2.0d0)/nmc)/(nmc-1.0d0)/nmc)*domainsize)/D1D2
                error_D2D1 = 100.0d0*(dsqrt((D2D1_sq-((D2D1/domainsize)**2.0d0)/nmc)/(nmc-1.0d0)/nmc)*domainsize)/D2D1
                if (D1D2.eq.0.0d0) then
                        write(501,*)
                        write(501,*) "Result D1 D2= 0 with nmc = ",nmc
                        nmc = nmc*10.0d0
                else if (D2D1.eq.0.0d0) then
                        write(501,*)
                        write(501,*) "Result D2 D1= 0 with nmc = ",nmc
                        nmc = nmc*10.0d0
                else
                        write(501,*) 
                        write(501,*) "With nmc = ",nmc
                        write(501,*) "Result D1 D2 = ",D1D2,"+-",error_D1D2
                        write(501,*) "Result D2 D1 = ",D2D1,"+-",error_D2D1
                        if ((error_D1D2.le.1.0d0).and.(error_D2D1.le.1.0d0)) then
                                ms = .false.
                        else
                               nmc = nmc*10.0d0
                        end if
                end if
                if (nmc.ge.10000000000000000) ms = .false.
        end do        
        close(501)

end program tension
