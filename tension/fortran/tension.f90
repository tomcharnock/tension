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

        function erfi(y)
        implicit none
        real(DP), parameter :: qa = 9.16461398268964d-01 
        real(DP), parameter :: qb = 2.31729200323405d-01
        real(DP), parameter :: qc = 4.88826640273108d-01
        real(DP), parameter :: qd = 1.24610454613712d-01
        real(DP), parameter :: q0 = 4.99999303439796d-01
        real(DP), parameter :: q1 = 1.16065025341614d-01
        real(DP), parameter :: q2 = 1.50689047360223d-01
        real(DP), parameter :: q3 = 2.69999308670029d-01
        real(DP), parameter :: q4 = -7.28846765585675d-02
        real(DP), parameter :: pa = 3.97886080735226000d+00
        real(DP), parameter :: pb = 1.20782237635245222d-01
        real(DP), parameter :: p0 = 2.44044510593190935d-01
        real(DP), parameter :: p1 = 4.34397492331430115d-01
        real(DP), parameter :: p2 = 6.86265948274097816d-01
        real(DP), parameter :: p3 = 9.56464974744799006d-01
        real(DP), parameter :: p4 = 1.16374581931560831d+00
        real(DP), parameter :: p5 = 1.21448730779995237d+00
        real(DP), parameter :: p6 = 1.05375024970847138d+00
        real(DP), parameter :: p7 = 7.13657635868730364d-01
        real(DP), parameter :: p8 = 3.16847638520135944d-01
        real(DP), parameter :: p9 = 1.47297938331485121d-02
        real(DP), parameter :: p10 = -1.05872177941595488d-01
        real(DP), parameter :: p11 = -7.43424357241784861d-02
        real(DP), parameter :: p12 = 2.20995927012179067d-03
        real(DP), parameter :: p13 = 3.46494207789099922d-02
        real(DP), parameter :: p14 = 1.42961988697898018d-02
        real(DP), parameter :: p15 = -1.18598117047771104d-02
        real(DP), parameter :: p16 = -1.12749169332504870d-02
        real(DP), parameter :: p17 = 3.39721910367775861d-03
        real(DP), parameter :: p18 = 6.85649426074558612d-03
        real(DP), parameter :: p19 = -7.71708358954120939d-04
        real(DP), parameter :: p20 = -3.51287146129100025d-03
        real(DP), parameter :: p21 = 1.05739299623423047d-04
        real(DP), parameter :: p22 = 1.12648096188977922d-03
        real(DP) erfi
        real(DP) s,t,u,w,x,y,z
        
        z = y
        if (y .gt. 1.0d0) z = 2.0d0 - y
        w = qa - dlog(z)
        u = dsqrt(w)
        s = (qc + dlog(u)) / w
        t = 1.0d0 / (u + qb)
        x = u * (1.0d0 - s * (0.5d0 + s * qd)) - ((((q4 * t + q3) * t + q2) * t + q1) * t + q0) * t
        t = pa / (pa + x)
        u = t - 0.5d0
        s = (((((((((p22 * u + p21) * u + p20) * u + p19) * u + p18) * u + p17) * u + p16) * u + p15) * u + p14) * u + p13) * u + p12 
        s = ((((((((((((s * u + p11) * u + p10) * u + p9) * u + p8) * u + p7) * u + p6) * u + p5) * u + p4) * u + p3) * u + p2) * u + p1) * u + p0) * t - z * dexp(x * x - pb)
        x = x + s * (1.0d0 + x * s)
        if (y .gt. 1.0d0) x = -x
        erfi = x
        return
        end function erfi

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

end module tension_mod

program tension
        use tension_mod
        use omp_lib
        implicit none
        character (len=50) :: filename
        real(DP) :: CMB_mean(5), CMB_cov(5,5),CMB_invcov(5,5),CMB_norm
        real(DP) :: LSS_mean(5), LSS_cov(5,5),LSS_invcov(5,5),LSS_norm
        real(DP) :: up(5),down(5),domainsize
        real(DP) :: results, error,percenterror
        real(DP) :: rand(5),x(5),normal,normal_CMB,normal_LSS,normal_sq,normal_sq_CMB,normal_sq_LSS
        real(DP) :: normal_res,normal_res_CMB,normal_res_LSS,normalise
        real(DP), allocatable :: normal_temp(:),normal_temp_CMB(:),normal_temp_LSS(:)
        real(DP), allocatable :: normal_sq_temp(:),normal_sq_temp_CMB(:),normal_sq_temp_LSS(:)
        integer(kind=8) :: i,nmc,j
        logical :: needsmorepoints
        integer :: threads,id
        integer :: jsr
        integer, allocatable :: seed(:)
        
        call getarg(1,filename)!filename = 'CMB/CMB'
        call read_mean_cov(filename,CMB_mean,CMB_cov)
        call inv(CMB_cov,CMB_invcov)
        CMB_norm = ((2.0d0*pi)**-2.5d0)*(det(CMB_cov)**-0.5d0)

        call getarg(2,filename)!filename = 'Strong/Strong_L'
        call read_mean_cov(filename,LSS_mean,LSS_cov)
        call inv(LSS_cov,LSS_invcov)
        LSS_norm = ((2.0d0*pi)**-2.5d0)*(det(LSS_cov)**-0.5d0)

        domainsize = 1
        do i=1,5
                up(i) = max(CMB_mean(i)+5.0d0*dsqrt(CMB_cov(i,i)),LSS_mean(i)+5.0d0*dsqrt(LSS_cov(i,i)))
                down(i) = min(CMB_mean(i)-5.0d0*dsqrt(CMB_cov(i,i)),LSS_mean(i)-5.0d0*dsqrt(LSS_cov(i,i)))
                domainsize = domainsize*(up(i)-down(i))
        end do
        !$omp parallel
        threads = omp_get_num_threads()
        !$omp end parallel
        call omp_set_num_threads(threads)
        jsr = 123456789
        allocate(seed(0:threads-1))
        allocate(normal_temp(0:threads-1),normal_sq_temp(0:threads-1))
        allocate(normal_temp_CMB(0:threads-1),normal_sq_temp_CMB(0:threads-1))
        allocate(normal_temp_LSS(0:threads-1),normal_sq_temp_LSS(0:threads-1))

        do i = 0,threads-1
                seed(i) = shr3(jsr)
        end do

        call getarg(3,filename)
        open(unit=501,file="output/"//trim(filename)//".txt")

        write(501,*) "Number of threads = ",threads
        nmc = 10
        needsmorepoints = .true.

        do while (needsmorepoints)        
                normal_temp(:) = 0.0d0
                normal_temp_CMB(:) = 0.0d0
                normal_temp_LSS(:) = 0.0d0
                normal_sq_temp(:) = 0.0d0
                normal_sq_temp_CMB(:) = 0.0d0
                normal_sq_temp_LSS(:) = 0.0d0
                !$omp parallel shared(seed,normal_temp,normal_temp_CMB,normal_temp_LSS,normal_sq_temp,normal_sq_temp_CMB,normal_sq_temp_LSS),private(x,rand,jsr,id,i,j,normal_res,normal_res_CMB,normal_res_LSS)
                !$omp do schedule(guided)
                do i = 1,nmc
                        id = omp_get_thread_num()
                        jsr = seed(id)
                        do j = 1,5
                                rand(j) = r4_uni(jsr)
                        end do
                        x = rand*(up-down)+down
                        normal_res_CMB = multivariate_normal(x,CMB_mean,CMB_norm,CMB_invcov)
                        normal_res_LSS = multivariate_normal(x,LSS_mean,LSS_norm,LSS_invcov)
                        normal_res = normal_res_CMB*normal_res_LSS
                        seed(id) = jsr
                        normal_sq_temp(id) = normal_sq_temp(id) + normal_res**2.0d0
                        normal_sq_temp_CMB(id) = normal_sq_temp_CMB(id) + normal_res_CMB**2.0d0
                        normal_sq_temp_LSS(id) = normal_sq_temp_LSS(id) + normal_res_LSS**2.0d0
                        normal_temp(id) = normal_temp(id) + normal_res
                        normal_temp_CMB(id) = normal_temp_CMB(id) + normal_res_CMB
                        normal_temp_LSS(id) = normal_temp_LSS(id) + normal_res_LSS
                end do
                !$omp end do
                !$omp end parallel
                normal = sum(normal_temp)
                normal_CMB = sum(normal_temp_CMB)
                normal_LSS = sum(normal_temp_LSS)
                normal_sq = sum(normal_sq_temp)
                normal_sq_CMB = sum(normal_sq_temp_CMB)
                normal_sq_LSS = sum(normal_sq_temp_LSS)
                normalise = (domainsize/nmc)**2.0d0*normal_LSS*normal_CMB
                results = domainsize*normal/nmc!*normalise
                error = domainsize*dsqrt(((normal_sq-((normal/nmc)**2.0d0)/nmc)/(nmc-1.0d0))/nmc)
                if (results.eq.0.0d0) then
                        write(501,*) 
                        write(501,*) "Result = 0 with nmc = ",nmc
                        nmc = nmc*10.0d0
                else
                        percenterror = 100.0d0*error/results
                
                        write(501,*) 
                        write(501,*) "With nmc = ",nmc
                        write(501,*) "Result = ",dsqrt(2.0d0)*erfi(results),"+-",percenterror
                        write(501,*) "Result = ",results
                        write(501,*) "Normalised result = ",dsqrt(2.0d0)*erfi(results*normalise)
                        write(501,*) "Normalisation = ",normalise
                        if (percenterror.le.1.0d-10) then
                                needsmorepoints = .false.
                        else
                               nmc = nmc*10.0d0
                        end if
                end if
                if (nmc.ge.10000000000000000) stop
        end do        
        close(501)

end program tension
