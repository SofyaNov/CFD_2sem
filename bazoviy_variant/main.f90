program navie
    implicit none
    INTEGER, PARAMETER :: IO = 12, IO_RES=13
    INTEGER :: I, J, NI, NJ, NITER
    REAL(8) :: L, H, U0, MU, NU, R0, P0, L1
    REAL(8) :: dx, dy, CFL, Uref, EPS, delta_t, A
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: X_Node, Y_Node
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: X_Cell, Y_Cell
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: P_c, U_c, V_c

    call readinput(IO,'Input.txt',L,H,NI,NJ,NITER,CFL,Uref,EPS,U0,MU,R0,P0,NU,L1)

    ALLOCATE(X_Node(NI, NJ), Y_Node(NI, NJ))
    ALLOCATE(X_Cell(0:NI, 0:NJ), Y_Cell(0:NI, 0:NJ))

    ALLOCATE(U_c(0:NI, 0:NJ), V_c(0:NI, 0:NJ), P_c(0:NI, 0:NJ))

    call setka(NI,NJ,L,H,dx,dy,X_Node,Y_Node,X_Cell,Y_Cell)

    WRITE(*,*) 'L = ', L, '  NI = ', NI, '  dx = ', dx
    WRITE(*,*) 'H = ', H, '  NJ = ', NJ, '  dy = ', dy
    WRITE(*,*) 'Re = ', U0 * L / NU
    
    call initial_field(NI,NJ,dx,dy,CFL,U0,P0,U_c,V_c,P_c,A,delta_t)
    call NavierStokes(NI, NJ, Niter, eps, dx, dy, NU, delta_t, a, U0, P0, &
                         U_c, V_c, P_c, IO_RES,L1)
    call OutputFields_Cell(IO, NI, NJ, X_Node, Y_Node, U_c, V_c, P_c)
    
    DEALLOCATE(X_Node, Y_Node, X_Cell, Y_Cell)
    DEALLOCATE(U_c, V_c, P_c)

end 



SUBROUTINE NavierStokes(NI, NJ, N_iter, eps, dx, dy, NU, delta_t, a, U0, P0, &
                         U, V, P, IO_RES,L1)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NI, NJ, N_iter, IO_RES
    REAL(8), INTENT(IN) :: eps, dx, dy, NU, delta_t, a, U0, P0,L1
    REAL(8), DIMENSION(0:NI, 0:NJ), INTENT(INOUT) :: U, V, P
    REAL(8), DIMENSION(0:NI, 0:NJ) :: U_new, V_new, P_new

    REAL(8), DIMENSION(0:NI,0:NJ) :: conv_u, conv_v, conv_p
    REAL(8), DIMENSION(0:NI,0:NJ) :: f_u_x, f_v_x, f_u_y, f_v_y, f_p_x, f_p_y
    real(8) :: u_bound, v_bound

    INTEGER :: i, j, t
    REAL(8) :: u_half, v_half, p_half, dudx, dvdx, dudy, dvdy
    REAL(8) :: ERR_U, ERR_V, ERR_P

    OPEN(IO_RES, FILE='residuals.plt',STATUS='REPLACE')

    WRITE(IO_RES,*) 't ERR_U ERR_V ERR_P'

    U_new = U
    V_new = V
    P_new = P
    f_u_x = 0.0D0; f_v_x = 0.0D0
    f_u_y = 0.0D0; f_v_y = 0.0D0

    DO t = 1, N_iter

        call bc(NI,NJ,L1,dx,P0,U0,U,V,P)
        
        do j=1,NJ-1
            do i=0,NI-1
                u_bound = 0.5D0*(U(i+1,j) + U(i,j))
                if (u_bound>=0) then 
                    u_half= U(i,j)
                    v_half= V(i,j)
                    p_half = P(i+1,j)
                else 
                    u_half= U(i+1,j)
                    v_half= V(i+1,j)
                    p_half = P(i,j)
                end if
                !p_half = 0.5D0*(P(i,j) + P(i+1,j))
                dudx = (U(i+1,j) - U(i,j)) / dx
                dvdx = (V(i+1,j) - V(i,j)) / dx

                f_u_x(i,j) = (u_half**2) / dx + p_half / dx - (NU * dudx) / dx
                f_v_x(i,j) = (u_half * v_half) / dx - (NU * dvdx) / dx
                f_p_x(i,j)=u_half/dx
            end do
        end do

        
        do j=0,NJ-1
            do i=1,NI-1
                v_bound = 0.5D0*(V(i,j+1) + V(i,j))
                if (v_bound>=0) then 
                    u_half= U(i,j)
                    v_half= V(i,j)
                    p_half = P(i,j+1)
                else 
                    u_half= U(i,j+1)
                    v_half= V(i,j+1)
                    p_half = P(i,j)
                end if
                !коррекция
                if (j == 0) then
                    u_half = (U(i,j+1)+U(i,j))/2.0
                    v_half = v_bound
                end if
                !p_half = 0.5D0*(P(i,j) + P(i,j+1)) 
                dudy = (U(i,j+1) - U(i,j)) / dy
                dvdy = (V(i,j+1) - V(i,j)) / dy

                f_u_y(i,j) = (u_half * v_half) / dy - (NU * dudy) / dy
                f_v_y(i,j) = (v_half**2) / dy + p_half / dy - (NU * dvdy) / dy
                f_p_y(i,j)=v_half/dy
            end do
        end do

        DO j = 1, NJ
            DO i = 1, NI

                conv_u(i,j) = f_u_x(i,j) - f_u_x(i-1,j) + f_u_y(i,j) - f_u_y(i,j-1)
                conv_v(i,j) = f_v_x(i,j) - f_v_x(i-1,j) + f_v_y(i,j) - f_v_y(i,j-1)
                conv_p(i,j) = f_p_x(i,j) - f_p_x(i-1,j) + f_p_y(i,j) - f_p_y(i,j-1)

                U_new(i,j) = U(i,j) - delta_t * conv_u(i,j)
                V_new(i,j) = V(i,j) - delta_t * conv_v(i,j)
                P_new(i,j) = P(i,j) -delta_t/A* conv_p(i,j)

            END DO
        END DO

        ERR_U = MAXVAL(ABS(U_new(1:NI-1, 1:NJ-1) - U(1:NI-1, 1:NJ-1)))
        ERR_V = MAXVAL(ABS(V_new(1:NI-1, 1:NJ-1) - V(1:NI-1, 1:NJ-1)))
        ERR_P = MAXVAL(ABS(P_new(1:NI-1, 1:NJ-1) - P(1:NI-1, 1:NJ-1)))

        WRITE(IO_RES,'(I8,3E16.8)') t, ERR_U, ERR_V, ERR_P
        FLUSH(IO_RES)

        U = U_new
        V = V_new
        P = P_new

        IF (MOD(t, 500) == 0) THEN
            WRITE(*,'(A,I6,A,3E12.4)') 'iter: ', t, ' | RES U, V, P = ', ERR_U, ERR_V, ERR_P
        END IF

        IF (ERR_U < eps .AND. ERR_V < eps .AND. ERR_P < eps) THEN
            WRITE(*,*) N_iter, t
            EXIT
        END IF

    END DO

    call bc(NI,NJ,L1,dx,P0,U0,U,V,P)

    CLOSE(IO_RES)

END SUBROUTINE NavierStokes


    
                





















SUBROUTINE OutputFields_Cell(IO, NI, NJ, X, Y, U, V, P)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NI, NJ, IO
    REAL(8), DIMENSION(NI, NJ), INTENT(IN) :: X, Y
    REAL(8), DIMENSION(0:NI, 0:NJ), INTENT(IN) :: U, V, P

    OPEN(IO, FILE='data.plt')

    WRITE(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"'
    WRITE(IO,*) 'ZONE I=', NI, ', J=', NJ, ', DATAPACKING=BLOCK, VARLOCATION=([3-5]=CELLCENTERED)'

    WRITE(IO,'(100E25.16)') X(1:NI, 1:NJ)
    WRITE(IO,'(100E25.16)') Y(1:NI, 1:NJ)
    WRITE(IO,'(100E25.16)') U(1:NI-1, 1:NJ-1)
    WRITE(IO,'(100E25.16)') V(1:NI-1, 1:NJ-1)
    WRITE(IO,'(100E25.16)') P(1:NI-1, 1:NJ-1)
    
    close(IO)

END SUBROUTINE OutputFields_Cell




subroutine bc(NI,NJ,L1,dx,P0,U0,U,V,P)
    implicit none
    integer, intent(in) :: NI,NJ
    real(8),intent(in) :: P0,U0,L1,dx
    real(8),dimension(0:NI,0:NJ), intent(inout) :: U,V,P
    integer :: i,ni1
    
    ni1=int(L1/dx+1)
    
    
    U(0,:)=2.0D0 * U0 - U(1,:)
    V(0, :) = -V(1, :)
    P(0, :) = P(1, :)
    U(NI, :) = U(NI-1, :)
    V(NI, :) = V(NI-1, :)
    P(NI, :) = P0

    do i=0,ni1
        U(i, 0) = U(i, 1)
        V(i, 0) = -V(i, 1)
        P(i, 0) = P0
    end do
    do i=ni1+1,NI
        U(i, 0) = -U(i, 1)
        V(i, 0) = -V(i, 1)
        P(i, 0) = P(i, 1)
    end do

    do i=1,NI
        if (V(i,NJ)>0) then
            U(i, NJ) = U(i,NJ-1)!2.0*U0-U(:, NJ-1)
            V(i, NJ) = V(i, NJ-1)
            P(i, NJ) = -P(i, NJ-1)
        else 
            U(i,NJ)=U0
            P(i, NJ) = P(i, NJ-1)
            V(i, NJ) = V(i, NJ-1)
        end if
    end do 
    
    
end subroutine bc




subroutine readinput(IO,filename,L,H,NI,NJ,NITER,CFL,Uref,EPS,U0,MU,R0,P0,NU,L1)
    implicit none
    INTEGER, intent(in) :: IO
    CHARACTER(*), intent(in) :: filename
    INTEGER, intent(out) :: NI,NJ,NITER
    REAL(8), intent(out) :: L,H, CFL, Uref,EPS,U0,MU,R0,P0,NU,L1

    OPEN(IO, FILE=filename)
        READ(IO, *) L
        READ(IO, *) H
        READ(IO, *) NI
        READ(IO, *) NJ
        READ(IO, *) NITER
        READ(IO, *) CFL
        READ(IO, *) Uref
        READ(IO, *) EPS
        READ(IO, *) U0
        READ(IO, *) MU
        READ(IO, *) R0
        READ(IO, *) P0
        READ(IO, *) L1
    CLOSE(IO)
    NU = MU / R0
END subroutine readinput




subroutine setka(NI,NJ,L,H,dx,dy,X_Node,Y_Node,X_Cell,Y_Cell)
    implicit none
    integer, intent(in) :: NI,NJ
    real(8), intent(in) :: L,H 
    real(8), intent(out) :: dx,dy
    real(8),dimension(NI,NJ), intent(out) :: X_Node,Y_Node
    real(8),dimension(0:NI,0:NJ), intent(out) :: X_Cell,Y_Cell 
    integer :: i,j

    dx = L / (NI - 1)
    dy = H / (NJ - 1)

    DO J = 1, NJ
        DO I = 1, NI
            X_Node(I,J) = (I - 1) * dx
            Y_Node(I,J) = (J - 1) * dy
        END DO
    END DO

    X_Cell(0, 1:NJ) = -dx / 2.0D0
    Y_Cell(0, 1:NJ) = Y_Node(1, 1:NJ) + dy / 2.0D0
    X_Cell(1:NI, 0) = X_Node(1:NI, 1) + dx / 2.0D0
    Y_Cell(1:NI, 0) = -dy / 2.0D0

    DO J = 1, NJ
        DO I = 1, NI
            X_Cell(I,J) = X_Node(I,J) + dx / 2.0D0
            Y_Cell(I,J) = Y_Node(I,J) + dy / 2.0D0
        END DO
    END DO
    X_Cell(0,0) = -dx / 2.0d0
    Y_Cell(0,0) = -dy / 2.0d0
end subroutine setka 




subroutine initial_field(NI,NJ,dx,dy,CFL,U0,P0,U_c,V_c,P_c,A,delta_t)
    implicit none
    integer,intent(in) :: NI,NJ
    real(8),intent(in) :: U0,P0,dx,dy,CFL
    real(8),dimension(0:NI,0:NJ), intent(out) :: U_c,V_c,P_c 
    real(8), intent(out) :: A, delta_t
    U_c = U0
    V_c = 0.0D0
    P_c = P0
    A=U0**2
    !delta_t = CFL * MIN(dx, dy) / (U0 + 1.0D-5)
    delta_t = CFL * MIN(dx, dy) / (U0 + SQRT(U0**2 + A))
end subroutine initial_field









