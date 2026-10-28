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
