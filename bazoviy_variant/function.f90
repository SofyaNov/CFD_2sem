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



