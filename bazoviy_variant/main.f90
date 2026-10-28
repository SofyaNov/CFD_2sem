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
