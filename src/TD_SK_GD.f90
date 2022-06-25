! update kernel weights using gradient descent

subroutine sk_update(N, p, z, rho, sigma, lam1, lam2, dist, param, verbose,&
    partol, abstol, maxit, K, Ka, wt, wtac, wtac_size, fn, conv)
    implicit none
    ! input
    integer, intent(in) :: N   ! sample size
    integer, intent(in) :: p   ! dimension
    double precision, dimension (N), intent(in) :: z   ! outcome variable
    double precision, intent(in) :: rho    ! TD index parameter
    double precision, intent(in) :: sigma  ! RBF kernel parameter exp[-rho||x_i-x_j||^2]
    double precision, intent(in) :: lam1, lam2   ! penalty, 1-coefficient, 2-weight
    double precision, dimension (N, N, p), intent(in) :: dist  ! distance array
    double precision, dimension (N), intent(in) :: param   ! coefficients
    integer, intent(in) :: verbose    ! flag for displaying messages
    ! double precision, intent(in) :: ftol ! threshold for function value convergence
    double precision, intent(in) :: partol ! threshold for parameter value convergence
    double precision, intent(in) :: abstol ! threshold for absolute convergence - stops when fn < abstol
    integer, intent(in) :: maxit   ! maximum iterations
    
    ! inout
    double precision, dimension (N, N), intent(inout) :: K   ! kernel matrix
    double precision, dimension (N), intent(inout) :: Ka   ! K * param
    double precision, dimension (p), intent(inout) :: wt   ! sparse kernel weight
    logical, dimension (p), intent(inout) :: wtac  ! active set of weight
    integer, intent(inout) :: wtac_size    ! size of active set of weight
    double precision, intent(inout) :: fn  ! objective function value

    ! output
    integer, intent(out) :: conv  ! convergence situation
    ! 2 - gradient descent for weights converges
    ! 7 - GD for weights goes over inner maxit
    ! 50 - all weights are 0
    ! 97 - should not happen: fn increased in GD for weight update

    ! local
    integer :: i, j  ! indices
    double precision, dimension (p) :: gradw, gradw_old ! gradient wrt w, old gradient
    logical :: accpt    ! acceptable point in line search
    double precision :: gradwsq ! l2 norm of gradw
    double precision, dimension (N) :: Kexp1, Kexp2 ! exp[-(rho-1) K*alpha], exp[(2-rho) K*alpha]
    double precision :: am2 ! RHS of Armijo-Goldstein condition
    double precision :: t   ! step size
    integer :: ls_iter  ! line search iteration
    double precision :: fn_old  ! old function value
    double precision :: lambda  ! locally replace lam1 - penalty for coefficients
    integer, dimension (:), allocatable :: wtac_id ! numeric vector for indexing active weights
    double precision, dimension (p) :: wt_old   ! old weights
    logical, dimension (p) :: wtac_old ! old weight flag
    double precision :: fnorm ! absolute difference between current and previous function values
    double precision :: parnorm ! sum of absolute difference between current and previous parameters
    integer :: iter  ! number of gradient descent iterations
    integer :: print_interval   ! interval for print information during running
    logical :: iterating    ! flag for convergence
    

    ! content
    ! set up message display interval - print_interval
    if (verbose > 100) then
        print_interval = mod(verbose, 100)
        if (print_interval == 0) print_interval = 99999
    else
        print_interval = 99999 ! to avoid warning message...
    end if

    iterating  = .true.

    iter = 1

    ! Initial allocation of wtac_id
    allocate(wtac_id(wtac_size))

    ! gd main loop
    do while (iterating) ! .and. wtac_size>0, exit in loop

        call grad_wt(N, p, Ka, z, param, dist, rho, wt, wtac, wtac_size, sigma, lam1, lam2, gradw)

        fn_old = fn
        wt_old = wt
        wtac_old = wtac
        gradw_old = gradw

        ! line search init
        t = 1.0d0
        ls_iter = 1
        accpt = .false.
        
        ! convert logical wtac to numeric wtac_id
        i = 0
        do j = 1, p
            if (wtac(j)) then
                i = i + 1
                wtac_id(i) = j
            end if
        end do
        gradwsq = sum( gradw(wtac_id)**2 )  ! dot_product(gradw, gradw)

        ! line search loop
        do while (.not. accpt)
            ! update wt
            do j = 1, p
                if ( wtac_old(j) ) then
                    wt(j) = wt_old(j) - t * gradw(j)
                    if ( wt(j)<0.0d0 ) then ! negative weights
                        wt(j) = 0.0d0   ! set to zero
                        wtac(j) = .false.   ! remove from active set
                        ! wtac_size = wtac_size - 1   ! lose one more in the active set - size minus 1
                    else
                        wtac(j) = .true. 
                        if ( wt(j)>1 ) then    ! weights greater than 1
                            wt(j) = 1.0d0   ! set to 1
                        end if
                    end if
                ! else
                !     wt(j) = 0.0d0
                !     wtac(j) = .false.
                end if
            end do

            ! get size of active set
            wtac_size = count(wtac)

            ! get id of active set
            if (allocated(wtac_id)) deallocate(wtac_id)
            allocate(wtac_id(wtac_size))

            ! convert logical wtac to numeric wtac_id
            i = 0
            do j = 1, p
                if (wtac(j)) then
                    i = i + 1
                    wtac_id(i) = j
                end if
            end do

            ! update kernel matrix
            do i = 1, N
                do j = i, N
                    K(i,j) = exp( -sigma * sum( (dist(i,j,wtac_id) * wt(wtac_id))**2 ) )
                    K(j,i) = K(i,j)
                end do
            end do

            ! calculate function value for line search
            lambda = lam1    ! for calling compute_fn locally
            call fn_grad_prep(N, K, param, rho, Ka, Kexp1, Kexp2)   ! fn_grad(N, K, z, param, lambda, rho, fn, gradf)
            call compute_fn(N, Ka, Kexp1, Kexp2, z, param, rho, lambda, fn)
            ! fn_grad does not include penalty on the weights
            ! so added here
            fn = fn + lam2 * sum( wt(wtac_id) )
            am2 = fn_old - 0.5d0 * t * gradwsq
            accpt = (fn <= am2)

            ! Armijo-Goldstein
            if (.not. accpt) then
                t = t * 0.9d0
            end if ! (fn <= am2)
            ls_iter = ls_iter + 1

            if (verbose > 100) then
                if ( mod(ls_iter, print_interval) == 0 ) then
                    call intpr("---Line search iteration---", -1, ls_iter, 1)
                    call dblepr("Gradient descent step size:", -1, t, 1)
                    call dblepr("Objective function:", -1, fn, 1)
                end if
            end if
        end do ! (.not. accpt)
        ! end of line search loop

        ! if wtac_size == 0, exit the weight update loop
        if (wtac_size == 0) then
            ! invalid update
            conv = 50
            call intpr("WARNING: All weights are zero in weight update iteration:", -1, iter, 1)
            exit
        end if

        if (iter>maxit) then
            conv = 7
            call intpr("Failure to converge before maximum iteration:", -1, maxit, 1)
            exit
        end if

        fnorm = fn_old - fn

        if (fnorm < 0) then
            call intpr("Fn increased at weight update iteration:", -1, iter, 1)
            conv = 97
            exit
        end if

        parnorm = sum(abs(wt_old - wt))

        iterating = parnorm>partol .and. fn>abstol

        if (verbose > 0) then
            call intpr("-----Weight update iteration-----", -1, iter, 1)
            call intpr("Total line search iterations =", -1, ls_iter, 1)
            if (verbose>=100) call dblepr("Gradient descent step size =", -1, t, 1)
            call intpr("Active weight count =", -1, wtac_size, 1)
            call dblepr("Objective fn =", -1, fn, 1)
            call dblepr("Decrease in fn =", -1, fnorm, 1)
            call dblepr("Mean of absolute parameter update =", -1, parnorm, 1)
        end if

        if (.not. iterating) then
            if (verbose>0) then
                call intpr("Gradient descent for weights converged at iteration:", -1, iter, 1)
            end if
            conv = 2
        end if

        iter = iter + 1
        
        ! write(*,*) "wt parnorm", parnorm

    end do  ! end of gd main loop
    ! prepare inout variables for output
    ! fn, K, Ka, wt, wtac, wtac_size
    ! if (verbose > -1) call intpr("Weight update iterations:", -1, iter-1, 1)
    ! if (verbose == 0) then
    !     if (conv == 2) iter = iter - 1
    !     call intpr("GD - Weight update iterations:", -1, iter, 1)
    ! end if

end subroutine sk_update