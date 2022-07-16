subroutine td_bfgs(N, K, y, lamreg, nlam, sk, initp, rho, ftol, partol, abstol,&
     maxit, verbose, resfn, resgradf, resparam, resKa, conv)
    implicit none
    ! Input variables
    integer, intent(in) :: N ! sample size
    double precision, dimension (N, N), intent(in) :: K ! kernel matrix
    double precision, dimension (N), intent(in) :: y ! outcome variable
    integer, intent(in) :: nlam ! number of lambdas - length of lamreg
    double precision, dimension (nlam), intent(in) :: lamreg    ! equals 1 when called by td_sk
    logical, intent(in) :: sk  ! flag for whether td_bfgs is used for td_sk
    double precision, dimension (N), intent(in) :: initp   ! initial value of param
    double precision, intent(in) :: rho ! index parameter in tweedie distribution
    double precision, intent(in) :: ftol ! threshold for function value convergence
    double precision, intent(in) :: partol ! threshold for parameter value convergence
    double precision, intent(in) :: abstol ! threshold for absolute convergence - stops when fn < abstol
    ! double precision, intent(in) :: p2 ! penalty for weights, used in TD_SK for printing results. 0 otherwise
    integer, intent(in) :: maxit ! maximum iterations
    logical, intent(in) :: verbose ! flag for displaying messages

    ! Output variables
    double precision, dimension (nlam), intent(out) :: resfn ! result function value
    double precision, dimension (N, nlam), intent(out) :: resgradf ! result gradient value
    double precision, dimension (N, nlam), intent(out) :: resparam ! result parameter value
    double precision, dimension (N, nlam), intent(out) :: resKa   ! output Ka: K * param
    integer, dimension (nlam), intent(out) :: conv  ! convergence situation
    ! 0 - BFGS converges
    ! 1 - (Switch to) GD, converges
    ! 8 - BFGS goes over inner maxit
    ! 9 - GD goes over inner maxit
    ! 98 - should not happen: fn increased in BFGS
    ! 99 - should not happen: fn increased in GD
    ! 999 - Initial. should be changed to something else

    ! local variables
    integer :: qn_iter ! number of quasi-newton iterations
    integer :: gd_iter ! number of gradient descent iterations
    integer :: gd ! flag for whether to use gradient descent; 0 - qn, 1 - gd
    double precision :: fn_init ! initial function value for param = 0
    double precision, dimension (N) :: gradf_init  ! initial gradient wrt param
    double precision, dimension (N, N) :: B_init    ! initial inverse hessian diag(N)
    double precision :: fn_old ! old function value
    double precision :: fn ! function value
    double precision :: lambda ! current lambda
    double precision, dimension (N, N) :: B ! estimated inverse hessian in BFGS algorithm
    double precision, dimension (N) :: pvec ! p = -B * gradf_old
    double precision :: gradproj ! gradf_old * pvec
    double precision :: t ! step size
    double precision :: t_old ! previous step size
    double precision :: lb, ub ! lower bound and upper bound in line search for step size
    double precision :: infty ! infinity
    double precision :: wp1, wp2, wp3, wp4 ! quantities used in wolfe-powell condition
    double precision :: denom ! yvec * svec, denominator in the update of B
    double precision :: fnorm ! absolute difference between current and previous function values
    double precision :: parnorm ! sum of absolute difference between current and previous parameters
    double precision :: gradfsq ! L2 norm of gradf_old
    double precision :: am2 ! RHS of Armijo-Goldstein condition
    double precision, dimension (N) :: gradf_old, gradf ! old gradient, gradient
    double precision, dimension (N) :: param_old, param ! old parameter estimates, estimates
    double precision, dimension (N) :: svec ! param - param_old = t * pvec, the amount of update
    double precision, dimension (N) :: yvec ! gradf - gradf_old
    double precision, dimension (N) :: s_By ! svec - B * yvec
    double precision, dimension (N) :: Ka ! K * param_old
    double precision, dimension (N) :: Kexp1, Kexp2 ! exp[-(rho-1) K*alpha], exp[(2-rho) K*alpha]
    double precision, dimension (N, N) :: tm1, tm2 ! Bnew = B + tm1 - tm2
    integer :: i, j, l ! indices, l for lambda
    logical :: phaseA ! flag for phase A in line search
    logical :: accpt ! acceptable point in line search
    ! integer :: ls_interval  ! interval for print line search information during running
    ! integer :: lp_interval  ! interval for print loop information
    integer :: ls_iter  ! line search iteration
    integer, dimension (2) :: qn_gd_iter    ! 2d array for printing qn and gd iter
    logical :: iterating    ! flag for convergence
    integer :: convergence ! flag for convergence

    ! set up message display interval -  lp_interval, ls_interval
    ! if (verbose > 100) then
    !     ls_interval = mod(verbose, 100)
    !     lp_interval = (verbose-ls_interval)/100
    !     if (ls_interval == 0) ls_interval = 99999
    ! else
    !     ls_interval = 99999 ! to avoid warning message...
    !     lp_interval = verbose
    ! end if

    ! get fn_init, gradf_init
    ! if (.not. sk) then
    fn_init = 1.0d0/N * ( sum(y)/(rho - 1.0d0) + N/(2.0d0 - rho) )
    gradf_init = 1.0d0/N * ( reshape(matmul(K, reshape(-y,(/N,1/))), (/N/)) + sum(K, dim = 1) )
    ! end if

    ! initialize B matrix
    do i = 1, N
        B_init(i, i) = 1.0d0
        do j = (i+1), N
            B_init(i,j) = 0.0d0
            B_init(j,i) = 0.0d0
        end do
    end do

    ! lambda loop
    do l = 1, nlam
        lambda = lamreg(l)
        ! Quasi-Newton iterations
        qn_iter = 1

        ! Initialize B matrix
        B = B_init

        ! Initialize first objective function and gradient calculation
        ! warm start feature for non-sk BFGS
        if (sk) then
            param = initp
            ! calculate initial fn and gradf
            call fn_grad(N, K, y, param, lambda, rho, Ka, fn, gradf)
        else
            ! sort lambda in decreasing order (in R wrapper)
            ! starting from the second lambda, use the estimated values from previous lambda as init.
            if (l == 1) then
                param = 0.0d0
                fn = fn_init
                gradf = gradf_init
            end if
        end if

        ! ! Initialize norms so 1st iteration runs
        ! fnorm = 1.0d5
        ! parnorm = 1.0d5

        ! gradient descent set to 0 - BFGS
        gd = 0
        infty = HUGE(fnorm)

        iterating = .true.
        convergence = 999

        ! BFGS main loop
        do while (iterating)
            ! save previous values of param, fn and gradf in *_old
            param_old = param
            fn_old = fn
            gradf_old = gradf

            ! update
            pvec = matmul(-B, gradf_old)
            gradproj = dot_product(gradf_old, pvec)
            wp4 = 0.9d0 * gradproj
            t = 1.0d0
            ! line search setup
            phaseA = .true.
            lb = 0.0d0
            ! acceptable point
            accpt = .false.
            ls_iter = 0
            
            ! line search loop phase A
            do while ( (.not. accpt) .and. phaseA )
                svec = t * pvec
                ! update param, fn, gradf
                param = param_old + svec
                call fn_grad(N, K, y, param, lambda, rho, Ka, fn, gradf)
                ! write(*,*) N
                wp1 = fn
                wp2 = fn_old + 1.0d-4 * t * gradproj
                wp3 = dot_product(gradf, pvec)
                
                ! prevent infinity
                if (wp1>infty .or. wp2>infty .or. wp3>infty) then
                    t = t * 0.2d0
                else
                    if (wp1<=wp2) then
                        if (wp3>=wp4) then
                            accpt = .true.
                        else
                            t = 2.0d0 * t
                        end if ! wp3>wp4
                    else
                        phaseA = .false.
                    end if ! wp1<=wp2
                end if ! (wp1>infty .or. wp2>infty .or. wp3>infty) ! prevent infty
                ls_iter = ls_iter + 1
                ! if (verbose > 100) then
                !     if ((mod(ls_iter, ls_interval) == 0) .and. (mod(qn_iter, lp_interval) == 0)) then
                !         call intpr("---Line search iteration---", -1, ls_iter, 1)
                !         call dblepr("Phase A, step size:", -1, t, 1)
                !         call dblepr("Objective function:", -1, fn+p2, 1)
                !     end if
                ! end if
            end do ! ((.not. accept) .and. phaseA)
            ! end of line search loop phase A
            
            ! line search loop phase B
            if (.not. phaseA) then
                ub = t
                t_old = 1.0d5
                do while (.not. accpt)
                    t_old = t
                    t = (lb + ub) * 0.5d0
                    
                    ! exit loop phase B and move on to gradient descent if line search fails
                    if (t == t_old) then
                        gd = 1
                        exit
                    end if ! t == t_old

                    svec = t * pvec
                    ! update param, fn, gradf
                    param = param_old + svec
                    call fn_grad(N, K, y, param, lambda, rho, Ka, fn, gradf)

                    wp1 = fn
                    wp2 = fn_old + 1.0d-4 * t * gradproj
                    wp3 = dot_product(gradf, pvec)
                    
                    if (wp1<=wp2) then
                        if (wp3>=wp4) then
                            accpt = .true.
                        else
                            lb = t
                        end if ! (wp3>=wp4)
                    else
                        ub = t
                    end if ! (wp1<=wp2)
                    ls_iter = ls_iter + 1
                    ! if (verbose > 100) then
                    !     if ((mod(ls_iter, ls_interval) == 0) .and. (mod(qn_iter, lp_interval) == 0)) then
                    !         call intpr("---Line search iteration---", -1, ls_iter, 1)
                    !         call dblepr("Phase B, step size:", -1, t, 1)
                    !         call dblepr("Objective function:", -1, fn+p2, 1)
                    !     end if
                    ! end if
                end do ! (.not. accpt)
                
                ! exit BFGS and move on to gradient descent if line search fails
                if (gd == 1) then
                    if (verbose) then
                        call intpr("-----Parameter update iteration-----", -1, qn_iter, 1)
                        call intpr("Switch to gradient descent. BFGS aborted at iteration:",-1,qn_iter,1)
                    end if
                    exit
                end if

            end if ! (.not. phaseA)
            ! end of line search loop phase B

            ! inverse BFGS udpate
            yvec = gradf - gradf_old
            denom = dot_product(yvec, svec)
            s_By = svec - matmul(B, yvec)
            do i = 1, N - 1
                tm1(i,i) = 2.0d0 * s_By(i) * svec(i)
                do j = i + 1, N
                    tm1(i, j) = s_By(i) * svec(j) + svec(i) * s_By(j)
                    tm1(j, i) = tm1(i, j)
                end do
            end do
            tm1(N, N) = 2.0d0 * s_By(N) * svec(N)
            tm1 = tm1/denom

            do i = 1, N - 1
                tm2(i,i) = svec(i) * svec(i)
                do j = i + 1, N
                    tm2(i, j) = svec(i) * svec(j)
                    tm2(j, i) = tm2(i, j)
                end do
            end do
            tm2(N, N) = svec(N) * svec(N)
            tm2 = dot_product(s_By, yvec) * tm2/(denom**2)
            B = B + tm1 - tm2


            fnorm = fn_old - fn
            parnorm = sum(abs(svec))/N
            ! write(*,*) "par", parnorm<=partol, "fn", fn<=abstol, "fnorm", fnorm, "gd", gd

            ! Precision loss can cause fn to increase in BFGS updates
            ! Should be the result of line search failure (t = t old)
            ! BFGS will abort and transition to GD at the same iteration
            ! This test stays just to make sure.
            ! If correct, line search failure will change gd to 1 and exit BFGS loop before reaching here.
            if (fnorm < 0) then
                convergence = 98
                call intpr( "Fn increased at BFGS iteration:", -1, qn_iter, 1)
                call intpr("gd", -1, gd, 1)
                exit
            end if

            ! check convergence
            if (qn_iter>maxit) then
                call intpr("BFGS fails to converge before maximum iteration", -1, maxit, 1)
                convergence = 8
                exit
            end if

            if (sk) then
                ! stop when parnorm and fn both converge
                iterating = .not.(parnorm<=partol .and. fn<=abstol)
            else
                ! stop when any one of fnorm, parnorm, fn converges
                iterating = (fnorm>ftol) .and. (parnorm>partol) .and. (fn>abstol)
            end if

            ! if (verbose>0 .and. gd==0) then
            !     if (mod(qn_iter,  lp_interval) == 0) then
            !         call intpr("-----Parameter update iteration-----", -1, qn_iter, 1)
            !         call intpr("Total BFGS line search iterations =", -1, ls_iter, 1)
            !         if (verbose>=100) call dblepr("BFGS step size =", -1, t, 1)
            !         call dblepr("Objective fn =", -1, fn+p2, 1)
            !         call dblepr("Decrease in fn =", -1, fnorm, 1)
            !         call dblepr("Mean of absolute parameter update =", -1, parnorm, 1)
            !     end if
            ! end if

            if (.not. iterating) then
                convergence = 0
                if (verbose) then
                    call intpr("BFGS converged at iteration:", -1, qn_iter, 1)
                end if
            end if

            qn_iter = qn_iter + 1

        end do ! while ((fnorm>ftol .or. parnorm>partol) .and. qn_iter<maxit)
        ! end of BFGS main loop
        
        ! Transition to gradient descent when line search in BFGS fails
        ! No feasible fn and gradf was generated from BFGS loop that failed
        ! so GD loop picks up from param_old, fn_old and gradf_old

        gd_iter = 1
        if (gd == 1) then
            ! Gradient descent interations
            iterating = .true.
            fn = fn_old
            
            ! gradient descent main loop
            do while (iterating)
                t = 1.0d0
                accpt = .false.
                gradfsq = dot_product(gradf_old, gradf_old)
                ls_iter = 0
                ! backtracking line search loop
                do while (.not. accpt)
                    svec = - t * gradf_old
                    param = param_old + svec
                    ! param = param_old - t * gradf_old
                    ! no gradf calculation here
                    call fn_grad_prep(N, K, param, rho, Ka, Kexp1, Kexp2)
                    call compute_fn(N, Ka, Kexp1, Kexp2, y, param, rho, lambda, fn)
                    am2 = fn_old - 0.5d0 * t * gradfsq
                    
                    ! Armijo-Goldstein
                    if (fn<=am2) then
                        accpt = .true.
                    else
                        t = t * 0.9d0
                    end if ! (fn<am2)

                    ls_iter = ls_iter + 1

                    ! if (verbose > 100) then
                    !     if ((mod(ls_iter, ls_interval)==0).and.(mod((qn_iter+gd_iter), lp_interval)==0)) then
                    !         call intpr("---Line search iteration---", -1, ls_iter, 1)
                    !         call dblepr("Gradient descent, step size:", -1, t, 1)
                    !         call dblepr("Objective function:", -1, fn+p2, 1)
                    !     end if
                    ! end if
                end do ! (.not. accpt)
                ! end of line search loop

                call compute_grad(N, K, Ka, Kexp1, Kexp2, y, lambda, gradf)
                fnorm = fn_old - fn

                ! should not happen
                if (fnorm<0) then
                    call intpr("Fn increased at GD iteration:", -1, qn_iter+gd_iter, 1)
                    convergence = 99
                    exit
                end if

                parnorm = sum(abs(svec))/N
                ! move to the end later
                param_old = param
                fn_old = fn
                gradf_old = gradf

                ! check maxit
                if (gd_iter+qn_iter>maxit) then
                    call intpr("Failure to converge before maximum iteration", -1, maxit, 1)
                    convergence = 9
                    exit
                end if

                ! check convergence
                if (sk) then
                    ! stop when parnorm and fn both converge
                    iterating = .not.(parnorm<=partol .and. fn<=abstol)
                else
                    ! stop when any one of fnorm, parnorm, fn converges
                    iterating = (fnorm>ftol) .and. (parnorm>partol) .and. (fn>abstol)
                end if

                ! if (verbose > 0) then
                !     if (mod((qn_iter+gd_iter), lp_interval) == 0) then
                !         call intpr("-----Parameter update iteration-----", -1, qn_iter+gd_iter, 1)
                !         call intpr("Total GD line search iterations =", -1, ls_iter, 1)
                !         if (verbose>=100) call dblepr("GD step size =", -1, t, 1)
                !         call dblepr("Objective fn =", -1, fn+p2, 1)
                !         call dblepr("Decrease in fn =", -1, fnorm, 1)
                !         call dblepr("Mean of absolute parameter update =", -1, parnorm, 1)
                !     end if
                ! end if

                if (.not. iterating) then
                    if (verbose) then
                        call intpr("Gradient descent converged at iteration:", -1, qn_iter+gd_iter, 1)
                    end if
                    convergence = 1
                end if

                gd_iter = gd_iter + 1

            end do ! ((fnorm>ftol .or. parnorm>partol) .and. gd_iter<(maxit-qn_iter))
            ! end of gradient descent main loop
        end if !(gd)

        if (verbose) then
            if (convergence == 0) qn_iter = qn_iter - 1
            qn_gd_iter(1) = qn_iter
            qn_gd_iter(2) = gd_iter-1
            ! call intpr("BFGS GD - Parameter update iterations:", -1, qn_gd_iter, 2)
        end if

        resparam(:, l) = param
        resfn(l) = fn
        resgradf(:, l) = gradf
        resKa(:, l) = Ka
        conv(l) = convergence
    end do ! end of lambda loop
end subroutine