! TD estimation with sparse kernel (RBF) feature

subroutine td_sk(N, p, x, y, lambda1, lambda2, sigmas, nval, rho,&
    fn_tol, param_tol, abs_tol, inner_ftol, inner_partol, max_it, inner_maxit,&
    verbose, fn_final, gradf_final, param_final, wt_final, convergence)
 
    ! input
    integer, intent(in) :: N   ! sample size
    integer, intent(in) :: p   ! dimension
    double precision, dimension (N, p), intent(in) :: x  ! covariate matrix
    double precision, dimension (N), intent(in) :: y   ! outcome variable
    double precision, dimension (nval), intent(in) :: lambda1, lambda2   ! penalty, 1-coefficient, 2-weight
    integer, intent(in) :: nval    ! number of penalty parameters - set to be the same for random 2d search
    double precision, dimension (nval), intent(in) :: sigmas  ! RBF kernel parameter exp[-rho||x_i-x_j||^2]
    double precision, intent(in) :: rho ! index parameter in tweedie distribution
    double precision, intent(in) :: fn_tol    ! threshold for function value convergence
    double precision, intent(in) :: param_tol  ! threshold for parameter value convergence
    double precision, intent(in) :: abs_tol  ! threshold for absolute convergence - stops when fn_old < abstol
    double precision, intent(in) :: inner_ftol, inner_partol   ! threshold for parameter/weight loop convergence
    integer, intent(in) :: max_it ! maximum iterations
    integer, intent(in) :: inner_maxit ! maximum iterations for parameter/weight loop
    logical, intent(in) :: verbose ! flag for displaying messages

    ! output
    double precision, dimension (nval), intent(out) :: fn_final   ! result function value
    double precision, dimension (N, nval), intent(out) :: gradf_final ! result gradient value
    double precision, dimension (N, nval), intent(out) :: param_final ! result parameter value
    double precision, dimension (p, nval), intent(out) :: wt_final  ! result weight
    integer, dimension (nval), intent(out) :: convergence  ! convergence

    ! local
    double precision, dimension (N, N, p) :: dist  ! distance array
    integer :: i, j ! indices
    integer :: lid   ! indices for lambda1 and lambda2
    double precision :: lam1, lam2, sigma   ! individual lambda1, lambda2 and sigma in each lam1-lam2-sigma loop
    double precision, dimension (N, N) :: K ! kernel matrix
    double precision, dimension (N) :: param    ! parameter estimates
    double precision :: fn  ! objective function value
    double precision, dimension (N) :: gradf    ! gradient wrt coefficient
    double precision :: lamreg  ! lambda1 within each lambda loop
    integer :: nlam = 1 ! number of lambda1's in calling coefficient update, set to 1
    double precision :: ftol, partol, abstol    ! tolerances for coefficients and weight update
    double precision :: fnorm   ! absolute difference between current and previous function values
    double precision :: parnorm ! sum of absolute difference between current and previous parameters
    double precision, dimension (N) :: param_old ! old parameter value
    double precision :: fn_old  ! old function value
    ! double precision, dimension (N) :: gradf_old ! old gradient value
    double precision :: resfn ! result function value
    double precision, dimension (N, 1) :: resgradf ! result gradient value
    double precision, dimension (N, 1) :: resparam ! result parameter value
    double precision, dimension (N, 1) :: resKa ! result K * param
    double precision, dimension (N) :: Ka   ! K * param, calculated in coefficient loop, used and updated by wt loop
    double precision, dimension (p) :: wt   ! weights
    logical, dimension (p) :: wtac  ! active set of weight
    integer :: wtac_size    ! size of active set of weight
    integer :: iter ! number of main iterations
    integer :: maxit    ! max iteration for parameter and weight loop
    double precision :: p2  ! penalty for weights, lam2*sum(weights)
    double precision :: fn_mid  ! fn after weight update before parameter update, for printing purpose
    logical :: sk   ! flag used to tell td_bfgs to initiate from previous param
    double precision, dimension (N) :: initp    ! initial parameter fed into td_bfgs
    integer :: conv ! convergence situation
    ! 10 - main loop goes over max_it
    ! 100 - fn increased after main loop iteration

    ! test
    ! double precision, dimension (p) :: gradw
    ! integer :: iter
    ! double precision :: gradwsq
    ! double precision :: wt_old (p)
    ! double precision :: gradw_old (p)
    ! double precision :: t=1.0d0

    ! content
    ftol = inner_ftol
    partol = inner_partol

    ! initialize abstol for parameter update
    ! fn after parameter update with new weights are not guaranteed to be smaller than that of the previous iteration
    abstol = 0.0d0

    ! maxit is max iteration for inner loops locally
    maxit = inner_maxit

    ! create dist array
    do i = 1, N
        do j = i, N
            dist(i,j,:) = x(i,:) - x(j,:)
            dist(j,i,:) = - dist(i,j,:)
        end do
    end do

    ! lambda loop
    do lid = 1, nval
        lam1 = lambda1(lid)
        lam2 = lambda2(lid)
        sigma = sigmas(lid)
        lamreg = lam1

        ! create initial K
        do i = 1, N
            do j = i, N
                K(i,j) = exp( -sigma * dot_product(dist(i,j,:), dist(i,j,:)) )  ! sum(dist(i,j,:)**2)
                K(j,i) = K(i,j)
            end do
        end do

        ! initialize weights
        wt = 1.0d0
        wtac = .true.
        wtac_size = p

        ! ! initialize a big fnorm and parnorm so first iteration runs
        ! fn = 1.0d5
        fnorm = 1.0d5
        parnorm = 1.0d5

        ! print initialization
        if (verbose) then
            call dblepr("lambda 1 = ", -1, lam1, 1)
            call dblepr("lambda 2 = ", -1, lam2, 1)
            call dblepr("sigma = ", -1, sigma, 1)
        end if

        ! Initial parameter update
        sk = .false.
        initp = 0.0d0
        p2 = lam2 * p
        call td_bfgs(N,K,y,lamreg,nlam,sk,initp,rho,ftol,partol,abstol,maxit,.false.,resfn,resgradf,resparam,resKa,conv)
        ! if (verbose > -1) call intpr("Convergence situation", -1, conv, 1)
        if (conv>=8) then
            call intpr("Program terminated at iteration:", -1, iter, 1)
            exit
        end if
        fn = resfn + p2
        gradf = resgradf(:,1)
        param = resparam(:,1)
        Ka = resKa(:,1)
        sk = .true.

        if (verbose) then
            call intpr("-----Summary of main loop iteration-----", -1, 0, 1)
            call dblepr("Objective function:", -1, fn, 1)
        end if

        iter = 1

        ! main loop - one weight loop + one param loop
        do while (fnorm > fn_tol .and. parnorm > param_tol .and. fn > abs_tol)

            fn_old = fn
            param_old = param

            ! if (verbose > -1) then
            !     call intpr("-----------Main loop iteration-----------", -1, iter, 1)
            !     call dblepr("-----Initialization for lambda 1-----", -1, lam1, 1)
            !     call dblepr("-----Initialization for lambda 2-----", -1, lam2, 1)
            !     call dblepr("-----Initialization for sigma-----", -1, sigma, 1)
            ! end if

            ! call weight update loop
            abstol = 0.0d0
            call sk_update(N,p,y,rho,sigma,lam1,lam2,dist,param,.false.,partol,abstol,maxit,K,Ka,wt,wtac,wtac_size,fn,conv)
            ! if (verbose > -1) call intpr("Convergence situation", -1, conv, 1)

            if (conv>=7) then
                ! call intpr("Inner weight loop terminated at max iteration:", -1, iter, 1)
                exit
            end if
            ! if (wtac_size == 0) then
            !     call intpr("All weights are zero. Program terminates at iteration:", -1, iter, 1)
            !     conv = 50
            !     exit
            ! end if

            fn_mid = fn

            ! call parameter update loop (BFGS)
            p2 = lam2 * sum(pack(wt, wtac))
            abstol = fn - p2
            initp = param_old
            call td_bfgs(N,K,y,lamreg,nlam,sk,initp,rho,ftol,partol,abstol,maxit,.false.,resfn,resgradf,resparam,resKa,conv)
            ! if (verbose > -1) call intpr("Convergence situation", -1, conv, 1)
            if (conv>=8) then
                call intpr("Inner alpha loop terminated at max iteration:", -1, iter, 1)
                exit
            end if
            fn = resfn + p2
            gradf = resgradf(:,1)
            param = resparam(:,1)
            Ka = resKa(:,1)

            fnorm = fn_old - fn
            parnorm = sum(abs(param_old - param))/N

            if (verbose) then
                call intpr("-----Summary of main loop iteration-----", -1, iter, 1)
                call intpr("Active weight count:", -1, wtac_size, 1)
                call dblepr("Objective function after weight update:", -1, fn_mid, 1)
                call dblepr("Objective function after parameter update:", -1, fn, 1)
                call dblepr("Decrease in objective function:", -1, fnorm, 1)
                call dblepr("Mean of absolute parameter update", -1, parnorm, 1)
            end if

            if (fnorm<0) then
                call intpr("Error: Fn increased after update iteration:", -1, iter, 1)
                conv = 100
                exit
            end if
            iter = iter + 1

            if (iter>max_it) then
                conv = 10
                call intpr("Failure to converge before main loop maximum iteration:", -1, max_it, 1)
                exit
            end if

        end do  ! end of main loop

        ! store output in *_final
        convergence(lid) = conv
        wt_final(:,lid) = wt
        fn_final(lid) = fn   ! need to change to lam2 * wtac_size later
        gradf_final(:,lid) = resgradf(:,1)
        param_final(:,lid) = resparam(:,1)
    end do  ! lid: lambda loop end



end subroutine td_sk

