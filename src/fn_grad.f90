! Gradient of fn (without penalty on w) wrt param
! preparation - imput K, return Ka, Kexp1, Kexp2
subroutine fn_grad_prep(N, K, param, rho, Ka, Kexp1, Kexp2)
    implicit none
    ! Input variables
    integer, intent(in) :: N ! sample size
    double precision, dimension (N, N), intent(in) :: K ! Kernel matrix
    double precision, dimension (N), intent(in) :: param ! parameter estimate
    double precision, intent(in) :: rho ! index parameter
    ! Output variables
    double precision, dimension (N), intent(out) :: Ka ! K * param
    double precision, dimension (N), intent(out) :: Kexp1, Kexp2 ! exp[-(rho-1) K*alpha], exp[(2-rho) K*alpha]
    ! Local variables 
    double precision, dimension (N) :: ExpKa ! exp(Ka)

    ! Content
    Ka = matmul(K, param)
    ExpKa = exp(Ka)
    Kexp1 = ExpKa**(1.0d0 - rho)
    Kexp2 = ExpKa**(2.0d0 - rho)
end subroutine fn_grad_prep

! Calculate fn following fn_grad_prep
subroutine compute_fn(N, Ka, Kexp1, Kexp2, y, param, rho, lambda, fn)
    implicit none
    ! input variables
    integer, intent(in) :: N ! sample size
    ! double precision, dimension (N, N) :: K ! Kernel matrix
    double precision, dimension (N), intent(in) :: Ka ! K * param
    double precision, dimension (N), intent(in) :: Kexp1, Kexp2 ! exp[-(rho-1) K*alpha], exp[(2-rho) K*alpha]
    double precision, dimension (N), intent(in) :: y ! outcome variable
    double precision, dimension (N), intent(in) :: param ! parameter estimate
    double precision, intent(in) :: rho
    double precision, intent(in) :: lambda
    ! output variable
    double precision, intent(out) :: fn
    ! content
    fn = 1.0d0/N * (dot_product(y, Kexp1)/(rho - 1.0d0) + sum(Kexp2)/(2.0d0 - rho)) + lambda * dot_product(param, Ka)
end subroutine compute_fn

! calculate gradf following compute_grad
subroutine compute_grad(N, K, Ka, Kexp1, Kexp2, y, lambda, gradf)
    implicit none
    ! input variables
    integer :: N ! sample size
    double precision, dimension (N, N), intent(in) :: K ! Kernel matrix
    double precision, dimension (N), intent(in) :: Ka ! K * param
    double precision, dimension (N), intent(in) :: Kexp1, Kexp2 ! exp[-(rho-1) K*alpha], exp[(2-rho) K*alpha]
    double precision, dimension (N), intent(in) :: y ! outcome variable
    double precision, intent(in) :: lambda
    ! output variable
    double precision, dimension (N), intent(out) :: gradf
    ! content
    gradf = 1.0d0/N * reshape(matmul(K, (-y * Kexp1)) + matmul(K, Kexp2), (/N/)) + 2.0d0 * lambda * Ka
end subroutine compute_grad

! compute both fn and gradf using call to fn_grad_prep
subroutine fn_grad(N, K, y, param, lambda, rho, Ka, fn, gradf)
    implicit none
    ! input
    integer :: N ! sample size
    double precision, dimension (N, N), intent(in) :: K ! Kernel matrix
    double precision, dimension (N), intent(in) :: y ! outcome variable
    double precision, dimension (N), intent(in) :: param    ! parameter estimate
    double precision, intent(in) :: lambda  ! penalty parameter for coefficient
    double precision, intent(in) :: rho ! index parameter
    ! output
    double precision, intent(out) :: fn  ! function value
    double precision, dimension (N), intent(out) :: gradf    ! gradient for coefficient
    double precision, dimension (N), intent(out) :: Ka ! K * param
    ! inout
    ! local
    double precision, dimension (N) :: Kexp1, Kexp2 ! exp[-(rho-1) K*alpha], exp[(2-rho) K*alpha]

    ! Content
    call fn_grad_prep(N, K, param, rho, Ka, Kexp1, Kexp2)
    fn = 1.0d0/N * (dot_product(y, Kexp1)/(rho - 1.0d0) + sum(Kexp2)/(2.0d0 - rho)) + lambda * dot_product(param, Ka)
    gradf = 1.0d0/N * (matmul(K, (-y * Kexp1)) + matmul(K, Kexp2)) + 2.0d0 * lambda * Ka
end subroutine fn_grad

! ----------------------------------------------------------------
! Gradient of fn (with penalty on w) wrt w, only the partial derivatives of w in the active set are calculated.
! Output reverted to vector of length p
! RBF kernel only (with parameter sigma)
subroutine grad_wt(N, p, Ka, y, param, dist, rho, wt, wtac, wtac_size, sigma, lam1, lam2, gradw)
    implicit none
    ! input
    integer, intent(in) :: N   ! sample size
    integer, intent(in) :: p   ! dimension
    double precision, dimension (N), intent(in) :: Ka    ! K * param
    double precision, dimension (N), intent(in) :: y   ! outcome variable
    double precision, dimension (N), intent(in) :: param ! coefficients
    double precision, dimension (N, N, p) :: dist   ! x_i - x_j
    double precision, intent(in) :: rho    ! tweedie index parameter
    double precision, dimension (p), intent(in) :: wt  ! kernel weights
    logical, dimension (p), intent(in) :: wtac    ! active set of weight
    integer, intent(in) :: wtac_size    ! size of active weight
    double precision, intent(in) :: sigma  ! RBF kernel parameter
    double precision, intent(in) :: lam1    ! parameter penalty coefficient
    double precision, intent(in) :: lam2    ! weight penalty coefficient

    ! output
    double precision, dimension (p), intent(out) :: gradw  ! gradient of fn wrt w

    ! local
    integer :: i, j
    double precision, dimension (N, N, wtac_size) :: dist_ac    ! x_i - x_j with positions corresponding to wtac
    double precision, dimension (N, N, wtac_size) :: weighted_dist  ! x_i - x_j weighted by wt with 0's deleted
    double precision, dimension (N, N) :: cmat  ! c matrix: c_ij = -2sigma exp(-sigma ||wx_i - wx_j||^2)
    double precision, dimension (N, N, wtac_size) :: dKdw   ! derive K  wrt w
    double precision, dimension (N, N, wtac_size) :: dp1dw  ! par_i * par_j * dKdw_ij
    double precision, dimension (N, N) :: dl1dK, dl2dK  ! derivatives of l1 and l2 wrt K
    double precision, dimension (N, wtac_size) :: dl1, dl2, dp1 ! derivatives of l1, l2, p1 wrt w for each i
    double precision, dimension (wtac_size) :: gradw_tmp    ! temporary gradw that does not have 0's
    integer, dimension (wtac_size) :: wtac_id ! numeric vector for indexing active weights


    ! convert logical wtac to numeric wtac_id
    j = 0
    do i = 1, p
        if (wtac(i)) then
        j = j + 1
        wtac_id(j) = i
        end if
    end do

    do i = 1, N
        do j = i, N
            dist_ac(i,j,:) = dist(i,j,wtac_id)
            dist_ac(j,i,:) = - dist_ac(i,j,:)

            weighted_dist(i,j,:) = dist_ac(i,j,:) * wt(wtac_id)
            weighted_dist(j,i,:) = - weighted_dist(i,j,:)

            cmat(i,j) = -2.0d0 * sigma * exp( - sigma * sum( weighted_dist(i,j,:)**2 ) )
            cmat(j,i) = cmat(i,j)

            dKdw(i,j,:) = cmat(i,j) * dist_ac(i,j,:) * weighted_dist(i,j,:)
            dKdw(j,i,:) = dKdw(i,j,:)

            dp1dw(i,j,:) = param(i) * param(j) * dKdw(i,j,:)
            dp1dw(j,i,:) = dp1dw(i,j,:)
        end do
        dl1dK(i,:) = -y(i) * exp((1.0d0-rho) * Ka(i)) * param
        dl1(i,:) = matmul(dl1dK(i,:), dKdw(i,:,:))

        dl2dK(i,:) = exp((2.0d0-rho) * Ka(i)) * param
        dl2(i,:) = matmul(dl2dK(i,:), dKdw(i,:,:))
    end do

    dp1 = sum(dp1dw, dim = 2)
    gradw_tmp = 1.0d0/N * sum(dl1, dim = 1) + 1.0d0/N * sum(dl2, dim = 1) + lam1 * sum(dp1, dim = 1) + lam2

    gradw = 0.0d0
    do i = 1, wtac_size
        gradw(wtac_id(i)) = gradw_tmp(i)
    end do

end subroutine grad_wt