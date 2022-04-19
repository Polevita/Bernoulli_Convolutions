module singlelam_subs
    implicit none 

    integer, parameter :: DRK = selected_real_kind(14)
    integer, parameter :: qp = selected_real_kind(33, 4931)
contains 

    function rightind(x,n,r) result(kr) 
        ! returns the number of right end point of the interval of the partition containing x
        implicit none 
        real(qp), intent(in) :: x,r
        integer :: kr 
        integer, intent(in) :: n

        kr = ceiling(x*n/r)+n+1
    end function rightind

    function leftind(x,n,r) result(kl) 
        ! returns the number of the left end point of the interval of the partition containing x
        implicit none 
        real(qp), intent(in) :: x,r
        integer :: kl 
        integer, intent(in) :: n

        kl = floor(x*n/r)+n+1
    end function leftind 

    function point(k,n,r) result(x)
        ! returns point of the partition number k 
        implicit none 
        real(qp) :: x,r
        integer, intent(in) :: k,n

        x = r*(k-n-1)/real(n,qp)
    end function point 

    function philamlama(k,n,r,phi,lam1,lam2,jmin,jmax,a) result(phik)
        implicit none 
        integer, intent(in) :: k,n,a,jmin,jmax
        real(qp), intent(in) :: r, phi(2*n), lam1, lam2
        real(qp) :: xnew, ynew, phik
        integer :: j0,j1

        xnew = min((point(k,n,r)+a)/lam1,(point(k,n,r)+a)/lam2)
        ynew = max((point(k+1,n,r)+a)/lam1,(point(k+1,n,r)+a)/lam2)

        phik = 0.0_qp
        j0 = max(jmin,floor(xnew*n/r)+n+1)
        j1 = min(jmax,ceiling(ynew*n/r)+n)

        if (j0<=j1) then 
            phik = maxval(phi(j0:j1)) !
        endif 
    end function philamlama

end module singlelam_subs  

program singlelam

    use singlelam_subs

    implicit none

    real(qp) :: dlam, delta, theta
    integer :: nintl, nit, nlam, lammin, lammax

    integer :: j, jj, jjj, jjjj, jmin, jmax, testvar, itnum
    real(qp) :: r, lmin, lmax, ep, w, dimguess, lambda

    real(qp), allocatable :: phi(:), phinew(:), dimdata(:,:), tphi(:)
    ! integer, allocatable :: itnum(:) 

    character(100) :: infile, filename
    character(8) :: myfile, file_ints

    logical :: stopp

    real*4 :: sec,sec1
 
   print *, 'Please give the theta > 10^{-24}'
   read *, theta
   
   if ( theta < 10.0_qp**(-24))  then
        print *, 'Theta too small. Stopping'
        stop
   endif
 
   print *, 'Please enter the size of neighbourhood of Lambda'
   read *, dlam
   if ( dlam < 10.0_qp**(-24))  then
        print *, 'Lambda is too small. Stopping.'
        stop
    endif
   
    print *, 'Please enter the value of lambda > 0.49'
    read *, lambda 

    if ( lambda < 0.490_qp)  then
        print *, 'Lambda is too small. Stopping.'
        stop
    endif

    print *, 'Please enter the your guess for the regularity exponent'
    read *, dimguess 

    if ( dimguess < 10.0_qp**(-16) )  then
        print *, 'The guess seems too small. Stopping.'
        stop
    endif

    print *, 'Please enter the number of partition intervals for test function'
    read *, nintl 

    if ( nintl < 1)  then
        print *, 'The number of the partition intervals is too small. Stopping.'
        stop
    endif

    print *, 'Please enter the number of iterations for the diffusion operator'
    read *, nit
    
    if ( nit < 1)  then
        print *, 'The number of the iterations intervals is too small. Stopping.'
        stop
    endif

    print *, 'Please enter filename for the function data (8 characters)'
    read (5,*) file_ints

    allocate(phi(2*nintl))
    allocate(phinew(2*nintl))
    phi = 0.0_qp

    print *, 'The length of lambda-intervals is', 1.20*dlam 
    itnum = 0
    sec = secnds(0.0)
    lmin = lambda - 0.50_qp*dlam
    lmax = lmin + dlam
    r = (1.0_qp/(1.0_qp - lmax)+2.10_qp)/lmin 

    jmin = floor(-(1.0_qp/(1.0_qp - lmax)+1.0)*nintl/r)+nintl+1
    jmax = ceiling( (1.0_qp/(1.0_qp - lmax)+1.0)*nintl/r)+nintl
    stopp = .false.
   
    print *, 'J_lambda:', point(jmin,nintl,r), point(jmax,nintl,r) 

    w = 0.50_qp*lmin**(-dimguess) 

    phi = 0.0_qp
    phi(jmin:jmax) = 1.0_qp 

    jj = 1
    itloop: do jj = 1, nit 
        phinew = phi
        do jjj = jmin-1, jmax+1 ! 
            ep = philamlama(jjj,nintl,r,phi,lmin,lmax,jmin,jmax,0) + philamlama(jjj,nintl,r,phi,lmin,lmax,jmin,jmax,-1)
            !  philamlama(jjj,nintl,r,phi,lmin,lmax,jmin,jmax,1)
            ! phi (jjj) = min(phi(jjj),ep*w+theta)
            ep = ep*w + theta
            if (ep < phi(jjj)) phinew(jjj) = ep 
        enddo
        phi = phinew
        phinew = 0.0_qp 

        if (maxval(phi)<0.99900_qp) then 
            print *, 'Correct. Confirmed with', jj, ' iterations' 
            stopp = .true.
            exit itloop
        endif
    enddo itloop
  
    if (stopp.eqv..false.) then 
        print *, 'I cannot validate your guess with these parameters. Sorry.'
    endif
    
    write(filename,'(a,a)') file_ints,'.dat'
    open(10,file=filename,status='replace')
    write(10,'("lambda = ",F16.12)') lambda
    write(10,'("Number of iterations: ", I6 )') jj
    write(10,'("Dimension ", F16.12 )') dimguess 
    do jjj = 1,2*nintl
        write(10,'(F16.12 F16.12)')  point(jjj,nintl,r), phi(jjj)
    enddo
    close(10)

    sec1 = secnds(sec)
    print *,j, ':', sec1, 'seconds to run'

    call sleep(3)

    deallocate(phi)
    deallocate(phinew)

end program singlelam
