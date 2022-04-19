module dimcorsingle_subs
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

    function philamlama(k,n,r,phi,lam1,lam2,a) result(phik)
        implicit none 
        integer, intent(in) :: k,n,a
        real(qp), intent(in) :: r, phi(2*n), lam1, lam2
        real(qp) :: xnew, ynew, phik

        xnew = min((point(k,n,r)+a)/lam1,(point(k,n,r)+a)/lam2)
        ynew = max((point(k+1,n,r)+a)/lam1,(point(k+1,n,r)+a)/lam2)

        !        phik = maxval(phi(ceiling(xnew*n/r)+n+1:floor(ynew*n/r)+n))
        phik = maxval(phi(floor(xnew*n/r)+n+1:ceiling(ynew*n/r)+n))

    end function philamlama

end module dimcorsingle_subs  

program dimcorsingle 

    use dimcorsingle_subs

    implicit none

    real(qp) :: dlam, delta !, salem(18)
    integer :: nintl, nit, nlam, lammin, lammax

    integer :: j, jj, jjj, testvar, jmin, jmax
    real(qp) :: r, lmin, lmax, ep, w, lambda, dimguess, theta

    real(qp), allocatable :: phi(:), phinew(:), dimdata(:,:)

    character(100) :: outfile
    character(8) :: myfile, file_ints

    logical :: stopp
 
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

    print *, 'Please enter the dimension guess'
    read *, dimguess 

    if ( dimguess < 10.0_qp**(-16) )  then
        print *, 'dimension guess seems too small. Stopping.'
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

    print *,'lambda_max', lambda+0.50_qp*dlam
    print *,'dimension guess', dimguess !-real(delta,qp) !  min(dimdata(j,2)-0.00001,dimdata(j,2)-0.150_qp*(1.0_qp - dimdata(j,2)))
    lmin = lambda - 0.50_qp*dlam ! lamleft+real(j-1,qp)*(lamright-lamleft)/nlam
    lmax = lmin+dlam
    w = 0.250_qp*lmin**(-dimguess+theta) ! (-min(dimdata(j,2)-0.00001,dimdata(j,2)-0.150_qp*(1.0_qp - dimdata(j,2))))    !   delta)

    r = (1.0_qp/(1.0_qp - lmax)+2.10_qp)/lmin 

    jmin = floor(-(1.0_qp/(1.0_qp - lmax)+1.0)*nintl/r)+nintl+1
    jmax = ceiling( (1.0_qp/(1.0_qp - lmax)+1.0)*nintl/r)+nintl
    print *, 'J_lambda', point(jmin,nintl,r), point(jmax,nintl,r)
    phi = 0.0_qp

    phi(jmin:jmax) = 1.0_qp !  floor(-(1.0_qp/(1.0_qp - lmax)+1.0)*nintl/r)+nintl+1:ceiling( (1.0_qp/(1.0_qp - lmax)+1.0)*nintl/r)+nintl) = 1.0_qp

    stopp = .false.
    jj = 1

    itloop: do jj = 1, nit 
        ep = 0.0_qp
        phinew = phi
        do jjj = jmin-1, jmax+1 !floor(-(1.0_qp/(1.0_qp - lmax)+1.0)*nintl/r)+nintl, ceiling( (1.0_qp/(1.0_qp - lmax)+1.0)*nintl/r)+nintl+1 
            ep =  w*(2.0_qp*philamlama(jjj,nintl,r,phi,lmin,lmax,0) + philamlama(jjj,nintl,r,phi,lmin,lmax,-1) + &
                philamlama(jjj,nintl,r,phi,lmin,lmax,1))
            !ep = w*ep  
            if (phi(jjj)>ep+theta) phinew(jjj) = ep  
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

    write(outfile,'(a,a)') file_ints,'.dat'
    open(10,file=outfile,status='replace')
    write(10,'("lambda = ",F16.12)') lambda
    write(10,'("dimension guess = ",F16.12)') dimguess
    write(10,'("Number of iterations: ", I6 )') jj
    write(10,'("Success? ", L4 )') stopp
    do jjj = 1,2*nintl
        write(10,'(F16.12 F16.12)')  point(jjj,nintl,r), phi(jjj)
    enddo
    close(10)


    deallocate(phi)
    deallocate(phinew)

end program dimcorsingle
