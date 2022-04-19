!Computing uniform lower bound on correlation dimension, \lambda > 1/3
module deletedig_subs 
    USE OMP_LIB
    implicit none

    integer, parameter :: DRK = selected_real_kind(14)
    integer, parameter :: qp = selected_real_kind(33,4931)

    type MAX_TREE
        type(MAX_TREE), pointer :: next
        integer  :: n
        real(qp), allocatable :: x(:)
    end type MAX_TREE
contains
    recursive subroutine mtree_create( mtree, n )
        type(MAX_TREE), pointer  :: mtree
        integer, intent(in) :: n
        allocate( mtree )    
        mtree%n  =  n
        allocate(mtree%x(n+1))
        mtree%x  =  0.0_qp
        if (n==1) then
            mtree%next  => null()
        else
            call mtree_create(mtree%next, (n+1)/2)       
        end if
    end subroutine mtree_create
    recursive subroutine mtree_destroy( mtree )
        type(MAX_TREE), pointer  :: mtree
        type(MAX_TREE), pointer  :: next
        next  => mtree%next
        if ( associated(next) ) then
            call mtree_destroy(next)
        endif
        deallocate( mtree%x)
        deallocate( mtree )
    end subroutine mtree_destroy
    subroutine mtree_append( mtree, x, n)
        type(MAX_TREE), pointer  :: mtree
        integer, intent(in) :: n
        real(qp), intent(in) :: x(n)
        mtree%x = x
        call mtree_update(mtree, 1, n)
    end subroutine mtree_append
    subroutine mtree_append2( mtree, x, n, j0, j1)
        type(MAX_TREE), pointer  :: mtree
        integer, intent(in) :: n, j0, j1
        real(qp), intent(in) :: x(n)
        mtree%x(j0:j1) = x(j0:j1)
        call mtree_update(mtree, j0, j1)
    end subroutine mtree_append2
    recursive subroutine mtree_update(mtree, j0, j1)
        type(MAX_TREE), pointer  :: mtree
        integer, intent(in) :: j0, j1
        integer :: n, i, j00, j11
        real(qp):: m
        type(MAX_TREE), pointer  :: next
        if (mtree%n > 1) then
            if (j1 > j0) then
                !    if ( associated(next) ) then
                j00 = 0
                j11 = 0
                next  => mtree%next
                do i = (j0+1)/2,(j1+1)/2
                    m = max(mtree%x(i*2),mtree%x(i*2-1))
                    if (m /= next%x(i)) then
                        next%x(i) = m
                        if (j00==0) j00=i
                        j11 = i
                    end if
                end do
                if (j00 /= 0) call mtree_update(next,j00,j11)
            else
                if (j1 == j0) then
                    next  => mtree%next
                    i = (j0+1)/2
                    m = max(mtree%x(i*2),mtree%x(i*2-1))
                    if (m /= next%x(i)) then
                        next%x(i) = m
                        call mtree_update(next,i,i)
                    end if
                end if
            end if
        end if
    end subroutine mtree_update
    recursive function mtree_get(mtree,j0,j1) result(m)
        type(MAX_TREE), pointer  :: mtree
        integer, intent(in) :: j0,j1
        integer :: j00,j11
        real(qp) :: m
        m = 0.0_qp
        if (j1 >= j0) then
            !    if ( associated(next) ) then
            if (mtree%n > 1) then
                m = max(mtree%x(j0),mtree%x(j1))
                j00=j0/2+1
                j11=j1/2
                if (j11 > j00) then
                    m = max(m,mtree_get(mtree%next,j00,j11))
                else
                    m = max(m,mtree%next%x(j00))
                end if
            else
                if (j0==1) m = mtree%x(1)
            endif    
        end if
    end function mtree_get

    function finv(x,a,lam) result(y)
        ! returns inverse image under f_{k,\lambda}(x) = \lambda x + a
        implicit none
        real(qp), intent(in) :: x,lam
        integer, intent(in) :: a
        real(qp) :: y 

        y = (x-a)/lam

    end function finv

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

    subroutine onelam(jjjj,initvals,dlam,theta,vareps,dimdrop,nm,nintl,nit,itnum, stepc, tphi)
        implicit none
        ! f_{\lambda,k} = \lambda x + k 
        ! r = the length of the interval = (-4/(1-\lambda),4/(1-\lambda)) 
        ! epsilon = the gap 
        ! alpha = dimension guess 
        ! nintl = the number of partition intervals
        ! nit = the number of iterations
        ! ep = end points of partition intervals
        ! nm = number of different maps
        real(qp), intent(in) :: initvals(3),dlam,theta, vareps,dimdrop
        real(qp) :: lambda, alpha, lammin, lammax
        integer :: nit, nintl, nlam, nalph, nlam0
        integer, allocatable :: indl(:,:), indr(:,:)
        real(qp), allocatable :: phi(:), phinew(:), ep(:)
        real(qp) :: lmin, lmax, w, x, r, almax, almin
        integer :: j, jj, jjj, jjjj, j5, j0, j1 , nm, j6
        real(qp), intent(inout) :: tphi
        integer, intent(inout) :: itnum, stepc
        type(MAX_TREE), pointer  :: mtree
        logical :: stopp 
        real*4 sec, sec1

        allocate(phi(2*nintl))
        !    allocate(phinew(2*nintl))
        allocate(indl(nm,2*nintl))
        allocate(indr(nm,2*nintl))
        call mtree_create( mtree, 2*nintl)
        phi = 0.0_qp
        tphi = 0.0_qp
        itnum = 0
        stepc = 0

        sec = secnds(0.0)

        lambda = initvals(1) ! lammin+(jjjj-1)*dlam

        lmin = lambda - 0.60_qp*dlam
        lmax = lambda + 0.60_qp*dlam

        r = (3.0_qp/(1.0_qp-lmax)+4.0_qp)!/lambda

        j0 = floor(-(3.0_qp/(1.0_qp-lmax)+1.0_qp)*nintl/r)+nintl
        j1 = ceiling((3.0_qp/(1.0_qp-lmax)+1.0_qp)*nintl/r)+nintl

        phi(j0:j1) = 1.0_qp

        do j = 1,nintl*2
            do jj = 1,nm
                indl(jj,j) = min(leftind(finv(point(j,nintl,r),jj-4,lmin),nintl,r), &
                    leftind(finv(point(j,nintl,r),jj-4,lmax),nintl,r))
                if (indl(jj,j) < 1)   indl(jj,j) = 1
                if (indl(jj,j) > 2*nintl)   indl(jj,j) = 2*nintl


                indr(jj,j) = max(rightind(finv(point(j,nintl,r),jj-4,lmin),nintl,r), &
                    rightind(finv(point(j,nintl,r),jj-4,lmax),nintl,r))  
                if (indr(jj,j) < 1)   indr(jj,j) = 1
                if (indr(jj,j) > 2*nintl)   indr(jj,j) = 2*nintl

            enddo
        enddo

        do j=1,nintl*2-1
            do jj=1,nm 
                if (indl(jj,j) > indr(jj,j+1)) then 
                    print *,jjjj, 'error: reversed image'
                endif
            enddo
        enddo

        almin = initvals(2)-dimdrop
        almax = initvals(3)
        nalph = ceiling(sqrt((almax-almin)/vareps))+1
        print *, jjjj, 'lambda=', lambda, 'taking partition into', nalph, 'intervals'

        do while (almax - almin>vareps)
            partloop: do j5 = 1,nalph
                alpha = almin+j5*(almax-almin)/real(nalph,qp)
                phi = 0.0_qp
                phi(j0:j1) = 1.0_qp
                stopp = .false.
                stepc = stepc + 1

                w = lmin**(-alpha)/9.0_qp
                print *,jjjj, 'Scaling factor: ', w
                j = 1

                !                phinew = phi 
                call mtree_append( mtree, phi, 2*nintl)
                itloop: do j = 1,nit 
                    do jj = 1,2*nintl-1
                        x = 0.0_qp
                        do jjj=1,3
                            !                            x = x + maxval(phi(indl(jjj,jj):indr(jjj,jj+1)))+maxval(phi(indl(jjj+4,jj):indr(jjj+4,jj+1)))+ &
                            !                                maxval(phi(indl(4,jj):indr(4,jj+1)))
                            x = x + mtree_get(mtree,indl(jjj,jj),indr(jjj,jj+1)) + mtree_get(mtree,indl(jjj+4,jj),indr(jjj+4,jj+1))& 
                                +mtree_get(mtree,indl(4,jj),indr(4,jj+1))
                        enddo
                        if ( x*w + theta < phi(jj)) then 
                            !                            phinew(jj) = x*w + eps 
                            phi(jj) = x*w + theta 
                        endif
                    enddo
                    call mtree_append( mtree, phi, 2*nintl)    
                    !                    phi = phinew

                    !                    if (maxval(phi)<0.9850_qp) then 
                    if (mtree_get(mtree,0,2*nintl)<0.9990_qp) then 
                        print *, jjjj, 'win', j, lambda, alpha
                        tphi = alpha 
                        itnum = j  
                       
                        if (almax - alpha < vareps ) then
                            print *, jjjj, 'last refinement was useful: ', j, stepc, tphi
                            goto 111
                        endif
                      
                        stopp = .true.
                        exit itloop
                    endif
                enddo itloop
                if (stopp.eqv..false.) then 
                    almin = alpha - (almax-almin)/nalph
                    almax = alpha 
                    print *, jjjj, almin, almax, j
                    exit partloop
                endif
            enddo partloop
        enddo ! while
        if (stopp.eqv..false.) then 
            print *, jjjj, 'left while loop:', itnum, stepc, tphi
        endif
        111     continue

        sec1 = secnds(sec)
        print *, jjjj, stopp, ':', sec1, 'seconds to run'
        deallocate(phi)
        deallocate(indl)
        deallocate(indr)
        !        deallocate(phinew)        
        call mtree_destroy( mtree )
    end subroutine onelam

end module deletedig_subs



program deleteddigit
    USE OMP_LIB
    use deletedig_subs
    implicit none
    ! f_{\lambda,k} = \lambda x + k 
    ! r = the length of the interval = (-4/(1-\lambda),4/(1-\lambda)) 
    ! theta = the gap 
    ! alpha = dimension guess 
    ! nintl = the number of partition intervals
    ! nit = the number of iterations
    ! ep = end points of partition intervals
    ! nm = number of maps
    real(qp) :: lambda, alpha, theta, lammin, lammax, dimdrop, vareps, tphi1
    integer :: nit, nintl, total, nlam, nalph, nlam0, itnum1, stepc1
    integer, allocatable :: itnum(:), stepc(:)
    real(qp), allocatable :: ep(:), tphi(:), initvals(:,:)
    real(qp) :: lmin, lmax, w, x, r, dlam, almax, almin
    integer :: j, jj, jjj, jjjj, j5, j0, j1 , nm

    logical :: stopp 
    character(100) :: infile, outfile
    real*4 sec, sec1

   print *, 'Please give the  theta > 10^{-24} '
   read *, theta
   
   if ( theta < 10.0_qp**(-24))  then
        print *, 'Theta too small. Stopping'
        stop
   endif
   
   print *, 'Please give the total number of intervals in partition Lambda provided'
   read *, total
   
   if ( total < 2)  then
        print *, 'At least one interval is necessary. Stopping.'
        stop
   endif

   print *, 'Please give the number of the first end point. '
   read *, nlam0

   if ( nlam0 > total)  then
        print *, 'The number of the first value is too large. Stopping.'
        stop
   endif

   print *, 'Please give the number of Lambda-intervals to consider'
   read *, nlam
 
   if ( nlam > total)  then
        print *, 'The number of the intervals is too large. Stopping.'
        stop
   endif

   print *, 'Please give the number of partition intervals for test function'
   read *, nintl     
 
   if ( nintl < 1)  then
        print *, 'The number of the partition intervals is too small. Stopping.'
        stop
   endif

   print *, 'Please give the number of iterations for the diffusion operator'
   read *, nit

   if ( nit < 1)  then
        print *, 'The number of the iterations intervals is too small. Stopping.'
        stop
   endif

   print *, 'Please give the refinement parameter  10^{-24} < epsilon < 1'
   read *, vareps
 
   if ( vareps < 10**(-24))  then
        print *, 'The refinement parameter is too small. Stopping.'
        stop
   endif
     
   print *, 'Please enter filename for the output data'
   read (5,*) outfile    
   print *, 'The result will be written into the file ', outfile    
   
   print *, 'Please enter filename for the initial data '
   read (5,*) infile
   print *,'I will read the data from the file ', infile    


    dimdrop = 0.005
    
    allocate(initvals(total,3))

    open(12,file = trim(infile)) 
    read(12,*) ((initvals(j,jj), jj=1,3), j=1,total)
    close(12)

    dlam = 0.50_qp*(initvals(1,1) - initvals(3,1))

    open(10,file=outfile,status='replace')
    write(10, '(" lambda ", 5X, " number of iterations ", 5X, " alpha ", 5X,  " result ")')
    write(10,'("number of iterations = ",I6)') nit
    write(10,'("partition intervals = ",I6)') nintl
    write(10,'("the length of lambda-intervals is = ", F16.12 )') 1.20_qp*dlam
    write(10, '("#", 5X, " lambda ", 5X, " number of iterations ", 5X, " dimension ")')    

    allocate(stepc(nlam+1))
    allocate(tphi(nlam+1))
    allocate(itnum(nlam+1))
  
    nm = 7
    
    stepc = 0
    itnum = 0 
    tphi = 0
    !$OMP PARALLEL SHARED(itnum, stepc, tphi) PRIVATE(j,itnum1,stepc1,tphi1)
    !    print *,'Threads=',OMP_GET_MAX_THREADS()
    !$OMP DO SCHEDULE (DYNAMIC, 1)
    do j = nlam0,nlam
        print *,'Thread #',1+OMP_GET_THREAD_NUM(),'/',OMP_GET_MAX_THREADS(),'j=',j
        call onelam(j,initvals(j,1:3),dlam,theta,vareps,dimdrop,nm,nintl,nit,itnum1,stepc1,tphi1)
        !        print *,'Output',j,itnum1,stepc1,tphi1
        itnum(j-nlam0+1)=itnum1
        stepc(j-nlam0+1)=stepc1
        tphi(j-nlam0+1)=tphi1
        !write(10,17) j, initvals(j,1), itnum1, tphi1
    enddo
    !$OMP END PARALLEL
    17     format(i5,4x,f16.12,4x,i6,4x,f16.12)   
    do j = 1,nlam
        write(10,17) j, initvals(j,1), itnum(j), tphi(j)
    end do
    close(10)
    deallocate(tphi)
    deallocate(itnum)
    deallocate(stepc)
    deallocate(initvals)

end program deleteddigit



