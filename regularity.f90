module regularity_subs
    implicit none 

    integer, parameter :: DRK = selected_real_kind(14)
    integer, parameter :: qp = selected_real_kind(33, 4931)
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
    real(qp) :: m
    m = 0.0_qp
    if (j1 >= j0) then
!    if ( associated(next) ) then
    if (mtree%n > 1) then
        m = mtree_get(mtree%next,j0/2+1,j1/2)
        if (mod(j0,2) /= 1) m = max(mtree%x(j0),m)
        if (mod(j1,2) == 1) m = max(mtree%x(j1),m)
    else
        if (j0==1) m = mtree%x(1)
    endif    
    end if
end function mtree_get

    
    function rightind(x,n,r,lmin) result(kr) 
        ! returns the number of right end point of the interval of the partition containing x
        implicit none 
        real(qp), intent(in) :: x,r,lmin
        integer :: kr 
        integer, intent(in) :: n

        !kr = ceiling(x*n/r)+n+1
        kr = ceiling(real(n-1,qp)*(x+2/lmin)/r) + 1
    end function rightind

    function leftind(x,n,r,lmin) result(kl) 
        ! returns the number of the left end point of the interval of the partition containing x
        implicit none 
        real(qp), intent(in) :: x,r,lmin
        integer :: kl 
        integer, intent(in) :: n

        !kl = floor(x*n/r)+n+1
        kl = floor(real(n-1,qp)*x/r + 2.0_qp*real(n-1,qp)/(lmin*r)) + 1
    end function leftind 

    function point(k,n,r,lmin) result(x)
        ! returns point of the partition number k 
        implicit none 
        real(qp) :: x,r,lmin
        integer, intent(in) :: k,n

        !x = r*(k-n-1)/real(n,qp)
        x = real(k-1,qp)*r/real(n-1,qp)-2.0_qp/lmin
    end function point 

    function philamlama(k,n,r,mtree,lam1,lam2,jmin,jmax,a) result(phik)
        implicit none 
        integer :: j00,j11
        integer, intent(in) :: k,n,a,jmin,jmax
        real(qp), intent(in) :: r, lam1, lam2
        real(qp) :: xnew, ynew, phik
        type(MAX_TREE), pointer, intent(in)  :: mtree

        xnew = min((point(k,n,r,lam1)+a)/lam1,(point(k,n,r,lam1)+a)/lam2)
        ynew = max((point(k+1,n,r,lam1)+a)/lam1,(point(k+1,n,r,lam1)+a)/lam2)

        phik = 0.0_qp
        j00 = max(jmin, leftind(xnew,n,r,lam1))
        
        j11 = min(jmax, rightind(ynew,n,r,lam1))
        
        phik = mtree_get(mtree,j00,j11)
    end function philamlama
    
    function philamlama2(k,n,r,mtree,lmin,lmax,j0,j1) result(ep)
        implicit none 
        integer, intent(in) :: k,n,j0,j1
        real(qp), intent(in) :: r, lmin,lmax
        type(MAX_TREE), pointer, intent(in)  :: mtree
        real(qp) :: xnew, ynew, ep
        ep = philamlama(k,n,r,mtree,lmin,lmax,j0,j1,0) + philamlama(k,n,r,mtree,lmin,lmax,j0,j1,-1) 
        !+ &   !philamlama(k,n,r,mtree,lmin,lmax,j0,j1,1)
    end function philamlama2
    
    subroutine onelam(j,initvals,dlam,theta,vareps,nintl,nit,itnum, stepc, tphi)
    implicit none
    integer :: jj, jjj, jjjj, j0, j1, nalph 
    real(qp) :: r, lmin, lmax, ep, w, almin, almax, alpha, phik, dimdrop 
    real(qp), allocatable :: phi(:)
    real(qp), intent(in) :: initvals(2), dlam,  theta, vareps
    real(qp), intent(inout) :: tphi
    integer, intent(inout) :: itnum, stepc
    integer, intent(in) :: nintl, nit, j
    type(MAX_TREE), pointer  :: mtree
    logical :: stopp
    real*4 sec, sec1
    allocate(phi(nintl))
    call mtree_create( mtree, nintl)
    phi = 0.0_qp
    tphi = 0.0_qp
    itnum = 0
    stepc = 0
    sec = secnds(0.0)
        print *,j,'lambda_max', initvals(1)
        lmin = initvals(1) - 0.60_qp*dlam
        lmax = lmin + 1.20_qp*dlam

        almin = initvals(2)
        almax = 1.0_qp
        nalph = ceiling(sqrt((almax - almin)/vareps))+1
        
        print *,j,'taking partition into', nalph, 'intervals'
        do while (almax - almin > vareps)
            partloop: do jjjj = 1,nalph
                alpha = almin+jjjj*(almax-almin)/real(nalph,qp) 
                stopp = .false.
                stepc = stepc + 1
                w = 0.50_qp*lmin**(-alpha) 

                r = (1.0_qp/(1.0_qp - lmax)+4.00_qp)/lmin
                
                j0 = floor((nintl-1)*(2.0_qp/lmin - 1.0_qp)/r)+1  !floor(-(1.0_qp/(1.0_qp - lmax)+1.0)*nintl/r)+nintl+1
                j1 = ceiling((nintl-1)*(2.0_qp/lmin + 1.0_qp/(1.0_qp-lmax) + 1.0_qp)/r)+1  !ceiling( (1.0_qp/(1.0_qp - lmax)+1.0)*nintl/r)+nintl

               ! print *, 'bounds:', j0+1,j1
                phi = 0.0_qp
                phi(j0+1:j1) = 1.0_qp
                call mtree_append( mtree, phi, nintl)
                jj = 1
                itloop: do jj = 1, nit !while ((jj < nit).and.(stopp.eqv..false.)) 
!                    ep = 0.0_qp
                    do jjj = j0,j1+1 
                        phik = philamlama2(jjj,nintl,r,mtree,lmin,lmax,j0,j1+1)*w+theta
                        if (phik<phi(jjj)) then 
                        phi(jjj)=phik
!                        call mtree_append2( mtree, phi, 2*nintl, jjj, jjj)
                        end if
                    enddo
                    call mtree_append2( mtree, phi, nintl, j0, j1+1)
                    phik = mtree_get(mtree,j0,j1+1)
                    if (phik < 0.99900_qp) then 
                        itnum = jj
                        tphi = alpha
                        print *,j,'win', itnum, tphi
                        if (almax-alpha< vareps) then
                            print *,j, 'last step was useful:', itnum, stepc, tphi, phik
                            goto 111
                        endif
                        stopp = .true.
                        exit itloop
                    endif
                enddo itloop
                if (stopp.eqv..false.) then
                    almin = alpha - (almax-almin)/nalph 
                    almax = alpha
                    print *,j,almin, almax, jj
                    exit partloop
                endif
            enddo partloop
            enddo ! while 
            if (stopp.eqv..false.) then
                print *,j,'Left while loop:', itnum, stepc, tphi, maxval(phi)
            endif
111     continue            
        deallocate(phi)
        call mtree_destroy( mtree )
        sec1 = secnds(sec)
        print *,j,itnum,stepc,tphi, ':', sec1, 'seconds to run'
    end subroutine onelam

end module regularity_subs  

program regularity
USE OMP_LIB

    use regularity_subs
    implicit none
    real(qp) :: dlam, theta,tphi1,vareps !, salem(18)
    integer :: nintl, nit, nlam, nlam0 = 1, nalph,j, jj,itnum1, stepc1
    real(qp), allocatable :: tphi(:), initvals(:,:)
    integer, allocatable :: itnum(:), stepc(:)
    character(100) :: infile, filename, myfile
   
   print *, 'Please give the  theta > 10^{-24} '
   read *, theta
   
   if ( theta < 10.0_qp**(-24))  then
        print *, 'Theta too small. Stopping'
        stop
   endif
   
   print *, 'Please give the number of Lambda-intervals to consider'
   read *, nlam
 
   if ( nlam < 1)  then
        print *, 'The number of the lambda-values should be positive. Stopping.'
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

   print *, 'Please give the size of the Lambda-neighbourhoods'
   read *, dlam

   if ( dlam < 10.0_qp**(-24) )  then
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
   read (5,*) myfile    
   print *, 'The result will be written into the file ', myfile    
   
   print *, 'Please enter filename for the initial data '
   read (5,*) infile
   print *,'I will read the data from the file ', infile    

    print *,'theta=',theta,'nlam=',nlam,'nintl=',nintl,'nit=',nit

    allocate(tphi(nlam))
    allocate(itnum(nlam))
    allocate(stepc(nlam))
    tphi = 0.0_qp !.false.
    itnum = 0
    stepc = 0 
    
    allocate(initvals(nlam,2))
    open(12, file=trim(infile))
    read(12,*) ((initvals(j,jj), jj=1,2), j=1,nlam)
    close(12) 
    
    print *, 'The length of lambda-intervals is', 1.20*dlam !0.150_qp/nlam

    !$OMP PARALLEL SHARED(itnum, stepc, tphi) PRIVATE(j,itnum1,stepc1,tphi1)
!    print *,'Threads=',OMP_GET_MAX_THREADS()
    !$OMP DO SCHEDULE (DYNAMIC, 1)
    do j = nlam0,nlam
        print *,'Thread #',1+OMP_GET_THREAD_NUM(),'/',OMP_GET_MAX_THREADS(),'j=',j
        call onelam(j,initvals(j,1:2),dlam,theta,vareps,nintl,nit,itnum1,stepc1,tphi1)
!        print *,'Output',j,itnum1,stepc1,tphi1
        itnum(j)=itnum1
        stepc(j)=stepc1
        tphi(j)=tphi1
    enddo
    !$OMP END PARALLEL

        open(7,file=trim(myfile),status='replace')
        !write(7,'("Dimension Guess = ",F16.12)'), dimguess
        do j = nlam0,nlam
            write(7,17) j,initvals(j,1), stepc(j), itnum(j), tphi(j)
17     format(i5,1x,f16.12,1x,2i5,1x,f16.12)
        enddo
        close(7)
        call sleep(3)
        deallocate(tphi)    
        deallocate(itnum)
        deallocate(stepc)
    end program regularity
