! COPYRIGHT (C) 2014 The Science and Technology Facilities Council (STFC)
! Authors:  Mario Arioli      STFC
!           Jonathan Hogg     STFC
!
! [x, info, rinfo] = ma57_solve(A, b, struct, nitre)
! A is the sparse matrix
! b is an n x nrhs vector of rhs
! struct is the factors returned by ma57_factor()
! nitre is OPTIONAL #steps of IR
!
! x is n x nrhs soln vector
! info is OPTIONAL statistics (integer)
! rinfo is OPTIONAL statistics (real)
!
SUBROUTINE mexFunction( nlhs, plhs, nrhs, prhs )
      use hsl_matlab
      implicit none

      INTEGER * 4 :: nlhs, nrhs
      integer(mwp_) :: plhs(*), prhs(*)

      ! Argument positions
      integer, parameter :: A_in       = 1, &
                            b_in       = 2, &
                            struct_in  = 3, &
                            nitre_in   = 4
      integer, parameter :: x_out      = 1, &
                            info_out   = 2, &
                            rinfo_out  = 3

      ! Interfaces
      interface
         double precision function mxGetNaN()
         end function mxGetNaN
      end interface

      ! Local variables

      INTEGER, DIMENSION(:), ALLOCATABLE :: ip_in,irn_in,jcn,irn
      INTEGER, DIMENSION(:), ALLOCATABLE :: ifact,icntl,iwork,info

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fact,work,val
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xxx,rhs
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: val_in,resid
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: cntl,rinfo

      INTEGER :: ne,err,nitre
      integer :: n, m, n_rhs
      integer :: lwork

      integer(mws_) :: mwm, mwn, mwnrhs

      DOUBLE PRECISION :: start,finish

      ! Check number of arguments
      if ( nrhs<3 .or. nrhs>4 ) &
        CALL matlab_error( " Wrong # of input arguments " )
      if ( (nlhs > 3) ) &
        CALL matlab_error( " Wrong # of output arguments " )

      ! Read in arguments
      call matlab_to_fortran(prhs(A_in), mwm, mwn, ip_in, irn_in, val_in, 'A')
      if (mwm.ne.mwn) &
        CALL mexErrMsgTxt("The matrix must be square")
      m = mwm
      n = mwn
      call matlab_to_fortran(prhs(b_in), xxx, mwn, mwnrhs, 'b')
      n_rhs = mwnrhs
      if(mwn.ne.n) &
        CALL mexErrMsgTxt( &
          "The RHS has dimensions incompatible with the input matrix")
      call read_struct(prhs(struct_in), fact, ifact, cntl, icntl, info, rinfo)
      if(info(1).lt.0) &
         call matlab_error('Previous factorization failed')

      ! Set desired number of IR steps
      if (nrhs .ge. nitre_in) then
        call matlab_to_fortran(prhs(nitre_in),nitre,'nitre')
        icntl(9) = nitre
      end if
      icntl(9) = max(0, icntl(9))
      if (icntl(9) .le. 0) &
        CALL matlab_warning("Iterative refinement is not executed")

      ! Allocate workspace
      lwork = max(n*n_rhs, 4*n)
      ALLOCATE( iwork(n), work(lwork), stat=err)
      IF (err .ne. 0) &
        CALL matlab_error( " Out of memory" )

      if(icntl(9).gt.0) then
         ! Take a copy of RHS for IR
         allocate(rhs(n,n_rhs), resid(n), stat=err)
         IF (err .ne. 0) &
            call matlab_error( " Out of memory" )
         rhs(:,:) = xxx(:,:)
      endif

      call cpu_time(start)
      print *, "try", icntl(1)
      call ma57cd(1, n, fact, size(fact), ifact, size(ifact), n_rhs, xxx, &
         size(xxx,1), work, size(work), iwork, icntl, info)
      call cpu_time(finish)
      rinfo(23) = finish-start

      if(icntl(9).gt.0) then
         ! Apply IR
         call cpu_time(start)
         call csr_to_coord(n, ip_in, irn_in, val_in, ne, irn, jcn, val)
         call ma57dd(2, n, ne, val, irn, jcn, fact, size(fact), ifact, &
            size(ifact), rhs, xxx, resid, work, iwork, icntl, cntl, info, rinfo)
         call cpu_time(finish)
         rinfo(24) = finish-start

         IF (info(1) == -8) &
            CALL matlab_warning("Iterative refinement did not converge in &
               &specified num. of iters")
         IF (info(1) == -14) &
            CALL matlab_warning( &
               "Failure in the estimates of the conditions numbers")
      endif

      ! Copy soln back to matlab
      plhs(x_out) = fortran_to_matlab(xxx)
      if(nlhs .ge. info_out) plhs(info_out) = fortran_to_matlab(info)
      if(nlhs .ge. rinfo_out) plhs(rinfo_out) = fortran_to_matlab(rinfo)

contains
   subroutine read_struct(mwptr, fact, ifact, cntl, icntl, info, rinfo)
      integer(mwp_) :: mwptr
      integer, dimension(:), allocatable :: ifact
      double precision, dimension(:), allocatable :: fact
      integer, dimension(:), allocatable :: icntl
      double precision, dimension(:), allocatable :: cntl
      integer, dimension(:), allocatable :: info
      double precision, dimension(:), allocatable :: rinfo

      integer(mwp_) :: fptr
      integer(mws_) :: mwlen

      if(.not.matlab_is_structure(mwptr)) &
         call mexErrMsgTxt('The third entry is not a structure')

      ! Read info and rinfo
      fptr = matlab_get_field_expert(mwptr, 1_mwi_, 'info ')
      call matlab_to_fortran(fptr, info, mwn, 'struct.info')
      fptr = matlab_get_field_expert(mwptr, 1_mwi_, 'rinfo')
      call matlab_to_fortran(fptr, rinfo, mwn, 'struct.rinfo')

      ! Check if MA57 completed sucessfully. If not, exit immediately
      ! If there are warnings, print them
      if(info(1).lt.0) return
      if(info(1).eq.4) &
          CALL matlab_warning("The matrix is not full rank: look at info(25)")

      ! If completion was sucessful, read factor and control data
      fptr = matlab_get_field_expert(mwptr, 1_mwi_, 'fact')
      call matlab_to_fortran(fptr, fact, mwlen, 'struct.fact')
      fptr = matlab_get_field_expert(mwptr, 1_mwi_, 'ifact')
      call matlab_to_fortran(fptr, ifact, mwlen, 'struct.ifact')
      fptr = matlab_get_field_expert(mwptr, 1_mwi_, 'cntl')
      call matlab_to_fortran(fptr, cntl, mwlen, 'struct.cntl')
      fptr = matlab_get_field_expert(mwptr, 1_mwi_, 'icntl')
      call matlab_to_fortran(fptr, icntl, mwlen, 'struct.icntl')
end subroutine read_struct

subroutine csr_to_coord(n, ip_in, irn_in, val_in, ne, irn, jcn, val)
   integer, intent(in) :: n
   integer, dimension(*), intent(in) :: ip_in, irn_in
   double precision, dimension(*), intent(in) :: val_in
   integer, intent(out) :: ne
   integer, dimension(:), allocatable, intent(out) :: irn, jcn
   double precision, dimension(:), allocatable, intent(out) :: val

   integer :: i, j, ne_in

   ! Allocate the space for the lower triangular part of the input matrix
   ne = ip_in(n+1)
   ALLOCATE( jcn(ne), irn(ne), val(ne), stat=err )
   IF (err .ne. 0) &
     CALL matlab_error( " Out of memory" )

   ne_in = 0
   DO i=1,n
     DO j=ip_in(i),ip_in(i+1)-1
        IF ((irn_in(j)) .ge. i) THEN
          ne_in = ne_in +1
          jcn(ne_in) = i
          irn(ne_in) = irn_in(j)
          val(ne_in) = val_in(j)
        END IF
     END DO
   END DO
   ne = ne_in
end subroutine csr_to_coord

END SUBROUTINE mexFunction

