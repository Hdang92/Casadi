!#define A_IN       (prhs(1))
!#define CONTROL    (prhs(2))

!#define REV_OUT    (plhs(1))
!#define INFO_OUT   (plhs(2))
!#define RINFO_OUT  (plhs(3))



! The gateway subroutine
SUBROUTINE mexFunction( nlhs, plhs, nrhs, prhs )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                              STATEMENT                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use hsl_matlab

      IMPLICIT NONE

! External mex-functions

      INTEGER * 4 :: nlhs, nrhs
      integer(mwp_) :: plhs( * ), prhs( * )

      integer(mwp_) :: fptr
      integer(mwi_) :: ind_i
      LOGICAL :: isset

! pointer for input/output parameters and local variables

      INTEGER, DIMENSION(:), ALLOCATABLE :: ip_in,irn_in,jcn,irn
      INTEGER, DIMENSION(:), ALLOCATABLE :: keep,iwork,ifact,newifc
      INTEGER, DIMENSION(:), ALLOCATABLE :: icntl,info

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fact,val,val_in,newfac
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: cntl,rinfo

      INTEGER, PARAMETER :: lcntl = 5,licntl = 20,linfo = 40,lrinfo = 24
      integer(mws_) :: mwm, mwn
      INTEGER :: i,j,n,loopf,loopc,order
      INTEGER :: liwork,lkeep,lfact,lifact
      INTEGER :: ne_in,ne,ic,lnew,linew
      INTEGER :: err

      DOUBLE PRECISION :: start, finish,count_time_a,count_time_f

      CHARACTER(len=132) msg
      CHARACTER (50) value
      CHARACTER(LEN=6) f1,f2,f3,f4,f5,f6,f7,f8
      CHARACTER(LEN=6), DIMENSION(8):: fieldnames
      CHARACTER(LEN=6) c1,c2,c3
      CHARACTER(LEN=6), DIMENSION(3):: controlnames

      INTRINSIC cpu_time
! The fieldsnames holds the Structure fields for passing the output
! of ma57_factor to the solve
      f1 = 'lfact '
      f2 = 'fact  '
      f3 = 'lifact'
      f4 = 'ifact '
      f5 = 'cntl  '
      f6 = 'icntl '
      f7 = 'info  '
      f8 = 'rinfo '
      fieldnames = (/f1,f2,f3,f4,f5,f6,f7,f8/)
! The controlnames hols the CONTROL structure fields in input
      c1 = 'order '
      c2 = 'thres '
      c3 = 'stpv  '
      controlnames = (/c1,c2,c3/)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                               CODE                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  CNTL and ICNTL defaults

      ALLOCATE( cntl( lcntl ), STAT = err )
      IF (err .ne. 0) THEN
        CALL mexErrMsgTxt( " Out of memory " )
      END IF
      ALLOCATE( icntl( licntl ), STAT = err )
      IF (err .ne. 0) THEN
        CALL mexErrMsgTxt( " Out of memory " )
      END IF
      ALLOCATE( info( linfo ), STAT = err )
      IF (err .ne. 0) THEN
        CALL mexErrMsgTxt( " Out of memory " )
      END IF

      CALL ma57id(cntl,icntl)
      icntl(1) = 0
      icntl(2) = 0
      icntl(3) = 0
      icntl(5) = -1
      icntl(6) = 4
      icntl(8) = 1

! Check of input arguments
! MATLAB_get_ptr(pm) determines the starting address of the
! real data in the mxArray that pm points to
!
!!
!!

      IF (nrhs > 2 .or. nrhs < 1 .or. nlhs < 1 .or. nlhs > 3) THEN
        CALL mexErrMsgTxt( " Wrong # of arguments " )
      END IF

! [m,n] = size(Matrix used)

      mwn = MATLAB_get_n(prhs(1))
      mwm = MATLAB_get_m(prhs(1))

!check if A sparse and square

      IF ((mwn .ne. mwm) .or. (.not.MATLAB_is_sparse(prhs(1)))) THEN
        CALL mexErrMsgTxt("The input matrix A must be square and sparse")
      END IF
!read CONTROL structure or set the defaults
      IF ( nrhs == 1 ) THEN
        order = 0
        cntl(1) = 0.01D+0
        cntl(4) = 0.0D+0
      END IF

      IF ( nrhs == 2 ) THEN
! TEST STRUCTURE
        IF (.not. MATLAB_is_structure(prhs(2))) THEN
          call mexErrMsgTxt('The second entry is not a structure')
        END IF
!Read control structure
        call matlab_get_field(prhs(2),c1,order,'control',isset,0)
        if (.not. isset) call MATLAB_warning( &
          'control.order is not assigned and it is set to default value 0')

        call matlab_get_field(prhs(2),c2,cntl(1),'control',isset, 0.01d0)
        if (.not. isset) call MATLAB_warning( &
           "control.thres is not assigned and it is set to default value 0.01")

        call matlab_get_field(prhs(2),c3,cntl(4),'control',isset, 0.0d0)
        if (.not. isset) call MATLAB_warning( &
          "control.stpv is not assigned and it is set to default value 0.0")

      END IF

      IF ( order <= 0 ) THEN
        icntl(6) = 2
      END IF
! this section will become active when the new nested dissection will be available
!      IF ( order == 1 ) THEN
!        icntl(6) = 4
!      END IF

! read A
      call matlab_to_fortran(prhs(1),mwm,mwn,ip_in,irn_in,val_in,"A")
      n = mwn

! ne = number of nonzero
      ne = ip_in(n+1)

!
! Allocate the space for the lower triangular part of the input matrix
!
      ALLOCATE( jcn( ne ), STAT = err )
      IF (err .ne. 0) THEN
        CALL mexErrMsgTxt( " Out of memory " )
      END IF
      ALLOCATE( irn( ne ), STAT = err )
      IF (err .ne. 0) THEN
        CALL mexErrMsgTxt( " Out of memory " )
      END IF
      ALLOCATE( val( ne ), STAT = err )
      IF (err .ne. 0) THEN
        CALL mexErrMsgTxt( " Out of memory " )
      END IF

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


      IF ( ALLOCATED( ip_in  ) ) DEALLOCATE( ip_in,  STAT = err)
      IF ( ALLOCATED( irn_in ) ) DEALLOCATE( irn_in, STAT = err)
      IF ( ALLOCATED( val_in ) ) DEALLOCATE( val_in, STAT = err)


! Fix others input parameters for analysis

      IF ( cntl(4) > 0 ) THEN
        cntl(5) = 0.0D+0
      END IF

      lkeep = 7*n+2*ne+42

      ALLOCATE( keep( lkeep ), STAT = err )
      IF (err .ne. 0) THEN
        CALL mexErrMsgTxt( " Out of memory " )
      END IF
      DO i=1,lkeep
        keep(i) = 0
      END DO
      IF ( order > 0 ) THEN
        icntl(6) = 1
        ! fix the input permutation
        DO i=1,n
          keep(i) = i
        END DO
      END IF

      liwork = 5*n
      ALLOCATE( iwork( liwork ), STAT = err )
      IF (err .ne. 0) THEN
        CALL mexErrMsgTxt( " Out of memory " )
      END IF
      ALLOCATE( rinfo( lrinfo ), STAT = err )
      IF (err .ne. 0) THEN
        CALL mexErrMsgTxt( " Out of memory " )
      END IF

      DO i=1,liwork
        iwork(i) = 0
      END DO
      DO i=1,linfo
        info(i) = 0
      END DO
      DO i=1,lrinfo
        rinfo(i) = 0
      END DO

! Do analysis

      count_time_a = 0

      CALL cpu_time (start)
      CALL ma57ad(n,ne,irn,jcn,lkeep,keep,iwork,icntl,info,rinfo)
      CALL cpu_time (finish)
      count_time_a = count_time_a + (finish - start)
      rinfo(21) = count_time_a


! ERROR MSG when MeTis is not present

      IF ( info(1) == -18 ) THEN
        CALL mexErrMsgTxt( &
          'The Library METIS does not exist but user set order = 1')
      END IF

! fix others input parameters for factorization

      lfact = 2*info(9)
!      lfact = 1*info(9)
      ALLOCATE( fact( lfact ), STAT = err )
      IF (err .ne. 0) THEN
        CALL mexErrMsgTxt( " Out of memory : malloc failed to allocate  &
        &the real array\n" )
      END IF

      lifact = 2*info(10)
!      lifact = 1*info(10)
      ALLOCATE( ifact( lifact ), STAT = err )
      IF (err .ne. 0) THEN
        CALL mexErrMsgTxt( " Out of memory : malloc failed to allocate  &
        &the integer array\n" )
      END IF

      count_time_f = 0
      loopf = -1
      loopc = 0

      DO WHILE (loopf < 0 )


        !  do factorization
        CALL cpu_time (start)
        CALL ma57bd(n,ne,val,fact,lfact,ifact,lifact,lkeep,keep,iwork,icntl, &
          cntl,info,rinfo)
        CALL cpu_time (finish)
        count_time_f = count_time_f + (finish - start)
        rinfo(22) = count_time_f

        loopc = loopc + 1

        IF (info(1) == 4) THEN
!          CALL mexPrintf('\n\n')
          CALL mexWarnMsgTxt("The matrix is not full rank: look at info(25)")
!           CALL mexWarnMsgTxt("MA57BD will produce a factorization &
!           &where D is singular \n\n")
        END IF

        IF (info(1) == 0 .OR. info(1) == 4) THEN
          loopf = 0
        END IF
        IF (info(1) < 0) THEN
          loopf = -1
        END IF

        IF (info(1) == -3) THEN
          lfact = info(17)
          IF ( ALLOCATED( fact ) ) DEALLOCATE( fact, STAT = err)
          ALLOCATE(fact(lfact), STAT = err)
          IF (err .ne. 0) THEN
            CALL mexErrMsgTxt( " Out of memory : malloc failed to allocate  &
            &the real array\n" )
          END IF
          loopf = -1
        END IF
        IF (info(1) == -4) THEN
          lifact = info(18)
          IF ( ALLOCATED( ifact ) ) DEALLOCATE( ifact, STAT = err)
          ALLOCATE(ifact(lifact), STAT = err)
          IF (err .ne. 0) THEN
            CALL mexErrMsgTxt( " Out of memory : malloc failed to allocate  &
            &the integer array\n" )
          END IF
          loopf = -1
        END IF

        IF (info(1) == 10) THEN
          ic = 0
          lnew = 2*lfact
          ALLOCATE( newfac( lnew ), STAT = err )
          IF (err .ne. 0) THEN
            CALL mexErrMsgTxt( " Out of memory " )
          END IF

          CALL ma57ed(n,ic,keep,fact,lfact,newfac,lnew,ifact,lifact,ifact, &
            lifact,info)

          IF ( ALLOCATED( fact ) ) DEALLOCATE( fact, STAT = err)
          ALLOCATE( fact( lnew ), STAT = err )
          IF (err .ne. 0) THEN
            CALL mexErrMsgTxt( " Out of memory " )
          END IF


          DO i =1,lnew
            fact(i) = newfac(i)
          END DO

          lfact = lnew
          IF ( ALLOCATED( newfac ) ) DEALLOCATE( newfac, STAT = err)

        END IF

        IF (info(1) == 11) THEN
          ic = 1
          linew = 2*lifact

          ALLOCATE( newifc( linew ), STAT = err )
          IF (err .ne. 0) THEN
            CALL mexErrMsgTxt( " Out of memory " )
          END IF

          CALL ma57ed(n,ic,keep,fact,lfact,fact,lfact,ifact,lifact,newifc, &
            linew,info)

          IF ( ALLOCATED( ifact ) ) DEALLOCATE( ifact, STAT = err)
          ALLOCATE(ifact(linew), STAT = err )
          IF (err .ne. 0) THEN
            CALL mexErrMsgTxt( " Out of memory " )
          END IF

          DO i=1,linew
            ifact(i) = newifc(i)
          END DO
          lifact = linew
          IF ( ALLOCATED( newifc ) ) DEALLOCATE( newifc, STAT = err)
        END IF

      END DO

! Free workspace

      IF ( ALLOCATED( val ) ) DEALLOCATE( val, STAT = err)
      IF ( ALLOCATED( jcn ) ) DEALLOCATE( jcn, STAT = err)
      IF ( ALLOCATED( irn ) ) DEALLOCATE( irn, STAT = err)
      IF ( ALLOCATED( keep ) ) DEALLOCATE( keep, STAT = err)
      IF ( ALLOCATED( iwork ) ) DEALLOCATE( iwork, STAT = err)


! Fix output arguments

! INFO_OUT = plhs(5)
      plhs(2) =  fortran_to_matlab(info)
! RINFO_OUT = plhs(6)
      plhs(3) = fortran_to_matlab(rinfo)
! REV_OUT = plhs(1)

      plhs(1) = matlab_create_structure(fieldnames)

      ind_i = 1
      call matlab_set_field(plhs(1),f1,lfact)
!       call matlab_set_field(plhs(1),f2,fact)
      call matlab_set_field(plhs(1),f3,lifact)
!       call matlab_set_field(plhs(1),f4,ifact)
!       call matlab_set_field(plhs(1),f5,cntl)
!       call matlab_set_field(plhs(1),f6,icntl)
!       call matlab_set_field(plhs(1),f7,info)
!       call matlab_set_field(plhs(1),f8,rinfo)

      fptr = fortran_to_matlab(fact)
      call mxSetField( plhs(1), ind_i, f2, fptr )

      fptr = fortran_to_matlab(ifact)
      call mxSetField( plhs(1), ind_i, f4, fptr )

      fptr = fortran_to_matlab(cntl)
      call mxSetField( plhs(1), ind_i, f5, fptr )

      fptr = fortran_to_matlab(icntl)
      call mxSetField( plhs(1), ind_i, f6, fptr )

      fptr = fortran_to_matlab(info)
      call mxSetField( plhs(1), ind_i, f7, fptr )

      fptr = fortran_to_matlab(rinfo)
      call mxSetField( plhs(1), ind_i, f8, fptr )
! Deallocations

      IF ( ALLOCATED( fact ) ) DEALLOCATE( fact, STAT = err)
      IF ( ALLOCATED( ifact ) ) DEALLOCATE( ifact, STAT = err)
      IF ( ALLOCATED( cntl ) ) DEALLOCATE( cntl, STAT = err)
      IF ( ALLOCATED( icntl ) ) DEALLOCATE( icntl, STAT = err)
      IF ( ALLOCATED( info ) ) DEALLOCATE( info, STAT = err)
      IF ( ALLOCATED( rinfo ) ) DEALLOCATE( rinfo, STAT = err)

END SUBROUTINE mexFunction
