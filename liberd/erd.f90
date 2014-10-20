module erd

  implicit none
  private

  public :: erd_batch

  public :: pack_coeffs ! for testing only

  integer, save                 :: INTSIZE = -1
  integer, save                 :: DBLSIZE = -1

  integer,          allocatable :: iwork(:) ! (INTSIZE)
  double precision, allocatable :: dwork(:) ! (DBLSIZE)

  contains

      subroutine erd_batch(la,    lb,    lc,    ld     &
                          , avect, bvect, cvect, dvect &
                          , aexps, bexps, cexps, dexps &
                          , acnts, bcnts, ccnts, dcnts &
                          , icont                      &
                          , batch, nints               &
                          )
        !
        !   Convenience wrapper for ERD__GENER_ERI_BATCH from the
        !   ERD integral package.
        !
        implicit none
   
        ! angular momenta:
        integer,          intent(in)  :: la, lb, lc, ld

        ! locations of the centers:   
        double precision, intent(in)  :: avect(3), bvect(3)
        double precision, intent(in)  :: cvect(3), dvect(3)
                                      
        ! exponents:                  
        double precision, intent(in)  :: aexps(:), bexps(:)
        double precision, intent(in)  :: cexps(:), dexps(:)
                                      
        ! contraction coeffs:         
        double precision, intent(in)  :: acnts(:,:), bcnts(:,:)
        double precision, intent(in)  :: ccnts(:,:), dcnts(:,:)

        ! tells what kind of contration to employ: 1 -- do, 0 -- dont:
        integer,          intent(in)  :: icont

        ! fake result:
        ! double precision, intent(out) :: res
        double precision, intent(out) :: batch(:) ! output integrals

        ! number of integrals filled in batch(*),
        ! zero if integrals vanish by selection rules:
        integer,          intent(out) :: nints
        ! *** end of interface ***

        ! total count of packed exponents and contraction coeffs:
        integer          :: nalpha, ncoeff

        ! total count of contractions:
        integer          :: ncsum

        ! shell indices as packed into linear arrays,
        ! so far quad(:) == [1,2,3,4] == [a,b,c,d]
        integer          :: quad(4)

        ! number of contractions and primitives IN THAT ORDER:
        integer          :: ncfps(4), npfps(4)

        ! linearly packed exponents and contraction coeffs:
        double precision :: alpha_pack(size(aexps)+size(bexps)+size(cexps)+size(dexps))

        ! for the uncontracted case (contractions == unity matrices) allocate more:
        double precision :: &
        coeff_pack(size(aexps)**2+size(bexps)**2+size(cexps)**2+size(dexps)**2)
!       double precision :: coeff_pack(size(acnts)+size(bcnts)+size(ccnts)+size(dcnts))

        ! starting and ending indices of exponents for segmented contractions,
        ! for generalized contractions assumed here these are alwyas 1 and
        ! the corresponding npfps(ishell). For the "uncontracted" case reserve
        ! more ranges than there are real contractions, one range per exponent:
        integer :: ccbeg_pack(size(aexps)+size(bexps)+size(cexps)+size(dexps))
        integer :: ccend_pack(size(aexps)+size(bexps)+size(cexps)+size(dexps))
!       integer :: ccbeg_pack(size(acnts,2)+size(bcnts,2)+size(ccnts,2)+size(dcnts,2))
!       integer :: ccend_pack(size(acnts,2)+size(bcnts,2)+size(ccnts,2)+size(dcnts,2))

        integer :: lvals(4)
        integer :: i
        integer :: imin, zmin
        integer :: iblk, zblk
        integer :: imax, zmax
        integer :: nfirst !, nints
        logical :: spherical = .true.

        ! for testing only, fill arrays with invalid entries:
        ccbeg_pack = -1
        ccend_pack = -1

        ! for convenience only:
        lvals(1) = la
        lvals(2) = lb
        lvals(3) = lc
        lvals(4) = ld

!       print *,'erd: call pack_coeffs()'
        call pack_coeffs( aexps, bexps, cexps, dexps    &
                        , acnts, bcnts, ccnts, dcnts    &
                        , quad, ncfps, npfps            &
                        , alpha_pack, nalpha            &
                        , coeff_pack, ncoeff            &
                        , ccbeg_pack, ccend_pack, ncsum &
                        , icont                         &
                        )

        ! maybe make pack_coeff() accept given order of shells:
        do i=1,4
          if( quad(i) /= i ) stop "FIXME: pack_coeffs reordered shells!"
        enddo

!       print *,'erd: call ERD__MEMORY_ERI_BATCH()'
        ! estimate the memory requirements for this integral batch:
        imax = 0
        zmax = 0
        call ERD__MEMORY_ERI_BATCH( nalpha, ncoeff                         &
                                  , ncfps(1), ncfps(2), ncfps(3), ncfps(4) &
                                  , npfps(1), npfps(2), npfps(3), npfps(4) &
                                  , lvals(1), lvals(2), lvals(3), lvals(4) &
                                  , avect(1), avect(2), avect(3)           &
                                  , bvect(1), bvect(2), bvect(3)           &
                                  , cvect(1), cvect(2), cvect(3)           &
                                  , dvect(1), dvect(2), dvect(3)           &
                                  , alpha_pack, coeff_pack                 &
                                  , spherical                              &
                                  , imin, iblk, zmin, zblk                 &
                                  )
        imax = max(imax, imin)
        zmax = max(zmax, zmin)

!       print *,'erd: ERD__MEMORY_ERI_BATCH returned imax=',imax,'zmax=',zmax

        ! realocate iwork, and dwork if necessary:
        if( imax > INTSIZE .or. zmax > DBLSIZE )then
          call reallocate_scratch(imax, zmax)
        endif

!       print *,'erd: call ERD__GENER_ERI_BATCH(',la, lb, lc, ld,')'

        ! compute the integral batch:
        call ERD__GENER_ERI_BATCH( INTSIZE, DBLSIZE                       &
                                 , nalpha, ncoeff, ncsum                  &
                                 , ncfps(1), ncfps(2), ncfps(3), ncfps(4) &
                                 , npfps(1), npfps(2), npfps(3), npfps(4) &
                                 , la, lb, lc, ld                         &
                                 , avect(1), avect(2), avect(3)           &
                                 , bvect(1), bvect(2), bvect(3)           &
                                 , cvect(1), cvect(2), cvect(3)           &
                                 , dvect(1), dvect(2), dvect(3)           &
                                 , alpha_pack, coeff_pack                 &
                                 , ccbeg_pack, ccend_pack                 &
                                 , spherical                              &
                                 , .true.                                 &
                                 , iwork                                  &
                                 , nints                                  &
                                 , nfirst                                 &
                                 , dwork                                  &
                                 )
!       print *,'erd: ERD__GENER_ERI_BATCH returned nints=',nints,'at nfirst=',nfirst

        ! FIXME: dont do this in production:
        if( nints == 0 )then
          print *,'erd: all integrals vanish'
          batch(:) = 0.0
        endif

        if( nints > size(batch) )then
          stop 'Error: integrals dont fit into output array!'
        endif

        ! copy result from dwork(nfirst:nfirst+nints-1):
        if( nints > 0 )then
!         call display( 2*lvals+1, ncfps, dwork(nfirst) )
          call copy(    2*lvals+1, ncfps, dwork(nfirst), batch )
        endif

      contains

        subroutine display( nm, ne, ints)
          implicit none
          integer, intent(in)           :: nm(4), ne(4)
          double precision, intent(in)  :: &
          ints( nm(1) * ne(1), nm(2) * ne(2), nm(3) * ne(3), nm(4) * ne(4) )
          ! *** end of interface ***

          integer :: i, j, k, l, ijkl

          ijkl = 0
          do l=1, nm(4) * ne(4)
          do k=1, nm(3) * ne(3)
          do j=1, nm(2) * ne(2)
          do i=1, nm(1) * ne(1)
            ijkl = ijkl + 1

            write(*,1000) ' [',i,j,'|',k,l,'] = ', ints(i,j,k,l)

          enddo
          enddo
          enddo
          enddo
          if( ijkl /= nints )then
            stop 'Error: integral count wrong!'
          endif
          1000 format (A2,2I3,A1,2I3,A4,F20.10)
        end subroutine display

        subroutine copy(nm, ne, ints, batch)
          !
          ! Piece of code for eathier addressing into batch(*)
          !
          implicit none
          integer, intent(in)   :: nm(4), ne(4)
          double precision, intent(in)  :: &
          ints(  nm(1)*ne(1), nm(2)*ne(2), nm(3)*ne(3), nm(4)*ne(4) )
          double precision, intent(out) :: &
          batch( nm(1)*ne(1), nm(2)*ne(2), nm(3)*ne(3), nm(4)*ne(4) )
          ! *** end of interface ***

          ! permutation of m-indices:
          integer :: ia(nm(1))
          integer :: ib(nm(2))
          integer :: ic(nm(3))
          integer :: id(nm(4))

          integer :: ls(4)
          integer :: ma,mb,mc,md
          integer :: i, j, k, l

          ! L-values:
          ls(:) = ( nm(:) - 1 ) / 2

          ! compute the permutations for magnetic shell indices:
          ia(:) = perm(ls(1))
          ib(:) = perm(ls(2))
          ic(:) = perm(ls(3))
          id(:) = perm(ls(4))

          ! unpack integrals in (1,2,3,4) (ang,exp) order:
          do md=1,nm(4)
          do mc=1,nm(3)
          do mb=1,nm(2)
          do ma=1,nm(1)

          do l=1,ne(4)
          do k=1,ne(3)
          do j=1,ne(2)
          do i=1,ne(1)

              batch( ia(ma) + (i - 1) * nm(1) &
                   , ib(mb) + (j - 1) * nm(2) &
                   , ic(mc) + (k - 1) * nm(3) &
                   , id(md) + (l - 1) * nm(4) &
                   ) =                        &
              batch( ia(ma) + (i - 1) * nm(1) &
                   , ib(mb) + (j - 1) * nm(2) &
                   , ic(mc) + (k - 1) * nm(3) &
                   , id(md) + (l - 1) * nm(4) &
                   )                          &
              + ints(   ma  + (i - 1) * nm(1) &
                    ,   mb  + (j - 1) * nm(2) &
                    ,   mc  + (k - 1) * nm(3) &
                    ,   md  + (l - 1) * nm(4) &
                    )
          enddo
          enddo
          enddo
          enddo

          enddo
          enddo
          enddo
          enddo
        end subroutine copy

        function perm(L) result(p)
          implicit none
          integer, intent(in) :: L
          integer             :: p(2*L+1)
          ! *** end of interface ***

          integer :: m

          select case ( L )

          case ( 0 ) ! s-shell
            p(1) = 1

          case ( 1 ) ! p-shell
            p(1) = 2 ! x goes to 2
            p(2) = 3 ! y goes to 3
            p(3) = 1 ! z goes to 1

!         case ( 2 ) ! d-shell
!           ! voodoo magic:
!           p(1) = 4 ! ? goes to 4
!           p(2) = 2 ! ? goes to 2
!           p(3) = 1 ! ? goes to 1
!           p(4) = 3 ! ? goes to 3
!           p(5) = 5 ! ? goes to 5

!         case ( 3 ) ! f-shell
!           ! voodoo magic:
!           p(1) = 6 ! ? goes to 6
!           p(2) = 4 ! ? goes to 4
!           p(3) = 2 ! ? goes to 2
!           p(4) = 1 ! ? goes to 1
!           p(5) = 3 ! ? goes to 3
!           p(6) = 5 ! ? goes to 5
!           p(7) = 7 ! ? goes to 7

!         case ( 4 ) ! g-shell
!           ! voodoo magic:
!           p(1) = 8 ! ? goes to 8
!           p(2) = 6 ! ? goes to 6
!           p(3) = 4 ! ? goes to 4
!           p(4) = 2 ! ? goes to 2
!           p(5) = 1 ! ? goes to 1
!           p(6) = 3 ! ? goes to 3
!           p(7) = 5 ! ? goes to 5
!           p(8) = 7 ! ? goes to 7
!           p(9) = 9 ! ? goes to 9

!         case ( 5 ) ! h-shell
!           ! voodoo magic:
!           p( 1) = 10 ! ? goes to 10
!           p( 2) =  8 ! ? goes to 8
!           p( 3) =  6 ! ? goes to 6
!           p( 4) =  4 ! ? goes to 4
!           p( 5) =  2 ! ? goes to 2
!           p( 6) =  1 ! ? goes to 1
!           p( 7) =  3 ! ? goes to 3
!           p( 8) =  5 ! ? goes to 5
!           p( 9) =  7 ! ? goes to 7
!           p(10) =  9 ! ? goes to 9
!           p(11) = 11 ! ? goes to 11

            ! there seems to be a rule emerging...
          case default
!           stop "Error: not yet implemented!"
            p(L + 1) = 1
            do m=1, L
              p(L + 1 - m) = 2 * m
              p(L + 1 + m) = 2 * m + 1
            enddo
          end select
        end function perm

      end subroutine erd_batch

      subroutine pack_coeffs( aexps, bexps, cexps, dexps    &
                            , acnts, bcnts, ccnts, dcnts    &
                            , quad, ncfps, npfps            &
                            , alpha_pack, nalpha            &
                            , coeff_pack, ncoeff            &
                            , ccbeg_pack, ccend_pack, ncsum &
                            , icont                         &
                            )
        !---------------------------------------------------------------------------
        !   Formats the integral exponents and contraction coefficients for use
        !   in the ERD integral package.
        !---------------------------------------------------------------------------
        implicit none
   
        ! exponents:
        double precision, intent(in) :: aexps(:), bexps(:)
        double precision, intent(in) :: cexps(:), dexps(:)

        ! contraction coeffs:
        double precision, intent(in) :: acnts(:,:), bcnts(:,:)
        double precision, intent(in) :: ccnts(:,:), dcnts(:,:)

        ! shell indices as packed into linear arrays,
        ! so far quad(:) == [1,2,3,4] == [a,b,c,d]
        integer, intent(out)          :: quad(4)

        ! number of contractions and primitives IN THAT ORDER:
        integer, intent(out)          :: ncfps(4), npfps(4)

        ! linearly packed exponents and contraction coeffs:
        double precision, intent(out) :: alpha_pack(:), coeff_pack(:)

        ! total count of packed exponents and contraction coeffs:
        integer,          intent(out) :: nalpha, ncoeff

        ! total count of contractions:
        integer,          intent(out) :: ncsum

        ! starting and ending indices of exponents for segmented contractions,
        ! for generalized contractions assumed here these are alwyas 1 and
        ! the corresponding npfps(ishell):
        integer,          intent(out) :: ccbeg_pack(:), ccend_pack(:)

        ! tells what kind of contration to employ: 1 -- do, 0 -- dont:
        integer,          intent(in)  :: icont
        ! *** end of interface ***

        integer i, j, k
        integer ishell

        ! shell indices:
        quad(1) = 1 ! a
        quad(2) = 2 ! b
        quad(3) = 3 ! c
        quad(4) = 4 ! d

        ! number of primitives:
        npfps(1) = size(aexps)
        npfps(2) = size(bexps)
        npfps(3) = size(cexps)
        npfps(4) = size(dexps)

        ! total exponent count:
        nalpha = sum(npfps)

        if( nalpha > size(alpha_pack) )then
           stop 'Error: alpha overflow in pack_coeffs'
        endif

        ! pack exponents into alpha_pack(:):
        k = 0
        do i=1,size(aexps)
          k = k + 1
          alpha_pack(k) = aexps(i)
        enddo
        do i=1,size(bexps)
          k = k + 1
          alpha_pack(k) = bexps(i)
        enddo
        do i=1,size(cexps)
          k = k + 1
          alpha_pack(k) = cexps(i)
        enddo
        do i=1,size(dexps)
          k = k + 1
          alpha_pack(k) = dexps(i)
        enddo
        if( k /= nalpha ) stop 'Error: total conunt of exponents wrong!'

        !
        ! pack contraction coefficients, as 
        !
        ! generalized contractions, (icont == 1)
        ! segmented contraction or  (icont == 2)
        ! simulated un-contractions (icont == 0)
        !
        select case ( icont )

        case ( 1 ) ! generalized contractions:
          !
          ! pack all contraction coeffs into coeff_pack(:),
          ! even if come of them may be treated as segmented:
          !

          ! number of (generalized) contractions:
          ncfps(1) = size(acnts,2)
          ncfps(2) = size(bcnts,2)
          ncfps(3) = size(ccnts,2)
          ncfps(4) = size(dcnts,2)

          ! total count of contraction coefficients:
          ncoeff = sum(npfps * ncfps)

          if( ncoeff > size(coeff_pack) )then
            stop 'Error: coeff overflow in pack_coeffs, case 1'
          endif

          ! total count of contractions:
          ncsum = sum(ncfps)

          if( ncsum > size(ccbeg_pack) .or. ncsum > size(ccend_pack) )then
            stop 'Error: ccbeg_pack overflow in pack_coeffs, case 1'
          endif

          k = 0
          do j=1,size(acnts,2)
            do i=1,size(acnts,1)
              k = k + 1
              coeff_pack(k) = acnts(i,j)
            enddo
          enddo
          do j=1,size(bcnts,2)
            do i=1,size(bcnts,1)
              k = k + 1
              coeff_pack(k) = bcnts(i,j)
            enddo
          enddo
          do j=1,size(ccnts,2)
            do i=1,size(ccnts,1)
              k = k + 1
              coeff_pack(k) = ccnts(i,j)
            enddo
          enddo
          do j=1,size(dcnts,2)
            do i=1,size(dcnts,1)
              k = k + 1
              coeff_pack(k) = dcnts(i,j)
            enddo
          enddo
          if( k /= ncoeff ) stop 'Error: total count of contraction coeffs wrong!'

          ! fill in contraction ranges (CCBEG:CCEND) for all
          ! contracted functions of all 4 shells:
          k = 0
          do ishell = 1, 4
             do j = 1, ncfps(ishell)
                k = k + 1
                ! for generalized contraction all ranges are the same:
                ccbeg_pack(k) = 1
                ccend_pack(k) = npfps(ishell)
                ! FIXME: for segmented contractions it may be different:
             enddo
          enddo
          if( k /= ncsum ) stop 'Error: total count of contractions wrong!'

        case ( 0 ) ! uncontracted:
          !
          ! pack fake unity contraction matrices into coeff_pack(:),
          ! first redefine necessary counts:
          !

          ! number of "contractions" == number of exponents:
          ncfps(1:4) = npfps(1:4)

          ! total number of contractions:
          ncsum = sum(ncfps)

          if( ncsum > size(ccbeg_pack) .or. ncsum > size(ccend_pack) )then
            stop 'Error: ccbeg_pack overflow in pack_coeffs, case 0'
          endif

          ! store all of unity matrices, sum of squares::
          ncoeff = sum( ncfps * npfps )

          if( ncoeff > size(coeff_pack) )then
            stop 'Error: coeff overflow in pack_coeffs, case 0'
          endif

          ! store "unity matrices":
          k = 0
          do ishell=1,4
            do j=1,ncfps(ishell)
              do i=1,npfps(ishell)
                k = k + 1
                ! delta(i,j):
                if( i /= j )then
                  coeff_pack(k) = 0.0D0
                else
                  coeff_pack(k) = 1.0D0
                endif
              enddo
            enddo
          enddo
          if( k /= ncoeff ) stop 'Error: total count of contraction coeffs wrong!'

          ! fill in contraction ranges (CCBEG:CCEND) for all
          ! contracted functions of all 4 shells:
          k = 0
          do ishell = 1, 4
             do j = 1, ncfps(ishell)
                k = k + 1
                ! for uncontracted functions each range is exactly one num:
                ccbeg_pack(k) = j ! not 1
                ccend_pack(k) = j ! not npfps(ishell)
             enddo
          enddo
          if( k /= ncsum ) stop 'Error: total count of contractions wrong!'

        case default
          stop "no such case"
        end select

!       write(*,1001) "quad", quad
!       write(*,1001) "npfps", npfps
!       write(*,1001) "ncfps", ncfps
!       write(*,1001) "ncsum", ncsum
!       write(*,1001) "ccbeg", ccbeg_pack
!       write(*,1001) "ccend", ccend_pack
!       write(*,1000) "coeff", coeff_pack

!       1000 format(A8,5F16.8,/:(8X,5F16.8))
!       1001 format(A8,10I3,/:(8X,10I3))
      end subroutine pack_coeffs

      subroutine reallocate_scratch(imax, zmax)
        !
        ! Reallocate iwork and dwork if necessary.
        !
        implicit none
        integer, intent(in) :: imax, zmax
        ! *** end of interface ***

        integer :: memstat

        if( imax > INTSIZE )then

          print *,'erd: reallocate_scratch resizing iwork from',INTSIZE,'to', imax

          if( allocated(iwork) )then
            deallocate(iwork, stat=memstat)
            if( memstat /= 0 ) stop "dealloc iwork failed"
          endif

          allocate(iwork(imax), stat=memstat)
          if( memstat /= 0 ) stop "alloc iwork failed"

          INTSIZE = imax
        endif

        if( zmax > DBLSIZE )then

          print *,'erd: reallocate_scratch resizing dwork from',DBLSIZE,'to', zmax

          if( allocated(dwork) )then
            deallocate(dwork, stat=memstat)
            if( memstat /= 0 ) stop "dealloc dwork failed"
          endif

          allocate(dwork(zmax), stat=memstat)
          if( memstat /= 0 ) stop "alloc dwork failed"

          DBLSIZE = zmax
        endif
      end subroutine reallocate_scratch

end module erd
