!
! ParaGauss, a program package for high-performance computations
! of molecular systems
! Copyright (C) 2014
! T. Belling, T. Grauschopf, S. Krüger, F. Nörtemann, M. Staufer,
! M. Mayer, V. A. Nasluzov, U. Birkenheuer, A. Hu, A. V. Matveev,
! A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman, D. I. Ganyushin,
! T. Kerdcharoen, A. Woiterski, A. B. Gordienko, S. Majumder,
! M. H. i Rotllant, R. Ramakrishnan, G. Dixit, A. Nikodem, T. Soini,
! M. Roderus, N. Rösch
!
! This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License version 2 as published
! by the Free Software Foundation [1].
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! [1] http://www.gnu.org/licenses/gpl-2.0.html
!
! Please see the accompanying LICENSE file for further information.
!
!===============================================================
! Public interface of module
!===============================================================
module relgrads_store
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !    the integral pices are mapped onto the ``hash'' values
  !    that are also used for filenames.
  !    Currently hash consists of 9 digits
  !             SSSXXXYYY
  !    each three stayuing for
  !     - symmetry irrep       (SSS)
  !     - 1st derivative index (XXX)
  !     - 2nd derivative index (YYY)
  !
  !  Reads/writes to disk (on low memory):
  !
  !     TMP_DIR/PPPSSSXXXYYY.swp
  !
  !  in addition to the hash the file names have:
  !
  !     - three digit processor index prefix (PPP)
  !     - file extension suffix              (.swp)
  !
  ! Copyright (c) Alexei Matveev
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------

!define FPP_TIMERS 1
# include "def.h"
  use type_module, only: &
       IK=>i4_kind, &
       RK=>r8_kind ! type specification parameters
  use error_module, only: MyID
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  integer(IK), parameter, public :: &
         RGNUCL =  1 &
       , RGPVSP =  2 &
       , RGKNTC =  3 &
       , RGOVRL =  4 &
       , RGPSEU =  5 &
       , RGOTHR =  6

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  public :: rg_open
  public :: rg_close
  public :: rg_put
  public :: rg_gr_put
  public :: rg_sd_put
  public :: rg_get
  public :: rg_gr_get
  public :: rg_sd_get

  public :: rg_bcast
  public :: rg_distr

  public :: rg_xx_get
  public :: rg_flag

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  integer(IK), parameter, private :: MAXHASH=9 ! 'SSSXXXYYY'
  character(len=MAXHASH), private :: ANYHASH=    '?????????'

  type, private :: inode
    type(block), pointer :: head => NULL()

    ! for statistics:
    integer(IK)          :: blk_count = 0
    integer(IK)          :: size      = 0 ! sum of all block-sizes

    ! for file backend:
    character(len=MAXHASH):: cur_hash = 'NULL'
    ! %cur_hash is set to non-NULL ``hash'' only after swapInOne(ino,hash)

    ! list of hashes/files that have already been created on disk:
    type(swap_t), pointer :: swp_list => NULL()
  end type

  integer(IK), parameter, private :: ITY=1 &
                                   , ISY=2 &
                                   , IUA=3 &
                                   , ILA=4 &
                                   , IUB=5 &
                                   , ILB=6 &
                                   , IXI=7 &
                                   , IXJ=8
  integer(IK), parameter, private :: MAXMETA=8

  ! to mark blocks set/clear these bits in %bits:
  integer(IK), parameter, private :: IDIRTY = 2**0 &
                                   , IMARK  = 2**1

  type, private :: block
    real(RK), allocatable :: sto(:)
    integer(IK)          :: meta(MAXMETA) = 0 ! typ,ua,la,ub,lb,i,j
    type(block), pointer :: prev   => NULL()
    type(block), pointer :: next   => NULL()
    integer(IK)          :: bits   =  IDIRTY  ! create new blocks as DIRTY
    ! for debug:
    integer(IK)          :: blk_number = 0
  end type

  type, private :: swap_t
    character(len=MAXHASH) :: hash = 'NULL' ! part of the filename
    ! hash will be derived from metadata such like derivative indices (x,y) and symmetry.
    type(swap_t), pointer  :: next => NULL()
  end type

  !------------ Declaration of constants and variables ----

  integer(IK), parameter   :: MAXORD =  2

  ! virtual filesystem:
  type(inode), private :: FS(0:MAXORD)

  integer(IK) :: mem=0
  integer(IK) :: mem0=0
  integer(IK) :: mem_swp=0

  ! after reaching this value start swaping:
  integer(IK), parameter :: MAXMEM=200 *1024**2/8 ! ... Megabytes
! integer(IK), parameter :: MAXMEM=HUGE(1)/1024
  real(RK)   , parameter :: MINLEV=0.6           ! reduce inode memory to this level in swapOut()

  integer(IK) :: flag(RGOTHR) = 0

  FPP_TIMER_DECL(tot)
  FPP_TIMER_DECL(disk)
  FPP_TIMER_DECL(bcst)
  FPP_TIMER_DECL(dist)

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  !*************************************************************
  subroutine rg_put(typ,irr,ua,la,ub,lb,ints)
    implicit none
    integer(IK), intent(in) :: typ
    integer(IK), intent(in) :: irr
    integer(IK), intent(in) :: ua,la
    integer(IK), intent(in) :: ub,lb
    real(RK)   , intent(in) :: ints(:,:,:,:) ! (NEA,NEB,NIA,NIB)
    ! *** end of interface ***

    call rg_xx_put(0,meta(typ,irr,ua,la,ub,lb,0,0),ints)
  end subroutine rg_put
  !*************************************************************

  !*************************************************************
  subroutine rg_gr_put(typ,irr,ua,la,ub,lb,i,ints)
    implicit none
    integer(IK), intent(in) :: typ
    integer(IK), intent(in) :: irr
    integer(IK), intent(in) :: ua,la
    integer(IK), intent(in) :: ub,lb
    integer(IK), intent(in) :: i
    real(RK)   , intent(in) :: ints(:,:,:,:) ! (NEA,NEB,NIA,NIB)
    ! *** end of interface ***

    call rg_xx_put(1,meta(typ,irr,ua,la,ub,lb,i,0),ints)
  end subroutine rg_gr_put
  !*************************************************************

  !*************************************************************
  subroutine rg_sd_put(typ,irr,ua,la,ub,lb,i,j,ints)
    implicit none
    integer(IK), intent(in) :: typ
    integer(IK), intent(in) :: irr
    integer(IK), intent(in) :: ua,la
    integer(IK), intent(in) :: ub,lb
    integer(IK), intent(in) :: i,j
    real(RK)   , intent(in) :: ints(:,:,:,:) ! (NEA,NEB,NIA,NIB)
    ! *** end of interface ***

    call rg_xx_put(2,meta(typ,irr,ua,la,ub,lb,i,j),ints)
  end subroutine rg_sd_put
  !*************************************************************

  !*************************************************************
  subroutine rg_get(irr,S,T,V,O)
    implicit none
    integer(IK), intent(in)  :: irr
    real(RK)   , intent(out) :: S(:,:) ! (N,N)
    real(RK)   , intent(out) :: T(:,:) ! (N,N)
    real(RK)   , intent(out) :: V(:,:) ! (N,N)
    real(RK)   , intent(out) :: O(:,:) ! (N,N)
    ! *** end of interface ***

    call rg_xx_get(0,irr,0,0,S,T,V,O)
  end subroutine rg_get
  !*************************************************************

  !*************************************************************
  subroutine rg_gr_get(irr,i,S,T,V,O)
    implicit none
    integer(IK), intent(in)  :: irr
    integer(IK), intent(in)  :: i
    real(RK)   , intent(out) :: S(:,:) ! (N,N)
    real(RK)   , intent(out) :: T(:,:) ! (N,N)
    real(RK)   , intent(out) :: V(:,:) ! (N,N)
    real(RK)   , intent(out) :: O(:,:) ! (N,N)
    ! *** end of interface ***

    call rg_xx_get(1,irr,i,0,S,T,V,O)
  end subroutine rg_gr_get
  !*************************************************************

  !*************************************************************
  subroutine rg_sd_get(irr,i,j,S,T,V,O)
    implicit none
    integer(IK), intent(in)  :: irr
    integer(IK), intent(in)  :: i,j
    real(RK)   , intent(out) :: S(:,:) ! (N,N)
    real(RK)   , intent(out) :: T(:,:) ! (N,N)
    real(RK)   , intent(out) :: V(:,:) ! (N,N)
    real(RK)   , intent(out) :: O(:,:) ! (N,N)
    ! *** end of interface ***

    call rg_xx_get(2,irr,i,j,S,T,V,O)
  end subroutine rg_sd_get
  !*************************************************************

  function meta(typ,sym,ua,la,ub,lb,i,j) result(m)
    implicit none
    integer(IK), intent(in) :: typ,sym,ua,la,ub,lb,i,j
    integer(IK)             :: m(MAXMETA)
    ! *** end of interface ***

    m(ITY) = typ
    m(ISY) = sym
    m(IUA) = ua
    m(ILA) = la
    m(IUB) = ub
    m(ILB) = lb
    m(IXI) = i
    m(IXJ) = j
  end function meta

  !*************************************************************
  subroutine rg_xx_put(der,met,ints)
    implicit none
    integer(IK), intent(in) :: der
    integer(IK), intent(in) :: met(:)
    real(RK)   , intent(in) :: ints(:,:,:,:) ! (NEA,NEB,NIA,NIB)
    ! *** end of interface ***

    if(size(ints)==0) RETURN ! dont write empty records!
    FPP_TIMER_START(tot)

    call append(met,ints,FS(der))

    call oom()

    ! mark stored integral types:
    flag(met(ITY)) = flag(met(ITY)) + 1
    FPP_TIMER_STOP(tot)
  end subroutine rg_xx_put

  subroutine oom()
    ! checks/prints memory usage
    ! calls out of memory hadler
    implicit none
    ! *** end of interface ***

    integer(IK) :: x(1),i
    integer(IK) :: lowmem

#ifdef FPP_DEBUG
    ! to print stats:
    call oom_stat()
#endif

    ! check if it is time for swapping:
    lowmem = mem - MAXMEM

    ! FIXME: do while ( lowmem > 0 ) ???
    if( lowmem > 0 )then
       ! swap to disk the largest inode.
       ! for that first find a memory hog:
       x = maxloc( FS(:)%size )

       ! FS allocated as FS(0:MAXDER):
       i = x(1) - 1

       ! try to free some memory from the largest inode (i==2?):
       if( i /= 0 ) call oom_inode(FS(i))
       ! FIXME: integrals (i=0) are bcast`ed
       !        with mask=ANYAHSH, which fails if FS(0) is partly in swap
    endif
  end subroutine oom

#ifdef FPP_DEBUG
  subroutine oom_stat()
    !
    ! Prints memory usage
    !
    implicit none
    ! *** end of interface ***

    integer(IK) :: lowmem
    character(len=5) :: ok
    real, parameter :: M = 1024**2/8.0

    ! check if it is time for swapping:
    lowmem = mem - MAXMEM

    ok = "ok..."
    if( lowmem > 0 ) ok = "OOM!!"

    ! print usage every .1M:
    if( abs(mem - mem0) > 10*M )then
       mem0 = mem
       write(*,'(X,A4,"RGS: MEM phys=",F7.1,"M swap=",F7.1,"M ino=",3F7.1,A6)') &
            MyID         &
            , mem/M        &
            , mem_swp/M    &
            , FS(:)%size/M &
            , ok
    endif
  end subroutine oom_stat
#endif

  subroutine oom_inode(ino)
    ! out of memory hadler
    implicit none
    type(inode), intent(inout) :: ino
    ! *** end of interface ***

    integer(IK) :: siz

    ! On exit the inode should have the size scaled down by MINLEV:
    siz = int (ino % size * MINLEV, kind = IK)

    !print *,MyID,'RGS oom_inode: ino%cur_hash=>'//ino%cur_hash//'<'
    if( ino%cur_hash /= 'NULL' )then
      ! some memory is already in swap, one may just drop some!
      call swapOutOne(ino,ino%cur_hash)
      ! (will also WRITE the dirty blocks with the same hash)

      ! FIXME: maybe swapping cur_hash out is not optimal?
      !        maybe one should track the least-recently-used ones?

      ino%cur_hash = 'NULL'
    endif

    if( ino%size <= siz ) RETURN

    ! continue swapping out dirty memory until the size reduces to siz:
    call swapOut(ino,siz)
  end subroutine oom_inode
  !*************************************************************

  !*************************************************************
  subroutine append(met,dat,ino)
    USE_MEMLOG
    implicit none
    integer(IK), intent(in) :: met(:)        ! metadata
    real(RK)   , intent(in) :: dat(:,:,:,:) ! (NEA,NEB,NIA,NIB)
    type(inode), intent(inout) :: ino
    ! *** end of interface ***

    type(block), pointer :: blk

    ! create new block:
    allocate(blk)

    ! copy meta-data:
    blk%meta = met

    ! allocate storage in a new block:
    allocate(blk%sto(size(dat)))

    ! update global memory counters:
    MEMLOG(+size_of (blk))
    mem = mem + size_of (blk)

    ! append new block to the list, update inode size:
    call insblk(ino,blk)

    ! unwindes the legacy index order in ``dat'':
    call cp_data_rec(dat,blk%sto)
  end subroutine append

  subroutine insblk(ino, blk)
    ! appends blk to ino%tail, updates inode size
    implicit none
    type(inode), intent(inout) :: ino
    type(block), pointer       :: blk
    ! *** end of interface **

ASSERT(associated(blk))
ASSERT(.not.associated(blk%prev))
ASSERT(.not.associated(blk%next))
ASSERT(allocated(blk%sto))

    ! attach current block list to the CDR field of the block:
    blk%next => ino%head

    ! block with old list in CDR is the new list:
    ino%head => blk

    ! for a dbouble linked list we need to update the CAR field of the
    ! old list, if it was not empty:
    if ( associated(blk%next) ) then
ASSERT(.not.associated(blk%next%prev))
       blk%next%prev => blk
    endif

    ! update inode size:
    ino%size = ino%size + size_of (blk)

    ! for debug:
    ino%blk_count  = ino%blk_count + 1
    blk%blk_number = ino%blk_count
  end subroutine insblk

  subroutine delblk(ino, blk)
    USE_MEMLOG
    implicit none
    type(inode), intent(inout) :: ino
    type(block), pointer       :: blk ! returns blk => blk%next
    ! *** end of interface **

    type(block), pointer :: prev
    type(block), pointer :: next ! return value

    ASSERT(associated(blk))

    ! save %prev, %next before deallocation:
    prev => blk%prev
    next => blk%next

    ! connect previous (-1) with next (+1):
    if ( associated(prev) ) then
      prev%next => next
    else
      ! it is the head that will be deleted, move it:
      ino%head => next
    endif

    ! connect next (+1) with previous (-1):
    if ( associated(next) ) then
      next%prev => prev
    endif

    ! delete the block in the middle (0):

    ! update global memory counters:
    MEMLOG(-size_of (blk))
    mem = mem - size_of (blk)

    ! update inode size:
    ino%size = ino%size - size_of (blk)

    ! actually deallocate:
    ASSERT(allocated(blk%sto))
    deallocate(blk)
    ! DONT nullify(blk)

    ! return blk(out) => blk(in)%next
    blk => next

    ! for debug:
    ino%blk_count  = ino%blk_count - 1
    !lk%blk_number = ino%blk_count
  end subroutine delblk

  function firstblk(ino,mask,mode) result(blk)
    ! returns first blk of ino that matches mask
    ! swaps the file in, if necessary
    implicit none
    type(inode)           , intent(inout) :: ino
    type(block), pointer                  :: blk
    integer(IK)           , intent(in)    :: mask(:) ! (MAXMETA)
    character(len=1)      , intent(in)    :: mode ! 'k' or 'd'
    optional :: mask, mode
    ! *** end of interface **

    character(len=MAXHASH) :: hash

    ! 1) if swap is in use make sure to bring relevant data
    ! into memory:
    if( associated(ino%swp_list) )then ! i.e. swap in use
      ! FIXME: cannot yet iterate over swap files with mask=* ...
      ASSERT(present(mask))
      ASSERT(size(mask)==MAXMETA)
      hash = mkhash(mask)
      ASSERT(hash/=ANYHASH)

      if( hash /= ino%cur_hash )then
        ! the pices we need to return maybe not in phys mem., therefore:

        ! A) drop the current contents:
        if( ino%cur_hash /= 'NULL' )then
          !print *,MyID,'RGS: SWAP OUT ',ino%cur_hash
          call swapOutOne(ino,ino%cur_hash)
          ino%cur_hash = 'NULL'
        endif

        ! B) swap-in the needed contents:
        !print *,MyID,'RGS: SWAP IN  ',hash,' need',mask
        call swapInOne(ino,hash,mode)
        ! keep or delete the file according to ``mode''
      endif
    endif

    ! 2) find first matching block:
    blk => ino%head
    do while ( associated(blk) )
      if ( match_meta(blk%meta, mask) ) EXIT
      blk => blk%next
    enddo
  end function firstblk

  subroutine nextblk(blk, mask)
    implicit none
    type(block), pointer :: blk ! returns blk => blk%next
    integer(IK), intent(in) :: mask(:) ! (MAXMETA)
    optional :: mask
    ! *** end of interface **

    ASSERT(associated(blk))
    ! one step forward:
    blk => blk%next

    ! more steps forward if no match:
    do while ( associated(blk) )
      if ( match_meta(blk%meta, mask) ) EXIT ! also matches any .not.present(mask)
      blk => blk%next
    enddo
  end subroutine nextblk

  subroutine delswp(swp)
    ! delete the whole list starting at swp
    implicit none
    type(swap_t), pointer              :: swp
    ! *** end of interface **

    type(swap_t), pointer :: nxt

    do while( associated(swp) )
      nxt => swp%next
      deallocate(swp)
      swp => nxt
    enddo
  end subroutine delswp

  subroutine delswp1(swp,hash)
    ! delete the list entry with hash, return (new) head
    implicit none
    type(swap_t), pointer              :: swp
    character(len=MAXHASH), intent(in) :: hash
    ! *** end of interface **

    type(swap_t), pointer :: pre,cur

    ASSERT(associated(swp))

    ! find the entry with the hash:
    pre => NULL()
    cur => swp
    do while( associated(cur) )
      if( cur%hash == hash ) exit
      pre => cur
      cur => cur%next
    enddo
    ! was it found?:
    ASSERT(associated(cur))

    ! connect previous with next:
    if( associated(pre) )then
      pre%next => cur%next
    else
      ! it is the head that will be deleted, save new head:
      swp => cur%next
    endif

    ! actually deallocate:
    deallocate(cur)
    ! DONT nullify(cur)
  end subroutine delswp1

  function size_of (blk) result(siz)
    ! returns size of a block
    implicit none
    type(block), intent(in) :: blk
    integer(IK)             :: siz
    ! *** end of interface **

    ASSERT(allocated(blk%sto))
    siz = size(blk%sto)
  end function size_of
  !*************************************************************

  !*************************************************************
  subroutine mem_unlink(ino,mask)
    implicit none
    type(inode)           , intent(inout) :: ino
    integer(IK)           , intent(in)    :: mask(:) ! (MAXMETA)
    optional :: mask ! if not present, matches anything
    ! *** end of interface ***

    type(block), pointer :: blk
    integer(IK) :: cblk,cdel,cdir
    logical     :: check

    cblk = 0
    cdel = 0
    cdir = 0
    blk => ino%head
    do while( associated(blk) )
      if ( match_meta(blk%meta, mask) )then ! match_meta() matches always if .not.present(mask)
        ! count:
        cdel = cdel + 1
        if( dirty(blk) ) cdir = cdir + 1
        call delblk(ino,blk) ! also iterates blk => blk%next
      else
        blk => blk%next
      endif
      cblk = cblk + 1
    enddo
    !WRITE(*,'(X,A4,"RGS: unlink: del",I6," blk: dirty=",I6," left=",I6," size=",I6)') MyID,cdel,cdir,cblk-cdel,ino%size
    ASSERT(ino%blk_count==cblk-cdel)

    check = .not.present(mask) ! same as ANYHASH
    if( present(mask) )then
      ASSERT(size(mask)==MAXMETA)
      check = ( mkhash(mask) == ANYHASH )
    endif
    if( check )then
      ! these must have been pointed to => NULL():
      ASSERT(.not.associated(ino%head))
      ASSERT(ino%blk_count==0)
    endif
  end subroutine mem_unlink

  subroutine unlink(ino)
    implicit none
    type(inode), intent(inout) :: ino
    ! *** end of interface ***

    call mem_unlink(ino) ! unlinks all

    ! delete the swap list if any:
    call delswp(ino%swp_list)
  end subroutine unlink
  !*************************************************************

  !*************************************************************
  subroutine rg_xx_get(der,irr,i,j,S,T,V,O,typ,mode)
    use dimensions, only: dimens, dimoff
    USE DEBUG, only : show
    implicit none
    integer(IK)     , intent(in)  :: der
    integer(IK)     , intent(in)  :: irr
    integer(IK)     , intent(in)  :: i,j
    real(RK)        , intent(out) :: S(:,:) ! (N,N)
    real(RK)        , intent(out) :: T(:,:) ! (N,N)
    real(RK)        , intent(out) :: V(:,:) ! (N,N)
    real(RK)        , intent(out) :: O(:,:) ! (N,N)
    integer(IK)     , intent(in)  :: typ
    character(len=1), intent(in)  :: mode ! 'k'eep or 'd'elete
    optional :: T,V,O,typ
    optional :: mode
    ! *** end of interface ***

    integer(IK), parameter :: UNCONTRACTED=0 ! to pass to dimens(...)
    integer(IK) :: N
    integer(IK) :: tp!,sy,ii,jj
    integer(IK) :: ua,la,ub,lb
    integer(IK) :: da,oa,db,ob
    type(block), pointer :: blk
    integer(IK) :: mask(MAXMETA)
    character(len=MAXHASH) :: hash


    FPP_TIMER_START(tot)
    N = dimens(irr,UNCONTRACTED)

    if(present(T))then
    ASSERT(present(V))
    ASSERT(present(O))
    ASSERT(N==size(S,1))
    ASSERT(N==size(S,2))
    ASSERT(N==size(T,1))
    ASSERT(N==size(T,2))
    ASSERT(N==size(V,1))
    ASSERT(N==size(V,2))
    ASSERT(N==size(O,1))
    ASSERT(N==size(O,2))
    else
    ASSERT(present(typ))
    endif

    ! some quadrupels may be missing (=zero)?
    S = 0.0_rk
    if(present(T)) T = 0.0_rk
    if(present(V)) V = 0.0_rk
    if(present(O)) O = 0.0_rk

    ! ``mask'' is a metadata "pattern" of interesting blocks,
    ! only derivative coordinates (and symmetry) do matter:
    mask = meta(-1,irr,-1,-1,-1,-1,i,j)

    ! initiate: blk => FS(der)%head
    blk => firstblk(FS(der),mask,mode)
    ! FIXME: dont hurry with mode=d when hash is not uniqe!

    DO WHILE( associated(blk) )
ASSERT(match_meta(blk%meta,mask))

          ! continue with blocks of
          ! - given symmetry irr, and
          ! - derivative indices i and j

          tp  = blk%meta(ITY) !typ
          ua  = blk%meta(IUA) !ua
          la  = blk%meta(ILA) !la
          ub  = blk%meta(IUB) !ub
          lb  = blk%meta(ILB) !lb

          da = dimens(irr,ua,la,UNCONTRACTED)
          oa = dimoff(irr,ua,la,UNCONTRACTED)
          db = dimens(irr,ub,lb,UNCONTRACTED)
          ob = dimoff(irr,ub,lb,UNCONTRACTED)

          ASSERT(size(blk%sto)==da*db)

          if( present(typ) )then
            if ( tp == typ )then
              call cp_data(oa,da,ob,db,blk%sto,S)
            endif
          else
          select case ( tp )
          case ( RGOVRL )
            call cp_data(oa,da,ob,db,blk%sto,S)
          case ( RGKNTC )
            call cp_data(oa,da,ob,db,blk%sto,T)
          case ( RGNUCL )
            call cp_data(oa,da,ob,db,blk%sto,V)
          case ( RGPVSP )
            call cp_data(oa,da,ob,db,blk%sto,O)
          case ( RGPSEU, RGOTHR )
            ! skip silently
          case default
            ABORT('no such typ!')
          end select
          endif

       ! iterate blk => blk%next:
       call nextblk(blk, mask)
    ENDDO

    ! free memory after reading:
    if( present(mode) )then
      select case(mode)
      case ('d')

        ! free consumed data (prevent swapping out):
        call mem_unlink(FS(der),mask)

        ! on disk file may have been deleted after swapInOne, check first:
        if( swp_inlist(FS(der)%swp_list,hash) )then
!         call swp_unlink(hash)
          WARN('delete '//hash)
        endif
        ! FIXME: when several derivatives result in the same hash!

      case ('k')
        ! do nothing
      case default
        ABORT('no such case')
      end select
    endif

    call rg_close('k') ! keep
!   call show('SXX',S)
!   call show('TXX',T)
    FPP_TIMER_STOP(tot)
  contains
    subroutine cp_data(oa,da,ob,db,sto,ints)
      implicit none
      integer(IK), intent(in)    :: oa,da,ob,db
      real(RK)   , intent(in)    :: sto(da,db)
      real(RK)   , intent(inout) :: ints(:,:)
      ! *** end of interface ***

      !all r_data_rec(ints(1+oa:da+oa,1+ob:db+ob))

      ints(1+oa:da+oa,1+ob:db+ob) = ints(1+oa:da+oa,1+ob:db+ob) &
                                  + sto(:,:)

      ! symmetrize:
      if( oa /= ob )then
        ints(1+ob:db+ob,1+oa:da+oa) = transpose(ints(1+oa:da+oa,1+ob:db+ob))
      else
        ASSERT(da==db)
      endif
!     print *,'ggg<',typ, irr, ii,'[',ua,la,'x',ub,lb,']<->(',1+oa,':',da+oa,'x',1+ob,':',db+ob,')'&
!            , sum(ints(1+oa:da+oa,1+ob:db+ob)) &
!            , sum(ints(1+ob:db+ob,1+oa:da+oa))
    end subroutine cp_data
  end subroutine rg_xx_get
  !*************************************************************

  !*************************************************************
  subroutine cp_data_rec(ints,sto)
    implicit none
    real(RK)   , intent(in)  :: ints(:,:,:,:) ! (NEA,NEB,NIA,NIB)
    real(RK)   , intent(out) :: sto(size(ints))
    ! integrals in a legacy storage format, where
    ! NEA, NEB = number of exponents at (UA,LA) and (UB,LB)
    ! NIA, NIB = number of symmetry independent functions ...
    ! *** end of interface ***

    integer(IK)            :: ea,eb,ia,ib,iab

    ! (NEA,NEB,NIA,NIB) -> (NEA,NIA,NEB,NIB)
    ! ( D1, D2, D3, D4) -> ( D1, D3, D2, D4) 
    iab = 0
    do ib=1,size(ints,4)
    do eb=1,size(ints,2)
    do ia=1,size(ints,3)
    do ea=1,size(ints,1)
      iab = iab + 1
      sto(iab) = ints(ea,eb,ia,ib)
    enddo
    enddo
    enddo
    enddo
  end subroutine cp_data_rec
  !*************************************************************

  subroutine rg_bcast(der)
    ! broadcast data over processors all-to-all
    use comm, only: comm_parallel, comm_bcast, comm_reduce
    use gradient_data_module, only: gradient_rel_index!(n_gr,n_irr)
    implicit none
    integer(IK), intent(in) :: der
    ! *** end of interface ***

    integer(IK) :: n_ir,s,n_gr,i

    if( .not. comm_parallel() ) RETURN

    FPP_TIMER_START(tot)
    FPP_TIMER_START(bcst)

    !print*, MyID,'RGS: rg_bcast enter'
    SELECT CASE(der)

    CASE (0) ! broadcast/distribute integrals between processors

      !print*, MyID,'RGS: rg_bcast: bcast all INTS'
      call bcast_inode(FS(0),meta(-1,-1,-1,-1,-1,-1,-1,-1))
      ! will break at first at firstblk(...,ANYHASH)
      ! that results from this pattern if parts of FS(0) are in swap

      call oom()

    CASE (1) ! broadcast/distribute gradients between processors

      n_ir = size(gradient_rel_index,2) ! NEW SHAPE
      n_gr = size(gradient_rel_index,1) ! NEW SHAPE

      !print*, MyID,'RGS: rg_bcast: bcast all ',n_gr,'*',n_ir,'GRADS'
      do s=1,n_ir
        do i=1,n_gr

          !print*, MyID,'RGS: rg_bcast: bcast GRADS',i,'irr',s
          call bcast_inode(FS(1),meta(-1,s,-1,-1,-1,-1,i,0))

          call oom()
        enddo
      enddo
    CASE DEFAULT
      ABORT('not implemented')
    END SELECT

    ! collect info on integral types (doesnt hurt to do many times):
    call comm_reduce(flag) ! sum all on 1
    call comm_bcast(flag)  ! bcast sum to all

    FPP_TIMER_STOP(bcst)
    FPP_TIMER_STOP(tot)
    !print*, MyID,'RGS: rg_bcast exit, time=',FPP_TIMER_VALUE(bcst)
  end subroutine rg_bcast

  subroutine bcast_inode(ino, mask)
    use comm
    implicit none
    type(inode), intent(inout)   :: ino
    integer(IK), intent(in)      :: mask(:) ! (MAXMETA)
    ! *** end of interface ***

    type(block), pointer     :: blk
    integer(IK)              :: p, np, rank

    ASSERT(size(mask)==MAXMETA)

    if(.not. comm_parallel()) return

    call comm_barrier() ! paranoya

    np   = comm_size()
    rank = comm_rank()

    DO p = 0, np - 1 ! each processor sends its data on its turn, other receive.
      if (rank == p) then ! my turn to send data

        ! initiate: blk => head
        blk => firstblk (ino, mask)

        do while (associated (blk))

          ! loop over matching blocks, skip marked (bcasted) blocks:
          if (match_meta (blk % meta, mask) .and. .not. marked (blk)) then
            !WRITE(*,'(X,A4,"RGS: bcast: SEND blk",8I3," from",I3," to all")') MyID, blk%meta, rank

            ! on sender side, blk is intent(in):
            call bcast_block (blk, root=p)
          endif

          call nextblk (blk, mask)
        enddo
        !WRITE(*,'(X,A4,"RGS: bcast: SEND NULL blk from",I3," to all")') MyID, rank
        ! send NULL blk to terminate the data stream:
        call bcast_block (blk, root=p)

      else ! we are receiving
        ASSERT(rank/=p)

        do ! until receiving NULL blk...

          ! on receiver side, blk is intent(out), so this is redundant:
          blk => NULL()

          ! receive (maybe NULL) blk from ``p'':
          call bcast_block (blk, root=p)

          ! exit on NULL blk:
          if (.not. associated (blk)) EXIT ! the do-loop

          !WRITE(*,'(X,A4,"RGS: bcast: RECV blk",8I3," from",I3)') MyID, blk%meta, p

          ! mark the block as "bcasted" in order not to send it back again:
          call mark (blk)

          ! append marked blk to inode, update inode size:
          call insblk (ino, blk)
        enddo
        !WRITE(*,'(X,A4,"RGS: bcast: RECV NULL blk from",I3)') MyID, p
      endif

      ! not to mix up data streams:
      call comm_barrier()
    ENDDO ! over processors
    call comm_barrier() ! paranoya
  end subroutine bcast_inode

  subroutine rg_distr(der)
    ! distribute (second) derivatives over processors
    ! as they will be required in modules/relgrads
    use comm, only: comm_parallel, comm_barrier
    use gradient_data_module, only: &
        gradient_der_index!(n_gr,n_gr,n_irr)
    implicit none
    integer(IK), intent(in) :: der
    ! *** end of interface ***

    integer(IK) :: s,i,j,n_ir,n_gr
    integer(IK) :: dst
    integer(IK) :: patt(MAXMETA)

    if( .not. comm_parallel() ) RETURN

    FPP_TIMER_START(tot)
    FPP_TIMER_START(dist)

    SELECT CASE(der)
    CASE (2) ! split second derivatives between the processors:
      n_ir = size(gradient_der_index,3)
      n_gr = size(gradient_der_index,1)
ASSERT(n_gr==size(gradient_der_index,2))

      do s=1,n_ir
        do j=1,n_gr
          do i=1,j
            ! destination to collect data on, convert to base-0:
            dst = gradient_der_index(i,j,s) - 1

            ! ``patt'' is a metadata "pattern" of interesting blocks
            patt = meta(-1,s,-1,-1,-1,-1,i,j)

            !print*, MyID,'RGS: rg_distr: gather',i,j,'on',dst
            call gather_inode(FS(2), patt, dst, mode='d')
            ! delete files after slurping 

            call oom()
          enddo
        enddo
      enddo
    CASE DEFAULT
      ABORT('not implemented')
    END SELECT
    FPP_TIMER_STOP(dist)
    FPP_TIMER_STOP(tot)
  end subroutine rg_distr

  subroutine gather_inode(ino, mask, dst, mode)
    ! gather data matching ``mask'' on ``dst''
    use comm
    implicit none
    type(inode), intent(inout)   :: ino
    integer(IK), intent(in)      :: mask(:) ! (MAXMETA)
    integer(IK), intent(in)      :: dst ! destination base-0
    character(len=1), intent(in) :: mode ! 'k' or 'd'
    ! *** end of interface ***

    type(block), pointer :: blk
    integer(IK)          :: np, rank

    if(.not. comm_parallel()) return

    np   = comm_size()
    rank = comm_rank()

    IF( rank /= dst )THEN ! WE ARE SENDING DATA ...
      ! initiate: blk => head
      blk => firstblk(ino, mask, mode)

      do while( associated(blk) )
         ! loop over matching blocks:
         if( .not. match_meta(blk%meta, mask) )then ! NEXT and CYCLE
           call nextblk(blk, mask)
           CYCLE
         endif

         !WRITE(*,'(X,A4,"RGS: gather: SEND blk",8I3," from",I3," to",I3)') MyID,blk%meta,rank,dst
         call send_block(blk, dst)

         ! destroy and proceed:
         call delblk(ino, blk)
         ! NOTE: delblk() also iterates blk => blk%next
         ! FIXME: %next does not necessarily match (and may be NULL)!
         !        The NULL-case must be handeled separately when one
         !        wants to iterate over SEVERAL swap files!
      enddo
      !WRITE(*,'(X,A4,"RGS: gather: SEND NULL blk from",I3," to",I3)') MyID,rank,dst
      ! send NULL blk to terminate the data stream:
      call send_block(blk, dst)

    ELSE ! I AM RECEIVING
      ASSERT(rank==dst)

      np = np - 1 ! number of procs to receive data from
      do while ( np > 0 )

        ! receive (maybe NULL) blk from anywhere:
        call recv_block(blk, COMM_ANY_SOURCE)

        if( associated(blk) )then
          !WRITE(*,'(X,A4,"RGS: gather: RECV blk",8I3)') MyID,blk%meta
          ! append new blk, update inode size:
          call insblk(ino, blk)
        else
          !WRITE(*,'(X,A4,"RGS: gather: RECV NULL blk")') MyID
          ! one more processor signalled the end of data:
          np = np - 1
        endif
      enddo
    ENDIF
    ! not to mix up the inter-process data flow:
    call comm_barrier()
  end subroutine gather_inode

  subroutine bcast_block (blk, root)
    !
    ! Can send blk => NULL()
    !
    use comm, only: comm_rank, comm_bcast, comm_parallel
    implicit none
    type (block), pointer :: blk
    integer (IK), intent (in) :: root
    ! *** end of interface ***

    integer :: siz

    if (.not. comm_parallel()) return

    if (comm_rank() == root) then
       siz = -1
       if (associated (blk)) then
          ASSERT(allocated(blk%sto))
          siz = size (blk % sto)
       endif

       call comm_bcast (siz, root)
       if (siz == -1) RETURN

       call comm_bcast (blk % sto, root)
       call comm_bcast (blk % meta, root)
    else
       blk => NULL()

       call comm_bcast (siz, root)
       if (siz == -1) RETURN

       allocate (blk)
       allocate (blk % sto(siz))

       call comm_bcast (blk % sto, root)
       call comm_bcast (blk % meta, root)

       ! update global memory counters:
       MEMLOG(+size_of (blk))
       mem = mem + size_of (blk)
    endif
  end subroutine bcast_block

  subroutine send_block(blk, dst)
    ! can send blk => NULL()
    use comm
    implicit none
    type(block), pointer    :: blk
    integer(IK), intent(in) :: dst
    ! *** end of interface ***

    integer(IK) :: siz

    if( .not. comm_parallel() ) RETURN

    siz = -1
    if( associated(blk) ) siz = size(blk%sto)

    call comm_send(siz, dst)

    if( siz /= -1 )then
      call comm_send(blk%sto, dst)
      call comm_send(blk%meta, dst)
    endif
  end subroutine send_block

  subroutine recv_block(blk, src)
    use comm
    USE_MEMLOG
    implicit none
    type(block), pointer       :: blk
    integer(IK), intent(in)    :: src ! could be ANY_SOURCE
    ! *** end of interface ***

    integer(IK) :: siz, source

    if( .not. comm_parallel() ) RETURN

    call comm_recv(siz, src)
    source = comm_source() ! in case src == ANY_SOURCE

    blk => NULL()
    if( siz == -1 ) RETURN

    allocate(blk)
    allocate(blk%sto(siz))

    ! recv from the actual source, not ANY_SOURCE:
    call comm_recv(blk%sto, source)
    call comm_recv(blk%meta, source)

    ! update global memory counters:
    MEMLOG(+size_of (blk))
    mem = mem + size_of (blk)
  end subroutine recv_block

  !*************************************************************
  subroutine rg_open(mode)
    implicit none
    character(len=1), intent(in)  :: mode ! w, r, rw
    ! *** end of interface ***

  end subroutine rg_open
  !*************************************************************

  !*************************************************************
  subroutine rg_close(mode)
    ! called as rg_close('k'eep)   from integral_shutdown()
    ! and    as rg_close('d'elete) at the end of integral_trafo()
    implicit none
    character(len=1), intent(in)    :: mode ! k, d
    ! *** end of interface ***

    integer(IK) :: der
    integer(IK) :: chk

    select case( mode )
    case ('k')
      !print *,MyID,'RGS: rg_close(k): '
#ifdef FPP_DEBUG
      call oom_stat()
#endif
    case ('d')
      !print *,MyID,'RGS: unlink(',mode,'): entered'
      !print *,MyID,'RGS: memory used:',mem*8.0/1024**2,'M (before deallocation)'
      chk = mem
      do der=0,MAXORD
        !print *,MyID,'RGS:   |->FS(',der,')=',FS(der)%size*8.0/1024**2,'M'
        chk = chk - FS(der)%size
        call unlink(FS(der))
      enddo
      !print *,MyID,'RGS: memory check:',chk*8.0/1024**2,'M (after deallocation)'
      ASSERT(mem==0)
      !print *,MyID,'RGS: unlink(',mode,'): MEM=',mem

      !print *,MyID,"RGS: timings, total    =",FPP_TIMER_VALUE(tot)
      !print *,MyID,"RGS:           |- disk =",FPP_TIMER_VALUE(disk)
      !print *,MyID,"RGS:           |- bcst =",FPP_TIMER_VALUE(bcst)
      !print *,MyID,"RGS:           `- dist =",FPP_TIMER_VALUE(dist)
!     print *,MyID,"disk__=", disk__, "diskr_=", diskr_

      ! reset flags:
      flag = 0
    case default
      ABORT('rg_close: mode?')
    end select
  end subroutine rg_close
  !*************************************************************

  function rg_flag(typ) result(yes)
    ! indicates presense of this integral type
    implicit none
    integer(IK), intent(in) :: typ
    logical                 :: yes ! result
    ! *** end of interface ***

    ASSERT(typ>=1)
    ASSERT(typ<=RGOTHR)
    yes = ( flag(typ) > 0 )
  end function rg_flag

  subroutine swapOut(ino,siz)
    ! Swaps out and deletes blocks
    ! until the inode size is below ``siz''
    ! ( set siz=0 if you dont want to stop! )
    implicit none
    type(inode), intent(inout) :: ino
    integer(IK), intent(in)    :: siz
    ! *** end of interface ***

!   ! first swap out blocks into already existing files:
!   call swapOutOld(ino,siz)

    ! continue swapping out other blocks:
    call swapOutAny(ino,siz)
  end subroutine swapOut

  subroutine swapOutAny(ino,siz)
    ! Swaps out and deletes blocks (matching ino%head)
    ! until the inode size is below ``siz''
    ! ( set siz=0 if you dont want to stop! )
    implicit none
    type(inode), intent(inout) :: ino
    integer(IK), intent(in)    :: siz
    ! *** end of interface ***

    type(block), pointer   :: blk
    character(len=MAXHASH) :: hash

    ! loop until size reduces to siz:
    do while ( associated(ino%head) .and. ino%size > siz )

      ! take the first block:
      blk => ino%head

      ! create its hash signature:
      hash = mkhash(blk%meta)

      ! swap out all blocks with the same hash:
      call swapOutOne(ino,hash)
    enddo
  end subroutine swapOutAny

! subroutine swapOutOld(ino,siz)
!   ! Swaps out blocks that match files already existing
!   ! on disk. Delete clean blocks after writing.
!   ! STOP SWAPPING if the inode size is below ``siz''
!   ! ( set siz=0 if you dont want to stop! )
!   implicit none
!   type(inode), intent(inout) :: ino
!   integer(IK), intent(in)    :: siz
!   ! *** end of interface ***

!   type(swap_t), pointer   :: swp

!   ! loop over swap entries:
!   swp => ino%swp_list
!   do while( associated(swp) .and. ino%size > siz )
!     ! swap out blocks for this swap entry:
!     call swapOutOne(ino,swp%hash)
!     ! proceed to next swap entry:
!     swp => swp%next
!   enddo
! end subroutine swapOutOld

  subroutine swapOutOne(ino,hash)
    ! Writes the ``dirty'' blocks matching ``hash'' into file
    ! (appends if the file already exists).
    ! After being written the blocks are deleted from
    ! memory. The ``clean'' matching blocks (that have already
    ! been written or fetched from disk) are also deleted.
    implicit none
    type(inode)           , intent(inout) :: ino
    character(len=MAXHASH), intent(in)    :: hash
    ! *** end of interface ***

    type(block), pointer   :: blk
    character(len=1)       :: fmode
    integer(IK)            :: iou = -1
    integer(IK)            :: siz0,siz1

    ! check if this file/hash was already created on disk:
    if( swp_inlist(ino%swp_list,hash) )then
      ! open EXISTING file for APPEND:
      fmode = 'a'
      !print *,MyID,'RGS: APPENDING to swap file: ',hash
    else
      ! open NEW file for WRITE:
      fmode = 'w'
      ! .. and append this hash to the list:
      call swp_insert(ino%swp_list,hash)
      !print *,MyID,'RGS: CREATING the swap file: ',hash
    endif

    ! remeber old size:
    siz0 = ino%size

    blk => ino%head
    ! loop over the blocks with hash == hash:
    do while ( associated(blk) )
      ! filter blocks by hash:
      if( .not. match_hash(mkhash(blk%meta),hash) )then
        ! cycle to next block:
        blk => blk%next
        CYCLE
      endif

      if( dirty(blk) )then

        ! before first write, open the file for WRITE or APPEND:
        if( iou < 0 )then
          call swp_open(iou,hash,fmode) ! sets iou > 0
        endif

        ! clear dirty bit as it will be saved on disk soon:
        call clear_dirty(blk)

        ! append this block to the file:
        call swp_write(iou,blk)
      endif

      ! delete just written or .not.dirty() blocks:
      call delblk(ino,blk) ! and iterate blk=>blk%next too!
    enddo

    ! close the swap file:
    if( iou > 0 )then
      call swp_close(iou,'k')
    endif

    ! get the new size:
    siz1 = ino%size
    !write(*,'(X,A4,"RGS: FLUSH",I8," of ",A9," from",I8," to",I8," (",A,")")') &
         !MyID,siz0-siz1,hash,siz0,siz1,fmode
  end subroutine swapOutOne

  subroutine swapInOne(ino,hash,mode)
    implicit none
    type(inode)           , intent(inout) :: ino
    character(len=MAXHASH), intent(in)    :: hash
    character(len=1), intent(in)          :: mode ! k or d
    optional :: mode
    ! *** end of interface ***

    type(block), pointer   :: blk
    integer(IK)            :: iou = -1
    character(len=1)       :: xmode ! k or d
    integer(IK)            :: siz0,siz1

!   ASSERT(.not.associated(ino%head))

    ASSERT(hash/=ANYHASH)
    ! so far only one file maybe swapped in:
    if( ino%cur_hash/='NULL' )then
      ABORT('swapInOne: hash='//ino%cur_hash)
    endif
    ASSERT(ino%cur_hash=='NULL')

    if( .not.swp_inlist(ino%swp_list, hash) )then
      RETURN
    endif

    ! keep the disk copy by default:
    xmode = 'k'
    if( present(mode) ) xmode = mode

    ! remember initial size:
    siz0 = ino%size

    ! open swap file:
    call swp_open(iou,hash,'r')

    do ! read until NULL blk returned on EOF

      ! read/allocate and return NULL if EOF:
      call swp_read(iou,blk)
      if( .not.associated(blk) ) EXIT ! the loop

      ! blocks on disk cannot be dirty:
      ASSERT(.not.dirty(blk))

      ! set dirty bit if deleting swap:
      if( xmode == 'd' ) call mark_dirty(blk)

      ! append new block to the list, update inode size:
      call insblk(ino,blk)
    enddo

    ! close file (keep or delete according to xmode):
    call swp_close(iou,xmode)

    select case ( xmode )
    case ('k')
      ! mark the current content of the memory:
      ino%cur_hash = hash
    case ('d')
      ! mark memory dirty:
      ino%cur_hash = 'NULL'
      ! find and delete the list entry:
      call delswp1(ino%swp_list,hash)
      !write(*,'(X,A4,"RGS: swapInOne: DELETED ",A9)') MyID,hash
    case default
      ABORT('swapIn action?')
    end select

    ! get the new size:
    siz1 = ino%size
    !write(*,'(X,A4,"RGS: SLURP",I8," of ",A9," from",I8," to",I8," (",A,")")') &
         !MyID,siz1-siz0,hash,siz0,siz1,xmode
  end subroutine swapInOne

  function swp_inlist(swp_list,hash) result(yes)
    implicit none
    type(swap_t), pointer              :: swp_list
    character(len=MAXHASH), intent(in) :: hash
    logical                            :: yes ! result
    ! *** end of interface ***

    type(swap_t), pointer :: swp

    yes = .false.

    swp => swp_list
    do while ( associated( swp ) )
      yes = ( hash == swp%hash )
      if( yes ) EXIT!from loop
      swp => swp%next
    enddo
  end function swp_inlist

  subroutine swp_insert(swp_list,hash)
    implicit none
    type(swap_t), pointer              :: swp_list
    character(len=MAXHASH), intent(in) :: hash
    ! *** end of interface ***

    type(swap_t), pointer :: new,last

    ! construct new list entry:
    allocate( new )
    new%hash = hash

    if( .not. associated(swp_list) )then
      ! point list head to new element:
      swp_list => new
    else
      ! find the last list element:
      last => swp_list

      do while ( associated( last%next ) )
        last => last%next
      enddo
      ! append new entry:
      last%next => new
    endif
  end subroutine swp_insert

  subroutine swp_write(iou,blk)
    use io, only: file_write
    implicit none
    integer(IK), intent(inout)         :: iou
    type(block), intent(in)            :: blk
    ! *** end of interface ***

    ASSERT(iou>0)
    FPP_TIMER_START(disk)

    call file_write(iou,size(blk%sto))
    call file_write(iou,blk%meta)
    call file_write(iou,blk%bits)
    call file_write(iou,blk%sto)

    mem_swp = mem_swp + size_of (blk)
    FPP_TIMER_STOP(disk)
  end subroutine swp_write

  subroutine swp_read(iou,blk)
    USE_MEMLOG
    use io, only: file_read
    implicit none
    integer(IK), intent(inout)         :: iou
    type(block), pointer               :: blk ! returns NULL on EOF
    ! *** end of interface ***

    integer(IK) :: siz
    integer(IK) :: stat

    ASSERT(iou>0)
    FPP_TIMER_START(disk)

    blk => NULL()

    call file_read(iou,siz,stat)
    if( stat/=0 ) goto 999 ! NULL blk

    ! create new blk:
    allocate(blk)
    allocate( blk%sto(siz) )
    MEMLOG(+size(blk%sto))
    mem = mem + size(blk%sto)

    call file_read(iou,blk%meta)
    call file_read(iou,blk%bits)
    call file_read(iou,blk%sto)
999 CONTINUE
    FPP_TIMER_STOP(disk)
  end subroutine swp_read

  subroutine swp_open(iou,fname,mode)
    use comm, only: comm_rank
    use filename_module, only: tmpfile
    use io, only: file_open
    implicit none
    integer(IK), intent(inout)   :: iou
    character(len=*), intent(in) :: fname
    character(len=1), intent(in) :: mode ! r,w,a
    ! *** end of interface ***

    character(len=LEN(fname)+7) :: filename

    ASSERT(iou==-1)
    FPP_TIMER_START(disk)

    ! prepend proc index to filenames:
    filename = itoa(comm_rank())//fname//'.swp'

    call file_open(iou, tmpfile(filename), mode)

    ASSERT(iou>0)
    FPP_TIMER_STOP(disk)
  end subroutine swp_open

  subroutine swp_close(iou,mode)
    use io, only: file_close
    implicit none
    integer(IK), intent(inout)         :: iou
    character(len=1), intent(in)       :: mode
    ! *** end of interface ***

    if( iou < 0 ) RETURN
    FPP_TIMER_START(disk)
    call file_close(iou,mode)
    iou = -1
    FPP_TIMER_STOP(disk)
  end subroutine swp_close

  subroutine swp_unlink(fname)
    implicit none
    character(len=*), intent(in) :: fname
    ! *** end of interface ***

    integer(IK) :: iou=-1

    call swp_open(iou,fname,'r')
    call swp_close(iou,'d')
  end subroutine swp_unlink

  function mkhash(meta) result(hash)
    ! creates a string of characters from
    ! symmetry- and derivative indices
    implicit none
    integer(IK), intent(in) :: meta(:) ! (MAXMETA)
    character(len=MAXHASH)  :: hash
    ! *** end of interface ***

    integer(IK) :: i,j

    ASSERT(size(meta)==MAXMETA)
    ASSERT(MAXHASH>=3*3)
#if 0
    ! each second derivative into its own file:
    hash = itoa(meta(ISY))//itoa(meta(IXI))//itoa(meta(IXJ))
    ! WARNING: results into too many files for large systems
#else
    ! group second derivatives by their position in hessian matrix (i,j):
    i = meta(IXI)
    j = meta(IXJ)

    ! fixed block size 12!
    if(i>0) i = 1+i/12 ! zero stays zero, negative stays negative
    if(j>0) j = 1+j/12 ! zero stays zero, negative stays negative

    hash =  itoa(meta(ISY))//itoa(i)//itoa(j)

    ! even better make hash by the access patern/sequence:
    !
    ! ij = 0
    ! do j=1,maxj
    !   do i=1,j
    !     ij = ij + 1
    !   enddo
    ! enddo
    !
    ! for that one needs an estimate for maxj!
#endif
  end function mkhash

  function match_hash(h1,h2) result(yes)
    ! understands ``?'' wildcard
    ! if .not.present(h2) then always return true
    implicit none
    character(len=MAXHASH), intent(in)  :: h1,h2
    logical                             :: yes ! reult
    optional :: h2
    ! *** end of interface ***

    integer(IK) :: i

    yes = .true.
    if( .not. present(h2) ) RETURN ! all match by default

    do i=1,MAXHASH
      if( h1(i:i) == '?' ) CYCLE
      if( h2(i:i) == '?' ) CYCLE
      yes = yes .and. ( h1(i:i) == h2(i:i) )
      if(.not.yes) EXIT
    enddo
  end function match_hash

  function match_meta(m1, m2) result(yes)
    ! understands ``?'' wildcard as -1
    implicit none
    integer(IK), intent(in)  :: m1(:) ! (MAXMETA)
    integer(IK), intent(in)  :: m2(:) ! (MAXMETA)
    logical                  :: yes ! reult
    optional :: m2
    ! *** end of interface ***

    integer(IK) :: i

    ASSERT(size(m1)==MAXMETA)

    yes = .true.
    if( .not. present(m2) ) RETURN ! all match by default

    ASSERT(size(m2)==MAXMETA)

    do i=1,MAXMETA
      if( m1(i) == -1 ) CYCLE
      if( m2(i) == -1 ) CYCLE
      yes = yes .and. ( m1(i) == m2(i) )
      if(.not.yes) EXIT
    enddo
  end function match_meta

  function itoa(i) result(a)
    implicit none
    integer(IK), intent(in) :: i
    character(len=3)        :: a
    ! *** end of interface ***

    if(i==-1)then
      a = '???'
      RETURN
    endif
    ASSERT(i>=0)
    ASSERT(i<1000)
    write(a,'(i3.3)') i
  end function itoa

  function dirty(blk) result(yes)
    implicit none
    type(block), intent(in)    :: blk
    logical                    :: yes ! result
    ! *** end of interface **

    yes = IAND( blk%bits, IDIRTY ) /= 0
  end function dirty

  subroutine mark_dirty(blk)
    implicit none
    type(block), intent(inout) :: blk
    ! *** end of interface **

    blk%bits = IOR( blk%bits, IDIRTY )
  end subroutine mark_dirty

  subroutine clear_dirty(blk)
    implicit none
    type(block), intent(inout) :: blk
    ! *** end of interface **

    blk%bits = IAND( blk%bits, NOT(IDIRTY) )
  end subroutine clear_dirty

  function marked(blk) result(yes)
    implicit none
    type(block), intent(in)    :: blk
    logical                    :: yes ! result
    ! *** end of interface **

    yes = IAND( blk%bits, IMARK ) /= 0
  end function marked

  subroutine mark(blk)
    implicit none
    type(block), intent(inout) :: blk
    ! *** end of interface **

    blk%bits = IOR( blk%bits, IMARK )
  end subroutine mark

  subroutine clear_mark(blk)
    implicit none
    type(block), intent(inout) :: blk
    ! *** end of interface **

    blk%bits = IAND( blk%bits, NOT(IMARK) )
  end subroutine clear_mark

  !--------------- End of module ----------------------------------
end module relgrads_store
