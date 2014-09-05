!
! ParaGauss,  a program package  for high-performance  computations of
! molecular systems
!
! Copyright (C) 2014     T. Belling,     T. Grauschopf,     S. Krüger,
! F. Nörtemann, M. Staufer,  M. Mayer, V. A. Nasluzov, U. Birkenheuer,
! A. Hu, A. V. Matveev, A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman,
! D. I. Ganyushin,   T. Kerdcharoen,   A. Woiterski,  A. B. Gordienko,
! S. Majumder,     M. H. i Rotllant,     R. Ramakrishnan,    G. Dixit,
! A. Nikodem, T. Soini, M. Roderus, N. Rösch
!
! This program is free software; you can redistribute it and/or modify
! it under  the terms of the  GNU General Public License  version 2 as
! published by the Free Software Foundation [1].
!
! This program is distributed in the  hope that it will be useful, but
! WITHOUT  ANY   WARRANTY;  without  even  the   implied  warranty  of
! MERCHANTABILITY  or FITNESS FOR  A PARTICULAR  PURPOSE. See  the GNU
! General Public License for more details.
!
! [1] http://www.gnu.org/licenses/gpl-2.0.html
!
! Please see the accompanying LICENSE file for further information.
!
!=====================================================================
! Public interface of module
!=====================================================================
module fit_trafo_tapes
  !-------------------------------------------------------------------
  !
  ! Copyright (c) Alexei Matveev
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

# include "def.h"
  use type_module, only: IK => i4_kind, RK => r8_kind
  use filename_module, only: strlen => filename_namelengthmax
  use readwriteblocked_module, only: TapeHandle=>readwriteblocked_tapehandle
#ifdef FPP_DEBUG
  use error_module, only: MyID
#endif
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------


  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------

  public init_fit_trafo_tapes,&
       & done_fit_trafo_tapes,&
       & Split3cInt,Merge3cInt,&
       & OpenTapes, CloseTapes,&
       & ReadTape,  WriteTape

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !
  ! FIXME: unifiy with istore.f90 functionality?
  !
  type, private :: FileHandle
     character(len=strlen) :: name
     integer(IK)           :: length, n_ff
  end type FileHandle

  !------------ Declaration of constants and variables ---------------

!!$  integer(IK),parameter :: MAX_ALLOC_REALS = 16*1048576_ik
!!$  integer(IK),parameter :: NDAMP  =   8_ik
  integer(IK),parameter :: NMAXRW = 256_ik

!!$  logical,parameter     :: verbose    = .true.
  logical,parameter     :: keep_files = .false.


  integer(IK),parameter    ::&
       & RE      = 1,&
       & IM      = 2

  integer(IK), parameter, public :: &
       TH_KIJ = 312, &
       TH_IJK = 123 ! storage modes for 3c-ints

  integer(IK),parameter,public ::&
       & TH_COUL = 1,&
       & TH_SS   = 2,&
       & TH_SX   = 3,&
       & TH_RS   = 4,&
       & TH_RX   = 5,&
       & TH_MAX  = 5

  integer(IK)              :: &
       & TH_NUM  = 1 ! or 5

  type(FileHandle)         ::&
       & i_Coul(2,TH_MAX),&
       & o_Coul(2)
  type(FileHandle),allocatable ::&
       & i_Tapes(:,:,:),&
       & o_Tapes(:,:),&
       & Tmp(:,:,:)

  type(TapeHandle),pointer ::&
       & i_Th(:,:,:),&
       & o_Th(:,:)

  integer(IK),parameter    :: NO_IRREP = -1

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine init_fit_trafo_tapes()
    use filename_module, only: tmpfile
    use quadrupel_fname, qfn_init => init
    use dimensions, only: u => IrrUBasDimSpor, c => IrrBasDimSpor
    use spin_orbit_module, only: is_on, op_RelFit
    implicit none
    ! *** end of interface ***

    integer(IK) :: memstat
    integer(IK) :: i, r, irr, n_irr

    ! number of irreps:
    ASSERT(size(u)==size(c))
    n_irr = size(c)

    call qfn_init()

    if ( .not. is_on(op_RelFit) )then
       TH_NUM = 1
    else
       TH_NUM = 5
    endif

    allocate(&
         & i_Tapes(n_irr, RE:IM, TH_NUM), &
         & o_Tapes(n_irr, RE:IM), &
         & Tmp(n_irr, RE:IM, TH_NUM), &
         & STAT=memstat)
    ASSERT(memstat==0)

    allocate(&
         & i_Th(n_irr, RE:IM, TH_NUM), &
         & o_Th(n_irr, RE:IM), &
         & STAT=memstat)
    ASSERT(memstat==0)

    if(keep_files)then
       o_Coul(RE)%name = tmpfile("coul_real.DAT")
       o_Coul(IM)%name = tmpfile("coul_imag.DAT")
    else
       o_Coul(RE)%name = tmpfile("coul_real.dat")
       o_Coul(IM)%name = tmpfile("coul_imag.dat")
    endif

    !
    ! Names and sizes of "coul_real.dat" and similar:
    !
    do i = TH_COUL, TH_NUM ! 1:1 or 1:5
       do r = RE, IM ! 1:2
          i_Coul(r, i)%name = fn3(i, r, "dat")
          i_Coul(r, i)%n_ff = nf(i)
          i_Coul(r, i)%length = nf(i) * sum(u * (u + 1) / 2)
       enddo
    enddo

    !
    ! Names and sizes of "coul_real_1.t01" and similar:
    !
    do irr = 1, n_irr
       do r = RE, IM
          o_Tapes(irr, r)%name = fn4(TH_COUL, r, irr, "t03")
          o_Tapes(irr, r)%length = nf(TH_COUL) * c(irr) * (c(irr) + 1) / 2 ! contracted!

          do i = TH_COUL, TH_NUM
             Tmp(irr, r, i)%name = fn4(i, r, irr, "t01")
             Tmp(irr, r, i)%length = nf(i) * u(irr) * (u(irr) + 1) / 2

             i_Tapes(irr, r, i)%name = fn4(i, r, irr, "t02")
             i_Tapes(irr, r, i)%length = nf(i) * u(irr) * (u(irr) + 1) / 2
          enddo
       enddo
    enddo

    contains

      function nf(t) result(n)
        !
        ! Number of fit Functions
        !
        use int_send_2cob3c_spor, only: IX_CHFIT_S, IX_CHFIT_R2, MaxFitIndex
        implicit none
        integer(IK), intent(in) :: t
        integer(IK) :: n
        ! *** end of interface ***

        select case (t)
        case (TH_COUL)
           n = MaxFitIndex()
        case (TH_SS, TH_SX)
           n = MaxFitIndex(IX_CHFIT_S)
        case (TH_RS, TH_RX)
           n = MaxFitIndex(IX_CHFIT_R2)
        case default
           n = -1
           ABORT("no such case")
        end select
      end function nf

      function pf(i) result(s)
        !
        ! Prefix
        !
        implicit none
        integer(IK), intent(in) :: i
        character(len=10) :: s
        ! *** end of interface ***

        select case (i)
        case (TH_COUL)
           s = "coul"
        case (TH_SS)
           s = "rcoul_pvsp"
        case (TH_SX)
           s = "rcoul_pvxp"
        case (TH_RS)
           s = "r2_pvsp"
        case (TH_RX)
           s = "r2_pvxp"
        case default
           s = "xxxxxxxxxx"
           ABORT("no such case")
        end select
      end function pf

      function fn3(i, r, ext) result(file)
        !
        ! File Name such as coul_real.dat
        !
        implicit none
        integer(IK), intent(in) :: i, r
        character(len=*), intent(in) :: ext
        character(len=strlen) :: file
        ! *** end of interface ***

        file = tmpfile(trim(pf(i)) // "_" // reim(r) // "." // ext)
      end function fn3

      function fn4(i, r, irr, ext) result(file)
        !
        ! File Name, such as coul_real_1.dat
        !
        implicit none
        integer(IK), intent(in) :: i, r, irr
        character(len=*), intent(in) :: ext
        character(len=strlen) :: file
        ! *** end of interface ***

        file = qfilename(trim(pf(i)) // "_" // reim(r), irr, ext)
      end function fn4

      function reim(r) result(s)
        !
        ! real or imag
        !
        implicit none
        integer(IK), intent(in) :: r
        character(len=4) :: s
        ! *** end of interface ***

        select case (r)
        case (RE)
           s = "real"
        case (IM)
           s = "imag"
        case default
           s = "xxxx"
           ABORT("no such case")
        end select
      end function reim
  end subroutine init_fit_trafo_tapes

  subroutine done_fit_trafo_tapes()
    use error_module, only: error
    use quadrupel_fname, qfn_done=>done
    implicit none
    ! *** end of interface ***

    integer(IK) :: memstat

    deallocate(i_Tapes,o_Tapes,Tmp,STAT=memstat)
    call error(memstat,"ftt/done_fit_trafo_tapes: dealloc failed 1")

    deallocate(i_Th,o_Th,STAT=memstat)
    call error(memstat,"ftt/done_fit_trafo_tapes: dealloc failed 2")

    call qfn_done()
  end subroutine done_fit_trafo_tapes

  subroutine ReadTape(th, irr, buf_re, buf_im)
    use readwriteblocked_module
    implicit none
    integer(IK),intent(in) :: th
    integer(IK),intent(in) :: irr
    real(RK), intent(out) :: buf_re(:)
    real(RK), intent(out) :: buf_im(:)
    ! *** end of interface ***

    ASSERT(th>=TH_COUL)
    ASSERT(th<=TH_NUM)
    if(irr.eq.NO_IRREP)return

    call readwriteblocked_read(buf_re, i_Th(irr, RE, th))
    call readwriteblocked_read(buf_im, i_Th(irr, IM, th))
  end subroutine ReadTape

  subroutine WriteTape(irr, buf_re, buf_im)
    use readwriteblocked_module
    implicit none
    integer(IK),intent(in) :: irr
    real(RK), intent(in) :: buf_re(:)
    real(RK), intent(in) :: buf_im(:)
    ! *** end of interface ***

    if(irr.eq.NO_IRREP)return

    call readwriteblocked_write(buf_re, o_Th(irr, RE))
    call readwriteblocked_write(buf_im, o_Th(irr, IM))
  end subroutine WriteTape

  subroutine OpenTapes(mode,irr)
    use readwriteblocked_module
    implicit none
    integer(IK),intent(in) :: mode
    integer(IK),intent(in) :: irr
    optional :: irr
    ! *** end of interface ***

    integer(IK) :: BlockLength,r,i

    if(present(irr))then
       if(irr.eq.NO_IRREP)return
    endif

    BlockLength = readwriteblocked_blocklength()

    select case(mode)
    case (TH_IJK)
    ASSERT(present(irr))
    do i=1,TH_NUM
       do r=RE,IM
          call readwriteblocked_startread(&
               & trim(i_Tapes(irr,r,i)%name),&
               & i_Th(irr,r,i),&
               & blocklength=BlockLength,&
               & total_length=i_Tapes(irr,r,i)%length)
       enddo
    enddo

    do r=RE,IM
       call readwriteblocked_startwrite(&
            & trim(o_Tapes(irr,r)%name),&
            & o_Th(irr,r),&
            & blocklength=BlockLength)
    enddo
    case (TH_KIJ)
       ASSERT(size(i_Th,1)>0)
       print *,'OpenTapes: TH_NUM=',TH_NUM
       do i=1,TH_NUM
          do r=RE,IM
             print *,'opening KJI ',trim(i_Coul(r,i)%name)
             call readwriteblocked_startread(&
                  & trim(i_Coul(r,i)%name),&
                  & i_Th(1,r,i),&
                  & blocklength=BlockLength,&
                  & total_length=i_Coul(r,i)%length)
          enddo
       enddo
       ! not going to write anything
    case default
      ABORT('no such mode')
    end select
  end subroutine OpenTapes

  subroutine CloseTapes(mode,irr)
    use readwriteblocked_module
    implicit none
    integer(IK),intent(in) :: mode
    integer(IK),intent(in) :: irr
    optional :: irr
    ! *** end of interface ***

    integer(IK) :: r,chklength,i
    integer(IK) :: ii
    character(len=6) :: keep="keep"

    select case(mode)
    case (TH_IJK)
    do r=RE,IM
       ASSERT(present(irr))
       if(irr.eq.NO_IRREP)return
       ii = irr
       call readwriteblocked_stopwrite(o_Th(ii,r),total_length=chklength)
       DPRINT MyID,'ftt/CloseTapes: chklength=',chklength
       o_Tapes(ii,r)%length = chklength
    enddo

    if(keep_files)then
       keep="keep"
    else
       keep="delete"
    endif

    case default
       ii = 1
       print *,' didnt write anything'
    end select

    do i=1,TH_NUM
       do r=RE,IM
          call readwriteblocked_stopread(i_Th(ii,r,i),status=keep)
       enddo
    enddo

  end subroutine CloseTapes

  subroutine Split3cInt()
    use error_module
    use quadrupel_fname
    use dimensions,&
         & IrrUDim => IrrUBasDimSpor,&
         & n_irr  => number_of_irreps
    implicit none
    ! *** end of interface ***

    integer(IK) :: i,irr,r,n_ob,n_ij,n_ff

    DPRINT MyID,"ftt/Split3cInt: entered"

    !  n_ff = MaxFitIndex()

    do i=1,TH_NUM
       do r=RE,IM
          ! split tape:
          call split_file(i_Coul(r,i),Tmp(1:n_irr,r,i))

          ! transpose tapes:
          do irr=1,n_irr

             n_ob = IrrUDim(irr)
             n_ij = (n_ob * (n_ob + 1))/2

             n_ff = i_Coul(r,i)%n_ff

             call transpose_file(n_ff,n_ij,Tmp(irr,r,i),i_tapes(irr,r,i))
          enddo
       enddo
    enddo
  end subroutine Split3cInt

  subroutine Merge3cInt()
    use error_module
    use quadrupel_fname
    use int_send_2cob3c_spor, only: MaxFitIndex
    use dimensions,&
         & IrrCDim => IrrBasDimSpor,&
         & n_irr  => number_of_irreps
    implicit none
    ! *** end of interface ***

    integer(IK) :: irr,r,n_ff,n_ob,n_ij

    DPRINT MyID,"ftt/Merge3cInt: entered"

    n_ff = MaxFitIndex()

    do r=RE,IM
       ! transpose files:
       do irr=1,n_irr

          n_ob = IrrCDim(irr)
          n_ij = (n_ob * (n_ob + 1))/2

          if(n_ff*n_ij/=o_Tapes(irr,r)%length)&
               & call error("ftt/Merge3cInt: o_Tapes length?")

          call transpose_file(n_ij,n_ff,o_Tapes(irr,r),Tmp(irr,r,TH_COUL))
!!!!!! Alexei, what should I do with this Tmp file
       enddo

       ! concat files:
       call concat_files(Tmp(:,r,TH_COUL),o_Coul(r))
    enddo
  end subroutine Merge3cInt

  subroutine split_file(from,to)
    use error_module
    use readwriteblocked_module
    implicit none
    type(FileHandle),intent(in) :: from,to(:)
    ! ..., to(n_files)
    ! *** end of interface ***

    integer(IK)      :: memstat
    type(TapeHandle) :: th_from,th_to
    integer(IK)      :: flength,MaxFLength,n_files,f,rest,n
    integer(IK)      :: IBlockLength,OBlockLength,BufSize
    integer(IK)      :: chklength
    real(RK),pointer :: buf(:)

    DPRINT MyID,"ftt/split_file: entered"

    if(from%length/=sum(to(:)%length))&
         & call error("ftt/split_file: length?")

    n_files    = size(to)
    MaxFLength = maxval(to(:)%length)

    flength = from%length

    IBlockLength = readwriteblocked_blocklength()

    BufSize = MaxAllocSize(MaxFLength)

    DPRINT MyID,'ftt/split_file: BufSize=',BufSize

    allocate(buf(BufSize),STAT=memstat)
    call error(memstat,"ftt/split_file: alloc buf failed")

!!$    DPRINT MyID,'ftt/split_file: open from ...'

    call readwriteblocked_startread(&
         & trim(from%name),&
         & th_from,&
         & blocklength=IBlockLength,&
         & total_length=from%length)

    do f=1,n_files

       OBlockLength = min(IBlockLength,to(f)%length)

!!$       DPRINT MyID,'ftt/split_file: open to(',f,'), OBlockLength=',OBlockLength

       call readwriteblocked_startwrite(&
            & trim(to(f)%name),&
            & th_to,&
            & blocklength=OBlockLength)

       ! rewriting:
       rest = to(f)%length
       do while(rest>0)
          n = min(rest,BufSize)
!!$          DPRINT MyID,'ftt/split_file: rest=',rest,'read in n=',n
          call readwriteblocked_read(buf(1:n), th_from)
          call readwriteblocked_write(buf(1:n),th_to)
          rest = rest - n
       enddo

       call readwriteblocked_stopwrite(th_to,total_length=chklength)
       if(to(f)%length/=chklength)call error("ftt/split_file: chklength?")
    enddo

    if(keep_files)then
       call readwriteblocked_stopread(th_from,status="keep")
    else
       call readwriteblocked_stopread(th_from,status="delete")
    endif

    deallocate(buf,STAT=memstat)
    call error(memstat,"ftt/split_file: dealloc buf failed")
  end subroutine split_file

  subroutine concat_files(from,to)
    use error_module
    use readwriteblocked_module
    implicit none
    type(FileHandle),intent(in)    :: from(:)
    type(FileHandle),intent(inout) :: to
    ! ..., from(n_files)
    ! *** end of interface ***

    integer(IK)      :: memstat
    type(TapeHandle) :: th_from,th_to
    integer(IK)      :: flength,MaxFLength,n_files,f,rest,n
    integer(IK)      :: IBlockLength,OBlockLength,BufSize
    integer(IK)      :: chklength
    real(RK),pointer :: buf(:)

    DPRINT MyID,"ftt/concat_files: entered"

    n_files    = size(from)
    MaxFLength = maxval(from(:)%length)

    flength   = sum(from(:)%length)
    to%length = flength

    OBlockLength = readwriteblocked_blocklength()

    BufSize = MaxAllocSize(MaxFLength)

    DPRINT MyID,'ftt/concat_files: BufSize=',BufSize

    allocate(buf(BufSize),STAT=memstat)
    call error(memstat,"ftt/concat_files: alloc buf failed")

!!$    DPRINT MyID,'ftt/concat_files: open to ...'

    call readwriteblocked_startwrite(&
         & trim(to%name),&
         & th_to,&
         & blocklength=OBlockLength)

    do f=1,n_files

       IBlockLength = min(OBlockLength,from(f)%length)

!!$       DPRINT MyID,'ftt/concat_files: open from(',f,'), IBlockLength=',IBlockLength

       call readwriteblocked_startread(&
            & trim(from(f)%name),&
            & th_from,&
            & blocklength=IBlockLength,&
            & total_length=from(f)%length)

       ! rewriting:
       rest = from(f)%length
       do while(rest>0)
          n = min(rest,BufSize)
!!$          DPRINT MyID,'ftt/concat_files: rest=',rest,'read in n=',n
          call readwriteblocked_read(buf(1:n), th_from)
          call readwriteblocked_write(buf(1:n),th_to)
          rest = rest - n
       enddo

!!$       if(keep_files)then
!!$          call readwriteblocked_stopread(th_from,status='keep')
!!$       else
!!$          call readwriteblocked_stopread(th_from,status='delete')
!!$       endif
       call readwriteblocked_stopread(th_from,status='delete')
    enddo

    call readwriteblocked_stopwrite(th_to,total_length=chklength)
    call error(chklength/=to%length,"ftt/concat_files: chklength?")

    DPRINT MyID,'ftt/concat_files: chklength=',chklength

    deallocate(buf,STAT=memstat)
    call error(memstat,"ftt/split_file: dealloc buf failed")
  end subroutine concat_files

  subroutine transpose_file(n,m,in,out)
    use error_module
    use readwriteblocked_module
    implicit none
    integer(IK),intent(in)         :: n,m
    type(FileHandle),intent(in)    :: in
    type(FileHandle),intent(inout) :: out
    ! *** end of interface ***

    type(TapeHandle) :: th_in,th_out
    integer(IK)      :: memstat
    integer(IK)      :: flength,BlockLength,BufSize,&
         & n_lim,m_lim,n_rest,m_rest,n_io,m_io,n_s,n_e,m_s,m_e
    real(RK),pointer :: abuf(:,:),bbuf(:,:)

    DPRINT MyID,"ftt/transpose_file: entered"

    call error(in%length/=n*m,"ftt/transpose_file: length?")

    flength    = in%length
    out%length = in%length

    if( n*m == 0 )then
       WARN('transpose_file: zero length, short-curcuit')
       goto 999 ! short-curcuit for zero-length files
    endif

    BlockLength = readwriteblocked_blocklength()

!!$    DPRINT MyID,'ftt/transpose_file: BlockLength(0)=',BlockLength

    BlockLength = min(BlockLength,flength)

!!$    DPRINT MyID,'ftt/transpose_file: BlockLength(1)=',BlockLength

    BufSize = MaxAllocSize(flength)

    DPRINT MyID,'ftt/transpose_file: BufSize(0)=',BufSize

    m_lim = BufSize/n
    n_lim = BufSize/m

    if(n_lim*m_lim < n*m)then
       call warn("ftt/transpose_file: WARNING-WARNING-WARNING-WARNING")
       call warn("ftt/transpose_file: low memory case")
       call warn("ftt/transpose_file: WARNING-WARNING-WARNING-WARNING")
    endif

    DPRINT MyID,'ftt/transpose_file: n_lim=',n_lim,'m_lim=',m_lim

    if(n_lim * m_lim == 0)&
         & call error("ftt/transpose_file: low of memory 0")

    BufSize = n * m_lim !<<< round it

    DPRINT MyID,'ftt/transpose_file: BufSize(2)=',BufSize

    if(BufSize.eq.0)&
         & call error("ftt/transpose_file: low of memory 1")

    if(flength/BufSize>NMAXRW)&
         & call error("ftt/transpose_file: low of memory 2")

    allocate(&
         & abuf(n,m_lim),&
         & bbuf(n_lim,m),&
         & STAT=memstat)
    call error(memstat,"ftt/transpose_file: alloc failed")

    ! rewriting:
    call OpenOut()

    n_rest = n
    n_s    = 1
    do while(n_rest>0)

!!$       DPRINT MyID,'ftt/transpose_file: n_rest=',n_rest

       call OpenIn()

       n_io = min(n_lim,n_rest)
       n_e  = n_s + n_io - 1

!!$       DPRINT MyID,'ftt/transpose_file: n_io=',n_io

       m_rest = m
       m_s    = 1
       do while(m_rest>0)

!!$          DPRINT MyID,'ftt/transpose_file: m_rest=',m_rest

          m_io = min(m_lim,m_rest)
          m_e  = m_s + m_io - 1

!!$          DPRINT MyID,'ftt/transpose_file: m_io=',m_io

          call ReadIn(abuf(1:n,1:m_io))

!!!bug:          bbuf(1:n_io,1:m_io) = abuf(n_s:n_e,m_s:m_e)
          bbuf(1:n_io,m_s:m_e) = abuf(n_s:n_e,1:m_io)

!!$          DPRINT MyID,'ftt/transpose_file: extracting (',n_s,':',n_e,',',m_s,':',m_e,')'

          m_rest = m_rest - m_io
          m_s    = m_s    + m_io
       enddo

       n_rest = n_rest - n_io
       n_s    = n_s    + n_io

       call CloseIn(DEL=(n_rest.eq.0))

       call WriteOut(bbuf(1:n_io,1:m))
    enddo

    call CloseOut()

    deallocate(abuf,bbuf,STAT=memstat)
    call error(memstat,"ftt/transpose_file: dealloc failed")

    RETURN   ! *** RETURN POINT FOR NORMAL CASE        ***

999 continue ! *** SHORT-CIRCUIT FOR ZERO-LENGTH FILES ***
    call OpenOut()
    call OpenIn()
    call CloseIn(DEL=.true.)
    call CloseOut()
    RETURN   ! *** RETURN POINT FOR ZERO-LENGTH FILES  ***

  contains
    !
    ! pieces of code, not subs:
    !
    subroutine OpenIn()
      implicit none
      ! *** end of interface ***

!!$      DPRINT MyID,'ftt/transpose_file: open in  ...'
      call readwriteblocked_startread(&
           & trim(in%name),&
           & th_in,&
           & blocklength=BlockLength,&
           & total_length=flength)
    end subroutine OpenIn

    subroutine CloseIn(del)
      implicit none
      logical,optional,intent(in) :: del
      ! *** end of interface ***

      logical :: delete

      delete = .not. keep_files
      if(present(del)) delete = del

!!$      DPRINT MyID,'ftt/transpose_file: close in  ...'

!!$      if(keep_files)then
!!$         call readwriteblocked_stopread(th_in,status="keep")
!!$      else
!!$         call readwriteblocked_stopread(th_in,status="delete")
!!$      endif
      if(delete)then
         call readwriteblocked_stopread(th_in,status="delete")
      else
         call readwriteblocked_stopread(th_in,status="keep")
      endif
    end subroutine CloseIn

    subroutine OpenOut()
      implicit none
      ! *** end of interface ***

!!$      DPRINT MyID,'ftt/transpose_file: open out  ...'
      call readwriteblocked_startwrite(&
           & trim(out%name),&
           & th_out,&
           & blocklength=BlockLength)
    end subroutine OpenOut

    subroutine CloseOut()
      use error_module
      implicit none
      ! *** end of interface ***

      integer(IK) :: chklength

!!$      DPRINT MyID,'ftt/transpose_file: close out  ...'
      call readwriteblocked_stopwrite(th_out,total_length=chklength)

      if(chklength/=flength)call error("ftt/transpose_file: chklength?")
    end subroutine CloseOut

    subroutine ReadIn(a2d)
      implicit none
      real(RK),intent(out) :: a2d(:,:)
      ! *** end of interface ***

      real(RK),dimension(size(a2d)) :: a1d

!!$      DPRINT MyID,'ftt/transpose_file: read in  ...'
      call readwriteblocked_read(a1d,th_in)

      a2d = reshape(a1d,shape(a2d))
    end subroutine ReadIn

    subroutine WriteOut(a2d)
      implicit none
      real(RK),intent(in) :: a2d(:,:)
      ! *** end of interface ***

      real(RK),dimension(size(a2d)) :: a1d

      a1d = reshape(transpose(a2d),(/size(a2d)/))

!!$      DPRINT MyID,'ftt/transpose_file: write out  ...'
      call readwriteblocked_write(a1d,th_out)
    end subroutine WriteOut
  end subroutine transpose_file

  function MaxAllocSize(need) result(n)
    use error_module
    use spin_orbit_module, only: max_alloc_reals
    implicit none
    integer(IK),intent(in) :: need
    integer(IK)            :: n !<<< result
    ! *** end of interface ***

!!$    integer(IK)      :: memstat
!!$    real(RK),pointer :: buf(:)

    integer(IK) :: max_size

    max_size  = 8*1048576_ik

    if(max_alloc_reals().ne.-1)then
       max_size = max_alloc_reals()
    endif

    if( need > max_size )then
       n = max_size
       call warn("ftt/MaxAllocSize: renice mem")
    else
       n = need
    endif

!!$    allocate(buf(n),STAT=memstat)
!!$
!!$    do while(memstat/=0 .and. n>0)
!!$
!!$       n = n/NDAMP
!!$       allocate(buf(n),STAT=memstat)
!!$    enddo
!!$
!!$    deallocate(buf,STAT=memstat)
!!$    call error(memstat,"ftt/MaxAllocSize: dealloc failed")
!!$
!!$    DPRINT MyID,'ftt/MaxAllocSize: down',need/n,'times, n=',n
!!$
!!$    call error(n.eq.0,"ftt/MaxAllocSize: low of memory")
  end function MaxAllocSize

  !--------------- End of module -------------------------------------
end module fit_trafo_tapes
