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
module scm
!
! Copyright (c) 2011-2013 Alexei Matveev
!
use iso_c_binding
implicit none
private

!
! Short structs  with size matching that  of some integer  type may be
! passed/returned  in  registers. That  depends  on  the platform  ABI
! convention.   BIND(C)  on  the  type(scm_t)  and  relevant  function
! interfaces requests the Fortran compiler to emulate the behaviour of
! (some) C compiler.
!
! http://gcc.gnu.org/onlinedocs/gcc/Code-Gen-Options.html
! www.x86-64.org/documentation/abi.pdf
!
! FIXME: Intel  compiler versions 11.1 and  12.1 break ABI  by using a
! custom convention of returning a struct. In this case one might want
! to replace all instances of type(scm_t) by integer(c_intptr_t).
!
type, public, bind(c) :: scm_t
   private
   integer (c_intptr_t) :: ptr
end type scm_t
!define SCM_T integer (c_intptr_t)
!define SCM_T type (scm_t)
!define SCM_T type (c_ptr)

interface
   !
   ! No object satises more than one of the following predicates:
   !
   ! boolean?  pair?  symbol?  number?  char?  string?  vector?
   ! bytevector?  port? procedure?  null?
   !
   function scm_boolean_p (object) result (yes) bind (c)
     !
     ! SCM scm_boolean_p (SCM obj)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: object
     type (scm_t) :: yes
   end function scm_boolean_p

   function scm_symbol_p (object) result (yes) bind (c)
     !
     ! SCM scm_symbol_p (SCM obj)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: object
     type (scm_t) :: yes
   end function scm_symbol_p

   function scm_string_p (object) result (yes) bind (c)
     !
     ! SCM scm_string_p (SCM obj)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: object
     type (scm_t) :: yes
   end function scm_string_p

   function scm_number_p (object) result (yes) bind (c)
     !
     ! SCM scm_number_p (SCM obj)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: object
     type (scm_t) :: yes
   end function scm_number_p

   function scm_pair_p (object) result (yes) bind (c)
     !
     ! SCM scm_pair_p (SCM obj)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: object
     type (scm_t) :: yes
   end function scm_pair_p

   function scm_null_p (object) result (yes) bind (c)
     !
     ! SCM scm_null_p (SCM x)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: object
     type (scm_t) :: yes
   end function scm_null_p

   !
   ! These may be usefull to tell appart the instances of of number?:
   !
   function scm_exact_p (object) result (yes) bind (c)
     !
     ! SCM scm_exact_p (SCM x)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: object
     type (scm_t) :: yes
   end function scm_exact_p

   function scm_inexact_p (object) result (yes) bind (c)
     !
     ! SCM scm_inexact_p (SCM x)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: object
     type (scm_t) :: yes
   end function scm_inexact_p

   function scm_list_p (obj) result (yes) bind (c)
     !
     ! SCM scm_list_p (SCM obj)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: obj
     type (scm_t) :: yes
   end function scm_list_p

   !
   ! Type constructors/accessors, integers, doubles, strings,
   ! logicals, lists:
   !

   !
   ! FIXME: Guile API exports a macro for some of these, we added a
   ! function, see guile-api.c:
   !
   function guile_macro_scm_is_true (object) result (yes) bind (c, name="scm_is_true")
     !
     ! int scm_is_true (SCM obj)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: object
     integer (c_int) :: yes
   end function guile_macro_scm_is_true

   function scm_to_int32 (exact) result (i) bind(c)
     !
     ! int scm_to_int32 (SCM exact)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: exact
     integer (c_int32_t) :: i
   end function scm_to_int32

   function scm_to_int64 (exact) result (i) bind(c)
     !
     ! int64 scm_to_int64 (SCM exact)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: exact
     integer (c_int64_t) :: i
   end function scm_to_int64

   function scm_from_pointer (ptr, finalizer) result (obj) bind(c)
     !
     ! SCM scm_from_pointer (void *ptr, void (*finalizer) (void*))
     !
     import
     implicit none
     type (c_ptr), intent (in), value :: ptr
     type (c_funptr), intent (in), value :: finalizer
     type (scm_t) :: obj
   end function scm_from_pointer

   !
   ! FIXME: Guile API exports a macro for this, see guile-api.c:
   !
   function scm_to_int (exact) result (i) bind (c, name="guile_macro_scm_to_int")
     !
     ! int scm_to_int (SCM exact)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: exact
     integer (c_int) :: i
   end function scm_to_int

   function scm_to_double (inexact) result (d) bind (c)
     !
     ! double scm_to_double (SCM inexact)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: inexact
     real (c_double) :: d
   end function scm_to_double

   function scm_from_latin1_stringn (str, len) result (string) bind (c)
     !
     ! SCM scm_from_latin1_stringn (const char *str, size_t len)
     !
     import
     implicit none
     character (kind=c_char), intent (in) :: str(*)
     integer (c_size_t), intent (in), value :: len
     type (scm_t) :: string
   end function scm_from_latin1_stringn

   function scm_from_locale_stringn (str, len) result (string) bind (c)
     !
     ! SCM scm_from_locale_stringn (const char *str, size_t len)
     !
     import
     implicit none
     character (kind=c_char), intent (in) :: str(*)
     integer (c_size_t), intent (in), value :: len
     type (scm_t) :: string
   end function scm_from_locale_stringn

   function scm_to_locale_stringbuf (str, buf, max_len) result (length) bind (c)
     !
     ! size_t scm_to_locale_stringbuf (SCM str, char *buf, size_t max_len)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: str
     character (kind=c_char) :: buf(*)
     integer (c_size_t), intent (in), value :: max_len
     integer (c_size_t) :: length
   end function scm_to_locale_stringbuf

   function scm_cons (car, cdr) result (pair) bind (c)
     !
     ! SCM scm_cons (SCM car, SCM cdr)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: car, cdr
     type (scm_t) :: pair
   end function scm_cons

   function scm_car (pair) result (car) bind (c)
     !
     ! SCM scm_car (SCM pair)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: pair
     type (scm_t) :: car
   end function scm_car

   function scm_cdr (pair) result (cdr) bind (c)
     !
     ! SCM scm_cdr (SCM pair)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: pair
     type (scm_t) :: cdr
   end function scm_cdr

   function scm_length (list) result (length) bind (c)
     !
     ! SCM scm_length (SCM list)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: list
     type (scm_t) :: length
   end function scm_length

   function scm_string_to_symbol (string) result (symbol) bind (c)
     !
     ! SCM scm_string_to_symbol (SCM string)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: string
     type (scm_t) :: symbol
   end function scm_string_to_symbol

   function scm_symbol_to_string (symbol) result (string) bind (c)
     !
     ! SCM scm_symbol_to_string (SCM symbol)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: symbol
     type (scm_t) :: string
   end function scm_symbol_to_string

   function scm_variable_bound_p (var) result (yes) bind (c)
     !
     ! SCM scm_variable_bound_p (SCM var)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: var
     type (scm_t) :: yes
   end function scm_variable_bound_p

   function scm_variable_ref (variable) result (value) bind (c)
     !
     ! SCM scm_variable_ref (SCM var)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: variable
     type (scm_t) :: value
   end function scm_variable_ref

   function scm_undefined () result (undef) bind (c)
     !
     ! SCM SCM_UNDEFINED
     !
     import
     implicit none
     type (scm_t) :: undef
   end function scm_undefined

   function scm_define (name, val) result (var) bind (c)
     !
     ! SCM scm_define (SCM name, SCM val)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: name, val
     type (scm_t) :: var
   end function scm_define

   function scm_c_define_gsubr (name, req, opt, rst, fcn) result (proc) bind (c)
     !
     ! SCM scm_c_define_gsubr (const char *name, int req, int opt, int rst, fcn)
     !
     import
     implicit none
     character (kind=c_char), intent (in) :: name(*) ! null-terminated, of course
     integer (c_int), intent (in), value :: req, opt, rst
     type (c_funptr), intent (in), value :: fcn
     type (scm_t) :: proc
   end function scm_c_define_gsubr

   subroutine scm_c_export_1 (name) bind (c)
     !
     ! void scm_c_export_1 (const char *name)
     !
     import
     implicit none
     character (kind=c_char), intent (in) :: name(*) ! null-terminated, of course
   end subroutine scm_c_export_1

   function scm_fluid_p (obj) result (yes) bind (c)
     !
     ! SCM scm_fluid_p (SCM)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: obj
     type (scm_t) :: yes
   end function scm_fluid_p

   function scm_fluid_ref (obj) result (maybe) bind (c)
     !
     ! SCM scm_fluid_ref (SCM)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: obj
     type (scm_t) :: maybe
   end function scm_fluid_ref
end interface


interface scm_lookup
   function scm_lookup (name) result (variable) bind (c)
     !
     ! SCM scm_lookup (SCM name)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: name
     type (scm_t) :: variable
   end function scm_lookup

   function scm_module_lookup (mod, name) result (variable) bind (c)
     !
     ! SCM scm_module_lookup (SCM module, SCM name)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: mod, name
     type (scm_t) :: variable
   end function scm_module_lookup
end interface scm_lookup


! FIXME: yet better add the  fortran wrappers for _c_-functions to not
! require   null-ternimated  strings   and  put   them   into  generic
! functions. See e.g. scm_resolve_module() generic interface.
interface
   function scm_c_lookup (name) result (variable) bind (c)
     !
     ! SCM scm_c_lookup (const char *name)
     !
     import
     implicit none
     character (kind=c_char), intent (in) :: name(*) ! null-terminated
     type (scm_t) :: variable
   end function scm_c_lookup

   function scm_c_module_lookup (mod, name) result (variable) bind (c)
     !
     ! SCM scm_c_module_lookup (SCM module, const char *name)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: mod
     character (kind=c_char), intent (in) :: name(*) ! null-terminated
     type (scm_t) :: variable
   end function scm_c_module_lookup
end interface

interface scm_resolve_module
   function scm_resolve_module (name) result (mod) bind (c)
     !
     ! SCM scm_resolve_module (SCM name)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: name ! list of symbols
     type (scm_t) :: mod
   end function scm_resolve_module

   ! Wrapper  for  scm_c_resolve_module()  that accepts  Fortran  (not
   ! terminated) strings.
   module procedure scm_f_resolve_module
end interface scm_resolve_module

interface
   function scm_c_resolve_module (name) result (mod) bind (c)
     !
     ! SCM scm_c_resolve_module (const char *name)
     !
     import
     implicit none
     character (kind=c_char), intent (in) :: name(*) ! null-terminated
     type (scm_t) :: mod
   end function scm_c_resolve_module
end interface

!
! Libguile  defines separate  scm_from_double(),  scm_from_int32() and
! scm_from_int64() with  a macro, scm_from_int() that  resolves to one
! of those.  For Fortran we  define a generic scm_from() that compiler
! will resolve to the correct function depending on the type, kind and
! rank (TKR) the argument:
!
interface scm_from
   function scm_from_int32 (i) result (exact) bind (c)
     !
     ! SCM scm_from_int32 (int  i), resolves from macro scm_from_int()
     ! on 32 bit.
     !
     import
     implicit none
     integer (c_int32_t), intent (in), value :: i
     type (scm_t) :: exact
   end function scm_from_int32

   function scm_from_int64 (i) result (exact) bind (c)
     !
     ! SCM scm_from_int64 (int  i), resolves from macro scm_from_int()
     ! on 64 bit.
     !
     import
     implicit none
     integer (c_int64_t), intent (in), value :: i
     type (scm_t) :: exact
   end function scm_from_int64

   function scm_from_double (d) result (inexact) bind (c)
     !
     ! SCM scm_from_double (double d)
     !
     import
     implicit none
     real (c_double), intent (in), value :: d
     type (scm_t) :: inexact
   end function scm_from_double

   module procedure scm_from_logical  ! logical -> SCM bool
   module procedure scm_from_string   ! string -> SCM string
   module procedure scm_from_c_funptr ! c_funptr -> SCM exact
end interface

!
! Defined assignment can be made  generic as the compiler can dispatch
! based on TKR of both arguments:
!
interface assignment(=)
   module procedure assign_to_int32
   module procedure assign_to_int64
   module procedure assign_to_double
   module procedure assign_to_string
   module procedure assign_to_logical
   module procedure assign_from_int32
   module procedure assign_from_int64
   module procedure assign_from_double
   module procedure assign_from_string
   module procedure assign_from_logical

   module procedure assign_to_c_ptr
   module procedure assign_from_c_funptr
end interface

interface scm_list
   !
   ! SCM scm_eol (); // returns an empty list
   ! SCM scm_list_1 (SCM e1);
   ! SCM scm_list_2 (SCM e1, SCM e2);
   ! SCM scm_list_3 (SCM e1, SCM e2, SCM e3);
   ! SCM scm_list_4 (SCM e1, SCM e2, SCM e3, SCM e4);
   ! SCM scm_list_5 (SCM e1, SCM e2, SCM e3, SCM e4, SCM e5);
   !
   pure function scm_eol () result (empty) bind (c)
     !
     ! SCM SCM_EOL; // macro or constant in C
     !
     import
     implicit none
     type (scm_t) :: empty
   end function scm_eol

   function scm_list_1 (e1) result (list) bind (c)
     import
     implicit none
     type (scm_t), intent (in), value :: e1
     type (scm_t) :: list
   end function scm_list_1

   function scm_list_2 (e1, e2) result (list) bind (c)
     import
     implicit none
     type (scm_t), intent (in), value :: e1, e2
     type (scm_t) :: list
   end function scm_list_2

   function scm_list_3 (e1, e2, e3) result (list) bind (c)
     import
     implicit none
     type (scm_t), intent (in), value :: e1, e2, e3
     type (scm_t) :: list
   end function scm_list_3

   function scm_list_4 (e1, e2, e3, e4) result (list) bind (c)
     import
     implicit none
     type (scm_t), intent (in), value :: e1, e2, e3, e4
     type (scm_t) :: list
   end function scm_list_4

   function scm_list_5 (e1, e2, e3, e4, e5) result (list) bind (c)
     import
     implicit none
     type (scm_t), intent (in), value :: e1, e2, e3, e4, e5
     type (scm_t) :: list
   end function scm_list_5
end interface

interface scm_call
   function scm_call_0 (proc) result (res) bind (c)
     !
     ! SCM scm_call_0 (SCM proc)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: proc
     type (scm_t) :: res
   end function scm_call_0

   function scm_call_1 (proc, arg1) result (res) bind (c)
     !
     ! SCM scm_call_1 (SCM proc, SCM arg1)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: proc, arg1
     type (scm_t) :: res
   end function scm_call_1

   function scm_call_2 (proc, arg1, arg2) result (res) bind (c)
     !
     ! SCM scm_call_2 (SCM proc, SCM arg1, SCM arg2)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: proc, arg1, arg2
     type (scm_t) :: res
   end function scm_call_2

   function scm_call_3 (proc, arg1, arg2, arg3) result (res) bind (c)
     !
     ! SCM scm_call_3 (SCM proc, SCM arg1, SCM arg2, SCM arg3)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: proc, arg1, arg2, arg3
     type (scm_t) :: res
   end function scm_call_3

   function scm_call_4 (proc, arg1, arg2, arg3, arg4) result (res) bind (c)
     !
     ! SCM scm_call_4 (SCM proc, SCM arg1, SCM arg2, SCM arg3, SCM arg4)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: proc, arg1, arg2, arg3, arg4
     type (scm_t) :: res
   end function scm_call_4

   function scm_call_5 (proc, arg1, arg2, arg3, arg4, arg5) result (res) bind (c)
     !
     ! SCM scm_call_4 (SCM proc, SCM arg1, SCM arg2, SCM arg3, SCM arg4, SCM arg5)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: proc, arg1, arg2, arg3, arg4, arg5
     type (scm_t) :: res
   end function scm_call_5

   !
   ! FIXME: extend up to 9 if necessary.
   !
end interface

interface scm_apply
   function scm_apply_0 (proc, arglst) result (res) bind (c)
     !
     ! SCM scm_apply_0 (SCM proc, SCM arglst)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: proc, arglst
     type (scm_t) :: res
   end function scm_apply_0

   ! FIXME: libguile  functions scm_apply() and  scm_apply_1() seem to
   ! have the same interface, but different meaning.
end interface

!
! Predicates:
!
public :: scm_is_symbol         ! SCM obj -> logical
public :: scm_is_string         ! SCM obj -> logical
public :: scm_is_number         ! SCM obj -> logical
public :: scm_is_pair           ! SCM obj -> logical
public :: scm_is_null           ! SCM obj -> logical

public :: scm_is_true           ! SCM obj -> logical
public :: scm_is_exact          ! SCM obj -> logical
public :: scm_is_inexact        ! SCM obj -> logical

public :: scm_list_p            ! SCM obj -> SCM logical

!
! Accessors:
!
public :: scm_to_int            ! SCM int -> integer
public :: scm_to_double         ! SCM double -> double
public :: scm_to_stringbuf      ! SCM string -> (Maybe string, integer)

public :: scm_car               ! SCM pair -> SCM car
public :: scm_cdr               ! SCM pair -> SCM cdr
public :: scm_length            ! SCM list -> SCM int

!
! Constructors:
!
public :: scm_from              ! any -> SCM any
public :: scm_string_to_symbol  ! SCM string -> SCM symbol
public :: assignment(=)         ! any = SCM any, SCM any = any

public :: scm_cons              ! SCM car -> SCM cdr -> SCM pair
public :: scm_list              ! ... -> SCM list

!
! Procedures:
!
public :: scm_call              ! SCM proc -> ... -> SCM res
public :: scm_apply             ! SCM proc -> SCM arglst -> SCM res

!
! Setting/quering environment:
!
interface scm_defined_p
   ! One argument version from this module:
   module procedure scm_defined_p_1 ! SCM symbol -> SCM bool

   ! Two argument version from libguile, SCM symbol -> SCM scope -> SCM bool
   function scm_defined_p_2 (symbol, scope) result (yes) bind (c, name="scm_defined_p")
     !
     ! SCM scm_defined_p (SCM symbol, SCM scope)
     !
     import
     implicit none
     type (scm_t), intent (in), value :: symbol, scope
     type (scm_t) :: yes
   end function scm_defined_p_2
end interface scm_defined_p

public :: scm_define         ! SCM symbol -> SCM value -> SCM variable
public :: scm_defined_p      ! SCM symbol -> [SCM scope] -> SCM bool
public :: scm_lookup ! SCM symbol -> SCM variable or (SCM module, SCM symbol) -> variable
public :: scm_variable_bound_p  ! SCM variable -> SCM bool
public :: scm_variable_ref      ! SCM variable -> SCM value
public :: scm_fluid_ref         ! SCM fluid -> SCM value or #f
public :: scm_undefined         ! () -> SCM_UNDEFINED
public :: scm_define_gsubr      ! (character, integer, integer, integer, c_funptr) -> SCM proc
public :: scm_export            ! (character)
public :: scm_resolve_module ! SCM name -> SCM module or string -> SCM module

contains

  function scm_f_resolve_module (name) result (mod)
    !
    ! Wrapper  for scm_c_resolve_module()  that  accepts Fortran  (not
    ! terminated) strings.
    !
    use iso_c_binding, only: C_NULL_CHAR
    implicit none
    character (len=*), intent (in) :: name
    type (scm_t) :: mod
    ! *** end of interface ***

    mod = scm_c_resolve_module (trim (name) // C_NULL_CHAR)
  end function scm_f_resolve_module

   function scm_is_true (object) result (yes)
     implicit none
     type (scm_t), intent (in), value :: object
     logical :: yes
     ! *** end of interface ***

     yes = guile_macro_scm_is_true (object) /= 0
   end function scm_is_true

   function scm_is_symbol (object) result (yes)
     implicit none
     type (scm_t), intent (in), value :: object
     logical :: yes
     ! *** end of interface ***

     yes = scm_is_true (scm_symbol_p (object))
   end function scm_is_symbol

   function scm_is_string (object) result (yes)
     implicit none
     type (scm_t), intent (in), value :: object
     logical :: yes
     ! *** end of interface ***

     yes = scm_is_true (scm_string_p (object))
   end function scm_is_string

   function scm_is_number (object) result (yes)
     implicit none
     type (scm_t), intent (in), value :: object
     logical :: yes
     ! *** end of interface ***

     yes = scm_is_true (scm_number_p (object))
   end function scm_is_number

   function scm_is_pair (object) result (yes)
     implicit none
     type (scm_t), intent (in), value :: object
     logical :: yes
     ! *** end of interface ***

     yes = scm_is_true (scm_pair_p (object))
   end function scm_is_pair

   function scm_is_null (object) result (yes)
     implicit none
     type (scm_t), intent (in), value :: object
     logical :: yes
     ! *** end of interface ***

     yes = scm_is_true (scm_null_p (object))
   end function scm_is_null

   function scm_is_exact (object) result (yes)
     implicit none
     type (scm_t), intent (in), value :: object
     logical :: yes
     ! *** end of interface ***

     yes = scm_is_true (scm_exact_p (object))
   end function scm_is_exact

   function scm_is_inexact (object) result (yes)
     implicit none
     type (scm_t), intent (in), value :: object
     logical :: yes
     ! *** end of interface ***

     yes = scm_is_true (scm_inexact_p (object))
   end function scm_is_inexact

  function scm_from_string (for_string) result (scm_string)
    !
    ! Fortran string -> SCM string.
    !
    implicit none
    character (len=*), intent (in) :: for_string
    type (scm_t) :: scm_string
    ! *** end of interface **

    integer (c_size_t) :: slen

    slen = len (for_string)
    scm_string = scm_from_locale_stringn (for_string, slen)
  end function scm_from_string

  subroutine scm_to_stringbuf (string, buf, length)
    !
    ! SCM string -> (Maybe (Fortran string), length)
    !
    ! The output string is only complete if length <= len(buf) on
    ! output.
    !
    type (scm_t) :: string
    character (len=*), intent (out) :: buf
    integer, intent (out) :: length ! default kind, not an integer (c_size_t)
    ! *** end of interface **

    integer (c_size_t) :: max_len, size_t_length

    max_len = len (buf)
    size_t_length = scm_to_locale_stringbuf (string, buf, max_len)

    ! if (size_t_length > max_len) then
    !    print *, "scm_to_stringbuf: WARNING! String too long:", size_t_length, ">", max_len
    ! endif

    ! Convert to plain integer, FIXME: overflow?
    length = int (size_t_length, kind = kind (length))

    ! clear trailing garbage, FIXME: should we?
    buf(length+1:max_len) = " "
  end subroutine scm_to_stringbuf

  function scm_define_gsubr (name, req, opt, rst, fcn) result (proc)
    !
    ! Fortranish wrapper (appends C_NULL_CHAR to name) for
    !
    ! SCM scm_c_define_gsubr (const char *name, int req, int opt, int rst, fcn)
    !
    implicit none
    character (len=*) :: name
    integer (c_int), intent (in) :: req, opt, rst
    type (c_funptr), intent (in) :: fcn
    type (scm_t) :: proc
    ! *** end of interface ***

    proc = scm_c_define_gsubr (name // C_NULL_CHAR, req, opt, rst, fcn)
  end function scm_define_gsubr

  subroutine scm_export (name)
    !
    ! Fortranish wrapper (appends C_NULL_CHAR to name) for
    !
    ! void scm_c_export_1 (const char *name)
    !
    implicit none
    character (len=*) :: name
    ! *** end of interface ***

    call scm_c_export_1 (name // C_NULL_CHAR)
  end subroutine scm_export

  function define (key, val) result (var)
    implicit none
    character (len=*), intent (in) :: key
    type (scm_t), intent (in), value :: val
    type (scm_t) :: var
    ! *** end of interface ***

    var = scm_define (scm_string_to_symbol (scm_from (key)), val)
  end function define

  function lookup (key) result (value)
    character (len=*), intent (in) :: key
    type (scm_t) :: value
    ! *** end of interface **

    type (scm_t) :: string, symbol, variable

    string = scm_from (key)
    symbol = scm_string_to_symbol (string)
    variable = scm_lookup (symbol)
    value = scm_variable_ref (variable)
  end function lookup

  function scm_defined_p_1 (symbol) result (yes)
    !
    ! This is one-argument version of scm_defined_p() from libguile.so
    ! offered  not  to prolifirate  use  of scm_undefined()  function.
    ! FIXME:  offer public SCM_UNDEFINED/SCM_UNSPECIFIED  constants to
    ! the users of this module.
    !
    implicit none
    type (scm_t), intent (in) :: symbol ! SCM symbol
    type (scm_t) :: yes                ! SCM bool
    ! *** end of interface ***

    ! FIXME: or should it be *unspecified*?
    yes = scm_defined_p (symbol, scm_undefined ())
  end function scm_defined_p_1

  recursive subroutine display (obj)
    implicit none
    type (scm_t), intent (in) :: obj
    ! *** end of interface **

    character (len=128) :: buf
    integer :: slen

    if (scm_is_number (obj)) then

       if (scm_is_exact (obj)) then
          write (*, '(I6)', advance='no') scm_to_int (obj)
       else
          write (*, '(F12.6)', advance='no') scm_to_double (obj)
       endif

    else if (scm_is_string (obj)) then

       call scm_to_stringbuf (obj, buf, slen)

       if (slen > len (buf)) then
          stop "display: ERROR: string too long!"
       endif

       write (*, "('""', A, '""')", advance='no') buf(1:slen)

    else if (scm_is_symbol (obj)) then

       call scm_to_stringbuf (scm_symbol_to_string (obj), buf, slen)

       if (slen > len (buf)) then
          stop "display: ERROR: symbol too long!"
       endif

       write (*, "('''', A)", advance='no') buf(1:slen)

    else if (scm_is_null (obj)) then

       write (*, '("()")', advance='no')

    else if (scm_is_pair (obj)) then

       write (*, '("(")', advance='no')
       call display (scm_car (obj))
       write (*, '(" . ")', advance='no')
       call display (scm_cdr (obj))
       write (*, '(")")', advance='no')

    else
       write (*, '("???")', advance='no')
    endif
  end subroutine display

  function test (symbol, object) result (out) bind (c)
    !
    ! FIXME:  Note that  even though  this funciton  is  BIND(C), when
    ! compiled  by  Intel  the   resulting  interface  is  NOT  binary
    ! compatible to
    !
    !     intptr_t test (intptr_t symbol, intptr_t object)
    !
    ! so that  calling it from C/Guile  will fail. The  call will even
    ! fail  from Intel,  as  the caller  assumes  the correct  (quoted
    ! above) type signature.
    !
    ! In  effect,  it  appears  that extending  Guile  interpreter  by
    ! fortran funcitons returning SCM  objects is impossible, since it
    ! would not work at least  with one compiler. It works though with
    ! Gfortran.
    !
    implicit none
    type (scm_t), intent (in), value :: symbol, object
    type (scm_t) :: out
    ! *** end of interface ***

    character (len=16) :: buf
    integer :: slen
    type (scm_t) :: var

    call display (symbol)
    write (*, *) ! newline
    call display (object)
    write (*, *) ! newline

    call scm_to_stringbuf (scm_symbol_to_string (symbol), buf, slen)

    if (slen > len (buf)) then
       stop "test: ERROR: symbol too long!"
    end if

    out = scm_cons (scm_eol(), scm_cons (lookup (buf(1:slen)), object))
    call display (out)
    write (*, *)

    var = define ("*string*", scm_from ("hello, world!"))
    var = define ("*double*", scm_from (0.3d0))
    var = define ("*integer*", scm_from (7))
  end function test

  !
  ! Defined  assignment  can  be  made  generic as  the  compiler  can
  ! dispatch based on TKR of both arguments:
  !
  subroutine assign_from_int32 (a, b)
    !
    ! SCM exact = int32
    !
    implicit none
    type (scm_t), intent (out) :: a
    integer (c_int32_t), intent (in) :: b
    ! *** end of interface ***

    a = scm_from_int32 (b)
  end subroutine assign_from_int32

  subroutine assign_to_int32 (a, b)
    !
    ! int32 = SCM exact
    !
    implicit none
    integer (c_int32_t), intent (out) :: a
    type (scm_t), intent (in) :: b
    ! *** end of interface ***

    a = scm_to_int32 (b)
  end subroutine assign_to_int32

  subroutine assign_from_int64 (a, b)
    !
    ! SCM exact = int64
    !
    implicit none
    type (scm_t), intent (out) :: a
    integer (c_int64_t), intent (in) :: b
    ! *** end of interface ***

    a = scm_from_int64 (b)
  end subroutine assign_from_int64

  subroutine assign_to_int64 (a, b)
    !
    ! int64 = SCM exact
    !
    implicit none
    integer (c_int64_t), intent (out) :: a
    type (scm_t), intent (in) :: b
    ! *** end of interface ***

    a = scm_to_int64 (b)
  end subroutine assign_to_int64

  subroutine assign_to_logical (a, b)
    !
    ! logical = SCM bool
    !
    implicit none
    logical, intent (out) :: a
    type (scm_t), intent (in) :: b
    ! *** end of interface ***

    ! scm_is_true() would  be a  lossy conversion. Abort  if b  is not
    ! boolean:
    if (.not. scm_is_true (scm_boolean_p (b))) then
       stop "Not boolean!"
    endif
    a = scm_is_true (b)
  end subroutine assign_to_logical

  subroutine assign_from_logical (a, b)
    !
    ! SCM bool = logical
    !
    implicit none
    type (scm_t), intent (out) :: a
    logical, intent (in) :: b
    ! *** end of interface ***

    a = scm_from_logical (b)
  end subroutine assign_from_logical

  function scm_from_logical (b) result (a)
    implicit none
    logical, intent (in) :: b
    type (scm_t) :: a
    ! *** end of interface ***

    ! FIXME: #t and #f as literals?
    if (b) then
       a = scm_null_p (scm_list ()) ! -> SCM #t
    else
       a = scm_pair_p (scm_list ()) ! -> SCM #f
    endif
  end function scm_from_logical

  subroutine assign_from_double (a, b)
    !
    ! SCM inexact = double
    !
    implicit none
    type (scm_t), intent (out) :: a
    real (c_double), intent (in) :: b
    ! *** end of interface ***

    a = scm_from_double (b)
  end subroutine assign_from_double

  subroutine assign_to_double (a, b)
    !
    ! double = SCM inexact
    !
    implicit none
    real (c_double), intent (out) :: a
    type (scm_t), intent (in) :: b
    ! *** end of interface ***

    a = scm_to_double (b)
  end subroutine assign_to_double

  subroutine assign_from_string (a, b)
    !
    ! SCM string = string
    !
    implicit none
    type (scm_t), intent (out) :: a
    character (len=*), intent (in) :: b
    ! *** end of interface ***

    a = scm_from_string (b)
  end subroutine assign_from_string

  subroutine assign_to_string (a, b)
    !
    ! string = SCM string
    !
    implicit none
    character (len=*), intent (out) :: a
    type (scm_t), intent (in) :: b
    ! *** end of interface ***

    integer (c_size_t) :: lena, lenb

    lena = len (a)
    lenb = scm_to_locale_stringbuf (b, a, lena)

    if (lenb > lena) then
       print *, "assign_to_string: string too long:", lenb, ">", lena
       stop 1
    endif
  end subroutine assign_to_string

  function scm_to_c_ptr (b) result(a)
    !
    ! SCM exact = c_funptr
    !
    use iso_c_binding, only: c_funptr, c_intptr_t
    implicit none
    type (scm_t), intent (in) :: b
    type (c_ptr) :: a
    ! *** end of interface ***

    integer (c_intptr_t) :: c

    ! Defined assignment of an SCM exact to intptr_t:
    c = b

    ! FIXME:  potentially non-portable  coversion of  an integer  to a
    ! (data) pointer:
    a = transfer (c, a)
  end function scm_to_c_ptr

  subroutine assign_to_c_ptr (a, b)
    !
    ! SCM exact = c_funptr
    !
    use iso_c_binding, only: c_ptr
    implicit none
    type (c_ptr), intent (out) :: a
    type (scm_t), intent (in) :: b
    ! *** end of interface ***

    a = scm_to_c_ptr (b)
  end subroutine assign_to_c_ptr

  function scm_from_c_funptr (b) result(a)
    !
    ! Conver type (c_funptr) to SCM pointer object.
    !
    use iso_c_binding, only: c_funptr, c_sizeof, c_null_funptr
    implicit none
    type (c_funptr), intent (in) :: b
    type (scm_t) :: a
    ! *** end of interface ***

    type (c_ptr) :: c

    ! Use F2008  c_sizeof().  It is however not  availbale in Gfortran
    ! 4.3 on Debian Lenny.  On the other hand sizeof() is available as
    ! an extension  in both gfortran and  ifort.  It is  used here for
    ! sanity check only.
    if (c_sizeof (b) /= c_sizeof (c)) then
       print *, "scm_from_c_funptr: size mismatch:", c_sizeof (b), "/=", c_sizeof (c)
       stop 1
    endif

    ! FIXME: potentially non-portable  coversion of a function pointer
    ! to a data pointer:
    c = transfer (b, c)

    ! Guile 2.0 function here:
    a = scm_from_pointer (c, c_null_funptr)
  end function scm_from_c_funptr

  subroutine assign_from_c_funptr (a, b)
    !
    ! SCM pointer object = type (c_funptr)
    !
    use iso_c_binding, only: c_funptr
    implicit none
    type (scm_t), intent (out) :: a
    type (c_funptr), intent (in) :: b
    ! *** end of interface ***

    a = scm_from_c_funptr (b)
  end subroutine assign_from_c_funptr

end module scm
