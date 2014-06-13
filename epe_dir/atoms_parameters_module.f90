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
module atoms_parameters_module

  use type_module
  use epecom_module

  implicit none
  save
  private

  type, public::par_item
     integer(kind=i4_kind)::nlth
     character(len=3),dimension(3):: atom_name
     real(kind=r8_kind),dimension(3,3)::bi,roi,ci,di,cutoff
     real(kind=r8_kind),dimension(3)::q_zi,pki,q_ni,q_si
  end type par_item

  type, public:: epe_par 
     type(par_item) :: item
     type(epe_par),pointer::next
  end type epe_par

  type(par_item),public :: parv
  type(epe_par),public, pointer ::p,pp,top

  integer(kind=i4_kind),parameter,public :: nlth=3
  character(len=3),public :: atom_name(nlth)=(/'O  ','Alo','Alt'/)

  real(kind=r8_kind),public :: bi(1:nlth,1:nlth)=reshape( &
                                            (/22764.3,1114.9,1012.6, &
                                             1114.9,0.00000,0.00000, &
                                             1012.6,0.00000,0.00000/),(/3,3/) )

  real(kind=r8_kind),public :: Roi(nlth,nlth)=reshape( &
                                            (/0.14900,0.3118000,0.3118000, &
                                              0.3118000,1.000000,1.000000, &
                                              0.3118000,1.000000,1.000000/),&
                                              (/3,3/))

  real(kind=r8_kind),public :: ci(nlth,nlth)=reshape( &
                                            (/20.3700,0.00000,0.00000, &
                                             0.00000,0.00000,0.00000, &
                                             0.00000,0.00000,0.00000/),&
                                                (/3,3/))

  real(kind=r8_kind),public :: di(nlth,nlth)=reshape( &
                                            (/0.00000,0.00000,0.00000, &
                                             0.00000,0.00000,0.00000, &
                                             0.00000,0.00000,0.00000/),&
                                                (/3,3/))
  real(kind=r8_kind),public :: cutoff(nlth,nlth)=reshape( &
            (/3.5,2.5,2.5, 2.5,2.5,2.5, 2.5,2.5,2.5/), (/3,3/))

real(kind=r8_kind),public :: q_zi(nlth)=(/-2.0,3.0,3.0/)
real(kind=r8_kind),public :: pki(nlth)=(/27.29,99999.0,99999.0/)
real(kind=r8_kind),public :: q_ni(nlth)=(/-2.0,3.0,3.0/)
real(kind=r8_kind),public :: q_si(nlth)=(/-2.21,3.0,3.0/)
        public init_epe_par,set_default_parv
contains
        subroutine set_default_parv
        parv=par_item(3,(/ 'Z  ','Z  ','Z  ' /),&
        reshape((/0.,0.,0.,0.,0.,0., 0., 0., 0./),(/3,3/)),&
        reshape((/1.,1.,1.,1.,1.,1., 1., 1., 1./),(/3,3/)),&
        reshape((/0.,0.,0.,0.,0.,0., 0., 0., 0./),(/3,3/)),&
        reshape((/0.,0.,0.,0.,0.,0., 0., 0., 0./),(/3,3/)),&
        reshape((/3.4,2.4,2.4,2.4,2.4,2.4, 2.4, 2.4, 2.4/),(/3,3/)),&
        (/ 0.0,0.0,0.0/),(/ 1.0,1.0,1.0/), &
        (/ 0.0,0.0,0.0/),(/ 0.0,0.0,0.0/))
        end subroutine set_default_parv
  subroutine init_epe_par
    nullify(pp)
    allocate(p)
    parv=par_item(3,(/'O  ','Alo','Alt'/),&
         reshape((/22764.3,1114.9,1012.6, &
                   1114.9,0.00000,0.00000, &
                   1012.6,0.00000,0.00000/),(/3,3/)),&
         reshape((/0.14900,0.3118000,0.3118000, &
                   0.3118000,1.000000,1.000000, &
                   0.3118000,1.000000,1.000000/),(/3,3/)),&
         reshape((/20.3700,0.00000,0.00000, &
                   0.00000,0.00000,0.00000, &
                   0.00000,0.00000,0.00000/),(/3,3/)),&
         reshape((/0.00000,0.00000,0.00000, &
                   0.00000,0.00000,0.00000, &
                   0.00000,0.00000,0.00000/),(/3,3/)),&
         reshape((/3.40000,2.40000,2.40000, &
                   2.40000,2.40000,2.40000, &
                   2.40000,2.40000,2.40000/),(/3,3/)),&
                (/-2.0,3.0,3.0/),(/27.29,99999.0,99999.0/),&
                (/-2.0,3.0,3.0/),(/-2.21,3.0,3.0/))
    p=epe_par(parv,pp)
    pp=>p
    allocate(p)
    parv=par_item(3,(/'O  ','Alo','Alt'/),&
         reshape((/22764.3,1114.9,1012.6, &
                   1114.9,0.00000,0.00000, &
                   1012.6,0.00000,0.00000/),(/3,3/)),&
         reshape((/0.14900,0.3118000,0.3118000, &
                   0.3118000,1.000000,1.000000, &
                   0.3118000,1.000000,1.000000/),(/3,3/)),&
         reshape((/20.3700,0.00000,0.00000, &
                   0.00000,0.00000,0.00000, &
                   0.00000,0.00000,0.00000/),(/3,3/)),&
         reshape((/0.00000,0.00000,0.00000, &
                   0.00000,0.00000,0.00000, &
                   0.00000,0.00000,0.00000/),(/3,3/)),&
         reshape((/3.40000,2.40000,2.40000, &
                   2.40000,2.40000,2.40000, &
                   2.40000,2.40000,2.40000/),(/3,3/)),&
                (/-2.0,3.0,3.0/),(/27.29,99999.0,99999.0/),&
                (/-2.0,3.0,3.0/),(/-2.21,3.0,3.0/))
    p=epe_par(parv,pp)
    top=>p
  end subroutine init_epe_par
end module atoms_parameters_module
