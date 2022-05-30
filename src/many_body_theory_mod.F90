!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!                                                                                                                                  !
!  Copyright (2003). See the README file in the top-level directory.                                                               !
!  This software is distributed with the GNU General Public License.                                                               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

module many_body_theory_mod

   use config_sc_mod
   use kind_mod

   implicit none ; private

   !*

   public :: many_body_theory

   !*

   interface many_body_theory
      module procedure many_body_theory_
   end interface

contains

   function many_body_theory_(cfg) result(changed)

      type(config_sc_obj), intent(inout) :: cfg
      logical :: changed

      changed = .false.

   end function

end module many_body_theory_mod
