module problem_mod
    use, intrinsic :: iso_c_binding
    implicit none

    procedure(setdim_abstract), pointer :: setdim
    procedure(startp_abstract), pointer :: startp
    procedure(setbounds_abstract), pointer :: setbounds
    procedure(functs_abstract), pointer :: functs
    procedure(fconstriq_abstract), pointer :: fconstriq

    interface
        subroutine setdim_abstract(n, m, q) bind(C)
            use, intrinsic :: iso_c_binding
            implicit none

            ! take references to integer values so that we can modify them
            integer(c_int) :: n  ! variables
            integer(c_int) :: m  ! constraints
            integer(c_int) :: q  ! objectives
        end subroutine

        subroutine startp_abstract(n, x) bind(C)
            use, intrinsic :: iso_c_binding
            implicit none
            
            ! take `n` by value so that we don't change it
            integer(c_int), value :: n
            ! but `x` is a pointer to pre-allocated buffer which we set
            real(c_double) :: x(n)
        end subroutine

        subroutine setbounds_abstract(n, lb, ub) bind(C)
            use, intrinsic :: iso_c_binding
            implicit none

            integer(c_int), value :: n
            real(c_double) :: lb(n)
            real(c_double) :: ub(n)
        end subroutine

        subroutine functs_abstract(n, x, q, f) bind(C)
            use, intrinsic :: iso_c_binding
            implicit none
            
            integer(c_int), value :: n  ! vars
            integer(c_int), value :: q  ! objectives
            real(c_double) :: x(n)  ! input
            real(c_double) :: f(q)  ! output
        end subroutine

        subroutine fconstriq_abstract(n, m, x, ciq) bind(C)
            use, intrinsic :: iso_c_binding
            implicit none
            
            integer(c_int), value :: n  ! vars
            integer(c_int), value :: m  ! constraints
            real(c_double) :: x(n)  ! input
            real(c_double) :: ciq(m)  ! output
        end subroutine
 
    end interface

    contains

        subroutine set_setdim_ptr(ptr_to_be_wrapped) bind(C, name="set_setdim_ptr")
            use, intrinsic :: iso_c_binding
            implicit none
            
            type(c_funptr), value :: ptr_to_be_wrapped

            call c_f_procpointer(ptr_to_be_wrapped, setdim)
        end subroutine

        subroutine set_startp_ptr(ptr_to_be_wrapped) bind(C, name="set_startp_ptr")
            use, intrinsic :: iso_c_binding
            implicit none
            
            type(c_funptr), value :: ptr_to_be_wrapped

            call c_f_procpointer(ptr_to_be_wrapped, startp)
        end subroutine

        subroutine set_setbounds_ptr(ptr_to_be_wrapped) bind(C, name="set_setbounds_ptr")
            use, intrinsic :: iso_c_binding
            implicit none
            
            type(c_funptr), value :: ptr_to_be_wrapped

            call c_f_procpointer(ptr_to_be_wrapped, setbounds)
        end subroutine

        subroutine set_functs_ptr(ptr_to_be_wrapped) bind(C, name="set_functs_ptr")
            use, intrinsic :: iso_c_binding
            implicit none
            
            type(c_funptr), value :: ptr_to_be_wrapped

            call c_f_procpointer(ptr_to_be_wrapped, functs)
        end subroutine

        subroutine set_fconstriq_ptr(ptr_to_be_wrapped) bind(C, name="set_fconstriq_ptr")
            use, intrinsic :: iso_c_binding
            implicit none
            
            type(c_funptr), value :: ptr_to_be_wrapped

            call c_f_procpointer(ptr_to_be_wrapped, fconstriq)
        end subroutine

end module problem_mod