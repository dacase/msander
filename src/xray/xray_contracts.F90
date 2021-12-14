module xray_contracts_module

  implicit none
  
contains

pure subroutine check_requirement(condition, message)
  implicit none
  logical, intent(in) :: condition
  character(len = *), optional, intent(in) :: message

  if (.not. condition) then
    if (present(message)) then
      error stop "Unmet requirement: " // message
    else
      error stop "Unmet requirement."
    end if
  endif

end subroutine check_requirement

pure subroutine check_assertion(condition, message)
  implicit none
  logical, intent(in) :: condition
  character(len = *), optional, intent(in) :: message

#ifndef NDEBUG
  if (.not. condition) then
    if (present(message)) then
      error stop "Assertion failed: " // message
    else
      error stop "Assertion failed."
    end if
  endif
#endif

end subroutine check_assertion

pure subroutine check_precondition(condition, message)
  implicit none
  logical, intent(in) :: condition
  character(len=*), optional, intent(in) :: message

#ifndef NDEBUG
  if (.not. condition) then
    if (present(message)) then
      error stop "Precondition failed: " // message
    else
      error stop "Precondition failed."
    end if
  endif
#endif

end subroutine check_precondition

pure subroutine check_postcondition(condition, message)
  implicit none
  logical, intent(in) :: condition
  character(len=*), optional, intent(in) :: message

#ifndef NDEBUG
  if (.not. condition) then
    if (present(message)) then
      error stop "Postcondition failed: " // message
    else
      error stop "Postcondition failed."
    end if
  endif
#endif

end subroutine check_postcondition

end module xray_contracts_module