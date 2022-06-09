module xray_bulk_mask_module
  
  use xray_pure_utils, only: real_kind

#ifdef CUDA
  use xray_bulk_mask_impl_gpu_module
#else
  use xray_bulk_mask_impl_cpu_module
#endif
  
  implicit none
  private
  
  public :: update_f_bulk
  public :: finalize
  public :: init
 
end module xray_bulk_mask_module