//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

#ifndef __GPU_H__
#define __GPU_H__

#include "gpuContext.h"
#include "gputypes.h"
// F95 interface
extern "C" void gpu_init_(void);
extern "C" void gpu_shutdown_(void);
extern "C" void gpu_upload_crd_(double atm_crd[][3]);
extern "C" void gpu_download_crd_(double atm_crd[][3]);
extern "C" void gpu_upload_frc_(double atm_frc[][3]);
#endif
