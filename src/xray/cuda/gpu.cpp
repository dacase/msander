//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------
#include <stdio.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>
#include <string.h>
#include <algorithm>

/*  #include "gpuContext.h"  */
#include "gpu.h"

using namespace std;

/*  #include "matrix.h"  */
/*  #include "bondRemap.h" */

//---------------------------------------------------------------------------------------------
// gpu_startup_: establish the identity of the GPU to use
extern "C" void gpu_startup_(void)
{
  gpu = new _gpuContext;
  PRINTMETHOD("gpu_startup"); 
}

//---------------------------------------------------------------------------------------------
// gpu_set_device_: set the device to use in the calculation
//
// Arguments:
//   device:   ID number of the GPU (CUDA_VISIBLE_DEVICES environment variable)
//---------------------------------------------------------------------------------------------
extern "C" void gpu_set_device_(int* device)
{
  PRINTMETHOD("gpu_set_device"); 
  gpu->gpu_device_id = *device;
}

//---------------------------------------------------------------------------------------------
// gpu_init_: master function for initializing the device.  This calls setup for all types of
//            simulations and detects essential information about the GPU, like Streaming
//            Multiprocessor (SM) version.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_init_(void)
{
  PRINTMETHOD("gpu_init");
  int LRFSize = 0;
  int SMCount = 0;
  int SMMajor = 0;
  int SMMinor = 0;
    
  int device = -1;
  int gpuCount = 0;
  cudaError_t status;
  cudaDeviceProp deviceProp;
  status = cudaGetDeviceCount(&gpuCount);
  RTERROR(status, "cudaGetDeviceCount failed");
  if (gpuCount == 0) {
    printf("No CUDA-capable devices found, exiting.\n");
    cudaDeviceReset();
    exit(-1);
  }

  // Activate zero-copy
  cudaSetDeviceFlags(cudaDeviceMapHost);

  // If the device id is -1 this means it was left at the default so we 
  // or the user specified -1 on the command line meaning autoselect GPU.
  // choose CUDA device with the most memory (flesh out later for multi-gpu)
  // otherwise we use the device id specified on the command line.
  if (gpu->gpu_device_id == -1) {
	  
    // Generate list of compatible GPUs scored by GPU revision first and total memory second
    int* pGPUList = new int[gpuCount];
    unsigned int* pGPUScore = new unsigned int[gpuCount];
    int gpus = 0;
    for (int i = 0; i < gpuCount; i++) {
      cudaGetDeviceProperties(&deviceProp, i);
      if (((deviceProp.major >= 3) || 
           ((deviceProp.major == 1) && (deviceProp.minor == 3)))) {
        pGPUList[gpus] = i;
        pGPUScore[gpus] = (deviceProp.major << 24) + (deviceProp.totalGlobalMem >> 20);
        gpus += (deviceProp.major >= 3);
      }
    }
    if (gpus == 0) {
      printf("Error searching for compatible GPU");
      exit(-1);
    }
    
    // Select best GPU according to score
    if (gpus > 0) {

      // Bubble sort (go bubblesort!) device list by score
      bool done = true;
      do {
        done = true;
        for (int i = 0; i < gpus - 1; i++) {
          if (pGPUScore[i] < pGPUScore[i + 1]) {
            done = false;
            int gpu = pGPUList[i];
            unsigned int score = pGPUScore[i];
            pGPUList[i] = pGPUList[i + 1];
            pGPUScore[i] = pGPUScore[i + 1];
            pGPUList[i + 1] = gpu;
            pGPUScore[i + 1] = score;
          }
        }
      } while (!done);
    }
            
    // Let CUDA select any device from this list
    status = cudaSetValidDevices(pGPUList, gpus);
    delete[] pGPUList;
    delete[] pGPUScore;
    RTERROR(status, "Error searching for compatible GPU");
    // Trick driver into creating a context on an available and valid GPU
    status = cudaFree(0);
    RTERROR(status, "Error selecting compatible GPU");

    // Get device
    status = cudaGetDevice(&device);
    RTERROR(status, "Error fetching current GPU");
  } else {
    cudaGetDeviceProperties(&deviceProp, gpu->gpu_device_id);
    if (deviceProp.major >= 3) {
      device = gpu->gpu_device_id;
    }
    else {
      printf("Selected GPU lacks SM 3.0 or better support, exiting.\n");
      cudaDeviceReset();
      exit(-1);
    }
  }

  if (device == -1) {
    printf("No double-precision capable gpu located, exiting.\n");
    cudaDeviceReset();
    exit(-1);
  }
    
  // Finally set CUDA device
  status = cudaSetDevice(device); 
  RTERROR(status, "Error setting CUDA device");  
  cudaDeviceSynchronize();

  // Determine kernel call configuration and grab desired additional GPU properties
  cudaGetDeviceProperties(&deviceProp, device);
  //gpu->bECCSupport = deviceProp.ECCEnabled || deviceProp.tccDriver ||
  //                   (strcasestr(deviceProp.name, "tesla") != NULL);
  // TL: comment out "telsa"--ECC support can be on non-tesla GPUs
  gpu->bECCSupport = deviceProp.ECCEnabled || deviceProp.tccDriver ;
  gpu->bCanMapHostMemory = deviceProp.canMapHostMemory;
  gpu->totalMemory = deviceProp.totalGlobalMem;

#ifdef GVERBOSE
  double memsize = (double)deviceProp.totalGlobalMem / (1024.0 * 1024.0);
  printf("Using GPU %d, %s, SM %d.%d, %.1f MBytes of memory\n", device, deviceProp.name,
         deviceProp.major, deviceProp.minor, memsize);
#endif

  // Store GPU Device ID for later use
  gpu->gpu_device_id = device;
  gpu->blocks = deviceProp.multiProcessorCount;
  gpu->GBBornRadiiBlocks = gpu->blocks;
  gpu->GBNonbondEnergy1Blocks = gpu->blocks;
  gpu->GBNonbondEnergy2Blocks = gpu->blocks;
  gpu->BNLBlocks = gpu->blocks;

  // Determine SM version
  unsigned int blocksPerSM;
  gpu->sm_version = SM_3X;
  gpu->threadsPerBlock = THREADS_PER_BLOCK;
  gpu->NLCalculateOffsetsThreadsPerBlock = NLCALCULATE_OFFSETS_THREADS_PER_BLOCK;

  return;
}

//---------------------------------------------------------------------------------------------
// gpu_shutdown_: function for shutting down GPUs and deleting memory allocated on the host
//                associated with them.
//---------------------------------------------------------------------------------------------
extern "C" void gpu_shutdown_(void)
{
  PRINTMETHOD("gpu_shutdown");
  delete gpu;
  cudaDeviceReset();
  return;
}

//---------------------------------------------------------------------------------------------
// gpu_get_device_info_: obtains information on the GPUs installed and what they are called.
//
// Arguments:
//   gpu_dev_count:  the number of GPU devices
//   gpu_dev_id:     identifications of the GPU device for which to obtain information
//   gpu_dev_mem:    total device memory of this GPU
//   gpu_num_proc:   the number of CUDA cores in this GPU
//   gpu_core_freq:  cycles/second of each core (i.e. 1.4GHz)
//   gpu_name:       the name of this GPU device
//   name_len:       length of the device name
//---------------------------------------------------------------------------------------------
extern "C" void gpu_get_device_info_(int* gpu_dev_count, int* gpu_dev_id, int* gpu_dev_mem,
                                     int* gpu_num_proc, double* gpu_core_freq, char* gpu_name,
                                     int* name_len)
{
  PRINTMETHOD("gpu_get_device_info");
  cudaError_t status;
  cudaDeviceProp deviceProp;
  size_t device_mem;

  status = cudaGetDeviceCount(gpu_dev_count);
  RTERROR(status, "cudaGetDeviceCount failed");

  *gpu_dev_id = gpu->gpu_device_id;
  cudaGetDeviceProperties(&deviceProp, *gpu_dev_id);

  device_mem = deviceProp.totalGlobalMem/(1024*1024);
  *gpu_dev_mem = (int )device_mem;

  *gpu_num_proc = (int )deviceProp.multiProcessorCount;
  *gpu_core_freq = (double )(deviceProp.clockRate * 1e-6f);

  strcpy(gpu_name,deviceProp.name);
  *name_len=strlen(deviceProp.name);
}

//---------------------------------------------------------------------------------------------
// gpu_check_titan_: function to check Amber compatibility with NVIDIA GTX TITAN or 780
//---------------------------------------------------------------------------------------------
extern "C" void gpu_check_titan_()
{
  PRINTMETHOD("gpu_check_titan");

  // Variables to get GPU name
  int gpu_dev_id = 0;
  cudaDeviceProp deviceProp;
  char* gpu_name;
  gpu_name = (char *) malloc(128 * sizeof(char));

  // Variables to get driver version
  float driver_version = 0.0f;
  FILE* nvidia_version_file;
  char* version_line;
  version_line = (char *) malloc(128 * sizeof(char));

  // Get GPU name
  gpu_dev_id = gpu->gpu_device_id;
  cudaGetDeviceProperties(&deviceProp, gpu_dev_id);
  strcpy(gpu_name,deviceProp.name);

  // Check for GeForce GTX TITAN or 780
  if (strcmp(gpu_name, "GeForce GTX TITAN") == 0 || strcmp(gpu_name, "GeForce GTX 780") == 0) {

    // Get driver version from /proc/driver/nvidia/version
    // Note that this may be linux-specific
    nvidia_version_file = fopen("/proc/driver/nvidia/version", "r");
    if (nvidia_version_file != NULL) {
      char* throwaway = fgets(version_line, 128, nvidia_version_file);
      if (version_line != NULL) {
        sscanf(version_line, "%*s %*s %*s %*s %*s %*s %*s %f", &driver_version);
      }
    }
    fclose(nvidia_version_file);

    // Check for NVIDIA driver version 325.15.
    // Note on NVIDIA drivers: driver_branch.build_in_branch (not necessarily ordered in time).
    // This works for now, but because of the way NVIDIA drivers are organized, 
    // this may not always be correct.
    if (driver_version < 325.15f) {
      printf("Error: NVIDIA graphics card and NVIDIA driver version "
             "incompatible with Amber\n");
      printf("GeForce GTX TITAN or GeForce GTX 780 need at least driver "
             "version 325.15 for Amber\n");
      exit(-1);
    }
  }
  free(gpu_name);
  free(version_line);
}

//---------------------------------------------------------------------------------------------
// gpu_get_memory_info_: returns KB of memory in use on CPU and GPU
//
// Arguments:
//   gpumemory: KB of memory on the GPU device
//   cpumemory: KB of memory in the host RAM
//---------------------------------------------------------------------------------------------
extern "C" void gpu_get_memory_info_(int* gpumemory, int* cpumemory)
{
  // The "ll" suffix denotes that the integer be formulated as a long long int
  //*gpumemory  = (int)(gpu->totalGPUMemory / 1024ll);
  //*cpumemory  = (int)(gpu->totalCPUMemory / 1024ll);

  *gpumemory = gpuMemoryInfo::Instance().totalGPUMemory / 1024ll;
  *cpumemory = gpuMemoryInfo::Instance().totalCPUMemory / 1024ll;
  return;
}

//---------------------------------------------------------------------------------------------
// gpu_upload_crd_: upload coordinates to the device.  This routine is called by all Fortran
//                  routines seeking to interface with the CUDA code.
//
// Arguments:
//    atm_crd:  the atomic coordinates
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_crd_(double atm_crd[][3])
{
  PRINTMETHOD("gpu_upload_crd");
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    gpu->pbImage->Download();
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    double *pCrd = gpu->pbImage->_pSysData;
    if (gpu->sim.pImageX != gpu->pbImage->_pDevData) {
      pCrd = gpu->pbImage->_pSysData + gpu->sim.stride3;
    }
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      pCrd[i1] = atm_crd[i][0];
      pCrd[i1 + gpu->sim.stride] = atm_crd[i][1];
      pCrd[i1 + gpu->sim.stride2] = atm_crd[i][2];
    }
    gpu->pbImage->Upload();
  }
  else {
    for (int i = 0; i < gpu->sim.atoms; i++) {
      gpu->pbAtom->_pSysData[i] = atm_crd[i][0];
      gpu->pbAtom->_pSysData[i + gpu->sim.stride] = atm_crd[i][1];
      gpu->pbAtom->_pSysData[i + gpu->sim.stride2] = atm_crd[i][2];
      gpu->pbAtomXYSP->_pSysData[i].x = atm_crd[i][0];
      gpu->pbAtomXYSP->_pSysData[i].y = atm_crd[i][1];
      gpu->pbAtomZSP->_pSysData[i] = atm_crd[i][2];
    }
    for (int i = gpu->sim.atoms; i < gpu->sim.stride; i++) {
      gpu->pbAtom->_pSysData[i] = 9999990000.0 + i*2000.0;
      gpu->pbAtom->_pSysData[i + gpu->sim.stride] = 9999990000.0 + i*2000.0;
      gpu->pbAtom->_pSysData[i + gpu->sim.stride2] = 9999990000.0 + i*2000.0;
      gpu->pbAtomXYSP->_pSysData[i].x = 9999990000.0f + i*2000.0;
      gpu->pbAtomXYSP->_pSysData[i].y = 9999990000.0f + i*2000.0;
      gpu->pbAtomZSP->_pSysData[i] = 9999990000.0f + i*2000.0;
    }
    gpu->pbAtom->Upload();
    gpu->pbAtomXYSP->Upload();
    gpu->pbAtomZSP->Upload();   
  }
}


//---------------------------------------------------------------------------------------------
// gpu_download_crd_: download coordinates from the device.  This routine is called by all
//                    routines seeking to interface with the GPU code.
//
// Arguments:
//   atm_crd:  the atomic coordinates
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_crd_(double atm_crd[][3])
{
  PRINTMETHOD("gpu_download_crd");
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    gpu->pbImage->Download();
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    double *pCrd = gpu->pbImage->_pSysData;
    if (gpu->sim.pImageX != gpu->pbImage->_pDevData) {
      pCrd = gpu->pbImage->_pSysData + gpu->sim.stride3;
    }
    int i;
    for (i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      atm_crd[i][0] = pCrd[i1];
      atm_crd[i][1] = pCrd[i1 + gpu->sim.stride];
      atm_crd[i][2] = pCrd[i1 + gpu->sim.stride2];
    }
  }
  else {
    gpu->pbAtom->Download();
    for (int i = 0; i < gpu->sim.atoms; i++) {
      atm_crd[i][0] = gpu->pbAtom->_pSysData[i];
      atm_crd[i][1] = gpu->pbAtom->_pSysData[i + gpu->sim.stride];
      atm_crd[i][2] = gpu->pbAtom->_pSysData[i + gpu->sim.stride2];
    }
  }     
}

//---------------------------------------------------------------------------------------------
// gpu_upload_frc_: upload the forces on all atoms to the device.  While it is anathema to
//                  transfer so much data between the host and the device, this is needed in
//                  Monte Carlo barostating (even when a change in box volume succeeds, the
//                  old forces must still be reloaded) and in REMD (when a move fails, the old
//                  forces again must be reloaded).
//
// Arguments:
//   atm_frc:    the forces on all atoms
//---------------------------------------------------------------------------------------------
extern "C" void gpu_upload_frc_(double atm_frc[][3])
{
  PRINTMETHOD("gpu_upload_frc");
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
    PMEAccumulator *pForce = gpu->pbForceAccumulator->_pSysData;
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1 = pImageAtomLookup[i];
      pForce[i1] = (PMEAccumulator)(FORCESCALE * atm_frc[i][0]);
      pForce[i1 + gpu->sim.stride] = (PMEAccumulator)(FORCESCALE * atm_frc[i][1]);
      pForce[i1 + gpu->sim.stride2] = (PMEAccumulator)(FORCESCALE * atm_frc[i][2]);
    }
    gpu->pbForceAccumulator->Upload();
  }
  else {
    for (int i = 0; i < gpu->sim.atoms; i++) {
      gpu->pbForceAccumulator->_pSysData[i] = (PMEAccumulator)(FORCESCALE * atm_frc[i][0]);
      gpu->pbForceAccumulator->_pSysData[i + gpu->sim.stride] = 
        (PMEAccumulator)(FORCESCALE * atm_frc[i][1]);
      gpu->pbForceAccumulator->_pSysData[i + gpu->sim.stride2] =
        (PMEAccumulator)(FORCESCALE * atm_frc[i][2]);
    }
    gpu->pbForceAccumulator->Upload();
  }
}

//---------------------------------------------------------------------------------------------
// gpu_download_frc_: download forces from the GPU device.  Generally, this is no longer used
//                    when updating forces on the GPU with contributions calculated on the
//                    host--those contributions are simply uploaded.  However, it is called
//                    in REMD and MC barostating, when failed moves mean that system must
//                    revert to its old state and use the old forces.
//
// Arguments:
//   atm_frc:   the forces on all atoms
//---------------------------------------------------------------------------------------------
extern "C" void gpu_download_frc_(double atm_frc[][3])
{
  PRINTMETHOD("gpu_download_frc");
  gpu->pbForceAccumulator->Download();
  PMEAccumulator *pForce  = gpu->pbForceAccumulator->_pSysData;
  if (gpu->bNeighborList && (gpu->pbImageIndex != NULL)) {
    if (gpu->bNewNeighborList) {
      gpu->pbImageIndex->Download();
      gpu->bNewNeighborList = false;
    }
    unsigned int* pImageAtomLookup = &(gpu->pbImageIndex->_pSysData[gpu->sim.imageStride * 2]);
#ifndef MPI
    if ((gpu->sim.ntp > 0) && (gpu->sim.barostat == 1)) {
      PMEAccumulator *pNBForce  = gpu->pbForceAccumulator->_pSysData + gpu->sim.stride3;
      for (int i = 0; i < gpu->sim.atoms; i++) {
        int i1 = pImageAtomLookup[i];
        atm_frc[i][0] = (double)(pForce[i1] + pNBForce[i1]) * (double)ONEOVERFORCESCALE;
        atm_frc[i][1] = (double)(pForce[i1 + gpu->sim.stride] +
                                 pNBForce[i1 + gpu->sim.stride]) * (double)ONEOVERFORCESCALE;
        atm_frc[i][2] = (double)(pForce[i1 + gpu->sim.stride2] +
                                 pNBForce[i1 + gpu->sim.stride2]) * (double)ONEOVERFORCESCALE;
#ifdef use_DPFP
        if (gpu->imin != 0)
        {
          atm_frc[i][0]    += (double)pNBForce[i1 + gpu->sim.stride3];
          atm_frc[i][1]    += (double)pNBForce[i1 + gpu->sim.stride4];
          atm_frc[i][2]    += (double)pNBForce[i1 + gpu->sim.stride5];        
        }
#endif
      }
    }
    else {
#endif
    for (int i = 0; i < gpu->sim.atoms; i++) {
      int i1              = pImageAtomLookup[i];
      atm_frc[i][0]       = (double)pForce[i1] * (double)ONEOVERFORCESCALE;
      atm_frc[i][1]       = (double)pForce[i1 + gpu->sim.stride] * (double)ONEOVERFORCESCALE;
      atm_frc[i][2]       = (double)pForce[i1 + gpu->sim.stride2] * (double)ONEOVERFORCESCALE;
#ifdef use_DPFP
      if (gpu->imin != 0)
      {
        atm_frc[i][0]    += (double)pForce[i1 + gpu->sim.stride3];
        atm_frc[i][1]    += (double)pForce[i1 + gpu->sim.stride4];
        atm_frc[i][2]    += (double)pForce[i1 + gpu->sim.stride5];        
      }
#endif
    }
#ifndef MPI
    }
#endif
  }
  else {
    for (int i = 0; i < gpu->sim.atoms; i++) {
      atm_frc[i][0] = (double)pForce[i] * (double)ONEOVERFORCESCALE;
      atm_frc[i][1] = (double)pForce[i + gpu->sim.stride] * (double)ONEOVERFORCESCALE;
      atm_frc[i][2] = (double)pForce[i + gpu->sim.stride2] * (double)ONEOVERFORCESCALE;
#ifdef use_DPFP
      if (gpu->imin != 0)
      {
        atm_frc[i][0]    += (double)pForce[i + gpu->sim.stride3];
        atm_frc[i][1]    += (double)pForce[i + gpu->sim.stride4];
        atm_frc[i][2]    += (double)pForce[i + gpu->sim.stride5];        
      }
#endif
    } 
  }
}


