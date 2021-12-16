#ifndef _GPU_BUFFER
#define _GPU_BUFFER

#include "gpuMemoryInfo.h"

//---------------------------------------------------------------------------------------------
// Template for all the various GpuBuffer types, of which _gpuContext uses so many
//---------------------------------------------------------------------------------------------
template <typename T> class GpuBuffer {

public:

  size_t _length;                  // The overall length of the array to allocate (total memory
                                   //   allocation will depend on the size of the data type)
  bool _bSysMem;                   // Flag to indicate that host (CPU) memory should be
                                   //   allocated to "shadow" what happens on the device (GPU)
  bool _bPinned;                   // Flag to indicate that the array held in host RAM will be
                                   //   page locked, or "pinned." Such data is more rapidly
                                   //   accessible for transfer to the GPU, as the GPU cannot
                                   //   directly access pageable host memory.  If allocated as
                                   //   pinned memory, _pSysData becomes a pre-allocated
                                   //   staging area for transferring memory to the GPU.
  T* _pSysData;                    // System data, held in CPU RAM
  T* _pDevData;                    // Device data, can be shadowed by _pSysData

  // Overloading, to make three constructor methods and allow size_t, ints and unsigned ints to
  // serve as the length argument for allocating a GpuBuffer struct.  In all cases, the
  // default behavior is to allocate memory on the host to shadow the device memory, and to
  // have that host memory be allocated as pageable (NOT pinned).
  GpuBuffer(int length, bool bSysMem = true, bool bPinned = false);
  GpuBuffer(unsigned int length, bool bSysMem = true, bool bPinned = false);
  GpuBuffer(size_t length, bool bSysMem = true, bool bPinned = false);

  // Having a virtual destructor ensures that there are no
  // memory leaks for derived GpuBuffer objects.
  virtual ~GpuBuffer();

  // Methods
  void Allocate();
  void Deallocate();
  void Upload(T* pBuff = NULL);
  void Download(T* pBuff = NULL);
};

//---------------------------------------------------------------------------------------------
// First of three constructors for the GpuBuffer data type
//---------------------------------------------------------------------------------------------
template <typename T> GpuBuffer<T>::GpuBuffer(unsigned int length, bool bSysMem,
                                              bool bPinned) :
_length(length), _bSysMem(bSysMem), _bPinned(bPinned), _pSysData(NULL), _pDevData(NULL)
{
  Allocate();
}

//---------------------------------------------------------------------------------------------
// Second of three constructors for the GpuBuffer data type
//---------------------------------------------------------------------------------------------
template <typename T> GpuBuffer<T>::GpuBuffer(int length, bool bSysMem, bool bPinned) :
_length(length), _bSysMem(bSysMem), _bPinned(bPinned), _pSysData(NULL), _pDevData(NULL)
{
  Allocate();
}

//---------------------------------------------------------------------------------------------
// Second of three constructors for the GpuBuffer data type
//---------------------------------------------------------------------------------------------
template <typename T> GpuBuffer<T>::GpuBuffer(size_t length, bool bSysMem, bool bPinned) :
_length(length), _bSysMem(bSysMem), _bPinned(bPinned), _pSysData(NULL), _pDevData(NULL)
{
  Allocate();
}

//---------------------------------------------------------------------------------------------
// Destructor for the GpuBuffer data type
//---------------------------------------------------------------------------------------------
template <typename T> GpuBuffer<T>::~GpuBuffer()
{
  Deallocate();
}

//---------------------------------------------------------------------------------------------
// Allocate: allocates memory on the host and on the device for a GpuBuffer of type T
//---------------------------------------------------------------------------------------------
template <typename T> void GpuBuffer<T>::Allocate()
{
  if (_length == 0) {
    return;
  }
  size_t sizeT = sizeof(T);
  
  cudaError_t status;

  // If pinned memory is to be allocated, pinned memory is something that lives in CPU RAM,
  // so that implies host memory WILL be allocated to shadow the device, and that it will
  // be un-pageable, unlike usual host memory allocations.  Pinned memory on the host will
  // correspond to a block of memory on the device, so allocating memory in this way implies
  // roping off an equivalent amount of memory on the device, to which we simply need to get
  // a pointer.
  if (_bPinned) {
    status = cudaHostAlloc((void **)&_pSysData, _length * sizeT, cudaHostAllocMapped);
    RTERROR(status, "cudaHostalloc GpuBuffer::Allocate failed");
    gpuMemoryInfo::Instance().totalCPUMemory += _length * sizeT;
    gpuMemoryInfo::Instance().totalGPUMemory += _length * sizeT;
    status = cudaHostGetDevicePointer((void **)&_pDevData, (void *)_pSysData, 0);
    RTERROR(status, "cudaGetDevicePointer GpuBuffer::failed to get device pointer");
    memset(_pSysData, 0, _length * sizeT);
  }
  else {

    // Check to see whether (pageable) host memory has been requested to shadow the device.
    // This memory will be allocated by the more familiar C++ new operator rather than
    // cudaHostAlloc().  An equivalent amount of memory will be allocated on the device,
    // not mapped to memory on the host.
    if (_bSysMem) {
      _pSysData = new T[_length];
      gpuMemoryInfo::Instance().totalCPUMemory += _length * sizeT;
      memset(_pSysData, 0, _length * sizeT);
    }
    status = cudaMalloc((void **)&_pDevData, _length * sizeT);
    gpuMemoryInfo::Instance().totalGPUMemory += _length * sizeT;
    RTERROR(status, "cudaMalloc GpuBuffer::Allocate failed");
    status = cudaMemset((void *)_pDevData, 0, _length * sizeT);
    RTERROR(status, "cudaMemset GpuBuffer::Allocate failed");
  }
}

//---------------------------------------------------------------------------------------------
// Deallocate: free memory associated with a GpuBuffer on the host and on the device
//---------------------------------------------------------------------------------------------
template <typename T> void GpuBuffer<T>::Deallocate()
{
  if (_length == 0) {
    return;
  }
  size_t sizeT = sizeof(T);

  cudaError_t status = cudaSuccess;

  // Free pinned memory allocated by cudaHostAlloc(),
  // or pageable memory allocated by the new operator
  try {
    if (_bPinned) {
      status = cudaFreeHost(_pSysData);
      gpuMemoryInfo::Instance().totalCPUMemory -= _length * sizeT;
      gpuMemoryInfo::Instance().totalGPUMemory -= _length * sizeT;
    } else {
      if (_bSysMem) {

        if (_pSysData) {
         delete[] _pSysData;
        gpuMemoryInfo::Instance().totalCPUMemory -= _length * sizeT;
      }
      }
      status = cudaFree(_pDevData);
      gpuMemoryInfo::Instance().totalGPUMemory -= _length * sizeT;
    }
  } 
  catch (...) {
    RTERROR(status, "cudaFree GpuBuffer::Deallocate failed");
  }
  _pSysData = NULL;
  _pDevData = NULL;
}

//---------------------------------------------------------------------------------------------
// Upload: send information in an array of type T from the host (held in CPU RAM) to the
//         device (GPU).
//
// Arguments:
//   pBuff:  if not NULL, the data will be uploaded from this array rather than _pSysData
//---------------------------------------------------------------------------------------------
template <typename T> void GpuBuffer<T>::Upload(T* pBuff)
{
  if (_length == 0) {
    return;
  }
  size_t sizeT = sizeof(T);

  if (pBuff) {
    cudaError_t status;
    status = cudaMemcpy(_pDevData, pBuff, _length * sizeT, cudaMemcpyHostToDevice);
    RTERROR(status, "cudaMemcpy GpuBuffer::Upload failed");
  }
  else if (_bSysMem) {
    cudaError_t status;
    status = cudaMemcpy(_pDevData, _pSysData, _length * sizeT, cudaMemcpyHostToDevice);
    RTERROR(status, "cudaMemcpy GpuBuffer::Upload failed");
  }
}

//---------------------------------------------------------------------------------------------
// Download: retrieve information in a device array of type T and write it to a host array
//           in CPU RAM.  The data can be copied into a temporary array or the mirror array
//           contained in the same GpuBuffer.
//
// Arguments:
//   pBuff:  if not NULL, the data will be downloaded into this array rather than _pSysData
//---------------------------------------------------------------------------------------------
template <typename T> void GpuBuffer<T>::Download(T* pBuff)
{
  if (_length == 0) {
    return;
  }
  size_t sizeT = sizeof(T);

  if (pBuff) {
    cudaError_t status;
    status = cudaMemcpy(pBuff, _pDevData, _length * sizeT, cudaMemcpyDeviceToHost);
    RTERROR(status, "cudaMemcpy GpuBuffer::Download failed");
  }
  else if (_bSysMem) {
    cudaError_t status;
    status = cudaMemcpy(_pSysData, _pDevData, _length * sizeT, cudaMemcpyDeviceToHost);
    RTERROR(status, "cudaMemcpy GpuBuffer::Download failed");
  }
}

#endif // _GPU_BUFFER
