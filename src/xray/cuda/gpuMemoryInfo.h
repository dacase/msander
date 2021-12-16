#ifndef _GPU_MEMORY_INFO
#define _GPU_MEMORY_INFO
#include "gputypes.h"

//Singleton pattern for the global instance 
class gpuMemoryInfo
{
public:

  // Memory parameters: mainly used by GpuBuffer templates
  aligned_lli totalCPUMemory;             // Approximate total allocated CPU memory
  aligned_lli totalGPUMemory;             // Approximate total allocated CPU memory

  static gpuMemoryInfo& Instance();

protected:

  gpuMemoryInfo(void);
  virtual ~gpuMemoryInfo(void);
};

#endif //_GPU_MEMORY_INFO
