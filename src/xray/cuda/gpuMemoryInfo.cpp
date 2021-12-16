#include "gpuMemoryInfo.h"

gpuMemoryInfo& gpuMemoryInfo::Instance()
{
  static gpuMemoryInfo thisInstance;
  return thisInstance;
};

gpuMemoryInfo::gpuMemoryInfo(void)
{
  totalCPUMemory = 0;
  totalGPUMemory = 0;
};

gpuMemoryInfo::~gpuMemoryInfo(void)
{
  // Blank destructor function
};
