#ifndef _GPU_CONTEXT
#define _GPU_CONTEXT

#include "base_gpuContext.h"
  typedef base_gpuContext* gpuContext;
  typedef base_gpuContext _gpuContext;

// Singleton pattern for the global instance 
class theGPUContext {

public:
  static gpuContext thisInstance;
  static bool init;
  static bool shutdown;

  static gpuContext GetPointer(void) {
    return (thisInstance);
  };
  
  static gpuContext Initialize(void)
  {
    if (init)
    {
      printf("CUDA GPU context already initialized, you're up to no good, shutting down.\n");
      Shutdown();
      exit(-1);
    }
    else
    {
      thisInstance = new(_gpuContext);
      init = true;
    }
    return thisInstance;
  }
  
  static void Shutdown(void)
  {
    if (shutdown)
    {
      printf("CUDA GPU context does not exist anymore, no need to try again.\n");
    }    
    else if (init)
    {
      delete thisInstance;
    }
    else
    {
      printf("CUDA GPU context does not exist yet\n");
    }
  }
  

protected:

  theGPUContext(void);
  virtual ~theGPUContext(void);
};

#endif //_GPU_CONTEXT
