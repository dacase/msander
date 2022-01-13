include ../../../config.h

.SUFFIXES: .cpp .cu .o

#CUDA Specific build flags  -- use nvcc in the PATH for now; 
#                              may want to set with configure later
NVCC=nvcc \
   -gencode arch=compute_52,code=sm_52 \
   -gencode arch=compute_53,code=sm_53 \
   -gencode arch=compute_60,code=sm_60 \
   -gencode arch=compute_60,code=sm_70 \
   -gencode arch=compute_60,code=sm_75 \
   -Wno-deprecated-declarations -use_fast_math -O3 \
   -I/usr/local/cuda/include -IB40C $(CUDAFLAGS) --std c++14

CPPOBJS = BulkMask.o BulkMaskCPU.o DPartial.o DPartialCPU.o NonBulk.o \
          NonBulkCPU.o xray_bulk_mask.o xray_dpartial.o xray_non_bulk.o \
          xray_unit_cell.o

CUOBJS = BulkMaskGPU.o	DPartialGPU.o	NonBulkGPU.o

.cpp.o:
	@echo "[XRAY] nvcc $<"
	$(VB)$(NVCC) -c $(CFLAGS) -o $@ $<

.cu.o:
	@echo "[XRAY] nvcc $<"
	$(VB)$(NVCC) -c $(CFLAGS) -o $@ $<

xray_cuda.o: $(CUOBJS) $(CPPOBJS)
	nvcc -dlink -o $@ $(CUOBJS) $(CPPOBJS)

libxray_cuda: $(LIBDIR)/libxray_cuda.a

$(LIBDIR)/libxray_cuda.a: xray_cuda.o
	ar rvs $@ $< $(CUOBJS) $(CPPOBJS)

clean:
	/bin/rm -f *.o
