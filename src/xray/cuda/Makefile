include ../../../config.h

.SUFFIXES: .cpp .cu .o

#CUDA Specific build flags  -- use nvcc in the PATH for now; 
#                              may want to set with configure later
NVCC=nvcc \
   -Wno-deprecated-declarations -use_fast_math -O3 \
   -I/usr/local/cuda/include -IB40C $(CUDAFLAGS) --std c++14

CPPOBJS = BulkMask.o BulkMaskCPU.o DPartial.o DPartialCPU.o NonBulk.o \
          NonBulkCPU.o xray_bulk_mask.o xray_dpartial.o xray_non_bulk.o \
          xray_unit_cell.o

CUOBJS = BulkMaskGPU.o	DPartialGPU.o	NonBulkGPU.o

.cpp.o:
	@echo "[XRAY] nvcc $<"
	$(VB)$(NVCC) -c -I$(INCDIR) -o $@ $<

.cu.o:
	@echo "[XRAY] nvcc $<"
	$(VB)$(NVCC) -c -o $@ $<

xray_cuda.o: $(CUOBJS) $(CPPOBJS)
	@echo "[XRAY] nvcc -dlink $@"
	$(VB)$(NVCC) -dlink -o $@ $(CUOBJS) $(CPPOBJS)

libxray_cuda: $(LIBDIR)/libxray_cuda.a

$(LIBDIR)/libxray_cuda.a: xray_cuda.o
	ar rvs $@ $< $(CUOBJS) $(CPPOBJS)

clean:
	/bin/rm -f *.o
