include ../../../config.h

.SUFFIXES: .cpp .cu .o

#CUDA Specific build flags
NVCC=/usr/local/cuda/bin/nvcc -gencode arch=compute_30,code=sm_30 -gencode \
arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode \
arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode \
arch=compute_53,code=sm_53 -gencode arch=compute_60,code=sm_60 -gencode \
arch=compute_61,code=sm_61 -gencode arch=compute_60,code=sm_70 -gencode \
arch=compute_60,code=sm_75 -Wno-deprecated-declarations -use_fast_math -O3 \
-I/usr/local/cuda//include -IB40C -DCUDA --std c++11

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
