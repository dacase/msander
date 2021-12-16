include ../../../config.h

.SUFFIXES: .cpp .cu .o

CPPOBJS = BulkMask.o BulkMaskCPU.o DPartial.o DPartialCPU.o NonBulk.o \
          NonBulkCPU.o xray_bulk_mask.o xray_dpartial.o xray_non_bulk.o \
          xray_unit_cell.o

CUOBJS = BulkMaskGPU.o	DPartialGPU.o	NonBulkGPU.o

.cpp.o:
	@echo "[XRAY] nvcc $<"
	$(VB)nvcc -c $(CFLAGS) -o $@ $<

.cu.o:
	@echo "[XRAY] nvcc $<"
	$(VB)nvcc -c $(CFLAGS) -o $@ $<

xray_cuda.o: $(CUOBJS) $(CPPOBJS)
	nvcc -dlink -o $@ $(CUOBJS) $(CPPOBJS)

libxray_cuda: $(LIBDIR)/libxray_cuda.a

$(LIBDIR)/libxray_cuda.a: xray_cuda.o
	ar rvs $@ $< $(CUOBJS) $(CPPOBJS)

clean:
	/bin/rm *.o *.a
