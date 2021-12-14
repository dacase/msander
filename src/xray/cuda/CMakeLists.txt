
cuda_add_library(
        pmemd_xray_cuda
    STATIC
        src/xray/BulkMask.cpp
        src/xray/BulkMask.h
        src/xray/BulkMaskCPU.cpp
        src/xray/BulkMaskCPU.h
        src/xray/BulkMaskGPU.cu
        src/xray/BulkMaskGPU.h
        src/xray/DPartial.cpp
        src/xray/DPartial.h
        src/xray/DPartialCPU.cpp
        src/xray/DPartialCPU.h
        src/xray/DPartialGPU.cu
        src/xray/DPartialGPU.h
        src/xray/NonBulk.cpp
        src/xray/NonBulk.h
        src/xray/NonBulkCPU.cpp
        src/xray/NonBulkCPU.h
        src/xray/NonBulkGPU.cu
        src/xray/NonBulkGPU.h
        src/xray/xray_bulk_mask.cpp
        src/xray/xray_bulk_mask.h
        src/xray/xray_common.h
        src/xray/xray_dpartial.cpp
        src/xray/xray_dpartial.h
        src/xray/xray_non_bulk.cpp
        src/xray/xray_non_bulk.h
        src/xray/xray_unit_cell.cpp
        src/xray/xray_unit_cell.h
)

target_include_directories(
        pmemd_xray_cuda
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

target_link_libraries(
        pmemd_xray_cuda
        stdc++
        fftw # used in preliminary CPU scaffolding and benchmarking
        pthread # fixme: should be provided by CUDA::cudart
)

# TODO: FindCUDA is deprecated since CMake-3.10:
# Switch to FindCUDAToolkit: https://cmake.org/cmake/help/latest/module/FindCUDAToolkit.html
cuda_add_cufft_to_target(pmemd_xray_cuda)

add_library(PMEMD::xray_cuda ALIAS pmemd_xray_cuda)

target_compile_definitions(
        pmemd_xray_cuda
    PUBLIC
        "-DCUDA"
)