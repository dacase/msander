cmake_minimum_required(VERSION 3.1)
project(fftw)

set(fftw_MAJOR_VERSION 3)
set(fftw_MINOR_VERSION 3)
set(fftw_PATCH_VERSION 0)

# Amber CMake build file for internal FFTW3
# File containing some little functions which wrap things which are very verbose in CMake in something a little easier to type
# Since these functions are used all over the place in the build system, this must be included after Utils.cmake but before anything else.

# allow this build file to work on its own, albeit still inside the Amber source tree
if(INSIDE_AMBER)
	# Define to disable Fortran wrappers.
	set(DISABLE_FORTRAN FALSE)
else()
	option(DISABLE_FORTRAN "Set to true to disable Fortran" FALSE)
	
	include(../../../cmake/AmberBuildSystemInit.cmake)
	
	enable_language(C)
	
	if(NOT DISABLE_FORTRAN)
		enable_language(Fortran)
	endif()
	
	include(AmberBuildSystem2ndInit)
	set(NEEDED_3RDPARTY_TOOLS libm)
	include(3rdPartyTools)
	include(CompilerFlags)
endif()

#------------------------------------------------------------------------
# Figure out what processor optimizations to use
#------------------------------------------------------------------------

#enable optimizations for machines with fused multiply-add
test(HAVE_FMA ${TARGET_ARCH} STREQUAL "ia64" OR ${TARGET_ARCH} STREQUAL "mips" OR ${TARGET_ARCH} STREQUAL "ppc")

# Use our own aligned malloc routine; mainly helpful for Windows systems lacking aligned allocation system-library routines.
test(WITH_OUR_MALLOC "${CMAKE_SYSTEM_NAME}" STREQUAL Windows)

# Use common Windows Fortran mangling styles for the Fortran interfaces.
test(WINDOWS_F77_MANGLING "${CMAKE_SYSTEM_NAME}" STREQUAL Windows)

#------------------------------------------------------------------------
# Hardcoded settings
#------------------------------------------------------------------------
# Define to enable extra FFTW debugging code.
test(FFTW_DEBUG ${CMAKE_BUILD_TYPE} STREQUAL Debug)

# Define to enable debugging malloc.
set(FFTW_DEBUG_MALLOC ${FFTW_DEBUG})

# Define to enable alignment debugging hacks.
set(FFTW_DEBUG_ALIGNMENT ${FFTW_DEBUG})

# Define to enable pseudorandom estimate planning for debugging.
set(FFTW_RANDOM_ESTIMATOR FALSE)

# Define to compile in single precision.
set(FFTW_SINGLE FALSE)
set(BENCHFFT_SINGLE FALSE)

# Define to compile in long-double precision.
set(FFTW_LDOUBLE FALSE)
set(BENCHFFT_LDOUBLE FALSE)

# Define to compile in quad precision.
set(FFTW_QUAD FALSE)
set(BENCHFFT_QUAD FALSE)

# Define to enable SSE/SSE2 optimizations.
set(HAVE_SSE2 ${SSE})

# Define to enable AVX optimizations.
set(HAVE_AVX FALSE)

# Define to enable Altivec optimizations.
set(HAVE_ALTIVEC FALSE)

# Define to enable MIPS paired-single optimizations.
set(HAVE_MIPS_PS FALSE)

# Use low-precision timers, making planner very slow
set(WITH_SLOW_TIMER FALSE)

# Define to enable use of MIPS ZBus cycle-counter.
set(HAVE_MIPS_ZBUS_TIMER FALSE)

# Assume that stack is aligned to (1<<X) bytes
#do not set(INCOMING_STACK_BOUNDARY )

# use OpenMP directives for parallelism
set(HAVE_OPENMP ${OPENMP}) #from Amber configuration code

set(ENABLE_MPI ${MPI})

# compile FFTW SMP threads library
test(ENABLE_THREADS NOT ${CMAKE_THREAD_LIBS_INIT} MATCHES ".*-NOTFOUND") #Amber code finds the Threads package

test(SMP OPENMP OR ENABLE_THREADS)

# combine threads into main libfftw3
# currently COMBINED_THREADS=FALSE can't work with shared libraries on Windows.  
# The functions used by the threading library don't have the required DLL exports and imports.
set(COMBINED_THREADS TRUE) 

#prefer OpenMP combined threads
test(COMBINE_OPENMP_THREADS HAVE_OPENMP AND COMBINED_THREADS)
test(COMBINE_REGULAR_THREADS COMBINED_THREADS AND ENABLE_THREADS AND (NOT HAVE_OPENMP))


#------------------------------------------------------------------------
# Compiler flags
#------------------------------------------------------------------------

if(${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
	# Stop icc from defining __GNUC__, except on MacOS where this fails
	if(NOT APPLE)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -no-gcc")
	endif()
	
elseif(${CMAKE_C_COMPILER_ID} STREQUAL "HP")
	
	check_c_compiler_flag("-Wp,-H128000" SUPPORTS_HP_CPP_FLAG)
	# must (sometimes) manually increase cpp limits to handle fftw3.h
	
	if(SUPPORTS_HP_CPP_FLAG)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wp,-H128000")
	endif()
	
elseif(${CMAKE_C_COMPILER_ID} STREQUAL "PGI")
	# -Masmkeyword required for asm("") cycle counters
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Masmkeyword")
	
	# 12.6 < PGI version < 13.6 will not compile FFTW3 with SSE
	if(USE_SYSTEM_FFTW AND HAVE_SSE2)
		if(${CMAKE_C_COMPILER_VERSION} VERSION_GREATER 12.6 AND ${CMAKE_C_COMPILER_VERSION} VERSION_LESS 13.6)
			message(FATAL_ERROR "PGI versions in-between 12.6 and 13.6 cannot compile FFTW3 with SSE. FFTW3 will be built without SSE.")
		endif()
		
		set(HAVE_SSE2 FALSE)
		
	endif()
	
endif()

#SSE flags
if(${CMAKE_C_COMPILER_ID} STREQUAL "Intel" OR ${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
	# SSE/SSE2
	if(HAVE_SSE2)
	    if(FFTW_SINGLE)
	    	set(SSE_FLAG -msse)
    	else()
    		set(SSE_FLAG -msse2)
    	endif()
    	
    	check_c_compiler_flag(${SSE_FLAG} SUPPORTS_SSE_FLAG)
    	
    	if(NOT SUPPORTS_SSE_FLAG)
    		message(FATAL_ERROR "You requested SSE support, but the sse compiler flag does not work")
		endif()
		
		set(SSE2_CFLAGS ${SSE_FLAG})
		
	endif()

	# AVX
	if(HAVE_AVX)
    	check_c_compiler_flag(-mavx SUPPORTS_AVX_FLAG)
    	
    	if(NOT SUPPORTS_AVX_FLAG)
    		message(FATAL_ERROR "You requested AVX support, but the -mavx compiler flag does not work")
		endif()
		
		set(AVX_CFLAGS -mavx)
	endif()

	#Altivec
	if(HAVE_ALTIVEC)
	    check_c_compiler_flag(-faltivec SUPPORTS_FALTIVEC_FLAG)
	    
	    check_c_compiler_flag("-maltivec -mabi=altivec" SUPPORTS_MALTIVEC_FLAG)
	    
	    if(SUPPORTS_FALTIVEC_FLAG)
	    	set(ALTIVEC_CFLAGS -faltivec)
    	elseif(SUPPORTS_MALTIVEC_FLAG)
	    	set(ALTIVEC_CFLAGS -maltivec -mabi=altivec)
	    else()
    		message(FATAL_ERROR "You requested Altivec support, but no altivec compiler flag was found.")
		endif()
	else()
		set(ALTIVEC_CFLAGS "")
	endif()
elseif(${CMAKE_C_COMPILER} STREQUAL PGI)
	# SSE/SSE2
	if(HAVE_SSE2)
	    set(SSE_FLAG -fastsse)

    	check_c_compiler_flag(${SSE_FLAG} SUPPORTS_SSE_FLAG)
    	
    	if(NOT SUPPORTS_SSE_FLAG)
    		message(FATAL_ERROR "You requested SSE support, but the sse compiler flag does not work")
		endif()
		
		set(SSE2_CFLAGS ${SSE_FLAG})
		
	endif()
endif()

# add stack alignment CFLAGS if so requested
if(DEFINED INCOMING_STACK_BOUNDARY)
	set(ISB_FLAG -mincoming-stack-boundary=${INCOMING_STACK_BOUNDARY})

	check_c_compiler_flag(${ISB_FLAG} SUPPORTS_ISB_FLAG)
    	
	if(NOT SUPPORTS_ISB_FLAG)
		message(SEND_ERROR "You requested to set the incoming stack boundary, but the compiler flag to do that doesn't work")
	endif()

	set(STACK_ALIGN_CFLAGS ${ISB_FLAG})
else()
	set(STACK_ALIGN_CFLAGS "")
endif()

#------------------------------------------------------------------------
# system introspection checks
#------------------------------------------------------------------------
# header checks
check_all_includes(C alloca.h
sys/types.h
sys/time.h
dlfcn.h
libintl.h
malloc.h
stddef.h
string.h
strings.h
unistd.h
limits.h
c_asm.h # Header file for enabling asm() on Digital Unix  
stdlib.h
stdarg.h
float.h
intrinsics.h # cray unicos
stdint.h
mach/mach_time.h
sys/sysctl.h)

check_include_file(altivec.h HAVE_ALTIVEC_H ${ALTIVEC_CFLAGS}) # MacOS X altivec detection

if(HAVE_STRING_H AND HAVE_STDLIB_H AND HAVE_STDARG_H AND HAVE_FLOAT_H)
	set(STDC_HEADERS TRUE)
endif()

# function checks

# the Automake version of this check goes all Master Prepper and builds its own version of alloca if it is not found
# for now, we'll take the simpler approach and just disable it if we don't find it.

# it seems to be defined as a macro on some platforms, so we use check_symbol_exists() since we know the header
if(HAVE_ALLOCA_H)
	check_symbol_exists(alloca alloca.h HAVE_ALLOCA)
else()
	check_function_exists(alloca HAVE_ALLOCA)
endif()

test(FFTW_ENABLE_ALLOCA HAVE_ALLOCA AND NOT (${TARGET_ARCH} STREQUAL "x86_64" AND MINGW))

check_function_exists(vprintf HAVE_VPRINTF)

if(FFTW_QUAD)
	if(${CMAKE_C_COMPILER_ID} STREQUAL GNU AND ${CMAKE_C_COMPILER_VERSION} VERSION_LESS 4.6.0)
		message(FATAL_ERROR "gcc 4.6 or later is required for quad precision support!")
	endif()
	
	set(CMAKE_REQUIRED_LIBRARIES quadmath)
	
	check_function_exists(sinq HAVE_DECL_SINQ)
	check_function_exists(cosq HAVE_DECL_COSQ)
	
	if(NOT HAVE_DECL_SINQ)
		message(FATAL_ERROR "quad precision requires libquadmath for quad-precision trigonometric routines!")
	endif()
	
endif()

# the main AMBER buildfiles already find libm,so we pass it down to the compile tests from there
if(libm_ENABLED)
	set(CMAKE_REQUIRED_LIBRARIES ${LIBM})
endif()
 
# check for the following functions in the math library and the standard library                  
check_all_functions(BSDgettimeofday
gettimeofday
gethrtime
read_real_time
time_base_to_time
sqrt
memset
_mm_malloc
_mm_free
clock_gettime
mach_absolute_time
sysctl
abort
snprintf
drand48
memalign 
posix_memalign
sinl
cosl
tanl) 

# check function definitions
check_symbol_exists(drand48 stdlib.h HAVE_DECL_DRAND48)
check_symbol_exists(srand48 stdlib.h HAVE_DECL_SRAND48)
check_symbol_exists(memalign malloc.h HAVE_DECL_MEMALIGN)
check_symbol_exists(posix_memalign stdlib.h HAVE_DECL_POSIX_MEMALIGN)
check_symbol_exists(sinl math.h HAVE_DECL_SINL)
check_symbol_exists(cosl math.h HAVE_DECL_COSL)

if(HAVE_INTRINSICS_H)
	check_symbol_exists(_rtc intrinsics.h HAVE__RTC)
endif()

if(FFTW_LDOUBLE AND NOT (HAVE_SINL AND HAVE_COSL))
	message(FATAL_ERROR "long-double precision requires long-double trigonometric routines!")
endif()

check_symbol_exists(isnan math.h HAVE_ISNAN)

# Check type sizes
check_all_types(double float int long "long long" "long double" "unsigned int" "unsigned long" "unsigned long long" ptrdiff_t size_t uintptr_t)

set(CMAKE_EXTRA_INCLUDE_FILES sys/time.h)
check_type_size(hrtime_t HRTIME_T)

if(MPI)
	set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
	set(CMAKE_EXTRA_INCLUDE_FILES mpi.h)
	
	check_type_size(MPI_Fint SIZEOF_MPI_FINT)
endif()

# why are we checking the size of something inside fftw, you ask?
# It's an enum, and the developers weren't sure that its size would be consistent.
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/api)
set(CMAKE_EXTRA_INCLUDE_FILES fftw3.h)
check_type_size(fftw_r2r_kind SIZEOF_FFTW_R2R_KIND)


#figure out GCC stack alignment
set(CC_HAS_STACK_ALIGNMENT FALSE)
check_c_compiler_flag(-mpreferred-stack-boundary=4 SUPPORTS_PREFERRED_STACK_BOUNDARY)
if(SUPPORTS_PREFERRED_STACK_BOUNDARY)
	check_c_compiler_flag(-malign-double SUPPORTS_MALIGN_DOUBLE)
	if(SUPPORTS_MALIGN_DOUBLE)
		set(TEST_CFLAGS "-O -malign-double")
	else()
		set(TEST_CFLAGS -O)
	endif()
	try_run(ALIGNMENT_RESULT_CODE ALIGNMENT_COMPILE_SUCCEEDED ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_gcc_aligns_stack.c COMPILE_DEFINITIONS "-DCMAKE_C_FLAGS=${TEST_CFLAGS}") 
		
	if(ALIGNMENT_COMPILE_SUCCEEDED AND "${ALIGNMENT_RESULT_CODE}" EQUAL 0)
		set(CC_HAS_STACK_ALIGNMENT TRUE)
	else()
		set(CC_HAS_STACK_ALIGNMENT FALSE)
	endif()
	
endif()
		
message(STATUS "FFTW: compiler supports stack alignment: ${CC_HAS_STACK_ALIGNMENT}")

#change name of pthread variable.... twice
set(HAVE_PTHREAD ${CMAKE_USE_PTHREADS_INIT})
set(USING_POSIX_THREADS ${CMAKE_USE_PTHREADS_INIT})

#------------------------------------------------------------------------
# Win32 DLL
#------------------------------------------------------------------------

if(TARGET_WINDOWS AND MPI)
	# fftw_mpi uses private fftw API which does not have dll exports and imports.  If we use MPI on Windows, we cannot build fftw as shared.
	set(BUILD_SHARED_LIBS FALSE)
endif()

test(WINDOWS_DLL TARGET_WINDOWS AND BUILD_SHARED_LIBS)

# set up the definitions to make FFTW use imports and exports
if(WINDOWS_DLL)
	# NOTE: each of these are used in different places.  We need to supply both (I think).
	set(DLL_EXPORT_DEFINITIONS -DDLL_EXPORT -DFFTW_DLL -DCOMPILING_FFTW)
	
	set(DLL_IMPORT_DEFINITIONS -DDLL_IMPORT -DFFTW_DLL -DCALLING_FFTW)
else()

	set(DLL_EXPORT_DEFINITIONS "")
	set(DLL_IMPORT_DEFINITIONS "")	
endif()

#------------------------------------------------------------------------
# Set up cycle counter
#------------------------------------------------------------------------
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/kernel)
check_c_source_compiles(
"#include \"cycle.h\"
#ifndef HAVE_TICK_COUNTER
#  error No cycle counter
#endif

void main(){}" HAVE_TICK_COUNTER)

if((NOT WITH_SLOW_TIMER) AND (NOT HAVE_TICK_COUNTER))
	unset(HAVE_TICK_COUNTER CACHE)
	message(WARNING "No cycle counter found.  FFTW will use ESTIMATE mode for all plans.  See the manual for more information.")
endif()

#------------------------------------------------------------------------
# Set up Fortran
#------------------------------------------------------------------------
if(NOT DISABLE_FORTRAN)
	include(FortranCInterface)
		
	FortranCInterface_VERIFY() #check that Fortran can link to C
	
	FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/fortran-mangling.h MACRO_NAMESPACE FORTRAN_MANGLE)
	
	if(NOT FortranCInterface_GLOBAL_FOUND)
		message(FATAL_ERROR "Failed to detect the Fortran mangling!  Internal fftw3 cannot build!")
	endif()
	
	# The fftw headers include the API headers twice, once with the underscore function and one without
	# This causes multiple definition errors if the functions are the same
	# So we check for that here.
	test(F77_FUNC_EQUIV ${FortranCInterface_GLOBAL_MACRO} STREQUAL ${FortranCInterface_GLOBAL__MACRO})
	
	
	test(WITH_G77_WRAPPERS ${CMAKE_C_COMPILER_ID} STREQUAL GNU)
endif()

#------------------------------------------------------------------------
# configuration file
#------------------------------------------------------------------------

#embed compiler command line into library
set(FFTW_CC "${CMAKE_C_COMPILER} ${CMAKE_C_FLAGS}")

set(PACKAGE_VERSION "${fftw_MAJOR_VERSION}.${fftw_MINOR_VERSION}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

#------------------------------------------------------------------------
# function for building generated SIMD codelets
#------------------------------------------------------------------------

#list of built codelet libraries to link into the final library
#this must be a cache variable because it needs to be set in a lower subdirectory
set(LIBFFTW_SIMD_OBJS "" CACHE INTERNAL "Stores list of SIMD object libraries to link into the final FFTW library")

function(generate_and_build_fft_subdir LIBRARY_NAME COMMON_DIR OPTIMIZATION_NAME OPTIMIZATION_HEADER) #CODELETS...	

	set(GENERATED_CODELET_DIR ${CMAKE_CURRENT_BINARY_DIR}/${OPTIMIZATION_NAME})
	file(MAKE_DIRECTORY ${GENERATED_CODELET_DIR})
	
	#lit of all the files we need to compile
	set(CODELET_FILES "")
	
	foreach(CODELET ${ARGN})

		set(CODELET_FILE ${GENERATED_CODELET_DIR}/${CODELET})

		# we use file(GENERATE) because it doesn't update the file if the content doesn't change, so we don't trigger unnecessary rebuilds		
		file(GENERATE OUTPUT ${CODELET_FILE} 
			CONTENT
"/* Generated automatically.  DO NOT EDIT! */
#define SIMD_HEADER \"${OPTIMIZATION_HEADER}\"
#include \"${COMMON_DIR}/${CODELET}\"")
		
		list(APPEND CODELET_FILES ${CODELET_FILE})
	endforeach()
	
	add_library(${LIBRARY_NAME} OBJECT ${CODELET_FILES})
	target_include_directories(${LIBRARY_NAME} PRIVATE ${COMMON_DIR}/..)
	make_pic_if_needed(${LIBRARY_NAME})
	
	set(LIBFFTW_SIMD_OBJS ${LIBFFTW_SIMD_OBJS} ${LIBRARY_NAME} CACHE INTERNAL "" FORCE)	
	
endfunction(generate_and_build_fft_subdir)

#------------------------------------------------------------------------
# recurse into subdirectories
#------------------------------------------------------------------------
include_directories(simd-support kernel api dft rdft reodft libbench2 tests)

add_subdirectory(api)
add_subdirectory(dft)
add_subdirectory(rdft)
add_subdirectory(kernel)
add_subdirectory(libbench2)
add_subdirectory(reodft)
add_subdirectory(simd-support)

if(SMP)
	add_subdirectory(threads)
endif()

#------------------------------------------------------------------------
# create the actual library
#------------------------------------------------------------------------
if(COMBINE_OPENMP_THREADS)
	set(THREAD_OBJECTS fftw3_omp_obj)
elseif(COMBINE_REGULAR_THREADS)
	set(THREAD_OBJECTS fftw3_threads_obj)
elseif()
	set(THREAD_OBJECTS "")
endif()

set(ALL_FFTW_OBJ_LIBS dft dft_scalar dft_scalar_codelets fftw_kernel fftw_api reodft rdft rdft_scalar rdft_scalar_r2cb rdft_scalar_r2cf rdft_scalar_r2r simd_support ${THREAD_OBJECTS} ${LIBFFTW_SIMD_OBJS})

# put the generator expression wrapper around each library name
set(ALL_FFTW_OBJ_LIBS_GENEX "")
foreach(OBJ_LIB ${ALL_FFTW_OBJ_LIBS})
	# add DLL export definitions
	target_compile_definitions(${OBJ_LIB} PRIVATE ${DLL_EXPORT_DEFINITIONS})

	list(APPEND ALL_FFTW_OBJ_LIBS_GENEX $<TARGET_OBJECTS:${OBJ_LIB}>)
endforeach()



add_library(fftw ${ALL_FFTW_OBJ_LIBS_GENEX})
	
# set interface include directories properly for libraries linking to fftw
target_include_directories(fftw INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/api ${CMAKE_CURRENT_BINARY_DIR}/api)	

if(HAVE_OPENMP AND COMBINE_OPENMP_THREADS)
	# link the OpenMP libraries
	set_property(TARGET fftw PROPERTY LINK_FLAGS ${OpenMP_C_FLAGS})
endif()

if(ENABLE_THREADS AND COMBINE_REGULAR_THREADS)
	target_link_libraries(fftw Threads::Threads)
endif()

if(libm_ENABLED)
	target_link_libraries(fftw ${LIBM})
endif()

#library is called "fftw" in CMake, but needs to be named "fftw3"
set_property(TARGET fftw PROPERTY OUTPUT_NAME "fftw3")

install_libraries(fftw)

#------------------------------------------------------------------------
# build remaining programs and libraries
#------------------------------------------------------------------------

# CALLING_FFTW changes the dll exports to imports
if(WINDOWS_DLL)
	target_compile_definitions(fftw INTERFACE CALLING_FFTW)
endif()


# Figure out which thread lib for supporting programs to link
# prefer openmp
set(SUPPORT_PROGRAM_THREAD_LIB "")
if(NOT COMBINED_THREADS)
	if(HAVE_OPENMP)
		set(SUPPORT_PROGRAM_THREAD_LIB fftw3_omp)
	elseif(ENABLE_THREADS)
		set(SUPPORT_PROGRAM_THREAD_LIB fftw3_threads)
	endif()
endif()

add_subdirectory(tests)
add_subdirectory(tools)

if(ENABLE_MPI)
	add_subdirectory(mpi)
endif()