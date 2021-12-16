//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

#ifndef __GPUTYPES_H__
#define __GPUTYPES_H__
#include <stdio.h>
#include <cassert>
#if !defined(_WIN32)
#include <unistd.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#endif
#include <sys/types.h>
#include <cuda.h>
#include <cufft.h>
#include <curand.h>
#include <vector_functions.h>
#include <cuda_runtime_api.h>
#include <builtin_types.h>
#include <cstring>


//---------------------------------------------------------------------------------------------
// Enforce use of CUDA 5.0 due to GK110 issues with 4.2
//---------------------------------------------------------------------------------------------
#if defined(CUDA_VERSION) && (CUDA_VERSION < 5000)
#error "CUDA support requires the use of a 5.0 or later CUDA toolkit.  "
#error "Aborting compilation."
#endif

//---------------------------------------------------------------------------------------------
// Control of single and double precision use. If use_DPFP is not defined then the code runs
// in the default SPFP mode which makes mixed use of single and double precision as needed.
// Defining use_DPFP makes it use double precision throughout while defining use_SPFP makes it
// use single precision for nonbonded force components, double-precision for bonded force
// components and SHAKE, and 24.40 bit fixed-point for force accumulation.
//
// Note using Double Precision throughout will cause, in most cases, large performance
// degradation.
//
// This option is now set in the configure script and the original defines are left here just
// for reference.
//---------------------------------------------------------------------------------------------
//#define use_DPFP
//#define use_SPFP

//TL: enforce correct usage of const in different env.(CPU or GPU)
#ifdef GL_CONST
#undef GL_CONST
#endif
#ifdef  __CUDA_ARCH__
#define GL_CONST  static const __constant__
#else
#define GL_CONST  static const
#endif

//---------------------------------------------------------------------------------------------
// Enforce definition of one and only one precision mode
//---------------------------------------------------------------------------------------------
/*  dac kludge, 12/2021  */
#define use_SPFP 1

#if !(defined(use_DPFP) && !defined(use_SPFP)) && \
    !(defined(use_SPFP) && !defined(use_DPFP))
#error "You must define one and only one precision mode (use_SPFP, "
#error "or use_DPFP) to build pmemd.cuda. Aborting compilation."
#endif

#if defined(_MSC_VER)
#define __align(_boundary_size) __declspec(align(_boundary_size))
#else
#define __align(_boundary_size) __attribute__((aligned(_boundary_size)))
#endif

//---------------------------------------------------------------------------------------------
// Type definitions to declare aliases for various P.O.D. types
//---------------------------------------------------------------------------------------------
typedef double                 __align(8) aligned_double;
typedef unsigned long int      __align(8) aligned_uli;
typedef long long int          __align(8) aligned_lli;
typedef unsigned long long int __align(8) PMEUllInt;
typedef unsigned int                      uint;

typedef long long int PMEAccumulator;
typedef long long int PMEForceAccumulator;
typedef long long int PMEEnergyAccumulator;
typedef long long int PMEVirialAccumulator;

//---------------------------------------------------------------------------------------------
// Additional, data type aliasing dependent on the precision mode
//---------------------------------------------------------------------------------------------
typedef double             PMEDouble;
#if defined(use_DPFP)
typedef double             PMEFloat;
typedef double             PMEForce;
typedef double             PMEEnergy;
typedef double             PMEVirial;
typedef double2            PMEDouble2;
typedef double3            PMEDouble3;
typedef double4            PMEDouble4;
typedef double2            PMEFloat2;
typedef double3            PMEFloat3;
typedef double4            PMEFloat4;
typedef cufftDoubleComplex PMEComplex;
#  ifdef MPI
static const MPI_Datatype MPI_PMEDOUBLE = MPI_DOUBLE_PRECISION;
static const MPI_Datatype MPI_PMEFLOAT = MPI_DOUBLE_PRECISION;
static const MPI_Datatype MPI_PMEACCUMULATOR = MPI_LONG_LONG_INT;
#  endif
#elif defined(use_SPFP)
typedef long long int      PMEForce;
typedef long long int      PMEEnergy;
typedef long long int      PMEVirial;
typedef float              PMEFloat;
typedef double2            PMEDouble2;
typedef double3            PMEDouble3;
typedef double4            PMEDouble4;
typedef float2             PMEFloat2;
typedef float3             PMEFloat3;
typedef float4             PMEFloat4;
typedef cufftComplex PMEComplex;
#  ifdef MPI
static const MPI_Datatype MPI_PMEDOUBLE = MPI_DOUBLE_PRECISION;
static const MPI_Datatype MPI_PMEFLOAT = MPI_FLOAT;
static const MPI_Datatype MPI_PMEACCUMULATOR = MPI_LONG_LONG_INT;
#  endif
#endif
// End branch for aliasing based on the precision mode

enum SM_VERSION {
  SM_10,
  SM_11,
  SM_12,
  SM_13,
  SM_2X,
  SM_3X,
};

// Miscellaneous constants
enum {
  GRID                  = 32,
  GRID_BITS             = 5,
  GRID_BITS_MASK        = (GRID - 1),
  GRID_PADDING_MASK     = 0xffffffff - (GRID - 1),
  CMAP_RESOLUTION       = 24,
  CMAP_DIMENSION        = CMAP_RESOLUTION + 4,
  CMAP_STEP_SIZE        = 15,
  NEIGHBOR_CELLS        = 14,
  CELL_IDX_BITS         = 10,
  CELL_IDY_BITS         = 10,
  CELL_IDZ_BITS         = 10,
  CELL_IDY_SHIFT        = CELL_IDX_BITS,
  CELL_IDZ_SHIFT        = CELL_IDX_BITS + CELL_IDY_BITS,
  PENCIL_IDY_SHIFT      = 14,
  PENCIL_IDZ_SHIFT      = 22,
  CELL_ID_MASK          = 0x0000003f,
  CELL_HASHX_BITS       = 2,
  CELL_HASHY_BITS       = 2,
  CELL_HASHZ_BITS       = 2,
  CELL_HASH_BITS        = (CELL_HASHX_BITS + CELL_HASHY_BITS + CELL_HASHZ_BITS),
  CELL_HASHX            = (1 << CELL_HASHX_BITS),
  CELL_HASHY            = (1 << CELL_HASHY_BITS),
  CELL_HASHZ            = (1 << CELL_HASHZ_BITS),
  CELL_HASHXY           = CELL_HASHX * CELL_HASHY,
  CELL_HASH_CELLS       = (1 << CELL_HASH_BITS),
  CELL_XC_SHIFT         = 0,
  CELL_YC_SHIFT         = 10,
  CELL_ZC_SHIFT         = 20,

  NLRECORD_CELL_COUNT_BITS = 4,
  NLRECORD_CELL_COUNT_MASK = ((1 << NLRECORD_CELL_COUNT_BITS) - 1),
  NLRECORD_YOFFSET_SHIFT   = NLRECORD_CELL_COUNT_BITS,
  NLRECORD_CELL_SHIFT      = 4,
  NLRECORD_CELL_TYPE_MASK  = ((1 << NLRECORD_CELL_SHIFT) - 1),

  NLENTRY_YMAX_SHIFT       = 3,
  NLENTRY_HOME_CELL_MASK   = 1,

  NLATOM_CELL_SHIFT        = 4,
  NLATOM_CELL_TYPE_MASK    = ((1 << NLATOM_CELL_SHIFT) - 1),

  NLEXCLUSION_SHIFT        = 8,
  NLEXCLUSION_ATOM_MASK    = ((1 << NLEXCLUSION_SHIFT) - 1),
  AMD_E_DIHEDRAL_OFFSET    = 19,
  GAMD_E_DIHEDRAL_OFFSET   = 20,
  NFX_OFFSET               = 21,
  NFY_OFFSET               = 22,
  NFZ_OFFSET               = 23,
  NFX1_OFFSET              = 21,
  NFY1_OFFSET              = 22,
  NFZ1_OFFSET              = 23,
  NFX2_OFFSET              = 24,
  NFY2_OFFSET              = 25,
  NFZ2_OFFSET              = 26,
  VIRIAL_OFFSET            = 27,
  ENERGY_TERMS             = 33,
  EVDWE_OFFSET             = 33,
  EEDE_OFFSET              = 34,
  EELTE_OFFSET             = 35,
  VIRIALE_OFFSET           = 36,
  EXTENDED_ENERGY_TERMS    = 42,
  AFE_TERMS                = 25,
  PADDING                  = 16,
  PADDING_MASK             = 0xfffffff0,
  TEST_DENSITY             = 37,
#ifdef AWSMM
  PP_ENERGIES           = 1,
  PP_VELOCITIES         = 2,
  PP_RMSD               = 4,
#endif
};

//---------------------------------------------------------------------------------------------
// NLRecordStruct: the "Neighbor List Record" struct keeps track of the non-bonded pair list.
//---------------------------------------------------------------------------------------------
struct NLRecordStruct {
  unsigned int neighborCell[NEIGHBOR_CELLS];  // Bits 0:3 give the cell type,
                                              //   bits 4:31 give the cell ID
  unsigned int homeCell;                      // Homecell ID
  unsigned int neighborCells;                 // Bits 0:3 give the number of
                                              //   neighbor cells, bits 4:31
                                              //   give the y offset
};

//---------------------------------------------------------------------------------------------
// NLRecord: this union makes it possible to read each NLRecordStruct as an array of unsigned
//           integers, called "array" no less.
//---------------------------------------------------------------------------------------------
union NLRecord {
  NLRecordStruct NL;
  unsigned int array[NEIGHBOR_CELLS + 2];
};

//---------------------------------------------------------------------------------------------
// NLEntryStruct: the neighbor list entry
//---------------------------------------------------------------------------------------------
struct NLEntryStruct {
  unsigned int ypos;       // First y atom
  unsigned int ymax;       // Bit 0 is the home cell flag.
                           // Bits 3:31 are the number of y atoms (1-32)
  unsigned int xatoms;     // Number of x atoms
  unsigned int offset;     // Offset to entry's atom list
};

//---------------------------------------------------------------------------------------------
// NLEntry: this union makes it possible to read the NLEntry struct as an array of unsigned
//          integers.  Note that NLRecord has an eponymous attribute NL but the two should not
//          be confused when dealing with each of these data types.
//---------------------------------------------------------------------------------------------
union NLEntry {
  NLEntryStruct NL;
  unsigned int array[4];
};

//---------------------------------------------------------------------------------------------
// FloatShift: a union to enable bit shifts on floating point numbers.  With the right steup
//             (e.g. add 2.0) a floating point number greater than zero can be converted to an
//             unsigned integer indexing into a lookup table by taking the last few bits of
//             the exponent (it will be > 10000000 for 2.0, and for floating point numbers of
//             limited size (i.e. less than the squared direct space cutoff only the last few
//             bits will ever be nonzero) and the first five bits of the mantissa (note the
//             frequent use of 32 intervals in the erfc coefficients table).
//---------------------------------------------------------------------------------------------
union FloatShift {
  PMEFloat f;
#ifndef use_DPFP
  int i;
  unsigned int ui;
#else
  long long int i;
  unsigned long long int ui;
#endif
};

//---------------------------------------------------------------------------------------------
// Kernel dimensions - we only support SM 3.0 or better due to double-precision and atomic
// operation requirements.
//---------------------------------------------------------------------------------------------
static const int THREADS_PER_BLOCK                        = 1024;
static const int NLCALCULATE_OFFSETS_THREADS_PER_BLOCK    = 720;
static const int NLBUILD_NEIGHBORLIST8_THREADS_PER_BLOCK  = 128;
static const int NLBUILD_NEIGHBORLIST16_THREADS_PER_BLOCK = 128;
static const int NLBUILD_NEIGHBORLIST32_THREADS_PER_BLOCK = 128;
#if (__CUDA_ARCH__ == 750)
static const int NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER   = 8;
#else
static const int NLBUILD_NEIGHBORLIST_BLOCKS_MULTIPLIER   = 10;
#endif
static const int LOCALFORCES_BLOCKS                       = 16;
static const int LOCALFORCES_THREADS_PER_BLOCK            = 64;
static const int AFE_EXCHANGE_THREADS_PER_BLOCK           = 256;
static const int CHARMMFORCES_BLOCKS                      = 16;
static const int CHARMMFORCES_THREADS_PER_BLOCK           = 64;
static const int NMRFORCES_BLOCKS                         = 16;
static const int NMRFORCES_THREADS_PER_BLOCK              = 64;
static const int CLEARFORCES_THREADS_PER_BLOCK            = 1024;
static const int NLCLEARFORCES_THREADS_PER_BLOCK          = 1024;
static const int NLREDUCEFORCES_THREADS_PER_BLOCK         = 1024;
static const int REDUCEFORCES_THREADS_PER_BLOCK           = 1024;
static const int REDUCEBUFFER_THREADS_PER_BLOCK           = 1024;
static const int SHAKE_THREADS_PER_BLOCK                  = 64;
static const int SHAKE_BLOCKS                             = 16;
static const int UPDATE_THREADS_PER_BLOCK                 = 256;
static const int GENERAL_THREADS_PER_BLOCK                = 1024;
static const int P2M_THREADS_PER_BLOCK                    = 96;
#if (__CUDA_ARCH__ == 750)
static const int P2M_BLOCKS_MULTIPLIER                    = 10;
#else
static const int P2M_BLOCKS_MULTIPLIER                    = 16;
#endif
static const int FILTER_IMAGE_THREADS_PER_BLOCK           = 128;
static const int FILTER_IMAGE_WARPS_PER_BLOCK             = 4;
static const int VS_THREADS_PER_BLOCK                     = 64;
static const int VS_BLOCKS_MULTIPLIER                     = 16;
static const int VS_BLOCKS_MINIMUM                        = 12;
#ifdef use_DPFP
static const int GBBORNRADII_THREADS_PER_BLOCK            = 1024;
static const int GBBORNRADII_BLOCKS_MULTIPLIER            = 1;
static const int GBNONBONDENERGY1_THREADS_PER_BLOCK       = 1024;
static const int GBNONBONDENERGY1_BLOCKS_MULTIPLIER       = 1;
static const int GBNONBONDENERGY2_THREADS_PER_BLOCK       = 512;
static const int GBNONBONDENERGY2_BLOCKS_MULTIPLIER       = 1;
static const int GBNONBONDENERGY1_BLOCKS_LOAD             = 1;
static const int GBNONBONDENERGY2_BLOCKS_LOAD             = 1;
static const int PMENONBONDENERGY_THREADS_PER_BLOCK       = 1024;
static const int PMENONBONDFORCES_THREADS_PER_BLOCK       = 1024;
static const int PMENONBONDENERGY_BLOCKS_MULTIPLIER       = 1;
static const int IPSNONBONDENERGY_THREADS_PER_BLOCK       = 1024;
static const int IPSNONBONDFORCES_THREADS_PER_BLOCK       = 1024;
static const int IPSNONBONDENERGY_BLOCKS_MULTIPLIER       = 1;
static const int TRANSPOSE_QMESH_THREADS_PER_BLOCK        = 64;
#else
static const int PMENONBONDFORCES_THREADS_PER_BLOCK       = 128;
static const int PMENONBONDENERGY_THREADS_PER_BLOCK       = 128;
static const int GBBORNRADII_THREADS_PER_BLOCK            = 64;
static const int GBNONBONDENERGY1_THREADS_PER_BLOCK       = 64;
static const int GBNONBONDENERGY2_THREADS_PER_BLOCK       = 64;
#if ((__CUDA_ARCH__ == 750) || (__CUDA_ARCH__ < 500))
static const int GBBORNRADII_BLOCKS_MULTIPLIER            = 16;
static const int GBNONBONDENERGY1_BLOCKS_MULTIPLIER       = 16;
static const int GBNONBONDENERGY2_BLOCKS_MULTIPLIER       = 16;
static const int GBNONBONDENERGY1_BLOCKS_LOAD             = 16;
static const int GBNONBONDENERGY2_BLOCKS_LOAD             = 16;
static const int PMENONBONDENERGY_BLOCKS_MULTIPLIER       = 8;
#else
static const int GBBORNRADII_BLOCKS_MULTIPLIER            = 24;
static const int GBNONBONDENERGY1_BLOCKS_MULTIPLIER       = 24;
static const int GBNONBONDENERGY2_BLOCKS_MULTIPLIER       = 24;
static const int GBNONBONDENERGY1_BLOCKS_LOAD             = 24;
static const int GBNONBONDENERGY2_BLOCKS_LOAD             = 24;
static const int PMENONBONDENERGY_BLOCKS_MULTIPLIER       = 10;
#endif
static const int IPSNONBONDENERGY_THREADS_PER_BLOCK       = 512;
static const int IPSNONBONDFORCES_THREADS_PER_BLOCK       = 512;
static const int IPSNONBONDENERGY_BLOCKS_MULTIPLIER       = 2;
static const int TRANSPOSE_QMESH_THREADS_PER_BLOCK        = 128;
#endif

static const int MAXMOLECULES                             = 1024;
static const int MAXPSMOLECULES                           = 2040;
static const int READ_SIZE                                = 128;
static const int GB_TEXTURE_WIDTH                         = 512;
#ifdef DPFP
static const int MAX_RANDOM_STEPS                         = 64;
static const int SMEM_MESH_SIZE                           = 5400;
#else
static const int MAX_RANDOM_STEPS                         = 128;
#endif

//---------------------------------------------------------------------------------------------
// PME data
//---------------------------------------------------------------------------------------------
#define PI_VAL 3.1415926535897932384626433832795
GL_CONST PMEDouble PI   = (PMEDouble)PI_VAL;
GL_CONST PMEFloat  PI_F = (PMEFloat)PI_VAL;

// Used for bounding allowable contents of nonbond cells
#define AVERAGE_ATOM_VOLUME_VAL 9.83
static const PMEDouble AVERAGE_ATOM_VOLUME = (PMEDouble)AVERAGE_ATOM_VOLUME_VAL;
static const PMEFloat AVERAGE_ATOM_VOLUMEF = (PMEFloat)AVERAGE_ATOM_VOLUME_VAL;

// Boltzmann's constant in internal units
#define KB_VAL (1.380658 * 6.0221367) / (4.184 * 1000.0)
GL_CONST PMEDouble KB   = (PMEDouble) KB_VAL;
GL_CONST PMEFloat  KB_F = (PMEFloat)  KB_VAL;

// Compilation flag: turns off sysmem shadowing of really large buffers
static const bool bShadowedOutputBuffers = false;

//---------------------------------------------------------------------------------------------
// Topological and neighbor list data types
//---------------------------------------------------------------------------------------------
struct listdata_rec {
  int offset;     //
  int cnt;        //
};

struct bond_rec {
  int atm_i;      // The i atom in the bond, per Amber conventions
  int atm_j;      // The j atom in the bond
  int parm_idx;   // Index in to the bonded term parameter arrays
  int bat_sc;     // Index for bat sc term (gpu unused)
};

struct angle_rec {
  int atm_i;      // The i atom in the angle, per Amber nomenclature
  int atm_j;      // The j atom in the angle
  int atm_k;      // The k atom in the angle
  int parm_idx;   // Index into the angle term parameter arrays
  int bat_sc;     // Index for bat sc term (gpu unused)
};

struct dihed_rec {
  int atm_i;      // The i atom in the dihedral angle, per Amber nomenclature
  int atm_j;      // The j atom in the diehdral
  int atm_k;      // The k atom in the dihedral
  int atm_l;      // The l atom in the dihedral
  int parm_idx;   // Index into the dihedral parameter arrays
  int bat_sc;     // Index for bat sc term (gpu unused)
};

struct angle_ub_rec {
  int atm_i;      // The i atom of a Urey-Bradley potential (there are only two, i and j)
  int atm_j;      // The j atom of a Urey-Bradley potential (arrangement is i-X-j)
  int parm_idx;   // Index into the Urey-Bradley parameter arrays
};

struct dihed_imp_rec {
  int atm_i;      // The i atom in the improper torsion angle, per Amber nomenclature
  int atm_j;      // The j atom in the diehdral
  int atm_k;      // The k atom in the dihedral
  int atm_l;      // The l atom in the dihedral
  int parm_idx;   // Index into the improper torsion parameter arrays
};
struct cmap_rec {
  int atm_i;      // The i atom in the CMAP pair of torsion angles
  int atm_j;      // The j atom in the CMAP
  int atm_k;      // The k atom in the CMAP
  int atm_l;      // The l atom in the CMAP
  int atm_m;      // The m atom in the CMAP
  int parm_idx;   // Index into the CMAP parameter arrays
};

struct shake_bond_rec {
  int atm_i;      // Atom i of the bond being SHAKEn
  int atm_j;      // Atom j of the bond being SHAKEn
  double parm;    // Target length of the bond being SHAKEn (just the length
                  //   is needed, so no need to index a parameter set)
};

struct gb_pot_ene_rec {
  double total;       // Total (potential) energy
  double vdw_tot;     // Total van-der Waals (that means Lennard-Jones) energy
  double elec_tot;    // Total electrostatic energy
  double gb;          // Generalized Born energy
  double surf;        // Energy related to surface area
  double bond;        // Bond strain energy
  double angle;       // Angle strain energy
  double dihedral;    // Dihedral / torsion energy (does not include improper torsions)
  double vdw_14;      // 1-4 Lennard-Jones energy
  double elec_14;     // 1-4 electrostatic energy
  double restraint;   // Restraint strain energy
  double angle_ub;    // Urey-Bradley angle energy
  double imp;         // Inproper torsion energy
  double cmap;        // CMAP total energy
  double dvdl;        // dV/dL calculated in Thermodynamic Integration
  double amd_boost;   // Energy owing to accelerated MD boost potentials
  double gamd_boost;  // Energy owning to Gaussian-accelerated MD boost potentials
  double emap;        // E-map restraint energy
  double nfe;         // NFE energy contribution
};

struct pme_pot_ene_rec {
  double total;          // Total potential energy of the system
  double vdw_tot;        // Total of direct and reciprocal space Lennard-Jones energies
  double vdw_dir;        // Direct-space Lennard-Jones energies
  double vdw_recip;      // Reciprocal space computed Lennard-Jones energies (really this is
                         //   the energy resulting from the homogeneity approximation, NOT
                         //   Lennard-Jones inverse r6 interactions being performed on a mesh,
                         //   and it is calculated by the Fortran code as part of the setup on
                         //   the CPU, then merely remembered by the GPU)
  double elec_tot;       // Total of direct space, reciprocal space, nb_adjust adjustment
                         //   energy, and self energy of each charge on the mesh
  double elec_dir;       // Direct space electrostatic energy
  double elec_recip;     // Reciprocal space electrostatic energy
  double elec_nb_adjust; // Non-bonded adjustments needed to remove 1-2 and 1-3 exclusions from
                         //   the electrostatics, particularly the energy coming off the mesh
  double elec_self;      // Self energy of each charge in the PME sum
  double hbond;          // Energy of hydrogen-bonding terms
  double bond;           // Bond strain energy
  double angle;          // Angle strain energy
  double dihedral;       // Dihedral / torsion energy
  double vdw_14;         // 1-4 Lennard-Jones interaction energy
  double elec_14;        // 1-4 electrostatic interaction energy
  double restraint;      // Restraint energy
  double angle_ub;       // Urey-Bradley angle energy
  double imp;            // Improper torsion energy
  double cmap;           // CMAP total energy
  double amd_boost;      // Energy due to accelerated MD boost potentials
  double gamd_boost;     // Energy due to Gaussian-aMD boost potentials
  double emap;
  double efield;
  double nfe;
};

struct afe_gpu_sc_ene_rec {
    double dvdl;
    double bond_R1;
    double bond_R2;
    double angle_R1;
    double angle_R2;
    double dihedral_R1;
    double dihedral_R2;
    double sc_res_dist_R1;
    double sc_res_dist_R2;
    double sc_res_ang_R1;
    double sc_res_ang_R2;
    double sc_res_tors_R1;
    double sc_res_tors_R2;
    double vdw_dir_R1;
    double vdw_dir_R2;
    double elec_dir_R1;
    double elec_dir_R2;
    double vdw_14_R1;
    double vdw_14_R2;
    double elec_14_R1;
    double elec_14_R2;
    double vdw_der_R1;
    double vdw_der_R2;
    double elec_der_R1;
    double elec_der_R2;
};

struct NTPData {
  double last_recip[9];  // The previous matrix for taking coordinates into fractional space
  double recip[9];       // The matrix to take coordinates into fractional space
  double ucell[9];       // The matrix to transform fractional coordinates into real space
  PMEFloat recipf[9];    // Single-precision variants of
  PMEFloat ucellf[9];    //   recip and ucell

  // Shortcuts for comparisons with the cutoff + skin boundary
  PMEFloat one_half_nonbond_skin_squared;
  PMEFloat cutPlusSkin2;
};

struct ep_frame_rec {
  int extra_pnt[2];  // Atom indices of the extra points in the topology
  int ep_cnt;        // Number of extra points associated with this record
  int type;          // Frame type for the extra point
  int parent_atm;    // Parent atom to which the extra point attaches
  int frame_atm1;    //
  int frame_atm2;    // Frame atoms that determine the extra point location(s)
  int frame_atm3;    //
};

#ifdef use_DPFP
#  define LSCALE  (1ll << 50)
#else
#  define LSCALE  (1 << 24)
#endif
#define XSCALE  (1ll << 12)
#define ESCALE  (1ll << 30)
#define VSCALE  (1ll << 40)
#define FSCALE  (1ll << 40)
#define NBSCALE_1 (1 << 20)    // Used in SPFP mode: floating point forces of magnitgude at
                               //   least 1 kcal/mol-A will see their mantissas fully
                               //   incorporated into the resulting int.
#define NBSHIFT_2  20          // Shift the accumulated integer representation of non-bonded
                               //   direct space forces to the left by this amount to get a
                               //   result that is consistent with units of SPFP forces
                               //   elsewhere in the code.
#define NRGSHIFT_2 10          // Shift the accumulated and reduced energy to the left by this
                               //   amount to put it on the same scale as other energies in the
                               //   code.
#define DFSCALE (1ll << 44)
GL_CONST PMEDouble LATTICESCALE               = (PMEDouble)LSCALE;
GL_CONST PMEFloat  LATTICESCALEF              = (PMEFloat)LSCALE;
GL_CONST PMEDouble ONEOVERLATTICESCALE        = (PMEDouble)1.0 / (PMEDouble)LSCALE;
GL_CONST PMEFloat  ONEOVERLATTICESCALEF       = (PMEDouble)1.0 / (PMEDouble)LSCALE;
GL_CONST PMEDouble PMEGRID_SCALE              = (PMEDouble)(1ll << 54);
GL_CONST PMEFloat  PMEGRID_SCALEF             = (PMEFloat)(1 << 22);
GL_CONST PMEDouble INV_PMEGRID_SCALE          = (PMEDouble)1.0 / (PMEDouble)(1ll << 54);
GL_CONST PMEFloat  INV_PMEGRID_SCALEF         = (PMEFloat)1.0 / (PMEFloat)(1 << 22);
GL_CONST PMEDouble ENERGYSCALE                = (PMEDouble)ESCALE;
GL_CONST PMEFloat  ENERGYSCALEF               = (PMEFloat)ESCALE;
GL_CONST PMEDouble XRAYSCALE                  = (PMEDouble)XSCALE;
GL_CONST PMEFloat  XRAYSCALEF                 = (PMEFloat)XSCALE;
GL_CONST PMEDouble ONEOVERXRAYSCALE           = (PMEDouble)1.0 / (PMEDouble)XSCALE;
GL_CONST PMEFloat  HALF_ENERGYSCALEF          = (PMEFloat)ESCALE * (PMEFloat)0.5;
GL_CONST PMEDouble ONEOVERENERGYSCALE         = (PMEDouble)1.0 / (PMEDouble)ESCALE;
GL_CONST PMEDouble ONEOVERENERGYSCALESQUARED  = (PMEDouble)1.0 / ((PMEDouble)ESCALE *
                                                                      (PMEDouble)ESCALE);
GL_CONST PMEDouble FORCESCALE                 = (PMEDouble)FSCALE;
GL_CONST PMEFloat  FORCESCALEF                = (PMEFloat)FSCALE;
GL_CONST PMEFloat  HALF_FORCESCALEF           = (PMEFloat)FSCALE * (PMEFloat)0.5;
GL_CONST PMEFloat  FORCESCALEF_NB1            = (PMEFloat)NBSCALE_1;
GL_CONST PMEFloat  ENERGYSCALEF_NB1           = (PMEFloat)NBSCALE_1;
GL_CONST PMEFloat  HALF_ENERGYSCALEF_NB1      = (PMEFloat)NBSCALE_1 * (PMEFloat)0.5;
GL_CONST PMEDouble ONEOVERFORCESCALE          = (PMEDouble)1.0 / (PMEDouble)FSCALE;
GL_CONST PMEFloat  ONEOVERFORCESCALEF         = (PMEDouble)1.0 / (PMEDouble)FSCALE;
GL_CONST PMEDouble ONEOVERFORCESCALESQUARED   = (PMEDouble)1.0 / ((PMEDouble)FSCALE *
                                                                      (PMEDouble)FSCALE);
GL_CONST PMEDouble DIHEDRALFORCESCALE         = (PMEDouble)DFSCALE;
GL_CONST PMEFloat  DIHEDRALFORCESCALEF        = (PMEFloat)DFSCALE;
GL_CONST PMEDouble ONEOVERDIHEDRALFORCESCALE  = (PMEDouble)1.0 / (PMEDouble)DFSCALE;
GL_CONST PMEFloat  ONEOVERDIHEDRALFORCESCALEF = (PMEDouble)1.0 / (PMEDouble)DFSCALE;

GL_CONST PMEDouble VIRIALSCALE = (PMEDouble)VSCALE;
GL_CONST PMEFloat VIRIALSCALEF = (PMEFloat)VSCALE;
GL_CONST PMEDouble ONEOVERVIRIALSCALE = (PMEDouble)1.0 / (PMEDouble)VSCALE;
GL_CONST PMEFloat  ONEOVERVIRIALSCALEF = (PMEDouble)1.0 / (PMEDouble)VSCALE;

GL_CONST unsigned int REAL32_NONZERO_MASK           = ((1u << 31) - 1);
GL_CONST unsigned long long int REAL64_NONZERO_MASK = ((1ull << 63) - 1);

#ifdef use_DPFP
#define eScale ENERGYSCALE 
#define fScale FORCESCALE 
#define vScale VIRIALSCALE 
#define ftoi llrint
#else
#define eScale ENERGYSCALEF 
#define fScale FORCESCALEF 
#define vScale VIRIALSCALEF 
#define ftoi fast_llrintf
#endif
struct KineticEnergyRecord {
  PMEFloat EKE;
  PMEFloat EKPH;
  PMEFloat EKPBS;
};

union KineticEnergy {
  struct KineticEnergyRecord KE;  // This repeated use of unions to make things readable as
  PMEFloat array[3];              //   arrays is sensible and conducive to clean code.
};

// Alchemical free energy kinetic energy storage
struct AFEKineticEnergyRecord
{
  PMEFloat TI_EKER1;
  PMEFloat TI_EKER2;
  PMEFloat TI_EKPBSR1;
  PMEFloat TI_EKPBSR2;
  PMEFloat TI_EKPHR1;
  PMEFloat TI_EKPHR2;
  PMEFloat TI_SC_EKER1;
  PMEFloat TI_SC_EKER2;
};

union AFEKineticEnergy
{
  // It might be possible to conditionally define this to only have 4 terms for Langevin
  // dynamics.  However, this should only happen if it's a huge performance boost.  That
  // would just amount to making two unions, and only using the one we want for the given
  // conditions.
  struct AFEKineticEnergyRecord AFEKE;
  PMEFloat array[8];
};

//---------------------------------------------------------------------------------------------
// Error and operation reporting macros.  Debugging flags below
// are commented out for normal compilation.
//---------------------------------------------------------------------------------------------
//#define GVERBOSE
//#define MEMTRACKING
//#define SYNCHRONOUS
#ifdef GVERBOSE
#ifndef MEMTRACKING
#define MEMTRACKING
#endif
#ifdef MPI
#define PRINTMETHOD(name) \
{ \
    printf("Method: %s on node %d\n", name, gpu->gpuID); \
    fflush(stdout); \
}

#ifdef SYNCHRONOUS
#define LAUNCHERROR(s) \
    { \
        printf("Launched %s on node %d\n", s, gpu->gpuID); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
        cudaDeviceSynchronize(); \
    }
#else
#define LAUNCHERROR(s) \
    { \
        printf("Launched %s on node %d\n", s, gpu->gpuID); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
    }
#endif
#define LAUNCHERROR_BLOCKING(s) \
    { \
        printf("Launched %s on node %d\n", s, gpu->gpuID); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
        cudaDeviceSynchronize(); \
    }
#define LAUNCHERROR_NONBLOCKING(s) \
    { \
        printf("Launched %s on node %d\n", s, gpu->gpuID); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
    }
#else
#define PRINTMETHOD(name) \
{ \
    printf("Method: %s\n", name); \
    fflush(stdout); \
}
#ifdef SYNCHRONOUS
#define LAUNCHERROR(s) \
    { \
        printf("Launched %s\n", s); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
        cudaDeviceSynchronize(); \
    }
#else
#define LAUNCHERROR(s) \
    { \
        printf("Launched %s\n", s); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
    }
#endif
#define LAUNCHERROR_BLOCKING(s) \
    { \
        printf("Launched %s\n", s); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
        cudaDeviceSynchronize(); \
    }
#define LAUNCHERROR_NONBLOCKING(s) \
    { \
        printf("Launched %s\n", s); \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
    }
#endif
#else
#define PRINTMETHOD(name)
#ifdef SYNCHRONOUS
#define LAUNCHERROR(s) \
    { \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
        cudaDeviceSynchronize(); \
    }
#else
#define LAUNCHERROR(s) \
    { \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
    }
#endif
#define LAUNCHERROR_BLOCKING(s) \
    { \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
        cudaDeviceSynchronize(); \
    }
#define LAUNCHERROR_NONBLOCKING(s) \
    { \
        cudaError_t status = cudaGetLastError(); \
        if (status != cudaSuccess) { \
            printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
            gpu_shutdown_(); \
            exit(-1); \
        } \
    }
#endif

#define RTERROR(status, s) \
    if (status != cudaSuccess) { \
        printf("%s %s\n", s, cudaGetErrorString(status)); \
        cudaDeviceReset(); \
        exit(-1); \
    }


#endif
