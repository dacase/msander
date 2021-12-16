#ifndef _BASE_GPU_CONTEXT
#define _BASE_GPU_CONTEXT

#include "gputypes.h"
#include "gpuBuffer.h"
/*  #include "simulationConst.h"  */
#include "gpuMemoryInfo.h"

class base_gpuContext {

public:
  virtual void Init();

  base_gpuContext();
  virtual ~base_gpuContext();

public:

  simulationConst sim;         // All simulation data destined for GPU constant RAM
  
  SM_VERSION sm_version;

  // Memory parameters
  bool bECCSupport;            // Flag for ECC support to detect Tesla versus consumer GPU
  bool bCanMapHostMemory;      // Flag for pinned memory support
  aligned_lli totalMemory;     // Total memory on GPU
  aligned_lli totalCPUMemory;  // Approximate total allocated CPU memory
  aligned_lli totalGPUMemory;  // Approximate total allocated CPU memory

  // AMBER parameters
  int ntt;              // AMBER Thermostat setting
  int ntb;              // AMBER PBC setting
  int ips;              // AMBER IPS setting
  int ntc;              // AMBER SHAKE setting
  int imin;             // AMBER minimization versus MD setting
  int ntf;              // AMBER force field setting
  int ntpr;             // AMBER status output interval
  int ntwe;             // AMBER energy output interval
  int ntr;              // AMBER Restraint flags
  int gbsa;             // AMBER GBSA active flag
  int step;             // Current step
  bool bCPURandoms;     // Generate random numbers on CPU instead of GPU

  // Atom stuff
  GpuBuffer<double>*         pbAtom;              // Atom coordinates
  GpuBuffer<PMEFloat2>*      pbAtomXYSP;          // Single Precision Atom X and Y coordinates
  GpuBuffer<PMEFloat>*       pbAtomZSP;           // Single Precision Atom Z coordinate
  GpuBuffer<PMEFloat2>*      pbAtomSigEps;        // Atom nonbond parameters
  GpuBuffer<unsigned int>*   pbAtomLJID;          // Atom Lennard Jone index
  GpuBuffer<PMEFloat>*       pbAtomS;             // Atom scaled Born Radius
  GpuBuffer<PMEFloat>*       pbAtomRBorn;         // Atom Born Radius
  GpuBuffer<double>*         pbAtomCharge;        // Atom charges
  GpuBuffer<PMEFloat>*       pbAtomChargeSP;      // Single precision atom charges
  GpuBuffer<PMEFloat2>*      pbAtomChargeSPLJID;  // Single precision atom charges and
                                                  //   Lennard Jones atom types
  GpuBuffer<double>*         pbAtomMass;          // Atom masses
  GpuBuffer<PMEDouble>*      pbReff;              // Effective Born Radius
  GpuBuffer<PMEFloat>*       pbReffSP;            // Single Precision Effective Born Radius
  GpuBuffer<PMEFloat>*       pbTemp7;             // Single Precision Born Force
  GpuBuffer<PMEFloat>*       pbPsi;               // GB intermediate psi value
  GpuBuffer<int>*            pbTIRegion;          // Softcore TI Region
  GpuBuffer<int>*            pbTILinearAtmID;     // AFE linear atom IDs for use in vector
                                                  //   exchange
  GpuBuffer<PMEDouble>*       pbBarLambda;        // Bar Lambda values for AFE
  GpuBuffer<unsigned long long int>* pbBarTot;    // Holds final mc bar values for AFE
  bool bCalculateLocalForces;                     // Flags to signal calculating local forces,
  bool bCalculateReciprocalSum;                   //   reciprocal sum, and direct sum forces
  bool bCalculateDirectSum;                       //   on this process
  GpuBuffer<PMEAccumulator>* pbPeerAccumulator;   // Remote data accumulator buffer
  PMEAccumulator** pPeerAccumulatorList;          // List of accumulator buffer addresses
  cudaIpcMemHandle_t* pPeerAccumulatorMemHandle;  // Memory handles for each peer accumulator
  GpuBuffer<double>*         pbVel;               // Atom velocities
  GpuBuffer<double>*         pbLVel;              // Atom last velocities
  GpuBuffer<PMEFloat>*       pbCenter;            // Atom recentering extrema

  // Lennard-Jones Stuff
  GpuBuffer<PMEFloat2>*      pbLJTerm;            // Pointer to LJ terms buffer
  GpuBuffer<int2>*           pbR4sources;         // Integers for atom ID numbers of 1/r4
                                                  //   sources followed by Lennard-Jones type
                                                  //   numbers for indexing potentials

  // GB Data
  GpuBuffer<PMEFloat2>*      pbNeckMaxValPos;     // GB Born Radii and Energy correction data
  GpuBuffer<PMEFloat>*       pbGBAlphaBetaGamma;  // GB per atom Alpha, Beta, and Gamma

  // GBSA3 Data pwsasa
  GpuBuffer<PMEFloat>*       pbGBSASigEpsiRadiMax; // GBSA 3 parameters

  // PME data
  cufftHandle forwardPlan;                        // PME FFT plan for cuFFT
  cufftHandle backwardPlan;                       // PME IFFT plan for cuFFT
#ifdef use_DPFP
  GpuBuffer<long long int>*  pblliXYZ_q;          // DPFP long long int PME charge grid buffer
#else
  GpuBuffer<int>*            pblliXYZ_q;          // SPFP int PME charge grid buffer
#endif
  GpuBuffer<PMEFloat>*       pbXYZ_q;             // PME charge grid/buffer
  GpuBuffer<PMEComplex>*     pbXYZ_qt;            // FFTed PME charge grid
  GpuBuffer<PMEFloat>*       pbPrefac;            // PME FFT pre-factors
  GpuBuffer<PMEFloat>*       pbFract;             // PME fractional coordinates
  GpuBuffer<int2>*           pbTileBoundary;      // PME charge interpolation box boundaries
                                                  //   (min/max)
#ifndef use_DPFP
  GpuBuffer<PMEFloat4>*      ErfcCoeffsTable;     // d/dr [erfc(a * r)/r] spline coefficients
  GpuBuffer<PMEFloat4>*      pbErfcCoeffsTable;   // d/dr [erfc(a * r)/r] spline coefficients
#endif

  // Self energy values, calculated on host
  PMEDouble ee_plasma;        // Component of PME self energy
  PMEDouble self_energy;      // Total self energy
  PMEDouble vdw_recip;        // Reciprocal vdw correction energy

  // Neighbor list stuff
  bool bNeighborList;             // Neighbor list activation flag
  bool bNewNeighborList;          // Newly generated neighbor list activation flag
  bool bNeedNewNeighborList;      // Need newly generated neighbor list activation flag
  bool bOddNLCells;               // Flag to determine whether charge grid needs 8 or 16
                                  //   output buffers
  bool bSmallBox;                 // Flag to determine whether any nonbond cell count is
                                  //   1 or 2 on any axis
  bool bPencilYShift;             // Flag to indicate that the pencils are staggered along
                                  //   the Y-axis as a function of Z, or if the natural tilt
                                  //   of the grid lines up with the unit cell
  unsigned int neighborListBits;  // Number of bits in cell list
  size_t sortTempBytes;           // Temporary bytes needed for CUB sort
  GpuBuffer<char>* pbSortTemp;    // CUB GPU temporary data buffer
  PMEFloat nonbond_skin;          // Nonbond skin thickness
  GpuBuffer<unsigned int>* pbImageIndex;      // Image indexing data
  GpuBuffer<unsigned int>* pbSubImageLookup;  // Sub-image indexing for all atoms.  This
                                              //   indexing references data in pbSubImageSP,
                                              //   NOT the pointers pQQSubImageSP and
                                              //   pLJSubImageSP (see above, in cudaSimulation)
                                              //   that point to each half of it!
  GpuBuffer<PMEFloat2>*    pbAtomXYSaveSP;    // Saved neighbor list coordinates
  GpuBuffer<PMEFloat>*     pbAtomZSaveSP;     // Saved neighbor list coordinates
  GpuBuffer<double>*       pbImage;           // Image coordinates
  GpuBuffer<double>*       pbImageVel;        // Image velocities
  GpuBuffer<double>*       pbImageLVel;       // Image last velocities
  GpuBuffer<double>*       pbImageMass;       // Image masses
  GpuBuffer<double>*       pbImageCharge;     // Image charges
  GpuBuffer<PMEFloat2>*    pbImageSigEps;     // Image sigma/epsilon data for nonbond
                                              //   interactions
  GpuBuffer<unsigned int>* pbImageLJID;       // Image Lennard Jones index for nonbond
                                              //   interactions
  GpuBuffer<int>*          pbImageTIRegion;          // Softcore TI region
  GpuBuffer<int>*          pbImageTILinearAtmID;     // atom indexes for linear scaling atoms
                                                     // for syncing coords between steps
  GpuBuffer<unsigned int>* pbUpdateIndex;            // index update during neighborlist update
                                                     // for use with syncing linear vectors
  GpuBuffer<unsigned int>* pbImageCellID;            // Image per atom nonbond cell ID
  GpuBuffer<unsigned int>* pbNLExclusionList;        // Raw exclusion list
  GpuBuffer<uint2>*        pbNLExclusionStartCount;  // Per atom exclusion bounds
  GpuBuffer<uint2>*        pbNLNonbondCellStartEnd;  // Nonbond cells start and end
  GpuBuffer<bool>*         pbNLbSkinTestFail;        // Skin test result buffer
  GpuBuffer<unsigned int>* pbNLCellHash;             // Spatial ordering hash for within cells
  GpuBuffer<NLRecord>*     pbNLRecord;               // Pointer to neighbor list records
  GpuBuffer<NLEntry>*      pbNLEntry;                // Active neighbor list
  GpuBuffer<unsigned int>* pbNLAtomList;             // Pointer to atom list and exclusion data
  GpuBuffer<unsigned int>* pbNLTotalOffset;          // Current atom list offset
  GpuBuffer<unsigned int>* pbFrcBlkCounters;         // Progress counters for various force
                                                     //   and neighbor list calculations
  GpuBuffer<unsigned int>* pbNLEntries;              // Total nonbond entries
  GpuBuffer<unsigned int>* pbBNLExclusionBuffer;     // Per-warp GMEM exclusion buffer

  // Bonded term work units
  bondwork* bondWorkUnitRecord;                 // A record of all bond work units in their
                                                //   verbose form (prior to conversion into
                                                //   the carefully ordered arrays kept by the
                                                //   GPU) made available so that the work units
                                                //   can be edited during the simulation.
  GpuBuffer<unsigned int>* pbBondWorkUnitUINT;  // Unsigned integer data for bonded work units
  GpuBuffer<PMEDouble2>* pbBondWorkUnitDBL2;    // Double2 data for bonded work units
  GpuBuffer<PMEFloat>* pbBondWorkUnitPFLOAT;    // PMEFloat data for bonded work units (like
                                                //   all PMEFloat data, this will be float32
                                                //   in SP*P but float64 in DPFP mode)
  GpuBuffer<PMEFloat2>* pbBondWorkUnitPFLOAT2;  // PMEFloat2 data for bonded work units
                                                //   (this holds dihedral and non-bonded 1:4
                                                //   interaction parameters)
  unsigned int *ptrBwuInstructions;             //
  unsigned int *ptrBwuBondID;                   //
  unsigned int *ptrBwuAnglID;                   //
  unsigned int *ptrBwuDiheID;                   //
  unsigned int *ptrBwuCmapID;                   //
  unsigned int *ptrBwuQQxcID;                   // 
  unsigned int *ptrBwuNB14ID;                   //
  unsigned int *ptrBwuNMR2ID;                   //
  unsigned int *ptrBwuNMR3ID;                   //
  unsigned int *ptrBwuNMR4ID;                   //
  unsigned int *ptrBwuUreyID;                   //
  unsigned int *ptrBwuCimpID;                   //
  unsigned int *ptrBwuCnstID;                   //
  unsigned int *ptrBwuCnstUpdateIdx;            // Pointers into the GpuBuffers listed above,
  unsigned int *ptrBwuBondStatus;               //   allowing easy access to host data in the
  unsigned int *ptrBwuAnglStatus;               //   event that it needs to be edited.  See
  unsigned int *ptrBwuDiheStatus;               //   also the pointers with similar names in
  unsigned int *ptrBwuCmapStatus;               //   the cudaSimulation struct, referencing
  unsigned int *ptrBwuQQxcStatus;               //   the equivalent data on the device.
  unsigned int *ptrBwuNB14Status;               //   Uploads and downloads cannot proceed with
  unsigned int *ptrBwuNMR2Status;               //   these pointers alone if the GpuBuffer
  unsigned int *ptrBwuNMR3Status;               //   method is to be used; the entire array
  unsigned int *ptrBwuNMR4Status;               //   must make the trip.  May be worthwhile to
  unsigned int *ptrBwuUreyStatus;               //   write data shuttle routines making use of
  unsigned int *ptrBwuCimpStatus;               //   these arrays.
  unsigned int *ptrBwuCnstStatus;               //
  PMEDouble2 *ptrBwuBond;                       //
  PMEDouble2 *ptrBwuAngl;                       //
  PMEDouble2 *ptrBwuNMR2;                       //
  PMEDouble2 *ptrBwuNMR3;                       //
  PMEDouble2 *ptrBwuNMR4;                       //
  PMEDouble2 *ptrBwuUrey;                       //
  PMEDouble2 *ptrBwuCimp;                       //
  PMEDouble2 *ptrBwuCnst;                       //
  PMEFloat *ptrBwuCharges;                      //
  PMEFloat *ptrBwuDihe3;                        //
  PMEFloat *ptrBwuEEnb14;                       //
  PMEFloat2 *ptrBwuDihe12;                      //
  PMEFloat2 *ptrBwuLJnb14;                      //

  // Atom shuttling
  GpuBuffer<double>* pbDataShuttle;             // Data transfer buffer for abstracted
                                                //   atom information (coordinates,
                                                //   forces, etc.)
  GpuBuffer<int>* pbShuttleTickets;             // List of atoms that will be on the
                                                //   shuttle to and from the CPU

  // NEB info
  GpuBuffer<int>* pbRMSMask;                    // Mask of atoms used for NEB force calculation.
  GpuBuffer<int>* pbFitMask;                    // Mask of atoms used for NEB structure fitting.
  GpuBuffer<int>* pbAtmIdx;
  GpuBuffer<PMEDouble>* pbNEBEnergyAll;         // Energy of the replicas in NEB.
  GpuBuffer<double>* pbTangents;                // Tangents to the path.
  GpuBuffer<double>* pbSpringForce;             // Artificial spring forces.
  GpuBuffer<double>* pbNEBForce;                // NEB forces.
  GpuBuffer<double>* pbNextDataShuttle;         // Next neighboring replica's coordinate.
  GpuBuffer<double>* pbPrevDataShuttle;         // Previous neighboring replica's coordinate.
  GpuBuffer<double>* pbKabschCOM;
  GpuBuffer<double>* pbtotFitMass;
  GpuBuffer<double>* pbDataSPR;
  GpuBuffer<double>* pbRotAtm;

  // Remapped bonded interactions and Shake constraint atom IDs
  GpuBuffer<int2>*   pbImageNMRCOMDistanceID;        // Remapped NMR COM Distance i, j
  GpuBuffer<int2>*   pbImageNMRCOMDistanceCOM;       // Remapped NMR COM DistanceCOM range
                                                     //   .x to .y
  GpuBuffer<int2>*   pbImageNMRCOMDistanceCOMGrp;    // Remapped NMR COM DistanceCOMGrp range
                                                     //   .x to .y
  GpuBuffer<int2>*   pbImageNMRr6avDistanceID;       // Remapped NMR r6av Distance i, j
  GpuBuffer<int2>*   pbImageNMRr6avDistancer6av;     // Remapped NMR r6av Distancer6av range
                                                     //   .x to .y
  GpuBuffer<int2>*   pbImageNMRr6avDistancer6avGrp;  // Remapped NMR r6av Distancer6avGrp
                                                     //   range .x to .y
  GpuBuffer<int2>*   pbImageNMRCOMAngleID1;          // Remapped NMR COM Angle i, j
  GpuBuffer<int>*    pbImageNMRCOMAngleID2;          // Remapped NMR COM Angle k
  GpuBuffer<int2>*   pbImageNMRCOMAngleCOM;          // Remapped NMR COM AngleCOM range
                                                     //   .x to .y
  GpuBuffer<int2>*   pbImageNMRCOMAngleCOMGrp;       // Remapped NMR COM AngleCOMGrp range
                                                     //   .x to .y
  GpuBuffer<int4>*   pbImageNMRCOMTorsionID1;        // Remapped NMR COM Torsion i, j, k , l
  GpuBuffer<int2>*   pbImageNMRCOMTorsionCOM;        // Remapped NMR COM TorsionCOM range
                                                     //   .x to .y
  GpuBuffer<int2>*   pbImageNMRCOMTorsionCOMGrp;     // Remapped NMR COM TorsionCOMGrp range
                                                     //   .x to .y
 
  GpuBuffer<int4>*   pbImageShakeID;                 // Remapped traditional SHAKE IDs
  GpuBuffer<int4>*   pbImageFastShakeID;             // Remapped H2O (Fast) SHAKE IDs
  GpuBuffer<int>*    pbImageSlowShakeID1;            // Remapped Central atom of XH4 (Slow)
                                                     //   Shake constraint
  GpuBuffer<int4>*   pbImageSlowShakeID2;            // Remapped XH4 (Slow) SHAKE constraint
                                                     //   hydrogens
  GpuBuffer<int4>*   pbImageSolventAtomID;           // Remapped solvent molecules/ions
  GpuBuffer<int>*    pbImageSoluteAtomID;            // Remapped solute atoms

  // NTP molecule data
  int maxSoluteMolecules;     // Maximum solute molecules for most NTP kernels
  int maxPSSoluteMolecules;   // Maximum solute molecules for NTP pressure scaling kernels
  GpuBuffer<int>*             pbSoluteAtomID;   // Solute per atom ID
  GpuBuffer<PMEDouble>*       pbSoluteAtomMass; // Solute atom masses
  GpuBuffer<PMEDouble>*       pbSolute;         // Current and last centers of mass, kinetic
                                                //   energy and inverse mass for each solute
                                                //   molecule
  GpuBuffer<PMEUllInt>*       pbUllSolute;      // Current COM and EKCOM in integer form
  GpuBuffer<int4>*            pbSolventAtomID;  // List of solvent molecules/ions of 4 or
                                                //   fewer atoms
  GpuBuffer<PMEDouble>*       pbSolvent;        // Last centers of mass, atom masses and
                                                //   inverse mass for each solvent molecule
  GpuBuffer<NTPData>*         pbNTPData;        // NTP mutable values
  GpuBuffer<int>*             pbAFEMolType;     // Molecule type for COM calculation for
                                                //   alchemical free energy
  GpuBuffer<int>*             pbAFEMolPartner;  // Partner molecule ID for COM calculation
                                                //   for alchemical free energy

  // NTP constraint molecule data
  GpuBuffer<PMEDouble>*    pbConstraintAtomX;   // Original atom x for constraint
  GpuBuffer<PMEDouble>*    pbConstraintAtomY;   // Original atom y for constraint
  GpuBuffer<PMEDouble>*    pbConstraintAtomZ;   // Original atom z for constraint
  GpuBuffer<PMEDouble>*    pbConstraintCOMX;    // Original X, Y, and Z COM for constraints in
  GpuBuffer<PMEDouble>*    pbConstraintCOMY;    //   fractional coordinates
  GpuBuffer<PMEDouble>*    pbConstraintCOMZ;    //

  // Radix sort data
  GpuBuffer<unsigned int>* pbRadixBlockSum;     // Radix sums per thread group
  GpuBuffer<unsigned int>* pbRadixCounter;      // Radix counters per thread group for scatter
  GpuBuffer<unsigned int>* pbRadixSum;          // Pointer to indivudal thread block radix sums

  // Random number stuff
  unsigned int randomCounter;   // Counter for triggering RNG
  curandGenerator_t RNG;        // CURAND RNG
  GpuBuffer<double>* pbRandom;  // Pointer to overall RNG buffer

  // Bonded interactions
  bool bLocalInteractions;   // Flag indicating presence or absence of local interactions
  bool bCharmmInteractions;  // Flag indicating presence or absence of CHARMM interactions
  GpuBuffer<PMEDouble2>* pbBond;           // Bond Kr, Req
  GpuBuffer<int2>*       pbBondID;         // Bond i, j
  int2*                  pQQxcID;          // Electrostatic exclusion i, j
  GpuBuffer<PMEDouble2>* pbBondAngle;      // Bond Angle Ka, Aeq
  GpuBuffer<int2>*       pbBondAngleID1;   // Bond Angle i, j
  GpuBuffer<int>*        pbBondAngleID2;   // Bond Angle k;
  GpuBuffer<PMEFloat2>*  pbDihedral1;      // Dihedral gmul, pn
  GpuBuffer<PMEFloat2>*  pbDihedral2;      // Dihedral pk, gamc
  GpuBuffer<PMEFloat>*   pbDihedral3;      // Dihedral gams
  GpuBuffer<int4>*       pbDihedralID1;    // Dihedral i, j, k, l
  GpuBuffer<PMEDouble2>* pbNb141;          // 1-4 nonbond scee, scnb0
  GpuBuffer<PMEDouble2>* pbNb142;          // 1-4 nonbond cn1, cn2
  GpuBuffer<int2>*       pbNb14ID;         // 1-4 nonbond i, j
  GpuBuffer<PMEDouble2>* pbConstraint1;    // Constraint weight and xc
  GpuBuffer<PMEDouble2>* pbConstraint2;    // Constraint yc and zc
  GpuBuffer<int>*        pbConstraintID;   // Atom constraint ID
  GpuBuffer<PMEDouble2>* pbUBAngle;        // Urey Bradley Angle rk, r0
  GpuBuffer<int2>*       pbUBAngleID;      // Urey Bradley Angle i, j
  GpuBuffer<PMEDouble2>* pbImpDihedral;    // Improper Dihedral pk, phase
  GpuBuffer<int4>*       pbImpDihedralID1; // Improper Dihedral i, j, k, l
  GpuBuffer<int4>*       pbCmapID1;        // Cmap i, j, k, l
  GpuBuffer<int>*        pbCmapID2;        // Cmap m
  GpuBuffer<int>*        pbCmapType;       // Cmap type
  GpuBuffer<PMEFloat4>*  pbCmapEnergy;     // Pointer to Cmap LUT data (E, dPhi, dPsi,
                                           //   dPhi_dPsi)

  // NMR stuff
  bool bNMRInteractions;  // Flag indicating presence or absence of NMR interactions
  int  NMRnstep;          // Imported NMR variable for time-dependent restraints
  GpuBuffer<double>*     pbNMRJarData;             // Jarzynski accumulated work data
  GpuBuffer<int2>*       pbNMRDistanceID;          // NMR distance i, j
  GpuBuffer<PMEDouble2>* pbNMRDistanceR1R2;        // NMR distance computed r1, r2
  GpuBuffer<PMEDouble2>* pbNMRDistanceR3R4;        // NMR distance computed r3, r4
  GpuBuffer<PMEDouble2>* pbNMRDistanceK2K3;        // NMR distance computed k2, k3
  GpuBuffer<PMEDouble>*  pbNMRDistanceK4;          // NMR distance computed k4
  GpuBuffer<PMEDouble>*  pbNMRDistanceAve;         // NMR distance restraint linear and
                                                   //   exponential averaged value
  GpuBuffer<PMEDouble2>* pbNMRDistanceTgtVal;      // NMR distance target and actual value for
                                                   //   current step
  GpuBuffer<int2>*       pbNMRDistanceStep;        // NMR distance first and last step for
                                                   //   application of restraint
  GpuBuffer<int>*        pbNMRDistanceInc;         // NMR distance increment for step weighting
  GpuBuffer<PMEDouble2>* pbNMRDistanceR1R2Slp;     // NMR distance r1, r2 slope
  GpuBuffer<PMEDouble2>* pbNMRDistanceR3R4Slp;     // NMR distance r3, r4 slope
  GpuBuffer<PMEDouble2>* pbNMRDistanceK2K3Slp;     // NMR distance k2, k3 slope
  GpuBuffer<PMEDouble>*  pbNMRDistanceK4Slp;       // NMR distance k4 slope
  GpuBuffer<PMEDouble2>* pbNMRDistanceR1R2Int;     // NMR distance r1, r2 intercept
  GpuBuffer<PMEDouble2>* pbNMRDistanceR3R4Int;     // NMR distance r3, r4 intercept
  GpuBuffer<PMEDouble2>* pbNMRDistanceK2K3Int;     // NMR distance k2, k3 intercept
  GpuBuffer<PMEDouble>*  pbNMRDistanceK4Int;       // NMR distance k4 intercept
  GpuBuffer<int2>*       pbNMRCOMDistanceID;       // NMR COM distance i, j
  GpuBuffer<int2>*       pbNMRCOMDistanceCOM;      // NMR COM distance COM ranges from .x to .y
  GpuBuffer<int2>*       pbNMRCOMDistanceCOMGrp;   // NMR COM distance COM indexing for no. of
                                                   //   atom ranges in a COM group from .x
                                                   //   to .y
  GpuBuffer<PMEDouble2>* pbNMRCOMDistanceR1R2;     // NMR COM distance computed r1, r2
  GpuBuffer<PMEDouble2>* pbNMRCOMDistanceR3R4;     // NMR COM distance computed r3, r4
  GpuBuffer<PMEDouble2>* pbNMRCOMDistanceK2K3;     // NMR COM distance computed k2, k3
  GpuBuffer<PMEDouble>*  pbNMRCOMDistanceK4;       // NMR COM distance computed k4
  GpuBuffer<PMEDouble>*  pbNMRCOMDistanceAve;      // NMR COM distance restraint linear and
                                                   //   exponential averaged value
  GpuBuffer<PMEDouble2>* pbNMRCOMDistanceTgtVal;   // NMR COM distance target and actual value
                                                   //   for current step
  GpuBuffer<int2>*       pbNMRCOMDistanceStep;     // NMR COM distance first and last step for
                                                   //   application of restraint
  GpuBuffer<int>*        pbNMRCOMDistanceInc;      // NMR COM distance increment for step
                                                   //   weighting
  GpuBuffer<PMEDouble2>* pbNMRCOMDistanceR1R2Slp;  // NMR COM distance r1, r2 slope
  GpuBuffer<PMEDouble2>* pbNMRCOMDistanceR3R4Slp;  // NMR COM distance r3, r4 slope
  GpuBuffer<PMEDouble2>* pbNMRCOMDistanceK2K3Slp;  // NMR COM distance k2, k3 slope
  GpuBuffer<PMEDouble>*  pbNMRCOMDistanceK4Slp;    // NMR COM distance k4 slope
  GpuBuffer<PMEDouble2>* pbNMRCOMDistanceR1R2Int;  // NMR COM distance r1, r2 intercept
  GpuBuffer<PMEDouble2>* pbNMRCOMDistanceR3R4Int;  // NMR COM distance r3, r4 intercept
  GpuBuffer<PMEDouble2>* pbNMRCOMDistanceK2K3Int;  // NMR COM distance k2, k3 intercept
  GpuBuffer<PMEDouble>*  pbNMRCOMDistanceK4Int;    // NMR COM distance k4 intercept
  GpuBuffer<int>*        pbNMRCOMDistanceWeights;  // NMR COM distance weights
  GpuBuffer<PMEDouble>*  pbNMRCOMDistanceXYZ;      // NMR COM X,Y,Z distance components
  GpuBuffer<int2>*       pbNMRr6avDistanceID;      // NMR r6av distance i, j
  GpuBuffer<int2>*       pbNMRr6avDistancer6av;    // NMR r6av distance r6av ranges from
                                                   //   .x to .y
  GpuBuffer<int2>*       pbNMRr6avDistancer6avGrp; // NMR r6av distance r6av indexing for no.
                                                   //   of atom ranges in a r6av group from
                                                   //   .x to .y
  GpuBuffer<PMEDouble2>* pbNMRr6avDistanceR1R2;    // NMR r6av distance computed r1, r2
  GpuBuffer<PMEDouble2>* pbNMRr6avDistanceR3R4;    // NMR r6av distance computed r3, r4
  GpuBuffer<PMEDouble2>* pbNMRr6avDistanceK2K3;    // NMR r6av distance computed k2, k3
  GpuBuffer<PMEDouble>*  pbNMRr6avDistanceK4;      // NMR r6av distance computed k4
  GpuBuffer<PMEDouble>*  pbNMRr6avDistanceAve;     // NMR r6av distance restraint linear and
                                                   //   exponential averaged value
  GpuBuffer<PMEDouble2>* pbNMRr6avDistanceTgtVal;  // NMR r6av distance target and actual
                                                   //   value for current step
  GpuBuffer<int2>*       pbNMRr6avDistanceStep;    // NMR r6av distance first and last step
                                                   //   for application of restraint
  GpuBuffer<int>*        pbNMRr6avDistanceInc;     // NMR r6av distance increment for step
                                                   //   weighting
  GpuBuffer<PMEDouble2>* pbNMRr6avDistanceR1R2Slp; // NMR r6av distance r1, r2 slope
  GpuBuffer<PMEDouble2>* pbNMRr6avDistanceR3R4Slp; // NMR r6av distance r3, r4 slope
  GpuBuffer<PMEDouble2>* pbNMRr6avDistanceK2K3Slp; // NMR r6av distance k2, k3 slope
  GpuBuffer<PMEDouble>*  pbNMRr6avDistanceK4Slp;   // NMR r6av distance k4 slope
  GpuBuffer<PMEDouble2>* pbNMRr6avDistanceR1R2Int; // NMR r6av distance r1, r2 intercept
  GpuBuffer<PMEDouble2>* pbNMRr6avDistanceR3R4Int; // NMR r6av distance r3, r4 intercept
  GpuBuffer<PMEDouble2>* pbNMRr6avDistanceK2K3Int; // NMR r6av distance k2, k3 intercept
  GpuBuffer<PMEDouble>*  pbNMRr6avDistanceK4Int;   // NMR r6av distance k4 intercept
  GpuBuffer<int>*        pbNMRr6avDistanceWeights;  // NMR r6av distance weights
  GpuBuffer<PMEDouble>*  pbNMRr6avDistanceXYZ;      // NMR r6av X,Y,Z distance components
 
  GpuBuffer<int2>*       pbNMRAngleID1;            // NMR angle i, j
  GpuBuffer<int>*        pbNMRAngleID2;            // NMR angle k
  GpuBuffer<PMEDouble2>* pbNMRAngleR1R2;           // NMR angle computed r1, r2
  GpuBuffer<PMEDouble2>* pbNMRAngleR3R4;           // NMR angle computed r3, r4
  GpuBuffer<PMEDouble2>* pbNMRAngleK2K3;           // NMR angle computed k2, k3
  GpuBuffer<PMEDouble>*  pbNMRAngleK4;             // NMR angle computed k4
  GpuBuffer<PMEDouble>*  pbNMRAngleAve;            // NMR angle restraint linear and
                                                   //   exponential averaged value
  GpuBuffer<PMEDouble2>* pbNMRAngleTgtVal;         // NMR angle target and actual value for
                                                   //   current step
  GpuBuffer<int2>*       pbNMRAngleStep;           // NMR angle first and last step for
                                                   //   application of restraint
  GpuBuffer<int>*        pbNMRAngleInc;            // NMR angle increment for step weighting
  GpuBuffer<PMEDouble2>* pbNMRAngleR1R2Slp;        // NMR angle r1, r2 slope
  GpuBuffer<PMEDouble2>* pbNMRAngleR3R4Slp;        // NMR angle r3, r4 slope
  GpuBuffer<PMEDouble2>* pbNMRAngleK2K3Slp;        // NMR angle k2, k3 slope
  GpuBuffer<PMEDouble>*  pbNMRAngleK4Slp;          // NMR angle k4 slope
  GpuBuffer<PMEDouble2>* pbNMRAngleR1R2Int;        // NMR angle r1, r2 intercept
  GpuBuffer<PMEDouble2>* pbNMRAngleR3R4Int;        // NMR angle r3, r4 intercept
  GpuBuffer<PMEDouble2>* pbNMRAngleK2K3Int;        // NMR angle k2, k3 intercept
  GpuBuffer<PMEDouble>*  pbNMRAngleK4Int;          // NMR angle k4 intercept
  //COM angle 
  GpuBuffer<int2>*       pbNMRCOMAngleID1;       // NMR COM angle i, j
  GpuBuffer<int>*        pbNMRCOMAngleID2;      // NMR angle k
  GpuBuffer<int2>*       pbNMRCOMAngleCOM;      // NMR COM angle COM ranges from .x to .y
  GpuBuffer<int2>*       pbNMRCOMAngleCOMGrp;   // NMR COM angle COM indexing for no. of
                                                   //   atom ranges in a COM group from .x
                                                   //   to .y
  GpuBuffer<PMEDouble2>* pbNMRCOMAngleR1R2;     // NMR COM angle computed r1, r2
  GpuBuffer<PMEDouble2>* pbNMRCOMAngleR3R4;     // NMR COM angle computed r3, r4
  GpuBuffer<PMEDouble2>* pbNMRCOMAngleK2K3;     // NMR COM angle computed k2, k3
  GpuBuffer<PMEDouble>*  pbNMRCOMAngleK4;       // NMR COM angle computed k4
  GpuBuffer<PMEDouble>*  pbNMRCOMAngleAve;      // NMR COM angle restraint linear and
                                                   //   exponential averaged value
  GpuBuffer<PMEDouble2>* pbNMRCOMAngleTgtVal;   // NMR COM angle target and actual value
                                                   //   for current step
  GpuBuffer<int2>*       pbNMRCOMAngleStep;     // NMR COM angle first and last step for
                                                   //   application of restraint
  GpuBuffer<int>*        pbNMRCOMAngleInc;      // NMR COM angle increment for step
                                                   //   weighting
  GpuBuffer<PMEDouble2>* pbNMRCOMAngleR1R2Slp;  // NMR COM angle r1, r2 slope
  GpuBuffer<PMEDouble2>* pbNMRCOMAngleR3R4Slp;  // NMR COM angle r3, r4 slope
  GpuBuffer<PMEDouble2>* pbNMRCOMAngleK2K3Slp;  // NMR COM angle k2, k3 slope
  GpuBuffer<PMEDouble>*  pbNMRCOMAngleK4Slp;    // NMR COM angle k4 slope
  GpuBuffer<PMEDouble2>* pbNMRCOMAngleR1R2Int;  // NMR COM angle r1, r2 intercept
  GpuBuffer<PMEDouble2>* pbNMRCOMAngleR3R4Int;  // NMR COM angle r3, r4 intercept
  GpuBuffer<PMEDouble2>* pbNMRCOMAngleK2K3Int;  // NMR COM angle k2, k3 intercept
  GpuBuffer<PMEDouble>*  pbNMRCOMAngleK4Int;    // NMR COM angle k4 intercept
 
  GpuBuffer<int4>*       pbNMRTorsionID1;          // NMR torsion i, j, k, l
  GpuBuffer<PMEDouble2>* pbNMRTorsionR1R2;         // NMR torsion computed r1, r2
  GpuBuffer<PMEDouble2>* pbNMRTorsionR3R4;         // NMR torsion computed r3, r4
  GpuBuffer<PMEDouble2>* pbNMRTorsionK2K3;         // NMR torsion computed k2, k3
  GpuBuffer<PMEDouble>*  pbNMRTorsionK4;           // NMR torsion computed k4
  GpuBuffer<PMEDouble>*  pbNMRTorsionAve1;         // NMR torsion restraint linear and
                                                   //   exponential averaged value
  GpuBuffer<PMEDouble>*  pbNMRTorsionAve2;         // NMR torsion restraint linear and
                                                   //   exponential averaged value
  GpuBuffer<PMEDouble2>* pbNMRTorsionTgtVal;       // NMR torsion target and actual value for
                                                   //   current step
  GpuBuffer<int2>*       pbNMRTorsionStep;         // NMR torsion first and last step for
                                                   //   application of restraint
  GpuBuffer<int>*        pbNMRTorsionInc;          // NMR torsion increment for step weighting
  GpuBuffer<PMEDouble2>* pbNMRTorsionR1R2Slp;      // NMR torsion r1, r2 slope
  GpuBuffer<PMEDouble2>* pbNMRTorsionR3R4Slp;      // NMR torsion r3, r4 slope
  GpuBuffer<PMEDouble2>* pbNMRTorsionK2K3Slp;      // NMR torsion k2, k3 slope
  GpuBuffer<PMEDouble>*  pbNMRTorsionK4Slp;        // NMR torsion k4 slope
  GpuBuffer<PMEDouble2>* pbNMRTorsionR1R2Int;      // NMR torsion r1, r2 intercept
  GpuBuffer<PMEDouble2>* pbNMRTorsionR3R4Int;      // NMR torsion r3, r4 intercept
  GpuBuffer<PMEDouble2>* pbNMRTorsionK2K3Int;      // NMR torsion k2, k3 intercept
  GpuBuffer<PMEDouble>*  pbNMRTorsionK4Int;        // NMR torsion k4 intercept
  // COM torsion 
  GpuBuffer<int4>*       pbNMRCOMTorsionID1;       // NMR COM torsion i, j, k, l
  GpuBuffer<int2>*       pbNMRCOMTorsionCOM;      // NMR COM torsion COM ranges from .x to .y
  GpuBuffer<int2>*       pbNMRCOMTorsionCOMGrp;   // NMR COM torsion COM indexing for no. of
                                                   //   atom ranges in a COM group from .x
                                                   //   to .y
  GpuBuffer<PMEDouble2>* pbNMRCOMTorsionR1R2;     // NMR COM torsion computed r1, r2
  GpuBuffer<PMEDouble2>* pbNMRCOMTorsionR3R4;     // NMR COM torsion computed r3, r4
  GpuBuffer<PMEDouble2>* pbNMRCOMTorsionK2K3;     // NMR COM torsion computed k2, k3
  GpuBuffer<PMEDouble>*  pbNMRCOMTorsionK4;       // NMR COM torsion computed k4
  GpuBuffer<PMEDouble>*  pbNMRCOMTorsionAve;      // NMR COM torsion restraint linear and
                                                   //   exponential averaged value
  GpuBuffer<PMEDouble2>* pbNMRCOMTorsionTgtVal;   // NMR COM torsion target and actual value
                                                   //   for current step
  GpuBuffer<int2>*       pbNMRCOMTorsionStep;     // NMR COM torsion first and last step for
                                                   //   application of restraint
  GpuBuffer<int>*        pbNMRCOMTorsionInc;      // NMR COM torsion increment for step
                                                   //   weighting
  GpuBuffer<PMEDouble2>* pbNMRCOMTorsionR1R2Slp;  // NMR COM torsion r1, r2 slope
  GpuBuffer<PMEDouble2>* pbNMRCOMTorsionR3R4Slp;  // NMR COM torsion r3, r4 slope
  GpuBuffer<PMEDouble2>* pbNMRCOMTorsionK2K3Slp;  // NMR COM torsion k2, k3 slope
  GpuBuffer<PMEDouble>*  pbNMRCOMTorsionK4Slp;    // NMR COM torsion k4 slope
  GpuBuffer<PMEDouble2>* pbNMRCOMTorsionR1R2Int;  // NMR COM torsion r1, r2 intercept
  GpuBuffer<PMEDouble2>* pbNMRCOMTorsionR3R4Int;  // NMR COM torsion r3, r4 intercept
  GpuBuffer<PMEDouble2>* pbNMRCOMTorsionK2K3Int;  // NMR COM torsion k2, k3 intercept
  GpuBuffer<PMEDouble>*  pbNMRCOMTorsionK4Int;    // NMR COM torsion k4 intercept
 

  // Force accumulators for atomic ops
  GpuBuffer<PMEAccumulator>* pbReffAccumulator;          // Effective Born Radius accumulator
  GpuBuffer<PMEAccumulator>* pbSumdeijdaAccumulator;     // Sumdeijda accumulator
  GpuBuffer<PMEAccumulator>* pbForceAccumulator;         // Force accumulators
  GpuBuffer<KineticEnergy>*  pbKineticEnergyBuffer;      // Kinetic energy accumulation buffer
  GpuBuffer<AFEKineticEnergy>* pbAFEKineticEnergyBuffer; // Extra terms for Alchemical Free
                                                         //   Energy kinetic energy
  GpuBuffer<unsigned long long int>* pbEnergyBuffer;     // Energy accumulation buffer
  GpuBuffer<unsigned long long int>* pbAFEBuffer;        // Alchemical free energy accumulation
                                                         //   buffer

  // MCWAT Buffers
  GpuBuffer<PMEFloat>* pbStericRadius;       // Radius for each atom to generate steric grid
  GpuBuffer<int>* pbVoxelOffset;       // Radius for each atom to generate steric grid
  GpuBuffer<unsigned int>* pbStericGrid;     // Steric grid
  GpuBuffer<int>* pbStericEmptyVoxels;       // Empty voxels in steric grid
  GpuBuffer<int>* pbStericThreadCount;       // Counts for each thread in steric grid
  GpuBuffer<int>* pbCoarseGrid;              // Coarse Grid
  GpuBuffer<int>* pbRandomWater;             // Random water atom ids for coarse grid
  GpuBuffer<PMEFloat>* pbCoarseTrans;        // Water translation for coarse grid
  GpuBuffer<PMEFloat>* pbCoarseGridEne;      // Coarse grid energy for a water move

  // AMD buffers
  PMEDouble* pAmdWeightsAndEnergy;           // AMD
  GpuBuffer<PMEDouble>* pbAMDfwgtd;

  // GaMD buffers
  PMEDouble* pGaMDWeightsAndEnergy;          // GaMD
  GpuBuffer<PMEDouble>* pbGaMDfwgtd;

  // Constant pH stuff
  GpuBuffer<double>* pbChargeRefreshBuffer;  // Pointer to new charges buffer

  // Extra points data
  GpuBuffer<int4>*   pbExtraPoint11Frame;        // Type 1 single extra point tuple:
                                                 //   parent_atm, atm1, atm2, atm3
  GpuBuffer<int>*    pbExtraPoint11Index;        // Type 1 single extra point EP
  GpuBuffer<double>* pbExtraPoint11;             // Type 1 single extra point local coordinates
  GpuBuffer<int4>*   pbExtraPoint12Frame;        // Type 2 single extra point tuple:
                                                 //   parent_atm, atm1, atm2, atm3
  GpuBuffer<int>*    pbExtraPoint12Index;        // Type 2 single extra point EP
  GpuBuffer<double>* pbExtraPoint12;             // Type 2 single extra point local coordinates
  GpuBuffer<int4>*   pbExtraPoint21Frame;        // Type 1 dual extra point tuple:
                                                 //   parent_atm, atm1, atm2, atm3
  GpuBuffer<int2>*   pbExtraPoint21Index;        // Type 1 dual extra point EP1, EP2
  GpuBuffer<double>* pbExtraPoint21;             // Type 1 dual extra point local coordinates
  GpuBuffer<int4>*   pbExtraPoint22Frame;        // Type 2 dual extra point tuple:
                                                 //   parent_atm, atm1, atm2, atm3
  GpuBuffer<int2>*   pbExtraPoint22Index;        // Type 2 dual extra point EP1. EP2
  GpuBuffer<double>* pbExtraPoint22;             // Type 2 dual extra point local coordinates
  GpuBuffer<int4>*   pbExtraPointCustomFrame;    // Remapped custom EP frames
  GpuBuffer<int>*    pbExtraPointCustomIndex;    // Remapped custom extra points
  GpuBuffer<double>* pbExtraPointCustom;         // Custom extra point frame measurements (this
                                                 //   is equivalent to the local coordinates of
                                                 //   other classes of EPs)
  GpuBuffer<int4>*   pbImageExtraPoint11Frame;   // Remapped Type 1 single extra point tuple:
                                                 //   parent_atm, atm1, atm2, atm3
  GpuBuffer<int>*    pbImageExtraPoint11Index;   // Remapped Type 1 single extra point EP
  GpuBuffer<int4>*   pbImageExtraPoint12Frame;   // Remapped Type 2 single extra point tuple:
                                                 //   parent_atm, atm1, atm2, atm3
  GpuBuffer<int>*    pbImageExtraPoint12Index;   // Remapped Type 2 single extra point EP
  GpuBuffer<int4>*   pbImageExtraPoint21Frame;   // Remapped Type 1 dual extra point tuple:
                                                 //   parent_atm, atm1, atm2, atm3
  GpuBuffer<int2>*   pbImageExtraPoint21Index;   // Remapped type 1 dual extra point EP1, EP2
  GpuBuffer<int4>*   pbImageExtraPoint22Frame;   // Remapped Type 2 dual extra point tuple:
                                                 //   parent_atm, atm1, atm2, atm3
  GpuBuffer<int2>*   pbImageExtraPoint22Index;   // Remapped Type 2 dual extra point EP1. EP2
  GpuBuffer<int>*    pbImageEPCustomIndex;       // Remapped custom extra points, neighbor list
  GpuBuffer<int4>*   pbImageEPCustomFrame;       // Remapped custom EP frames, neighbor list
  
  // SHAKE constraints
  bool bUseHMR;                             // Use HMR for SHAKE flag
  GpuBuffer<int4>*    pbShakeID;            // SHAKE central atom plus up to 3 hydrogens
  GpuBuffer<double2>* pbShakeParm;          // SHAKE central atom mass and equilibrium bond
                                            //   length
  GpuBuffer<double>*  pbShakeInvMassH;      // SHAKE HMR inverse hydrogen mass
  GpuBuffer<int4>*    pbFastShakeID;        // H2O oxygen plus two hydrogens atom ID
  GpuBuffer<int>*     pbSlowShakeID1;       // Central atom of XH4 Shake constraint
  GpuBuffer<int4>*    pbSlowShakeID2;       // XH4 SHAKE constraint hydrogens
  GpuBuffer<double2>* pbSlowShakeParm;      // XH4 SHAKE central atom mass and equilibrium
                                            //   bond length
  GpuBuffer<double>*  pbSlowShakeInvMassH;  // XH4 SHAKE HMR inverse hydrogen mass

  // Launch parameters: these set the numbers of threads for each block relating to
  // GPU multithreading.  Allow -1 for gpu_device_id default --> choose gpu based on memory.
  int gpu_device_id;
  unsigned int blocks;
  unsigned int threadsPerBlock;
  unsigned int clearForcesThreadsPerBlock;
  unsigned int NLClearForcesThreadsPerBlock;
  unsigned int reduceForcesThreadsPerBlock;
  unsigned int NLReduceForcesThreadsPerBlock;
  unsigned int reduceBufferThreadsPerBlock;
  unsigned int localForcesBlocks;
  unsigned int localForcesThreadsPerBlock;
  unsigned int AFEExchangeBlocks;
  unsigned int AFEExchangeThreadsPerBlock;
  unsigned int CHARMMForcesBlocks;
  unsigned int CHARMMForcesThreadsPerBlock;
  unsigned int NMRForcesBlocks;
  unsigned int NMRForcesThreadsPerBlock;
  unsigned int GBBornRadiiThreadsPerBlock;
  unsigned int GBBornRadiiIGB78ThreadsPerBlock;
  unsigned int GBNonbondEnergy1ThreadsPerBlock;
  unsigned int GBNonbondEnergy2ThreadsPerBlock;
  unsigned int GBNonbondEnergy2IGB78ThreadsPerBlock;
  unsigned int PMENonbondEnergyThreadsPerBlock;
  unsigned int PMENonbondForcesThreadsPerBlock;
  unsigned int IPSNonbondEnergyThreadsPerBlock;
  unsigned int IPSNonbondForcesThreadsPerBlock;
  unsigned int updateThreadsPerBlock;
  unsigned int shakeThreadsPerBlock;
  unsigned int generalThreadsPerBlock;
  unsigned int GBBornRadiiBlocks;
  unsigned int GBNonbondEnergy1Blocks;
  unsigned int GBNonbondEnergy2Blocks;
  unsigned int PMENonbondBlocks;
  unsigned int IPSNonbondBlocks;
  unsigned int BNLBlocks;
  unsigned int NLCalculateOffsetsThreadsPerBlock;
  unsigned int NLBuildNeighborList32ThreadsPerBlock;
  unsigned int NLBuildNeighborList16ThreadsPerBlock;
  unsigned int NLBuildNeighborList8ThreadsPerBlock;
  unsigned int updateBlocks;
  unsigned int bondWorkBlocks;
  int          ExpansionBlocks;
  unsigned int readSize;

  // Kludge for HW bug on Ge Force GPUs.  The Ge Force GPUs
  // freak out when a kernel uses lots of DP and Texture
  bool bNoDPTexture;

  // Nonbond Kernel stuff
  GpuBuffer<unsigned int>* pbExclusion;   // Exclusion masks
  GpuBuffer<unsigned int>* pbWorkUnit;    // Work unit list
  GpuBuffer<unsigned int>* pbGBPosition;  // Nonbond kernel worunit positions

  // All simulation data is destined for GPU constant RAM
  // Single-node multi-gpu parameters
  bool bSingleNode;         // Flag to indicate MPI run is all on one node
  bool bP2P;                // Flag to indicate P2P connectivity between all processes
  bool bNCCL;               // Flag to indicate NCCL is used

  // Multi-gpu parameters
  int nGpus;                // Number of GPUs involved in calculation
  int gpuID;                // Local GPU ID
  int minLocalCell;         // First nonbond cell owned by local GPU
  int maxLocalCell;         // End of nonbond cells owned by local GPU
  int* pMinLocalCell;       // First cell owned by each node
  int* pMaxLocalCell;       // Last cell owned by each node
  int* pMinLocalAtom;       // First atom owned by each node
  int* pMaxLocalAtom;       // Last atom owned by each node
#ifdef MPI
  MPI_Comm comm;            // MPI Communicator for all collective MPI calls
#ifdef NCCL
  ncclComm_t ncclComm;      // NCCL Communicator  
#endif
#endif
#ifdef AWSMM
  // AWSMM stuff
  int postProcessingFlags;  // List of desired post-processing steps
  int nAlphaCarbons;        // Number of alpha carbons in system
  int* pAlphaCarbonIndex;   // Alpha carbon indices
  double* pRefAlphaCarbon;  // Alpha carbon initial coordinates
  FILE* pPPRMSD;            // Pointer to RMSD output file
  FILE* pPPEnergy;          // Pointer to energy output file
  FILE* pPPVelocity;        // Pointer to velocity output file
#endif
};
#endif //_BASE_GPU_CONTEXT
