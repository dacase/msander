#include "base_gpuContext.h"
#include "bondRemap.h"

void base_gpuContext::Init()
{
  sim.InitSimulationConst();
}

base_gpuContext::base_gpuContext() :

  bECCSupport(false), bCanMapHostMemory(false), totalCPUMemory(0), totalGPUMemory(0),
  totalMemory(0), ntt(0), ntpr(0), ntwe(0),
#ifdef MPI
  nGpus(1), gpuID(0),
#endif
  sm_version(SM_3X), bNoDPTexture(false), imin(0), ntb(1), ips(0), ntc(1), ntf(1), ntr(0),
  bCPURandoms(false), gbsa(0), step(0), forwardPlan(0), backwardPlan(0), bPencilYShift(false),
  nonbond_skin(2.0), sortTempBytes(0),

  // Atom data, in single and double precision
  pbAtom(NULL), pbAtomXYSP(NULL), pbAtomZSP(NULL), pbAtomSigEps(NULL), pbAtomLJID(NULL),
  pbAtomCharge(NULL), pbAtomChargeSP(NULL), pbAtomChargeSPLJID(NULL), pbAtomRBorn(NULL),
  pbAtomS(NULL), pbAtomMass(NULL), pbCenter(NULL), pbR4sources(NULL),
#ifdef MPI
  pMinLocalCell(NULL), pMaxLocalCell(NULL), pMinLocalAtom(NULL), pMaxLocalAtom(NULL),
  pbPeerAccumulator(NULL), pPeerAccumulatorList(NULL), pPeerAccumulatorMemHandle(NULL),
  bCalculateLocalForces(false), bCalculateDirectSum(false), bCalculateReciprocalSum(false),
  bSingleNode(false), bP2P(false),
#endif
  pbVel(NULL), pbLVel(NULL), bLocalInteractions(true), bCharmmInteractions(false),
  bNMRInteractions(false), NMRnstep(0), pbBond(NULL), pbBondID(NULL), pbBondAngle(NULL),
  pbAFEMolType(NULL), pbAFEMolPartner(NULL), pbBondAngleID1(NULL), pbBondAngleID2(NULL),
  pbDihedral1(NULL), pbDihedral2(NULL), pbDihedral3(NULL), pbDihedralID1(NULL), pbNb141(NULL),
  pbNb142(NULL), pbNb14ID(NULL), pbConstraint1(NULL), pbConstraint2(NULL),
  pbConstraintID(NULL), pbUBAngle(NULL), pbUBAngleID(NULL), pbImpDihedral(NULL),
  pbImpDihedralID1(NULL), pbCmapID1(NULL), pbCmapID2(NULL), pbCmapType(NULL),
  pbCmapEnergy(NULL),

  // NMR-related data
  pbNMRJarData(NULL), pbNMRDistanceID(NULL), pbNMRDistanceR1R2(NULL), pbNMRDistanceR3R4(NULL),
  pbNMRDistanceK2K3(NULL), pbNMRDistanceK4(NULL), pbNMRDistanceAve(NULL),
  pbNMRDistanceTgtVal(NULL), pbNMRDistanceStep(NULL), pbNMRDistanceInc(NULL),
  pbNMRDistanceR1R2Slp(NULL), pbNMRDistanceR3R4Slp(NULL), pbNMRDistanceK2K3Slp(NULL),
  pbNMRDistanceK4Slp(NULL), pbNMRDistanceR1R2Int(NULL), pbNMRDistanceR3R4Int(NULL),
  pbNMRDistanceK2K3Int(NULL), pbNMRDistanceK4Int(NULL), pbNMRCOMDistanceID(NULL),
  pbNMRCOMDistanceCOM(NULL), pbNMRCOMDistanceCOMGrp(NULL), pbNMRCOMDistanceR1R2(NULL),
  pbNMRCOMDistanceR3R4(NULL), pbNMRCOMDistanceK2K3(NULL), pbNMRCOMDistanceK4(NULL),
  pbNMRCOMDistanceAve(NULL), pbNMRCOMDistanceTgtVal(NULL), pbNMRCOMDistanceStep(NULL),
  pbNMRCOMDistanceInc(NULL), pbNMRCOMDistanceR1R2Slp(NULL), pbNMRCOMDistanceR3R4Slp(NULL),
  pbNMRCOMDistanceK2K3Slp(NULL), pbNMRCOMDistanceK4Slp(NULL), pbNMRCOMDistanceR1R2Int(NULL),
  pbNMRCOMDistanceR3R4Int(NULL), pbNMRCOMDistanceK2K3Int(NULL), pbNMRCOMDistanceK4Int(NULL),
  pbNMRCOMDistanceWeights(NULL), pbNMRCOMDistanceXYZ(NULL), pbNMRr6avDistanceID(NULL),
  pbNMRr6avDistancer6av(NULL), pbNMRr6avDistancer6avGrp(NULL), pbNMRr6avDistanceR1R2(NULL),
  pbNMRr6avDistanceR3R4(NULL), pbNMRr6avDistanceK2K3(NULL), pbNMRr6avDistanceK4(NULL),
  pbNMRr6avDistanceAve(NULL), pbNMRr6avDistanceTgtVal(NULL), pbNMRr6avDistanceStep(NULL),
  pbNMRr6avDistanceInc(NULL), pbNMRr6avDistanceR1R2Slp(NULL), pbNMRr6avDistanceR3R4Slp(NULL),
  pbNMRr6avDistanceK2K3Slp(NULL), pbNMRr6avDistanceK4Slp(NULL), pbNMRr6avDistanceR1R2Int(NULL),
  pbNMRr6avDistanceR3R4Int(NULL), pbNMRr6avDistanceK2K3Int(NULL), pbNMRr6avDistanceK4Int(NULL),
  pbNMRr6avDistanceWeights(NULL), pbNMRr6avDistanceXYZ(NULL),

  pbNMRCOMAngleID1(NULL), pbNMRCOMAngleID2(NULL), pbNMRCOMAngleCOM(NULL), 
  pbNMRCOMAngleCOMGrp(NULL), pbNMRCOMAngleR1R2(NULL), pbNMRCOMAngleR3R4(NULL), 
  pbNMRCOMAngleK2K3(NULL), pbNMRCOMAngleK4(NULL), pbNMRCOMAngleAve(NULL), 
  pbNMRCOMAngleTgtVal(NULL), pbNMRCOMAngleStep(NULL), pbNMRCOMAngleInc(NULL), 
  pbNMRCOMAngleR1R2Slp(NULL), pbNMRCOMAngleR3R4Slp(NULL), pbNMRCOMAngleK2K3Slp(NULL),
  pbNMRCOMAngleK4Slp(NULL), pbNMRCOMAngleR1R2Int(NULL), pbNMRCOMAngleR3R4Int(NULL), 
  pbNMRCOMAngleK2K3Int(NULL), pbNMRCOMAngleK4Int(NULL),

  pbNMRAngleID1(NULL), pbNMRAngleID2(NULL), pbNMRAngleR1R2(NULL), pbNMRAngleR3R4(NULL),
  pbNMRAngleK2K3(NULL), pbNMRAngleK4(NULL), pbNMRAngleAve(NULL), pbNMRAngleTgtVal(NULL),
  pbNMRAngleStep(NULL), pbNMRAngleInc(NULL), pbNMRAngleR1R2Slp(NULL), pbNMRAngleR3R4Slp(NULL),
  pbNMRAngleK2K3Slp(NULL), pbNMRAngleK4Slp(NULL), pbNMRAngleR1R2Int(NULL),
  pbNMRAngleR3R4Int(NULL), pbNMRAngleK2K3Int(NULL), pbNMRAngleK4Int(NULL),
  pbNMRTorsionID1(NULL),
  pbNMRTorsionR1R2(NULL), pbNMRTorsionR3R4(NULL), pbNMRTorsionK2K3(NULL), pbNMRTorsionK4(NULL),
  pbNMRTorsionAve1(NULL), pbNMRTorsionAve2(NULL), pbNMRTorsionTgtVal(NULL),
  pbNMRTorsionStep(NULL), pbNMRTorsionInc(NULL), pbNMRTorsionR1R2Slp(NULL),
  pbNMRTorsionR3R4Slp(NULL), pbNMRTorsionK2K3Slp(NULL), pbNMRTorsionK4Slp(NULL),
  pbNMRTorsionR1R2Int(NULL), pbNMRTorsionR3R4Int(NULL), pbNMRTorsionK2K3Int(NULL),
  pbNMRTorsionK4Int(NULL),

  pbNMRCOMTorsionID1(NULL),
  pbNMRCOMTorsionCOM(NULL), pbNMRCOMTorsionCOMGrp(NULL), pbNMRCOMTorsionR1R2(NULL),
  pbNMRCOMTorsionR3R4(NULL), pbNMRCOMTorsionK2K3(NULL), pbNMRCOMTorsionK4(NULL),
  pbNMRCOMTorsionAve(NULL), pbNMRCOMTorsionTgtVal(NULL), pbNMRCOMTorsionStep(NULL),
  pbNMRCOMTorsionInc(NULL), pbNMRCOMTorsionR1R2Slp(NULL), pbNMRCOMTorsionR3R4Slp(NULL),
  pbNMRCOMTorsionK2K3Slp(NULL), pbNMRCOMTorsionK4Slp(NULL), pbNMRCOMTorsionR1R2Int(NULL),
  pbNMRCOMTorsionR3R4Int(NULL), pbNMRCOMTorsionK2K3Int(NULL), pbNMRCOMTorsionK4Int(NULL),
 
  // Other parameters
  pbReff(NULL), pbReffSP(NULL), pbPsi(NULL), pbTemp7(NULL), pbTIRegion(NULL),
  pbTILinearAtmID(NULL), pbBarLambda(NULL), pbBarTot(NULL),
  pbReffAccumulator(NULL), pbSumdeijdaAccumulator(NULL),
  pbForceAccumulator(NULL), pbEnergyBuffer(NULL), pbKineticEnergyBuffer(NULL),
  pbAFEBuffer(NULL), pbAFEKineticEnergyBuffer(NULL), pbWorkUnit(NULL),
  pbExclusion(NULL), pbGBPosition(NULL), pbNeckMaxValPos(NULL),
  pbGBAlphaBetaGamma(NULL), randomCounter(0), pbRandom(NULL), bUseHMR(false),
  pbShakeID(NULL), pbShakeParm(NULL), pbShakeInvMassH(NULL), pbFastShakeID(NULL),
  pbSlowShakeID1(NULL), pbSlowShakeID2(NULL), pbSlowShakeParm(NULL),
  pbSlowShakeInvMassH(NULL), pbUpdateIndex(NULL),

  // Extra points data
  pbExtraPoint11Frame(NULL), pbExtraPoint11Index(NULL),
  pbExtraPoint11(NULL), pbExtraPoint12Frame(NULL), pbExtraPoint12Index(NULL),
  pbExtraPoint12(NULL), pbExtraPoint21Frame(NULL), pbExtraPoint21Index(NULL),
  pbExtraPoint21(NULL), pbExtraPoint22Frame(NULL), pbExtraPoint22Index(NULL),
  pbExtraPoint22(NULL), pbImageExtraPoint11Frame(NULL), pbImageExtraPoint11Index(NULL),
  pbImageExtraPoint12Frame(NULL), pbImageExtraPoint12Index(NULL),
  pbImageExtraPoint21Frame(NULL), pbImageExtraPoint21Index(NULL),
  pbImageExtraPoint22Frame(NULL), pbImageExtraPoint22Index(NULL),

  // Other parameters
  pbChargeRefreshBuffer(NULL), pbPrefac(NULL), pbFract(NULL), pbTileBoundary(NULL),
#ifndef use_DPFP
  ErfcCoeffsTable(NULL), pbErfcCoeffsTable(NULL),
#endif
  pblliXYZ_q(NULL), pbXYZ_q(NULL), pbXYZ_qt(NULL),
  bNeighborList(false), bNeedNewNeighborList(true), bNewNeighborList(false), bSmallBox(false),
  bOddNLCells(false), neighborListBits(32), pbAtomXYSaveSP(NULL), pbAtomZSaveSP(NULL),
  pbImage(NULL), pbImageIndex(NULL), pbSortTemp(NULL), pbSubImageLookup(NULL), pbImageVel(NULL),
  pbImageLVel(NULL), pbImageMass(NULL), pbImageCharge(NULL), pbImageSigEps(NULL),
  pbImageLJID(NULL), pbImageCellID(NULL), pbImageTIRegion(NULL), pbImageTILinearAtmID(NULL),
  pbImageNMRCOMDistanceID(NULL), pbImageNMRCOMDistanceCOM(NULL),
  pbImageNMRCOMDistanceCOMGrp(NULL), pbImageNMRr6avDistanceID(NULL),
  pbImageNMRr6avDistancer6av(NULL), pbImageNMRr6avDistancer6avGrp(NULL),
  pbImageNMRCOMAngleID1(NULL), pbImageNMRCOMAngleID2(NULL), pbImageNMRCOMAngleCOM(NULL),
  pbImageNMRCOMAngleCOMGrp(NULL),
  pbImageNMRCOMTorsionID1(NULL), pbImageNMRCOMTorsionCOM(NULL), pbImageNMRCOMTorsionCOMGrp(NULL),
  pbImageShakeID(NULL), pbImageFastShakeID(NULL), pbImageSlowShakeID1(NULL),
  pbImageSlowShakeID2(NULL), pbImageSolventAtomID(NULL), pbImageSoluteAtomID(NULL),
  pbBNLExclusionBuffer(NULL), pbNLExclusionList(NULL), pbNLExclusionStartCount(NULL),
  pbNLAtomList(NULL), pbNLTotalOffset(NULL), pbFrcBlkCounters(NULL), pbNLEntries(NULL),
  pbNLNonbondCellStartEnd(NULL), pbNLbSkinTestFail(NULL), pbNLCellHash(NULL), pbNLRecord(NULL),
  pbNLEntry(NULL), maxSoluteMolecules(0), maxPSSoluteMolecules(0), pbSoluteAtomID(NULL),
  pbSoluteAtomMass(NULL), pbSolute(NULL), pbUllSolute(NULL), pbSolventAtomID(NULL),
  pbSolvent(NULL), pbNTPData(NULL), pbConstraintAtomX(NULL), pbConstraintAtomY(NULL),
  pbConstraintAtomZ(NULL), pbConstraintCOMX(NULL), pbConstraintCOMY(NULL),
  pbConstraintCOMZ(NULL), ee_plasma(0.0), self_energy(0.0), vdw_recip(0.0), pbLJTerm(NULL),
#ifdef AWSMM
  postProcessingFlags(0), nAlphaCarbons(0), pAlphaCarbonIndex(NULL), pRefAlphaCarbon(NULL),
  pPPRMSD(NULL), pPPEnergy(NULL), pPPVelocity(NULL),
#endif

  // Bonded term work units
  pbBondWorkUnitUINT(NULL), pbBondWorkUnitDBL2(NULL), pbBondWorkUnitPFLOAT(NULL),
  pbBondWorkUnitPFLOAT2(NULL),

  // Atom shuttling
  pbDataShuttle(NULL), pbShuttleTickets(NULL),

  //NEB buffers
  pbRMSMask(NULL), pbFitMask(NULL), pbAtmIdx(NULL), pbNEBEnergyAll(NULL),
  pbTangents(NULL), pbSpringForce(NULL), pbNEBForce(NULL),
  pbNextDataShuttle(NULL), pbPrevDataShuttle(NULL), pbKabschCOM(NULL),
  pbtotFitMass(NULL), pbDataSPR(NULL), pbRotAtm(NULL),

  //MCWat buffers
  pbStericRadius(NULL), pbStericGrid(NULL), pbStericThreadCount(NULL),
  pbCoarseGrid(NULL), pbRandomWater(NULL), pbCoarseTrans(NULL),
  pbCoarseGridEne(NULL), pbStericEmptyVoxels(NULL), pbVoxelOffset(NULL),

  // AMD buffers
  pbAMDfwgtd(NULL), pAmdWeightsAndEnergy(NULL),

  // GaMD buffers
  pbGaMDfwgtd(NULL), pGaMDWeightsAndEnergy(NULL)
{
  sim.InitSimulationConst();
}

//---------------------------------------------------------------------------------------------
// Destructor for the _gpuContext class
//---------------------------------------------------------------------------------------------
base_gpuContext::~base_gpuContext()
{
  int i;
  // Delete Atom data
  delete pbAtom;
  delete pbAtomXYSP;
  delete pbAtomZSP;
  delete pbAtomSigEps;
  delete pbAtomLJID;
  delete pbAtomRBorn;
  delete pbAtomS;
  delete pbAtomCharge;
  delete pbAtomChargeSP;
  delete pbAtomChargeSPLJID;
  delete pbAtomMass;
  delete pbReff;
  delete pbReffSP;
  delete pbPsi;
  delete pbTemp7;
  delete pbTIRegion;
  delete pbTILinearAtmID;
  delete pbBarLambda;
  delete pbBarTot;
#ifdef MPI
  delete[] pMinLocalCell;
  delete[] pMaxLocalCell;
  delete[] pMinLocalAtom;
  delete[] pMaxLocalAtom;
  delete pbPeerAccumulator;
  delete[] pPeerAccumulatorList;
  delete[] pPeerAccumulatorMemHandle;
#endif
  delete pbVel;
  delete pbLVel;
  delete pbCenter;
  delete pbLJTerm;
  if (sim.nR4sources > 0) {
    delete pbR4sources;
  }

  // Delete PME stuff
  delete pbPrefac;
  delete pbFract;
  delete pblliXYZ_q;
  delete pbXYZ_q;
  delete pbXYZ_qt;
#ifndef use_DPFP
  delete ErfcCoeffsTable;
  delete pbErfcCoeffsTable;
#endif
  cufftDestroy(forwardPlan);
  cufftDestroy(backwardPlan);

  // Delete neighbor list stuff
  delete pbAtomXYSaveSP;
  delete pbAtomZSaveSP;
  delete pbImage;
  delete pbSortTemp;
  delete pbImageIndex;
  delete pbImageVel;
  delete pbImageLVel;
  delete pbImageMass;
  delete pbImageCharge;
  delete pbImageSigEps;
  delete pbImageLJID;
  delete pbImageCellID;
  delete pbImageTIRegion;
  if (sim.ti_mode != 0 && sim.TIPaddedLinearAtmCnt > 0) {
    delete pbUpdateIndex;
    delete pbImageTILinearAtmID;
  }
  delete pbImageShakeID;
  delete pbImageFastShakeID;
  delete pbImageSlowShakeID1;
  delete pbImageSlowShakeID2;
  delete pbImageSolventAtomID;
  delete pbImageSoluteAtomID;
  delete pbNLNonbondCellStartEnd;
  delete pbBNLExclusionBuffer;
  delete pbNLExclusionList;
  delete pbNLExclusionStartCount;
  delete pbNLAtomList;
  delete pbNLTotalOffset;
  delete pbFrcBlkCounters;
  delete pbNLEntries;
  delete pbNLbSkinTestFail;
  delete pbNLCellHash;
  delete pbNLRecord;
  delete pbNLEntry;
  
  // Delete bond work units data
  delete pbBondWorkUnitUINT;
  delete pbBondWorkUnitDBL2;
  delete pbBondWorkUnitPFLOAT;
  delete pbBondWorkUnitPFLOAT2;
  for (i = 0; i < sim.bondWorkUnits; i++) {
    DestroyBondWorkUnit(&bondWorkUnitRecord[i]);
  }
  free(bondWorkUnitRecord);

  // Delete atom shuttling stuff
  delete pbDataShuttle;
  delete pbShuttleTickets;

  // Delete MC Wat
  delete pbStericRadius;
  delete pbVoxelOffset;
  delete pbStericEmptyVoxels;
  delete pbStericGrid;
  delete pbStericThreadCount;
  delete pbCoarseGrid;
  delete pbRandomWater;
  delete pbCoarseTrans;
  delete pbCoarseGridEne;

  // Delete NEB stuff
  delete pbRMSMask;
  delete pbFitMask;
  delete pbAtmIdx;
  delete pbNEBEnergyAll;
  delete pbTangents;
  delete pbSpringForce;
  delete pbNEBForce;
  delete pbNextDataShuttle;
  delete pbPrevDataShuttle;
  delete pbKabschCOM;
  delete pbtotFitMass;
  delete pbDataSPR;
  delete pbRotAtm;

  // Delete NTP stuff
  delete pbSoluteAtomID;
  delete pbSoluteAtomMass;
  delete pbSolute;
  delete pbUllSolute;
  delete pbSolventAtomID;
  delete pbSolvent;
  delete pbNTPData;
  delete pbConstraintAtomX;
  delete pbConstraintAtomY;
  delete pbConstraintAtomZ;
  delete pbConstraintCOMX;
  delete pbConstraintCOMY;
  delete pbConstraintCOMZ;

  // Delete GB stuff
  delete pbWorkUnit;
  delete pbExclusion;
  delete pbGBPosition;
  delete pbNeckMaxValPos;
  delete pbGBAlphaBetaGamma;

  // Delete random number stuff
  delete pbRandom;

  // Delete bonded parameter data
  delete pbBond;
  delete pbBondID;
  delete pbBondAngle;
  delete pbBondAngleID1;
  delete pbBondAngleID2;
  delete pbDihedral1;
  delete pbDihedral2;
  delete pbDihedral3;
  delete pbDihedralID1;
  delete pbNb141;
  delete pbNb142;
  delete pbNb14ID;
  delete pbConstraint1;
  delete pbConstraint2;
  delete pbConstraintID;
  delete pbUBAngle;
  delete pbUBAngleID;
  delete pbImpDihedral;
  delete pbImpDihedralID1;
  delete pbCmapID1;
  delete pbCmapID2;
  delete pbCmapType;
  delete pbCmapEnergy;

  // Delete NMR stuff
  delete pbNMRJarData;
  delete pbNMRDistanceID;
  delete pbNMRDistanceR1R2;
  delete pbNMRDistanceR3R4;
  delete pbNMRDistanceK2K3;
  delete pbNMRDistanceAve;
  delete pbNMRDistanceTgtVal;
  delete pbNMRDistanceStep;
  delete pbNMRDistanceInc;
  delete pbNMRDistanceR1R2Slp;
  delete pbNMRDistanceR3R4Slp;
  delete pbNMRDistanceK2K3Slp;
  delete pbNMRDistanceR1R2Int;
  delete pbNMRDistanceR3R4Int;
  delete pbNMRDistanceK2K3Int;
  delete pbNMRCOMDistanceID;
  delete pbNMRCOMDistanceCOM;
  delete pbNMRCOMDistanceCOMGrp;
  delete pbNMRCOMDistanceR1R2;
  delete pbNMRCOMDistanceR3R4;
  delete pbNMRCOMDistanceK2K3;
  delete pbNMRCOMDistanceAve;
  delete pbNMRCOMDistanceTgtVal;
  delete pbNMRCOMDistanceStep;
  delete pbNMRCOMDistanceInc;
  delete pbNMRCOMDistanceR1R2Slp;
  delete pbNMRCOMDistanceR3R4Slp;
  delete pbNMRCOMDistanceK2K3Slp;
  delete pbNMRCOMDistanceR1R2Int;
  delete pbNMRCOMDistanceR3R4Int;
  delete pbNMRCOMDistanceK2K3Int;
  delete pbNMRCOMDistanceWeights;
  delete pbNMRCOMDistanceXYZ;
  delete pbNMRr6avDistanceID;
  delete pbNMRr6avDistancer6av;
  delete pbNMRr6avDistancer6avGrp;
  delete pbNMRr6avDistanceR1R2;
  delete pbNMRr6avDistanceR3R4;
  delete pbNMRr6avDistanceK2K3;
  delete pbNMRr6avDistanceAve;
  delete pbNMRr6avDistanceTgtVal;
  delete pbNMRr6avDistanceStep;
  delete pbNMRr6avDistanceInc;
  delete pbNMRr6avDistanceR1R2Slp;
  delete pbNMRr6avDistanceR3R4Slp;
  delete pbNMRr6avDistanceK2K3Slp;
  delete pbNMRr6avDistanceR1R2Int;
  delete pbNMRr6avDistanceR3R4Int;
  delete pbNMRr6avDistanceK2K3Int;
  delete pbNMRAngleID1;
  delete pbNMRAngleID2;
  delete pbNMRAngleR1R2;
  delete pbNMRAngleR3R4;
  delete pbNMRAngleK2K3;
  delete pbNMRAngleAve;
  delete pbNMRAngleTgtVal;
  delete pbNMRAngleStep;
  delete pbNMRAngleInc;
  delete pbNMRAngleR1R2Slp;
  delete pbNMRAngleR3R4Slp;
  delete pbNMRAngleK2K3Slp;
  delete pbNMRAngleR1R2Int;
  delete pbNMRAngleR3R4Int;
  delete pbNMRAngleK2K3Int;

  delete pbNMRCOMAngleID1;
  delete pbNMRCOMAngleID2;
  delete pbNMRCOMAngleCOM;
  delete pbNMRCOMAngleCOMGrp;
  delete pbNMRCOMAngleR1R2;
  delete pbNMRCOMAngleR3R4;
  delete pbNMRCOMAngleK2K3;
  delete pbNMRCOMAngleAve;
  delete pbNMRCOMAngleTgtVal;
  delete pbNMRCOMAngleStep;
  delete pbNMRCOMAngleInc;
  delete pbNMRCOMAngleR1R2Slp;
  delete pbNMRCOMAngleR3R4Slp;
  delete pbNMRCOMAngleK2K3Slp;
  delete pbNMRCOMAngleR1R2Int;
  delete pbNMRCOMAngleR3R4Int;
  delete pbNMRCOMAngleK2K3Int;
 
  delete pbNMRTorsionID1;
  delete pbNMRTorsionR1R2;
  delete pbNMRTorsionR3R4;
  delete pbNMRTorsionK2K3;
  delete pbNMRTorsionAve1;
  delete pbNMRTorsionAve2;
  delete pbNMRTorsionTgtVal;
  delete pbNMRTorsionStep;
  delete pbNMRTorsionInc;
  delete pbNMRTorsionR1R2Slp;
  delete pbNMRTorsionR3R4Slp;
  delete pbNMRTorsionK2K3Slp;
  delete pbNMRTorsionR1R2Int;
  delete pbNMRTorsionR3R4Int;
  delete pbNMRTorsionK2K3Int;

  delete pbNMRCOMTorsionID1;
  delete pbNMRCOMTorsionCOM;
  delete pbNMRCOMTorsionCOMGrp;
  delete pbNMRCOMTorsionR1R2;
  delete pbNMRCOMTorsionR3R4;
  delete pbNMRCOMTorsionK2K3;
  delete pbNMRCOMTorsionAve;
  delete pbNMRCOMTorsionTgtVal;
  delete pbNMRCOMTorsionStep;
  delete pbNMRCOMTorsionInc;
  delete pbNMRCOMTorsionR1R2Slp;
  delete pbNMRCOMTorsionR3R4Slp;
  delete pbNMRCOMTorsionK2K3Slp;
  delete pbNMRCOMTorsionR1R2Int;
  delete pbNMRCOMTorsionR3R4Int;
  delete pbNMRCOMTorsionK2K3Int;

  delete pbImageNMRCOMDistanceID;
  delete pbImageNMRCOMDistanceCOM;
  delete pbImageNMRCOMDistanceCOMGrp;
  delete pbImageNMRr6avDistanceID;
  delete pbImageNMRr6avDistancer6av;
  delete pbImageNMRr6avDistancer6avGrp;
  delete pbImageNMRCOMAngleID1;
  delete pbImageNMRCOMAngleID2;
  delete pbImageNMRCOMAngleCOM;
  delete pbImageNMRCOMAngleCOMGrp;
  delete pbImageNMRCOMTorsionID1;
  delete pbImageNMRCOMTorsionCOM;
  delete pbImageNMRCOMTorsionCOMGrp;

  // Delete Shake constraint data
  delete pbShakeID;
  delete pbShakeParm;
  delete pbShakeInvMassH;
  delete pbFastShakeID;
  delete pbSlowShakeID1;
  delete pbSlowShakeID2;
  delete pbSlowShakeParm;
  delete pbSlowShakeInvMassH;

  // Delete extra points data
  delete pbExtraPoint11Frame;
  delete pbExtraPoint11Index;
  delete pbExtraPoint11;
  delete pbExtraPoint12Frame;
  delete pbExtraPoint12Index;
  delete pbExtraPoint12;
  delete pbExtraPoint21Frame;
  delete pbExtraPoint21Index;
  delete pbExtraPoint21;
  delete pbExtraPoint22Frame;
  delete pbExtraPoint22Index;
  delete pbExtraPoint22;
  delete pbImageExtraPoint11Frame;
  delete pbImageExtraPoint11Index;
  delete pbImageExtraPoint12Frame;
  delete pbImageExtraPoint12Index;
  delete pbImageExtraPoint21Frame;
  delete pbImageExtraPoint21Index;
  delete pbImageExtraPoint22Frame;
  delete pbImageExtraPoint22Index;

  // Delete output and/or accumulator buffers
  delete pbReffAccumulator;
  delete pbSumdeijdaAccumulator;
  delete pbForceAccumulator;
  delete pbEnergyBuffer;
  delete pbAFEBuffer;
  delete pbKineticEnergyBuffer;
  delete pbAFEKineticEnergyBuffer;

#ifdef AWSMM
  // Delete AWSMM data
  delete[] pAlphaCarbonIndex;
  delete[] pRefAlphaCarbon;
#endif

  // Delete AMD buffers
  delete[] pAmdWeightsAndEnergy;     // AMD
  delete pbAMDfwgtd;

  // Delete GaMD buffers
  delete[] pGaMDWeightsAndEnergy;     // GaMD
  delete pbGaMDfwgtd;

  // Delete constant pH stuff
  delete pbChargeRefreshBuffer;

  // delete texture object
  cudaDestroyTextureObject(sim.texImageX);
  cudaDestroyTextureObject(sim.texAtomX);
  cudaDestroyTextureObject(sim.texXYZ_q);
  cudaDestroyTextureObject(sim.texOldAtomX);
#if !defined(use_DPFP)
  cudaDestroyTextureObject(sim.texAtomXYSP);
  cudaDestroyTextureObject(sim.texAtomChargeSPLJID);
  cudaDestroyTextureObject(sim.texLJTerm);
  cudaDestroyTextureObject(sim.texAtomZSP);
  cudaDestroyTextureObject(sim.texErfcCoeffsTable);
#endif
}
