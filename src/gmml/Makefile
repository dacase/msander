#############################################################################
# Makefile for building: libgmml.so.1.0.0
#############################################################################
include ../../config.h

MAKEFILE      = Makefile

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = -D_REENTRANT -DQT_DEPRECATED_WARNINGS
CFLAGS        = -pipe -Wall -W -O2 -fPIC $(DEFINES)
CXXFLAGS      = -pipe -Wall -W -std=c++11 -O2 -fPIC $(DEFINES)
INCPATH       = -I. 
DEL_FILE      = rm -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p
COPY          = cp -f
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
INSTALL_FILE  = install -m 644 -p
INSTALL_PROGRAM = install -m 755 -p
INSTALL_DIR   = cp -f -R
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
TAR           = tar -cf
COMPRESS      = gzip -9f
DISTNAME      = gmml1.0.0
LINK          = g++
LFLAGS        = -Wl,-O1 -shared -Wl,-soname,libgmml.so.1
LIBS          = $(SUBLIBS)  
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = strip

####### Files

OBJECTS       = build/codetests.o \
		build/cell.o \
		build/coordinate.o \
		build/geometrytopology.o \
		build/grid.o \
		build/plane.o \
		build/rotation.o \
		build/glycosidiclinkage.o \
		build/monosaccharide.o \
		build/oligosaccharide.o \
		build/element.o \
		build/guesses.o \
		build/atom.o \
		build/atomnode.o \
		build/dockingatom.o \
		build/moleculardynamicatom.o \
		build/molecule.o \
		build/oligosaccharidedetectionatom.o \
		build/overlaps.o \
		build/quantommechanicatom.o \
		build/ring_shape_detection.o \
		build/superimposition.o \
		build/residue_linkage.o \
		build/rotatable_dihedral.o \
		build/ciffile.o \
		build/ciffileatom.o \
		build/ciffileprocessingexception.o \
		build/carbohydratebuilder.o \
		build/carbohydratebuilderjsonparser.o \
		build/condensedsequence.o \
		build/condensedsequenceglycam06residue.o \
		build/condensedsequenceprocessingexception.o \
		build/condensedsequenceresidue.o \
		build/sequencestring.o \
		build/coordinatefile.o \
		build/coordinatefileprocessingexception.o \
		build/pdbatomcard.o \
		build/pdbatomsection.o \
		build/pdbauthorsection.o \
		build/pdbcaveatsection.o \
		build/pdbcispeptidecard.o \
		build/pdbcispeptidesection.o \
		build/pdbcompoundsection.o \
		build/pdbcompoundspecification.o \
		build/pdbconnectsection.o \
		build/pdbcrystallographiccard.o \
		build/pdbdatabasereference.o \
		build/pdbdatabasereferencesection.o \
		build/pdbdisulfidebondsection.o \
		build/pdbdisulfideresidue.o \
		build/pdbdisulfideresiduebond.o \
		build/pdbexperimentaldatasection.o \
		build/pdbfile.o \
		build/pdbfileprocessingexception.o \
		build/pdbformulacard.o \
		build/pdbformulasection.o \
		build/pdbheadercard.o \
		build/pdbhelixcard.o \
		build/pdbhelixresidue.o \
		build/pdbhelixsection.o \
		build/pdbheterogenatomsection.o \
		build/pdbheterogencard.o \
		build/pdbheterogennamecard.o \
		build/pdbheterogennamesection.o \
		build/pdbheterogensection.o \
		build/pdbheterogensynonymcard.o \
		build/pdbheterogensynonymsection.o \
		build/pdbjournalsection.o \
		build/pdbkeywordssection.o \
		build/pdblinkcard.o \
		build/pdblinkcardresidue.o \
		build/pdblinksection.o \
		build/pdbmastercard.o \
		build/pdbmatrixncard.o \
		build/pdbmatrixnsection.o \
		build/pdbmodelcard.o \
		build/pdbmodelresidueset.o \
		build/pdbmodelsection.o \
		build/pdbmodeltypesection.o \
		build/pdbnummodelcard.o \
		build/pdbobsoletecard.o \
		build/pdbobsoletesection.o \
		build/pdboriginxncard.o \
		build/pdboriginxnsection.o \
		build/pdbremarksection.o \
		build/pdbresidue.o \
		build/pdbresiduemodificationcard.o \
		build/pdbresiduemodificationsection.o \
		build/pdbresiduesequencecard.o \
		build/pdbresiduesequencesection.o \
		build/pdbrevisiondatacard.o \
		build/pdbrevisiondatasection.o \
		build/pdbscalencard.o \
		build/pdbscalensection.o \
		build/pdbsequenceadvancedcard.o \
		build/pdbsequenceadvancedsection.o \
		build/pdbsheetcard.o \
		build/pdbsheetsection.o \
		build/pdbsheetstrand.o \
		build/pdbsheetstrandresidue.o \
		build/pdbsitecard.o \
		build/pdbsiteresidue.o \
		build/pdbsitesection.o \
		build/pdbsourcecard.o \
		build/pdbsourcesection.o \
		build/pdbsplitsection.o \
		build/pdbsupersededentriescard.o \
		build/pdbsupersededentriessection.o \
		build/pdbtitlesection.o \
		build/pdbqtatom.o \
		build/pdbqtatomcard.o \
		build/pdbqtbranchcard.o \
		build/pdbqtcompoundcard.o \
		build/pdbqtfile.o \
		build/pdbqtfileprocessingexception.o \
		build/pdbqtmodel.o \
		build/pdbqtmodelcard.o \
		build/pdbqtmodelresidueset.o \
		build/pdbqtremarkcard.o \
		build/pdbqtrootcard.o \
		build/pdbqttorsionaldofcard.o \
		build/topologyangle.o \
		build/topologyangletype.o \
		build/topologyassembly.o \
		build/topologyatom.o \
		build/topologyatompair.o \
		build/topologybond.o \
		build/topologybondtype.o \
		build/topologydihedral.o \
		build/topologydihedraltype.o \
		build/topologyfile.o \
		build/topologyfileprocessingexception.o \
		build/topologyresidue.o \
		build/response.o \
		build/amberatomtypeinfo.o \
		build/amberelements.o \
		build/bondlengthbytypepair.o \
		build/dihedralangledata.o \
		build/glycam06residuecodes.o \
		build/glycam06residueinfo.o \
		build/assembly.o \
		build/assemblybuilder.o \
		build/geometryoperation.o \
		build/glycamnaming.o \
		build/ionization.o \
		build/manipulationoperation.o \
		build/selection.o \
		build/solvation.o \
		build/structurebuilder.o \
		build/edge.o \
		build/graph.o \
		build/node.o \
		build/residue.o \
		build/residue3dtemplatelibrary.o \
		build/residuebuilder.o \
		build/residuenode.o \
		build/residueproperties.o \
		build/selections.o \
		build/libraryfile.o \
		build/libraryfileatom.o \
		build/libraryfileprocessingexception.o \
		build/libraryfileresidue.o \
		build/offfile.o \
		build/offfileatom.o \
		build/offfileresidue.o \
		build/parameterfile.o \
		build/parameterfileangle.o \
		build/parameterfileatom.o \
		build/parameterfilebond.o \
		build/parameterfiledihedral.o \
		build/parameterfiledihedralterm.o \
		build/parameterfileprocessingexception.o \
		build/prepfile.o \
		build/prepfileatom.o \
		build/prepfileprocessingexception.o \
		build/prepfileresidue.o \
		build/pdbpreprocessor.o \
		build/pdbpreprocessoralternateresidue.o \
		build/pdbpreprocessorchaintermination.o \
		build/pdbpreprocessordisulfidebond.o \
		build/pdbpreprocessorhistidinemapping.o \
		build/pdbpreprocessormissingresidue.o \
		build/pdbpreprocessorreplacedhydrogen.o \
		build/pdbpreprocessorresidueinfo.o \
		build/pdbpreprocessorunrecognizedheavyatom.o \
		build/pdbpreprocessorunrecognizedresidue.o \
		build/analysis.o \
		build/gf_query.o \
		build/population.o \
		build/query.o \
		build/coordinatebuilder.o \
		build/librarybuilder.o \
		build/pdbbuilder.o \
		build/pdbqtbuilder.o \
		build/prepbuilder.o \
		build/topologybuilder.o \
		build/oligosaccharidedetection.o \
		build/ringperception.o

TARGETA       = $(LIBDIR)/libgmml.a

####### Build rules

gmml: $(TARGETA) $(BINDIR)/pdb2glycam

$(BINDIR)/pdb2glycam: $(TARGETA) pdb2glycam.cc
	@echo "[GMML] CXX $@"
	$(VB)$(CXX) -std=c++0x -I./includes/ -L$(LIBDIR) pdb2glycam.cc \
        -lgmml -lpthread -o $@

$(TARGETA):  $(OBJECTS) $(OBJCOMP) 
	@$(DEL_FILE) $(TARGETA) 
	@echo "[GMML] AR $@"
	$(VB)$(AR) $(TARGETA) $(OBJECTS)

clean:
	$(VB)$(DEL_FILE) $(OBJECTS)

distclean: clean 
	$(VB)$(DEL_FILE) $(TARGETA)

####### Compile

build/codetests.o: src/CodeUtils/codetests.cc includes/CodeUtils/codetests.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/codetests.o src/CodeUtils/codetests.cc

build/cell.o: src/GeometryTopology/cell.cc includes/GeometryTopology/cell.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/cell.o src/GeometryTopology/cell.cc

build/coordinate.o: src/GeometryTopology/coordinate.cc includes/GeometryTopology/coordinate.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/coordinate.o src/GeometryTopology/coordinate.cc

build/geometrytopology.o: src/GeometryTopology/geometrytopology.cc includes/GeometryTopology/geometrytopology.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/External_Libraries/Eigen_Algebra_Template_Library/Geometry \
		includes/External_Libraries/Eigen_Algebra_Template_Library/Core \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/DisableStupidWarnings.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/Macros.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/MKL_support.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/misc/blas.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/Constants.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/Meta.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/ForwardDeclarations.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/StaticAssert.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/XprHelper.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/Memory.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/NumTraits.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/GenericPacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/MathFunctionsImpl.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/SSE/PacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AVX/PacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AVX512/PacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AVX512/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/SSE/Complex.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/SSE/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AVX/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AVX/Complex.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AVX/TypeCasting.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/SSE/TypeCasting.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AltiVec/PacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AltiVec/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AltiVec/Complex.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/NEON/PacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/NEON/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/NEON/Complex.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/ZVector/PacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/ZVector/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/ZVector/Complex.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/Half.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/PacketMathHalf.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/TypeCasting.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/PacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/Default/Settings.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/functors/TernaryFunctors.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/functors/BinaryFunctors.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/functors/UnaryFunctors.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/functors/NullaryFunctors.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/functors/StlFunctors.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/functors/AssignmentFunctors.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/Complex.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/IO.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/DenseCoeffsBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/DenseBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/plugins/BlockMethods.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/MatrixBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/plugins/CommonCwiseUnaryOps.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/plugins/CommonCwiseBinaryOps.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/plugins/MatrixCwiseUnaryOps.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/plugins/MatrixCwiseBinaryOps.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/EigenBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Product.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CoreEvaluators.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/AssignEvaluator.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Assign.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/ArrayBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/plugins/ArrayCwiseUnaryOps.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/plugins/ArrayCwiseBinaryOps.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/BlasUtil.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/DenseStorage.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/NestByValue.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/ReturnByValue.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/NoAlias.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/PlainObjectBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Matrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Array.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CwiseTernaryOp.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CwiseBinaryOp.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CwiseUnaryOp.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CwiseNullaryOp.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CwiseUnaryView.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/SelfCwiseBinaryOp.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Dot.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/StableNorm.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Stride.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/MapBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Map.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Ref.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Block.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/VectorBlock.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Transpose.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/DiagonalMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Diagonal.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/DiagonalProduct.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Redux.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Visitor.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Fuzzy.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Swap.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CommaInitializer.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/GeneralProduct.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Solve.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Inverse.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/SolverBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/PermutationMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Transpositions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/TriangularMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/SelfAdjointView.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/GeneralBlockPanelKernel.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/Parallelizer.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/ProductEvaluators.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixVector.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/SolveTriangular.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixMatrixTriangular.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointMatrixVector.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointMatrixMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointProduct.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointRank2Update.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/TriangularMatrixVector.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/TriangularMatrixMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/TriangularSolverMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/TriangularSolverVector.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/BandMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CoreIterators.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/ConditionEstimator.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/BooleanRedux.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Select.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/VectorwiseOp.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Random.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Replicate.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Reverse.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/ArrayWrapper.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixMatrix_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixVector_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixMatrixTriangular_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointMatrixMatrix_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointMatrixVector_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/TriangularMatrixMatrix_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/TriangularMatrixVector_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/TriangularSolverMatrix_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Assign_MKL.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/GlobalFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/ReenableStupidWarnings.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/SVD \
		includes/External_Libraries/Eigen_Algebra_Template_Library/QR \
		includes/External_Libraries/Eigen_Algebra_Template_Library/Cholesky \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Cholesky/LLT.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Cholesky/LDLT.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/misc/lapacke.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/misc/lapacke_mangling.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Cholesky/LLT_LAPACKE.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/Jacobi \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Jacobi/Jacobi.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/Householder \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Householder/Householder.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Householder/HouseholderSequence.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Householder/BlockHouseholder.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/QR/HouseholderQR.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/QR/FullPivHouseholderQR.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/QR/ColPivHouseholderQR.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/QR/CompleteOrthogonalDecomposition.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/QR/HouseholderQR_LAPACKE.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/QR/ColPivHouseholderQR_LAPACKE.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/misc/RealSvd2x2.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/SVD/UpperBidiagonalization.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/SVD/SVDBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/SVD/JacobiSVD.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/SVD/BDCSVD.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/SVD/JacobiSVD_LAPACKE.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/LU \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/misc/Kernel.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/misc/Image.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/LU/FullPivLU.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/LU/PartialPivLU.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/LU/PartialPivLU_LAPACKE.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/LU/Determinant.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/LU/InverseImpl.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/LU/arch/Inverse_SSE.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/OrthoMethods.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/EulerAngles.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Homogeneous.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/RotationBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Rotation2D.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Quaternion.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/AngleAxis.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Transform.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Translation.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Scaling.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Hyperplane.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/ParametrizedLine.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/AlignedBox.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Umeyama.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/arch/Geometry_SSE.h \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp
	@echo "[GMML] CXX $<"
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/geometrytopology.o src/GeometryTopology/geometrytopology.cc

build/grid.o: src/GeometryTopology/grid.cc includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/cell.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/grid.o src/GeometryTopology/grid.cc

build/plane.o: src/GeometryTopology/plane.cc includes/GeometryTopology/plane.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/plane.o src/GeometryTopology/plane.cc

build/rotation.o: src/GeometryTopology/rotation.cc includes/GeometryTopology/rotation.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/rotation.o src/GeometryTopology/rotation.cc

build/glycosidiclinkage.o: src/Glycan/glycosidiclinkage.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/glycosidiclinkage.o src/Glycan/glycosidiclinkage.cc

build/monosaccharide.o: src/Glycan/monosaccharide.cc includes/MolecularModeling/atom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/sugarname.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/generictypedefs.hpp \
		includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/MolecularModeling/ring_shape_detection.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/monosaccharide.o src/Glycan/monosaccharide.cc

build/oligosaccharide.o: src/Glycan/oligosaccharide.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/oligosaccharide.o src/Glycan/oligosaccharide.cc

build/element.o: src/MolecularMetadata/element.cc includes/MolecularMetadata/element.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/element.o src/MolecularMetadata/element.cc

build/guesses.o: src/MolecularMetadata/guesses.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/guesses.o src/MolecularMetadata/guesses.cc

build/atom.o: src/MolecularModeling/atom.cc includes/MolecularModeling/atom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/generictypedefs.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/common.hpp \
		includes/Glycan/sugarname.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/atom.o src/MolecularModeling/atom.cc

build/atomnode.o: src/MolecularModeling/atomnode.cc includes/MolecularModeling/atomnode.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/generictypedefs.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/common.hpp \
		includes/Glycan/sugarname.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/utils.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/atomnode.o src/MolecularModeling/atomnode.cc

build/dockingatom.o: src/MolecularModeling/dockingatom.cc includes/MolecularModeling/dockingatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/dockingatom.o src/MolecularModeling/dockingatom.cc

build/moleculardynamicatom.o: src/MolecularModeling/moleculardynamicatom.cc includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/moleculardynamicatom.o src/MolecularModeling/moleculardynamicatom.cc

build/molecule.o: src/MolecularModeling/molecule.cc includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/molecule.o src/MolecularModeling/molecule.cc

build/oligosaccharidedetectionatom.o: src/MolecularModeling/oligosaccharidedetectionatom.cc includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/oligosaccharidedetectionatom.o src/MolecularModeling/oligosaccharidedetectionatom.cc

build/overlaps.o: src/MolecularModeling/overlaps.cc includes/MolecularModeling/overlaps.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/generictypedefs.hpp \
		includes/common.hpp \
		includes/Glycan/sugarname.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/overlaps.o src/MolecularModeling/overlaps.cc

build/quantommechanicatom.o: src/MolecularModeling/quantommechanicatom.cc includes/MolecularModeling/quantommechanicatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/quantommechanicatom.o src/MolecularModeling/quantommechanicatom.cc

build/ring_shape_detection.o: src/MolecularModeling/ring_shape_detection.cc includes/MolecularModeling/ring_shape_detection.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/generictypedefs.hpp \
		includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/Glycan/sugarname.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/Glycan/note.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/ring_shape_detection.o src/MolecularModeling/ring_shape_detection.cc

build/superimposition.o: src/MolecularModeling/superimposition.cc includes/MolecularModeling/superimposition.hpp \
		includes/External_Libraries/Eigen_Algebra_Template_Library/Geometry \
		includes/External_Libraries/Eigen_Algebra_Template_Library/Core \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/DisableStupidWarnings.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/Macros.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/MKL_support.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/misc/blas.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/Constants.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/Meta.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/ForwardDeclarations.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/StaticAssert.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/XprHelper.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/Memory.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/NumTraits.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/GenericPacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/MathFunctionsImpl.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/SSE/PacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AVX/PacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AVX512/PacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AVX512/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/SSE/Complex.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/SSE/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AVX/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AVX/Complex.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AVX/TypeCasting.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/SSE/TypeCasting.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AltiVec/PacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AltiVec/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/AltiVec/Complex.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/NEON/PacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/NEON/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/NEON/Complex.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/ZVector/PacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/ZVector/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/ZVector/Complex.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/Half.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/PacketMathHalf.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/TypeCasting.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/PacketMath.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/MathFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/Default/Settings.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/functors/TernaryFunctors.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/functors/BinaryFunctors.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/functors/UnaryFunctors.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/functors/NullaryFunctors.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/functors/StlFunctors.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/functors/AssignmentFunctors.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/arch/CUDA/Complex.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/IO.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/DenseCoeffsBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/DenseBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/plugins/BlockMethods.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/MatrixBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/plugins/CommonCwiseUnaryOps.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/plugins/CommonCwiseBinaryOps.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/plugins/MatrixCwiseUnaryOps.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/plugins/MatrixCwiseBinaryOps.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/EigenBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Product.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CoreEvaluators.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/AssignEvaluator.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Assign.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/ArrayBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/plugins/ArrayCwiseUnaryOps.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/plugins/ArrayCwiseBinaryOps.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/BlasUtil.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/DenseStorage.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/NestByValue.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/ReturnByValue.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/NoAlias.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/PlainObjectBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Matrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Array.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CwiseTernaryOp.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CwiseBinaryOp.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CwiseUnaryOp.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CwiseNullaryOp.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CwiseUnaryView.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/SelfCwiseBinaryOp.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Dot.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/StableNorm.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Stride.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/MapBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Map.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Ref.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Block.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/VectorBlock.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Transpose.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/DiagonalMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Diagonal.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/DiagonalProduct.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Redux.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Visitor.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Fuzzy.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Swap.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CommaInitializer.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/GeneralProduct.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Solve.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Inverse.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/SolverBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/PermutationMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Transpositions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/TriangularMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/SelfAdjointView.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/GeneralBlockPanelKernel.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/Parallelizer.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/ProductEvaluators.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixVector.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/SolveTriangular.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixMatrixTriangular.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointMatrixVector.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointMatrixMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointProduct.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointRank2Update.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/TriangularMatrixVector.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/TriangularMatrixMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/TriangularSolverMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/TriangularSolverVector.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/BandMatrix.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/CoreIterators.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/ConditionEstimator.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/BooleanRedux.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Select.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/VectorwiseOp.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Random.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Replicate.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Reverse.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/ArrayWrapper.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixMatrix_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixVector_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/GeneralMatrixMatrixTriangular_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointMatrixMatrix_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/SelfadjointMatrixVector_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/TriangularMatrixMatrix_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/TriangularMatrixVector_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/products/TriangularSolverMatrix_BLAS.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/Assign_MKL.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/GlobalFunctions.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Core/util/ReenableStupidWarnings.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/SVD \
		includes/External_Libraries/Eigen_Algebra_Template_Library/QR \
		includes/External_Libraries/Eigen_Algebra_Template_Library/Cholesky \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Cholesky/LLT.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Cholesky/LDLT.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/misc/lapacke.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/misc/lapacke_mangling.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Cholesky/LLT_LAPACKE.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/Jacobi \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Jacobi/Jacobi.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/Householder \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Householder/Householder.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Householder/HouseholderSequence.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Householder/BlockHouseholder.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/QR/HouseholderQR.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/QR/FullPivHouseholderQR.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/QR/ColPivHouseholderQR.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/QR/CompleteOrthogonalDecomposition.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/QR/HouseholderQR_LAPACKE.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/QR/ColPivHouseholderQR_LAPACKE.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/misc/RealSvd2x2.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/SVD/UpperBidiagonalization.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/SVD/SVDBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/SVD/JacobiSVD.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/SVD/BDCSVD.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/SVD/JacobiSVD_LAPACKE.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/LU \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/misc/Kernel.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/misc/Image.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/LU/FullPivLU.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/LU/PartialPivLU.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/LU/PartialPivLU_LAPACKE.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/LU/Determinant.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/LU/InverseImpl.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/LU/arch/Inverse_SSE.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/OrthoMethods.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/EulerAngles.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Homogeneous.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/RotationBase.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Rotation2D.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Quaternion.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/AngleAxis.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Transform.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Translation.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Scaling.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Hyperplane.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/ParametrizedLine.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/AlignedBox.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/Umeyama.h \
		includes/External_Libraries/Eigen_Algebra_Template_Library/src/Geometry/arch/Geometry_SSE.h \
		includes/GeometryTopology/coordinate.hpp \
		includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/superimposition.o src/MolecularModeling/superimposition.cc

build/residue_linkage.o: src/GeometryTopology/ResidueLinkages/residue_linkage.cc includes/GeometryTopology/ResidueLinkages/residue_linkage.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/MolecularModeling/Selections/selections.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/GeometryTopology/ResidueLinkages/rotatable_dihedral.hpp \
		includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp \
		includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp \
		includes/MolecularModeling/overlaps.hpp \
		includes/External_Libraries/PCG/pcg_extras.hpp \
		includes/External_Libraries/PCG/pcg_random.hpp \
		includes/MolecularMetadata/GLYCAM/glycam06residuecodes.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/residue_linkage.o src/GeometryTopology/ResidueLinkages/residue_linkage.cc

build/rotatable_dihedral.o: src/GeometryTopology/ResidueLinkages/rotatable_dihedral.cc includes/GeometryTopology/ResidueLinkages/rotatable_dihedral.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp \
		includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/utils.hpp \
		includes/External_Libraries/PCG/pcg_random.hpp \
		includes/External_Libraries/PCG/pcg_extras.hpp \
		includes/GeometryTopology/geometrytopology.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/rotatable_dihedral.o src/GeometryTopology/ResidueLinkages/rotatable_dihedral.cc

build/ciffile.o: src/InputSet/CifFileSpace/ciffile.cc includes/InputSet/CifFileSpace/ciffile.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/InputSet/CifFileSpace/ciffileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/ciffile.o src/InputSet/CifFileSpace/ciffile.cc

build/ciffileatom.o: src/InputSet/CifFileSpace/ciffileatom.cc includes/InputSet/CifFileSpace/ciffileatom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/ciffileatom.o src/InputSet/CifFileSpace/ciffileatom.cc

build/ciffileprocessingexception.o: src/InputSet/CifFileSpace/ciffileprocessingexception.cc includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp \
		includes/InputSet/CifFileSpace/ciffileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/ciffileprocessingexception.o src/InputSet/CifFileSpace/ciffileprocessingexception.cc

build/carbohydratebuilder.o: src/InputSet/CondensedSequenceSpace/carbohydratebuilder.cc includes/InputSet/CondensedSequenceSpace/carbohydratebuilder.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/note.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/GeometryTopology/ResidueLinkages/residue_linkage.hpp \
		includes/MolecularModeling/Selections/selections.hpp \
		includes/GeometryTopology/ResidueLinkages/rotatable_dihedral.hpp \
		includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp \
		includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp \
		includes/External_Libraries/json.hpp \
		includes/ParameterSet/OffFileSpace/offfile.hpp \
		includes/MolecularMetadata/elementattributes.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/carbohydratebuilder.o src/InputSet/CondensedSequenceSpace/carbohydratebuilder.cc

build/carbohydratebuilderjsonparser.o: src/InputSet/CondensedSequenceSpace/carbohydratebuilderjsonparser.cc includes/InputSet/CondensedSequenceSpace/carbohydratebuilder.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/note.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/External_Libraries/json.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/carbohydratebuilderjsonparser.o src/InputSet/CondensedSequenceSpace/carbohydratebuilderjsonparser.cc

build/condensedsequence.o: src/InputSet/CondensedSequenceSpace/condensedsequence.cc includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/note.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.hpp \
		includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/condensedsequence.o src/InputSet/CondensedSequenceSpace/condensedsequence.cc

build/condensedsequenceglycam06residue.o: src/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.cc includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/condensedsequenceglycam06residue.o src/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.cc

build/condensedsequenceprocessingexception.o: src/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.cc includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/condensedsequenceprocessingexception.o src/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.cc

build/condensedsequenceresidue.o: src/InputSet/CondensedSequenceSpace/condensedsequenceresidue.cc includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/note.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.hpp \
		includes/utils.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/condensedsequenceresidue.o src/InputSet/CondensedSequenceSpace/condensedsequenceresidue.cc

build/sequencestring.o: src/InputSet/CondensedSequenceSpace/sequencestring.cpp includes/InputSet/CondensedSequenceSpace/sequencestring.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/sequencestring.o src/InputSet/CondensedSequenceSpace/sequencestring.cpp

build/coordinatefile.o: src/InputSet/CoordinateFileSpace/coordinatefile.cc includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/coordinatefile.o src/InputSet/CoordinateFileSpace/coordinatefile.cc

build/coordinatefileprocessingexception.o: src/InputSet/CoordinateFileSpace/coordinatefileprocessingexception.cc includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/coordinatefileprocessingexception.o src/InputSet/CoordinateFileSpace/coordinatefileprocessingexception.cc

build/pdbatomcard.o: src/InputSet/PdbFileSpace/pdbatomcard.cc includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbatomcard.o src/InputSet/PdbFileSpace/pdbatomcard.cc

build/pdbatomsection.o: src/InputSet/PdbFileSpace/pdbatomsection.cc includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbatomsection.o src/InputSet/PdbFileSpace/pdbatomsection.cc

build/pdbauthorsection.o: src/InputSet/PdbFileSpace/pdbauthorsection.cc includes/InputSet/PdbFileSpace/pdbauthorsection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbauthorsection.o src/InputSet/PdbFileSpace/pdbauthorsection.cc

build/pdbcaveatsection.o: src/InputSet/PdbFileSpace/pdbcaveatsection.cc includes/InputSet/PdbFileSpace/pdbcaveatsection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbcaveatsection.o src/InputSet/PdbFileSpace/pdbcaveatsection.cc

build/pdbcispeptidecard.o: src/InputSet/PdbFileSpace/pdbcispeptidecard.cc includes/InputSet/PdbFileSpace/pdbcispeptidecard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbcispeptidecard.o src/InputSet/PdbFileSpace/pdbcispeptidecard.cc

build/pdbcispeptidesection.o: src/InputSet/PdbFileSpace/pdbcispeptidesection.cc includes/InputSet/PdbFileSpace/pdbcispeptidesection.hpp \
		includes/InputSet/PdbFileSpace/pdbcispeptidecard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbcispeptidesection.o src/InputSet/PdbFileSpace/pdbcispeptidesection.cc

build/pdbcompoundsection.o: src/InputSet/PdbFileSpace/pdbcompoundsection.cc includes/InputSet/PdbFileSpace/pdbcompoundsection.hpp \
		includes/InputSet/PdbFileSpace/pdbcompoundspecification.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbcompoundsection.o src/InputSet/PdbFileSpace/pdbcompoundsection.cc

build/pdbcompoundspecification.o: src/InputSet/PdbFileSpace/pdbcompoundspecification.cc includes/InputSet/PdbFileSpace/pdbcompoundspecification.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbcompoundspecification.o src/InputSet/PdbFileSpace/pdbcompoundspecification.cc

build/pdbconnectsection.o: src/InputSet/PdbFileSpace/pdbconnectsection.cc includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbconnectsection.o src/InputSet/PdbFileSpace/pdbconnectsection.cc

build/pdbcrystallographiccard.o: src/InputSet/PdbFileSpace/pdbcrystallographiccard.cc includes/InputSet/PdbFileSpace/pdbcrystallographiccard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbcrystallographiccard.o src/InputSet/PdbFileSpace/pdbcrystallographiccard.cc

build/pdbdatabasereference.o: src/InputSet/PdbFileSpace/pdbdatabasereference.cc includes/InputSet/PdbFileSpace/pdbdatabasereference.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbdatabasereference.o src/InputSet/PdbFileSpace/pdbdatabasereference.cc

build/pdbdatabasereferencesection.o: src/InputSet/PdbFileSpace/pdbdatabasereferencesection.cc includes/InputSet/PdbFileSpace/pdbdatabasereferencesection.hpp \
		includes/InputSet/PdbFileSpace/pdbdatabasereference.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbdatabasereferencesection.o src/InputSet/PdbFileSpace/pdbdatabasereferencesection.cc

build/pdbdisulfidebondsection.o: src/InputSet/PdbFileSpace/pdbdisulfidebondsection.cc includes/InputSet/PdbFileSpace/pdbdisulfidebondsection.hpp \
		includes/InputSet/PdbFileSpace/pdbdisulfideresiduebond.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbdisulfidebondsection.o src/InputSet/PdbFileSpace/pdbdisulfidebondsection.cc

build/pdbdisulfideresidue.o: src/InputSet/PdbFileSpace/pdbdisulfideresidue.cc includes/InputSet/PdbFileSpace/pdbdisulfideresidue.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbdisulfideresidue.o src/InputSet/PdbFileSpace/pdbdisulfideresidue.cc

build/pdbdisulfideresiduebond.o: src/InputSet/PdbFileSpace/pdbdisulfideresiduebond.cc includes/InputSet/PdbFileSpace/pdbdisulfideresiduebond.hpp \
		includes/InputSet/PdbFileSpace/pdbdisulfideresidue.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbdisulfideresiduebond.o src/InputSet/PdbFileSpace/pdbdisulfideresiduebond.cc

build/pdbexperimentaldatasection.o: src/InputSet/PdbFileSpace/pdbexperimentaldatasection.cc includes/InputSet/PdbFileSpace/pdbexperimentaldatasection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbexperimentaldatasection.o src/InputSet/PdbFileSpace/pdbexperimentaldatasection.cc

build/pdbfile.o: src/InputSet/PdbFileSpace/pdbfile.cc includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbFileSpace/pdbheadercard.hpp \
		includes/InputSet/PdbFileSpace/pdbobsoletesection.hpp \
		includes/InputSet/PdbFileSpace/pdbobsoletecard.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbsplitsection.hpp \
		includes/InputSet/PdbFileSpace/pdbcaveatsection.hpp \
		includes/InputSet/PdbFileSpace/pdbcompoundsection.hpp \
		includes/InputSet/PdbFileSpace/pdbcompoundspecification.hpp \
		includes/InputSet/PdbFileSpace/pdbsourcesection.hpp \
		includes/InputSet/PdbFileSpace/pdbsourcecard.hpp \
		includes/InputSet/PdbFileSpace/pdbkeywordssection.hpp \
		includes/InputSet/PdbFileSpace/pdbexperimentaldatasection.hpp \
		includes/InputSet/PdbFileSpace/pdbnummodelcard.hpp \
		includes/InputSet/PdbFileSpace/pdbmodeltypesection.hpp \
		includes/InputSet/PdbFileSpace/pdbauthorsection.hpp \
		includes/InputSet/PdbFileSpace/pdbrevisiondatasection.hpp \
		includes/InputSet/PdbFileSpace/pdbrevisiondatacard.hpp \
		includes/InputSet/PdbFileSpace/pdbsupersededentriessection.hpp \
		includes/InputSet/PdbFileSpace/pdbsupersededentriescard.hpp \
		includes/InputSet/PdbFileSpace/pdbjournalsection.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbFileSpace/pdbdatabasereferencesection.hpp \
		includes/InputSet/PdbFileSpace/pdbdatabasereference.hpp \
		includes/InputSet/PdbFileSpace/pdbsequenceadvancedsection.hpp \
		includes/InputSet/PdbFileSpace/pdbsequenceadvancedcard.hpp \
		includes/InputSet/PdbFileSpace/pdbresiduesequencesection.hpp \
		includes/InputSet/PdbFileSpace/pdbresiduesequencecard.hpp \
		includes/InputSet/PdbFileSpace/pdbresiduemodificationsection.hpp \
		includes/InputSet/PdbFileSpace/pdbresiduemodificationcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogensection.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogencard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogennamesection.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogennamecard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogensynonymsection.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogensynonymcard.hpp \
		includes/InputSet/PdbFileSpace/pdbformulasection.hpp \
		includes/InputSet/PdbFileSpace/pdbformulacard.hpp \
		includes/InputSet/PdbFileSpace/pdbhelixsection.hpp \
		includes/InputSet/PdbFileSpace/pdbhelixcard.hpp \
		includes/InputSet/PdbFileSpace/pdbhelixresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbsheetsection.hpp \
		includes/InputSet/PdbFileSpace/pdbsheetcard.hpp \
		includes/InputSet/PdbFileSpace/pdbsheetstrand.hpp \
		includes/InputSet/PdbFileSpace/pdbsheetstrandresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbdisulfidebondsection.hpp \
		includes/InputSet/PdbFileSpace/pdbdisulfideresiduebond.hpp \
		includes/InputSet/PdbFileSpace/pdbdisulfideresidue.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbcispeptidesection.hpp \
		includes/InputSet/PdbFileSpace/pdbcispeptidecard.hpp \
		includes/InputSet/PdbFileSpace/pdbsitesection.hpp \
		includes/InputSet/PdbFileSpace/pdbsitecard.hpp \
		includes/InputSet/PdbFileSpace/pdbsiteresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbcrystallographiccard.hpp \
		includes/InputSet/PdbFileSpace/pdboriginxncard.hpp \
		includes/InputSet/PdbFileSpace/pdboriginxnsection.hpp \
		includes/InputSet/PdbFileSpace/pdbscalencard.hpp \
		includes/InputSet/PdbFileSpace/pdbscalensection.hpp \
		includes/InputSet/PdbFileSpace/pdbmatrixnsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmatrixncard.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmastercard.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbFileSpace/pdbresidue.hpp \
		includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbfile.o src/InputSet/PdbFileSpace/pdbfile.cc

build/pdbfileprocessingexception.o: src/InputSet/PdbFileSpace/pdbfileprocessingexception.cc includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbfileprocessingexception.o src/InputSet/PdbFileSpace/pdbfileprocessingexception.cc

build/pdbformulacard.o: src/InputSet/PdbFileSpace/pdbformulacard.cc includes/InputSet/PdbFileSpace/pdbformulacard.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbformulacard.o src/InputSet/PdbFileSpace/pdbformulacard.cc

build/pdbformulasection.o: src/InputSet/PdbFileSpace/pdbformulasection.cc includes/InputSet/PdbFileSpace/pdbformulacard.hpp \
		includes/InputSet/PdbFileSpace/pdbformulasection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbformulasection.o src/InputSet/PdbFileSpace/pdbformulasection.cc

build/pdbheadercard.o: src/InputSet/PdbFileSpace/pdbheadercard.cc includes/InputSet/PdbFileSpace/pdbheadercard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbheadercard.o src/InputSet/PdbFileSpace/pdbheadercard.cc

build/pdbhelixcard.o: src/InputSet/PdbFileSpace/pdbhelixcard.cc includes/InputSet/PdbFileSpace/pdbhelixcard.hpp \
		includes/InputSet/PdbFileSpace/pdbhelixresidue.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbhelixcard.o src/InputSet/PdbFileSpace/pdbhelixcard.cc

build/pdbhelixresidue.o: src/InputSet/PdbFileSpace/pdbhelixresidue.cc includes/InputSet/PdbFileSpace/pdbhelixresidue.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbhelixresidue.o src/InputSet/PdbFileSpace/pdbhelixresidue.cc

build/pdbhelixsection.o: src/InputSet/PdbFileSpace/pdbhelixsection.cc includes/InputSet/PdbFileSpace/pdbhelixcard.hpp \
		includes/InputSet/PdbFileSpace/pdbhelixsection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbhelixsection.o src/InputSet/PdbFileSpace/pdbhelixsection.cc

build/pdbheterogenatomsection.o: src/InputSet/PdbFileSpace/pdbheterogenatomsection.cc includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbheterogenatomsection.o src/InputSet/PdbFileSpace/pdbheterogenatomsection.cc

build/pdbheterogencard.o: src/InputSet/PdbFileSpace/pdbheterogencard.cc includes/InputSet/PdbFileSpace/pdbheterogencard.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbheterogencard.o src/InputSet/PdbFileSpace/pdbheterogencard.cc

build/pdbheterogennamecard.o: src/InputSet/PdbFileSpace/pdbheterogennamecard.cc includes/InputSet/PdbFileSpace/pdbheterogennamecard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbheterogennamecard.o src/InputSet/PdbFileSpace/pdbheterogennamecard.cc

build/pdbheterogennamesection.o: src/InputSet/PdbFileSpace/pdbheterogennamesection.cc includes/InputSet/PdbFileSpace/pdbheterogennamecard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogennamesection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbheterogennamesection.o src/InputSet/PdbFileSpace/pdbheterogennamesection.cc

build/pdbheterogensection.o: src/InputSet/PdbFileSpace/pdbheterogensection.cc includes/InputSet/PdbFileSpace/pdbheterogencard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogensection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbheterogensection.o src/InputSet/PdbFileSpace/pdbheterogensection.cc

build/pdbheterogensynonymcard.o: src/InputSet/PdbFileSpace/pdbheterogensynonymcard.cc includes/InputSet/PdbFileSpace/pdbheterogensynonymcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbheterogensynonymcard.o src/InputSet/PdbFileSpace/pdbheterogensynonymcard.cc

build/pdbheterogensynonymsection.o: src/InputSet/PdbFileSpace/pdbheterogensynonymsection.cc includes/InputSet/PdbFileSpace/pdbheterogensynonymcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogensynonymsection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbheterogensynonymsection.o src/InputSet/PdbFileSpace/pdbheterogensynonymsection.cc

build/pdbjournalsection.o: src/InputSet/PdbFileSpace/pdbjournalsection.cc includes/InputSet/PdbFileSpace/pdbjournalsection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbjournalsection.o src/InputSet/PdbFileSpace/pdbjournalsection.cc

build/pdbkeywordssection.o: src/InputSet/PdbFileSpace/pdbkeywordssection.cc includes/InputSet/PdbFileSpace/pdbkeywordssection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbkeywordssection.o src/InputSet/PdbFileSpace/pdbkeywordssection.cc

build/pdblinkcard.o: src/InputSet/PdbFileSpace/pdblinkcard.cc includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdblinkcard.o src/InputSet/PdbFileSpace/pdblinkcard.cc

build/pdblinkcardresidue.o: src/InputSet/PdbFileSpace/pdblinkcardresidue.cc includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdblinkcardresidue.o src/InputSet/PdbFileSpace/pdblinkcardresidue.cc

build/pdblinksection.o: src/InputSet/PdbFileSpace/pdblinksection.cc includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdblinksection.o src/InputSet/PdbFileSpace/pdblinksection.cc

build/pdbmastercard.o: src/InputSet/PdbFileSpace/pdbmastercard.cc includes/InputSet/PdbFileSpace/pdbmastercard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbmastercard.o src/InputSet/PdbFileSpace/pdbmastercard.cc

build/pdbmatrixncard.o: src/InputSet/PdbFileSpace/pdbmatrixncard.cc includes/InputSet/PdbFileSpace/pdbmatrixncard.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbmatrixncard.o src/InputSet/PdbFileSpace/pdbmatrixncard.cc

build/pdbmatrixnsection.o: src/InputSet/PdbFileSpace/pdbmatrixnsection.cc includes/InputSet/PdbFileSpace/pdbmatrixnsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmatrixncard.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbmatrixnsection.o src/InputSet/PdbFileSpace/pdbmatrixnsection.cc

build/pdbmodelcard.o: src/InputSet/PdbFileSpace/pdbmodelcard.cc includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbmodelcard.o src/InputSet/PdbFileSpace/pdbmodelcard.cc

build/pdbmodelresidueset.o: src/InputSet/PdbFileSpace/pdbmodelresidueset.cc includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbmodelresidueset.o src/InputSet/PdbFileSpace/pdbmodelresidueset.cc

build/pdbmodelsection.o: src/InputSet/PdbFileSpace/pdbmodelsection.cc includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbmodelsection.o src/InputSet/PdbFileSpace/pdbmodelsection.cc

build/pdbmodeltypesection.o: src/InputSet/PdbFileSpace/pdbmodeltypesection.cc includes/InputSet/PdbFileSpace/pdbmodeltypesection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbmodeltypesection.o src/InputSet/PdbFileSpace/pdbmodeltypesection.cc

build/pdbnummodelcard.o: src/InputSet/PdbFileSpace/pdbnummodelcard.cc includes/InputSet/PdbFileSpace/pdbnummodelcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbnummodelcard.o src/InputSet/PdbFileSpace/pdbnummodelcard.cc

build/pdbobsoletecard.o: src/InputSet/PdbFileSpace/pdbobsoletecard.cc includes/InputSet/PdbFileSpace/pdbobsoletecard.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbobsoletecard.o src/InputSet/PdbFileSpace/pdbobsoletecard.cc

build/pdbobsoletesection.o: src/InputSet/PdbFileSpace/pdbobsoletesection.cc includes/InputSet/PdbFileSpace/pdbobsoletecard.hpp \
		includes/InputSet/PdbFileSpace/pdbobsoletesection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbobsoletesection.o src/InputSet/PdbFileSpace/pdbobsoletesection.cc

build/pdboriginxncard.o: src/InputSet/PdbFileSpace/pdboriginxncard.cc includes/InputSet/PdbFileSpace/pdboriginxncard.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdboriginxncard.o src/InputSet/PdbFileSpace/pdboriginxncard.cc

build/pdboriginxnsection.o: src/InputSet/PdbFileSpace/pdboriginxnsection.cc includes/InputSet/PdbFileSpace/pdboriginxnsection.hpp \
		includes/InputSet/PdbFileSpace/pdboriginxncard.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdboriginxnsection.o src/InputSet/PdbFileSpace/pdboriginxnsection.cc

build/pdbremarksection.o: src/InputSet/PdbFileSpace/pdbremarksection.cc includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbremarksection.o src/InputSet/PdbFileSpace/pdbremarksection.cc

build/pdbresidue.o: src/InputSet/PdbFileSpace/pdbresidue.cc includes/InputSet/PdbFileSpace/pdbresidue.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbresidue.o src/InputSet/PdbFileSpace/pdbresidue.cc

build/pdbresiduemodificationcard.o: src/InputSet/PdbFileSpace/pdbresiduemodificationcard.cc includes/InputSet/PdbFileSpace/pdbresiduemodificationcard.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbresiduemodificationcard.o src/InputSet/PdbFileSpace/pdbresiduemodificationcard.cc

build/pdbresiduemodificationsection.o: src/InputSet/PdbFileSpace/pdbresiduemodificationsection.cc includes/InputSet/PdbFileSpace/pdbresiduemodificationcard.hpp \
		includes/InputSet/PdbFileSpace/pdbresiduemodificationsection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbresiduemodificationsection.o src/InputSet/PdbFileSpace/pdbresiduemodificationsection.cc

build/pdbresiduesequencecard.o: src/InputSet/PdbFileSpace/pdbresiduesequencecard.cc includes/InputSet/PdbFileSpace/pdbresiduesequencecard.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbresiduesequencecard.o src/InputSet/PdbFileSpace/pdbresiduesequencecard.cc

build/pdbresiduesequencesection.o: src/InputSet/PdbFileSpace/pdbresiduesequencesection.cc includes/InputSet/PdbFileSpace/pdbresiduesequencecard.hpp \
		includes/InputSet/PdbFileSpace/pdbresiduesequencesection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbresiduesequencesection.o src/InputSet/PdbFileSpace/pdbresiduesequencesection.cc

build/pdbrevisiondatacard.o: src/InputSet/PdbFileSpace/pdbrevisiondatacard.cc includes/InputSet/PdbFileSpace/pdbrevisiondatacard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbrevisiondatacard.o src/InputSet/PdbFileSpace/pdbrevisiondatacard.cc

build/pdbrevisiondatasection.o: src/InputSet/PdbFileSpace/pdbrevisiondatasection.cc includes/InputSet/PdbFileSpace/pdbrevisiondatasection.hpp \
		includes/InputSet/PdbFileSpace/pdbrevisiondatacard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbrevisiondatasection.o src/InputSet/PdbFileSpace/pdbrevisiondatasection.cc

build/pdbscalencard.o: src/InputSet/PdbFileSpace/pdbscalencard.cc includes/InputSet/PdbFileSpace/pdbscalencard.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbscalencard.o src/InputSet/PdbFileSpace/pdbscalencard.cc

build/pdbscalensection.o: src/InputSet/PdbFileSpace/pdbscalensection.cc includes/InputSet/PdbFileSpace/pdbscalensection.hpp \
		includes/InputSet/PdbFileSpace/pdbscalencard.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbscalensection.o src/InputSet/PdbFileSpace/pdbscalensection.cc

build/pdbsequenceadvancedcard.o: src/InputSet/PdbFileSpace/pdbsequenceadvancedcard.cc includes/InputSet/PdbFileSpace/pdbsequenceadvancedcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbsequenceadvancedcard.o src/InputSet/PdbFileSpace/pdbsequenceadvancedcard.cc

build/pdbsequenceadvancedsection.o: src/InputSet/PdbFileSpace/pdbsequenceadvancedsection.cc includes/InputSet/PdbFileSpace/pdbsequenceadvancedsection.hpp \
		includes/InputSet/PdbFileSpace/pdbsequenceadvancedcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbsequenceadvancedsection.o src/InputSet/PdbFileSpace/pdbsequenceadvancedsection.cc

build/pdbsheetcard.o: src/InputSet/PdbFileSpace/pdbsheetcard.cc includes/InputSet/PdbFileSpace/pdbsheetcard.hpp \
		includes/InputSet/PdbFileSpace/pdbsheetstrand.hpp \
		includes/InputSet/PdbFileSpace/pdbsheetstrandresidue.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbsheetcard.o src/InputSet/PdbFileSpace/pdbsheetcard.cc

build/pdbsheetsection.o: src/InputSet/PdbFileSpace/pdbsheetsection.cc includes/InputSet/PdbFileSpace/pdbsheetcard.hpp \
		includes/InputSet/PdbFileSpace/pdbsheetsection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbsheetsection.o src/InputSet/PdbFileSpace/pdbsheetsection.cc

build/pdbsheetstrand.o: src/InputSet/PdbFileSpace/pdbsheetstrand.cc includes/InputSet/PdbFileSpace/pdbsheetstrand.hpp \
		includes/InputSet/PdbFileSpace/pdbsheetstrandresidue.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbsheetstrand.o src/InputSet/PdbFileSpace/pdbsheetstrand.cc

build/pdbsheetstrandresidue.o: src/InputSet/PdbFileSpace/pdbsheetstrandresidue.cc includes/InputSet/PdbFileSpace/pdbsheetstrandresidue.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbsheetstrandresidue.o src/InputSet/PdbFileSpace/pdbsheetstrandresidue.cc

build/pdbsitecard.o: src/InputSet/PdbFileSpace/pdbsitecard.cc includes/InputSet/PdbFileSpace/pdbsitecard.hpp \
		includes/InputSet/PdbFileSpace/pdbsiteresidue.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbsitecard.o src/InputSet/PdbFileSpace/pdbsitecard.cc

build/pdbsiteresidue.o: src/InputSet/PdbFileSpace/pdbsiteresidue.cc includes/InputSet/PdbFileSpace/pdbsiteresidue.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbsiteresidue.o src/InputSet/PdbFileSpace/pdbsiteresidue.cc

build/pdbsitesection.o: src/InputSet/PdbFileSpace/pdbsitesection.cc includes/InputSet/PdbFileSpace/pdbsitesection.hpp \
		includes/InputSet/PdbFileSpace/pdbsitecard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbsitesection.o src/InputSet/PdbFileSpace/pdbsitesection.cc

build/pdbsourcecard.o: src/InputSet/PdbFileSpace/pdbsourcecard.cc includes/InputSet/PdbFileSpace/pdbsourcecard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbsourcecard.o src/InputSet/PdbFileSpace/pdbsourcecard.cc

build/pdbsourcesection.o: src/InputSet/PdbFileSpace/pdbsourcesection.cc includes/InputSet/PdbFileSpace/pdbsourcesection.hpp \
		includes/InputSet/PdbFileSpace/pdbsourcecard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbsourcesection.o src/InputSet/PdbFileSpace/pdbsourcesection.cc

build/pdbsplitsection.o: src/InputSet/PdbFileSpace/pdbsplitsection.cc includes/InputSet/PdbFileSpace/pdbsplitsection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbsplitsection.o src/InputSet/PdbFileSpace/pdbsplitsection.cc

build/pdbsupersededentriescard.o: src/InputSet/PdbFileSpace/pdbsupersededentriescard.cc includes/InputSet/PdbFileSpace/pdbsupersededentriescard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbsupersededentriescard.o src/InputSet/PdbFileSpace/pdbsupersededentriescard.cc

build/pdbsupersededentriessection.o: src/InputSet/PdbFileSpace/pdbsupersededentriessection.cc includes/InputSet/PdbFileSpace/pdbsupersededentriessection.hpp \
		includes/InputSet/PdbFileSpace/pdbsupersededentriescard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbsupersededentriessection.o src/InputSet/PdbFileSpace/pdbsupersededentriessection.cc

build/pdbtitlesection.o: src/InputSet/PdbFileSpace/pdbtitlesection.cc includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbtitlesection.o src/InputSet/PdbFileSpace/pdbtitlesection.cc

build/pdbqtatom.o: src/InputSet/PdbqtFileSpace/pdbqtatom.cc includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbqtatom.o src/InputSet/PdbqtFileSpace/pdbqtatom.cc

build/pdbqtatomcard.o: src/InputSet/PdbqtFileSpace/pdbqtatomcard.cc includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbqtatomcard.o src/InputSet/PdbqtFileSpace/pdbqtatomcard.cc

build/pdbqtbranchcard.o: src/InputSet/PdbqtFileSpace/pdbqtbranchcard.cc includes/InputSet/PdbqtFileSpace/pdbqtbranchcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbqtbranchcard.o src/InputSet/PdbqtFileSpace/pdbqtbranchcard.cc

build/pdbqtcompoundcard.o: src/InputSet/PdbqtFileSpace/pdbqtcompoundcard.cc includes/InputSet/PdbqtFileSpace/pdbqtcompoundcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbqtcompoundcard.o src/InputSet/PdbqtFileSpace/pdbqtcompoundcard.cc

build/pdbqtfile.o: src/InputSet/PdbqtFileSpace/pdbqtfile.cc includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtcompoundcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtremarkcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtrootcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtbranchcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqttorsionaldofcard.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbqtfile.o src/InputSet/PdbqtFileSpace/pdbqtfile.cc

build/pdbqtfileprocessingexception.o: src/InputSet/PdbqtFileSpace/pdbqtfileprocessingexception.cc includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbqtfileprocessingexception.o src/InputSet/PdbqtFileSpace/pdbqtfileprocessingexception.cc

build/pdbqtmodel.o: src/InputSet/PdbqtFileSpace/pdbqtmodel.cc includes/InputSet/PdbqtFileSpace/pdbqttorsionaldofcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtcompoundcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtremarkcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbqtmodel.o src/InputSet/PdbqtFileSpace/pdbqtmodel.cc

build/pdbqtmodelcard.o: src/InputSet/PdbqtFileSpace/pdbqtmodelcard.cc includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbqtmodelcard.o src/InputSet/PdbqtFileSpace/pdbqtmodelcard.cc

build/pdbqtmodelresidueset.o: src/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.cc includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtbranchcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtrootcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbqtmodelresidueset.o src/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.cc

build/pdbqtremarkcard.o: src/InputSet/PdbqtFileSpace/pdbqtremarkcard.cc includes/InputSet/PdbqtFileSpace/pdbqtremarkcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbqtremarkcard.o src/InputSet/PdbqtFileSpace/pdbqtremarkcard.cc

build/pdbqtrootcard.o: src/InputSet/PdbqtFileSpace/pdbqtrootcard.cc includes/InputSet/PdbqtFileSpace/pdbqtrootcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbqtrootcard.o src/InputSet/PdbqtFileSpace/pdbqtrootcard.cc

build/pdbqttorsionaldofcard.o: src/InputSet/PdbqtFileSpace/pdbqttorsionaldofcard.cc includes/InputSet/PdbqtFileSpace/pdbqttorsionaldofcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbqttorsionaldofcard.o src/InputSet/PdbqtFileSpace/pdbqttorsionaldofcard.cc

build/topologyangle.o: src/InputSet/TopologyFileSpace/topologyangle.cc includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/topologyangle.o src/InputSet/TopologyFileSpace/topologyangle.cc

build/topologyangletype.o: src/InputSet/TopologyFileSpace/topologyangletype.cc includes/InputSet/TopologyFileSpace/topologyangletype.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/topologyangletype.o src/InputSet/TopologyFileSpace/topologyangletype.cc

build/topologyassembly.o: src/InputSet/TopologyFileSpace/topologyassembly.cc includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/topologyassembly.o src/InputSet/TopologyFileSpace/topologyassembly.cc

build/topologyatom.o: src/InputSet/TopologyFileSpace/topologyatom.cc includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/topologyatom.o src/InputSet/TopologyFileSpace/topologyatom.cc

build/topologyatompair.o: src/InputSet/TopologyFileSpace/topologyatompair.cc includes/InputSet/TopologyFileSpace/topologyatompair.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/topologyatompair.o src/InputSet/TopologyFileSpace/topologyatompair.cc

build/topologybond.o: src/InputSet/TopologyFileSpace/topologybond.cc includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/topologybond.o src/InputSet/TopologyFileSpace/topologybond.cc

build/topologybondtype.o: src/InputSet/TopologyFileSpace/topologybondtype.cc includes/InputSet/TopologyFileSpace/topologybondtype.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/topologybondtype.o src/InputSet/TopologyFileSpace/topologybondtype.cc

build/topologydihedral.o: src/InputSet/TopologyFileSpace/topologydihedral.cc includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/topologydihedral.o src/InputSet/TopologyFileSpace/topologydihedral.cc

build/topologydihedraltype.o: src/InputSet/TopologyFileSpace/topologydihedraltype.cc includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/topologydihedraltype.o src/InputSet/TopologyFileSpace/topologydihedraltype.cc

build/topologyfile.o: src/InputSet/TopologyFileSpace/topologyfile.cc includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyfileprocessingexception.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/topologyfile.o src/InputSet/TopologyFileSpace/topologyfile.cc

build/topologyfileprocessingexception.o: src/InputSet/TopologyFileSpace/topologyfileprocessingexception.cc includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp \
		includes/InputSet/TopologyFileSpace/topologyfileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/topologyfileprocessingexception.o src/InputSet/TopologyFileSpace/topologyfileprocessingexception.cc

build/topologyresidue.o: src/InputSet/TopologyFileSpace/topologyresidue.cc includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/topologyresidue.o src/InputSet/TopologyFileSpace/topologyresidue.cc

build/response.o: src/InputSet/Utilities/response.cc includes/InputSet/Utilities/response.hpp \
		includes/Glycan/note.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/response.o src/InputSet/Utilities/response.cc

build/amberatomtypeinfo.o: src/MolecularMetadata/GLYCAM/amberatomtypeinfo.cpp includes/MolecularMetadata/GLYCAM/amberatomtypeinfo.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/amberatomtypeinfo.o src/MolecularMetadata/GLYCAM/amberatomtypeinfo.cpp

build/amberelements.o: src/MolecularMetadata/GLYCAM/amberelements.cpp includes/MolecularMetadata/AMBER/amberelements.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/amberelements.o src/MolecularMetadata/GLYCAM/amberelements.cpp

build/bondlengthbytypepair.o: src/MolecularMetadata/GLYCAM/bondlengthbytypepair.cpp includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/bondlengthbytypepair.o src/MolecularMetadata/GLYCAM/bondlengthbytypepair.cpp

build/dihedralangledata.o: src/MolecularMetadata/GLYCAM/dihedralangledata.cc includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/dihedralangledata.o src/MolecularMetadata/GLYCAM/dihedralangledata.cc

build/glycam06residuecodes.o: src/MolecularMetadata/GLYCAM/glycam06residuecodes.cc includes/MolecularMetadata/GLYCAM/glycam06residuecodes.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/glycam06residuecodes.o src/MolecularMetadata/GLYCAM/glycam06residuecodes.cc

build/glycam06residueinfo.o: src/MolecularMetadata/GLYCAM/glycam06residueinfo.cpp includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/glycam06residueinfo.o src/MolecularMetadata/GLYCAM/glycam06residueinfo.cpp

build/assembly.o: src/MolecularModeling/Assembly/assembly.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/OffFileSpace/offfile.hpp \
		includes/MolecularMetadata/elementattributes.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/assembly.o src/MolecularModeling/Assembly/assembly.cc

build/assemblybuilder.o: src/MolecularModeling/Assembly/assemblybuilder.cc includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/GeometryTopology/cell.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/GeometryTopology/rotation.hpp \
		includes/GeometryTopology/ResidueLinkages/residue_linkage.hpp \
		includes/MolecularModeling/Selections/selections.hpp \
		includes/GeometryTopology/ResidueLinkages/rotatable_dihedral.hpp \
		includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp \
		includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefileprocessingexception.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdbmastercard.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyfileprocessingexception.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/MolecularModeling/overlaps.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/assemblybuilder.o src/MolecularModeling/Assembly/assemblybuilder.cc

build/geometryoperation.o: src/MolecularModeling/Assembly/geometryoperation.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/MolecularModeling/overlaps.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp \
		includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp \
		includes/MolecularMetadata/GLYCAM/amberatomtypeinfo.hpp \
		includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/geometryoperation.o src/MolecularModeling/Assembly/geometryoperation.cc

build/glycamnaming.o: src/MolecularModeling/Assembly/glycamnaming.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/glycamnaming.o src/MolecularModeling/Assembly/glycamnaming.cc

build/ionization.o: src/MolecularModeling/Assembly/ionization.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/ionization.o src/MolecularModeling/Assembly/ionization.cc

build/manipulationoperation.o: src/MolecularModeling/Assembly/manipulationoperation.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/manipulationoperation.o src/MolecularModeling/Assembly/manipulationoperation.cc

build/selection.o: src/MolecularModeling/Assembly/selection.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/selection.o src/MolecularModeling/Assembly/selection.cc

build/solvation.o: src/MolecularModeling/Assembly/solvation.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/solvation.o src/MolecularModeling/Assembly/solvation.cc

build/structurebuilder.o: src/MolecularModeling/Assembly/structurebuilder.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/structurebuilder.o src/MolecularModeling/Assembly/structurebuilder.cc

build/edge.o: src/MolecularModeling/Graph/edge.cc includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/Graph/node.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/edge.o src/MolecularModeling/Graph/edge.cc

build/graph.o: src/MolecularModeling/Graph/graph.cc includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/graph.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/graph.o src/MolecularModeling/Graph/graph.cc

build/node.o: src/MolecularModeling/Graph/node.cc includes/MolecularModeling/Graph/node.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/node.o src/MolecularModeling/Graph/node.cc

build/residue.o: src/MolecularModeling/Residue/residue.cc includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/MolecularModeling/overlaps.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/residue.o src/MolecularModeling/Residue/residue.cc

build/residue3dtemplatelibrary.o: src/MolecularModeling/Residue/residue3dtemplatelibrary.cc 
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/residue3dtemplatelibrary.o src/MolecularModeling/Residue/residue3dtemplatelibrary.cc

build/residuebuilder.o: src/MolecularModeling/Residue/residuebuilder.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/GeometryTopology/geometrytopology.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp \
		includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp \
		includes/MolecularMetadata/GLYCAM/amberatomtypeinfo.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/residuebuilder.o src/MolecularModeling/Residue/residuebuilder.cc

build/residuenode.o: src/MolecularModeling/Residue/residuenode.cc includes/MolecularModeling/atom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/MolecularModeling/residuenode.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/residuenode.o src/MolecularModeling/Residue/residuenode.cc

build/residueproperties.o: src/MolecularModeling/Residue/residueproperties.cc includes/MolecularModeling/residueproperties.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/residueproperties.o src/MolecularModeling/Residue/residueproperties.cc

build/selections.o: src/MolecularModeling/Selections/selections.cc includes/MolecularModeling/Selections/selections.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/generictypedefs.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/common.hpp \
		includes/Glycan/sugarname.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/selections.o src/MolecularModeling/Selections/selections.cc

build/libraryfile.o: src/ParameterSet/LibraryFileSpace/libraryfile.cc includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/libraryfile.o src/ParameterSet/LibraryFileSpace/libraryfile.cc

build/libraryfileatom.o: src/ParameterSet/LibraryFileSpace/libraryfileatom.cc includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/GeometryTopology/coordinate.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/libraryfileatom.o src/ParameterSet/LibraryFileSpace/libraryfileatom.cc

build/libraryfileprocessingexception.o: src/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.cc includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/libraryfileprocessingexception.o src/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.cc

build/libraryfileresidue.o: src/ParameterSet/LibraryFileSpace/libraryfileresidue.cc includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/libraryfileresidue.o src/ParameterSet/LibraryFileSpace/libraryfileresidue.cc

build/offfile.o: src/ParameterSet/OffFileSpace/offfile.cc includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/ParameterSet/OffFileSpace/offfile.hpp \
		includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/MolecularMetadata/elementattributes.hpp \
		includes/ParameterSet/OffFileSpace/offfileresidue.hpp \
		includes/ParameterSet/OffFileSpace/offfileatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/offfile.o src/ParameterSet/OffFileSpace/offfile.cc

build/offfileatom.o: src/ParameterSet/OffFileSpace/offfileatom.cc includes/ParameterSet/OffFileSpace/offfileatom.hpp \
		includes/GeometryTopology/coordinate.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/offfileatom.o src/ParameterSet/OffFileSpace/offfileatom.cc

build/offfileresidue.o: src/ParameterSet/OffFileSpace/offfileresidue.cc includes/ParameterSet/OffFileSpace/offfileresidue.hpp \
		includes/ParameterSet/OffFileSpace/offfileatom.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/offfileresidue.o src/ParameterSet/OffFileSpace/offfileresidue.cc

build/parameterfile.o: src/ParameterSet/ParameterFileSpace/parameterfile.cc includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/parameterfile.o src/ParameterSet/ParameterFileSpace/parameterfile.cc

build/parameterfileangle.o: src/ParameterSet/ParameterFileSpace/parameterfileangle.cc includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/parameterfileangle.o src/ParameterSet/ParameterFileSpace/parameterfileangle.cc

build/parameterfileatom.o: src/ParameterSet/ParameterFileSpace/parameterfileatom.cc includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/parameterfileatom.o src/ParameterSet/ParameterFileSpace/parameterfileatom.cc

build/parameterfilebond.o: src/ParameterSet/ParameterFileSpace/parameterfilebond.cc includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/parameterfilebond.o src/ParameterSet/ParameterFileSpace/parameterfilebond.cc

build/parameterfiledihedral.o: src/ParameterSet/ParameterFileSpace/parameterfiledihedral.cc includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/parameterfiledihedral.o src/ParameterSet/ParameterFileSpace/parameterfiledihedral.cc

build/parameterfiledihedralterm.o: src/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.cc includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/parameterfiledihedralterm.o src/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.cc

build/parameterfileprocessingexception.o: src/ParameterSet/ParameterFileSpace/parameterfileprocessingexception.cc includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/parameterfileprocessingexception.o src/ParameterSet/ParameterFileSpace/parameterfileprocessingexception.cc

build/prepfile.o: src/ParameterSet/PrepFileSpace/prepfile.cc includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/prepfile.o src/ParameterSet/PrepFileSpace/prepfile.cc

build/prepfileatom.o: src/ParameterSet/PrepFileSpace/prepfileatom.cc includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/prepfileatom.o src/ParameterSet/PrepFileSpace/prepfileatom.cc

build/prepfileprocessingexception.o: src/ParameterSet/PrepFileSpace/prepfileprocessingexception.cc includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/prepfileprocessingexception.o src/ParameterSet/PrepFileSpace/prepfileprocessingexception.cc

build/prepfileresidue.o: src/ParameterSet/PrepFileSpace/prepfileresidue.cc includes/utils.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/prepfileresidue.o src/ParameterSet/PrepFileSpace/prepfileresidue.cc

build/pdbpreprocessor.o: src/Resolver/PdbPreprocessor/pdbpreprocessor.cc includes/Resolver/PdbPreprocessor/pdbpreprocessor.hpp \
		includes/InputSet/PdbFileSpace/pdbresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessorresidueinfo.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbpreprocessor.o src/Resolver/PdbPreprocessor/pdbpreprocessor.cc

build/pdbpreprocessoralternateresidue.o: src/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.cc includes/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbpreprocessoralternateresidue.o src/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.cc

build/pdbpreprocessorchaintermination.o: src/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.cc includes/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbpreprocessorchaintermination.o src/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.cc

build/pdbpreprocessordisulfidebond.o: src/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.cc includes/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbpreprocessordisulfidebond.o src/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.cc

build/pdbpreprocessorhistidinemapping.o: src/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.cc includes/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbpreprocessorhistidinemapping.o src/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.cc

build/pdbpreprocessormissingresidue.o: src/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.cc includes/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbpreprocessormissingresidue.o src/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.cc

build/pdbpreprocessorreplacedhydrogen.o: src/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.cc includes/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbpreprocessorreplacedhydrogen.o src/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.cc

build/pdbpreprocessorresidueinfo.o: src/Resolver/PdbPreprocessor/pdbpreprocessorresidueinfo.cc includes/Resolver/PdbPreprocessor/pdbpreprocessorresidueinfo.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbpreprocessorresidueinfo.o src/Resolver/PdbPreprocessor/pdbpreprocessorresidueinfo.cc

build/pdbpreprocessorunrecognizedheavyatom.o: src/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.cc includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbpreprocessorunrecognizedheavyatom.o src/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.cc

build/pdbpreprocessorunrecognizedresidue.o: src/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.cc includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbpreprocessorunrecognizedresidue.o src/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.cc

build/analysis.o: src/MolecularModeling/Assembly/Ontology/analysis.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/analysis.o src/MolecularModeling/Assembly/Ontology/analysis.cc

build/gf_query.o: src/MolecularModeling/Assembly/Ontology/gf_query.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/gf_query.o src/MolecularModeling/Assembly/Ontology/gf_query.cc

build/population.o: src/MolecularModeling/Assembly/Ontology/population.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/population.o src/MolecularModeling/Assembly/Ontology/population.cc

build/query.o: src/MolecularModeling/Assembly/Ontology/query.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/query.o src/MolecularModeling/Assembly/Ontology/query.cc

build/coordinatebuilder.o: src/MolecularModeling/Assembly/Output/coordinatebuilder.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/coordinatebuilder.o src/MolecularModeling/Assembly/Output/coordinatebuilder.cc

build/librarybuilder.o: src/MolecularModeling/Assembly/Output/librarybuilder.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/librarybuilder.o src/MolecularModeling/Assembly/Output/librarybuilder.cc

build/pdbbuilder.o: src/MolecularModeling/Assembly/Output/pdbbuilder.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbbuilder.o src/MolecularModeling/Assembly/Output/pdbbuilder.cc

build/pdbqtbuilder.o: src/MolecularModeling/Assembly/Output/pdbqtbuilder.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/pdbqtbuilder.o src/MolecularModeling/Assembly/Output/pdbqtbuilder.cc

build/prepbuilder.o: src/MolecularModeling/Assembly/Output/prepbuilder.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/prepbuilder.o src/MolecularModeling/Assembly/Output/prepbuilder.cc

build/topologybuilder.o: src/MolecularModeling/Assembly/Output/topologybuilder.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/topologybuilder.o src/MolecularModeling/Assembly/Output/topologybuilder.cc

build/oligosaccharidedetection.o: src/MolecularModeling/Assembly/SugarIdentification/oligosaccharidedetection.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp \
		includes/MolecularModeling/ring_shape_detection.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/oligosaccharidedetection.o src/MolecularModeling/Assembly/SugarIdentification/oligosaccharidedetection.cc

build/ringperception.o: src/MolecularModeling/Assembly/SugarIdentification/ringperception.cc includes/MolecularModeling/assembly.hpp \
		includes/GeometryTopology/coordinate.hpp \
		includes/GeometryTopology/plane.hpp \
		includes/common.hpp \
		includes/generictypedefs.hpp \
		includes/Glycan/sugarname.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/oligosaccharidedetectionatom.hpp \
		includes/Glycan/chemicalcode.hpp \
		includes/Glycan/monosaccharide.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/Glycan/note.hpp \
		includes/Glycan/ontologyvocabulary.hpp \
		includes/InputSet/PdbFileSpace/pdbfile.hpp \
		includes/InputSet/PdbFileSpace/inputfile.hpp \
		includes/utils.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp \
		includes/InputSet/TopologyFileSpace/topologyfile.hpp \
		includes/InputSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/InputSet/PdbFileSpace/pdbremarksection.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelsection.hpp \
		includes/InputSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/Glycan/oligosaccharide.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/residueproperties.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/Glycan/glycosidiclinkage.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp \
		includes/InputSet/Utilities/response.hpp \
		includes/MolecularModeling/molecule.hpp \
		includes/MolecularModeling/Graph/graph.hpp \
		includes/MolecularModeling/Graph/node.hpp \
		includes/MolecularModeling/Graph/edge.hpp \
		includes/MolecularModeling/residuenode.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp \
		includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp \
		includes/InputSet/TopologyFileSpace/topologyassembly.hpp \
		includes/InputSet/TopologyFileSpace/topologyresidue.hpp \
		includes/InputSet/TopologyFileSpace/topologyatom.hpp \
		includes/InputSet/TopologyFileSpace/topologybond.hpp \
		includes/InputSet/TopologyFileSpace/topologybondtype.hpp \
		includes/InputSet/TopologyFileSpace/topologyangle.hpp \
		includes/InputSet/TopologyFileSpace/topologyangletype.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedral.hpp \
		includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/InputSet/TopologyFileSpace/topologyatompair.hpp \
		includes/InputSet/PdbFileSpace/pdbtitlesection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomcard.hpp \
		includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbatomsection.hpp \
		includes/InputSet/PdbFileSpace/pdbconnectsection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcard.hpp \
		includes/InputSet/PdbFileSpace/pdblinksection.hpp \
		includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp \
		includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp \
		includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/GeometryTopology/grid.hpp \
		includes/GeometryTopology/cell.hpp
	@echo "[GMML] CXX $<"
	$(VB)$(CXX) -c $(CXXFLAGS) $(INCPATH) -o build/ringperception.o src/MolecularModeling/Assembly/SugarIdentification/ringperception.cc

