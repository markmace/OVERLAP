//SPHALERON TRANSITION RATE SIMULATION
#define METRIC_FLAG MINKOWSKI_FLAG
#define FERMION_FLAG OVERLAP_FLAG

#include <iostream>
#include <string>

/////////////////////////////////////////
// SU(Nc) GROUP AND ALGEBRA OPERATIONS //
/////////////////////////////////////////

//INCLUDE SU(Nc) ALGEBRA AND GROUP DEFINITION
#include "SUNc/Definitions.cpp"

/////////////////////////////
// RANDOM NUMBER GENERATOR //
/////////////////////////////

//INCLUDE RANDOM NUMBER GENERATOR
#include "MISC/RNG/GSLRandomNumberGenerator.cpp"

///////////////////////////
//  OUTPUT HANDLING      //
///////////////////////////

#include "IO/StringManipulation.cpp"
#include "IO/OutputManagement.cpp"

///////////////////////////
//   INPUT HANDLING      //
///////////////////////////
#include "IO/InputManagement.cpp"

/////////////////////////////////////////////
// DEFINITION OF LATTICE GRID AND INDEXING //
/////////////////////////////////////////////

#include "LATTICE/3DGrid.cpp"
#include "LATTICE/GaugeLinks.cpp"
#include "LATTICE/ElectricFields.cpp"
#include "LATTICE/GaugeTransformations.cpp"
#include "LATTICE/GLinksEFields.cpp"
#include "LATTICE/Indexing.cpp"
#include "LATTICE/Momenta.cpp"

///////////////////////////
//GENERAL PURPOSE MACROS //
///////////////////////////

//INCLUDE MACROS TO COMPUTE PLAQUETTES

#include "CALC/Plaquettes.cpp"
#include "CALC/AvgFieldStrength.cpp"

//INCLUDE MACROS TO DETERMINE AVERAGE FIELD STRENGTH
#include "CALC/AvgFieldStrength.cpp"

//INCLUDE MACROS TO COMPUTE GAUSS LAW VIOLATION
#include "CALC/GaussViolation.cpp"

//INCLUDE MACROS TO PERFORM GAUGE TRANSFORMATIONS
#include "CALC/GaugeTransformation.cpp"

//INCLUDE UNIMPROVED AND IMPROVED CHERN-SIMONS OPERATORS
#include "CALC/ChernSimonsOperators.cpp"

//INCLUDE ROUTINES TO MEASURE CHERN SIMONS NUMBER DERIVATIVE
#include "TOPOLOGY/ChernSimonsDerivative.cpp"

////////////////////////
// INITIAL CONDITIONS //
////////////////////////

DOUBLE Qs=1.0;
DOUBLE n0=1.0;

///////////////////////
//REAL-TIME DYNAMICS //
///////////////////////

//INCLUDE DYNAMICS
#include "DYNAMICS/Dynamics.cpp"

//////////////////////
//BASIC OBSERVABLES //
//////////////////////

//INCLUDE MEASUREMENT OF BULK OBSERVABLES
#include "OBSERVABLES/BulkObservables.cpp"

//INCLUDE MEASUREMENT OF GAUSS LAW VIOLATION
#include "OBSERVABLES/GaussLawViolation.cpp"

//INCLUDE MEASUREMENT OF UNITARITY VIOLATION
#include "OBSERVABLES/UnitarityViolation.cpp"

//INCLUDE MEASUREMENT OF ENERGY MOMENTUM TENSOR
#include "OBSERVABLES/EnergyMomentumTensor.cpp"

//INCLUDE MEASUREMENT OF HARD SCALES
#include "OBSERVABLES/HardScales.cpp"

//INCLUDE HISTOGRAM
#include "MISC/HISTOGRAM/Histogram.cpp"
#include "MISC/HISTOGRAM/MultiHistogram.cpp"

////////////////////////////////
//GAUSS LAW FIXING PROCEDURES //
////////////////////////////////
#include "GAUSSLAWFIXING/GaussRestore.cpp"

///////////////////////
//INITIAL CONDITIONS //
///////////////////////

#include "INITIALCONDITIONS/SetZero.cpp"
#include "INITIALCONDITIONS/SetRandomMatrices.cpp"
#include "INITIALCONDITIONS/SetQuasiParticles.cpp"

//////////////////////
// THERMAL DYNAMICS //
//////////////////////

#include "DYNAMICS/ThermalDynamics.cpp"

////////////////////////////
//GAUGE FIXING PROCEDURES //
////////////////////////////

//INCLUDE GAUGE TRANSFORMATION RULES AND BUFFERS
#include "GAUGETRANSFORMATION/GaugeTransformation.cpp"

//INCLUDE COULOMB GAUGE FIXING PROCEDURE
#include "GAUGETRANSFORMATION/CoulombGaugeFixing.cpp"

//INCLUDE MEASUREMENT OF SPECTRA
#include "OBSERVABLES/Spectra.cpp"

////////////////////////////////////////
//MEASUREMENTS OF CHERN-SIMONS NUMBER //
////////////////////////////////////////

//INCLUDE VACUUM TOPOLOGY ESTIMATOR
#include "TOPOLOGY/VacuumEstimator.cpp"

/////////////////////////////
//YANG-MILLS GRADIENT FLOW //
/////////////////////////////

//GRADIENT FLOW COOLING
#include "GRADIENTFLOW/GradientFlow.cpp"
#include "GRADIENTFLOW/Cooling.cpp"

////////////////////////////////////////
//MEASUREMENTS OF CHERN-SIMONS NUMBER //
////////////////////////////////////////

#include "TOPOLOGY/CoolingMethod.cpp"
#include "TOPOLOGY/WindingNumber.cpp"
#include "TOPOLOGY/IntegralEstimateMethod.cpp"


////////////////////////
// HANDMADE SPHALERON //
////////////////////////

// SPHALERON SIZE AND TRANSITION TIME IN LATTICE UNITS //
DOUBLE rSphaleron=6.0; DOUBLE tSphaleron=6.0;

#include "INITIALCONDITIONS/HandmadeSphaleron.cpp"

//////////////////////////
//INPUT/OUTPUT ROUTINES //
//////////////////////////

#include "IO/SaveConfiguration.cpp"
#include "IO/LoadConfiguration.cpp"

//////////
//  QED //
//////////

// INCLUDE QED //
#include "QED/QED.cpp"

/////////////////////
// WILSON FERMIONS //
/////////////////////

// DIRAC ALGEBRA //
#include "DiracAlgebra.cpp"

// FERMION MASS AND WILSON TERM //
DOUBLE mFermion=0.0; DOUBLE rWilson=1.0;

// U(1) MAGNETIC FIELD FLUX
INT BFlux=1.0;

// DOMAIN WALL HEIGHT //
DOUBLE MDWHeight=1.4;

// BASIC FERMION FIELD OBJECT //
#include "FERMIONS/FermionField.cpp"

// EVOLUTION
#if FERMION_FLAG==WILSON_FLAG
#include "FERMIONS/WILSON/EvolutionOperator.cpp"
#include "FERMIONS/WILSON/Update.cpp"
#include "FERMIONS/WILSON/Currents.cpp"
#endif

#if FERMION_FLAG==OVERLAP_FLAG
#include "FERMIONS/OVERLAP/OverlapEvolutionOperator.cpp"
#include "FERMIONS/OVERLAP/OverlapUpdate.cpp"
#include "FERMIONS/OVERLAP/OverlapCurrents.cpp"
#endif

// INITIAL CONDITIONS //
#include "INITIALCONDITIONS/FERMIONS/SetFreeFermions.cpp"

// INCLUDE DIAGONALIZATION ROUTINE //
//#include "MISC/Diagonalization.cpp"
#include "MISC/LapackDiagonalization.cpp"


INT INPUT_QED_EIGENFUNCTIONS=0; INT SAVE_QED_EIGENFUNCTIONS=1;

// MAGNETIC FIELDS //
//#include "INITIALCONDITIONS/FERMIONS/SetQEDEigenfunctions.cpp"

#include "INITIALCONDITIONS/FERMIONS/SetOverlapQEDEigenfunctions.cpp"



///////////////////////////
// SIMULATION PROCEDURE  //
///////////////////////////


namespace Simulation {
    
    // MPI RANDOM NUMBER SEED //
    INT MY_MPI_RNG_SEED;
    
    /////////////////////////////
    //INITIALIZATION AND EXIT  //
    /////////////////////////////
    
    INT CoolingFrequency,CalibFrequency;
    
    void Init(){
        
        ///////////////////////////////////
        //SET NCS MEASUREMENT PARAMETERS //
        ///////////////////////////////////
        
        // SET COOLING FREQUENCY //
        CoolingFrequency=5; ChernSimonsNumber::CoolingMethod::StandardCoolingMaxSteps=0;
        CalibFrequency=30*CoolingFrequency; ChernSimonsNumber::CoolingMethod::CalibrationCoolingMaxSteps=0;
        
        // SET BLOCKING PARAMETERS //
        Cooling::BlockFrequency=32;
        
        ChernSimonsNumber::CoolingMethod::StandardBlockingLevel=0;
        ChernSimonsNumber::CoolingMethod::CalibrationBlockingLevel=2;
        
        // CONSISTENCY CHECK FOR BLOCKING //
        if((ChernSimonsNumber::CoolingMethod::StandardBlockingLevel*Cooling::BlockFrequency*ChernSimonsNumber::CoolingMethod::StandardBlockingLevel>ChernSimonsNumber::CoolingMethod::StandardCoolingMaxSteps) || (ChernSimonsNumber::CoolingMethod::StandardBlockingLevel*Cooling::BlockFrequency*ChernSimonsNumber::CoolingMethod::CalibrationBlockingLevel>ChernSimonsNumber::CoolingMethod::CalibrationCoolingMaxSteps)){
            
            std::cerr << "#ERROR -- INCONSISTENT COOLING PARAMETERS" << std::endl;
            exit(0);
            
        }
        
        
        /////////////////////
        // INITIALIZATIONS //
        /////////////////////
        
        // INITIALIZE 3D LATTICE SETUP //
        Lattice::Init();
        
        // INITIALIZE GAUGE TRANSFORMATIONS //
        GaugeTransformation::Init();
        
        // INITIALIZE GAUGE TRANSFORMED FIELDS //
        GaugeFixedVariables::Init();
        
        // INITIALIZE FOURIER SPACE VARIABLES //
        FourierSpace::Init();
        
        // INITIALIZE COOLING METHOD //
        ChernSimonsNumber::CoolingMethod::Init();
        
        // INITIALIZE WINDING NUMBER MEASUREMENT //
        WindingNumber::Init();
        
        // INITIALIZE QED //
        QED::Init();
        
        // INITIALIZE FERMIONS //
        Fermions::Init();
        
        // INITIALIZE EVOLUTION //
        EvolutionOperator::Init();
        
    }
    
    /////////////////////
    //CREATE INFO FILE //
    /////////////////////
    
    void CreateInfoFile(){
        
        // CREATE INFO FILE //
        std::ofstream OutStream;
        OutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"InfoID",MY_MPI_RNG_SEED,".txt").c_str());
        
        // CREATE INFO FILE CONTAINING PARAMETERS //
        OutStream << "#LATTICE DATA" << std::endl;
        OutStream << "#Nx,Ny,Nz= " << GLinks::U->N[0] << " " << GLinks::U->N[1] << " " << GLinks::U->N[2] << std::endl;
        OutStream << "#ax,ay,az= " << GLinks::U->a[0] << " " << GLinks::U->a[1] << " " << GLinks::U->a[2] << std::endl;
        OutStream << "#dTau= " << Dynamics::dTau << std::endl;
        OutStream << std::endl;
        
        OutStream << "#FERMIONS" << std::endl;
        OutStream << "#mFermion= " << mFermion << std::endl;
        OutStream << "#rWilson= " << rWilson << std::endl;
        OutStream << "#MDWHeight= " << MDWHeight << std::endl;
        
        #if FERMION_FLAG==WILSON_FLAG
        OutStream << "#OperatorImprovement=WILSON" << std::endl;
        #endif
        #if FERMION_FLAG==OVERLAP_FLAG
        OutStream << "#OperatorImprovement=OVERLAP" << std::endl;
        #endif
        OutStream << std::endl;
        
        OutStream << "#SPHALERON" << std::endl;
        OutStream << "#rSphaleron=" << rSphaleron << std::endl;
        OutStream << "#tSphaleron=" << tSphaleron << std::endl;
        OutStream << std::endl;
        
        //OutStream << "#QED" << std::endl;
        //OutStream << "#BFlux= " << BFlux << std::endl;
        
        
        OutStream << "#COOLING PARAMETERS" << std::endl;
        OutStream << "#SqrtDcTauOverSpacing = " << GradientFlow::SqrtDcTauOverSpacing << std::endl;
        OutStream << "#CoolingFrequency = " << CoolingFrequency << std::endl;
        OutStream << "#StandardBlockingLevel= " << ChernSimonsNumber::CoolingMethod::StandardBlockingLevel << std::endl;
        OutStream << "#StandardCoolingMaxSteps= " << ChernSimonsNumber::CoolingMethod::StandardCoolingMaxSteps << std::endl;
        OutStream << "#BlockFrequency= " << Cooling::BlockFrequency << std::endl;
        OutStream << std::endl;
        
        OutStream << "#CALIBRATION PARAMETERS" << std::endl;
        OutStream << "#CalibFrequency= " << CalibFrequency << std::endl;
        OutStream << "#CalibrationCoolingMaxSteps= " << ChernSimonsNumber::CoolingMethod::CalibrationCoolingMaxSteps << std::endl;
        OutStream << "#CalibBlockingLevel= " << ChernSimonsNumber::CoolingMethod::CalibrationBlockingLevel << std::endl;
        OutStream << std::endl;
        
        
        // CLOSE OUPUT STREAM //
        OutStream.close();
        
    }
    
    ////////////////////////////////////
    // CLASSICAL YANG-MILLS EVOLUTION //
    ////////////////////////////////////
    
    void Evolve(DOUBLE MaxTime){
        
        
        ///////////////////////
        // OUTPUT MANAGEMENT //
        ///////////////////////
        
        std::ofstream EnergyOutStream,CoolNCsOutStream,CalibNCsOutStream;
        
        if(MPIBasic::ID==0){
            
            EnergyOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"EnergyID",MY_MPI_RNG_SEED,".txt").c_str());
            CoolNCsOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"CoolNCsID",MY_MPI_RNG_SEED,".txt").c_str());
            CalibNCsOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"CalibNCsID",MY_MPI_RNG_SEED,".txt").c_str());
            
        }
        
        ///////////////////////////////
        // INITIAL STATE OBSERVABLES //
        ///////////////////////////////
        
        if(MPIBasic::ID==0){
            
            // MEASURE BULK OBSERVABLES //
            Observables::Bulk::Update();
            
            // MEASURE HARD SCALES //
            Observables::HardScales::Update();
            
            // CREATE OUPUT //
            EnergyOutStream << Dynamics::Time() << " " << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << " " << Observables::Bulk::ELECTRIC() << " " << Observables::Bulk::MAGNETIC() << " " << Observables::HardScales::LambdaXX() << " " << Observables::HardScales::LambdaYY() << " " << Observables::HardScales::LambdaZZ() << std::endl;
            
        }
        
        /////////////////////////////////////
        // PREPARE GAUGE FIXING ALGORITHMS //
        /////////////////////////////////////
        
        // RESET COULOMB GAUGE FIXING ALGORITHMS //
        CoulombGaugeFixing::Reset();
        
        //////////////////////////////////////
        // INITIALIZE TOPOLOGY MEASUREMENTS //
        //////////////////////////////////////
        
        if(MPIBasic::ID==0){
            
            // RESET CHERN SIMONS-MEASUREMENTS //
            ChernSimonsNumber::DeltaNCsRealTime=DOUBLE(0.0);
            ChernSimonsNumber::DeltaNCsCoolRealTime=DOUBLE(0.0);
            ChernSimonsNumber::DeltaNCsCooling=DOUBLE(0.0);
            
            // START COOLING METHOD //
            ChernSimonsNumber::CoolingMethod::Start();
            
            // CREATE INITIAL COOLING OUTPUT //
            CoolNCsOutStream << Dynamics::Time() << " " << ChernSimonsNumber::DeltaNCsCoolRealTime << " " << ChernSimonsNumber::DeltaNCsCooling << std::endl;
            
            // CALIBRATE COOLING //
            ChernSimonsNumber::CoolingMethod::Calibrate();
            
            // CREATE INITIAL CALIBRATION OUTPUT //
            
            CalibNCsOutStream << Dynamics::Time() << " " << ChernSimonsNumber::DeltaNCsCalibration << " " << ChernSimonsNumber::DeltaNCsPreviousCalibration << " " << ChernSimonsNumber::DeltaNCsRealTimeCalibration << " " << ChernSimonsNumber::DeltaNCsRealTimePreviousCalibration << std::endl;
            
        }
        
        
        //////////////////////
        // MEASURE CURRENTS //
        //////////////////////
        
        if(MPIBasic::ID==0){
            std::cerr << "#COMPUTING CURRENTS AT T=" << Dynamics::Time() << std::endl;
        }
        
        Fermions::Observables::ComputeCurrents(GLinks::U,Fermions::Psi);
        
        // SYNCHRONIZE //
        MPI_Barrier(MPI_COMM_WORLD);
        
        
        ////////////////////
        // TIME EVOLUTION //
        ////////////////////
        
        // COMMANDLINE OUTPUT //
        if(MPIBasic::ID==0){
            std::cerr << "#STARTING TIME EVOLUTION" << std::endl;
            Timing::Reset();
            
        }
        
        // PERFORM EVOLUTION //
        while(Dynamics::Time()<MaxTime){
            
            
            // EVOLVE FERMIONS //
            Dynamics::Fermions::Update(GLinks::U,EFields::E,Fermions::Psi,Fermions::PsiMid);
            
            // COMPUTE FERMIONIC CURRENTS //
            if(Dynamics::tSteps%5==0){
                
                if(MPIBasic::ID==0){
                    std::cerr << "#COMPUTING CURRENTS AT T=" << Dynamics::Time() << " ELAPSED TIME " << Timing::Get() << " s" << std::endl;
                    Timing::Reset();
                }
                
                Fermions::Observables::ComputeCurrents(GLinks::U,Fermions::Psi);
                
                if(MPIBasic::ID==0){
                    std::cerr << "#CURRENTS COMPUTED IN " << Timing::Get() << " s" << std::endl;
                    Timing::Reset();
                }
                
            }
            
            
            // CHECK BULK OBSERVABLES //
            if(Dynamics::tSteps%20==0 && MPIBasic::ID==0){
                
                // MEASURE BULK OBSERVABLES //
                Observables::Bulk::Update();
                
                // MEASURE HARD SCALES //
                Observables::HardScales::Update();
                
                // CREATE OUPUT //
                EnergyOutStream << Dynamics::Time() << " " << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << " " << Observables::Bulk::ELECTRIC() << " " << Observables::Bulk::MAGNETIC() << " " << Observables::HardScales::LambdaXX() << " " << Observables::HardScales::LambdaYY() << " " << Observables::HardScales::LambdaZZ() << std::endl;
            }
            
            // COOLED CHERN SIMONS DERIVATIVE //
            if(Dynamics::tSteps%CoolingFrequency==0 && MPIBasic::ID==0){
                
                // UPDATE COOLING //
                ChernSimonsNumber::CoolingMethod::Update(1);
                
                // MEASURE DIFFERNCE IN CHERN SIMONS NUMBER  //
                CoolNCsOutStream << Dynamics::Time() << " " << ChernSimonsNumber::DeltaNCsCoolRealTime << " " << ChernSimonsNumber::DeltaNCsCooling << std::endl;
            }
            
            // COOLING CALIBRATION //
            if(Dynamics::tSteps%CalibFrequency==0 && MPIBasic::ID==0){
                
                // UPDATE COOLING //
                ChernSimonsNumber::CoolingMethod::Calibrate();
                
                // MEASURE OUTPUT //
                CalibNCsOutStream << Dynamics::Time() << " " << ChernSimonsNumber::DeltaNCsCalibration << " " << ChernSimonsNumber::DeltaNCsPreviousCalibration << " " << ChernSimonsNumber::DeltaNCsRealTimeCalibration << " " << ChernSimonsNumber::DeltaNCsRealTimePreviousCalibration << std::endl;
            }
            
            // PRECISION CHECKS //
            if(Dynamics::tSteps%1000==0){
                
                //CHECK GAUSS LAW VIOLATION
                Observables::GaussLaw::CheckViolation();
                
                //CHECK UNITARITY VIOLATION
                Observables::Unitarity::CheckViolation();
                
            }
            
            // SYNCRHONIZE //
            MPI_Barrier(MPI_COMM_WORLD);
            
        }
        
        // CLOSE OUTPUT STREAMS //
        if(MPIBasic::ID==0){
            EnergyOutStream.close();   CoolNCsOutStream.close();   CalibNCsOutStream.close();
            
        }
    }
    
    //////////////////////////
    //SIMULATION PROCDEDURE //
    //////////////////////////
    
    void Run(INT MPI_RNG_SEED){
        
        ///////////
        // SETUP //
        ///////////
        
        //SET SEED //
        MY_MPI_RNG_SEED=MPI_RNG_SEED;
        
        //INITIALIZE RANDOM NUMBER GENERATOR //
        RandomNumberGenerator::Init(MY_MPI_RNG_SEED);
        
        //INITIALIZE DYNAMICS //
        Dynamics::Reset();
        
        //////////////////////
        // CREATE INFO FILE //
        //////////////////////
        
        if(MPIBasic::ID==0){
            CreateInfoFile();
        }
        
        
        ////////////////////////////
        // SET INITIAL CONDITIONS //
        ////////////////////////////
        
        // SET SU(N) GAUGE LINKS //
        GaugeLinks *UFinal=new GaugeLinks(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2]);
        InitialConditions::HandmadeSphaleron::Setup(GLinks::U,UFinal,EFields::E,GaugeTransformation::G);
        
        // SET U(1) GAUGE LINKS //
        QED::SetConstantMagneticField(BFlux);
        
        
        /* // OPTION TO PERFORM RANDOM U(1) GAUGE TRANSFORMATION //
         COMPLEX *G_A=new COMPLEX[Lattice::N[0]*Lattice::N[1]];
         
         for(INT y=0;y<Lattice::N[1];y++){
         for(INT x=0;x<Lattice::N[0];x++){
         
         DOUBLE arg=2.0*M_PI*RandomNumberGenerator::rng();
         G_A[QED::Index2D(x,y)]=std::exp(ComplexI*arg);
         
         }
         }
         
         QED::PerformGaugeTranformation(G_A);
         // END OPTION // */
        
        // SYNCHRONIZE //
        MPI_Barrier(MPI_COMM_WORLD);
        
        // SET INITIAL CONDITIONS FOR FERMIONS //
        
        // OPTION TO USE FREE FERMION MODE FUNCTIONS //
        //InitialConditions::SetFreeFermions(Fermions::Psi,Fermions::PsiMid);
        // END OPTION //
        
        // OPTION TO USE QED EIGENFUNCTIONS //
        // WILSON
        //InitialConditions::SetQEDEigenfunctions(Fermions::Psi,Fermions::PsiMid);
        
        // OVERLAP
        InitialConditions::SetOverlapQEDEigenfunctions(Fermions::Psi,Fermions::PsiMid);
        // END OPTION //
        
        // OPTION TO SAVE INITIAL CONFIGURATION //
        IO::SaveConfiguration(StringManipulation::StringCast("UOutT",Dynamics::Time()).c_str(),StringManipulation::StringCast("EOutT",Dynamics::Time()).c_str());
        // END OPTION //
        
        ////////////////////
        // TIME EVOLUTION //
        ////////////////////
        
        Evolve(1.5*tSphaleron);
        
        
    }
    
    
}
