#ifndef __COULOMBGAUGEFIXING__CPP__
#define __COULOMBGAUGEFIXING__CPP__

#define GAUGE_FIXING_ALGORITHM_FOURIER_ACCELERATION

namespace CoulombGaugeFixing{
    
    //GAUGE VIOLATION MEASURES
    DOUBLE GlobalMaxDeviation;
    DOUBLE GlobalAvgDeviation;
    
    DOUBLE GlobalMaxStepSize;
    DOUBLE GlobalAvgStepSize;
    
    DOUBLE GlobalMinGaugeFunctional;
    DOUBLE GlobalSumGaugeFunctional;
    
    //INCLUDE GAUGE FIXING TOOL SETS
    #include "CoulombGaugeDeviation.cpp"
    
    //INCLUDE GAUGE FIXING ALGORITHMS
    #include "ALGORITHMS/LosAlamos.cpp"
    #include "ALGORITHMS/Jacobi.cpp"
    #include "ALGORITHMS/FourierAcceleration.cpp"
    
    void Reset(){
        
        GaugeTransformation::SetIdentity();
        
    }
    
    void SetCoulombGauge(std::string fname, INT TotalSteps){
        
        //COPY CURRENT FIELDS
        GaugeTransformation::Operations::CopyFields();
        
        //APPLY PREVIOUS GAUGE TRANSFORMATION
        GaugeTransformation::Operations::GaugeTransformLinks();
        
        //PERFORM GAUGE FIXING
        INT nSteps=0;
        
        if(MPIBasic::ID==0){
            std::cerr << "#STARTING COULOMB GAUGE FIXING" << std::endl;
        }
        
        //CREATE OUTPUT STREAM FOR LOG FILE //
        std::ofstream OutStream;
        
        if(MPIBasic::ID==0){
            std::string OutputFile=StringManipulation::StringCast(IO::OutputDirectory,fname,"ID",RandomNumberGenerator::MySEED,".log");        
            OutStream.open(OutputFile.c_str());
        }
        
        
        //RESET BREAK CONDITION //
        GlobalMaxDeviation=1;
        
        while(nSteps<TotalSteps && GlobalMaxDeviation>std::pow(10.0,-15)){
            
            // PERFORM UPDATE //
            #ifdef GAUGE_FIXING_ALGORITHM_LOS_ALAMOS
            LosAlamosAlgorithm::UpdateGaugeTransformation();
            #endif
            
            #ifdef GAUGE_FIXING_ALGORITHM_JACOBI
            JacobiAlgorithm::UpdateGaugeTransformation();
            #endif
            
            #ifdef GAUGE_FIXING_ALGORITHM_FOURIER_ACCELERATION
            FourierAccelerationAlgorithm::UpdateGaugeTransformation();
            #endif
            
            // MONITOR PRECISION //
            if(MPIBasic::ID==0){
                OutStream << nSteps << " " << GlobalMaxDeviation << " " << GlobalAvgDeviation << " " << GlobalMaxStepSize << " " <<  GlobalAvgStepSize << " " <<  GlobalMinGaugeFunctional << " " << GlobalSumGaugeFunctional << std::endl;
            }
            // INCREASE STEP COUNTER //
            nSteps++;
            
        }
        
        //COMMAND LINE OUTPUT
        if(MPIBasic::ID==0){
            std::cerr << "#GAUGE FIXING PRECISION AFTER FINAL ITERATION" << std::endl;
            std::cerr << nSteps << " " << GlobalMaxDeviation << " " << GlobalAvgDeviation << " " << GlobalMaxStepSize << " " <<  GlobalAvgStepSize << " " <<  GlobalMinGaugeFunctional << " " << GlobalSumGaugeFunctional << std::endl;
        }
        
        // CLOSE OUTSTREAM //
        if(MPIBasic::ID==0){
            OutStream.close();
        }
        
        //APPLY GAUGE TRANSFORMATION
        GaugeTransformation::Operations::GaugeTransformLinks();
        GaugeTransformation::Operations::GaugeTransformElectricFields();
    }
    
    //OVERLOADED
    void SetCoulombGauge(std::string fname){
        
        SetCoulombGauge(fname,1000);
    }
    
}


#endif
