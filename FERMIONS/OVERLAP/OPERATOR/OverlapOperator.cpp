//////////////////////
// COMPUTES -i h_ov //
//////////////////////
#include "RitzSolver.cpp"
#include "ZolotarevCoefficients.cpp"
#include "OverlapU.cpp"
#include "OverlapSignTest.cpp"


namespace OverlapOperator{

    // PARAMETERS FOR ZOLOTAREV COEFFICIENTS
    static const DOUBLE ZolotarevErrTol=std::pow(10.0,-12.0);
    static const INT NumofZolotarevCoeffs=25;
    
    // ERRORS FOR OVERLAP HAMILTONIAN MULTIPLE CONJUGATE GRID
    static const INT MaximumCGIter=1000;
    static const DOUBLE ConGradErr=std::pow(10.0,-16.0);
    
    // ZOLOTAREV COEFFICIENTS
    DOUBLE bz[NumofZolotarevCoeffs];
    DOUBLE dz[NumofZolotarevCoeffs];
    
    
    void Compute(u_cc *u,wcf *phi, wcf *psi,DOUBLE LambdaMin,DOUBLE LambdaMax){
        
        INT SizeOfVectorSpace=Lattice::N[0]*Lattice::N[1]*Lattice::N[2]*Nc*DiracAlgebra::SpinorComponents;
        
        COMPLEX* RandVector=new COMPLEX[SizeOfVectorSpace];
        
        DOUBLE Norm=0;
        
        for(INT i=0;i<SizeOfVectorSpace;i++){
            
            DOUBLE Re=RandomNumberGenerator::rng(); DOUBLE Im=RandomNumberGenerator::rng();
            
            Norm+=Re*Re+Im*Im;
            
            RandVector[i]=COMPLEX(Re,Im);
            
        }
        
        for(INT i=0;i<SizeOfVectorSpace;i++){
            
            RandVector[i]/=std::sqrt(Norm);
            
        }

        

        // CALCULATE ZOLOTAREV COEFFICIENTS //
        OverlapOperations::ZolotarevCoefficients(LambdaMin,LambdaMax,ZolotarevErrTol,NumofZolotarevCoeffs,bz,dz);
        
        // OPTION FOR  OVERLAP SIGN FUNCTION TEST //
        //EvolutionOperator::OverlapOperations::OverlapSignTest(u,bz,dz,RandVector,LambdaMin,LambdaMax,ConGradErr,NumofZolotarevCoeffs,SizeOfVectorSpace,MaximumCGIter);
        // END OPTION //

        // CALCULATE OVERLAP EVOLUTION OPERATOR U //
        OverlapOperations::OverlapU(u,phi,psi,bz,dz,NumofZolotarevCoeffs,SizeOfVectorSpace,LambdaMin,LambdaMax,MaximumCGIter,ConGradErr);
        
        // CLEAR MEMORY
        delete RandVector;
        
    }
    
}



