
namespace Dynamics{
    
    namespace Fermions{
        
        // GLOBAL EIGENVALUES FOR RITZ SOLVER //
        DOUBLE MinEValue,MaxEValue;
        
        void SingleOverlapModeUpdate(FermionField *Psi,FermionField* NewPsi,INT sIndex,wcf *DerivativePsi,DOUBLE LambdaMin,DOUBLE LambdaMax){
            
            // COMPUTE EVOLUTION OPERATOR ACTING ON PSI //
            EvolutionOperator::ComputeSingleOverlapMode(Psi,sIndex,DerivativePsi,LambdaMin,LambdaMax);
            
            // UPDATE NEW PSI //
            for(INT z=0;z<=Psi->N[2]-1;z++){
                for(INT y=0;y<=Psi->N[1]-1;y++){
                    for(INT x=0;x<=Psi->N[0]-1;x++){
                        
                        for(INT i=0;i<Nc;i++){
                            for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){
                                
                                NewPsi->Get(x,y,z,alpha,i,sIndex)[0]+=dTau*DerivativePsi[Psi->Index3D(x,y,z)].c[i].f[alpha];
                            }
                        }
                        
                    }
                }
            }
            
        }
        
        
        void PerformOverlapUpdate(FermionField *Psi,FermionField *NewPsi,DOUBLE LambdaMin,DOUBLE LambdaMax){
            
            // GET DERIVATIVE BUFFER //
            wcf *DerivativePsi=new wcf[Psi->Volume];                
            
            // UPDATE ALL MODES //
            for(INT s=0;s<Psi->NumberOfSamples;s++){
                
                SingleOverlapModeUpdate(Psi,NewPsi,s,DerivativePsi,LambdaMin,LambdaMax);
                
            }
            
            // CLEAR BUFFER //
            delete DerivativePsi;
            
        }
        
        // INCLUDE RITZ FOR EIGENVALUES
        #include "OPERATOR/RitzSolver.cpp"
        
        // PARAMETERS FOR RITZ COMPUTATION
        INT NumberOfRestarts=100;
        INT MaximumIterations=2000;
        DOUBLE NormRatio=std::pow(10.0,-16.0);
        
        // NUMBER OF EIGENVALUES FOR RITZ
        INT NumberOfEVals=1;
        
        void ComputeRitzEigenvalues(){
            
            // DEFINE SIZE OF VECTOR SPACE
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
            
            // CALCULATE MAXIMUM EIGENVALUE OF WILSON HAMILTONIAN //
            // ASSIGN VALUE -- NEEDS NEGATIVE SIGN TO BECAUSE USE OF INVERSE OF ALGORITHM //
            MaxEValue=-DOUBLE(1.05)*OverlapOperations::RitzSolver(EvolutionOperator::Ucc,RandVector,NumberOfEVals,NormRatio,MaximumIterations,-1.0,NumberOfRestarts,SizeOfVectorSpace);
            
            // CALCULATE MINIMUM EIGENVALUE OF WILSON HAMILTONIAN //
            // ASSIGN VALUE //
            MinEValue=DOUBLE(0.95)*OverlapOperations::RitzSolver(EvolutionOperator::Ucc,RandVector,NumberOfEVals,NormRatio,MaximumIterations,1.0,NumberOfRestarts,SizeOfVectorSpace);
  
            // CLEAR MEMORY
            delete RandVector;
            
        }
        
        void Update(GaugeLinks *U,ElectricFields *E,FermionField *Psi,FermionField *PsiMid){
            
            // SETUP EVOLUTION OPERATOR //
            EvolutionOperator::Setup(U);
            
            // CALCULATE EIGEVANLUES VIA RITZ
            ComputeRitzEigenvalues();
                        
            // EVOLVE FERMION FIELD AT HALF-INTEGER STEPS //
            PerformOverlapUpdate(PsiMid,Psi,MinEValue,MaxEValue);
            
            // EVOLVE GAUGE LINKS FOR HALF STEP //
            if(Dynamics::Time()<tSphaleron){
                Dynamics::UpdateGaugeLinks(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1,GLinks::U,EFields::E,0.5*dTau);
            }
            
            // SETUP EVOLUTION OPERATOR //
            EvolutionOperator::Setup(U);
            
            // CALCULATE EIGEVANLUES VIA RITZ
            ComputeRitzEigenvalues();
                        
            // EVOLVE FERMION FIELD AT INTEGER STEPS //
            PerformOverlapUpdate(Psi,PsiMid,MinEValue,MaxEValue);
            
            // EVOLVE GAUGE LINKS FOR HALF STEP //
            if(Dynamics::Time()<tSphaleron){
                Dynamics::UpdateGaugeLinks(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1,GLinks::U,EFields::E,0.5*dTau);
            }
            
            // INCREASE STEP COUNTER //
            tSteps++;
            
            //SYNCHRONIZE ALL MPI NODES
            MPI_Barrier(MPI_COMM_WORLD);
            
        }
        
    }
}