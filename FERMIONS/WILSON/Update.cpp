
namespace Dynamics{
    
    namespace Fermions{
        
        void SingleModeUpdate(FermionField *Psi,FermionField* NewPsi,INT sIndex,wcf *DerivativePsi){
            
            // COMPUTE EVOLUTION OPERATOR ACTING ON PSI //
            EvolutionOperator::ComputeSingleMode(Psi,sIndex,DerivativePsi);
            
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
        
        
        void PerformUpdate(FermionField *Psi,FermionField *NewPsi){
            
            // GET DERIVATIVE BUFFER //
            wcf *DerivativePsi=new wcf[Psi->Volume];                
            
            // UPDATE ALL MODES //
            for(INT s=0;s<Psi->NumberOfSamples;s++){
                
                SingleModeUpdate(Psi,NewPsi,s,DerivativePsi);
                
            }
            
            // CLEAR BUFFER //
            delete DerivativePsi;
            
        }
        
        
        void Update(GaugeLinks *U,ElectricFields *E,FermionField *Psi,FermionField *PsiMid){
            
            // SETUP EVOLUTION OPERATOR //
            EvolutionOperator::Setup(U);
            
            // EVOLVE FERMION FIELD AT HALF-INTEGER STEPS //
            PerformUpdate(PsiMid,Psi);
            
            // EVOLVE GAUGE LINKS FOR HALF STEP //
            if(Dynamics::Time()<tSphaleron){
                Dynamics::UpdateGaugeLinks(0,GLinks::U->N[0]-1,0,GLinks::U->N[1]-1,0,GLinks::U->N[2]-1,GLinks::U,EFields::E,0.5*dTau);
            }
            
            // SETUP EVOLUTION OPERATOR //
            EvolutionOperator::Setup(U);

            // EVOLVE FERMION FIELD AT INTEGER STEPS //
            PerformUpdate(Psi,PsiMid);
            
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