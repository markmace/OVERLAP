#include "BASIS/QEDEigenfunctions.cpp"

namespace InitialConditions{
    
    void SetQEDEigenfunctions(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,FermionField *Psi,FermionField *PsiMid){
        
        DOUBLE Normalization=std::sqrt(Psi->N[0]*Psi->N[1]);
        
        for(INT s=0;s<Psi->NumberOfSamples;s++){
            
            // DETERMINE WHICH MODE FUNCTION TO USE //
            INT LevelIndex,pZIndex,ColorIndex,ParticleIndex;
            Psi->GetLocalModeQuantumNumbers(s,LevelIndex,pZIndex,ColorIndex,ParticleIndex);
                        
            // GET CORRESPONDING EIGENFUNCTION // 
            INT BasisIndex;
            
            if(ParticleIndex==0){
                BasisIndex=(DiracAlgebra::SpinorComponents*Psi->N[0]*Psi->N[1])/2 -1 -LevelIndex;
            }
            else{
                BasisIndex=(DiracAlgebra::SpinorComponents*Psi->N[0]*Psi->N[1])/2 + LevelIndex;
            }
            
            DOUBLE Frequency; COMPLEX *eVec=new COMPLEX[DiracAlgebra::SpinorComponents*Psi->N[0]*Psi->N[1]];
            QEDEigenfunctions::EigenModes[pZIndex]->GetEigenMode(BasisIndex,Frequency,eVec);
            
            // CHECK SIGN OF FREQUENCY //
            if((ParticleIndex==0 && Frequency<0) || (ParticleIndex==1 && Frequency>0)){
                std::cerr << "#ERROR -- " << ParticleIndex << " " << Frequency << std::endl;
            }
            
            #pragma omp parallel for
            for(INT z=zLow;z<=zHigh;z++){
                
                COMPLEX ComplexPhase=exp(+2.0*M_PI*ComplexI*DOUBLE(pZIndex)*DOUBLE(z)/DOUBLE(Lattice::N[2]));
                
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){
                            
                            for(INT i=0;i<Nc;i++){
                                
                                if(i==ColorIndex){
                                    Psi->Get(x,y,z,alpha,i,s)[0]=ComplexPhase*Normalization*eVec[QEDEigenfunctions::VectorIndex(x,y,alpha)];
                                    PsiMid->Get(x,y,z,alpha,i,s)[0]=ComplexPhase*Normalization*eVec[QEDEigenfunctions::VectorIndex(x,y,alpha)]*std::exp(-0.5*ComplexI*Frequency*Dynamics::dTau);
                                }
                                else{
                                    Psi->Get(x,y,z,alpha,i,s)[0]=0.0;
                                    PsiMid->Get(x,y,z,alpha,i,s)[0]=0.0;
                                }
                                
                            }
                            
                        }
                        
                    }
                }
            }
            
            
            
        }
        
        
    }
    
    void SetQEDEigenfunctions(FermionField *Psi,FermionField *PsiMid){
        
        
        // COMMANDLINE OUTPUT //
        if(MPIBasic::ID==0){
            std::cerr << "#SETTING FERMION FIELDS IN BACKGROUND MAGNETIC FIELD" << std::endl;
        }
        
        // COMPUTE EIGENFUNCTIONS //
        QEDEigenfunctions::Setup(QED::U,INPUT_QED_EIGENFUNCTIONS,SAVE_QED_EIGENFUNCTIONS);
        
        // SYNCRHONIZE //
        MPI_Barrier(MPI_COMM_WORLD);
        
        // SET INITIAL FERMION FIELDS //
        SetQEDEigenfunctions(0,Psi->N[0]-1,0,Psi->N[1]-1,0,Psi->N[2]-1,Psi,PsiMid);
        
        // SYNCRHONIZE //
        MPI_Barrier(MPI_COMM_WORLD);
        
        // COMMANDLINE OUTPUT //
        if(MPIBasic::ID==0){
            std::cerr << "#DONE -- SETTING INITIAL CONDITIONS" << std::endl;
        }
    }

}