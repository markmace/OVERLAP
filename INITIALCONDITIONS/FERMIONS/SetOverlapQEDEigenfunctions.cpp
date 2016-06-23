#include "BASIS/OverlapQEDEigenfunctions.cpp"

namespace InitialConditions{
    
    void SetOverlapQEDEigenfunctions(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,FermionField *Psi,FermionField *PsiMid){
        
        // WARNING-- CHECK //
        DOUBLE Normalization=std::sqrt(Psi->N[0]*Psi->N[1]*Psi->N[2]);
        
        for(INT s=0;s<Psi->NumberOfSamples;s++){
            
            INT LevelIndex,ParticleIndex;
            Psi->GetLocalModeQuantumNumbers(s,LevelIndex,ParticleIndex);
            
            INT BasisIndex;
            
            if(ParticleIndex==0){
                BasisIndex=(DiracAlgebra::SpinorComponents*Psi->N[0]*Psi->N[1]*Psi->N[2]*Nc)/2 -1 -LevelIndex;
            }
            else{
                BasisIndex=(DiracAlgebra::SpinorComponents*Psi->N[0]*Psi->N[1]*Psi->N[2]*Nc)/2 + LevelIndex;
            }
            
            // GET CORRESPONDING EIGENFUNCTION //
            DOUBLE Frequency; COMPLEX *eVec=new COMPLEX[DiracAlgebra::SpinorComponents*Psi->N[0]*Psi->N[1]*Psi->N[2]*Nc];
            OverlapQEDEigenfunctions::EigenModes[0]->GetEigenMode(BasisIndex,Frequency,eVec);
            
            // CHECK SIGN OF FREQUENCY //
            if((ParticleIndex==0 && Frequency<0) || (ParticleIndex==1 && Frequency>0)){
                std::cerr << "#FREQ ERROR -- " << ParticleIndex << " " << BasisIndex << " " << Frequency << std::endl;
            }
            
            
            #pragma omp parallel for
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){
                            
                            for(INT i=0;i<Nc;i++){
                                
                                Psi->Get(x,y,z,alpha,i,s)[0]=Normalization*eVec[OverlapQEDEigenfunctions::VectorIndex(x,y,z,alpha,i)];
                                PsiMid->Get(x,y,z,alpha,i,s)[0]=Normalization*eVec[OverlapQEDEigenfunctions::VectorIndex(x,y,z,alpha,i)]*std::exp(-0.5*ComplexI*Frequency*Dynamics::dTau);
                                
                            }
                            
                        }
                        
                    }
                }
            }
            
            delete eVec;
            
            
            
        }
        
        
    }
    
    void SetOverlapQEDEigenfunctions(FermionField *Psi,FermionField *PsiMid){
        
        
        // COMMANDLINE OUTPUT //
        if(MPIBasic::ID==0){
            std::cerr << "#SETTING FERMION FIELDS IN BACKGROUND MAGNETIC FIELD" << std::endl;
        }
        
        // COMPUTE EIGENFUNCTIONS //
        OverlapQEDEigenfunctions::Setup(QED::U,INPUT_QED_EIGENFUNCTIONS,SAVE_QED_EIGENFUNCTIONS);
        
        // SYNCRHONIZE //
        MPI_Barrier(MPI_COMM_WORLD);
        
        // SET INITIAL FERMION FIELDS //
        SetOverlapQEDEigenfunctions(0,Psi->N[0]-1,0,Psi->N[1]-1,0,Psi->N[2]-1,Psi,PsiMid);
        
        // SYNCRHONIZE //
        MPI_Barrier(MPI_COMM_WORLD);
        
        // COMMANDLINE OUTPUT //
        if(MPIBasic::ID==0){
            std::cerr << "#DONE -- SETTING INITIAL CONDITIONS" << std::endl;
        }
    }
    
}