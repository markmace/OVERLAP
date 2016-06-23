namespace InitialConditions{
    
    void SetFreeFermions(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,FermionField *Psi,FermionField *PsiMid){
        
        #pragma omp parallel for
        for(INT z=zLow;z<=zHigh;z++){
            for(INT y=yLow;y<=yHigh;y++){
                for(INT x=xLow;x<=xHigh;x++){
                    
                    for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){
                        for(INT i=0;i<Nc;i++){
                            
                            for(INT s=0;s<Psi->NumberOfSamples;s++){
                                
                                // DETERMINE WHICH MODE FUNCTION TO USE //
                                INT pXIndex,pYIndex,pZIndex,SpinIndex,ColorIndex,ParticleIndex;
                                Psi->GetLocalModeQuantumNumbers(s,pXIndex,pYIndex,pZIndex,SpinIndex,ColorIndex,ParticleIndex);
                                
                                // GET DIRAC SPINOR AND COMPLEX PHASE FACTOR //
                                COMPLEX ComplexPhase; COMPLEX DiracSpinor[4]; DOUBLE Frequency;
                                #if FERMION_FLAG==WILSON_FLAG
                                    Psi->GetWilsonDiracSpinor(pXIndex,pYIndex,pZIndex,SpinIndex,ParticleIndex,DiracSpinor,Frequency);
                                    Psi->GetWilsonPhase(pXIndex,pYIndex,pZIndex,x,y,z,ParticleIndex,ComplexPhase);
                                #endif
                                
                                #if FERMION_FLAG==OVERLAP_FLAG
                                    Psi->GetOverlapDiracSpinor(pXIndex,pYIndex,pZIndex,SpinIndex,ParticleIndex,DiracSpinor,Frequency);
                                    Psi->GetOverlapPhase(pXIndex,pYIndex,pZIndex,x,y,z,ParticleIndex,ComplexPhase);
                                #endif
                                                                
                                if((ParticleIndex==0 && Frequency<0) || (ParticleIndex==1 && Frequency>0)){
                                    std::cerr << "#NOOOOOOOOO" << std::endl;
                                    exit(0);
                                }
                                
                                // SET INITIAL CONDITION //
                                if(i==ColorIndex){
                                    Psi->Get(x,y,z,alpha,i,s)[0]=ComplexPhase*DiracSpinor[alpha];
                                    PsiMid->Get(x,y,z,alpha,i,s)[0]=ComplexPhase*DiracSpinor[alpha]*std::exp(-0.5*ComplexI*Frequency*Dynamics::dTau);
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
    
    
    void SetFreeFermions(FermionField *Psi,FermionField *PsiMid){
        
        // COMMANDLINE OUTPUT //
        if(MPIBasic::ID==0){
            std::cerr << "#SETTING FREE FERMION FIELDS" << std::endl;
        }
        
        SetFreeFermions(0,Psi->N[0]-1,0,Psi->N[1]-1,0,Psi->N[2]-1,Psi,PsiMid);
        
        // SYNCRHONIZE //
        MPI_Barrier(MPI_COMM_WORLD);
        
        // COMMANDLINE OUTPUT //
        if(MPIBasic::ID==0){
            std::cerr << "#DONE -- SETTING INITIAL CONDITIONS" << std::endl;
        }
    }
    
}