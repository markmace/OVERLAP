#ifndef QED_EIGENFUNCTIONS_CPP
#define QED_EIGENFUNCTIONS_CPP

#define PERIODICDELTA(x,y,N) (((x)==(y) || (x)==(y)+(N) || (x)==(y)-(N))?1.0:0.0)

namespace QEDEigenfunctions{
    
    // EIGENFUNCTIONS IN THE PRESENCE OF CONSTANT MAGNETIC FIELD //
    Eigensystem **EigenModes;
    
    int VectorIndex(INT x,INT y,INT alpha){
        
        return alpha+DiracAlgebra::SpinorComponents*(x+Lattice::N[0]*y);
        
    }
    
    // INDEXING OF WILSON OPERATOR //
    int OperatorIndex(INT x1,INT y1,INT x2,INT y2,INT alpha,INT beta){
        
        return VectorIndex(x2,y2,beta)+DiracAlgebra::SpinorComponents*Lattice::N[0]*Lattice::N[1]*VectorIndex(x1,y1,alpha);
    }
        
    // COMPUTE EIGENFUNCTIONS //
    void Setup(COMPLEX *U_A,INT INPUT_FLAG,INT OUTPUT_FLAG){
        
        // COMMANDLINE OUTPUT //
        if(MPIBasic::ID==0){
            std::cerr << "#COMPUTING EIGENFUNCTIONS IN MAGNETIC FIELD" << std::endl;
        }
        
        
        
        // DIRAC SPINOR COUPLINGS //
        COMPLEX g0[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        
        COMPLEX gXu_1[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        COMPLEX gXd_1[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        
        COMPLEX gYu_1[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        COMPLEX gYd_1[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        
        COMPLEX gZu_1[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        COMPLEX gZd_1[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        
        
        COMPLEX gXu_2[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        COMPLEX gXd_2[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        
        COMPLEX gYu_2[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        COMPLEX gYd_2[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        
        COMPLEX gZu_2[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        COMPLEX gZd_2[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        
        
        COMPLEX gXu_3[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        COMPLEX gXd_3[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        
        COMPLEX gYu_3[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        COMPLEX gYd_3[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        
        COMPLEX gZu_3[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        COMPLEX gZd_3[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
        
        
        
        // SET DIRAC SPINOR COUPLINGS //
        for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){
            for(INT beta=0;beta<DiracAlgebra::SpinorComponents;beta++){
                
                using namespace WilsonCoefficients;
                
                // COEFFICIENTS FOR WILSON TERMS //
                DOUBLE cWilson=(WilsonCoefficients::C1+2.0*WilsonCoefficients::C2+3.0*WilsonCoefficients::C3);
                
                // FERMION MASS AND SCALAR WILSON TERM //
                g0[alpha][beta]=(mFermion+(cWilson*rWilson)/Lattice::a[0]+(cWilson*rWilson)/Lattice::a[1]+(cWilson*rWilson)/Lattice::a[2])*DiracAlgebra::Gamma0[alpha][beta];
                
                gXu_1[alpha][beta]=-0.5*(1.0*C1*rWilson)/Lattice::a[0]*DiracAlgebra::Gamma0[alpha][beta];
                gXd_1[alpha][beta]=-0.5*(1.0*C1*rWilson)/Lattice::a[0]*DiracAlgebra::Gamma0[alpha][beta];
                
                gYu_1[alpha][beta]=-0.5*(1.0*C1*rWilson)/Lattice::a[1]*DiracAlgebra::Gamma0[alpha][beta];
                gYd_1[alpha][beta]=-0.5*(1.0*C1*rWilson)/Lattice::a[1]*DiracAlgebra::Gamma0[alpha][beta];
                
                gZu_1[alpha][beta]=-0.5*(1.0*C1*rWilson)/Lattice::a[2]*DiracAlgebra::Gamma0[alpha][beta];
                gZd_1[alpha][beta]=-0.5*(1.0*C1*rWilson)/Lattice::a[2]*DiracAlgebra::Gamma0[alpha][beta];
                
                
                gXu_2[alpha][beta]=-0.5*(2.0*C2*rWilson)/Lattice::a[0]*DiracAlgebra::Gamma0[alpha][beta];
                gXd_2[alpha][beta]=-0.5*(2.0*C2*rWilson)/Lattice::a[0]*DiracAlgebra::Gamma0[alpha][beta];
                
                gYu_2[alpha][beta]=-0.5*(2.0*C2*rWilson)/Lattice::a[1]*DiracAlgebra::Gamma0[alpha][beta];
                gYd_2[alpha][beta]=-0.5*(2.0*C2*rWilson)/Lattice::a[1]*DiracAlgebra::Gamma0[alpha][beta];
                
                gZu_2[alpha][beta]=-0.5*(2.0*C2*rWilson)/Lattice::a[2]*DiracAlgebra::Gamma0[alpha][beta];
                gZd_2[alpha][beta]=-0.5*(2.0*C2*rWilson)/Lattice::a[2]*DiracAlgebra::Gamma0[alpha][beta];
                
                
                gXu_3[alpha][beta]=-0.5*(3.0*C3*rWilson)/Lattice::a[0]*DiracAlgebra::Gamma0[alpha][beta];
                gXd_3[alpha][beta]=-0.5*(3.0*C3*rWilson)/Lattice::a[0]*DiracAlgebra::Gamma0[alpha][beta];
                
                gYu_3[alpha][beta]=-0.5*(3.0*C3*rWilson)/Lattice::a[1]*DiracAlgebra::Gamma0[alpha][beta];
                gYd_3[alpha][beta]=-0.5*(3.0*C3*rWilson)/Lattice::a[1]*DiracAlgebra::Gamma0[alpha][beta];
                
                gZu_3[alpha][beta]=-0.5*(3.0*C3*rWilson)/Lattice::a[2]*DiracAlgebra::Gamma0[alpha][beta];
                gZd_3[alpha][beta]=-0.5*(3.0*C3*rWilson)/Lattice::a[2]*DiracAlgebra::Gamma0[alpha][beta];
                
                
                // SPATIAL DERIVATIVE TERMS //
                for(INT gamma=0;gamma<DiracAlgebra::SpinorComponents;gamma++){
                    
                    gXu_1[alpha][beta]-=0.5*ComplexI*C1/Lattice::a[0]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaX[gamma][beta];
                    gXd_1[alpha][beta]+=0.5*ComplexI*C1/Lattice::a[0]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaX[gamma][beta];
                    
                    gYu_1[alpha][beta]-=0.5*ComplexI*C1/Lattice::a[1]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaY[gamma][beta];
                    gYd_1[alpha][beta]+=0.5*ComplexI*C1/Lattice::a[1]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaY[gamma][beta];
                    
                    gZu_1[alpha][beta]-=0.5*ComplexI*C1/Lattice::a[2]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaZ[gamma][beta];
                    gZd_1[alpha][beta]+=0.5*ComplexI*C1/Lattice::a[2]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaZ[gamma][beta];
                    
                    
                    gXu_2[alpha][beta]-=0.5*ComplexI*C2/Lattice::a[0]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaX[gamma][beta];
                    gXd_2[alpha][beta]+=0.5*ComplexI*C2/Lattice::a[0]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaX[gamma][beta];
                    
                    gYu_2[alpha][beta]-=0.5*ComplexI*C2/Lattice::a[1]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaY[gamma][beta];
                    gYd_2[alpha][beta]+=0.5*ComplexI*C2/Lattice::a[1]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaY[gamma][beta];
                    
                    gZu_2[alpha][beta]-=0.5*ComplexI*C2/Lattice::a[2]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaZ[gamma][beta];
                    gZd_2[alpha][beta]+=0.5*ComplexI*C2/Lattice::a[2]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaZ[gamma][beta];
                    
                    
                    gXu_3[alpha][beta]-=0.5*ComplexI*C3/Lattice::a[0]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaX[gamma][beta];
                    gXd_3[alpha][beta]+=0.5*ComplexI*C3/Lattice::a[0]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaX[gamma][beta];
                    
                    gYu_3[alpha][beta]-=0.5*ComplexI*C3/Lattice::a[1]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaY[gamma][beta];
                    gYd_3[alpha][beta]+=0.5*ComplexI*C3/Lattice::a[1]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaY[gamma][beta];
                    
                    gZu_3[alpha][beta]-=0.5*ComplexI*C3/Lattice::a[2]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaZ[gamma][beta];
                    gZd_3[alpha][beta]+=0.5*ComplexI*C3/Lattice::a[2]*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaZ[gamma][beta];
                    
                }
                
                
            }
        }

        
        // COMPUTE COMPLETE SET OF EIGENFUNCTIONS //
        EigenModes=new Eigensystem*[Lattice::N[2]];
        
        for(INT pZIndex=0;pZIndex<Lattice::N[2];pZIndex++){
            
            
            /////////////////////
            // SET HAMILTONIAN //
            /////////////////////

            // SET TIMER //
            if(MPIBasic::ID==0){
                std::cerr << "#COMPUTING WILSON DIRAC OPERATOR FOR pZ= " << pZIndex << std::endl;
                Timing::Reset();
            }
            
            using namespace WilsonCoefficients;
            
            // GET MOMENTUM EIGENVALUES //
            COMPLEX pZu_1=exp(+2.0*M_PI*ComplexI*DOUBLE(pZIndex)/DOUBLE(Lattice::N[2]));
            COMPLEX pZd_1=exp(-2.0*M_PI*ComplexI*DOUBLE(pZIndex)/DOUBLE(Lattice::N[2]));
            
            COMPLEX pZu_2=exp(+4.0*M_PI*ComplexI*DOUBLE(pZIndex)/DOUBLE(Lattice::N[2]));
            COMPLEX pZd_2=exp(-4.0*M_PI*ComplexI*DOUBLE(pZIndex)/DOUBLE(Lattice::N[2]));
            
            COMPLEX pZu_3=exp(+6.0*M_PI*ComplexI*DOUBLE(pZIndex)/DOUBLE(Lattice::N[2]));
            COMPLEX pZd_3=exp(-6.0*M_PI*ComplexI*DOUBLE(pZIndex)/DOUBLE(Lattice::N[2]));
            
            // DIMENSION OF WILSON OPERATOR //
            INT OperatorDimension=Lattice::N[0]*Lattice::N[1]*DiracAlgebra::SpinorComponents;
            
            // CONSTRUCT WILSON OPERATOR //
            COMPLEX *WilsonOperator=new COMPLEX[SQR(OperatorDimension)];
            
            COMPLEX wOO,wDO,wUO,wOD,wOU; 
            COMPLEX w2DO,w2UO,wO2D,wO2U;
            COMPLEX w3DO,w3UO,wO3D,wO3U;
            
            for(INT x1=0;x1<Lattice::N[0];x1++){
                for(INT y1=0;y1<Lattice::N[1];y1++){
                    
                    for(INT x2=0;x2<Lattice::N[0];x2++){
                        for(INT y2=0;y2<Lattice::N[1];y2++){
                            
                            for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){
                                for(INT beta=0;beta<DiracAlgebra::SpinorComponents;beta++){
                                    
                                    // SET EFFECITVELY LOCAL TERM -- MASS TERM PLUS Z-DERIVATIVE //
                                    wOO=PERIODICDELTA(x1,x2,Lattice::N[0])*PERIODICDELTA(y1,y2,Lattice::N[1])*(g0[alpha][beta] + gZu_1[alpha][beta]*pZu_1+gZd_1[alpha][beta]*pZd_1 + gZu_2[alpha][beta]*pZu_2+gZd_2[alpha][beta]*pZd_2 + gZu_3[alpha][beta]*pZu_3+gZd_3[alpha][beta]*pZd_3);
                                    
                                    // SET x-DERIVATIVE //
                                    wDO =PERIODICDELTA(x1-1,x2,Lattice::N[0])*PERIODICDELTA(y1,y2,Lattice::N[1])*gXd_1[alpha][beta]*conj(U_A[QED::Index2D(x1-1,y1,0)]);
                                    w2DO=PERIODICDELTA(x1-2,x2,Lattice::N[0])*PERIODICDELTA(y1,y2,Lattice::N[1])*gXd_2[alpha][beta]*conj(U_A[QED::Index2D(x1-2,y1,0)]*U_A[QED::Index2D(x1-1,y1,0)]);
                                    w3DO=PERIODICDELTA(x1-3,x2,Lattice::N[0])*PERIODICDELTA(y1,y2,Lattice::N[1])*gXd_3[alpha][beta]*conj(U_A[QED::Index2D(x1-3,y1,0)]*U_A[QED::Index2D(x1-2,y1,0)]*U_A[QED::Index2D(x1-1,y1,0)]);
                                    
                                    wUO =PERIODICDELTA(x1+1,x2,Lattice::N[0])*PERIODICDELTA(y1,y2,Lattice::N[1])*gXu_1[alpha][beta]*(U_A[QED::Index2D(x1,y1,0)]);
                                    w2UO=PERIODICDELTA(x1+2,x2,Lattice::N[0])*PERIODICDELTA(y1,y2,Lattice::N[1])*gXu_2[alpha][beta]*(U_A[QED::Index2D(x1,y1,0)]*U_A[QED::Index2D(x1+1,y1,0)]);
                                    w3UO=PERIODICDELTA(x1+3,x2,Lattice::N[0])*PERIODICDELTA(y1,y2,Lattice::N[1])*gXu_3[alpha][beta]*(U_A[QED::Index2D(x1,y1,0)]*U_A[QED::Index2D(x1+1,y1,0)]*U_A[QED::Index2D(x1+2,y1,0)]);
                                    
                                    //SET y-DERIVATIVE //
                                    wOD =PERIODICDELTA(x1,x2,Lattice::N[0])*PERIODICDELTA(y1-1,y2,Lattice::N[1])*gYd_1[alpha][beta]*conj(U_A[QED::Index2D(x1,y1-1,1)]);
                                    wO2D=PERIODICDELTA(x1,x2,Lattice::N[0])*PERIODICDELTA(y1-2,y2,Lattice::N[1])*gYd_2[alpha][beta]*conj(U_A[QED::Index2D(x1,y1-2,1)]*U_A[QED::Index2D(x1,y1-1,1)]);
                                    wO3D=PERIODICDELTA(x1,x2,Lattice::N[0])*PERIODICDELTA(y1-3,y2,Lattice::N[1])*gYd_3[alpha][beta]*conj(U_A[QED::Index2D(x1,y1-3,1)]*U_A[QED::Index2D(x1,y1-2,1)]*U_A[QED::Index2D(x1,y1-1,1)]);
                                    
                                    
                                    wOU =PERIODICDELTA(x1,x2,Lattice::N[0])*PERIODICDELTA(y1+1,y2,Lattice::N[1])*gYu_1[alpha][beta]*(U_A[QED::Index2D(x1,y1,1)]);
                                    wO2U=PERIODICDELTA(x1,x2,Lattice::N[0])*PERIODICDELTA(y1+2,y2,Lattice::N[1])*gYu_2[alpha][beta]*(U_A[QED::Index2D(x1,y1,1)]*U_A[QED::Index2D(x1,y1+1,1)]);
                                    wO3U=PERIODICDELTA(x1,x2,Lattice::N[0])*PERIODICDELTA(y1+3,y2,Lattice::N[1])*gYu_3[alpha][beta]*(U_A[QED::Index2D(x1,y1,1)]*U_A[QED::Index2D(x1,y1+1,1)]*U_A[QED::Index2D(x1,y1+2,1)]);
                                    
                                    
                                    // SET MATRIX ENTRY //
                                    WilsonOperator[OperatorIndex(x1,y1,x2,y2,alpha,beta)]=wOO+wDO+w2DO+w3DO+wUO+w2UO+w3UO+wOD+wO2D+wO3D+wOU+wO2U+wO3U;
                                    
                                }
                                
                            }
                            
                        }
                    }
                    
                }
            }
            
            // TIMING //
            if(MPIBasic::ID==0){
                std::cerr << "#TIMING -- " << Timing::Get() << std::endl;
                Timing::Reset();
            }
            
            
            ///////////////////////
            // SETUP EIGENSYSTEM //
            ///////////////////////
            
            EigenModes[pZIndex]=new Eigensystem(WilsonOperator,OperatorDimension);            
            
            ////////////////////////////////////////////
            // COMPUTE EIGENSYSTEM BY DIAGONALIZATION //
            ////////////////////////////////////////////
            
            if(INPUT_FLAG==0){
            
                // COMMANDLINE OUTPUT //
                if(MPIBasic::ID==0){
                    std::cerr << "#DIAGONALIZING DIRAC OPERATOR FOR pZ= " << pZIndex << std::endl;
                    Timing::Reset();
                }
                
                // DIAGONALIZE //
                EigenModes[pZIndex]->Diagonalize();
                
                // TIMING //
                if(MPIBasic::ID==0){
                    std::cerr << "#TIMING -- " << Timing::Get() << std::endl;
                    Timing::Reset();
                }
                
            }
            
            /////////////////////////////////////
            // GET EIGENSYSTEM FROM INPUT FILE //
            /////////////////////////////////////
            
            if(INPUT_FLAG==1){
                  
                // COMMANDLINE OUTPUT //
                if(MPIBasic::ID==0){
                    std::cerr << "#IMPORTING EIGENMODES OF DIRAC OPERATOR FOR pZ= " << pZIndex << std::endl;
                    Timing::Reset();
                }
                
                // GET EIGENVALUE AND EIGENVECTOR //
                DOUBLE eVal; COMPLEX *eVec=new COMPLEX[OperatorDimension];
                
                // OPEN INPUT FILE ON MASTER //
                std::ifstream InStream;
                if(MPIBasic::ID==0){
                    InStream.open(StringManipulation::StringCast(IO::OutputDirectory,"DiracSpectrum_pZ",pZIndex,".txt").c_str());
                }
                
                // GET EIGENMODE //
                for(INT s=0;s<OperatorDimension;s++){
                    
                    // GET INPUT ON MASTER //
                    if(MPIBasic::ID==0){
                        
                        // GET EIGENVALUE //
                        std::string strBuffer,InputLine;
                        getline(InStream,InputLine);
                        
                        std::stringstream EValString(InputLine);
                        EValString >> strBuffer; EValString >> eVal;
                        
                        // GET EIGENVECTOR //
                        getline(InStream,InputLine);
                        std::stringstream EVecString(InputLine);
                        
                        DOUBLE Re,Im;
                        
                        for(INT k=0;k<OperatorDimension;k++){
                            
                            EVecString >> Re; EVecString >> Im;
                            eVec[k]=COMPLEX(Re,Im);
                        }
                        
                    }
                    
                    // BROADCAST //
                    MPI_Bcast(&eVal,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
                    MPI_Bcast(eVec,OperatorDimension,MPI::DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
                                        
                    // SET //
                    EigenModes[pZIndex]->SetEigenMode(s,eVal,eVec);
                    
                }
                
                
                // CLOSE INPUT FILE ON MASTER //
                if(MPIBasic::ID==0){
                    InStream.close();
                }
                
                // CLEAN-UP //
                delete eVec;
                
                // TIMING //
                if(MPIBasic::ID==0){
                    std::cerr << "#TIMING -- " << Timing::Get() << std::endl;
                    Timing::Reset();
                }
                
            }
            
            //////////////////////////////////////
            // CHECK CORRECTNESS OF EIGENSYSTEM //
            //////////////////////////////////////
            
            // COMMANDLINE OUTPUT //
            if(MPIBasic::ID==0){
                std::cerr << "#CHECKING EIGENSYSTEM FOR pZ= " << pZIndex << std::endl;
            }
            
            EigenModes[pZIndex]->Check();
            
            // TIMING //
            if(MPIBasic::ID==0){
                std::cerr << "#TIMING -- " << Timing::Get() << std::endl;
                Timing::Reset();
            }
            
            
            //////////////////////
            // SAVE EIGENSYSTEM //
            //////////////////////
            
            if(OUTPUT_FLAG==1){
            
                if(MPIBasic::ID==0){
                    
                    // COMMANDLINE OUTPUT //
                    std::cerr << "#SAVING SPECTRUM OF DIRAC OPERATOR FOR pZ= " << pZIndex << std::endl;
                    
                    // CREATE OUT-STREAM //
                    std::ofstream OutStream;
                    OutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"DiracSpectrum_pZ",pZIndex,".txt").c_str());
                    
                    OutStream.precision(OUTPUT_PRECISION);
                    
                    // ACCESS EIGENMODE //
                    DOUBLE eVal; COMPLEX *eVec=new COMPLEX[OperatorDimension];
                    
                    for(INT s=0;s<OperatorDimension;s++){
                        
                        // GET EIGENMODE //
                        EigenModes[pZIndex]->GetEigenMode(s,eVal,eVec);
                        
                        OutStream << "#EVal= " << eVal << std::endl;
                        
                        for(INT y=0;y<Lattice::N[1];y++){
                            for(INT x=0;x<Lattice::N[0];x++){
                                for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){
                                    
                                    OutStream << real(eVec[VectorIndex(x,y,alpha)]) << " " << imag(eVec[VectorIndex(x,y,alpha)]) << " ";
                                    
                                }
                                
                            }
                        }
                        
                        OutStream << std::endl;
                        
                    }
                    
                    // CLOSE OUT-STREAM //
                    OutStream.close();
                    
                    // CLEAN-UP //
                    delete eVec;
                    
                }
                
            }
            
            // CLEAN-UP //
            delete WilsonOperator;
            
        }
        
        
    }
    
    
    
    
}

#endif
