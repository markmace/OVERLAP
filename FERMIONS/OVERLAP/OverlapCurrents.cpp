namespace Fermions{
    
    namespace Observables{
        
        void ComputeCurrents(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,FermionField *Psi){
            
            // VECTOR AND AXIAL CURRENTS //
            DOUBLE *jV=new DOUBLE[4*Psi->Volume];
            DOUBLE *jA=new DOUBLE[4*Psi->Volume];
            
            ///////////////////////////////////
            // SET LATTICE CURRENT COUPLINGS //
            ///////////////////////////////////
            
            COMPLEX gVT[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gVX[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gVY[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gVZ[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            
            COMPLEX gAT[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gAX[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gAY[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gAZ[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            
            COMPLEX gVT2[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gVX2[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gVY2[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gVZ2[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            
            COMPLEX gAT2[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gAX2[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gAY2[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gAZ2[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            
            for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){
                for(INT beta=0;beta<DiracAlgebra::SpinorComponents;beta++){
                    
                    // RESET //
                    gVT[alpha][beta]=0.0; gVX[alpha][beta]=0.0; gVY[alpha][beta]=0.0; gVZ[alpha][beta]=0.0;
                    gAT[alpha][beta]=0.0; gAX[alpha][beta]=0.0; gAY[alpha][beta]=0.0; gAZ[alpha][beta]=0.0;
                    
                    gVT2[alpha][beta]=0.0; gVX2[alpha][beta]=0.0; gVY2[alpha][beta]=0.0; gVZ2[alpha][beta]=0.0;
                    gAT2[alpha][beta]=0.0; gAX2[alpha][beta]=0.0; gAY2[alpha][beta]=0.0; gAZ2[alpha][beta]=0.0;
                    
                    for(INT gamma=0;gamma<DiracAlgebra::SpinorComponents;gamma++){
                        
                        // STANDARD VECTOR TERMS //
                        // GAMMA_0 GAMMA_mu
                        
                        gVT[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::Gamma0[gamma][beta];
                        gVX[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaX[gamma][beta];
                        gVY[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaY[gamma][beta];
                        gVZ[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaZ[gamma][beta];
                      
                    
                        for(INT delta=0;delta<DiracAlgebra::SpinorComponents;delta++){
                            
                            // STANDARD AXIAL TERMS //
                            // GAMMA_0 GAMMA_mu GAMMA_5
                            gAT[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::Gamma0[gamma][delta]*DiracAlgebra::Gamma5[delta][beta];
                            gAX[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaX[gamma][delta]*DiracAlgebra::Gamma5[delta][beta];
                            gAY[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaY[gamma][delta]*DiracAlgebra::Gamma5[delta][beta];
                            gAZ[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaZ[gamma][delta]*DiracAlgebra::Gamma5[delta][beta];
                            
                            
                            // NON-STANDARD VECTOR TERMS //
                            // -i/2 GAMMA_0 GAMMMA_MU GAMMA_0
                            gVT2[alpha][beta]+=COMPLEX(0.0,-0.5)*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::Gamma0[gamma][delta]*DiracAlgebra::Gamma0[delta][beta];
                            gVX2[alpha][beta]+=COMPLEX(0.0,-0.5)*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaX[gamma][delta]*DiracAlgebra::Gamma0[delta][beta];
                            gVY2[alpha][beta]+=COMPLEX(0.0,-0.5)*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaY[gamma][delta]*DiracAlgebra::Gamma0[delta][beta];
                            gVZ2[alpha][beta]+=COMPLEX(0.0,-0.5)*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaZ[gamma][delta]*DiracAlgebra::Gamma0[delta][beta];
                             
                             
       
                            for(INT epsilon=0;epsilon<DiracAlgebra::SpinorComponents;epsilon++){
                                
                                // NON-STANDARD AXIAL TERMS //
                                // -i/2 GAMMA_0 GAMMMA_MU GAMMA_5 GAMMA_0
                                gAT2[alpha][beta]+=COMPLEX(0.0,-0.5)*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::Gamma0[gamma][delta]*DiracAlgebra::Gamma5[delta][epsilon]*DiracAlgebra::Gamma0[epsilon][beta];
                                gAX2[alpha][beta]+=COMPLEX(0.0,-0.5)*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaX[gamma][delta]*DiracAlgebra::Gamma5[delta][epsilon]*DiracAlgebra::Gamma0[epsilon][beta];
                                gAY2[alpha][beta]+=COMPLEX(0.0,-0.5)*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaY[gamma][delta]*DiracAlgebra::Gamma5[delta][epsilon]*DiracAlgebra::Gamma0[epsilon][beta];
                                gAZ2[alpha][beta]+=COMPLEX(0.0,-0.5)*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaZ[gamma][delta]*DiracAlgebra::Gamma5[delta][epsilon]*DiracAlgebra::Gamma0[epsilon][beta];
                                
                                
                            }
        
                            
                            
                            
                        }
                        
                        
                    }
                    
                }
                
            }
            
            
            
            ///////////////////////////
            // COMPUTE CONTRIBUTIONS //
            ///////////////////////////
            
            // RESET //
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        for(INT mu=0;mu<DiracAlgebra::SpinorComponents;mu++){
                            
                            jV[mu+4*Psi->Index3D(x,y,z)]=0.0;
                            jA[mu+4*Psi->Index3D(x,y,z)]=0.0;
                            
                        }
                    }
                }
            }
            
            
            // SETUP EVOLUTION OPERATOR //
            EvolutionOperator::Setup(U);
            
            // CALCULATE EIGEVANLUES VIA RITZ
            Dynamics::Fermions::ComputeRitzEigenvalues();
            
            // BUFFER FOR U^evo_ov Psi //
            wcf *UEvoOverlapPsi=new wcf[Lattice::N[0]*Lattice::N[1]*Lattice::N[2]];

            // COMPUTE CONTRIBUTION FROM EACH MODE
            for(INT s=0;s<Psi->NumberOfSamples;s++){
                
                // COMPUTE H_ov Psi FOR THIS MODE //
                EvolutionOperator::OverlapOperator::Compute(EvolutionOperator::Ucc,Psi->GetMode(s),UEvoOverlapPsi,Dynamics::Fermions::MinEValue,Dynamics::Fermions::MaxEValue);
                
                // COMPUTE AT EACH LOCATION //
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                            
                            // VECTOR AND AXIAL CURRENTS //
                            COMPLEX j0v=0.0; COMPLEX j1v=0.0; COMPLEX j2v=0.0; COMPLEX j3v=0.0;
                            COMPLEX j0a=0.0; COMPLEX j1a=0.0; COMPLEX j2a=0.0; COMPLEX j3a=0.0;
                            
                            //////////////////////
                            // COMPUTE CURRENTS //
                            //////////////////////
                            
                            for(INT i=0;i<Nc;i++){
                                for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){
                                    for(INT beta=0;beta<DiracAlgebra::SpinorComponents;beta++){
                                        
                                        
                                        j0v+=gVT[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s)*Psi->CommutatorExpectationValue(s);
                                        j0a+=gAT[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s)*Psi->CommutatorExpectationValue(s);
                                        
                                        j1v+=gVX[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s)*Psi->CommutatorExpectationValue(s);
                                        j1a+=gAX[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s)*Psi->CommutatorExpectationValue(s);
                                        
                                        j2v+=gVY[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s)*Psi->CommutatorExpectationValue(s);
                                        j2a+=gAY[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s)*Psi->CommutatorExpectationValue(s);
                                        
                                        j3v+=gVZ[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s)*Psi->CommutatorExpectationValue(s);
                                        j3a+=gAZ[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s)*Psi->CommutatorExpectationValue(s);
                                        
                                        j0v+=gVT2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta])*Psi->CommutatorExpectationValue(s);
                                        j0a+=gAT2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta])*Psi->CommutatorExpectationValue(s);
                                        
                                        j1v+=gVX2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta])*Psi->CommutatorExpectationValue(s);
                                        j1a+=gAX2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta])*Psi->CommutatorExpectationValue(s);
                                        
                                        j2v+=gVY2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta])*Psi->CommutatorExpectationValue(s);
                                        j2a+=gAY2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta])*Psi->CommutatorExpectationValue(s);
                                        
                                        j3v+=gVZ2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta])*Psi->CommutatorExpectationValue(s);
                                        j3a+=gAZ2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta])*Psi->CommutatorExpectationValue(s);
                                        
                                        /*
                                         // WITHOUT COMMUTATOR FOR SAYANTAN
                                        j0v+=gVT[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s);
                                        j0a+=gAT[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s);
                                        
                                        j1v+=gVX[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s);
                                        j1a+=gAX[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s);
                                        
                                        j2v+=gVY[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s);
                                        j2a+=gAY[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s);
                                        
                                        j3v+=gVZ[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s);
                                        j3a+=gAZ[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s);
                                      
                                        
                                        j0v+=gVT2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta]);
                                        j0a+=gAT2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta]);
                                        
                                        j1v+=gVX2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta]);
                                        j1a+=gAX2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta]);
                                        
                                        j2v+=gVY2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta]);
                                        j2a+=gAY2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta]);
                                        
                                        j3v+=gVZ2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta]);
                                        j3a+=gAZ2[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*(UEvoOverlapPsi[Psi->Index3D(x,y,z)].c[i].f[beta]);
                                        */
                                        
                                    }
                                }
                                
                            }
                            
                            
                            // SET PRE-FACTOR //
                            DOUBLE PreFactor=1.0/(2.0*Psi->N[0]*Psi->N[1]*Psi->N[2]);
                            

                            // WRITE VALUES TO GLOBAL ARRAY //
                            jV[0+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j0v);
                            jV[1+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j1v);
                            jV[2+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j2v);
                            jV[3+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j3v);
                            jA[0+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j0a);
                            jA[1+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j1a);
                            jA[2+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j2a);
                            jA[3+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j3a);
                            
                            
                        }
                    }
                }
                
            }
            
            
            // CLEAN-UP //
            delete UEvoOverlapPsi;
            
            // SYNCHRONIZE ALL MPI NODES //
            MPI_Barrier(MPI_COMM_WORLD);
            
            
            // PERFORM REDUCTION //
            DOUBLE *jVGlob=new DOUBLE[4*Psi->Volume];
            DOUBLE *jAGlob=new DOUBLE[4*Psi->Volume];
            
            MPI_Reduce(jV,jVGlob,4*Psi->Volume,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            MPI_Reduce(jA,jAGlob,4*Psi->Volume,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            
            
            // OUTPUT //
            if(MPIBasic::ID==0){
                
                std::ofstream VectorOutStream,AxialOutStream;
                
                VectorOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"VectorCurrentT",Dynamics::Time(),".txt").c_str());
                AxialOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"AxialCurrentT",Dynamics::Time(),".txt").c_str());
                
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                            VectorOutStream << Dynamics::Time() << " " << x << " " << y <<  " " << z << " " << jVGlob[0+4*Psi->Index3D(x,y,z)] << " " << jVGlob[1+4*Psi->Index3D(x,y,z)] << " " << jVGlob[2+4*Psi->Index3D(x,y,z)] << " " << jVGlob[3+4*Psi->Index3D(x,y,z)] << std::endl;
                            
                            AxialOutStream <<  Dynamics::Time() << " " << x << " " << y <<  " " << z << " " << jAGlob[0+4*Psi->Index3D(x,y,z)] << " " << jAGlob[1+4*Psi->Index3D(x,y,z)] << " " << jAGlob[2+4*Psi->Index3D(x,y,z)] << " " << jAGlob[3+4*Psi->Index3D(x,y,z)] << std::endl;
                            
                            
                        }
                        
                        VectorOutStream << std::endl;
                        AxialOutStream << std::endl;
                    }
                    
                    VectorOutStream << std::endl;
                    AxialOutStream << std::endl;
                }
                
                
                // DEBUG OUTPUT //
                std::cout << Dynamics::Time() << " " << 0 << " " << jVGlob[0+4*Psi->Index3D(1,3,2)] << "  " << jVGlob[1+4*Psi->Index3D(1,3,2)] << " " << jVGlob[2+4*Psi->Index3D(1,3,2)] << " " << jVGlob[3+4*Psi->Index3D(1,3,2)] << std::endl;
                
                std::cout << Dynamics::Time() << " " << 1 << " " << jAGlob[0+4*Psi->Index3D(1,3,2)] << "  " << jAGlob[1+4*Psi->Index3D(1,3,2)] << " " << jAGlob[2+4*Psi->Index3D(1,3,2)] << " " << jAGlob[3+4*Psi->Index3D(1,3,2)] << std::endl;
                
            }
            
            
            // CLEAN-UP //
            delete jV; delete jA;
            delete jVGlob; delete jAGlob;
            
        }
        
        void ComputeCurrents(GaugeLinks *U,FermionField *Psi){
            
            ComputeCurrents(0,Psi->N[0]-1,0,Psi->N[1]-1,0,Psi->N[2]-1,U,Psi);
            
        }
        
    }
    
    
}