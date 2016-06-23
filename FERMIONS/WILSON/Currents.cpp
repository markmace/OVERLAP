#include "LoopMacros.cpp"

namespace Fermions{
    
    namespace Observables{
        
        void ComputeCurrents(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,FermionField *Psi){
            
            // VECTOR AND AXIAL CURRENTS //
            DOUBLE *jV=new DOUBLE[4*Psi->Volume];   
            DOUBLE *jA=new DOUBLE[4*Psi->Volume];
            
            DOUBLE *jVw=new DOUBLE[4*Psi->Volume];
            DOUBLE *jAw=new DOUBLE[4*Psi->Volume];
            
            DOUBLE *PS=new DOUBLE[Psi->Volume];
            DOUBLE *WA=new DOUBLE[Psi->Volume];
            
            ///////////////////////////////////
            // SET LATTICE CURRENT COUPLINGS //
            ///////////////////////////////////
            
            COMPLEX gVT[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gVX[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gVY[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gVZ[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            
            COMPLEX gVwX[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gVwY[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gVwZ[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            
            COMPLEX gAT[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gAX[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gAY[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gAZ[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            
            COMPLEX gAwX[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gAwY[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            COMPLEX gAwZ[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            
            COMPLEX gPS[DiracAlgebra::SpinorComponents][DiracAlgebra::SpinorComponents];
            
            
            for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){
                for(INT beta=0;beta<DiracAlgebra::SpinorComponents;beta++){
                    
                    // RESET //
                    gVT[alpha][beta]=0.0; gVX[alpha][beta]=0.0; gVY[alpha][beta]=0.0; gVZ[alpha][beta]=0.0;
                    gAT[alpha][beta]=0.0; gAX[alpha][beta]=0.0; gAY[alpha][beta]=0.0; gAZ[alpha][beta]=0.0;
                    
                    gVwX[alpha][beta]=0.0; gVwY[alpha][beta]=0.0; gVwZ[alpha][beta]=0.0;
                    gAwX[alpha][beta]=0.0; gAwY[alpha][beta]=0.0; gAwZ[alpha][beta]=0.0;
                    
                    gPS[alpha][beta]=0.0;
                    
                    for(INT gamma=0;gamma<DiracAlgebra::SpinorComponents;gamma++){
                        
                        // STANDARD VECTOR TERMS //
                        gVT[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::Gamma0[gamma][beta];
                        gVX[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaX[gamma][beta];
                        gVY[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaY[gamma][beta];
                        gVZ[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaZ[gamma][beta];
                        
                        // PSEUDOSCALAR DENSITY TERMS //
                        gPS[alpha][beta]+=ComplexI*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::Gamma5[gamma][beta];
                        
                        for(INT delta=0;delta<DiracAlgebra::SpinorComponents;delta++){
                            
                            // STANDARD AXIAL TERMS //
                            gAT[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::Gamma0[gamma][delta]*DiracAlgebra::Gamma5[delta][beta];
                            gAX[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaX[gamma][delta]*DiracAlgebra::Gamma5[delta][beta];
                            gAY[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaY[gamma][delta]*DiracAlgebra::Gamma5[delta][beta];
                            gAZ[alpha][beta]+=DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::GammaZ[gamma][delta]*DiracAlgebra::Gamma5[delta][beta];
                            
                        }
                        
                        // AXIAL WILSON TERMS //
                        gAwX[alpha][beta]-=ComplexI*rWilson*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::Gamma5[gamma][beta];
                        gAwY[alpha][beta]-=ComplexI*rWilson*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::Gamma5[gamma][beta];
                        gAwZ[alpha][beta]-=ComplexI*rWilson*DiracAlgebra::Gamma0[alpha][gamma]*DiracAlgebra::Gamma5[gamma][beta];
                        
                    }
                    
                    // VECTOR WILSON TERMS //
                    gVwX[alpha][beta]-=ComplexI*rWilson*DiracAlgebra::Gamma0[alpha][beta];
                    gVwY[alpha][beta]-=ComplexI*rWilson*DiracAlgebra::Gamma0[alpha][beta];
                    gVwZ[alpha][beta]-=ComplexI*rWilson*DiracAlgebra::Gamma0[alpha][beta];
                    
                    
                }
                
            }
            
            
            
            ///////////////////////////
            // COMPUTE CONTRIBUTIONS //
            ///////////////////////////
            
            // RESET //
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        PS[Psi->Index3D(x,y,z)]=0.0;
                        WA[Psi->Index3D(x,y,z)]=0.0;
                        
                        for(INT mu=0;mu<4;mu++){
                            jV[mu+4*Psi->Index3D(x,y,z)]=0.0;
                            jA[mu+4*Psi->Index3D(x,y,z)]=0.0;
                            jVw[mu+4*Psi->Index3D(x,y,z)]=0.0;
                            jAw[mu+4*Psi->Index3D(x,y,z)]=0.0;
                            
                        }
                        
                        
                    }
                }
            }
            
            
            // COMPUTE CONTRIBUTION FROM EACH MODE
            for(INT s=0;s<Psi->NumberOfSamples;s++){
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                            // VECTOR AND AXIAL CURRENTS //
                            COMPLEX j0v=0.0; COMPLEX j1v=0.0; COMPLEX j2v=0.0; COMPLEX j3v=0.0;
                            COMPLEX j0a=0.0; COMPLEX j1a=0.0; COMPLEX j2a=0.0; COMPLEX j3a=0.0;
                            
                            // WILSON TERMS //
                            COMPLEX w1v=0.0; COMPLEX w2v=0.0; COMPLEX w3v=0.0;
                            COMPLEX w1a=0.0; COMPLEX w2a=0.0; COMPLEX w3a=0.0;
                            
                            // WILSON ANOMALY TERM //
                            COMPLEX WilsonAnomaly=0.0;
                            
                            // PSUEDOSCALE DENSITY //
                            COMPLEX psd=0.0;
                            
                            // SU(N) MATRIX BUFFERS //
                            SU_Nc_MATRIX_FORMAT Ux[Nc*Nc];
                            SU_Nc_MATRIX_FORMAT Uy[Nc*Nc];
                            SU_Nc_MATRIX_FORMAT Uz[Nc*Nc];
                            
                            SU_Nc_MATRIX_FORMAT UDxM[Nc*Nc];
                            SU_Nc_MATRIX_FORMAT UDyM[Nc*Nc];
                            SU_Nc_MATRIX_FORMAT UDzM[Nc*Nc];
                            
                            // PRE-FETCH GAUGE LINKS //
                            SUNcGroup::Operations::GetMatrix(U->Get(x,y,z,0),Ux);
                            SUNcGroup::Operations::GetMatrix(U->Get(x,y,z,1),Uy);
                            SUNcGroup::Operations::GetMatrix(U->Get(x,y,z,2),Uz);
                            
                            SUNcGroup::Operations::GetInverseMatrix(U->Get(x-1,y,z,0),UDxM);
                            SUNcGroup::Operations::GetInverseMatrix(U->Get(x,y-1,z,1),UDyM);
                            SUNcGroup::Operations::GetInverseMatrix(U->Get(x,y,z-1,2),UDzM);
                            
                            //////////////////////
                            // COMPUTE CURRENTS //
                            //////////////////////
                            
                            // j0v //
                            START_LOCAL_FERMION_LOOP(gVT)
                            j0v+=gVT[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s)*Psi->CommutatorExpectationValue(s);
                            END_LOCAL_FERMION_LOOP
                            
                            // j0a //
                            START_LOCAL_FERMION_LOOP(gAT)
                            j0a+=gAT[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s)*Psi->CommutatorExpectationValue(s);
                            END_LOCAL_FERMION_LOOP
                            
                            // pseudo-scalar //
                            START_LOCAL_FERMION_LOOP(gPS)
                            psd+=gPS[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s)*Psi->CommutatorExpectationValue(s);
                            END_LOCAL_FERMION_LOOP
                            
                            
                            // j1v //
                            START_FERMION_LOOP(gVX)
                            j1v+=gVX[alpha][beta]*Ux[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x+1,y,z,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            END_FERMION_LOOP
                            
                            // j1a //
                            START_FERMION_LOOP(gAX)
                            j1a+=gAX[alpha][beta]*Ux[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x+1,y,z,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            END_FERMION_LOOP
                            
                            // j2v //
                            START_FERMION_LOOP(gVY)
                            j2v+=gVY[alpha][beta]*Uy[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y+1,z,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            END_FERMION_LOOP
                            
                            // j2a //
                            START_FERMION_LOOP(gAY)
                            j2a+=gAY[alpha][beta]*Uy[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y+1,z,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            END_FERMION_LOOP
                            
                            // j3v //
                            START_FERMION_LOOP(gVZ)
                            j3v+=gVZ[alpha][beta]*Uz[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z+1,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            END_FERMION_LOOP
                            
                            // j3a //
                            START_FERMION_LOOP(gAZ)
                            j3a+=gAZ[alpha][beta]*Uz[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z+1,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            END_FERMION_LOOP
                            
                            //////////////////////////
                            // COMPUTE WILSON TERMS //
                            //////////////////////////
                            
                            
                            // LOCAL TERMS //
                            
                            START_LOCAL_FERMION_LOOP(gAwX)
                            WilsonAnomaly-=2.0*gAwX[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s)*Psi->CommutatorExpectationValue(s);
                            END_LOCAL_FERMION_LOOP
                            
                            START_LOCAL_FERMION_LOOP(gAwY)
                            WilsonAnomaly-=2.0*gAwY[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s)*Psi->CommutatorExpectationValue(s);
                            END_LOCAL_FERMION_LOOP
                            
                            START_LOCAL_FERMION_LOOP(gAwZ)
                            WilsonAnomaly-=2.0*gAwZ[alpha][beta]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z,beta,i,s)*Psi->CommutatorExpectationValue(s);
                            END_LOCAL_FERMION_LOOP
                            
                            
                            
                            // UPWARD DERIVATIVES //
                            
                            START_FERMION_LOOP(gVwX)
                            w1v+=gVwX[alpha][beta]*Ux[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x+1,y,z,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            END_FERMION_LOOP
                            
                            START_FERMION_LOOP(gAwX)
                            
                            w1a+=gAwX[alpha][beta]*Ux[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x+1,y,z,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            
                            WilsonAnomaly+=gAwX[alpha][beta]*Ux[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x+1,y,z,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            
                            END_FERMION_LOOP
                            
                            START_FERMION_LOOP(gVwY)
                            w2v+=gVwY[alpha][beta]*Uy[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y+1,z,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            END_FERMION_LOOP
                            
                            START_FERMION_LOOP(gAwY)
                            
                            w2a+=gAwY[alpha][beta]*Uy[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y+1,z,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            
                            WilsonAnomaly+=gAwY[alpha][beta]*Uy[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y+1,z,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            
                            END_FERMION_LOOP
                            
                            START_FERMION_LOOP(gVwZ)
                            w3v+=gVwZ[alpha][beta]*Uz[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z+1,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            END_FERMION_LOOP
                            
                            START_FERMION_LOOP(gAwZ)
                            
                            w3a+=gAwZ[alpha][beta]*Uz[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z+1,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            
                            WilsonAnomaly+=gAwZ[alpha][beta]*Uz[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z+1,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            
                            END_FERMION_LOOP
                            
                            
                            
                            // DOWNWARD DERIVATIVES //
                            
                            START_FERMION_LOOP(gAwX)
                            WilsonAnomaly+=gAwX[alpha][beta]*UDxM[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x-1,y,z,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            END_FERMION_LOOP
                            
                            START_FERMION_LOOP(gAwY)
                            WilsonAnomaly+=gAwY[alpha][beta]*UDyM[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y-1,z,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            END_FERMION_LOOP
                            
                            START_FERMION_LOOP(gAwZ)
                            WilsonAnomaly+=gAwZ[alpha][beta]*UDzM[i+Nc*j]*conj(Psi->GetValue(x,y,z,alpha,i,s))*Psi->GetValue(x,y,z-1,beta,j,s)*Psi->CommutatorExpectationValue(s);
                            END_FERMION_LOOP
                            
                            // -- //
                            
                            
                            // SET PRE-FACTOR //
                            DOUBLE PreFactor=1.0/(2.0*Psi->N[0]*Psi->N[1]*Psi->N[2]);
                            
                            // WRITE VALUES TO GLOBAL ARRAY //
                            jV[0+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j0v); jV[1+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j1v);  jV[2+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j2v);  jV[3+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j3v); 
                            jA[0+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j0a); jA[1+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j1a);  jA[2+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j2a);  jA[3+4*Psi->Index3D(x,y,z)]+=PreFactor*real(j3a);
                            
                            jVw[0+4*Psi->Index3D(x,y,z)]=0.0; jVw[1+4*Psi->Index3D(x,y,z)]+=PreFactor*real(w1v);  jVw[2+4*Psi->Index3D(x,y,z)]+=PreFactor*real(w2v);  jVw[3+4*Psi->Index3D(x,y,z)]+=PreFactor*real(w3v);
                            
                            jAw[0+4*Psi->Index3D(x,y,z)]=0.0; jAw[1+4*Psi->Index3D(x,y,z)]+=PreFactor*real(w1a);  jAw[2+4*Psi->Index3D(x,y,z)]+=PreFactor*real(w2a);  jAw[3+4*Psi->Index3D(x,y,z)]+=PreFactor*real(w3a);
                            
                            PS[Psi->Index3D(x,y,z)]+=PreFactor*real(psd);
                            
                            WA[Psi->Index3D(x,y,z)]+=PreFactor*real(WilsonAnomaly);
                            
                            
                            
                        }
                    }
                }
                
            }
            
            // SYNCHRONIZE ALL MPI NODES //
            MPI_Barrier(MPI_COMM_WORLD);
            
            
            // PERFORM REDUCTION //
            DOUBLE *jVGlob=new DOUBLE[4*Psi->Volume];   
            DOUBLE *jAGlob=new DOUBLE[4*Psi->Volume];
            
            DOUBLE *jVwGlob=new DOUBLE[4*Psi->Volume];
            DOUBLE *jAwGlob=new DOUBLE[4*Psi->Volume];
            
            DOUBLE *PSGlob=new DOUBLE[Psi->Volume];
            DOUBLE *WAGlob=new DOUBLE[Psi->Volume];
            
            
            MPI_Reduce(jV,jVGlob,4*Psi->Volume,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            MPI_Reduce(jA,jAGlob,4*Psi->Volume,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            
            MPI_Reduce(jVw,jVwGlob,4*Psi->Volume,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            MPI_Reduce(jAw,jAwGlob,4*Psi->Volume,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            
            MPI_Reduce(PS,PSGlob,Psi->Volume,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            MPI_Reduce(WA,WAGlob,Psi->Volume,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            
            
            
            // OUTPUT //
            if(MPIBasic::ID==0){
                
                std::ofstream VectorOutStream,AxialOutStream,PseudoScalarOutStream,WilsonAnomalyOutStream;
                
                VectorOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"VectorCurrentT",Dynamics::Time(),".txt").c_str());
                AxialOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"AxialCurrentT",Dynamics::Time(),".txt").c_str());
                PseudoScalarOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"PseudoScalarT",Dynamics::Time(),".txt").c_str());
                WilsonAnomalyOutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"WilsonAnomalyT",Dynamics::Time(),".txt").c_str());
                
                for(INT z=zLow;z<=zHigh;z++){
                    for(INT y=yLow;y<=yHigh;y++){
                        for(INT x=xLow;x<=xHigh;x++){
                            
                            VectorOutStream << Dynamics::Time() << " " << x << " " << y <<  " " << z << " " << jVGlob[0+4*Psi->Index3D(x,y,z)] << " " << jVGlob[1+4*Psi->Index3D(x,y,z)] << " " << jVGlob[2+4*Psi->Index3D(x,y,z)] << " " << jVGlob[3+4*Psi->Index3D(x,y,z)] << " " << jVwGlob[0+4*Psi->Index3D(x,y,z)] << " " << jVwGlob[1+4*Psi->Index3D(x,y,z)] << " " << jVwGlob[2+4*Psi->Index3D(x,y,z)] << " " << jVwGlob[3+4*Psi->Index3D(x,y,z)] << std::endl;
                            
                            AxialOutStream <<  Dynamics::Time() << " " << x << " " << y <<  " " << z << " " << jAGlob[0+4*Psi->Index3D(x,y,z)] << " " << jAGlob[1+4*Psi->Index3D(x,y,z)] << " " << jAGlob[2+4*Psi->Index3D(x,y,z)] << " " << jAGlob[3+4*Psi->Index3D(x,y,z)] << " " << jAwGlob[0+4*Psi->Index3D(x,y,z)] << " " << jAwGlob[1+4*Psi->Index3D(x,y,z)] << " " << jAwGlob[2+4*Psi->Index3D(x,y,z)] << " " << jAwGlob[3+4*Psi->Index3D(x,y,z)] << std::endl;
                            
                            PseudoScalarOutStream << Dynamics::Time() << " " << x << " " << y <<  " " << z << " " << PSGlob[Psi->Index3D(x,y,z)] << std::endl;
                            
                            WilsonAnomalyOutStream << Dynamics::Time() << " " << x << " " << y << " " << z << " " << WAGlob[Psi->Index3D(x,y,z)] << std::endl;
                            
                            
                        }
                        
                        VectorOutStream << std::endl;
                        AxialOutStream << std::endl;
                        PseudoScalarOutStream << std::endl;
                        WilsonAnomalyOutStream << std::endl;
                    }
                    
                    VectorOutStream << std::endl;
                    AxialOutStream << std::endl;
                    PseudoScalarOutStream << std::endl;
                    WilsonAnomalyOutStream << std::endl;
                }
                
                
                // DEBUG OUTPUT //
                std::cout << Dynamics::Time() << " " << 0 << " " << jVGlob[0+4*Psi->Index3D(1,3,2)] << "  " << jVGlob[1+4*Psi->Index3D(1,3,2)] << " " << jVGlob[2+4*Psi->Index3D(1,3,2)] << " " << jVGlob[3+4*Psi->Index3D(1,3,2)] << std::endl;
                
                std::cout << Dynamics::Time() << " " << 10 << " " << jVwGlob[0+4*Psi->Index3D(1,3,2)] << "  " << jVwGlob[1+4*Psi->Index3D(1,3,2)] << " " << jVwGlob[2+4*Psi->Index3D(1,3,2)] << " " << jVwGlob[3+4*Psi->Index3D(1,3,2)] << std::endl;
                
                std::cout << Dynamics::Time() << " " << 1 << " " << jAGlob[0+4*Psi->Index3D(1,3,2)] << "  " << jAGlob[1+4*Psi->Index3D(1,3,2)] << " " << jAGlob[2+4*Psi->Index3D(1,3,2)] << " " << jAGlob[3+4*Psi->Index3D(1,3,2)] << std::endl;
                
                std::cout << Dynamics::Time() << " " << 11 << " " << jAwGlob[0+4*Psi->Index3D(1,3,2)] << "  " << jAwGlob[1+4*Psi->Index3D(1,3,2)] << " " << jAwGlob[2+4*Psi->Index3D(1,3,2)] << " " << jAwGlob[3+4*Psi->Index3D(1,3,2)] << std::endl;
                
                std::cout << Dynamics::Time() << " " << 2 << " " << PSGlob[Psi->Index3D(1,3,2)] << std::endl;
                
                std::cout << Dynamics::Time() << " " << 12 << " " << WAGlob[Psi->Index3D(1,3,2)] << std::endl;
                
            }
            
            // CLEAN-UP //
            delete jV; delete jA; delete jVw; delete jAw; delete PS; delete WA;
            delete jVGlob; delete jAGlob; delete jVwGlob; delete jAwGlob; delete PSGlob; delete WAGlob;
            
        }
        
        void ComputeCurrents(GaugeLinks *U,FermionField *Psi){
            
            ComputeCurrents(0,Psi->N[0]-1,0,Psi->N[1]-1,0,Psi->N[2]-1,U,Psi);
            
        }
        
    }
    
    
}