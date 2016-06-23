#ifndef OVERLAP_QED_EIGENFUNCTIONS_CPP
#define OVERLAP_QED_EIGENFUNCTIONS_CPP

#define PERIODICDELTA(x,y,N) (((x)==(y) || (x)==(y)+(N) || (x)==(y)-(N))?1.0:0.0)

#include <iomanip>

namespace OverlapQEDEigenfunctions{
    
    // EIGENFUNCTIONS IN THE PRESENCE OF CONSTANT MAGNETIC FIELD //
    Eigensystem **EigenModes;
    
    int VectorIndex(INT x,INT y,INT z,INT alpha, INT i){
        
        return GLinks::U->Index3D(x,y,z)+Lattice::N[0]*Lattice::N[1]*Lattice::N[2]*(alpha+i*DiracAlgebra::SpinorComponents);

    }
    
    // CONVERTS U(1) GAUGE LINKS INTO UCC FORMAT //
    void ConvertU1GaugeLinks(COMPLEX *U,u_cc *UccLinks){
        
        for(INT z=0;z<GLinks::U->N[2];z++){
            for(INT y=0;y<GLinks::U->N[1];y++){
                for(INT x=0;x<GLinks::U->N[0];x++){
                    
                    for(INT mu=0;mu<GLinks::U->Dimension;mu++){
                        
                        INT CCIndex=GLinks::U->Index3D(x,y,z)+mu*(GLinks::U->N[0])*(GLinks::U->N[1])*(GLinks::U->N[2]);
                        
                        // SET ALL TO ZERO //
                        for(INT i=0;i<Nc;i++){
                            for(INT j=0;j<Nc;j++){
                                
                                // CHECK CONVENTION HERE //
                                UccLinks[CCIndex].c[i].c[j]=0.0;
                                
                            }
                        }
                        
                        // SET VALUE OF QED FIELD TO UccLink
                        UccLinks[CCIndex].c[0].c[0]=QED::GetUVal(x,y,z,mu);
                        
                    }
                    
                }
            }
        }
        
    }
    
    // NORMALIZE A SINGLE VECTOR WITH START AND ENDING INDEX
    void Normalize(COMPLEX *V,INT size,INT nMin,INT nMax){
        
        DOUBLE total=0.0;
        
        for(INT i=nMin;i<nMax;i++){
            total+=norm(V[i]);
        }
        
        for(INT i=nMin;i<nMax;i++){
            V[i]/=sqrt(total);
        }
        
    }
    
    // OVERLOAD
    // NORMALIZE A SINGLE VECTOR
    void Normalize(COMPLEX *V,INT size){
        
        Normalize(V,size,0,size);
        
    }
    
    // ORTHOGONALIZE v2 WITH RESPECT TO v1 //
    void Orthogonalize(COMPLEX *v1,COMPLEX *v2,INT Length,INT nEigen,INT nMin, INT nMax){
        
        COMPLEX val=0.0;
        
        for(INT i=0;i<nEigen-1;i++){
            
            val=0.0;
            for(INT j=nMin;j<nMax;j++){
                val+=v2[j]*conj(v1[j-nMin+Length*i]);
            }
            
            for(INT j=nMin;j<nMax;j++){
                v2[j]-=v1[j-nMin+Length*i]*val;
            }
        }
        
    }
    
    // CONVERTS GAUGE LINKS INTO UCC FORMAT //
    void ConvertBackGaugeLinks(u_cc *UccLinks,GaugeLinks *U){
        
        for(INT z=0;z<GLinks::U->N[2];z++){
            for(INT y=0;y<GLinks::U->N[1];y++){
                for(INT x=0;x<GLinks::U->N[0];x++){
                    
                    for(INT mu=0;mu<GLinks::U->Dimension;mu++){
                        
                        INT CCIndex=GLinks::U->Index3D(x,y,z)+mu*(U->N[0])*(U->N[1])*(U->N[2]);
                        
                        SU_Nc_MATRIX_FORMAT UMat[Nc*Nc];
                        
                        
                        for(INT i=0;i<Nc;i++){
                            for(INT j=0;j<Nc;j++){
                                
                                // CHECK CONVENTION HERE //
                                UMat[i+j*Nc]=UccLinks[CCIndex].c[i].c[j];
                                
                            }
                        }
                        
                        COPY_SUNcMatrix(UMat,U->Get(x,y,z,mu));
                        
                    }
                    
                }
            }
        }
        
    }
    
    // NOTE -- DEFINING ALL VECTORS AS COLUMNS //
    // ALL MATRICES ARE DEFINED AS FOLLOWS     //
    ///////////////////////////////////////////
    //                                       //
    //  A_0      A_n      ...       A_n(n-1) //
    //  A_1                         A_2n-1   //
    //  .                             .      //
    //  .                             .      //
    //  .                             .      //
    //  A_(n-1)  A_(2n-1) ...       A_n^2-1  //
    //                                       //
    ///////////////////////////////////////////
    
    
    // COMPUTE EIGENFUNCTIONS //
    void Setup(COMPLEX *U_A,INT INPUT_FLAG,INT OUTPUT_FLAG){
        
        // DIMENSION OF OVERLAP OPERATOR //
        INT OperatorDimension=Lattice::N[0]*Lattice::N[1]*Lattice::N[2]*Nc*DiracAlgebra::SpinorComponents;
        // STANDARD FREE GAUGE LINKS FOR DEBUGGING //
        //EvolutionOperator::ConvertGaugeLinks(GLinks::U,EvolutionOperator::Ucc);
        
        // SET U(1) GAUGE LINKS AS INITIAL UCC LINKS
        ConvertU1GaugeLinks(U_A,EvolutionOperator::Ucc);
        
        // OUTPUT GAUGE LINKS
        /*
        if(MPIBasic::ID==0){
        
            // COMMANDLINE OUTPUT //
            std::cerr << "#GAUGE LINKS SAVED AT T=" << Dynamics::Time() << std::endl;
            
            // OUTPUT STREAMS //
            std::ofstream UOutStream;
            
            // OUTPUT FILES //
            std::string UOutputFile=StringManipulation::StringCast(IO::OutputDirectory,"AnotherUOutT",Dynamics::Time(),"ID",RandomNumberGenerator::MySEED,".txt");
            
            // OPEN FILES //
            UOutStream.open(UOutputFile.c_str());
            
            // SET PRECISION //
            UOutStream.precision(OUTPUT_PRECISION);
            
            // CREATE OUTPUT //
            for(INT z=0;z<GLinks::U->N[2];z++){
                for(INT y=0;y<GLinks::U->N[1];y++){
                    for(INT x=0;x<GLinks::U->N[0];x++){
                        for(INT mu=0;mu<Lattice::Dimension;mu++){
                            
                            INT CCIndex=GLinks::U->Index3D(x,y,z)+mu*(GLinks::U->N[0])*(GLinks::U->N[1])*(GLinks::U->N[2]);
                            
                            SU_Nc_MATRIX_FORMAT UMat[Nc*Nc];
                            
                            for(INT i=0;i<Nc;i++){
                                for(INT j=0;j<Nc;j++){
                                    
                                    // CHECK CONVENTION HERE //
                                    UMat[i+j*Nc]=EvolutionOperator::Ucc[CCIndex].c[i].c[j];
                                    
                                }
                            }
                            
                            DOUBLE Uo[Nc*Nc];
                            
                            Uo[0]=real(COMPLEX(0.0,0.5)*(-UMat[1]-UMat[2]));
                            Uo[1]=real(COMPLEX(0.5,0.0)*(-UMat[1]+UMat[2]));
                            Uo[2]=real(COMPLEX(0.0,0.5)*( UMat[0]-UMat[3]));
                            Uo[3]=real(COMPLEX(0.5,0.0)*( UMat[0]+UMat[3]));
                            

                            // OUTPUT GAUGE LINKS //
                            UOutStream << std::fixed << std::setprecision(9) << x << " " << y << " " << z << " " << mu << " " << Uo[0] << " " << Uo[1] <<  " " << Uo[2] <<  " " << Uo[3] << std::endl;
                            
                        }
                        
                    }
                }
            }
            
            // CLOSE OUTPUT STREAM //
            UOutStream.close();

        
        }
        */
        
        // COMMANDLINE OUTPUT //
        if(MPIBasic::ID==0){
            std::cerr << "#COMPUTING OVERLAP EIGENFUNCTIONS IN MAGNETIC FIELD" << std::endl;
        }
        
        // COMPUTE COMPLETE SET OF EIGENFUNCTIONS //
        EigenModes=new Eigensystem*[1];
        
        ////////////////////////////////////////////////////////////////////////
        // SET RANDOM ORTHONORMAL VECTOR SPACE FOR HAMILTONIAN RECONSTRUCTION //
        ////////////////////////////////////////////////////////////////////////
        
        // DETERMINE ELEMENTARY COMPLEX VECTOR SPACE //
        COMPLEX *RandomVectorSpace=new COMPLEX[OperatorDimension*OperatorDimension];
        
        // SET UP RANDOM COMPLEX ORTHONORMAL VECTOR SPACE //
        COMPLEX *InitialRandomVectors=new COMPLEX[OperatorDimension*OperatorDimension];
        COMPLEX *TempVector=new COMPLEX[OperatorDimension];
        
        // COMMANDLINE OUTPUT //
        if(MPIBasic::ID==0){
            std::cerr << "INITIALIZING RANDOM VECTOR SPACE" << std::endl;
        }
        
        for(INT j=0;j<OperatorDimension;j++){
            for(INT i=0;i<OperatorDimension;i++){
                // SET INITIAL VECTOR SPACE TO COMPLEX RANDOM //
                // VECTORS ARE SET STORED AS COLUMNS
                InitialRandomVectors[i*OperatorDimension+j]=COMPLEX(RandomNumberGenerator::rng(),RandomNumberGenerator::rng());
                // SET FINAL VECTOR SPACE TO ZERO //
                RandomVectorSpace[i+OperatorDimension*j]=COMPLEX(0.0);
                
            }
            
        }
        
        // ZERO-TH VECTOR CASE -- NORMALIZATION
        for(INT j=0;j<OperatorDimension;j++){
            RandomVectorSpace[j]=InitialRandomVectors[j];
        }
        
        Normalize(RandomVectorSpace,OperatorDimension);
        
        // 1,...,OperatorDimension-1 VECTOR CASE //
        for(INT j=1;j<OperatorDimension;j++){
            
            // SET COLUMN VECTOR TO BE ORTHOGONALIZED //
            for(INT i=0;i<OperatorDimension;i++){
                TempVector[i]=InitialRandomVectors[i+j*OperatorDimension];
            }
            
            // ORTHOGONALIZE WRT ENTIRE PREVIOUS VECTOR SPACE //
            Orthogonalize(RandomVectorSpace,TempVector,OperatorDimension,j+1,0,OperatorDimension);
            
            // GET NORM OF NEWLY ORTHOGONAL VECTOR //
            Normalize(TempVector,OperatorDimension);
            
            // SAVE ORTHONORMAL VECTOR AS COLUMN //
            for(INT i=0;i<OperatorDimension;i++){
                
                RandomVectorSpace[i+j*OperatorDimension]=TempVector[i];
            }
            
        }
        
        // COMMANDLINE OUTPUT //
        if(MPIBasic::ID==0){
            std::cerr << "RANDOM VECTOR SPACE SET" << std::endl;
        }
        // CLEAN-UP FOR RANDOM VECTOR SPACE BUFFERS
        delete[] InitialRandomVectors;
        delete[] TempVector;
        
        // OUTPUT INITIAL VECTOR SPACE TO FILE FOR DEBU //
        /*
        if(MPIBasic::ID==0){
            
            // COMMANDLINE OUTPUT //
            std::cerr << "#SAVING INITIAL VECTOR SPACE" << std::endl;
            
            // CREATE OUT-STREAM //
            std::ofstream OutStream;
            OutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"InitialVectorSpace.txt").c_str());
            
            OutStream.precision(OUTPUT_PRECISION);
            
            for(INT i=0;i<OperatorDimension;i++){
                for(INT j=0;j<OperatorDimension;j++){
                    
                    OutStream << std::fixed  << std::setprecision(9) << std::real(RandomVectorSpace[i+j*OperatorDimension]) << "+I*" << std::imag(RandomVectorSpace[i+j*OperatorDimension]) << " ";
                    
                }
                OutStream << std::endl;
                
                
            }
            
            // CLOSE OUT-STREAM //
            OutStream.close();
            
        }
         */
        
        ////////////////////////////////////////////////////////////////////
        // DETERMINE AND SOLVE EIGENVALUE PROBLEM FOR OVERLAP HAMILTONIAN //
        ////////////////////////////////////////////////////////////////////
        
        // DEFINE EIGENVALUE TEST VECTOR IN COMPLEX AND WCF FORMATS //
        COMPLEX *TestVector=new COMPLEX[OperatorDimension];
        wcf *TV=new wcf[GLinks::U->N[0]*GLinks::U->N[1]*GLinks::U->N[2]];
        
        // DEFINE OVERLAP TRANSFORMATION OUTPUT
        wcf *OverlapOutput=new wcf[GLinks::U->N[0]*GLinks::U->N[1]*GLinks::U->N[2]];
        
        // DEFINE VECTORS AFTER OVERLAP UPDATE
        COMPLEX *HamiltonianTransformedVectorSpace=new COMPLEX[OperatorDimension*OperatorDimension];
        
        // DETERMINE VECTORS AFTER OVERLAP UPDATE //
        for(INT s=0;s<OperatorDimension;s++){
            
            for(INT i=0;i<Lattice::N[0]*Lattice::N[1]*Lattice::N[2];i++){
                for(INT k=0;k<DiracAlgebra::SpinorComponents;k++){
                    for(INT j=0;j<Nc;j++){
                        
                        // INDEX -- SAME AS VECTOR INDEX //
                        INT m=i+Lattice::N[0]*Lattice::N[1]*Lattice::N[2]*(k+j*DiracAlgebra::SpinorComponents);
                        
                        // SET TEST VECTOR IN TWO FORMATS //
                        TestVector[m]=RandomVectorSpace[m+s*OperatorDimension];
                        TV[i].c[j].f[k]=RandomVectorSpace[m+s*OperatorDimension];
                        
                        
                    }
                }
            }
            
            // CALCULATE MAXIMUM EIGENVALUE OF WILSON HAMILTONIAN //
            DOUBLE eMax=-DOUBLE(1.05)*EvolutionOperator::OverlapOperations::RitzSolver(EvolutionOperator::Ucc,TestVector,Dynamics::Fermions::NumberOfEVals,Dynamics::Fermions::NormRatio,Dynamics::Fermions::MaximumIterations,-1.0,Dynamics::Fermions::NumberOfRestarts,OperatorDimension);
            
            // CALCULATE MINIMUM EIGENVALUE OF WILSON HAMILTONIAN //
            DOUBLE eMin=DOUBLE(0.95)*EvolutionOperator::OverlapOperations::RitzSolver(EvolutionOperator::Ucc,TestVector,Dynamics::Fermions::NumberOfEVals,Dynamics::Fermions::NormRatio,Dynamics::Fermions::MaximumIterations,1.0,Dynamics::Fermions::NumberOfRestarts,OperatorDimension);
            /*
            if(MPIBasic::ID==0){
                std::cerr << "#eMax eMin " << eMax << " " << eMin << std::endl;
            }
            */
            
            /////////////////////
            // SET HAMILTONIAN //
            /////////////////////
            
            // SET TIMER //
            if(MPIBasic::ID==0){
                std::cerr << "#COMPUTING OVERLAP DIRAC OPERATOR FOR s= " << s << std::endl;
                Timing::Reset();
            }
            
            
            // COMPUTE OVERLAP OPERATOR WITH RANDOM ORTHOGONAL TEST VECTOR //
            EvolutionOperator::OverlapOperator::Compute(EvolutionOperator::Ucc,TV,OverlapOutput,eMin,eMax);
            // SAVE OUTPUT //
            for(INT i=0;i<Lattice::N[0]*Lattice::N[1]*Lattice::N[2];i++){
                for(INT k=0;k<DiracAlgebra::SpinorComponents;k++){
                    for(INT j=0;j<Nc;j++){
                        
                        // INDEX //
                        INT m=i+Lattice::N[0]*Lattice::N[1]*Lattice::N[2]*(k+j*DiracAlgebra::SpinorComponents);
                        // SAVED AS COLUMNS
                        HamiltonianTransformedVectorSpace[m+s*OperatorDimension]=OverlapOutput[i].c[j].f[k];
                        
                    }
                }
            }
            
            // TIMING //
            if(MPIBasic::ID==0){
                std::cerr << "#TIMING -- " << Timing::Get() << std::endl;
                Timing::Reset();
            }
            
        }
        
        // OUTPUT FINAL VECTOR SPACE TO FILE FOR DEBUG
        /*
        if(MPIBasic::ID==0){
            
            // COMMANDLINE OUTPUT //
            std::cerr << "#SAVING FINAL VECTOR SPACE" << std::endl;
            
            // CREATE OUT-STREAM //
            std::ofstream OutStream;
            OutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"FinalVectorSpace.txt").c_str());
            
            OutStream.precision(OUTPUT_PRECISION);
            
            for(INT i=0;i<OperatorDimension;i++){
                for(INT j=0;j<OperatorDimension;j++){
                    
                    OutStream << std::fixed  << std::setprecision(9) << std::real(HamiltonianTransformedVectorSpace[i+j*OperatorDimension]) << "+I*" << std::imag(HamiltonianTransformedVectorSpace[i+j*OperatorDimension]) << " ";
                    
                }
                OutStream << std::endl;
                
                
            }
            
            // CLOSE OUT-STREAM //
            OutStream.close();
            
        }
        */
        // RECONSTRUCTED OPERLAP OPERATOR //
        COMPLEX *ReconstructedOverlapHamiltonian=new COMPLEX[OperatorDimension*OperatorDimension];
        
        // RECONSTRUCT THE UPDATE OPERATOR FROM [HamiltonianTransformedVectorSpace].[RandomVectorSpace]^\dagger //
        for(INT i=0;i<OperatorDimension;i++){
            for(INT j=0;j<OperatorDimension;j++){
                
                ReconstructedOverlapHamiltonian[i+j*OperatorDimension]=0.0;
                
                for(INT k=0;k<OperatorDimension;k++){
                    
                    ReconstructedOverlapHamiltonian[i+j*OperatorDimension]+=HamiltonianTransformedVectorSpace[i+OperatorDimension*k]*conj(RandomVectorSpace[k*OperatorDimension+j]);
                    
                }
                // MULTIPLY BY I FOR HAMILTONIAN //
                ReconstructedOverlapHamiltonian[i+j*OperatorDimension]*=COMPLEX(0.0,1.0);
            }
        }
        
        // OUTPUT HAMILTONIAN TO FILE FOR DEBUG
        /*
        if(MPIBasic::ID==0){
            
            // COMMANDLINE OUTPUT //
            std::cerr << "#SAVING HAMILTONIAN" << std::endl;
            
            // CREATE OUT-STREAM //
            std::ofstream OutStream;
            OutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"Hamiltonian.txt").c_str());
            
            OutStream.precision(OUTPUT_PRECISION);
        
            for(INT i=0;i<OperatorDimension;i++){
                for(INT j=0;j<OperatorDimension;j++){
                    
                    OutStream << std::fixed  << std::setprecision(9) << std::real(ReconstructedOverlapHamiltonian[i+j*OperatorDimension]) << "+I*" << std::imag(ReconstructedOverlapHamiltonian[i+j*OperatorDimension]) << " ";
                    
                }
                OutStream << std::endl;

                
            }
            
            // CLOSE OUT-STREAM //
            OutStream.close();
            
        }
         */

        // CHECK HERMITICITY //
        if(MPIBasic::ID==0){
            
            for(INT i=0;i<OperatorDimension;i++){
                for(INT j=0;j<OperatorDimension;j++){
                    
                    COMPLEX diff=ReconstructedOverlapHamiltonian[i+j*OperatorDimension]-conj(ReconstructedOverlapHamiltonian[j+i*OperatorDimension]);
                    if(norm(diff)>std::pow(10,-12)){
                        
                        std::cerr << "#DIFF -- NOT HERMITIAN " << i << " "  << j << " " << norm(diff) << std::endl;
                    }
                    
                }
            }
        }
        ///////////////////////
        // SETUP EIGENSYSTEM //
        ///////////////////////
        
        EigenModes[0]=new Eigensystem(ReconstructedOverlapHamiltonian,OperatorDimension);
        
        
        ////////////////////////////////////////////
        // COMPUTE EIGENSYSTEM BY DIAGONALIZATION //
        ////////////////////////////////////////////
        
        if(INPUT_FLAG==0){
            
            // COMMANDLINE OUTPUT //
            if(MPIBasic::ID==0){
                std::cerr << "#DIAGONALIZING OVERLAP HAMILTONIAN OPERATOR" << std::endl;
                Timing::Reset();
            }
            
            // DIAGONALIZE //
            EigenModes[0]->Diagonalize();
            
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
            std::cerr << "#CHECKING EIGENSYSTEM" << std::endl;
            Timing::Reset();
        }
        
        EigenModes[0]->Check();
        
        // TIMING //
        if(MPIBasic::ID==0){
            std::cerr << "#TIMING -- " << Timing::Get() << std::endl;
            Timing::Reset();
        }
         
        if(MPIBasic::ID==0){
            std::cerr << "#DONE SETTING OVERLAP EIGENFUNCTIONS IN MAGNETIC FIELD" << std::endl;
        }
        /*
        // OUTPUT EIGENSYSTEM
        if(MPIBasic::ID==0){
            
            // COMMANDLINE OUTPUT //
            std::cerr << "#SAVING EIGENSYSTEM" << std::endl;
            
            // CREATE OUT-STREAM //
            std::ofstream OutStream,OutStream2;
            OutStream.open(StringManipulation::StringCast(IO::OutputDirectory,"Eigenvectors.txt").c_str());
            OutStream2.open(StringManipulation::StringCast(IO::OutputDirectory,"Eigenvalues.txt").c_str());

            
            OutStream.precision(OUTPUT_PRECISION);
            OutStream2.precision(OUTPUT_PRECISION);

            
            DOUBLE eVal; COMPLEX *eVec=new COMPLEX[OperatorDimension];
            
            for(INT s=0;s<OperatorDimension;s++){
                
                EigenModes[0]->GetEigenMode(s,eVal,eVec);
                
                OutStream2 << std::fixed  << std::setprecision(9) << eVal << " ";
                
                for(INT j=0;j<OperatorDimension;j++){
                    
                    OutStream << std::fixed  << std::setprecision(9) << std::real(eVec[j]) << "+I*" << std::imag(eVec[j]) << " ";
                    
                }
                OutStream << std::endl;
                
                
                
                
            }
            // CLOSE OUT-STREAM //
            OutStream.close();
            OutStream2.close();

            
        }
         */
        
        // CLEAN-UP
        delete[] RandomVectorSpace;
        delete[] TestVector;
        delete[] TV;
        delete[] HamiltonianTransformedVectorSpace;
        
        delete[] OverlapOutput;
        delete[] ReconstructedOverlapHamiltonian;
    }
    
    
    
    
}

#endif