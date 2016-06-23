// CONVERSION FROM C SOURCE CODE
// NOTE FOR NERSC NEEDS TO BE cheev_ INSTEAD OF cheev
extern "C" {
    void zheev(char const* jobz,char const* uplo,INT* n,COMPLEX* a,INT* lda, DOUBLE* w, COMPLEX* work, int* lwork, DOUBLE* rwork,INT* info );
}

////////////////////////////////////////
//                                     //
// NOTE: LAPACK AND THIS CODE DEFINES MATRICES AS    //
//                                     //
//          A_0 A_3 A_6                //
//          A_1 A_4 A_7                //
//          A_2 A_5 A_8                //
//                                     //
//              OR                     //
//                                     //
//          A_00 A_01 A_02             //
//          A_10 A_11 A_12             //
//          A_20 A_21 A_22             //
//                                     //
//    WHEREAS GSL DEFINES IT AS        //
//                                     //
//          A_0 A_1 A_1                //
//          A_3 A_4 A_5                //
//          A_6 A_7 A_8                //
//                                     //
//              OR                     //
//                                     //
//          A_00 A_01 A_02             //
//          A_10 A_11 A_12             //
//          A_20 A_21 A_22             //
//                                     //
//                                     //
/////////////////////////////////////////



class Eigensystem {
    
private:
    
    // LINEAR DIMENSION OF THE MATRIX //
    INT MatrixSize;
    
private:
    
    // MATRIX //
    std::complex<double> *HermitianMatrix;
    
    std::complex<double> *HermitianMatrixUnTransposed;
    
    
    // EIGENVECTOR AND EIGENVALUES //
    DOUBLE *EigenValues;
    COMPLEX *EigenVectors;
    
public:
    
    // ACCESS EIGENMODES //
    void GetEigenMode(INT s,DOUBLE &eVal,COMPLEX *eVec){
        
        eVal=EigenValues[s];
        
        for(INT i=0;i<MatrixSize;i++){
            eVec[i]=EigenVectors[i+s*MatrixSize];
        }
        
    }
    
    // ACCESS EIGENMODES //
    void SetEigenMode(INT s,DOUBLE eVal,COMPLEX *eVec){
        
        EigenValues[s]=eVal;
        
        for(INT i=0;i<MatrixSize;i++){
            EigenVectors[i+s*MatrixSize]=eVec[i];
        }
        
    }
    
public:
    
    // DIAGONALIZATION USING LAPACK LIBRARY //
    void Diagonalize(){
        
        // COPY MATRIX //
        std::complex<double> *CopyMatrix=new std::complex<double>[MatrixSize*MatrixSize];
        std::memcpy(CopyMatrix,this->HermitianMatrix,MatrixSize*MatrixSize*sizeof(std::complex<double>));
        
        DOUBLE *eVals=new DOUBLE[MatrixSize];
        COMPLEX *eVecs=new COMPLEX[MatrixSize*MatrixSize];

        // LOCAL BUFFERS //
        INT info, lwork;
        COMPLEX wkopt;
        // rwork dimension should be at least max(1,3*n-2) //
        DOUBLE *rwork=new DOUBLE[3*MatrixSize-2];

        // QUERY AND ALLOCATE WORKSPACE //
        lwork = -1;
        zheev( "V", "L", &MatrixSize, CopyMatrix, &MatrixSize, eVals, &wkopt, &lwork, rwork, &info );
        lwork = int(real(wkopt));
        // DEFINE WORKSPACE
        COMPLEX* work= new COMPLEX[lwork*sizeof(COMPLEX)];
        // SOLVE EIGENPROBLEM //
        zheev( "V", "L", &MatrixSize, CopyMatrix, &MatrixSize, eVals, work, &lwork, rwork, &info );

        // CHECK CONVERGENCE //
        if( info > 0 ) {
            std::cerr << "The algorithm failed to compute eigenvalues." << std::endl;
            exit( 1 );
        }

        // SET EIGENSYSTEM //
        // IN REVERSE ORDER AS LAPACK GIVES LOWEST TO HIGHEST //
        /*for(int s=0;s<MatrixSize;s++){

            
            // SET EIGENVALUE //
            this->EigenValues[s] = eVals[s];
            
            std::cout << this->EigenValues[s] << std::endl;
            
            // SET EIGENVECTOR //
            for(int j=0;j<MatrixSize;j++){
                
                this->EigenVectors[j+s*MatrixSize]=CopyMatrix[j+s*MatrixSize];
                
                std::cout << this->EigenVectors[j+s*MatrixSize] << " ";
                
            }
            
            std::cout << std::endl;
            
            
        }*/
        
        // IN REVERSE ORDER AS LAPACK GIVES LOWEST TO HIGHEST //
        for(int s=0;s<MatrixSize;s++){
            
            
            // SET EIGENVALUE //
            this->EigenValues[((MatrixSize-1)-s)] = eVals[s];
            
            //std::cout << this->EigenValues[((MatrixSize-1)-s)] << std::endl;
            
            // SET EIGENVECTOR //
            for(int j=0;j<MatrixSize;j++){
                
                this->EigenVectors[j+((MatrixSize-1)-s)*MatrixSize]=CopyMatrix[j+s*MatrixSize];
                
                //std::cout << this->EigenVectors[j+((MatrixSize-1)-s)*MatrixSize] << " ";
                
            }
            
            //std::cout << std::endl;
            
            
        }
        
        
        // FREE WORKSPACE
        delete[] work;
        delete[] rwork;
        
        delete[] CopyMatrix;
        
        // FREE TEMP EIGENSYSTEM
        delete[] eVals;
        delete[] eVecs;
        
        
    }
    
public:
    
    // CHECK THAT EIGENFUNCTIONS ARE CORRECT //
    void Check(){
        
        INT CHECKVALUE=1;
        
        for(INT s=0;s<MatrixSize;s++){
            
            for(INT i=0;i<MatrixSize;i++){
                
                COMPLEX qVal=0.0;
                
                for(INT j=0;j<MatrixSize;j++){
                    qVal+=HermitianMatrix[i+MatrixSize*j]*EigenVectors[j+s*MatrixSize];
                    
                }
                // WARNING -- CHANGED PRECISION -- LOWER THAN WE'D LIKE //
                if(std::abs(qVal-EigenValues[s]*EigenVectors[i+s*MatrixSize])>std::pow(10.0,-8)){
                    CHECKVALUE=0;
                }
                
            }
            
        }
        
        if(CHECKVALUE==1){
            if(MPIBasic::ID==0){
                std::cerr << "#DIAGONALIZATION SUCCESSFUL" << std::endl;
            }
        }
        else{
            if(MPIBasic::ID==0){
                std::cerr << "#DIAGONALIZATION FAILED" << std::endl;
            }
            exit(0);
        }
        
    }
    
public:
    
    void Output(std::string fname){
        
        std::ofstream OutStream;
        OutStream.open(fname.c_str());
        
        for(INT s=0;s<MatrixSize;s++){
            
            OutStream << "#EVal= " << EigenValues[s] << std::endl;
            for(INT i=0;i<MatrixSize;i++){
                OutStream << real(EigenVectors[i+s*MatrixSize]) << " " << imag(EigenVectors[i+s*MatrixSize]) << " ";
            }
            OutStream << std::endl;
            
        }
        
        OutStream.close();
        
    }
    
    
public:
    
    // CONSTRUCTOR //
    Eigensystem(std::complex<double> *InputMatrix,INT InputSize){
        
        // SETUP //
        this->MatrixSize=InputSize;
        
        EigenVectors=new COMPLEX[MatrixSize*MatrixSize];
        EigenValues=new DOUBLE[MatrixSize];
        /*
        // SET TEMP MATRIX TO BE TRANSPOSED //
        this->HermitianMatrixUnTransposed=new std::complex<double>[MatrixSize*MatrixSize];
        
        // COPY
        std::memcpy(this->HermitianMatrixUnTransposed,InputMatrix,MatrixSize*MatrixSize*sizeof(std::complex<double>));
        */
        // SET MATRIX //
        this->HermitianMatrix=new std::complex<double>[MatrixSize*MatrixSize];
        std::memcpy(this->HermitianMatrix,InputMatrix,MatrixSize*MatrixSize*sizeof(std::complex<double>));

        /*
        // TRANSFORM TO ROW DOMINANT FORMAT LIKE GSL //
        for(INT i=0;i<InputSize;i++){
            for(INT j=0;j<InputSize;j++){
                
                this->HermitianMatrix[i+InputSize*j]=this->HermitianMatrixUnTransposed[j+InputSize*i];
                
            }
        }
        
        delete[] this->HermitianMatrixUnTransposed;
        */
        
    }
    
    // DESTRUCTOR //
    ~Eigensystem(){
        
        delete[] EigenVectors;
        delete[] EigenValues;
        
        delete HermitianMatrix;
    }
    
};