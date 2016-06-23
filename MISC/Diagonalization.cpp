#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

class Eigensystem {
    
private:
    
    // LINEAR DIMENSION OF THE MATRIX //
    INT MatrixSize;
    
private:
    
    // MATRIX //
    std::complex<double> *HermitianMatrix;
    
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
    
    // DIAGONALIZATION USING GSL LIBRARY //
    void Diagonalize(){
        
        // COPY MATRIX //
        std::complex<double> *CopyMatrix=new std::complex<double>[MatrixSize*MatrixSize];
        std::memcpy(CopyMatrix,this->HermitianMatrix,MatrixSize*MatrixSize*sizeof(std::complex<double>));
        
        // SET GSL MATRIX //
        gsl_matrix_complex_view GSLMatrix = gsl_matrix_complex_view_array (reinterpret_cast<double*>(CopyMatrix),MatrixSize,MatrixSize);
        
        // GSL ALLOCATION //
        gsl_vector *eVals = gsl_vector_alloc (MatrixSize);
        gsl_matrix_complex *eVecs = gsl_matrix_complex_alloc (MatrixSize, MatrixSize);
        
        // DIAGONALIZE //
        gsl_eigen_hermv_workspace * GSLMatrixWorkspace = gsl_eigen_hermv_alloc (MatrixSize);
        gsl_eigen_hermv (&GSLMatrix.matrix, eVals, eVecs, GSLMatrixWorkspace);
        gsl_eigen_hermv_free (GSLMatrixWorkspace);
        
        // SORT EIGENVALUES //
        gsl_eigen_hermv_sort (eVals, eVecs, GSL_EIGEN_SORT_VAL_DESC);
        /*
        // SET EIGENSYSTEM //
        for (int s=0;s<MatrixSize;s++){
            
            // SET EIGENVALUE //
            this->EigenValues[s] = gsl_vector_get (eVals, s);
            
            std::cout << this->EigenValues[s] << std::endl;
            
            // SET EIGENVECTOR //
            for(int j=0;j<MatrixSize;j++){
                
                
                gsl_complex val=gsl_matrix_complex_get(eVecs,j,s);
                this->EigenVectors[j+s*MatrixSize]=COMPLEX(GSL_REAL(val),GSL_IMAG(val));
                
                std::cout << this->EigenVectors[j+s*MatrixSize] << " ";
                
            }
            
            std::cout << std::endl;
            
            
        }
         
         */
        // SET EIGENSYSTEM //
        for (int s=0;s<MatrixSize;s++){
            
            // SET EIGENVALUE //
            this->EigenValues[s] = gsl_vector_get (eVals, s);
            
            //std::cout << this->EigenValues[s] << std::endl;
            
            //std::cout << "{";
            
            // SET EIGENVECTOR //
            for(int j=0;j<MatrixSize;j++){
                
                gsl_complex val=gsl_matrix_complex_get(eVecs,j,s);

                this->EigenVectors[j+s*MatrixSize]=COMPLEX(GSL_REAL(val),GSL_IMAG(val));
                /*
                std::cout << "{" <<real(this->EigenVectors[j+s*MatrixSize]) << "+I*" << imag(this->EigenVectors[j+s*MatrixSize]) << "}";
                if(j!=MatrixSize-1){
                    std::cout << ",";
                }
                */
                
            }
            
            //std::cout << "}"<< std::endl;
            
            
        }
        //std::cout << std::endl;
        
        
        // CLEANUP //
        gsl_vector_free (eVals);
        gsl_matrix_complex_free (eVecs);
        
        delete CopyMatrix;
    }
    
public:
    
    // CHECK THAT EIGENFUNCTIONS ARE CORRECT //
    void Check(){
        
        INT CHECKVALUE=1;
        
        for(INT s=0;s<MatrixSize;s++){
            
            for(INT i=0;i<MatrixSize;i++){
                
                COMPLEX qVal=0.0;
                
                for(INT j=0;j<MatrixSize;j++){
                    qVal+=HermitianMatrix[i*MatrixSize+j]*EigenVectors[j+s*MatrixSize];
                }
                
                if(std::abs(qVal-EigenValues[s]*EigenVectors[i+s*MatrixSize])>std::pow(10.0,-10)){
                    CHECKVALUE=0;
                    
                    std::cerr << s <<  " " << std::abs(qVal-EigenValues[s]*EigenVectors[i+s*MatrixSize]) << std::endl;
                }
                
            }
            
        }
        
        if(CHECKVALUE==1){
            std::cerr << "#DIAGONALIZATION SUCCESSFUL" << std::endl;
        }
        else{
            std::cerr << "#DIAGONALIZATION FAILED" << std::endl;
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
        
        // SET MATRIX //
        this->HermitianMatrix=new std::complex<double>[MatrixSize*MatrixSize];
        std::memcpy(this->HermitianMatrix,InputMatrix,MatrixSize*MatrixSize*sizeof(std::complex<double>));
        
    }
    
    // DESTRUCTOR //
    ~Eigensystem(){
        
        delete EigenVectors;
        delete EigenValues;
        
        delete HermitianMatrix;
    }
    
};
