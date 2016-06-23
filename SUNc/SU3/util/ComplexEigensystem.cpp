#ifndef __COMPLEX__EIGENSYSTEM__CPP__
#define __COMPLEX__EIGENSYSTEM__CPP__

/* ZGEEV prototype */
extern "C" {
    void zgeev(char const* jobvl, char const* jobvr, int* n, COMPLEX* a, int* lda, COMPLEX* w, COMPLEX* vl, int* ldvl,COMPLEX* vr, int* ldvr, COMPLEX* work, int* lwork, double* rwork, int* info );
}

/* Parameters */
#define LDA Nc
#define LDVL Nc
#define LDVR Nc

////////////////////////////////////////
//                                    //
// NOTE: LAPACK DEFINES MATRICES AS   //
//                                    //
//          A_0 A_3 A_6               //
//          A_1 A_4 A_7               //
//          A_2 A_5 A_8               //
//                                    //
////////////////////////////////////////

namespace ComplexEigensystem{
    
    INT GetComplexEigensystem(COMPLEX A[9],COMPLEX EVals[3], COMPLEX EVecs[9]) {
        
        int n = Nc, lda = LDA, ldvl = LDVL, ldvr = LDVR, info, lwork;
        COMPLEX wkopt;
        COMPLEX* work;
        // LOCAL BUFFERS //
        // rwork dimension should be at least 2*n //
        double rwork[2*Nc];
        COMPLEX vl[1]; // LEFT EIGENVECTORS NEVER USED
        
        // QUERY AND ALLOCATE WORKSPACE //
        lwork = -1;
        zgeev( "N", "Vectors", &n, A, &lda, EVals, vl, &ldvl, EVecs, &ldvr, &wkopt, &lwork, rwork, &info );
        lwork = int(real(wkopt));
        work = new COMPLEX[lwork*sizeof(COMPLEX)];
        // SOLVE EIGENPROBLEM //
        zgeev( "N", "Vectors", &n, A, &lda, EVals, vl, &ldvl, EVecs, &ldvr,work, &lwork, rwork, &info );
        // CHECK CONVERGENCE //
        if( info > 0 ) {
            std::cerr << "The algorithm failed to compute eigenvalues." << std::endl;
            return -1;
            exit( 1 );
        }
        
        // FREE WORKSPACE
        delete[] work;
        
        return 0;
    }
}

#endif

