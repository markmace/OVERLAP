#ifndef GAUGEF_H_GUARD
#define GAUGEF_H_GUARD
typedef struct { COMPLEX c[2]; } u_c;
typedef struct { u_c  c[2]; } u_cc;
#endif

#include <time.h>

namespace EvolutionOperator{
    
    // GAUGE LINKS //
    u_cc *Ucc;
    
    time_t t0,t1;
    clock_t c0,c1;
    
    // INITIALIZE //
    void Init(){
        Ucc=new u_cc[Lattice::N[0]*Lattice::N[1]*Lattice::N[2]*Lattice::Dimension];
    }

    // CONVERTS GAUGE LINKS INTO UCC FORMAT //
    void ConvertGaugeLinks(GaugeLinks *U,u_cc *UccLinks){
        
        for(INT z=0;z<U->N[2];z++){
            for(INT y=0;y<U->N[1];y++){
                for(INT x=0;x<U->N[0];x++){
                    
                    for(INT mu=0;mu<U->Dimension;mu++){
                        
                        INT CCIndex=U->Index3D(x,y,z)+mu*(U->N[0])*(U->N[1])*(U->N[2]);
                        
                        SU_Nc_MATRIX_FORMAT UMat[Nc*Nc];
                        SUNcGroup::Operations::GetMatrix(U->Get(x,y,z,mu),UMat);
                        
                        for(INT i=0;i<Nc;i++){
                            for(INT j=0;j<Nc;j++){
                                
                                // CHECK CONVENTION HERE //
                                UccLinks[CCIndex].c[i].c[j]=UMat[i+j*Nc];
                                
                            }
                        }
                        
                    }
                    
                }
            }
        }
        
    }
    
    
    // SETUP //
    void Setup(GaugeLinks *U){
        ConvertGaugeLinks(U,Ucc);        
    }
    
    // INCLUDE OVERLAP OPERATOR //
    #include "OPERATOR/OverlapOperator.cpp"
    
    // COMPUTE \partial_{t} Psi= -i Gamma0 [ -i DSlash_s +m] Psi //
    void ComputeSingleOverlapMode(FermionField *Psi,INT sIndex,wcf *DerivativePsi,DOUBLE LambdaMin,DOUBLE LambdaMax){
        
        //c0=clock();
        
        OverlapOperator::Compute(Ucc,Psi->GetMode(sIndex),DerivativePsi,LambdaMin,LambdaMax);
        
        //c1=clock();
        
        //printf("the elapsed cpu time in calculating overlap vector: %f\n",(float)(c1-c0)/CLOCKS_PER_SEC);
        
    }
    
    
}
