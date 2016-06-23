#ifndef __FERMIONFIELDS__CPP__
#define __FERMIONFIELDS__CPP__

#include "../MPI/ModePartition.cpp"

#ifndef VECDEF_H_GUARD
#define VECDEF_H_GUARD
typedef struct { COMPLEX f[4]; } wf;
typedef struct { wf  c[2]; } wcf;
#endif

class FermionField{
    
public:
    
    static const INT Dimension=3;
    
    // DIMENSIONS //
    INT N[Dimension];
    
    // VOLUME //
    INT Volume;
    
    // LATTICE SPACING //
    DOUBLE a[Dimension];
    
    // NUMBER OF SAMPLES //
    INT NumberOfSamples;
    
    // 3D LATTICE INDEX //
    INT Index3D(INT x,INT y,INT z){
        return MOD(x,N[0])+N[0]*(MOD(y,N[1])+N[1]*MOD(z,N[2]));
    }
    
public:
    
    // GET MODE INDEX //
    INT GetModeIndex(INT pXIndex,INT pYIndex,INT pZIndex,INT SpinIndex,INT ColorIndex,INT ParticleIndex){
        
        return ParticleIndex+2*(ColorIndex+Nc*(SpinIndex+2*(pXIndex+N[0]*(pYIndex+N[1]*pZIndex))));
        
    }
    
    // GET MODE QUANTUM NUMBER //
    void GetLocalModeQuantumNumbers(INT sLoc,INT &pXIndex,INT &pYIndex,INT &pZIndex,INT &SpinIndex,INT &ColorIndex,INT &ParticleIndex){
        
        INT sGlobal=sLoc+MPIModePartition::FirstIndex;
        
        ParticleIndex=sGlobal%2; sGlobal=(sGlobal-ParticleIndex)/2;
        ColorIndex=sGlobal%Nc; sGlobal=(sGlobal-ColorIndex)/Nc;
        SpinIndex=sGlobal%2; sGlobal=(sGlobal-SpinIndex)/2;
        pXIndex=sGlobal%(N[0]); sGlobal=(sGlobal-pXIndex)/(N[0]);
        pYIndex=sGlobal%(N[1]); sGlobal=(sGlobal-pYIndex)/(N[1]);
        pZIndex=sGlobal%(N[2]); sGlobal=(sGlobal-pZIndex);
        
        if(sGlobal!=0){
            std::cerr << "##### ERRRROOOOOOORRRR IN INDEXING ####" << std::endl;
        }
        
    }
    
    // USED IN OVERLAP B-FIELD EFUNCTIONS
    void GetLocalModeQuantumNumbers(INT sLoc,INT &LevelIndex,INT &ParticleIndex){
        
        INT sGlobal=sLoc+MPIModePartition::FirstIndex;
        
        ParticleIndex=sGlobal%2; sGlobal=(sGlobal-ParticleIndex)/2;
        LevelIndex=sGlobal%(2*N[0]*N[1]*N[2]*Nc); sGlobal=(sGlobal-LevelIndex)/(2*N[0]*N[1]*N[2]*Nc);
        
        if(sGlobal!=0){
            std::cerr << "##### ERRRROOOOOOORRRR IN INDEXING ####" << std::endl;
        }
        
    }
    
    // FOR WILSON B-FIELD EFUNCTIONS
    void GetLocalModeQuantumNumbers(INT sLoc,INT &LevelIndex,INT &pZIndex,INT &ColorIndex,INT &ParticleIndex){
        
        INT sGlobal=sLoc+MPIModePartition::FirstIndex;
        
        ParticleIndex=sGlobal%2; sGlobal=(sGlobal-ParticleIndex)/2;
        ColorIndex=sGlobal%Nc; sGlobal=(sGlobal-ColorIndex)/Nc;
        LevelIndex=sGlobal%(2*N[0]*N[1]); sGlobal=(sGlobal-LevelIndex)/(2*N[0]*N[1]);
        pZIndex=sGlobal%(N[2]); sGlobal=(sGlobal-pZIndex);
        
        if(sGlobal!=0){
            std::cerr << "##### ERRRROOOOOOORRRR IN INDEXING ####" << std::endl;
        }
        
    }
    
    void GetGlobalModeQuantumNumbers(INT sGlobalInput,INT &pXIndex,INT &pYIndex,INT &pZIndex,INT &SpinIndex,INT &ColorIndex,INT &ParticleIndex){
        
        INT sGlobal=sGlobalInput;
        
        ParticleIndex=sGlobal%2; sGlobal=(sGlobal-ParticleIndex)/2;
        ColorIndex=sGlobal%Nc; sGlobal=(sGlobal-ColorIndex)/Nc;
        SpinIndex=sGlobal%2; sGlobal=(sGlobal-SpinIndex)/2;
        pXIndex=sGlobal%(N[0]); sGlobal=(sGlobal-pXIndex)/(N[0]);
        pYIndex=sGlobal%(N[1]); sGlobal=(sGlobal-pYIndex)/(N[1]);
        pZIndex=sGlobal%(N[2]); sGlobal=(sGlobal-pZIndex);
        
        if(sGlobal!=0){
            std::cerr << "##### ERRRROOOOOOORRRR IN INDEXING ####" << std::endl;
        }
        
    }
    
    DOUBLE CommutatorExpectationValue(INT sLoc){
        
        INT sGlobal=sLoc+MPIModePartition::FirstIndex;
        
        if(sGlobal%2==0){
            return 1;
        }
        
        else{
            return -1;
        }
        
        
    }
    
    DOUBLE VacuumExpectationValue(INT sLoc){
        
        INT sGlobal=sLoc+MPIModePartition::FirstIndex;
        
        if(sGlobal%2==0){
            return 0;
        }
        
        else{
            return 1;
        }
        
    }
    
#define HELICITY_POSITIVE_FLAG 0
#define HELICITY_NEGATIVE_FLAG 1
    
#define POSITIVE_ENERGY_FLAG 0
#define NEGATIVE_ENERGY_FLAG 1
    
    // GET FREE DIRAC SPINOR IN HELICITY UP/DOWN BASIS FOR WILSON FERMION //
    void GetWilsonDiracSpinor(INT pXIndex,INT pYIndex,INT pZIndex,INT SpinIndex,INT ParticleIndex,COMPLEX *u,DOUBLE &Frequency){
        
        // GET LATTICE MOMENTUM //
        DOUBLE ppx,ppy,ppz,pSqr,pAbs,mEff,wEff,HelicityEigenvalue;
        
        // PARTICLES GO AS EXP(-iWT+iPX) //
        if(ParticleIndex==POSITIVE_ENERGY_FLAG){
            
            ppx=sin(2.0*M_PI*pXIndex/DOUBLE(N[0]))/a[0];
            ppy=sin(2.0*M_PI*pYIndex/DOUBLE(N[1]))/a[1];
            ppz=sin(2.0*M_PI*pZIndex/DOUBLE(N[2]))/a[2];
            
            pSqr=SQR(ppx)+SQR(ppy)+SQR(ppz);    pAbs=sqrt(pSqr);
            
            mEff=mFermion+2.0*rWilson*(SQR(sin(M_PI*pXIndex/DOUBLE(N[0])))/a[0]+SQR(sin(M_PI*pYIndex/DOUBLE(N[1])))/a[1]+SQR(sin(M_PI*pZIndex/DOUBLE(N[2])))/a[2]);
            wEff=sqrt(pSqr+SQR(mEff));
            
            Frequency=wEff;
            
        }
        // ANTI-PARTICLES GO AS EXP(iWT-iPX) //
        else if(ParticleIndex==NEGATIVE_ENERGY_FLAG){
            
            ppx=sin(2.0*M_PI*(N[0]-pXIndex)/DOUBLE(N[0]))/a[0];
            ppy=sin(2.0*M_PI*(N[1]-pYIndex)/DOUBLE(N[1]))/a[1];
            ppz=sin(2.0*M_PI*(N[2]-pZIndex)/DOUBLE(N[2]))/a[2];
            
            pSqr=SQR(ppx)+SQR(ppy)+SQR(ppz); pAbs=sqrt(pSqr);
            
            mEff=mFermion+2.0*rWilson*(SQR(sin(M_PI*(N[0]-pXIndex)/DOUBLE(N[0])))/a[0]+SQR(sin(M_PI*(N[1]-pYIndex)/DOUBLE(N[1])))/a[1]+SQR(sin(M_PI*(N[2]-pZIndex)/DOUBLE(N[2])))/a[2]);
            wEff=sqrt(pSqr+SQR(mEff));
            
            Frequency=-wEff;
        }
        else{
            std::cerr << "ERROR" << std::endl;
            exit(0);
        }
        
        
        // COMPUTE WEYL SPINORS IN HELICITY BASIS //
        COMPLEX Phi[2];
        
        if( (pXIndex==0 || pXIndex==N[0]/2) && (pYIndex==0 || pYIndex==N[1]/2) ){
            
            if(SpinIndex==HELICITY_POSITIVE_FLAG){
                Phi[0]=1.0; Phi[1]=0.0;  HelicityEigenvalue=+ppz/std::abs(ppz);
            }
            else if(SpinIndex==HELICITY_NEGATIVE_FLAG){
                Phi[0]=0.0; Phi[1]=1.0;  HelicityEigenvalue=-ppz/std::abs(ppz);
            }
            else{
                std::cerr << "ERROR" << std::endl;
                exit(0);
            }
            
            // SPECIAL CASE //
            if((pZIndex==0 || pZIndex==N[2]/2)){
                HelicityEigenvalue=0;
            }
            
        }
        
        else{
            
            
            if(SpinIndex==HELICITY_POSITIVE_FLAG){
                HelicityEigenvalue=+1.0;
            }
            else if(SpinIndex==HELICITY_NEGATIVE_FLAG){
                HelicityEigenvalue=-1.0;
            }
            else{
                std::cerr << "ERROR" << std::endl;
                exit(0);
            }
            
            COMPLEX q=-(ppz-HelicityEigenvalue*pAbs)/(ppx-ComplexI*ppy); DOUBLE d=sqrt(1.0+SQR_ABS(q));
            
            Phi[0]=1.0/d; Phi[1]=q/d;
            
        }
        
        // DETERMINE DIRAC SPINOR //
        if((pXIndex==0 || pXIndex==N[0]/2) && (pYIndex==0 || pYIndex==N[1]/2) && (pZIndex==0 || pZIndex==N[2]/2)){
            
            // NOTE THAT WHEN rWilson IS NEGATIVE ONE HAS mEff<0 //
            if(ParticleIndex==POSITIVE_ENERGY_FLAG){
                
                if(mEff>=0){
                    u[0]=Phi[0]; u[1]=Phi[1]; u[2]=0.0; u[3]=0.0;
                }
                else{
                    u[0]=0.0; u[1]=0.0; u[2]=Phi[0]; u[3]=Phi[1];
                }
                
            }
            else if(ParticleIndex==NEGATIVE_ENERGY_FLAG){
                
                if(mEff>=0){
                    u[0]=0.0; u[1]=0.0; u[2]=Phi[0]; u[3]=Phi[1];
                }
                else{
                    u[0]=Phi[0]; u[1]=Phi[1]; u[2]=0.0; u[3]=0.0;
                }
            }
            else{
                std::cerr << "ERROR" << std::endl;
                exit(0);
            }
            
        }
        
        else{
            
            DOUBLE DiracNorm=sqrt(2.0*Frequency*(Frequency-mEff)/pSqr);
            
            u[0]=Phi[0]/DiracNorm;
            u[1]=Phi[1]/DiracNorm;
            
            u[2]=HelicityEigenvalue*(Frequency-mEff)/pAbs*(Phi[0]/DiracNorm);
            u[3]=HelicityEigenvalue*(Frequency-mEff)/pAbs*(Phi[1]/DiracNorm);
            
        }
        
        // CHECK CORRECT EIGENVALUE PROPERTIES //
        COMPLEX Hamiltonian[4][4];
        
        for(INT alpha=0;alpha<4;alpha++){
            for(INT beta=0;beta<4;beta++){
                
                using namespace DiracAlgebra;
                
                Hamiltonian[alpha][beta]=mEff*Gamma0[alpha][beta];
                
                for(INT gamma=0;gamma<4;gamma++){
                    Hamiltonian[alpha][beta]+=Gamma0[alpha][gamma]*(ppx*GammaX[gamma][beta]+ppy*GammaY[gamma][beta]+ppz*GammaZ[gamma][beta]); // Gamma_0 pSlash for Wilson
                }
                
            }
        }
        
        
        for(INT alpha=0;alpha<4;alpha++){
            
            COMPLEX Check1=Frequency*u[alpha];
            
            COMPLEX Check2=0.0;
            for(INT beta=0;beta<4;beta++){
                Check2+=Hamiltonian[alpha][beta]*u[beta];
            }
            
            if(std::abs(Check1-Check2)>std::pow(10.0,-6)){
                std::cerr << "ERR " << pXIndex << " " << pYIndex << " " << pZIndex << " " << SpinIndex << " " << ParticleIndex << " " << Check1 << " " << Check2 << " " << Check1/Check2 << std::endl;
                
                std::cerr << Phi[0] << " " << Phi[1] << " " << mEff << " " << Frequency << std::endl;
            }
            
        }
        
        
        // CHECK NORMALIZATION //
        
        DOUBLE CheckNorm=SQR_ABS(u[0])+SQR_ABS(u[1])+SQR_ABS(u[2])+SQR_ABS(u[3]);
        
        if(std::abs(CheckNorm-1.0)>std::pow(10.0,-6)){
            std::cerr << "ERR " << CheckNorm  << " " << SQR_ABS(Phi[0])+SQR_ABS(Phi[1]) << " " << ParticleIndex << " " << SpinIndex << std::endl;
        }
        
        
        
        
    }
    
    // GET FREE DIRAC PHASE FOR WILSON FERMION //
    void GetWilsonPhase(INT pXIndex,INT pYIndex,INT pZIndex,INT x,INT y,INT z,INT ParticleIndex,COMPLEX &Phase){
        
        if(ParticleIndex==0){
            Phase=exp(+2.0*M_PI*ComplexI*(pXIndex*x/DOUBLE(N[0])+pYIndex*y/DOUBLE(N[1])+pZIndex*z/DOUBLE(N[2])));
        }
        if(ParticleIndex==1){
            Phase=exp(-2.0*M_PI*ComplexI*(pXIndex*x/DOUBLE(N[0])+pYIndex*y/DOUBLE(N[1])+pZIndex*z/DOUBLE(N[2])));
        }
    }
    
    // GET FREE DIRAC SPINOR FOR OVERLAP FERMION //
    void GetOverlapDiracSpinor(INT pXIndex,INT pYIndex,INT pZIndex,INT SpinIndex,INT ParticleIndex,COMPLEX *u,DOUBLE &Frequency){
        
        // GET LATTICE MOMENTUM //
        DOUBLE hpx,hpy,hpz,hp5,sSqr,sAbs,ppx,ppy,ppz,pSqr,pAbs,mEff,wEff;
        COMPLEX HelicityEigenvalue;
        
        // PARTICLES GO AS EXP(-iWT+iPX) //
        if(ParticleIndex==POSITIVE_ENERGY_FLAG){
            
            hpx=sin(2.0*M_PI*pXIndex/DOUBLE(N[0]))/a[0];
            hpy=sin(2.0*M_PI*pYIndex/DOUBLE(N[1]))/a[1];
            hpz=sin(2.0*M_PI*pZIndex/DOUBLE(N[2]))/a[2];
            
            hp5=3.0-MDWHeight-cos(2.0*M_PI*pXIndex/DOUBLE(N[0]))/a[0]-cos(2.0*M_PI*pYIndex/DOUBLE(N[1]))/a[1]-cos(2.0*M_PI*pZIndex/DOUBLE(N[2]))/a[2];
            
            sSqr=SQR(hpx)+SQR(hpy)+SQR(hpz)+SQR(hp5);    sAbs=sqrt(sSqr);
            
            ppx=MDWHeight*hpx/sAbs;
            ppy=MDWHeight*hpy/sAbs;
            ppz=MDWHeight*hpz/sAbs;
            
            pSqr=SQR(ppx)+SQR(ppy)+SQR(ppz);    pAbs=sqrt(pSqr);
            
            mEff=MDWHeight*(1+hp5/sAbs);
            wEff=sqrt(pSqr+SQR(mEff));
            
            Frequency=wEff;
            
        }
        // ANTI-PARTICLES GO AS EXP(iWT-iPX) //
        else if(ParticleIndex==NEGATIVE_ENERGY_FLAG){
            
            hpx=sin(2.0*M_PI*(N[0]-pXIndex)/DOUBLE(N[0]))/a[0];
            hpy=sin(2.0*M_PI*(N[1]-pYIndex)/DOUBLE(N[1]))/a[1];
            hpz=sin(2.0*M_PI*(N[2]-pZIndex)/DOUBLE(N[2]))/a[2];
            
            hp5=3.0-MDWHeight-cos(2.0*M_PI*(N[0]-pXIndex)/DOUBLE(N[0]))/a[0]-cos(2.0*M_PI*(N[1]-pYIndex)/DOUBLE(N[1]))/a[1]-cos(2.0*M_PI*(N[2]-pZIndex)/DOUBLE(N[2]))/a[2];
            
            sSqr=SQR(hpx)+SQR(hpy)+SQR(hpz)+SQR(hp5);    sAbs=sqrt(sSqr);
            
            ppx=MDWHeight*hpx/sAbs;
            ppy=MDWHeight*hpy/sAbs;
            ppz=MDWHeight*hpz/sAbs;
            
            pSqr=SQR(ppx)+SQR(ppy)+SQR(ppz);    pAbs=sqrt(pSqr);
            
            mEff=MDWHeight*(1+hp5/sAbs);
            wEff=sqrt(pSqr+SQR(mEff));
            
            Frequency=-wEff;
        }
        else{
            std::cerr << "ERROR" << std::endl;
            exit(0);
        }
        
        // COMPUTE WEYL SPINORS IN HELICITY BASIS //
        COMPLEX Phi[2];
        
        if( (pXIndex==0 || pXIndex==N[0]/2) && (pYIndex==0 || pYIndex==N[1]/2) ){
            
            if(SpinIndex==HELICITY_POSITIVE_FLAG){
                Phi[0]=1.0; Phi[1]=0.0;  HelicityEigenvalue=+ppz/std::abs(ppz);
            }
            else if(SpinIndex==HELICITY_NEGATIVE_FLAG){
                Phi[0]=0.0; Phi[1]=1.0;  HelicityEigenvalue=-ppz/std::abs(ppz);
            }
            else{
                std::cerr << "ERROR" << std::endl;
                exit(0);
            }
            
            // SPECIAL CASE //
            if((pZIndex==0 || pZIndex==N[2]/2)){
                HelicityEigenvalue=0;
            }
            
        }
        
        else{
            
            
            if(SpinIndex==HELICITY_POSITIVE_FLAG){
                HelicityEigenvalue=+1.0;
            }
            else if(SpinIndex==HELICITY_NEGATIVE_FLAG){
                HelicityEigenvalue=-1.0;
            }
            else{
                std::cerr << "ERROR" << std::endl;
                exit(0);
            }
            
            COMPLEX q=-(ppz-HelicityEigenvalue*pAbs)/(ppx-ComplexI*ppy); DOUBLE d=sqrt(1.0+SQR_ABS(q));
            
            Phi[0]=1.0/d; Phi[1]=q/d;
            
        }
        
        // DETERMINE DIRAC SPINOR //
        if((pXIndex==0 || pXIndex==N[0]/2) && (pYIndex==0 || pYIndex==N[1]/2) && (pZIndex==0 || pZIndex==N[2]/2)){
            
            // NOTE THAT WHEN rWilson IS NEGATIVE ONE HAS mEff<0 //
            if(ParticleIndex==POSITIVE_ENERGY_FLAG){
                
                if(mEff>=0){
                    u[0]=Phi[0]; u[1]=Phi[1]; u[2]=0.0; u[3]=0.0;
                }
                else{
                    u[0]=0.0; u[1]=0.0; u[2]=Phi[0]; u[3]=Phi[1];
                }
                
            }
            else if(ParticleIndex==NEGATIVE_ENERGY_FLAG){
                
                if(mEff>=0){
                    u[0]=0.0; u[1]=0.0; u[2]=Phi[0]; u[3]=Phi[1];
                }
                else{
                    u[0]=Phi[0]; u[1]=Phi[1]; u[2]=0.0; u[3]=0.0;
                }
            }
            else{
                std::cerr << "ERROR" << std::endl;
                exit(0);
            }
            
        }
        
        else{
            
            DOUBLE DiracNorm=sqrt(2.0*Frequency*(Frequency-mEff)/pSqr);
            
            u[0]=Phi[0]/DiracNorm;
            u[1]=Phi[1]/DiracNorm;
            
            u[2]=HelicityEigenvalue*(Frequency-mEff)/pAbs*(Phi[0]/DiracNorm);
            u[3]=HelicityEigenvalue*(Frequency-mEff)/pAbs*(Phi[1]/DiracNorm);
            
        }
        
        // CHECK CORRECT EIGENVALUE PROPERTIES //
        COMPLEX Hamiltonian[4][4];
        
        for(INT alpha=0;alpha<4;alpha++){
            for(INT beta=0;beta<4;beta++){
                
                using namespace DiracAlgebra;
                
                Hamiltonian[alpha][beta]=0.0;
                
                for(INT gamma=0;gamma<4;gamma++){
                    Hamiltonian[alpha][beta]+=MDWHeight*Gamma0[alpha][gamma]*(IdentityMatrix4D[gamma][beta]+(hpx*GammaX[gamma][beta]+hpy*GammaY[gamma][beta]+hpz*GammaZ[gamma][beta]+hp5*IdentityMatrix4D[gamma][beta])/sAbs);
                }
                
            }
        }
        
        
        for(INT alpha=0;alpha<4;alpha++){
            
            COMPLEX Check1=Frequency*u[alpha];
            
            COMPLEX Check2=0.0;
            for(INT beta=0;beta<4;beta++){
                Check2+=Hamiltonian[alpha][beta]*u[beta];
            }
            
            if(std::abs(Check1-Check2)>std::pow(10.0,-6)){
                std::cerr << "ERR " << pXIndex << " " << pYIndex << " " << pZIndex << " " << SpinIndex << " " << ParticleIndex << " " << Check1 << " " << Check2 << " " << Check1/Check2 << std::endl;
                
                std::cerr << Phi[0] << " " << Phi[1] << " " << mEff << " " << Frequency << std::endl;
            }
            
        }
        
        
        // CHECK NORMALIZATION //
        
        DOUBLE CheckNorm=SQR_ABS(u[0])+SQR_ABS(u[1])+SQR_ABS(u[2])+SQR_ABS(u[3]);
        
        if(std::abs(CheckNorm-1.0)>std::pow(10.0,-6)){
            std::cerr << "ERR " << CheckNorm  << " " << SQR_ABS(Phi[0])+SQR_ABS(Phi[1]) << " " << ParticleIndex << " " << SpinIndex << std::endl;
        }
        
        
    }
    
    // GET PHASE FOR OVERLAP FERMION //
    void GetOverlapPhase(INT pXIndex,INT pYIndex,INT pZIndex,INT x,INT y,INT z,INT ParticleIndex,COMPLEX &Phase){
        
        if(ParticleIndex==0){
            Phase=exp(+2.0*M_PI*ComplexI*(pXIndex*x/DOUBLE(N[0])+pYIndex*y/DOUBLE(N[1])+pZIndex*z/DOUBLE(N[2])));
        }
        if(ParticleIndex==1){
            Phase=exp(-2.0*M_PI*ComplexI*(pXIndex*x/DOUBLE(N[0])+pYIndex*y/DOUBLE(N[1])+pZIndex*z/DOUBLE(N[2])));
        }
    }
    
private:
    
    // SPINOR FIELDS //
    wcf **SpinorField;
    
public:
    
    COMPLEX *Get(INT x,INT y,INT z,INT alpha,INT i,INT s){
        
        return &SpinorField[s][Index3D(x,y,z)].c[i].f[alpha];
    }
    
    COMPLEX GetValue(INT x,INT y,INT z,INT alpha,INT i,INT s){
        
        return SpinorField[s][Index3D(x,y,z)].c[i].f[alpha];
    }
    
    wcf* GetMode(INT sIndex){
        return SpinorField[sIndex];
    }
    
public:
    
    void SetZero(){
        
        for(INT s=0;s<NumberOfSamples;s++){
            
            for(INT z=0;z<N[2];z++){
                for(INT y=0;y<N[1];y++){
                    for(INT x=0;x<N[0];x++){
                        for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){
                            for(INT i=0;i<Nc;i++){
                                
                                this->Get(x,y,z,alpha,i,s)[0]=0.0;
                            }
                        }
                    }
                }
            }
            
        }
        
    }
    
public:
    
    std::string ToString(INT x,INT y,INT z,INT s){
        
        std::stringstream ss;
        
        for(INT i=0;i<Nc;i++){
            for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){
                ss << real(this->GetValue(x,y,z,alpha,i,s)) << " " << imag(this->GetValue(x,y,z,alpha,i,s)) << " ";
            }
        }
        
        return ss.str();
    }
    
public:
    
    // CONSTRUCTOR //
    FermionField(INT Nx,INT Ny,INT Nz,DOUBLE aX,DOUBLE aY,DOUBLE aZ,INT NumberOfSamples){
        
        // SET DIMENSIONS //
        this->N[0]=Nx; this->N[1]=Ny; this->N[2]=Nz;    this->Volume=Nx*Ny*Nz;
        
        // SET LATTICE SPACINGS //
        this->a[0]=aX; this->a[1]=aY; this->a[2]=aZ;
        
        // SET NUMBER OF SAMPLES //
        this->NumberOfSamples=NumberOfSamples;
        
        // ALLOCATE //
        SpinorField=new wcf*[NumberOfSamples];
        
        for(INT s=0;s<NumberOfSamples;s++){
            SpinorField[s]=new wcf[N[0]*N[1]*N[2]];
        }
        
    }
    
    // DESTRUCTOR //
    ~FermionField(){
        
        for(INT s=0;s<NumberOfSamples;s++){
            delete SpinorField[s];
        }
        
        delete SpinorField;
    }
    
    
};

namespace Fermions {
    
    // FERMION FIELDS //
    FermionField *Psi; FermionField *PsiMid;
    
    
    // INITIALIZE //
    void Init(){
        
        // INITIALIZE PARTITION //
        MPIModePartition::Init();
        
        // INITIALIZE FIELDS //
        Psi=new FermionField(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2],MPIModePartition::ModesPerNode);
        PsiMid=new FermionField(Lattice::N[0],Lattice::N[1],Lattice::N[2],Lattice::a[0],Lattice::a[1],Lattice::a[2],MPIModePartition::ModesPerNode);
        
    }
    
}

#endif