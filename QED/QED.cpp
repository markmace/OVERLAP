// DYNAMICAL QED FIELDS // 
namespace QED {
    
    COMPLEX *U; DOUBLE *E;
    
    INT Index2D(INT x,INT y){
        return MOD(x,Lattice::N[0])+Lattice::N[0]*MOD(y,Lattice::N[1]);
    }
    
    INT Index2D(INT x,INT y,INT mu){
        return mu+2*(MOD(x,Lattice::N[0])+Lattice::N[0]*MOD(y,Lattice::N[1]));
    }
    
    COMPLEX GetUVal(INT x,INT y,INT z,INT mu){
        
        if(mu<2){
            return QED::U[QED::Index2D(x,y,mu)];
        }
        else{
            return COMPLEX(1.0,0.0);
        }
        
    }
    
    void Init(){
        U=new COMPLEX[2*Lattice::N[0]*Lattice::N[1]];
        E=new DOUBLE[2*Lattice::N[0]*Lattice::N[1]];
    }
    
    void SetZero(){
        
        for(INT y=0;y<Lattice::N[1];y++){
            for(INT x=0;x<Lattice::N[0];x++){
                
                for(INT mu=0;mu<2;mu++){
                    
                    U[QED::Index2D(x,y,mu)]=COMPLEX(1.0,0.0);
                    E[QED::Index2D(x,y,mu)]=DOUBLE(0.0);
                    
                }
                
            }
        }
        
    }
    
    void PerformGaugeTranformation(COMPLEX *G){
        
        for(INT y=0;y<Lattice::N[1];y++){
            for(INT x=0;x<Lattice::N[0];x++){
                
                COMPLEX Ux=U[QED::Index2D(x,y,0)]; COMPLEX Uy=U[QED::Index2D(x,y,1)];
                
                U[QED::Index2D(x,y,0)]=G[QED::Index2D(x,y)]*Ux*conj(G[QED::Index2D(x+1,y)]);
                U[QED::Index2D(x,y,1)]=G[QED::Index2D(x,y)]*Uy*conj(G[QED::Index2D(x,y+1)]);
                
            }
        }
    }
    
    // SET CONSTANT MANGETIC FIELD -- c.f. 1111.4956 //
    void SetConstantMagneticField(INT Nb){
        
        for(INT y=0;y<Lattice::N[1];y++){
            for(INT x=0;x<Lattice::N[0];x++){
                
                U[QED::Index2D(x,y,1)]=std::exp(-2.0*M_PI*ComplexI*DOUBLE(Nb*x)/DOUBLE(Lattice::N[0]*Lattice::N[1]));
                
                if(x==Lattice::N[0]-1){
                    U[QED::Index2D(x,y,0)]=std::exp(2.0*M_PI*ComplexI*DOUBLE(Nb*y*Lattice::N[0])/DOUBLE(Lattice::N[0]*Lattice::N[1]));
                }
                else{
                    U[QED::Index2D(x,y,0)]=COMPLEX(1.0,0.0);
                }
                
                for(INT mu=0;mu<2;mu++){
                    E[QED::Index2D(x,y,mu)]=DOUBLE(0.0);
                }
                
                
            }
        }
        
    }
    
    DOUBLE ComputeMagneticField(INT x,INT y,INT z){
        
        COMPLEX UBox=GetUVal(x,y,z,0)*GetUVal(x+1,y,z,1)*conj(GetUVal(x,y+1,z,0))*conj(GetUVal(x,y,z,1));
        
        return -imag(UBox);
        
    }
    
    void SaveMagneticField(std::string fname){
        
        std::ofstream OutStream;
        OutStream.open(fname.c_str());
        
        for(INT z=0;z<Lattice::N[2];z++){
            for(INT y=0;y<Lattice::N[1];y++){
                for(INT x=0;x<Lattice::N[0];x++){
                    
                    OutStream << x << " " << y << " " << z << " " << ComputeMagneticField(x,y,z) << std::endl;
                }
            }
        }
        
        OutStream.close();
        
    }
    
}