namespace InitialConditions{
    
    void SetZero(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *U,ElectricFields *E){
        
        for(INT z=zLow;z<=zHigh;z++){
            for(INT y=yLow;y<=yHigh;y++){
                for(INT x=xLow;x<=xHigh;x++){
                    
                    //SET ELECTRIC FIELDS TO ZERO
                    for(int mu=0;mu<Lattice::Dimension;mu++){
                        for(int a=0;a<SUNcAlgebra::VectorSize;a++){
                            E->Get(x,y,z,mu,a)[0]=0.0;
                        }
                    }
                    
                    //SET GAUGE LINKS TO UNITY
                    for(int mu=0;mu<Lattice::Dimension;mu++){
                        COPY_SUNcMatrix(U->Get(x,y,z,mu),SUNcGroup::UnitMatrix);
                    }
                    
                }
            }
        }
    }
    
    void SetZero(GaugeLinks *U,ElectricFields *E){
        
        SetZero(0,U->N[0]-1,0,U->N[1]-1,0,U->N[2]-1,U,E);
        
    }
    
    
    void SetZero(){
                
        SetZero(GLinks::U,EFields::E);
        
    }
    
    
    
}
