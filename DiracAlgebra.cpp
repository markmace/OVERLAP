namespace DiracAlgebra{
    
    static const INT SpinorComponents=4;
    
    COMPLEX Gamma0[SpinorComponents][SpinorComponents]={{1.0,0.0,0.0,0.0},{0.0,1.0,0.0,0.0},{0.0,0.0,-1.0,0.0},{0.0,0.0,0.0,-1.0}};
    COMPLEX GammaX[SpinorComponents][SpinorComponents]={{0.0,0.0,0.0,1.0},{0.0,0.0,1.0,0.0},{0.0,-1.0,0.0,0.0},{-1.0,0.0,0.0,0.0}};
    COMPLEX GammaY[SpinorComponents][SpinorComponents]={{0.0,0.0,0.0,-ComplexI},{0.0,0.0,ComplexI,0.0},{0.0,ComplexI,0.0,0.0},{-ComplexI,0.0,0.0,0.0}};
    COMPLEX GammaZ[SpinorComponents][SpinorComponents]={{0.0,0.0,1.0,0.0},{0.0,0.0,0.0,-1.0},{-1.0,0.0,0.0,0.0},{0.0,1.0,0.0,0.0}};
    COMPLEX Gamma5[SpinorComponents][SpinorComponents]={{0.0,0.0,1.0,0.0},{0.0,0.0,0.0,1.0},{1.0,0.0,0.0,0.0},{0.0,1.0,0.0,0.0}};

    
    COMPLEX IdentityMatrix4D[SpinorComponents][SpinorComponents]={{1.0,0.0,0.0,0.0},{0.0,1.0,0.0,0.0},{0.0,0.0,1.0,0.0},{0.0,0.0,0.0,1.0}};
    
    
    void Output(){
        
        std::cerr << "#GAMMA 0" << std::endl;
        
        for(INT alpha=0;alpha<SpinorComponents;alpha++){
            for(INT beta=0;beta<SpinorComponents;beta++){
                
                std::cerr << Gamma0[alpha][beta] << " ";
            }
            std::cerr << std::endl;
        }
        
        
        std::cerr << "#GAMMA X" << std::endl;
        
        for(INT alpha=0;alpha<SpinorComponents;alpha++){
            for(INT beta=0;beta<SpinorComponents;beta++){
                
                std::cerr << GammaX[alpha][beta] << " ";
            }
            std::cerr << std::endl;
        }
        
        std::cerr << "#GAMMA Y" << std::endl;
        
        for(INT alpha=0;alpha<SpinorComponents;alpha++){
            for(INT beta=0;beta<SpinorComponents;beta++){
                
                std::cerr << GammaY[alpha][beta] << " ";
            }
            std::cerr << std::endl;
        }
        
        
        std::cerr << "#GAMMA Z" << std::endl;
        
        for(INT alpha=0;alpha<SpinorComponents;alpha++){
            for(INT beta=0;beta<SpinorComponents;beta++){
                
                std::cerr << GammaZ[alpha][beta] << " ";
            }
            std::cerr << std::endl;
        }
        
        std::cerr << "#GAMMA 5" << std::endl;
        
        for(INT alpha=0;alpha<SpinorComponents;alpha++){
            for(INT beta=0;beta<SpinorComponents;beta++){
                
                std::cerr << Gamma5[alpha][beta] << " ";
            }
            std::cerr << std::endl;
        }
    }
    
}