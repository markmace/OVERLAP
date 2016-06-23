namespace InitialConditions{
    
    
    namespace HandmadeSphaleron{
        
        // CENTER THE MAP ON THE LATTICE TORUS //
        void LocalizedLatticeMap(INT ix, INT iy, INT iz, DOUBLE &xT,DOUBLE &yT,DOUBLE &zT){
            
            // GET INVERSE SPHALERON WIDTH //
            DOUBLE bbx=1.0/rSphaleron; DOUBLE bby=1.0/rSphaleron; DOUBLE bbz=1.0/rSphaleron;
                    
            // GET COORDINATE MAP //
            xT = 2.*M_PI*( atan( bbx*DOUBLE(ix-Lattice::N[0]/2.0) ) - atan( bbx*DOUBLE(-Lattice::N[0]/2.0) ) ) / ( atan( bbx*DOUBLE(Lattice::N[0]/2.0) ) - atan( -bbx*DOUBLE(Lattice::N[0]/2.0) ) );
            yT = 2.*M_PI*( atan( bby*DOUBLE(iy-Lattice::N[1]/2.0) ) - atan( bby*DOUBLE(-Lattice::N[1]/2.0) ) ) / ( atan( bby*DOUBLE(Lattice::N[1]/2.0) ) - atan( -bby*DOUBLE(Lattice::N[1]/2.0) ) );
            zT = 2.*M_PI*( atan( bbz*DOUBLE(iz-Lattice::N[2]/2.0) ) - atan( bbz*DOUBLE(-Lattice::N[2]/2.0) ) ) / ( atan( bbz*DOUBLE(Lattice::N[2]/2.0) ) - atan( -bbz*DOUBLE(Lattice::N[2]/2.0) ) );
            
        }
        
        
        // PERFORM THE MAP FROM LATTICE TORUS TO ONE POINT COMPACTIFIED R^3 //
        void TorusProjection(DOUBLE xT,DOUBLE yT,DOUBLE zT, DOUBLE &x,DOUBLE &y, DOUBLE &z, bool &InfFlag){
            
            if( xT == 0.0 || xT == 2.0*M_PI || yT == 0.0 || yT == 2.0*M_PI || zT == 0.0 || zT == 2.0*M_PI  ){InfFlag = true;}
            else{
                InfFlag = false;    x=std::tan(0.5*(xT-M_PI)); y=std::tan(0.5*(yT-M_PI)); z=std::tan(0.5*(zT-M_PI));
            }
            
        }
        
        
        // PERFORM THE MAP FROM ONE POINT COMPACTIFIED R^3 TO S^3 EMBEDDED IN R^4 //
        void InverseStereographicProjection(DOUBLE x,DOUBLE y, DOUBLE z, bool InfFlag, SU_Nc_FUNDAMENTAL_FORMAT *G){
            
            DOUBLE xS,yS,zS,uS;
            
            if(InfFlag==false){
                
                DOUBLE s=SQR(x)+SQR(y)+SQR(z);
                
                uS=(s-1.0)/(s+1.0);
                
                xS=(1.0-uS)*x;
                yS=(1.0-uS)*y;
                zS=(1.0-uS)*z;
            }
            
            else{
                xS=0.0;
                yS=0.0;
                zS=0.0;
                uS=1.0;
            }
            
            
            G[0]=xS; G[1]=yS; G[2]=zS; G[3]=uS;
            
            if(SUNcGroup::Operations::UnitarityNorm(G)>std::pow(10.0,-10)){
                std::cerr << "#ERROR -- UNITARY VIOLATION " <<  SUNcGroup::Operations::UnitarityNorm(G) << std::endl;
            }
            
            
        }
        
        void GetGaugeTransformation(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeTransformations *G){
            
            DOUBLE xT,yT,zT; DOUBLE xR,yR,zR;
            
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        bool InfFlag;
                        
                        LocalizedLatticeMap(x,y,z,xT,yT,zT);
                        TorusProjection(xT,yT,zT,xR,yR,zR,InfFlag);
                        InverseStereographicProjection(xR,yR,zR,InfFlag,G->Get(x,y,z));
                        
                    }
                }
            }
            
        }
        
        void InterpolateElectricField(INT xLow,INT xHigh,INT yLow,INT yHigh,INT zLow,INT zHigh,GaugeLinks *UOld,GaugeLinks *UNew,ElectricFields *E){
            
            for(INT z=zLow;z<=zHigh;z++){
                for(INT y=yLow;y<=yHigh;y++){
                    for(INT x=xLow;x<=xHigh;x++){
                        
                        for(INT mu=0;mu<Lattice::Dimension;mu++){
                            
                            // GET INTERPOLATING FIELD //
                            SUNcGroup::AdvancedOperations::GeodesicInterpolation(UOld->Get(x,y,z,mu),UNew->Get(x,y,z,mu),E->Get(x,y,z,mu,0));
                            
                            // SET TRANSITION TIME SCALE //
                            for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                                E->Get(x,y,z,mu,a)[0]*=1.0/tSphaleron;
                            }
                            
                        }
                        
                    }
                }
            }
            
        }
        
        void InterpolateElectricField(GaugeLinks *UOld,GaugeLinks *UNew,ElectricFields *E){
            
            InterpolateElectricField(0,E->N[0]-1,0,E->N[1]-1,0,E->N[2]-1,UOld,UNew,E);
            
        }
        
        void Setup(GaugeLinks *UInitial,GaugeLinks *UFinal,ElectricFields *E,GaugeTransformations *G){
            
            if(MPIBasic::ID==0){
                std::cerr << "#SETTING FIELDS TO ZERO" << std::endl;
            }
            
            SetZero(UInitial,E);
            
            if(MPIBasic::ID==0){
                std::cerr << "#CONSTRUCTING GAUGE TRANSFORMATION" << std::endl;
            }
            GetGaugeTransformation(0,G->N[0]-1,0,G->N[1]-1,0,G->N[2]-1,G);
                        
            DOUBLE WindingNumberIntegralEstimate=WindingNumber::IntegralEstimate::Measure(G);
            DOUBLE WindingNumber=WindingNumber::Measure(G);  
            
            if(MPIBasic::ID==0){
                std::cerr << "#WINDING NUMBER (INTEGRAL ESTIMATE) -- Nw=" << WindingNumber << " " << WindingNumberIntegralEstimate << std::endl;
            }
            GaugeTransformation::Operations::GaugeTransformLinks(UInitial,UFinal,G);

            
            DOUBLE NCsInitial=ChernSimonsNumber::VacuumEstimator::NCS(UInitial);
            DOUBLE NCsFinal=ChernSimonsNumber::VacuumEstimator::NCS(UFinal);
            
            if(MPIBasic::ID==0){
                std::cerr << "#DELTA NCS (VACUUM ESTIMATOR)-- NCs(0)=" << NCsInitial << " NCs(1)=" << NCsFinal << std::endl;
            }
            
            InterpolateElectricField(UInitial,UFinal,E);

            
            ///////////////////////////////
            // CHECK INITIAL OBSERVABLES //
            ///////////////////////////////
            
            if(MPIBasic::ID==0){
                std::cerr << "#CHECKING INITIAL STATE OBSERVABLES" << std::endl;
            }
            
            //CHECK GAUSS LAW VIOLATION //
            Observables::GaussLaw::CheckViolation();
            
            //CHECK UNITARITY VIOLATION //
            Observables::Unitarity::CheckViolation();
            
            //CHECK INITIAL ENERGY DENSITY //
            Observables::Bulk::Update();
            
            if(MPIBasic::ID==0){
                
                std::cerr << "#INITIAL ENERGY DENSITIES AND PRESSURES -- T00 TXX TYY TZZ" << std::endl;
                std::cerr << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << std::endl;
            }
                        
            
        }
        
        
    }
    
}