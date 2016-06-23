namespace IO{
    
    void SaveConfiguration(std::string Ufname,std::string Efname){
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#CONFIGURATION SAVED AT T=" << Dynamics::Time() << std::endl;
        
        // OUTPUT STREAMS //
        std::ofstream UOutStream,EOutStream;
        
        // OUTPUT FILES //
        std::string UOutputFile=StringManipulation::StringCast(IO::OutputDirectory,Ufname,"ID",RandomNumberGenerator::MySEED,".txt");
        std::string EOutputFile=StringManipulation::StringCast(IO::OutputDirectory,Efname,"ID",RandomNumberGenerator::MySEED,".txt");

        // OPEN FILES //
        UOutStream.open(UOutputFile.c_str()); EOutStream.open(EOutputFile.c_str());
        
        // SET PRECISION //
        UOutStream.precision(OUTPUT_PRECISION); EOutStream.precision(OUTPUT_PRECISION);
        
        // CREATE OUTPUT //
        for(INT z=0;z<GLinks::U->N[2];z++){
            for(INT y=0;y<GLinks::U->N[1];y++){
                for(INT x=0;x<GLinks::U->N[0];x++){
                    for(INT mu=0;mu<Lattice::Dimension;mu++){
                    
                        // OUTPUT GAUGE LINKS //
                        // STANDARD OUTPUT //
                        UOutStream << x << " " << y << " " << z << " " << mu << " " << SUNcGroup::IO::MatrixToString(GLinks::U->Get(x,y,z,mu)) << std::endl;
                        // END STANDARD GAUGE LINK OUTPUT //
                        
                        // BEGIN SAYANTAN OUTPUT //
                        /*
                        COMPLEX UMat[Nc*Nc];
                        SUNcGroup::Operations::GetMatrix(GLinks::U->Get(x,y,z,mu),UMat);
                    
                        UOutStream << std::real(UMat[0]) <<  " " << std::imag(UMat[0]) << std::endl;
                        UOutStream << std::real(UMat[2]) <<  " " << std::imag(UMat[2]) << std::endl;
                        UOutStream << std::real(UMat[1]) <<  " " << std::imag(UMat[1]) << std::endl;
                        UOutStream << std::real(UMat[3]) <<  " " << std::imag(UMat[3]) << std::endl;
                        */
                        // END SAYANTAN OUTPUT //
                        

                        // OUTPUT ELECTRIC FIELDS //
                        EOutStream << x << " " << y << " " << z << " " << mu;
                        
                        for(INT a=0;a<SUNcAlgebra::VectorSize;a++){
                            EOutStream << " " << EFields::E->Get(x,y,z,mu,a)[0];
                        }
                        
                        EOutStream << std::endl;
                        
                    }
                    
                }
            }
        }
        
        // CLOSE OUTPUT STREAM //
        UOutStream.close(); EOutStream.close();
        
    }
    void SaveGaugeLinks(std::string Ufname){
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#GAUGE LINKS SAVED AT T=" << Dynamics::Time() << std::endl;
        
        // OUTPUT STREAMS //
        std::ofstream UOutStream;
        
        // OUTPUT FILES //
        std::string UOutputFile=StringManipulation::StringCast(IO::OutputDirectory,Ufname,"ID",RandomNumberGenerator::MySEED,".txt");
        
        // OPEN FILES //
        UOutStream.open(UOutputFile.c_str());
        
        // SET PRECISION //
        UOutStream.precision(OUTPUT_PRECISION);
        
        // CREATE OUTPUT //
        for(INT z=0;z<GLinks::U->N[2];z++){
            for(INT y=0;y<GLinks::U->N[1];y++){
                for(INT x=0;x<GLinks::U->N[0];x++){
                    for(INT mu=0;mu<Lattice::Dimension;mu++){
                        
                        // OUTPUT GAUGE LINKS //
                        // STANDARD OUTPUT //
                        UOutStream << x << " " << y << " " << z << " " << mu << " " << SUNcGroup::IO::MatrixToString(GLinks::U->Get(x,y,z,mu)) << std::endl;
                        // END STANDARD GAUGE LINK OUTPUT //
                        
                        // BEGIN SAYANTAN OUTPUT //
                        /*
                         COMPLEX UMat[Nc*Nc];
                         SUNcGroup::Operations::GetMatrix(GLinks::U->Get(x,y,z,mu),UMat);
                         
                         UOutStream << std::real(UMat[0]) <<  " " << std::imag(UMat[0]) << std::endl;
                         UOutStream << std::real(UMat[2]) <<  " " << std::imag(UMat[2]) << std::endl;
                         UOutStream << std::real(UMat[1]) <<  " " << std::imag(UMat[1]) << std::endl;
                         UOutStream << std::real(UMat[3]) <<  " " << std::imag(UMat[3]) << std::endl;
                         */
                        // END SAYANTAN OUTPUT //
                        
                        
                        
                    }
                    
                }
            }
        }
        
        // CLOSE OUTPUT STREAM //
        UOutStream.close();
        
    }

    
    
}
