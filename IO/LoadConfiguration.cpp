namespace IO{
    
    void LoadConfiguration(GaugeLinks *U,ElectricFields *E){
        
        //NUMBER OF LINES READ FROM FILE
        INT InputCount=0;
        
        // INPUT STREAMS //
        std::ifstream UInStream,EInStream;
        
        // INPUT FILES //
        std::string UInputFile=StringManipulation::StringCast(IO::InputDirectory,IO::InputUFile);
        std::string EInputFile=StringManipulation::StringCast(IO::InputDirectory,IO::InputEFile);
        
        // OPEN FILES //
        UInStream.open(UInputFile.c_str()); EInStream.open(EInputFile.c_str());
        
        // SET PRECISION //
        UInStream.precision(OUTPUT_PRECISION); EInStream.precision(OUTPUT_PRECISION);
        
        std::string ULink;
        
        // INPUT DATA //
        INT x; INT y; INT z; INT mu;
        
        // BUFFER FOR STRING TO MATRIX //
        std::string UStr;
        
        // BUFFER FOR ELECTRIC FIELD //
        DOUBLE E0; DOUBLE E1; DOUBLE E2;
        
        #if SU_Nc_FLAG==SU3_FLAG
        
        DOUBLE E3; DOUBLE E4; DOUBLE E5; DOUBLE E6; DOUBLE E7;
        
        #endif
        
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#READING IN GAUGE LINKS FROM FILE " << UInputFile << std::endl;
        
        // MONITOR PRECISION //
        //DOUBLE MaxUnitarityViolation=DOUBLE(0.0);
        
        //GET GAUGE LINK DATA FROM INPUT FILES
        while(UInStream.good()){
            
            //READ LINES
            getline(UInStream,ULink);
            
            //PROCESS FILE LINE BY LINE
            if(!(ULink.empty())){
                
                //STRING TOKEN
                std::stringstream ULinkValues(ULink);
                
                //GET POSITIONS IN FILE
                ULinkValues >> x; ULinkValues >> y; ULinkValues >> z; ULinkValues >> mu;
                
                std::stringstream UString;    std::string strBuff;
                
                while(ULinkValues >> strBuff){UString << strBuff << " ";}
                
                SUNcGroup::IO::StringToMatrix(UString.str(),U->Get(x,y,z,mu));

                
                //INCREASE POSITION COUNT
                InputCount++;
            }
        }
        
        if((InputCount)!=Lattice::N[0]*Lattice::N[1]*Lattice::N[2]*Lattice::Dimension){
            std::cerr << "#GAUGE LINK INPUT CAUSED FATAL ERROR -- GAUGE LINK NOT LOADED CORRECTLY" << std::endl;
            std::cerr << "# " << InputCount << " "  << Lattice::N[0]*Lattice::N[1]*Lattice::N[2]*Lattice::Dimension << std::endl;

            exit(0);
        }
        
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#FINISHED READING IN GAUGE LINKS FROM FILE " << UInputFile << std::endl;
        
        // RESET
        InputCount=0;
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#READING IN ELECTRIC FIELDS FROM FILE " << EInputFile << std::endl;
        
        // INPUT ELECTRIC FIELDS //
        if (EInStream.is_open()){
            while(!EInStream.eof()){
                
                #if SU_Nc_FLAG==SU2_FLAG
                // INPUT ELECTRIC FIELDS //
                EInStream >> x >> y >> z >> mu >> E0 >> E1 >> E2;
                
                E->Get(x,y,z,mu,0)[0]=E0;
                E->Get(x,y,z,mu,1)[0]=E1;
                E->Get(x,y,z,mu,2)[0]=E2;
                #endif
                
                #if SU_Nc_FLAG==SU3_FLAG
                // INPUT ELECTRIC FIELDS //
                EInStream >> x >> y >> z >> mu >> E0 >> E1 >> E2 >> E3 >> E4 >> E5 >> E6 >> E7;
                
                E->Get(x,y,z,mu,0)[0]=E0;
                E->Get(x,y,z,mu,1)[0]=E1;
                E->Get(x,y,z,mu,2)[0]=E2;
                E->Get(x,y,z,mu,3)[0]=E3;
                E->Get(x,y,z,mu,4)[0]=E4;
                E->Get(x,y,z,mu,5)[0]=E5;
                E->Get(x,y,z,mu,6)[0]=E6;
                E->Get(x,y,z,mu,7)[0]=E7;
                
                #endif
                
                //INCREASE POSITION COUNT
                InputCount++;
                
            }
            
        }
        
        if((InputCount-1)!=Lattice::N[0]*Lattice::N[1]*Lattice::N[2]*Lattice::Dimension){
            std::cerr << "#ELECTRIC FIELD INPUT CAUSED FATAL ERROR -- ELECTRIC FIELD NOT LOADED CORRECTLY" << std::endl;
            std::cerr << "# " << InputCount << " "  << Lattice::N[0]*Lattice::N[1]*Lattice::N[2]*Lattice::Dimension << std::endl;

            exit(0);
        }
        
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#FINISHED READING IN ELECTRIC FIELDS FROM FILE " << EInputFile << std::endl;
        
        // CLOSE INPUT STREAM //
        UInStream.close(); EInStream.close();
        
        // MEASURE BULK OBSERVABLES //
        Observables::Bulk::Update();
        // OUTPUT
        std::cerr << "#INITIAL ENERGY DENSITIES AND PRESSURES -- T00 TXX TYY TZZ" << std::endl;
        std::cerr << Observables::Bulk::T00() << " " << Observables::Bulk::TXX() << " " << Observables::Bulk::TYY() << " " << Observables::Bulk::TZZ() << std::endl;
        
    }

}


