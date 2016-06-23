namespace  IO {
    
    std::string OutputDirectory;
    
    template<typename GenericArgument>
    
    void SetOutputDirectory(GenericArgument x){
        
        OutputDirectory=StringManipulation::StringCast(x,"/");
        
        if(MPIBasic::ID==0){
            std::cerr << "#OUTPUT DIRECTORY IS " << OutputDirectory <<  std::endl;
        }
    }

}
