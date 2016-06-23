#ifndef __MODEPARTITION__CPP__
#define __MODEPARTITION__CPP__

namespace MPIModePartition{
    
    INT NumberOfModes;
    INT ModesPerNode;
    
    INT FirstIndex;
    
    // INITIALIZE PARTITION //
    void Init(){
        
        NumberOfModes=2*2*Nc*Lattice::N[0]*Lattice::N[1]*Lattice::N[2];
        
        ModesPerNode=NumberOfModes/MPIBasic::NumberOfNodes;
        
        if(ModesPerNode*MPIBasic::NumberOfNodes!=NumberOfModes){
            std::cerr << "#ERROR IN DISTRIBUTION OF " << NumberOfModes << " MODES ON " << MPIBasic::NumberOfNodes << " NODES" << std::endl;
            exit(0);
        }
        
        FirstIndex=MPIBasic::ID*ModesPerNode;
        
    }
    
}

#endif
