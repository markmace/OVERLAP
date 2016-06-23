namespace Timing{
    
    DOUBLE c0;
    DOUBLE c0W;
    
    void Reset(){
        c0=clock();
    }
    
    DOUBLE Get(){
        return (clock()-c0)/CLOCKS_PER_SEC;
    }
        
    void ResetWilsonClock(){
        c0W=clock();
    }
    
    DOUBLE GetWilsonClock(){
        return (clock()-c0W)/CLOCKS_PER_SEC;
    }
    
}