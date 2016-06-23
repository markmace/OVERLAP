namespace OverlapOperations{
    
    void OverlapSignTest(u_cc *Ug,DOUBLE *bz,DOUBLE *dz,COMPLEX *ini,DOUBLE amin,DOUBLE amax,DOUBLE cgerr,INT no,INT mvfc,INT nitcgm){
        //void overlapgwD(u_cc *Ug,double *bz,double *dz,std::complex<double> *fnal,std::complex<double> *ini,double amin,double amax,double cgerr,int no,int mvfc,int nitcgm){
        
        INT i;
        COMPLEX *chi2,*fnal2,*fnal3;
        wcf *ovi,*ovf;
        ovi= new wcf[Lattice::Volume];
        ovf=new wcf[Lattice::Volume];
        chi2=new COMPLEX[mvfc];
        fnal2=new COMPLEX[mvfc];
        fnal3=new COMPLEX[mvfc];
        
        /*
         
         
         zolotrr(Ug,bz,dz,chi2,ini,amin,amax,no,cgerr,mvfc,nitcgm);
         
         for(i=0;i<mv;i++)
         {
         for(k=0;k<nf;k++)
         {
         for(j=0;j<nc;j++)
         {
         m=i+mv*k+j*mv*nf;
         
         ovi[i].c[j].f[k]=chi2[m];
         
         
         
         }
         }
         }
         
         WilsonHamiltonianOperator::Compute(Ug,ovi,ovf);
         
         for(i=0;i<mv;i++)
         {
         for(k=0;k<nf;k++)
         {
         for(j=0;j<nc;j++)
         {
         m=i+mv*k+j*mv*nf;
         
         fnal[m]=ovf[i].c[j].f[k];
         
         }
         }
         }
         
         
         
         
         for(i=0;i<mvfc;i++)
         fnal[i]=ini[i]+fnal[i];
         
         
         
         // To check Ginsparg-Wilson relation I g0 D+D I g0=-D D
         
         
         
         
         
         //Make D g5|ini>
         
         zolotrr(Ug,bz,dz,chi2,ini,amin,amax,no,cgerr,mvfc,nitcgm);
         
         for(i=0;i<mv;i++)
         {
         for(k=0;k<nf;k++)
         {
         for(j=0;j<nc;j++)
         {
         m=i+mv*k+j*mv*nf;
         
         ovi[i].c[j].f[k]=chi2[m];
         
         
         
         }
         }
         }
         
         WilsonHamiltonianOperator::Compute(Ug,ovi,ovf);
         
         for(i=0;i<mv;i++)
         {
         for(k=0;k<nf;k++)
         {
         for(j=0;j<nc;j++)
         {
         m=i+mv*k+j*mv*nf;
         
         fnal2[m]=ovf[i].c[j].f[k];
         
         
         }
         }
         }
         
         
         
         for(i=0;i<mvfc;i++)
         fnal2[i]=ini[i]+fnal2[i];
         
         
         //To make D  D|ini>
         zolotrr(Ug,bz,dz,chi2,fnal,amin,amax,no,cgerr,mvfc,nitcgm);
         
         for(i=0;i<mv;i++)
         {
         for(k=0;k<nf;k++)
         {
         for(j=0;j<nc;j++)
         {
         m=i+mv*k+j*mv*nf;
         
         ovi[i].c[j].f[k]=chi2[m];
         
         
         
         }
         }
         }
         
         WilsonHamiltonianOperator::Compute(Ug,ovi,ovf);
         
         for(i=0;i<mv;i++)
         {
         for(k=0;k<nf;k++)
         {
         for(j=0;j<nc;j++)
         {
         m=i+mv*k+j*mv*nf;
         
         fnal3[m]=ovf[i].c[j].f[k];
         
         
         }
         }
         }
         
         
         
         for(i=0;i<mvfc;i++)
         {
         fnal3[i]=fnal[i]+fnal3[i];
         fnal3[i]=fnal[i]+fnal2[i]+fnal3[i];
         }
         */
        //printf("The deviation from GW relation is %0.9le\n",sqrtf(norm(fnal3,fnal3)));
        
        /* to calculate epsilon*/
        Zolotrr(Ug,bz,dz,fnal2,ini,amin,amax,no,cgerr,mvfc,nitcgm);
        Zolotrr(Ug,bz,dz,chi2,fnal2,amin,amax,no,cgerr,mvfc,nitcgm);
        WilsonUU(Ug,chi2,fnal3);
        for(i=0;i<mvfc;i++){
            fnal3[i]=fnal3[i]-ini[i];
        }
        if(MPIBasic::ID==0){
            
            std::cerr << "The value of epsilon is is " << sqrt(norm(fnal3,fnal3,mvfc)) << std::endl;;
        }
        delete[] chi2;
        delete[] ovi;
        delete[] ovf;
        delete[] fnal3;
        delete[] fnal2;
        
        
        
    }
}


/*namespace OverlapOperations{
 
 #include "../../WILSON/OPERATOR/WilsonOperator.cpp"
 
 void OverlapgSignFunctionTest(u_cc *u,COMPLEX *fnal,COMPLEX *ini,DOUBLE amin,DOUBLE amax,DOUBLE cgerr,DOUBLE fm,INT no,INT mvfc){
 
 COMPLEX *chi2,*fnal2,*fnal3;
 wcf *ovi,*ovf,*ovf1,*ovf2;
 ovi=new wcf[Lattice::N[0]*Lattice::N[1]*Lattice::N[2]];
 ovf=new wcf[Lattice::N[0]*Lattice::N[1]*Lattice::N[2]];
 chi2=new COMPLEX[mvfc];
 fnal2=new COMPLEX[mvfc];
 fnal3=new COMPLEX[mvfc];
 
 //std::cerr << "The initial norm is " << norm(ini,ini)) << std::endl; //

Zolotrr(chi2,ini,amin,amax,no,errtol,cgerr);

for(INT i=0;i<mv;i++){
    for(INT k=0;k<nf;k++){
        for(INT j=0;j<nc;j++){
            
            INT m=i+mv*k+j*mv*nf;
            
            ovi[i].c[j].f[k]=chi2[m];
            
            
        }
    }
}

WilsonEvoOperator::Compute(u,ovi,ovf);

for(INT i=0;i<mv;i++){
    for(k=0;k<nf;k++){
        for(j=0;j<nc;j++){
            
            INT m=i+mv*k+j*mv*nf;
            
            fnal[m]=ovf[i].c[j].f[k];
            
        }
    }
}




for(i=0;i<mvfc;i++){
    fnal[i]=ini[i]+fnal[i];
}



// To check Ginsparg-Wilson relation I g0 D+D I g0=-D D //


// Make D g5|ini> //

Zolotrr(u,bz,dz,chi2,ini,amin,amax,no,cgerr,mvfc,nitcgm);

for(INT i=0;i<mv;i++){
    for(INT k=0;k<nf;k++){
        for(INT j=0;j<nc;j++){
            
            INT m=i+mv*k+j*mv*nf;
            
            ovi[i].c[j].f[k]=chi2[m];
            
            
            
        }
    }
}

WilsonEvoOperator::Compute(u,ovi,ovf);

for(INT i=0;i<mv;i++){
    for(INT k=0;k<nf;k++){
        for(INT j=0;j<nc;j++){
            
            INT m=i+mv*k+j*mv*nf;
            
            fnal2[m]=ovf[i].c[j].f[k];
            
        }
    }
}



for(INT i=0;i<mvfc;i++){
    fnal2[i]=ini[i]+fnal2[i];
}


// To make D  D|ini> //
//        Zolotrr(u,bz,dz,chi2,ini,amin,amax,no,cgerr,mvfc,nitcgm);

Zolotrr(u,bz,dz,chi2,fnal,amin,amax,no,errtol,cgerr);

for(INT i=0;i<mv;i++){
    for(INT k=0;k<nf;k++){
        for(INT j=0;j<nc;j++){
            
            INT m=i+mv*k+j*mv*nf;
            
            ovi[i].c[j].f[k]=chi2[m];
            
            
        }
    }
}

WilsonEvoOperator::Compute(u,ovi,ovf);

for(INT i=0;i<mv;i++){
    for(INT k=0;k<nf;k++){
        for(INT j=0;j<nc;j++){
            
            INT m=i+mv*k+j*mv*nf;
            
            fnal3[m]=ovf[i].c[j].f[k];
            
            
        }
    }
}



for(INT i=0;i<mvfc;i++){
    
    fnal3[i]=fnal[i]+fnal3[i];
    fnal3[i]=fnal[i]+fnal2[i]+fnal3[i];
}

std::cerr << "The deviation from GW relation is " << sqrtf(norm(fnal3,fnal3))) std::endl;;

// to calculate epsilon //
Zolotrr(u,bz,dz,fnal2,ini,amin,amax,no,cgerr,mvfc,nitcgm);
Zolotrr(u,bz,dz,chi2,fnal2,amin,amax,no,cgerr,mvfc,nitcgm);
WilsonUU(u,chi2,fnal3);
for(INT i=0;i<mvfc;i++){
    fnal3[i]=fnal3[i]-ini[i];
}
std::cerr << "The value of epsilon is is " << std::sqrt(norm(fnal3,fnal3))) << std::endl;

delete[] chi2;
delete[] ovi;
delete[] ovf;
delete[] fnal3;
delete[] fnal2;



}
}
*/
