// CALCULATING -i h_W
namespace WilsonEvoOperator{
    
    // define function iup //
    INT iup(INT n1,INT n2){
        
        INT x,y,z;
        z=(INT)n1/(Lattice::N[1]*Lattice::N[0]);
        y=(INT)(n1-z*Lattice::N[1]*Lattice::N[0])/(Lattice::N[0]);
        x=n1-z*Lattice::N[1]*Lattice::N[0]-y*Lattice::N[0];
        if(n2==1) return(n1+(x+1)%(Lattice::N[0])-x);
        if(n2==2) return (n1+Lattice::N[0]*((y+1)%Lattice::N[1]-y));
        if(n2==3) return(n1+Lattice::N[0]*Lattice::N[1]*((z+1)%Lattice::N[2]-z));
        else{exit(0);}
    }
    
    // define function idn //
    INT idn(INT n1,INT n2){
        
        INT x,y,z;
        z=(INT)n1/(Lattice::N[1]*Lattice::N[0]);
        y=(INT)(n1-z*Lattice::N[1]*Lattice::N[0])/(Lattice::N[0]);
        x=n1-z*Lattice::N[1]*Lattice::N[0]-y*Lattice::N[0];
        if(n2==1) return(n1+(x-1+Lattice::N[0])%(Lattice::N[0])-x);
        if(n2==2) return (n1+Lattice::N[0]*((y-1+Lattice::N[1])%Lattice::N[1]-y));
        if(n2==3) return(n1+Lattice::N[0]*Lattice::N[1]*((z-1+Lattice::N[2])%Lattice::N[2]-z));
        else{exit(0);}
    }
    
    void Compute(u_cc *u,wcf *phi, wcf *psi){
        
        INT jf,ix,mu,jx,bc;
        DOUBLE amass;
        COMPLEX d[3][4][4];
        COMPLEX a1,a2,a3,b1,b2,b3;
        COMPLEX udag;
        
        /*Definition of mass */
        amass=6.0*rWilson+2.0*mFermion-2.0*MDWHeight;

        /* Initialization of matrices*/
        for (jf=0;jf<DiracAlgebra::SpinorComponents;jf++){
            for(ix=0;ix<Lattice::N[0]*Lattice::N[1]*Lattice::N[2];ix++){
                
                psi[ix].c[0].f[jf]=amass*phi[ix].c[0].f[jf];
                psi[ix].c[1].f[jf]=amass*phi[ix].c[1].f[jf];
                
            }
        }
        
    
        /* Constructing derivative vectors and Wilson term */
        for(ix=0;ix<Lattice::N[0]*Lattice::N[1]*Lattice::N[2];ix++){
            for(jf=0;jf<DiracAlgebra::SpinorComponents;jf++){
                for(mu=0;mu<Lattice::Dimension;mu++){
                    
                    jx=iup(ix,mu+1);
                    
                    /* bc is the bond coordinate here it is U(x)*/
                    bc=ix+mu*Lattice::N[0]*Lattice::N[1]*Lattice::N[2];
                    
                    a1=u[bc].c[0].c[0]*phi[jx].c[0].f[jf]+u[bc].c[0].c[1]*phi[jx].c[1].f[jf];
                    
                    a2=u[bc].c[1].c[0]*phi[jx].c[0].f[jf]+u[bc].c[1].c[1]*phi[jx].c[1].f[jf];
                    
                    
                    
                    /* for dagger*/
                    jx=idn(ix,mu+1);
                    /*Here it is U(x-mu)^dagger*/
                    
                    bc=jx+mu*Lattice::N[0]*Lattice::N[1]*Lattice::N[2];
                    
                    b1=conj(u[bc].c[0].c[0])*phi[jx].c[0].f[jf]+conj(u[bc].c[1].c[0])*phi[jx].c[1].f[jf];
                    
                    b2=conj(u[bc].c[0].c[1])*phi[jx].c[0].f[jf]+conj(u[bc].c[1].c[1])*phi[jx].c[1].f[jf];
                    
                    
                    d[0][jf][mu]=(a1-b1);
                    d[1][jf][mu]=(a2-b2);
                    
                    
                    psi[ix].c[0].f[jf]=psi[ix].c[0].f[jf]-a1-b1;
                    psi[ix].c[1].f[jf]=psi[ix].c[1].f[jf]-a2-b2;
                }
                
            }
            /* To get the various Dirac components of d which are the derivatives */
            /* Dirac gamma matrices are hermitian*/
            psi[ix].c[0].f[0]=psi[ix].c[0].f[0]-ComplexI*d[0][3][0]-1.0*d[0][3][1]-ComplexI*d[0][2][2];
            psi[ix].c[0].f[1]=psi[ix].c[0].f[1]-ComplexI*d[0][2][0]+1.0*d[0][2][1]+ComplexI*d[0][3][2];
            psi[ix].c[0].f[2]=psi[ix].c[0].f[2]+ComplexI*d[0][1][0]+1.0*d[0][1][1]+ComplexI*d[0][0][2];
            psi[ix].c[0].f[3]=psi[ix].c[0].f[3]+ComplexI*d[0][0][0]-1.0*d[0][0][1]-ComplexI*d[0][1][2];
            
            psi[ix].c[1].f[0]=psi[ix].c[1].f[0]-ComplexI*d[1][3][0]-1.0*d[1][3][1]-ComplexI*d[1][2][2];
            psi[ix].c[1].f[1]=psi[ix].c[1].f[1]-ComplexI*d[1][2][0]+1.0*d[1][2][1]+ComplexI*d[1][3][2];
            psi[ix].c[1].f[2]=psi[ix].c[1].f[2]+ComplexI*d[1][1][0]+1.0*d[1][1][1]+ComplexI*d[1][0][2];
            psi[ix].c[1].f[3]=psi[ix].c[1].f[3]+ComplexI*d[1][0][0]-1.0*d[1][0][1]-ComplexI*d[1][1][2];
            
            /* Multiply by -I gamma0 and halve the result*/
            psi[ix].c[0].f[0]=-ComplexI*psi[ix].c[0].f[0]*0.5;
            psi[ix].c[0].f[1]=-ComplexI*psi[ix].c[0].f[1]*0.5;
            psi[ix].c[0].f[2]=ComplexI*psi[ix].c[0].f[2]*0.5;
            psi[ix].c[0].f[3]=ComplexI*psi[ix].c[0].f[3]*0.5;
            
            psi[ix].c[1].f[0]=-ComplexI*psi[ix].c[1].f[0]*0.5;
            psi[ix].c[1].f[1]=-ComplexI*psi[ix].c[1].f[1]*0.5;
            psi[ix].c[1].f[2]=ComplexI*psi[ix].c[1].f[2]*0.5;
            psi[ix].c[1].f[3]=ComplexI*psi[ix].c[1].f[3]*0.5;
            
        }
        
        return;
    }

    
    
}


