
namespace OverlapOperations{
    
    // DETERMINE MINIMUM AND MAXIMUM EIGENVALUES OF HAMILTONIAN VIA RITZ ALGORITHM
    
    void orthog(COMPLEX *x1,COMPLEX *x2,INT length,INT neigen,INT nmin, INT nmax){
        /* orthogonalise x2 wrt x1*/
        INT i,j;
        COMPLEX val;
        
        for(i=0;i<neigen-1;i++){
            
            val=0.0;
            for(j=nmin-1;j<nmax;j++){
                val=val+x2[j]*conj(x1[j-nmin+1+length*i]);
            }
            
            for(j=nmin-1;j<nmax;j++){
                x2[j]=x2[j]-x1[j-nmin+1+length*i]*val;
            }
        }
        
    }
    
    COMPLEX complexnorm(COMPLEX *vector1,COMPLEX *vector2,INT mvfc){
        
        COMPLEX  tot=0.0;
        
        for(INT i=0;i<mvfc;i++){
            tot+=vector2[i]*conj(vector1[i]);
        }
        
        return tot;
    }
    
    DOUBLE norm(COMPLEX *vector1,COMPLEX *vector2,INT mvfc){
        
        COMPLEX tot=0.0;
        
        for(INT i=0;i<mvfc;i++){
            tot+=vector2[i]*conj(vector1[i]);
        }
        
        return real(tot);
    }
    
#include "../../WILSON/OPERATOR/WilsonOperator.cpp"
    
    void WilsonUU(u_cc *Ug,COMPLEX *a,COMPLEX *b){
        
        wcf *chi = new wcf[Lattice::N[0]*Lattice::N[1]*Lattice::N[2]];
        wcf *A = new wcf[Lattice::N[0]*Lattice::N[1]*Lattice::N[2]];
        wcf *B = new wcf[Lattice::N[0]*Lattice::N[1]*Lattice::N[2]];
        INT m;
        
        for(INT i=0;i<Lattice::N[0]*Lattice::N[1]*Lattice::N[2];i++){
            for(INT k=0;k<DiracAlgebra::SpinorComponents;k++){
                for(INT j=0;j<Nc;j++){
                    
                    m=i+Lattice::N[0]*Lattice::N[1]*Lattice::N[2]*(k+j*DiracAlgebra::SpinorComponents);
                    
                    A[i].c[j].f[k]=a[m];
                    
                }
            }
        }
        
        WilsonEvoOperator::Compute(Ug,A,chi);
        
        WilsonEvoOperator::Compute(Ug,chi,B);
        
        for(INT i=0;i<Lattice::N[0]*Lattice::N[1]*Lattice::N[2];i++){
            for(INT k=0;k<DiracAlgebra::SpinorComponents;k++){
                for(INT j=0;j<Nc;j++){
                    
                    m=i+Lattice::N[0]*Lattice::N[1]*Lattice::N[2]*(k+j*DiracAlgebra::SpinorComponents);
                    
                    // NEGATIVE SIGN BECAUSE ARGUEMENT OF SQRT IS -U_evo^Wilson ^2 //
                    b[m]=-B[i].c[j].f[k];
                    
                }
            }
        }
        
        delete chi;
        delete A;
        delete B;
        
    }
    
    
    DOUBLE RitzSolver(u_cc *u,COMPLEX *vector,INT neigen,DOUBLE delitr,INT nitlz,DOUBLE sign,INT nrestart,INT mvfc){
        
        
        DOUBLE omeg1;
        
        INT nmin,nmax;
        INT i1,i2,i3,i4,i8,itt,ittot,length,istart;
        DOUBLE xnorm,gnorm,pnorm,gnormim;
        COMPLEX xpnorm,xynorm;
        DOUBLE beta,ct,st;
        COMPLEX amb,oxynorm;
        COMPLEX *xr=new COMPLEX[mvfc*neigen];
        COMPLEX *yr=new COMPLEX[mvfc*neigen];
        COMPLEX *zr=new COMPLEX[mvfc*neigen];
        COMPLEX *gr=new COMPLEX[mvfc*neigen];
        COMPLEX *pr=new COMPLEX[mvfc*neigen];
        
        DOUBLE ccc,cox,cox2,sox,sqp,sqm,sig1,sig2;
        INT cond;
        
        // Assuming ichiral=1 //
        nmin=1;
        nmax=mvfc;
        length=nmax-nmin+1;
        
        for(i1=nmin-1;i1<nmax;i1++){
            xr[i1]=vector[i1-nmin+1+length*(neigen-1)];
        }
        
        
        orthog(vector,xr,length,neigen,nmin,nmax);
        xnorm=norm(xr,xr,mvfc);
        xnorm=1.0/sqrt(xnorm);
        
        
        for(i3=nmin-1;i3<nmax;i3++){
            xr[i3]=xr[i3]*xnorm;
        }

        ittot=0;
        istart=-1;
        cond=0;
        /* option for restart*/
        while(cond>=0){
            
            istart=istart+1;
            
            WilsonUU(u,xr,yr);
            
            for(i2=nmin-1;i2<nmax;i2++){
                yr[i2]=sign*yr[i2];
            }
            
            //To orthogonalize//
            orthog(vector,yr,length,neigen,nmin,nmax);
            
            xynorm=complexnorm(xr,yr,mvfc);
            
            //std::cerr << "#xynorm " << xynorm << std::endl;
            
            for(i2=nmin-1;i2<nmax;i2++){
                
                gr[i2]=yr[i2]-xynorm*xr[i2];
                pr[i2]=gr[i2];
                
            }
            
            gnormim=norm(gr,gr,mvfc);
            orthog(vector,pr,length,neigen,nmin,nmax);
            pnorm=sqrt(gnormim);
            
            
            for(i2=nmin-1;i2<nmax;i2++){
                pr[i2]=pr[i2]/pnorm;
            }
            
            for(itt=1;itt<=nrestart;itt++){
                
                ittot=ittot+1;
                //std::cerr << "#xynorm " << xynorm << std::endl;
                
                WilsonUU(u,pr,zr);
                
                for(i2=nmin-1;i2<nmax;i2++){
                    zr[i2]=sign*zr[i2];
                }
                
                orthog(vector,zr,length,neigen,nmin,nmax);
                amb=xynorm;
                amb=amb-norm(pr,zr,mvfc);
                ccc=real(complexnorm(xr,zr,mvfc));
                cox=2.0*ccc/real(amb);
                
                /*	    theta=0.5*atan(cox);
                 ct=cos(theta);
                 st=sin(theta);
                 if(ct*st*ccc>0)
                 {
                 theta=theta-0.5*PI;
                 ct=cos(theta);
                 st=sin(theta);
                 }*/
                sig1=-std::copysign(1.0,real(amb));
                sig2=std::copysign(1.0,cox)*sig1;
                cox2=cox*cox;
                if(cox2<=1.0e-3){
                    
                    sqp=1.0-0.125*cox2+0.0859375*cox2*cox2-0.0673828125*cox2*cox2*cox2+0.0562439*cox2*cox2*cox2*cox2+0.0487*cox2*cox2*cox2*cox2*cox2;
                    sqm=0.5*sig1*sig2*cox*(1.0-0.375*cox2+0.2421875*cox2*cox2-0.182617875*cox2*cox2*cox2+0.1482849*cox2*cox2*cox2*cox2-0.12575*cox2*cox2*cox2*cox2*cox2+0.1097*cox2*cox2*cox2*cox2*cox2*cox2);
                    
                    if(sig1<0.0){
                        ct=sqm;
                        st=sig2*sqp;
                    }
                    else{
                        ct=sqp;
                        st=sig2*sqm;
                    }
                }
                else{
                    sox=sig1/sqrt(1.0+cox2);
                    ct=sqrt((1.0+sox)*0.5);
                    st=sqrt((1.0-sox)*0.5)*sig2;
                }
                
                
                for(i8=nmin-1;i8<nmax;i8++)
                {
                    xr[i8]=xr[i8]*ct+pr[i8]*st;
                    yr[i8]=yr[i8]*ct+zr[i8]*st;
                }
                
                orthog(vector,xr,length,neigen,nmin,nmax);
                orthog(vector,yr,length,neigen,nmin,nmax);
                xnorm=norm(xr,xr,mvfc);
                
                oxynorm=xynorm;
                xynorm=complexnorm(xr,yr,mvfc);
                xnorm=1.0/sqrt(xnorm);
                xynorm=xynorm*xnorm;
                
                
                for(i4=nmin-1;i4<nmax;i4++){
                    
                    xr[i4]=xr[i4]*xnorm;
                    gr[i4]=yr[i4]-xynorm*xr[i4];
                }
                
                gnorm=norm(gr,gr,mvfc);
                
                
                xpnorm=complexnorm(xr,pr,mvfc);
                
                if(real(xynorm)>real(oxynorm)){
                    cond=-1;
                    break;
                }
                
                if(gnorm<delitr*sign*real(xynorm))
                {
                    cond=-1;
                    break;
                }
                
                beta=gnorm/gnormim*pnorm;
                /*beta=gnorm/gnormim*ct;
                 if(beta>2.0) beta=2.0;*/
                
                for(i4=nmin-1;i4<nmax;i4++){
                    pr[i4]=gr[i4]+beta*(pr[i4]-xr[i4]*xpnorm);
                }
                
                orthog(vector,pr,length,neigen,nmin,nmax);
                
                pnorm=sqrt(norm(pr,pr,mvfc));
                
                for(i4=nmin-1;i4<nmax;i4++){
                    pr[i4]=pr[i4]/pnorm;
                }
                
                gnormim=gnorm;
                
                if(ittot>=nitlz){
                    
                    std::cerr << "NUMBER OF ITERATIONS EXCEED THE MAXIMUM ALLOWED WITH " << ittot << " RESTARTS" << std::endl;
                    cond=-1;
                    break;
                }
                
            }
            
        }
        
        omeg1=real(xynorm);
        
        //std::cerr << "The value of eigenavlue is " << omeg1 << " restart is " << ittot << " " << istart << std::endl;

        delete[] xr;
        delete[] yr;
        delete[] zr;
        delete[] pr;
        delete[] gr;
        
        return omeg1;
        
    }
}


