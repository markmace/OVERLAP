// OVERLAP HAMILTONIAN
namespace OverlapOperations{
    
    // COMPUTES SQRT(-U_WILSON^2)
    void Zolotrr(u_cc *U,DOUBLE *bz,DOUBLE *dz,COMPLEX *chi1,COMPLEX *ini,DOUBLE amin,DOUBLE amax,INT no,DOUBLE cgerr,INT mvfc,INT nitcgm){
        
        INT i,j,k;
        DOUBLE dr;
        COMPLEX zm,nr;
        DOUBLE betam;
        DOUBLE rnorm,bnorm;
        
        DOUBLE *z = new DOUBLE[no];
        COMPLEX *x = new COMPLEX[mvfc*no];
        COMPLEX *r  = new COMPLEX[mvfc];
        COMPLEX *ps  = new COMPLEX[mvfc*no];
        COMPLEX *q  = new COMPLEX[mvfc];
        COMPLEX *spbuf  = new COMPLEX[mvfc];
        DOUBLE *alphas  = new DOUBLE[no];
        DOUBLE *betas  = new DOUBLE[no];
        DOUBLE *gcoeff  = new DOUBLE[no];
        
        
        /*Initialization of chi1*/
        
        for(j=0;j<mvfc;j++){
            chi1[j]=0.0;
        }
        
        /*Multicongrad iteration*/
        for(j=0;j<mvfc;j++){
            
            r[j]=ini[j];
            for(k=0;k<no;k++){
                ps[j+k*mvfc]=ini[j];
            }
        }
        
        rnorm=norm(r,r,mvfc);
        
        i=0;
        bnorm=sqrt(rnorm);
        zm=1.0;
        betam=1.0;
        alphas[0]=0.0;
        
        for(k=0;k<no;k++){
            
            gcoeff[k]=1.0;
            z[k]=1.0;
            for(j=0;j<mvfc;j++){
                x[j+k*mvfc]=0.0;
            }
        }
        
        while(rnorm>cgerr*bnorm){
            
            for(j=0;j<mvfc;j++){
                spbuf[j]=ps[j+0*mvfc];
            }
            
            WilsonUU(U,spbuf,q);
            
            for(j=0;j<mvfc;j++){
                q[j]=q[j]+ps[j+0*mvfc]*dz[0];
            }
            
            //std::cerr << "The starting value of pAp is " << norm(spbuf,q,mvfc) << std::endl;
            //std::cerr << "The starting value of dz[0] " << dz[0] << std::endl;
            
            betas[0]=rnorm/norm(spbuf,q,mvfc);
            
            //std::cerr << "The value of beta is " << betas[0] << std::endl;
            
            
            for(j=0;j<mvfc;j++){
                r[j]=r[j]-betas[0]*q[j];
            }
            
            //std::cerr << "The value of rnorm after beta[0] is " << norm(r,r,mvfc) << std::endl;
            
            for(k=1;k<no;k++){
                
                dr=betas[0]*alphas[0]*(1.0-z[k])+betam*(1.0+(dz[k]-dz[0])*betas[0]);
                z[k]=betam/dr;
                betas[k]=betas[0]*z[k];
                
            }
            
            alphas[0]=norm(r,r,mvfc)/rnorm;
            rnorm=norm(r,r,mvfc);
            
            //std::cerr << "The value of rnorm after alphas[0] is " << rnorm << std::endl;
            
            for(k=0;k<no;k++){
                
                for(j=0;j<mvfc;j++){
                    x[j+k*mvfc]=x[j+k*mvfc]+betas[k]*ps[j+k*mvfc];
                }
                
                alphas[k]=alphas[0]*z[k]*betas[k]/betas[0];
                gcoeff[k]=gcoeff[k]*z[k];
                
                for(j=0;j<mvfc;j++){
                    ps[j+k*mvfc]=ps[j+k*mvfc]*alphas[k]+r[j]*gcoeff[k];
                }
                
            }
            
            betam=betas[0];
            i++;
            
            if(i>=nitcgm){
                
                std::cerr << "Too many CG iteration steps, res=" << rnorm << std::endl;
                break;
            }
            
            
            
        }
        
        // OPTION TO OUTPUT THE TOTAL NUMBER OF CONJUGATE GRADIENT ITERATIONS //
        //if(MPIBasic::ID==0){
        //    std::cerr << "The total no of cg iterations are " << i << std::endl;
        //}
        // END OPTION
        
        for(k=0;k<no;k++){
            
            for(j=0;j<mvfc;j++)
                chi1[j]=chi1[j]+x[j+k*mvfc]*bz[k];
        }
        
        
        
        delete x;
        delete z;
        delete spbuf;
        delete ps;
        delete q;
        delete r;
        delete alphas;
        delete betas;
        delete gcoeff;
        
        return;
        
    }

    // COMPUTES U_OV
    void OverlapU(u_cc *u,wcf *phi,wcf *psi,DOUBLE *bz,DOUBLE *dz,INT no,INT mvfc,DOUBLE amin,DOUBLE amax,INT nitcgm,DOUBLE cgerr){
        
        INT i,k,j,m;
        COMPLEX *chi2=new COMPLEX[mvfc];
        COMPLEX *ini =new COMPLEX[mvfc];
        wcf *ovi= new wcf[Lattice::N[0]*Lattice::N[1]*Lattice::N[2]];
        
        for(i=0;i<Lattice::N[0]*Lattice::N[1]*Lattice::N[2];i++){
            for(k=0;k<DiracAlgebra::SpinorComponents;k++){
                for(j=0;j<Nc;j++){
                    
                    m=i+Lattice::N[0]*Lattice::N[1]*Lattice::N[2]*(k+j*DiracAlgebra::SpinorComponents);
                    
                    ini[m]=phi[i].c[j].f[k];
                    
                }
            }
        }
        
        // COMPUTES SQRT(-U_WILSON^2)
        Zolotrr(u,bz,dz,chi2,ini,amin,amax,no,cgerr,mvfc,nitcgm);
                
        for(i=0;i<Lattice::N[0]*Lattice::N[1]*Lattice::N[2];i++){
            for(k=0;k<DiracAlgebra::SpinorComponents;k++){
                for(j=0;j<Nc;j++){
                    
                    m=i+Lattice::N[0]*Lattice::N[1]*Lattice::N[2]*(k+j*DiracAlgebra::SpinorComponents);
                    
                    ovi[i].c[j].f[k]=chi2[m];
                }
            }
        }
        
        // COMPUTES U_W
        WilsonEvoOperator::Compute(u,ovi,psi);
        
        for(INT ix=0;ix<Lattice::N[0]*Lattice::N[1]*Lattice::N[2];ix++){
            
            // Multiply by -I gamma0 and halve the fermion mass //
            psi[ix].c[0].f[0]=MDWHeight*(-ComplexI*(1.0+0.5*mFermion)*phi[ix].c[0].f[0]+(1.0-0.5*mFermion)*psi[ix].c[0].f[0]);
            psi[ix].c[0].f[1]=MDWHeight*(-ComplexI*(1.0+0.5*mFermion)*phi[ix].c[0].f[1]+(1.0-0.5*mFermion)*psi[ix].c[0].f[1]);
            psi[ix].c[0].f[2]=MDWHeight*(ComplexI*(1.0+0.5*mFermion)*phi[ix].c[0].f[2]+(1.0-0.5*mFermion)*psi[ix].c[0].f[2]);
            psi[ix].c[0].f[3]=MDWHeight*(ComplexI*(1.0+0.5*mFermion)*phi[ix].c[0].f[3]+(1.0-0.5*mFermion)*psi[ix].c[0].f[3]);
            
            
            psi[ix].c[1].f[0]=MDWHeight*(-ComplexI*(1.0+0.5*mFermion)*phi[ix].c[1].f[0]+(1.0-0.5*mFermion)*psi[ix].c[1].f[0]);
            psi[ix].c[1].f[1]=MDWHeight*(-ComplexI*(1.0+0.5*mFermion)*phi[ix].c[1].f[1]+(1.0-0.5*mFermion)*psi[ix].c[1].f[1]);
            psi[ix].c[1].f[2]=MDWHeight*(ComplexI*(1.0+0.5*mFermion)*phi[ix].c[1].f[2]+(1.0-0.5*mFermion)*psi[ix].c[1].f[2]);
            psi[ix].c[1].f[3]=MDWHeight*(ComplexI*(1.0+0.5*mFermion)*phi[ix].c[1].f[3]+(1.0-0.5*mFermion)*psi[ix].c[1].f[3]);
            
        }
        
        delete ini;
        delete ovi;
        delete chi2;
        
    }
}


