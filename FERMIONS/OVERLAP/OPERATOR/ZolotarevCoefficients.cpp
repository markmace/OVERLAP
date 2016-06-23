
namespace OverlapOperations{
    
    // COMPUTES ZOLOTAREZ COEFFICIENTS //
    DOUBLE CompElliptic(DOUBLE Kappa,DOUBLE err){
        
        DOUBLE mu,Xn,Yn,Zn,xndev,yndev,zndev,epsilon;
        DOUBLE e2,e3,s,k,xnroot,ynroot,znroot,lambda;
        // SETTING DUMMY INITIAL VALUES
        Xn=0.0;
        Yn=1.0-Kappa;
        Zn=1.0;
        epsilon=err;
        
        mu=(Xn+Yn+Zn)/3.0;
        
        xndev=2.0-(mu+Xn)/mu;
        yndev=2.0-(mu+Yn)/mu;
        zndev=2.0-(mu+Zn)/mu;
        
        while(epsilon>=err){
            
            mu=(Xn+Yn+Zn)/3.0;
            
            xndev=2.0-(mu+Xn)/mu;
            yndev=2.0-(mu+Yn)/mu;
            zndev=2.0-(mu+Zn)/mu;
            
            epsilon=fabs(xndev);
            if(fabs(yndev)>epsilon) epsilon=fabs(yndev);
            if(fabs(zndev)>epsilon) epsilon=fabs(zndev);
            
            /* the implementation for  epsilon=max(abs(xndev),abs(yndev),abs(zndev));*/
            xnroot=sqrt(Xn);
            ynroot=sqrt(Yn);
            znroot=sqrt(Zn);
            lambda=xnroot*(ynroot+znroot)+ynroot*znroot;
            Xn=0.25*(Xn+lambda);
            Yn=0.25*(Yn+lambda);
            Zn=0.25*(Zn+lambda);
        }
        e2=xndev*yndev-zndev*zndev;
        e3=xndev*yndev*zndev;
        s=1.0+(e2/24.0-0.1-e3/44.0)*e2+e3/14.0;
        k=s/sqrt(mu);
        
        return k;
        
    }
    
    DOUBLE coeffb(DOUBLE *cval, INT n,INT norder){
        
        DOUBLE valb,prodnr;
        INT j=0;
        prodnr=1.0;
        
        for(INT i=0;i<norder-1;i++){
            prodnr=prodnr*(cval[2*i+1]-cval[2*n]);
            if(i==n){
                j=j+1;
            }
            prodnr=prodnr/(cval[2*j]-cval[2*n]);
            j++;
            
        }
        
        valb=prodnr;
        return valb;
    }
    
    DOUBLE ZolotarevCoefficients(DOUBLE amin,DOUBLE amax,DOUBLE errtol,INT no,DOUBLE *bz,DOUBLE *dz){
        
        DOUBLE xn,zn,tempbeta,betak,yn,anorm,fmin,fmax,phi;
        
        // WARNING -- WHY IS THIS SO BIG AND DOESNT CHANGE E.G. WITH LATTICE SIZE //
        INT maxzolot=4000; /*size of phi during iteration*/
        
        DOUBLE *alpha=new DOUBLE[maxzolot];
        DOUBLE *gamma=new DOUBLE[maxzolot];
        
        INT npoints;
        DOUBLE K,pow2;
        DOUBLE dx;
        DOUBLE *c= new DOUBLE[2*no];
        DOUBLE eps;
        DOUBLE kappa;
        
        kappa=sqrt(1.0-amin/amax);
        K=CompElliptic(kappa,errtol);
        
        zn=0.5*K/no;
        
        for(INT l=0;l<2*no-1;l++){
            
            xn=(l+1)*zn;
            alpha[0]=1.0;
            betak=sqrtf(1.0-kappa);
            gamma[0]=sqrt(kappa);
            INT i=1;
            pow2=1.0;
            
            while(fabs(gamma[i-1]/alpha[i-1])>=errtol){
                
                tempbeta=betak;
                alpha[i]=(alpha[i-1]+tempbeta)*0.5;
                betak=sqrtf(alpha[i-1]*tempbeta);
                gamma[i]=(alpha[i-1]-tempbeta)*0.5;
                i++;
                pow2=pow2*2.0;
                
            }
            i=i-1;
            phi=pow2*alpha[i]*xn;
            while(i>=1){
                
                phi=0.5*(phi+std::asin(gamma[i]*std::sin(phi)/alpha[i]));
                i=i-1;
            }
            //sn=std::sin(phi);
            c[l]=std::tan(phi)*std::tan(phi);
        }
        
        /*For very small condition number*/
        for(INT l=0;l<no;l++){
            
            dz[l]=c[2*l]*amin;
            bz[l]=coeffb(c,l,no);
        }
        
        npoints=5000;
        dx=(std::log(amax)-std::log(amin))/npoints;
        
        yn=std::log(amin)-dx;
        
        // STARTING VALEUS //
        fmin=std::pow(1.0,10.0);
        fmax=0.0;
        for(INT i=0;i<npoints;i++){
            
            yn=yn+dx;
            xn=exp(yn);
            zn=0.0;
            
            for(INT l=0;l<no;l++){
                zn=zn+bz[l]/(xn+dz[l]);
            }
            
            zn=zn*sqrtf(xn);
            if(zn>fmax){
                fmax=zn;
            }
            if(zn<fmin){
                fmin=zn;
            }
        }
        /*    printf("%le %le\n",fmax,fmin);*/
        anorm=2.0/(fmin+fmax);
        eps=(fmax-fmin)/(fmax+fmin);
        /*      printf("the value of eps is %le\n",eps);*/
        
        for(INT k=0;k<no;k++){
            bz[k]=anorm*bz[k];
        }
        /*
         for(l=0;l<no;l++){
         std::cerr << "#ZF l b[l] d[l] " << l << " "  << bz[l] << " " << dz[l] << std::endl;
         }
         */
        delete c;
        delete alpha;
        delete gamma;
        
        return eps;
        
    }
}



