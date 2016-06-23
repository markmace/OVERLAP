#ifndef __LOOP_MACROS_CPP__
#define __LOOP_MACROS_CPP__

#define START_FERMION_LOOP(gMatrixCoupling) \
\
for(INT i=0;i<Nc;i++){ \
for(INT j=0;j<Nc;j++){ \
\
for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){\
for(INT beta=0;beta<DiracAlgebra::SpinorComponents;beta++){ \
\
if(gMatrixCoupling[alpha][beta]!=0.0){ \



#define END_FERMION_LOOP \
}\
}\
}\
}\
}

#define START_FERMION_LOOP_2(gMatrixCoupling,gMatrixCoupling2) \
\
for(INT i=0;i<Nc;i++){ \
for(INT j=0;j<Nc;j++){ \
\
for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){\
for(INT beta=0;beta<DiracAlgebra::SpinorComponents;beta++){ \
\
if(gMatrixCoupling[alpha][beta]!=0.0 || gMatrixCoupling2[alpha][beta]!=0.0){ \



#define END_FERMION_LOOP \
}\
}\
}\
}\
}


#define START_LOCAL_FERMION_LOOP(gMatrixCoupling) \
\
for(INT i=0;i<Nc;i++){ \
\
for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){\
for(INT beta=0;beta<DiracAlgebra::SpinorComponents;beta++){ \
\
if(gMatrixCoupling[alpha][beta]!=0.0){ \



#define END_LOCAL_FERMION_LOOP \
}\
}\
}\
}

#define START_LOCAL_FERMION_LOOP_2(gMatrixCoupling,gMatrixCoupling2) \
\
for(INT i=0;i<Nc;i++){ \
\
for(INT alpha=0;alpha<DiracAlgebra::SpinorComponents;alpha++){\
for(INT beta=0;beta<DiracAlgebra::SpinorComponents;beta++){ \
\
if(gMatrixCoupling[alpha][beta]!=0.0 || gMatrixCoupling2[alpha][beta]!=0.0){ \



#define END_LOCAL_FERMION_LOOP \
}\
}\
}\
}

#endif