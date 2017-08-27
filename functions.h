#ifndef FUNCTIONS
#define FUNCTIONS
///// define LAPACK fucntion /////
#define F77NAME(x) x##_
extern "C" {
    void F77NAME(dgbsv)(const int& n, const int& kl, const int& ku, 
                    const int& nrhs, const double * A, const int& ldab, 
                    int * ipiv, double * B, const int& ldb, int& info);
    void F77NAME(dscal) (const int& n, const double& alpha, double * x, const int& incx);

    void F77NAME(daxpy) (const int& n, const double& alpha, const double* x,
    				const int& incx, double * y, const int& incy);

    void F77NAME(dcopy) (const int& n, const double* x, const int& incx, 
    				double* y, const int& incy);
    void F77NAME(dgbmv)(const char& trans, const int& m, const int& n,
                const int& kl, const int& ku,
                const double& alpha, const double* a, const int& lda,
                const double* x, const int& incx, const double& beta,
                double* y, const int& incy);

}


// function declarations used in question 1
double *createElStiff(double A, double E, double I, double l);
double *assembleGlobStiffB(double *arr, int rs, int Nx);
double *createElForce(double qx, double qy, double l);
double *assembleGlobForce(double *arr, int rs, int Nx, double fy);
double *solveStatic(double A, double E, double I, double l, int rs, int Nx, double qx, double qy, double fy);

// additional function declarations for question 2
double *assembleGlobStiffBw(double *arr, int rs, int Nx);
double *assembleMass(double alpha, double l, double rho, double A, int rs, int Nx);
double *assembleMassB(double alpha, double l, double rho, double A, int rs, int Nx);
double *additionOp(double *A, double *B, int rs, double mult);
double *buildUn(int rs, double *stiffMat, double *massMat, double *un, double delT);
double *buildUnm(int rs, double *massMat, double *unminus1);
double *buildFn(int rs, double qx, double qy, double fy, double l, int Nx, double delT2);
double *buildB(int rs, double *fn, double *un_term, double *unminus1_term);
double *solveDynExp(double A, double E, double I, double l, int rs, int Nx, double qx, double qy, 
							double fy, double T, int Nt, double alpha, double delT, double rho, double Tt);

// // additional function declarations for question 3
double *buildKeff(int rs, double beta, double delT, double *massMat, double *stiffMat);
double *buildS(int rs, double beta, double delT, double *un, double *unDot, double *unDotDot);
double *buildMs(int rs, double *massMat, double *S);
double *buildFnplus1(int rs, double qx, double qy, double fy, double l, int Nx, double delT2);
double *buildBQ3(int rs, double *Fnplus1, double *MS);
double *buildUnDotDotp1(int rs, double beta, double delT, double *unp1, double *un, double *unDot, double *unDotDot);
double *buildUnDotp1(int rs, double gamma, double delT, double *unDot, double *unDotDot, double *unDotDotp1);
double *additionOpQ3(double *A, double *B, int rs, double mult);
double *solveDynImp(double A, double E, double I, double l, int rs, int Nx, double qx, double qy, 
							double fy, double T, int Nt, double alpha, double delT, double rho, double Tt);

//additonal function declarations for question 4
double *assembleGlobStiffBP(double *arr, int nLoc, int rs, int rank);
double *split(double *arr, int rs, int rs_p, int rank);
double *combine(double *arr1, double *arr2, int rs_p, int rs);
double *solveDynExpPar(double A, double E, double I, double l, int rs, int Nx, double qx, double qy, 
                            double fy, double T, int Nt, double alpha, double delT, double rho, int rank, int size);

// additional function declarations for question 5
double *addColumnK(int rs, double *arr);
double *addColumnB(int rs, double *arr);
double *assembleGlobStiffBP(double *arr, int rs, int Nx, int p, int rank);
double *buildKeffP(int rs, double beta, double delT, double *massMat, double *stiffMat, int z, int rank);
double *additionOpP(double *A, double *B, int rs, double mult, int z);
double *addColumnKP(int rs, double *arr, int z, int rank);
double *combineP(double *arr1, double*arr2, int rs);
double *addZeros(double *arr, int rs, int Nx, int rank);
double *buildS0(int rs, double beta, double delT, double *un, double *unDot, double *unDotDot, int rank);
double *buildMs0(int rs, double *massMat, double *S, int rank);
double *buildFnplus10(int rs, double qx, double qy, double fy, double l, int Nx, double delT2, int rank);
double *buildB0(int rs, double *Fnplus1, double *MS, int rank);
double *addColumnB0(int rs, double *arr, int rank);
double *buildUnDotDotp10(int rs, double beta, double delT, double *unp1, double *un, double *unDot,
                     double *unDotDot, int rank);
double *buildUnDotp10(int rs, double gamma, double delT, double *unDot, double *unDotDot, double *unDotDotp1, int rank);
double *solveDynImpPar(double A, double E, double I, double l, int rs, int Nx, double qx, double qy, 
                            double fy, double T, int Nt, double alpha, double delT, double rho, int rank, int size);
#endif