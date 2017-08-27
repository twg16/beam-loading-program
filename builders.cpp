//////////////////////////////////////////////
/////////////////////////////////////////////

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <mpi.h>

///// define LAPACK fucntions /////
#define F77NAME(x) x##_
extern "C" {
	void F77NAME(dgbsv)(const int& n, const int& kl, const int& ku, 
                    const int& nrhs, const double * A, const int& ldab, 
                    int * ipiv, double * B, const int& ldb, int& info);
    void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                    const int& lda, int * ipiv, double * B, const int& ldb,
                    int& info);

    void F77NAME(dgemv) (const char& trans,  const int& m,
    const int& n,       const double& alpha,
    const double* a,    const int& lda,
    const double* x,    const int& incx,
    const double& beta, double* y, const int& incy);

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
    void F77NAME(pdgbsv)(const int& n, const int& kl, const int& ku, 
                const int& nrhs, const double * A, const int& ja,
                const int* desca, int * ipiv, double * B, const int& ib,
                const int* descb, double* work, const int& lwork, 
                int* info);
    void F77NAME(pdgbtrf) (const int& n, const int& bwl, const int& bwu,
                const double* a,     const int& ja,
                const int* desca,    int* ipiv,
                double* AF, const int& LAF, double* work, const int& LWORK,
                int& info);
    void F77NAME(pdgbtrs) (const char& trans, const int& n, const int& kl,
                const int &ku, const int& nrhs,   const double* a,
                const int& ja, const int* desca, const int* ipiv, double* b,
                const int& ib, const int* descb, double* AF, const int& LAF,
                double* work, const int& LWORK, int& info);

    void Cblacs_get(int, int, int*);
	void Cblacs_pinfo(int*, int*);
	void Cblacs_gridinit(int*, char*, int, int);
	void Cblacs_gridinfo(int, int*, int*, int*, int*);
	void Cblacs_exit(int);
	void Cblacs_gridexit(int);
}


///// function to create elemental stiffness matrix /////
double *createElStiff(double A, double E, double I, double l)
{
	double *arr = 0;
	arr = new double[36];
	arr[0] = (A*E)/l;
	arr[3] = -(A*E)/l;
	arr[7] = (12*E*I)/(l*l*l);
	arr[8] = (6*E*I)/(l*l);
	arr[10] = -(12*E*I)/(l*l*l);
	arr[11] = (6*E*I)/(l*l);
	arr[13] = (6*E*I)/(l*l);
	arr[14] = (4*E*I)/(l);
	arr[16] = -(6*E*I)/(l*l);
	arr[17] = (2*E*I)/l;
	arr[18] = -(A*E)/l;
	arr[21] = (A*E)/l;
	arr[25] = -(12*E*I)/(l*l*l);
	arr[26] = -(6*E*I)/(l*l);
	arr[28] = (12*E*I)/(l*l*l);
	arr[29] = -(6*E*I)/(l*l);
	arr[31] = (6*E*I)/(l*l);
	arr[32] = (2*E*I)/l;
	arr[34] = -(6*E*I)/(l*l);
	arr[35] = (4*E*I)/(l);

	//std::cout << "------- Finished creating elemental stiffness matrix -------" << std::endl;
	return arr;
}

double *createElForce(double qx, double qy, double l)
{
	// create array that  will be returned
	double *arr = 0;
	arr = new double[6];
	for (int i = 0; i < 6; ++i)
	{
		arr[i] = 0.0;
	}

	arr[0] = l*qx/2;
	arr[1] = l*qy/2;
	arr[2] = l*(qy*l)/12;
	arr[3] = l*qx/2;
	arr[4] = l*qy/2;
	arr[5] = -l*(qy*l)/12;

	//std::cout << "------- Finished creating elemental force matrix -------" << std::endl;
	return arr;
}

double *assembleGlobForce(double *arr, int rs, int Nx, double fy)
{
	// creating array that will be returned
	double *globArr = 0;
	globArr = new double[rs];
	for (int i = 0; i < rs; ++i)
	{
		globArr[i] = 0.0;
	}

	// create index array
	int *indexArr = 0;
	indexArr = new int[rs];
	for (int i = 0; i <Nx-3; ++i)
	{
		indexArr[i] = i*3;
	}

	for (int k = 0; k < rs-3; k+=3)
	{
		for (int i = k; i < k+6; ++i)
		{
			globArr[i] += arr[i-k];
		}

		if (k == 0)
		{
			for(int j = 0; j < 3; ++j)
			{
				globArr[j] += arr[j+3];
			}
		}
		if (k == rs-6)
		{
			for(int j = k+3; j < k+6; ++j)
			{
				globArr[j] += arr[(j%6)];
			}
		}
	}
	//apply the concenetrated force
	globArr[rs/2] = globArr[rs/2] + fy;
	//std::cout << "------- Finished assembling global matrix -------" << std::endl;
	return globArr;
}

// CONSTRUCT BANDED FORM
double *assembleGlobStiffB(double *arr, int rs, int Nx)
{
	double *globArr = 0;
	globArr = new double[rs*13];
	// making arrays to refer to
	double a[3] = {0, arr[11], 0};
	double b[3] = {arr[3], arr[10], arr[17]};
	double c[3] = {0, 0, arr[16]};
	double d[3] = {0, arr[8]+arr[29], 0};
	double e[3] = {arr[0]+arr[21], arr[7]+arr[28], arr[14]+arr[35]};
	double f[3] = {0, arr[13]+arr[34], 0};
	double g[3] = {0, 0, arr[26]};
	double h[3] = {arr[18], arr[25], arr[32]};
	double i[3] = {0, arr[31], 0};


	for (int n = 0; n < rs*13; ++n)
	{
		if (n%13 == 4 and n > 13*4)
		{
			globArr[n] = a[((n+1)%3)];
			// std::cout << ((n+1)%3) << std::endl;
		}

		if (n%13 == 5 and n > 13*3)
		{
			//std::cout << ((n+1)%3) << std::endl;
			globArr[n] = b[((n+1)%3)];
		}

		if (n%13 == 6 and n > 13*2)
		{
			//std::cout << ((n+1)%3) << std::endl;
			globArr[n] = c[((n+1)%3)];
		}

		if (n%13 == 7 and n > 13*1)
		{
			// std::cout << ((n+1)%3) << std::endl;
			globArr[n] = d[((n+1)%3)];
		}

		if (n%13 == 8 and n > 13*0)
		{
			// std::cout << ((n+1)%3) << std::endl;
			globArr[n] = e[((n+1)%3)];
		}

		if (n%13 == 9 and n < 13*(rs-1))
		{
			// std::cout << ((n)%3) << std::endl;
			globArr[n] = f[((n)%3)];
		}

		if (n%13 == 10 and n < 13*(rs-2))
		{
			// std::cout << ((n-1)%3) << std::endl;
			globArr[n] = g[((n-1)%3)];
		}

		if (n%13 == 11 and n < 13*(rs-3))
		{
			// std::cout << ((n-2)%3) << std::endl;
			globArr[n] = h[((n-2)%3)];
		}

		if (n%13 == 12 and n < 13*(rs-4))
		{
			// std::cout << ((n)%3) << std::endl;
			globArr[n] = i[((n)%3)];
		}
	}
	return globArr;

}
// ASSEMBLE BANDED STRUCTURE WITHOUT EXTRA ZEROS USING ARRAY WITH EXTRA ZEROS
double *assembleGlobStiffBw(double *arr, int rs, int Nx) //without extra lines above
{
	double *globArr = 0;
	globArr = new double[rs*9];
	double *indexArr = new double[rs];
	int q = 0;
	// make array with indices pointing to columns
	for (int k = 4; k < rs*13; k=k+13)
	{
		indexArr[q] = k;
		q = q + 1;
	}

	for (int w = 0; w < rs; ++w)
	{
		for (int y = 0; y < 9; ++y)
		{
			int z = indexArr[w]+y;
			globArr[(w*9)+y] = arr[z];
		}
	}

	delete[] indexArr;
	return globArr;
}
// QUESTION 2 FUNCTIONS //
double *assembleMass(double alpha, double l, double rho, double A, int rs, int Nx)
{
	//double alpha = 1.0/24.0;
	double *arr = 0;
	arr = new double[rs];
	double *ellArr = 0;
	ellArr = new double[6];
	double *massArr = 0;
	massArr = new double[rs*rs];

	// create index array
	int *indexArr = 0;
	indexArr = new int[rs];
	for (int i = 0; i <Nx-3; ++i)
	{
		indexArr[i] = i*3;
	}
	// create elemental matrix as array with diagonals
	ellArr[0] = (rho*A*l)*0.5;
	ellArr[1] = (rho*A*l)*0.5;
	ellArr[2] = (rho*A*l)*(alpha*l*l);
	ellArr[3] = (rho*A*l)*0.5;
	ellArr[4] = (rho*A*l)*0.5;
	ellArr[5] = (rho*A*l)*(alpha*l*l);
	for (int k = 0; k < rs-3; k+=3)
	{
		for (int i = k; i < k+6; ++i)
		{
			arr[i] += ellArr[i-k];
		}

		if (k == 0)
		{
			for(int j = 0; j < 3; ++j)
			{
				arr[j] += ellArr[j+3];
			}
		}
		if (k == rs-6)
		{
			for(int j = k+3; j < k+6; ++j)
			{
				arr[j] += ellArr[(j%6)];
			}
		}
	}
	// assemble the global mass matrix from elemental array
	int j = 0;
	for (int i = 0; i < (rs*rs); i +=rs+1)
	{
		massArr[i] = 2*ellArr[j%6];
		j++;
	}
	delete[] indexArr;
	delete[] ellArr;
	delete[] arr;
	//std::cout << "------- Finished assembling mass matrix -------" << std::endl;
	return massArr;
}

double *assembleMassB(double alpha, double l, double rho, double A, int rs, int Nx)
{
	// double *arr = 0;
	// arr = new double[rs];
	double *ellArr = 0;
	ellArr = new double[6];
	double *massArr = 0;
	massArr = new double[rs]; // had to change from rs*rs to rs. NOTE THIS

	// create index array
	int *indexArr = 0;
	indexArr = new int[rs];
	for (int i = 0; i <Nx-3; ++i)
	{
		indexArr[i] = i*3;
	}
	// create elemental matrix as array with diagonals
	ellArr[0] = 2*(rho*A*l)*0.5;
	ellArr[1] = 2*(rho*A*l)*0.5;
	ellArr[2] = 2*(rho*A*l)*(alpha*l*l);
	ellArr[3] = 2*(rho*A*l)*0.5;
	ellArr[4] = 2*(rho*A*l)*0.5;
	ellArr[5] = 2*(rho*A*l)*(alpha*l*l);
	
	for (int i = 0; i < rs; ++i)
	{
		massArr[i] = ellArr[i%6];
		// std::cout << massArr[i] << std::endl;
	}
	return massArr;
}

double *additionOp(double *A, double *B, int rs, double mult)
{
	// make a copy of A
	double *foo = new double[rs*9];
	for (int i = 0; i < rs*9; ++i)
	{
		foo[i] = A[i];
	}

	int q = 0;
	for (int i = 0; i < rs*9; ++i)
	{
		if (i%9 == 4)
		{
			foo[i] = foo[i] + (mult*B[q]);
			q = q + 1;
		}
	}

	return foo;
}


double *buildUn(int rs, double *stiffMat, double *massMat, double *un, double delT)
{
	double delT2 = delT * delT;

	double *var = 0;
	var = new double[rs];
	double *var1 = 0;
	var1 = new double[rs*9];
	double *var2 = new double[rs*9];
	double *un_term = 0;
	un_term = new double[rs];
	// double var2 = new double[rs]
	for (int i = 0; i < rs; ++i)
	{
		var[i] = 0.0;
		un_term[i] = 0.0;
	}
	for (int i = 0; i < rs*9; ++i)
	{
		var1[i] = 0.0;
		var2[i] = 0.0;
	}
	//F77NAME(dcopy) (rs*rs, massMat, 1, var, 1);
	F77NAME(dcopy) (rs*9, stiffMat, 1, var1, 1);
	//multiply stiffness by delT squared = var1
	F77NAME(dscal) (rs*9, delT2, var1, 1);
	// return var1;
	// calculating -2 * M = var
	F77NAME(daxpy) (rs, -2, massMat, 1, var, 1);
	// return var1;
	// calculating Un coefficient = var
	// F77NAME(daxpy) (rs*9, 1.0, var1, 1, var, 1);
	var2 = additionOp(var1, var, rs, 1.0);
	// // //calculating un term = un
	F77NAME(dgbmv) ('N', rs, rs, 4, 4, 1.0, var2, 9, un, 1, 0.0, un_term, 1);
	delete[] var1;
	delete[] var;
	delete[] var2;
	return un_term;
}

double *buildUnm(int rs, double *massMat, double *unminus1)
{
	double *unminus1_term = 0;
	unminus1_term = new double[rs];
	for (int i = 0; i < rs; ++i)
	{
		unminus1_term[i] = 0;
	}
	for (int i = 0; i < rs; ++i)
	{
		unminus1_term[i] = massMat[i] * unminus1[i];
	}
	// return unminus1_term;
	return unminus1_term;
}

double *buildFn(int rs, double qx, double qy, double fy, double l, int Nx, double delT2)
{
	double *arr = 0;
	arr = new double[rs];
	for (int i = 0; i < rs; ++i)
	{
		arr[i] = 0;
	}
	double *elForceMat = createElForce(qx, qy, l);
	double *globForceMat = assembleGlobForce(elForceMat, rs, Nx, fy);
	F77NAME(dscal) (rs, delT2, globForceMat, 1);
	arr = globForceMat;
	return arr;
}
double *buildB(int rs, double *fn, double *un_term, double *unminus1_term)
{
	double *arr = new double[rs];
	for (int i = 0; i < rs; ++i)
	{
		arr[i] = fn[i];
	}
	F77NAME(daxpy) (rs, -1, un_term, 1, arr, 1);
	F77NAME(daxpy) (rs, -1, unminus1_term, 1, arr, 1);
	//std::cout << "build 4 complete" << std::endl;
	return arr;
}

void solveDynamicExp(double alpha, double l, double A, double rho, int rs, int Nx, double qy, int Nt,
					double fy, double T, double delT, double E, double I, double qx)
{
	// create mass matrix
	double *massMat =  assembleMass(alpha, l, A, rho, rs, Nx);
	// instantiate necessary arrays to use throughout
	double *un = new double[rs];
	double *unminus1 = new double[rs];
	double *un_term = new double[rs];
	double *unminus1_term = new double[rs];
	for (int i = 0; i < rs; ++i)
	{
		un[i] = 0.0;
		unminus1[i] = 0.0;
		un_term[i] = 0.0;
		unminus1_term[i] = 0.0;
	}
	// define step in forces applied
	double delqy = qy/(Nt-1);
	double delfy = fy/(Nt-1);
	// reset forces before time loop starts
	qy = 0;
	fy = 0;
	// TIME LOOP //
	for (double n = 0; n <T; n+= delT)
	{
		// build global stiffness matrix
		double *elStiffMat = createElStiff(A, E, I, l);
		double *globStiffMat = assembleGlobStiffB(elStiffMat, rs, Nx);
		// std::cout << "qy = " << qy << std::endl;
		// std::cout << "fy = " << fy << std::endl;
		// build right hand side terms and put together using buildB()
		double *un_term = buildUn(rs, globStiffMat, massMat, un, delT);
		double *unminus1_term = buildUnm(rs, massMat, unminus1);
		double *fn = buildFn(rs, qx, qy, fy, l, Nx, delT*delT);
		double *unplus1 = buildB(rs, fn, un_term, unminus1_term);
		//solving the Ax=b problem
		const int nrhs = 1;
    	int info = 0;
    	int* ipiv = new int[rs];
	    F77NAME(dgesv) (rs, nrhs, massMat, rs, ipiv, unplus1, rs, info);
	    double *soln = new double [rs];
	    soln = unplus1;
	 //    for (int i = 0; i < rs; ++i)
		// {
		// 	std::cout << soln[i] << std::endl;
		// }
	    // using dcopy to set Un = Un=1 and Un+1 = Un
		F77NAME(dcopy) (rs, un, 1.0, unminus1, 1);
		F77NAME(dcopy) (rs, soln, 1.0, un, 1);
		//update forces applied to the beam
   		qy = qy + delqy;
		fy = fy + delfy;
		// solfile << un[(rs/2)];
		// solfile << "\n";
	}
	// solution written to file
	std::ofstream solfile;
	solfile.open("solution.txt");
	for (int i = 0; i < rs; ++i)
	{
		solfile << un[i];
		solfile << "\n";
	}
	solfile.close();
}

/////// Q3 functions ///////

double *additionOpQ3(double *A, double *B, int rs, double mult)
{
	// make a copy of A
	double *foo = new double[rs*13];
	for (int i = 0; i < rs*13; ++i)
	{
		foo[i] = A[i];
	}
	int q = 0;
	for (int i = 0; i < rs*13; ++i)
	{
		if (i%13 == 8)
		{
			foo[i] = foo[i] + (mult*B[q]);
			q = q + 1;
		}
	}
	return foo;
}

double *buildKeff(int rs, double beta, double delT, double *massMat, double *stiffMat)
{
	double *arr = new double[rs*13];
	double *arr1 = new double[rs*13];
	F77NAME(dcopy) (rs*13, stiffMat, 1.0, arr, 1);
	double mult = 1/(beta*delT*delT);
	arr1 = additionOpQ3(arr, massMat, rs, mult);
	delete[] arr;
	return arr1;
}

double *buildS(int rs, double beta, double delT, double *un, double *unDot, double *unDotDot)
{
	double mult1 = 1/(beta*delT*delT);
	double mult2 = 1/(beta*delT);
	double mult3 = (1/(2*beta))-1;

	double *var1 = new double[rs];
	F77NAME(dcopy) (rs, un, 1.0, var1, 1); // var1 = un
	double *var2 = new double[rs];
	F77NAME(dcopy) (rs, unDot, 1.0, var2, 1); // var2 = unDot
	double *var3 = new double[rs];
	F77NAME(dcopy) (rs, unDotDot, 1.0, var3, 1); // var3 = unDotDot

	// multiply un by mult1
	F77NAME(dscal) (rs, mult1, var1, 1);
	// multiply unDot by mult2
	F77NAME(dscal) (rs, mult2, var2, 1);
	// sum (mult1*un) and (mult2*unDot)
	F77NAME(daxpy) (rs, 1.0, var1, 1, var2, 1); // result = var2
	// add (mult3*unDotDot) and answer to line above
	F77NAME(daxpy) (rs, mult3, var3, 1, var2, 1); // result = var2
	return var2;
}

double *buildMs(int rs, double *massMat, double *S)
{
	// initialise array to return
	double *arr = new double[rs];
	for (int i = 0; i < rs; ++i)
	{
		arr[i] = 0;;
	}
	for (int i = 0; i < rs; ++i)
	{
		arr[i] = massMat[i] * S[i];
	}
	// mnultiply M by S
	//F77NAME(dgemv) ('N', rs, rs, 1.0, massMat, rs, S, 1, 0.0, arr, 1);
	return arr;
}

double *buildFnplus1(int rs, double qx, double qy, double fy, double l, int Nx, double delT2)
{
	// double *arr = 0;
	// arr = new double[rs];
	// for (int i = 0; i < rs; ++i)
	// {
	// 	arr[i] = 0;
	// }
	double *elForceMat = createElForce(qx, qy, l);
	double *globForceMat = assembleGlobForce(elForceMat, rs, Nx, fy);
	//std::cout << "GLOBFORCEMAT[1] = " << globForceMat[1] << std::endl;
	//std::cout << "DELT2 = " << delT2 << std::endl; 
	//F77NAME(dscal) (rs, delT2, globForceMat, 1);
	//arr = globForceMat;
	return globForceMat;
}

double *buildBQ3(int rs, double *Fnplus1, double *MS)
{
	// std::cout << "building B Q3 way" << std::endl; 
	double *var4 = new double[rs];
	F77NAME(dcopy) (rs, MS, 1.0, var4, 1); // var4 = MS

	F77NAME(daxpy) (rs, 1.0, Fnplus1, 1, var4, 1);
	return var4;
}
double *buildUnDotDotp1(int rs, double beta, double delT, double *unp1, double *un, double *unDot, double *unDotDot)
{
	double *var5 = new double[rs];
	F77NAME(dcopy) (rs, unp1, 1.0, var5, 1); // var5 = unp1
	double *var6 = new double[rs];
	F77NAME(dcopy) (rs, unDot, 1.0, var6, 1); // var6 = unDot
	double *var7 = new double[rs];
	F77NAME(dcopy) (rs, un, 1.0, var7, 1); // var7 = un
	double *var8 = new double[rs];
	F77NAME(dcopy) (rs, unDotDot, 1.0, var8, 1); // var8 = unDotDot

	double mult1 = 1/(beta*delT*delT);
	double mult2 = -1*(1/(beta*delT));
	double mult3 = -1*((1/(2*beta))-1);
	// calculate first term
	F77NAME(daxpy) (rs, -1.0, var7, 1, var5, 1); // result = var5
	F77NAME(dscal) (rs, mult1, var5, 1); // result = var5
	// calculate second term
	F77NAME(dscal) (rs, mult2, var6, 1); // result = var6
	// calculate third term
	F77NAME(dscal) (rs, mult3, var8, 1); // result = var8
	// sum terms held in un, unDot and unDotDot
	F77NAME(daxpy) (rs, 1.0, var5, 1, var6, 1); // result = var6
	F77NAME(daxpy) (rs, 1.0, var6, 1, var8, 1); // result = var8
	return var8;
}
double *buildUnDotp1(int rs, double gamma, double delT, double *unDot, double *unDotDot, double *unDotDotp1)
{
	double *var9 = new double[rs];
	F77NAME(dcopy) (rs, unDot, 1.0, var9, 1); // var9 = unDot

	double mult1 = delT*(1-gamma);
	double mult2 = delT*gamma;
	// sum first and second first terms
	F77NAME(daxpy) (rs, mult1, unDotDot, 1, var9, 1); // result = var9
	// add final term to first and second terms
	F77NAME(daxpy) (rs, mult2, unDotDotp1, 1, var9, 1); // result = var9
	return var9;
}  

double *assembleGlobStiffBP(double *arr, int nLoc, int rs, int rank)
{
	// need to change rs to n-m
	double *globArr = 0;
	globArr = new double[rs*13];
	// making arrays to refer to
	double a[3] = {0, arr[11], 0};
	double b[3] = {arr[3], arr[10], arr[17]};
	double c[3] = {0, 0, arr[16]};
	double d[3] = {0, arr[8]+arr[29], 0};
	double e[3] = {arr[0]+arr[21], arr[7]+arr[28], arr[14]+arr[35]};
	double f[3] = {0, arr[13]+arr[34], 0};
	double g[3] = {0, 0, arr[26]};
	double h[3] = {arr[18], arr[25], arr[32]};
	double i[3] = {0, arr[31], 0};

	// for (int n = start; n <finish*13; ++n)
	for (int n = 0; n < 13*rs; ++n)
	{
		if (n%13 == 4 and n > 13*4)
		{
			globArr[n] = a[((n+1)%3)];
			// std::cout << ((n+1)%3) << std::endl;
		}

		if (n%13 == 5 and n > 13*3)
		{
			//std::cout << ((n+1)%3) << std::endl;
			globArr[n] = b[((n+1)%3)];
		}

		if (n%13 == 6 and n > 13*2)
		{
			//std::cout << ((n+1)%3) << std::endl;
			globArr[n] = c[((n+1)%3)];
		}

		if (n%13 == 7 and n > 13*1)
		{
			// std::cout << ((n+1)%3) << std::endl;
			globArr[n] = d[((n+1)%3)];
		}

		if (n%13 == 8 and n > 13*0)
		{
			// std::cout << ((n+1)%3) << std::endl;
			globArr[n] = e[((n+1)%3)];
		}

		if (n%13 == 9 and n < 13*(rs-1))
		{
			// std::cout << ((n)%3) << std::endl;
			globArr[n] = f[((n)%3)];
		}

		if (n%13 == 10 and n < 13*(rs-2))
		{
			// std::cout << ((n-1)%3) << std::endl;
			globArr[n] = g[((n-1)%3)];
		}

		if (n%13 == 11 and n < 13*(rs-3))
		{
			// std::cout << ((n-2)%3) << std::endl;
			globArr[n] = h[((n-2)%3)];
		}

		if (n%13 == 12 and n < 13*(rs-4))
		{
			// std::cout << ((n)%3) << std::endl;
			globArr[n] = i[((n)%3)];
		}
	}

	return globArr;
}

double *split(double *arr, int rs, int rs_p, int rank)
{
	double *fArr = new double[rs_p];

	if (rank == 0)
	{
		for (int i = 0; i < rs_p; ++i)
		{
			fArr[i] = arr[i];
		}
	}
	if (rank == 1)
	{
		for (int i = 0; i < rs_p; ++i)
		{
			fArr[i] = arr[(rs-rs_p)+i];
		}
	}
	return fArr;
}
double *combine(double *arr1, double *arr2, int rs_p, int rs)
{
	double *arr = new double[rs];
	for (int i = 0; i < rs; ++i)
	{
		arr[i] = 0;
	}
	for (int i = 0; i < rs_p-3; ++i)
	{
		arr[i] = arr1[i];
	}
	int cnt = 5;
	for (int i = rs_p-4; i < rs; ++i)
	{
		arr[i] = arr2[cnt];
		cnt = cnt + 1; 
	}

	return arr;
}
//

double *addColumnK(int rs, double *arr)
{
	double *returnArr = new double[(rs+1)*13];
	for (int i = 0; i < rs*13; ++i)
	{
		returnArr[i] = arr[i];
	} 

	returnArr[((rs+1)*13)-5] = 1.0;
	return returnArr;
}
double *addColumnB(int rs, double *arr)
{
	double *returnArr = new double[(rs+1)];
	for (int i = 0; i < rs; ++i)
	{
		returnArr[i] = arr[i];
	} 
	return arr;
}

double *assembleGlobStiffBP(double *arr, int rs, int Nx, int p, int rank)
{
	if (rank ==0)
	{
		int k = 9+p;
		double *globArr = 0;
		globArr = new double[rs*k];
		// making arrays to refer to
		double a[3] = {0, arr[11], 0};
		double b[3] = {arr[3], arr[10], arr[17]};
		double c[3] = {0, 0, arr[16]};
		double d[3] = {0, arr[8]+arr[29], 0};
		double e[3] = {arr[0]+arr[21], arr[7]+arr[28], arr[14]+arr[35]};
		double f[3] = {0, arr[13]+arr[34], 0};
		double g[3] = {0, 0, arr[26]};
		double h[3] = {arr[18], arr[25], arr[32]};
		double i[3] = {0, arr[31], 0};

		for (int n = 0; n < rs*k; ++n)
		{
			if (n%k == p and n > k*4)
			{
				globArr[n] = a[((n+1)%3)];
				// std::cout << ((n+1)%3) << std::endl;
			}

			if (n%k == p+1 and n > k*3)
			{
				//std::cout << ((n+1)%3) << std::endl;
				globArr[n] = b[((n+1)%3)];
			}

			if (n%k == p+2 and n > k*2)
			{
				//std::cout << ((n+1)%3) << std::endl;
				globArr[n] = c[((n+1)%3)];
			}

			if (n%k == p+3 and n > k*1)
			{
				// std::cout << ((n+1)%3) << std::endl;
				globArr[n] = d[((n+1)%3)];
			}

			if (n%k == p+4 and n > k*0)
			{
				// std::cout << ((n+1)%3) << std::endl;
				globArr[n] = e[((n+1)%3)];

			}

			if (n%k == p+5 and n < k*(rs-1))
			{
				// std::cout << ((n)%3) << std::endl;
				globArr[n] = f[((n)%3)];
			}

			if (n%k == p+6 and n < k*(rs-2))
			{
				// std::cout << ((n-1)%3) << std::endl;
				globArr[n] = g[((n-1)%3)];
			}

			if (n%k == p+7 and n < k*(rs-3))
			{
				// std::cout << ((n-2)%3) << std::endl;
				globArr[n] = h[((n-2)%3)];
			}

			if (n%k == p+8 and n < k*(rs-4))
			{
				// std::cout << ((n)%3) << std::endl;
				globArr[n] = i[((n)%3)];
			}
		}
		return globArr;
	}
	else
	{
		return 0;
	}
}

double *addZeros(double *arr, int rs, int Nx, int rank) //without extra lines above
{
	if (rank == 0)
	{
		double *globArr = 0;
		globArr = new double[rs*17];
		double *indexArr = new double[rs];
		int q = 0;
		// make array with indices pointing to columns
		for (int k = 4; k < rs*13; k=k+13)
		{
			indexArr[q] = k;
			q = q + 1;
		}

		for (int w = 0; w < rs; ++w)
		{
			for (int y = 0; y < 9; ++y)
			{
				int z = indexArr[w]+y;
				globArr[(w*17)+y+8] = arr[z];
			}
		}

		delete[] indexArr;
		return globArr;
	}
	else
	{
		delete[] arr;
		return 0;
	}
}

double *additionOpP(double *A, double *B, int rs, double mult, int z)
{
	// make a copy of A
	double *foo = new double[rs*z];
	for (int i = 0; i < rs*z; ++i)
	{
		foo[i] = A[i];
	}
	int q = 0;
	for (int i = 0; i < rs*z; ++i)
	{
		if (i%z == 12)
		{
			foo[i] = foo[i] + (mult*B[q]);
			q = q + 1;
		}
	}
	return foo;
}

double *buildKeffP(int rs, double beta, double delT, double *massMat, double *stiffMat, int z, int rank)
{
	if (rank == 0)
	{
		double *arr = new double[(rs+1)*z];
		double *arr1 = new double[(rs+1)*z];
		F77NAME(dcopy) ((rs+1)*z, stiffMat, 1.0, arr, 1);
		double mult = 1/(beta*delT*delT);
		arr1 = additionOpP(arr, massMat, rs, mult, 17);
		delete[] arr;
		return arr1;
	}
	else
	{
		return 0;
	}
}

double *addColumnKP(int rs, double *arr, int z, int rank)
{
	if (rank == 0)
	{
		double *returnArr = new double[(rs+1)*z];
		for (int i = 0; i < rs*z; ++i)
		{
			returnArr[i] = arr[i];
		} 

		returnArr[((rs+1)*z)-5] = 1.0;
		delete[] arr;
		return returnArr;
	}
	else 
	{
		return 0;
	}
}

// double *combineP(double *arr1, double*arr2, int rs)
// {
// 	double *arr = new double[rs];
// 	for (int i = 0; i < rs; ++i)
// 	{
// 		if(i < (rs+1)/2)
// 		{
// 			arr[i] = arr1[i];
// 		}
// 		if(i > (rs+1)/2)
// 		{
// 			arr[i] = arr2[i];
// 		}
// 	}
// }

double *buildS0(int rs, double beta, double delT, double *un, double *unDot, double *unDotDot, int rank)
{
	if (rank == 0)
	{
		double mult1 = 1/(beta*delT*delT);
		double mult2 = 1/(beta*delT);
		double mult3 = (1/(2*beta))-1;

		double *var1 = new double[rs];
		F77NAME(dcopy) (rs, un, 1.0, var1, 1); // var1 = un
		double *var2 = new double[rs];
		F77NAME(dcopy) (rs, unDot, 1.0, var2, 1); // var2 = unDot
		double *var3 = new double[rs];
		F77NAME(dcopy) (rs, unDotDot, 1.0, var3, 1); // var3 = unDotDot

		// multiply un by mult1
		F77NAME(dscal) (rs, mult1, var1, 1);
		// multiply unDot by mult2
		F77NAME(dscal) (rs, mult2, var2, 1);
		// sum (mult1*un) and (mult2*unDot)
		F77NAME(daxpy) (rs, 1.0, var1, 1, var2, 1); // result = var2
		// add (mult3*unDotDot) and answer to line above
		F77NAME(daxpy) (rs, mult3, var3, 1, var2, 1); // result = var2
		return var2;
	}
	else{
		return 0;
	}
}

double *buildMs0(int rs, double *massMat, double *S, int rank)
{
	if (rank == 0)
	{
		// initialise array to return
		double *arr = new double[rs];
		for (int i = 0; i < rs; ++i)
		{
			arr[i] = 0;;
		}
		for (int i = 0; i < rs; ++i)
		{
			arr[i] = massMat[i] * S[i];
		}
		// mnultiply M by S
		//F77NAME(dgemv) ('N', rs, rs, 1.0, massMat, rs, S, 1, 0.0, arr, 1);
		return arr;
	}
	else{ 
		return 0;}
}

double *buildFnplus10(int rs, double qx, double qy, double fy, double l, int Nx, double delT2, int rank)
{
	if (rank == 0)
	{
		double *elForceMat = createElForce(qx, qy, l);
		double *globForceMat = assembleGlobForce(elForceMat, rs, Nx, fy);
		return globForceMat;
	}
	else{
		return 0;}
}

double *buildB0(int rs, double *Fnplus1, double *MS, int rank)
{
	if (rank == 0)
	{
		double *var4 = new double[rs];
		F77NAME(dcopy) (rs, MS, 1.0, var4, 1); // var4 = MS

		F77NAME(daxpy) (rs, 1.0, Fnplus1, 1, var4, 1);
		return var4;
	}
	else{
		return 0;}
}

double *addColumnB0(int rs, double *arr, int rank)
{
	if (rank == 0)
	{
		double *returnArr = new double[(rs+1)];
		for (int i = 0; i < rs; ++i)
		{
			returnArr[i] = arr[i];
		} 
		return arr;
	}
	else{	
		return 0;}
}

double *buildUnDotDotp10(int rs, double beta, double delT, double *unp1, double *un, double *unDot, double *unDotDot, int rank)
{
	if (rank == 0)
	{

		double *var5 = new double[rs];
		F77NAME(dcopy) (rs, unp1, 1.0, var5, 1); // var5 = unp1
		double *var6 = new double[rs];
		F77NAME(dcopy) (rs, unDot, 1.0, var6, 1); // var6 = unDot
		double *var7 = new double[rs];
		F77NAME(dcopy) (rs, un, 1.0, var7, 1); // var7 = un
		double *var8 = new double[rs];
		F77NAME(dcopy) (rs, unDotDot, 1.0, var8, 1); // var8 = unDotDot

		double mult1 = 1/(beta*delT*delT);
		double mult2 = -1*(1/(beta*delT));
		double mult3 = -1*((1/(2*beta))-1);
		// calculate first term
		F77NAME(daxpy) (rs, -1.0, var7, 1, var5, 1); // result = var5
		F77NAME(dscal) (rs, mult1, var5, 1); // result = var5
		// calculate second term
		F77NAME(dscal) (rs, mult2, var6, 1); // result = var6
		// calculate third term
		F77NAME(dscal) (rs, mult3, var8, 1); // result = var8
		// sum terms held in un, unDot and unDotDot
		F77NAME(daxpy) (rs, 1.0, var5, 1, var6, 1); // result = var6
		F77NAME(daxpy) (rs, 1.0, var6, 1, var8, 1); // result = var8
		return var8;
	}
	else{return 0;}
}
double *buildUnDotp10(int rs, double gamma, double delT, double *unDot, double *unDotDot, double *unDotDotp1, int rank)
{
	if (rank == 0)
	{
		double *var9 = new double[rs];
		F77NAME(dcopy) (rs, unDot, 1.0, var9, 1); // var9 = unDot

		double mult1 = delT*(1-gamma);
		double mult2 = delT*gamma;
		// sum first and second first terms
		F77NAME(daxpy) (rs, mult1, unDotDot, 1, var9, 1); // result = var9
		// add final term to first and second terms
		F77NAME(daxpy) (rs, mult2, unDotDotp1, 1, var9, 1); // result = var9
		return var9;
	}
	else{
		return 0;
	}
}


//////////////////////////////////////////////////
/////////////// SOLVERS //////////////////////////
//////////////////////////////////////////////////
double *solveStatic(double A, double E, double I, double l, int rs, int Nx, double qx, double qy, double fy)
{
	double *elStiffMat = createElStiff(A, E, I, l);
	double *globStiffMat = assembleGlobStiffB(elStiffMat, rs, Nx);
	double *elForceMat = createElForce(qx, qy, l);
	double *globForceMat = assembleGlobForce(elForceMat, rs, Nx, fy);
	const int nrhs = 1;
    int info = 0;
    int* ipiv = new int[rs];
    F77NAME(dgbsv) (rs, 4, 4, nrhs, globStiffMat, 13, ipiv, globForceMat, rs, info);
    // std::cout << info << std::endl;
    return globForceMat;
}

double *solveDynExp(double A, double E, double I, double l, int rs, int Nx, double qx, double qy, double fy, double T, int Nt, double alpha, double delT, double rho, double Tt)
{

	// instantiate necessary arrays to use throughout
	double *un = new double[rs]();
	double *unminus1 = new double[rs]();
	double *un_term = new double[rs]();
	double *unminus1_term = new double[rs]();
	double *fn = new double[rs]();
	double *unplus1 = new double[rs]();
	double *soln = new double[rs]();

	//parameter for altering loading speed
	// double Tt = T/4;
	//////////////////////////////////////
	double NTt = (Tt/delT);
	// std::cout << "Tt/Nt = " << NTt << std::endl;
	// define step in forces applied
	double delqy = qy/(NTt);
	double delfy = fy/(NTt);
	// reset forces before time loop starts
	qy = 0;
	fy = 0;

	double *elStiffMat = createElStiff(A, E, I, l);
	double *globStiffMatB = assembleGlobStiffB(elStiffMat, rs, Nx);
	double *globStiffMat = assembleGlobStiffBw(globStiffMatB, rs, Nx);
	double *massMat =  assembleMassB(alpha, l, A, rho, rs, Nx);

	// TIME LOOP //
	//iterator for writing centre points
	int i = 0;
	double *cenDisp = new double[Nt];
	for (double n = 0; n < T; n+= delT)
	{		
		// build right hand side terms and put together using buildB()
		un_term = buildUn(rs, globStiffMat, massMat, un, delT);
		unminus1_term = buildUnm(rs, massMat, unminus1);
		fn = buildFn(rs, qx, qy, fy, l, Nx, delT*delT);
		unplus1 = buildB(rs, fn, un_term, unminus1_term);

		//solving the Ax=b problem
    	int info = 0;
    	int* ipiv = new int[rs];

	    // F77NAME(dgesv) (rs, nrhs, massMat, rs, ipiv, unplus1, rs, info);
	    F77NAME(dgbsv) (rs, 0, 0, 1, massMat, 1, ipiv, unplus1, rs, info);

	    soln = unplus1;

	    // using dcopy to set Un = Un=1 and Un+1 = Un
		F77NAME(dcopy) (rs, un, 1.0, unminus1, 1);
		F77NAME(dcopy) (rs, soln, 1.0, un, 1);
	
		// update forces applied to the beam
		if (n < Tt)
		{
   			qy = qy + delqy;
			fy = fy + delfy;
		}
		else{}

		cenDisp[i] = soln[int((rs)/2)];
		i += 1;

	}
	std::ofstream solfile;
	solfile.open("solution_2.txt");
	for (int i = 0; i < rs; ++i)
		{
			solfile << soln[i];
			solfile << "\n";
		}
	solfile.close();

	// returning array of central node displacement
	return cenDisp;
}

double *solveDynImp(double A, double E, double I, double l, int rs, int Nx, double qx, double qy, 
							double fy, double T, int Nt, double alpha, double delT, double rho, double Tt)
{
	double beta = 0.25;
	double gamma = 0.5;
	double NTt = (Tt/delT);
	// define step in forces applied
	double delqy = qy/(NTt);
	double delfy = fy/(NTt);
	// reset forces before time loop starts
	qy = 0;
	fy = 0;

	double *un = new double[rs]();
	double *unCopy = new double[rs]();
	double *unDot = new double[rs]();
	double *unDotDot = new double[rs]();
	double *unDotDotCopy = new double[rs]();
	double *unp1 = new double[rs]();
	double *unDotp1 = new double[rs]();
	double *unDotDotp1 = new double[rs]();
	
	double *massMat =  assembleMassB(alpha, l, A, rho, rs, Nx);
	double *elStiffMat = createElStiff(A, E, I, l);
	double *globStiffMat = assembleGlobStiffB(elStiffMat, rs, Nx);

	int i = 0;
	double *cenDisp = new double[Nt];
	// TIME LOOP //
	for (double n = 0; n < T; n+= delT)
	{
		// build the coeffients matrix
		double *Keff = buildKeff(rs, beta, delT, massMat, globStiffMat);
 		// build the right hand side term in brackets
		double *S = buildS(rs, beta, delT, un, unDot, unDotDot);
		// carry out the matrix-vector multiplication 
		double *MS = buildMs(rs, massMat, S);
		// build Fn+1
		double *Fnplus1 = buildFnplus1(rs, qx, qy+delqy, fy+delfy, l, Nx, delT*delT);
		//add Fnplus1 and MS to complete the right hand side of the equation
		double *B = buildBQ3(rs, Fnplus1, MS);

		//solve the system
	    int info = 0;
	    int* ipiv = new int[rs];
		// Keff needs to have extra zeros
		F77NAME(dgbsv) (rs, 4, 4, 1, Keff, 13, ipiv, B, rs, info);
		// copy the solution into unp1
		F77NAME(dcopy) (rs, B, 1.0, unp1, 1);
		// copy un before overwriting with Un+1
		F77NAME(dcopy) (rs, un, 1.0, unCopy, 1);
		// set Un = Unplus1
		F77NAME(dcopy) (rs, unp1, 1.0, un, 1);
		// calculate unDotDotp1
		double *unDotDotp1 = buildUnDotDotp1(rs, beta, delT, un, unCopy, unDot, unDotDot);
		// copy unDotDot before overwriting with unDotDotp1
		F77NAME(dcopy) (rs, unDotDot, 1.0, unDotDotCopy, 1);
		// set unDotDot = unDotDotp1
		F77NAME(dcopy) (rs, unDotDotp1, 1.0, unDotDot, 1);
		// calculate unDotp1
		double *unDotp1 = buildUnDotp1(rs, gamma, delT, unDot, unDotDotCopy, unDotDot);
		F77NAME(dcopy) (rs, unDotp1, 1.0, unDot, 1);

		// update qy and fy
		if (n < Tt)
		{
   			qy = qy + delqy;
			fy = fy + delfy;
		}
		else{}

		cenDisp[i] = un[int((rs)/2)];
		i += 1;	
	}

	std::ofstream solfile;
	solfile.open("solution_3.txt");
	for (int i = 0; i < rs; ++i)
		{
			solfile << un[i];
			solfile << "\n";
		}
	solfile.close();

	// returning array of central node displacement
	return cenDisp;

}

double *solveDynExpPar(double A, double E, double I, double l, int rs, int Nx, double qx, double qy, double fy, double T, int Nt, double alpha, double delT, double rho, int rank, int size)
{
	int nLoc = (((3*Nx)-1)/2)-1;
	// std::cout << "rank = " << rank << " size = " << size << std::endl;

	int rs_p;
	if (size == 2)
	{
		rs_p = (rs/2)+5;
	}
	else if (size == 1)
	{
		rs_p = rs;
	}
	double *un = new double[rs]();
	double *unminus1 = new double[rs]();
	double *un_term = new double[rs]();
	double *unminus1_term = new double[rs]();
	// for (int i = 0; i < rs; ++i)
	// {
	// 	un[i] = 0.0;
	// 	unminus1[i] = 0.0;
	// 	un_term[i] = 0.0;
	// 	unminus1_term[i] = 0.0;
	// }
	// define step in forces applied
	double delqy = qy/(Nt-1);
	double delfy = fy/(Nt-1);
	// reset forces before time loop starts
	qy = 0;
	fy = 0;

	// build local and global stiffness matrix
	double *elStiffMat = createElStiff(A, E, I, l);
	double *globStiffMat = assembleGlobStiffB(elStiffMat, rs, Nx);
	double *locStiffMat = assembleGlobStiffBP(elStiffMat, nLoc, rs_p, rank); 
	locStiffMat = assembleGlobStiffBw(locStiffMat, rs_p, Nx);
	globStiffMat = assembleGlobStiffBw(globStiffMat, rs, Nx);
	double *locMassMat = assembleMassB(alpha, l, rho, A, rs_p, Nx);

	// created buffers for passing between processe
	double *soln = new double[rs_p];
	double *bufSol = new double[rs_p];
	double *bufUn = new double[rs_p];
	double *bufUnminus1 = new double[rs_p];
	double *unplus1 = new double[rs];

	// create local arrays for un and unminus1
	double *locUn = new double(rs_p);
	double *locUnminus1 = new double(rs_p);

	// TIME LOOP //
	for (double n = 0; n < T; n+= delT)
	{		
		// split up un and unminus1 into local values
		locUn = split(un, rs, rs_p, rank);
		locUnminus1 = split(unminus1, rs, rs_p, rank);

		// build the local un and unminus1 terms
		double *un_term = buildUn(rs_p, locStiffMat, locMassMat, locUn, delT);
		double *unminus1_term = buildUnm(rs_p, locMassMat, locUnminus1);

		// build global fn and split for local processes
		double *fn = buildFn(rs, qx, qy, fy, l, Nx, delT*delT);
		double *locForce = split(fn, rs, rs_p, rank);
		// build the right hand side of the equation
		double *b = buildB(rs_p, locForce, un_term, unminus1_term);

		// Solve the Ax = b problem taking advantage of the diagonal nature of M
	    for(int i = 0; i < rs_p; ++i)
	    {
	    	soln[i] = (1/locMassMat[i])*(b[i]);
	    }
	    // updating values for when there is only 1 process
	    if (rank == 0 and size == 1)
	    {
	   		unplus1 = soln;
	    	un = locUn;
	    	unminus1 = locUnminus1;

	    	F77NAME(dcopy) (rs, un, 1.0, unminus1, 1);
			F77NAME(dcopy) (rs, unplus1, 1.0, un, 1);
	    }

	    // sending solution from rank 1 to rank 0
	    if(rank == 1 and size > 1)
	    {
	    	bufSol = soln;
	    	bufUn = locUn;
	    	bufUnminus1 = locUnminus1;
	    	MPI_Send(&bufSol[0], rs_p, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	    	MPI_Send(&bufUn[0], rs_p, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
	    	MPI_Send(&bufUnminus1[0], rs_p, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
	    }
	    // recieving the solution at rank 0
	    else if (rank == 0 and size > 1)
	    {
	    	MPI_Recv(&bufSol[0], rs_p, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    	MPI_Recv(&bufUn[0], rs_p, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    	MPI_Recv(&bufUnminus1[0], rs_p, MPI_DOUBLE, 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	    MPI_Barrier(MPI_COMM_WORLD); 

	    if(rank == 0 and size >1)
	    {	
	    	// combine rank 0 and rank 1 solutions on rank 0
	    	unplus1 = combine(soln, bufSol, rs_p, rs);
	    	un = combine(locUn, bufUn, rs_p, rs);
	    	unminus1 = combine(locUnminus1, bufUnminus1, rs_p, rs);

	    	// update un and unminus1 values
	    	F77NAME(dcopy) (rs, un, 1.0, unminus1, 1);
			F77NAME(dcopy) (rs, unplus1, 1.0, un, 1);

			// send full solution from rank 0 to rank 1 for splitting at the top of time loop
	    	MPI_Send(&unplus1[0], rs, MPI_DOUBLE, 1, 4, MPI_COMM_WORLD);
	    	MPI_Send(&un[0], rs, MPI_DOUBLE, 1, 5, MPI_COMM_WORLD);
	    	MPI_Send(&unminus1[0], rs, MPI_DOUBLE, 1, 6, MPI_COMM_WORLD);
	    }
	    // recieve full solution at rank 1
	    else if(rank == 1 and size > 1)
	    {
	    	MPI_Recv(&unplus1[0], rs, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    	MPI_Recv(&un[0], rs, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    	MPI_Recv(&unminus1[0], rs, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);

		// // update forces applied to the beam
   		qy = qy + delqy;
		fy = fy + delfy;

		MPI_Barrier(MPI_COMM_WORLD);
 	}

	return un;   
}

double *solveDynImpPar(double A, double E, double I, double l, int rs, int Nx, double qx, double qy, 
                       double fy, double T, int Nt, double alpha, double delT, double rho, int rank, int size)
{
	double beta = 0.25;
	double gamma = 0.5;
	int globrs = rs+1;
	int locrs = globrs/size;
	if (globrs%4 != 0 and size == 4 )
	{
		double *errPtr = 0;
		if (rank == 0)
		{
			std::cerr <<"ERROR: global size [((Nx-1)*3)+1] of problem must be divisible by 4" << std::endl;
		}
		return errPtr;
	}

	// defining pointers for local matrices
	double *locKeff = new double[locrs*17];
	double *locKeffCopy = new double[locrs*17];
	double *locB = new double[locrs];
	double *soln = new double[globrs];

	double delqy = qy/(Nt-1);
	double delfy = fy/(Nt-1);
	qy = 0;
	fy = 0;

	double *un = new double[rs]();
	double *unCopy = new double[rs]();
	double *unDot = new double[rs]();
	double *unDotDot = new double[rs]();
	double *unDotDotCopy = new double[rs]();
	double *unp1 = new double[rs]();
	double *unDotp1 = new double[rs]();
	double *unDotDotp1 = new double[rs]();
	// set initial conditions to 0
	// for (int i = 0; i < rs; ++i)
	// {
	// 	un[i] = 0.0;
	// 	unCopy[i] = 0.0;
	// 	unDot[i] = 0;
	// 	unDotDot[i] = 0;
	// 	unDotDotCopy[i] = 0;
	// 	unp1[i] = 0.0;
	// 	unDotp1[i] = 0;
	// 	unDotDotp1[i] = 0;
	// }

	const int nb   = locrs;		// Block size
    const int kl   = 4;   			// Number of lower diagonals
    const int ku   = 4;   			// Number of upper diagonals
    const int lda  = 1 + 2*kl + 2*ku; // Leading dimension (num of rows)
    // const int ldb  = nb;
    const int nrhs = 1;
    const int ja   = 1;             // Offset index in global array (col)
    const int ib   = 1;             // Offset index in global array (row)
    const int lwork= (nb+ku)*(kl+ku)+6*(kl+ku)*(kl+2*ku) + std::max(nrhs*(nb+2*kl+4*ku), 1);
    int info = 0;

    int nrow = 1;
    int ncol = 4;
    char order = 'R';
	int ctx;
    int mype;
    int npe;
    int myrow;
    int mycol;
    Cblacs_pinfo(&mype, &npe);
    Cblacs_get( 0, 0, &ctx );
    Cblacs_gridinit( &ctx, &order, 1, npe );
    Cblacs_gridinfo( ctx, &nrow, &ncol, &myrow, &mycol);

	double* wk   = new double[lwork];   // Local workspace
    int* ipiv = new int[lda*nb];         // Local pivoting array

	int desca[7];
    desca[0] = 501;         // Type is a banded matrix 1-by-P
    desca[1] = ctx;         // Context
    desca[2] = globrs;           // Problem size
    desca[3] = locrs;          // Blocking
    desca[4] = 0;           // Process row/column
    desca[5] = lda;         // Local size
    desca[6] = 0;           // Reserved

    int descb[7];
    descb[0] = 502;         // Type is a banded matrix P-by-1 (RHS)
    descb[1] = ctx;         // Context
    descb[2] = globrs;           // Problem size
    descb[3] = locrs;          // Blocking
    descb[4] = 0;           // Process row/column
    descb[5] = locrs;          // Local size
    descb[6] = 0;           // Reserved

    // Build global mass matrix
    double *massMat =  assembleMassB(alpha, l, A, rho, rs, Nx);
	massMat = addColumnB(rs, massMat);

	// build global stiffness matrix
	double *elStiffMat = createElStiff(A, E, I, l);
	// arrays pass rank into them so they are only created on rank 0
	double *globStiffMat = assembleGlobStiffBP(elStiffMat, rs, Nx, 4, rank);
	globStiffMat = addZeros(globStiffMat, rs, Nx, rank);

	// build Keff matrix on rank 0
   	double *globKeff = buildKeffP(rs, beta, delT, massMat, globStiffMat, 17, rank);
	double *Keff = addColumnKP(rs, globKeff, 17, rank);

	// scatter global Keff between other processes
	MPI_Scatter(&Keff[0], locrs*17, MPI_DOUBLE, locKeff, locrs*17, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double *S = new double[rs];
	double *MS = new double[rs];
	double *Fnplus1 = new double[rs];
	double *B = new double[rs];

	for (double n = 0; n < T; n+= delT)
	{
		// reset local Keff as pdgbsv will alter the matrix after it finishes
		F77NAME(dcopy) (locrs*17, locKeff, 1.0, locKeffCopy, 1);

		// building the right hand side of the problem and scattering among processes
		S = buildS0(rs, beta, delT, un, unDot, unDotDot, rank);
		MS = buildMs0(rs, massMat, S, rank);
		Fnplus1 = buildFnplus10(rs, qx, qy+delqy, fy+delfy, l, Nx, delT*delT, rank);
		B = buildB0(rs, Fnplus1, MS, rank);
		B = addColumnB0(rs, B, rank);
		MPI_Scatter(&B[0], locrs, MPI_DOUBLE, locB, locrs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		// solving the Ax=b problem using pdgbsv
	    F77NAME(pdgbsv) (globrs, kl, ku, nrhs, locKeffCopy, ja, desca, ipiv, locB, ib, descb, wk, lwork, &info);

	    // gather solution into soln on process 0
		MPI_Gather(locB, locrs, MPI_DOUBLE, soln, locrs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		// copy un before overwriting with Un+1
		if (rank == 0)
		{
			F77NAME(dcopy) (rs, un, 1.0, unCopy, 1);

			// copy solution to un
			F77NAME(dcopy) (rs, soln, 1.0, un, 1);

			unDotDotp1 = buildUnDotDotp1(rs, beta, delT, un, unCopy, unDot, unDotDot);

			// copy unDotDot before overwriting with unDotDotp1
			F77NAME(dcopy) (rs, unDotDot, 1.0, unDotDotCopy, 1);

			// set unDotDot = unDotDotp1
			F77NAME(dcopy) (rs, unDotDotp1, 1.0, unDotDot, 1);

			// calculate unDotp1
			unDotp1 = buildUnDotp1(rs, gamma, delT, unDot, unDotDotCopy, unDotDot);

			// copy unDotp1 into unDot
			F77NAME(dcopy) (rs, unDotp1, 1.0, unDot, 1);
		}
		qy = qy + delqy;
		fy = fy + delfy;
		
	}

	// // Free CBLACS context
    Cblacs_gridexit( ctx );

	return un;
}
