#include <iostream>
#include "functions.h"
#include <stdlib.h>
#include <typeinfo>
#include <fstream>
#include <mpi.h>

/*
L = 10
Nx = 501
A = 0.012
I = 0.0000144
E = 210000000000
T = 1
Nt = 11
rho = 7850
*/


// Print a matrix for debugging //
void print_matrix(double *A, int SIZEx, int SIZEy) {
	for (int i = 0; i < SIZEy; i++) {
		for (int j = 0; j < SIZEx; j++) {
			std::cout << A[(j*SIZEy)+i] << "  ";
		}
		std::cout << std::endl;
	}
}

// main //
int main(int argc, char* argv[])
{
	// define parameters
	double L, A, I, E, l, qx, qy, fy, T, delT, alpha, rho;
	int Nx, size, rs, Question, Nt;
	Question = atof(argv[1]);
	L = atof(argv[2]);
	Nx = atoi(argv[3]);
	A = atof(argv[4]);
	I = atof(argv[5]);
	E = atof(argv[6]);
	T = atof(argv[7]);
	Nt = atoi(argv[8]);
	rho = atof(argv[9]);

	alpha = 1.0/24;
	delT = T/(Nt);
	l = L/(Nx-1);
	size = 3*(Nx);
	rs = size-6;
	qx = 0;
	qy = -1000;
	fy = -1000;

	/////////////////////////////////////////////////////////////////////////
	// validating input
	if (Nx%2 == 0)
	{
		std::cerr << "ERROR: Number of nodes (Nx) must be an odd number" << std::endl;
		return 2;
	}
	////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	// Running correct code dependent on the value of parameter 'Question'
	if (Question == 1)
	{	
		// solve problem
		double *soln = solveStatic(A, E, I, l, rs, Nx, qx, qy, fy);
		
		// write to file
		std::ofstream solfile;
		solfile.open("solution_1.txt");
		for (int i = 0; i < rs; ++i)
		{
			solfile << soln[i];
			solfile << "\n";
		}
		solfile.close();
		// std::cout << "Q1 finished" << std::endl;
		return 0;
	}
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	if (Question == 2)
	{
		int itr = 0;
		double *amps = new double[100];
		double *time = new double[100];
		double *oscArr = new double[Nt];
		double *soln = new double[Nt];
		double min;
		double max;
		double amp;
		int k;

		// running solver to return and write to file the displacements of the central node
		double Tt = T;
		double *cenDisp = solveDynExp(A, E, I, l, rs, Nx, qx, qy, fy, T, Nt, alpha, delT, rho, Tt);
		std::ofstream dispfile;
		dispfile.open("solution_2disp.txt");
		for (int i = 0; i < Nt; ++i)
		{
			dispfile << cenDisp[i];
			dispfile << "\n";
			// std::cout << un[i] << std::endl;
		}
		dispfile.close();
		///////////////////////////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////////////////////
		// COMMENTED OUT CODE BLOCK USED TO LOOP THROUGH DIFFERENT LOADING SPEEDS TO PRODUCE AMPLITUDE OF OSCILLATIONS
		///////////////////////////////////////////////////////////////////////////////////
		// for (float ldParam = 0.02; ldParam < 1.0; ldParam += 0.02)
		// {
		//	// Tt is the loading time
		// 	double Tt = T*ldParam;
		// 	soln = solveDynExp(A, E, I, l, rs, Nx, qx, qy, fy, T, Nt, alpha, delT, rho, Tt);
			
		//	// fill oscArr with the central displacement with values after full load has been applied	
		// 	k = 0;
		// 	for (int i = Nt*ldParam; i < Nt; ++i)
		// 	{
		// 		oscArr[k] = 0.0;
		// 		oscArr[k] = soln[i];
		// 		k += 1; 
		// 	}
		// 	min = 10.0;
		// 	max = -10.0;
		//  // find the maximum and minimum values to calculate amplitude
		// 	for (int p = 1; p < int(Nt-(Nt*ldParam)); ++p)
		// 	{
		// 		if (oscArr[p] > oscArr[p-1] and oscArr[p] > max)
		// 		{
		// 			max = oscArr[p];
		// 		}
		// 		if (oscArr[p] < oscArr[p-1] and oscArr[p] < min)
		// 		{
		// 			min = oscArr[p];
		// 		}
		// 	}
		// 	// append aplitude to amp array
		// 	amp = 0.5*(max-min);
		// 	amps[itr] = amp;
		// 	time[itr] = Tt;
		// 	itr += 1;
		// }
		// // writing amplitudes and times to file
		// std::ofstream ampfile;
		// ampfile.open("solution_2amp.txt");
		// for (int i = 0; i < 100; ++i)
		// {
		// 	ampfile << amps[i] << "\t" << time[i];
		// 	ampfile << "\n";
		// }
		// ampfile.close();
		///////////////////////////////////////////////////////////////////////////////////

		return 0;
	}

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	if (Question == 3)
	{
		int itr = 0;
		double *amps = new double[100];
		double *time = new double[100];
		double *oscArr = new double[Nt];
		double *soln = new double[Nt];
		double min;
		double max;
		double amp;
		int k;

		// running solver to return and write to file the displacements of the central node
		double Tt = T;
		double *cenDisp = solveDynImp(A, E, I, l, rs, Nx, qx, qy, fy, T, Nt, alpha, delT, rho, Tt);
		std::ofstream dispfile;
		dispfile.open("solution_3disp.txt");
		for (int i = 0; i < Nt; ++i)
		{
			dispfile << cenDisp[i];
			dispfile << "\n";
			// std::cout << un[i] << std::endl;
		}
		dispfile.close();
		///////////////////////////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////////////////////
		// COMMENTED OUT CODE BLOCK USED TO LOOP THROUGH DIFFERENT LOADING SPEEDS TO PRODUCE AMPLITUDE OF OSCILLATIONS
		//  note: takes roughly 1 minute to run
		///////////////////////////////////////////////////////////////////////////////////
		// for (float ldParam = 0.02; ldParam < 1.0; ldParam += 0.02)
		// {
		//	// Tt is the loading time
		// 	double Tt = T*ldParam;
		// 	soln = solveDynImp(A, E, I, l, rs, Nx, qx, qy, fy, T, Nt, alpha, delT, rho, Tt);
		
		//	// fill oscArr with the central displacement with values after full load has been applied	
		// 	k = 0;
		// 	for (int i = Nt*ldParam; i < Nt; ++i)
		// 	{
		// 		oscArr[k] = 0.0;
		// 		oscArr[k] = soln[i];
		// 		k += 1; 
		// 	}
		// 	min = 10.0;
		// 	max = -10.0;
		//  // find the maximum and minimum values to calculate amplitude
		// 	for (int p = 1; p < int(Nt-(Nt*ldParam)); ++p)
		// 	{
		// 		if (oscArr[p] > oscArr[p-1] and oscArr[p] > max)
		// 		{
		// 			max = oscArr[p];
		// 		}
		// 		if (oscArr[p] < oscArr[p-1] and oscArr[p] < min)
		// 		{
		// 			min = oscArr[p];
		// 		}
		// 	}
		// 	// append aplitude to amp array
		// 	amp = 0.5*(max-min);
		// 	amps[itr] = amp;
		// 	time[itr] = Tt;
		// 	itr += 1;
		// }

		// // write amplitude array to file
		// std::ofstream ampfile;
		// ampfile.open("solution_3amp.txt");
		// for (int i = 0; i < 100; ++i)
		// {
		// 	ampfile << amps[i] << "\t" << time[i];
		// 	ampfile << "\n";
		// }
		// ampfile.close();
		///////////////////////////////////////////////////////////////////////////////////

		return 0;
	}

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	if (Question == 4)
	{
		// solve the problem
		int rank;
		int size;
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (size > 2)
		{
			if (rank == 0)
			{
				std::cerr << "ERROR: number of processes (size) must be equal to 1 or 2." << std::endl;
			}
			MPI_Finalize();
			return 0;
		}
		
		double *soln = solveDynExpPar(A, E, I, l, rs, Nx, qx, qy, fy, T, Nt, alpha, delT, rho, rank, size);

		std::ofstream solfile;
		solfile.open("solution_4.txt");
		for (int i = 0; i < rs; ++i)
		{
			solfile << soln[i];
			solfile << "\n";
		}

		MPI_Finalize();
		return 0;
	}

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	if (Question == 5)
	{
		// initialise MPI
		int rank;
		int size;
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		double *soln = solveDynImpPar(A, E, I, l, rs, Nx, qx, qy, fy, T, Nt, alpha, delT, rho, rank, size);
		if (soln == 0)
			{
				MPI_Finalize();
				return 0;
			}

		std::ofstream solfile;
		if (rank == 0)
		{	
			solfile.open("solution_5.txt");
			for (int i = 0; i < rs; ++i)
			{
				solfile << soln[i];
				solfile << "\n";
			}
		}
		MPI_Finalize();
		return 0;

	}
}
