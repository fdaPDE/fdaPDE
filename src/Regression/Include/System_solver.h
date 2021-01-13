#ifndef __SYSTEM_SOLVER_H__
#define __SYSTEM_SOLVER_H__

#include <memory>
#include <type_traits>
#include "../../FdaPDE.h"

class SystemSolver
{
public:
	SpMat buildSystemMatrix(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE);
	MatrixXr system_solve(SpMat& M, MatrixXr& b);
	MatrixXr system_solve(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE, MatrixXr& b);
};

class MassLumping : public SystemSolver
{
	SpMat diag;
	bool isdiagInitialized = false;
	void builddiag(SpMat& M);

public:
	MassLumping(SpMat& M) { this->builddiag(M);};
	SpMat buildSystemMatrix(SpMat& M);
	SpMat buildSystemMatrix(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE);
	inline MatrixXr system_solve(MatrixXr& b) {return VectorXr(diag).cwiseInverse().asDiagonal()*b;}
	MatrixXr system_solve(SpMat& M, MatrixXr& b);
	MatrixXr system_solve(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE, MatrixXr& b);
};

class diagonalPreconditioner : public SystemSolver
{
	VectorXr diag;
	bool isdiagInitialized = false;

public:
	diagonalPreconditioner(VectorXr& prec) : diag(prec) { isdiagInitialized = true; };
	inline void compute(VectorXr& prec) { diag = prec; isdiagInitialized = true; }
	SpMat preconditionMatrix(SpMat& M);
	SpMat preconditionMatrix(SpMat& M, VectorXr& prec);
	inline SpMat buildSystemMatrix(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE) { return this->preconditionMatrix(SystemSolver::buildSystemMatrix(NW, SE, SW, NE)); }
	inline SpMat buildSystemMatrix(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE, VectorXr& prec)
	{
		return this->preconditionMatrix(SystemSolver::buildSystemMatrix(NW, SE, SW, NE),b);
	}
	MatrixXr system_solve(SpMat& M, MatrixXr& b, bool preconditioned);
	MatrixXr system_solve(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE, MatrixXr& b, bool preconditioned);
};

class JacobiPreconditioner : public diagonalPreconditioner
{
public:
	// eliminare il costruttore con SpVec e compute con VectorXr&
	void compute(SpMat& M);
	JacobiPreconditioner(SpMat& M) { compute(M); };
};

class blocksPreconditioner : public SystemSolver
{
	SpMat NWblock;
	SpMat SEblock;
	SpMat SWblock;
	SpMat NEblock;

	Eigen::SparseLU<SpMat> NWdec;
	Eigen::SparseLU<SpMat> SEdec;

	bool decomposed = false;

	UInt nnodes;

public:
	blocksPreconditioner(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE) : NWblock(NW), SEblock(SE), SWblock(SW), NEblock(NE)
	{
		nnodes = NW.outerSize();
		NWdec.compute(NW);
		SEdec.compute(SE);
		decomposed = true;
	};

	SpMat buildSystemMatrix();
	SpMat buildSystemMatrix(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE);
	MatrixXr preconditionRHS(MatrixXr& b);
	MatrixXr system_solve(SpMat& M, MatrixXr& b, bool preconditioned);
	MatrixXr system_solve(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE, MatrixXr& b, bool preconditioned);
};

#endif