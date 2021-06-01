#ifndef __SYSTEM_SOLVER_H__
#define __SYSTEM_SOLVER_H__

#include <memory>
#include <type_traits>
#include "../../FdaPDE.h"

class VirtualSolver
{
protected:
	Eigen::SparseLU<SpMat> Mdec;
	bool decomposed = false;
public:
	VirtualSolver() = default;
	SpMat buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE) = 0;
	MatrixXr system_solve(const SpMat& M, const MatrixXr& b) = 0;
	MatrixXr system_solve(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE, const MatrixXr& b) = 0;
	MatrixXr system_solve(const MatrixXr& b) = 0;
	SpMat preconditionMatrix(const SpMat& M) = 0;
	void compute(const SpMat& M) = 0;
	MatrixXr preconditionRHS(const MatrixXr& b) = 0;
};

class SystemSolver : public VirtualSolver
{
public:
	SystemSolver() = default;
	SystemSolver(const SpMat& M)
	{
		Mdec.compute(M);
		decomposed = true;
	};
	SystemSolver(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
	{
		Mdec.compute(buildSystemMatrix(NW, SE, SW, NE));
		decomposed = true;
	};
	SpMat buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE);
	MatrixXr system_solve(const SpMat& M, const MatrixXr& b);
	MatrixXr system_solve(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE, const MatrixXr& b);
	MatrixXr system_solve(const MatrixXr& b);
};

class MassLumping : public SystemSolver
{
	VectorXr diag;
	bool isdiagInitialized = false;
	void builddiag(const SpMat& M);

public:
	MassLumping(const SpMat& M) { builddiag(M); };
	MassLumping(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE) { this->buildSystemMatrix(NW, SW, SW, NE); };
	SpMat buildSystemMatrix(const SpMat& M);
	SpMat buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE);
	MatrixXr system_solve(const MatrixXr& b);
	MatrixXr system_solve(const SpMat& M, const MatrixXr& b);
	MatrixXr system_solve(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE, const MatrixXr& b);
};

class diagonalPreconditioner : public SystemSolver
{
protected:
	VectorXr diag;
	bool isdiagInitialized = false;

public:
	diagonalPreconditioner() {diag = VectorXr::Ones(1);};
	diagonalPreconditioner(const VectorXr& prec) : diag(prec) { isdiagInitialized = true; };
	//diagonalPreconditioner(SpMat& M, VecorXr& prec) : diag(prec) { isdiagInitialized = true; Mdec.compute(M); decomposed = true; };
	inline void preconditioner(const VectorXr& prec) { diag = prec; isdiagInitialized = true; }
	inline SpMat preconditioner() const { return SpMat(diag.asDiagonal()); }
	SpMat preconditionMatrix(const SpMat& M);
	SpMat preconditionMatrix(const SpMat& M, const VectorXr& prec);
	MatrixXr preconditionRHS(const MatrixXr& b);
	SpMat buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE);
	SpMat buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE, const VectorXr& prec);
	MatrixXr system_solve(const SpMat& M, const MatrixXr& b, bool preconditioned);
	MatrixXr system_solve(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE, const MatrixXr& b, bool preconditioned);
	MatrixXr system_solve(const MatrixXr& b, bool preconditioned);
};

class JacobiPreconditioner : public diagonalPreconditioner
{
public:
	// eliminare il costruttore con SpVec e compute con VectorXr&
	void compute(const SpMat& M);
	JacobiPreconditioner(const SpMat& M) { this->compute(M); };
};

class blocksPreconditioner : public SystemSolver
{
	SpMat NWblock;
	SpMat SEblock;
	SpMat SWblock;
	SpMat NEblock;

	Eigen::SparseLU<SpMat> NWdec;
	Eigen::SparseLU<SpMat> SEdec;

	bool blocksDecomposed = false;

	UInt nnodes;

public:
	blocksPreconditioner();
	blocksPreconditioner(const SpMat& M) { compute(M); };
	blocksPreconditioner(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE) : NWblock(NW), SEblock(SE), SWblock(SW), NEblock(NE)
	{
		nnodes = NW.outerSize();
		NWdec.compute(NW);
		SEdec.compute(SE);
		decomposed = true;
	};

	void compute(const SpMat& M);
	SpMat buildSystemMatrix();
	SpMat buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE);
	MatrixXr preconditionRHS(const MatrixXr& b);
	MatrixXr system_solve(const MatrixXr& b, bool preconditioned);
	MatrixXr system_solve(const SpMat& M, const MatrixXr& b, bool preconditioned);
	MatrixXr system_solve(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE, const MatrixXr& b, bool preconditioned);
};

#endif