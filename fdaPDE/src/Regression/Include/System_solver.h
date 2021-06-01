#ifndef __SYSTEM_SOLVER_H__
#define __SYSTEM_SOLVER_H__

#include <memory>
#include <type_traits>
#include "../../FdaPDE.h"

// Base class for a generic SparseLU solver of a system Mx=b
// with M square matrix composed of four blocks: M = | NW | NE |
//													 | SW | SE |
//
class BaseSolver
{
protected:
	Eigen::SparseLU<SpMat> Mdec;
	bool decomposed = false;
	inline void system_factorize(const SpMat& M) { Mdec.compute(M); decomposed = true; }
public:
	BaseSolver() = default;
	explicit BaseSolver(const SpMat& M) { system_factorize(M); };
	virtual SpMat buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE); // Assemble the system matrix starting from the four blocks
	virtual void compute(const SpMat& M) { system_factorize(M); }
	//Solve the system
	virtual MatrixXr system_solve(const MatrixXr& b) const;
	virtual MatrixXr system_solve(const SpMat& M, const MatrixXr& b);
	virtual void compute(const Real lambda_, const UInt nnodes) { return; }
};


// The bottom-right block of the system matrix (i.e. the FE mass matrix) is lumped
// Solving the system Ml x = b, with Ml = | NW  | NE  |
//										  | SW  | SE' |
// SE' = diag(m_ii), m_ii=sum_j (SE_ij)
//
class MassLumping : public BaseSolver
{
protected:
	VectorXr diag;    // Stores the diagonal of SE'
	void lumpMassMatrix(const SpMat& M);    // Compute the diagonal terms of SE'
	bool lumped = false;
public:
	MassLumping() : diag(VectorXr::Ones(1)) {};
	explicit MassLumping(const SpMat& M) { compute(M); };
	MassLumping(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE);
	SpMat buildSystemMatrix(const SpMat& M);	// Compute the system matrix Ml starting from M
	SpMat buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE) override; // Assemble the system matrix Ml starting from the four blocks
	inline void compute(const SpMat& M) override
	{
		if(!lumped)
			buildSystemMatrix(M);
		system_factorize(M);
	}
	using BaseSolver::system_solve;
	MatrixXr system_solve(const SpMat& M, const MatrixXr& b) override;
	using BaseSolver::compute;
};

class IdentityMass : public BaseSolver
{
	Real lambda;
public:
	IdentityMass() : lambda(1.0) {};
	explicit IdentityMass(const SpMat& M) : lambda(1.0) { compute(M); };
	SpMat buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE) override; // Assemble the system matrix starting from the four blocks
	using BaseSolver::system_solve;
	using BaseSolver::compute;
	void compute(const Real lambda_, const UInt nnodes) override { lambda = lambda_; }
};

class LumpedPreconditioner : public MassLumping
{
	//Real lambda;
public:
	LumpedPreconditioner() = default; //: lambda(1.0) {};
	explicit LumpedPreconditioner(const SpMat& M) { compute(M); }; // : lambda(1.0) { compute(M); };
	SpMat buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE) override; // Assemble the system matrix starting from the four blocks
	MatrixXr system_solve(const MatrixXr& b) const override;
	MatrixXr system_solve(const SpMat& M, const MatrixXr& b) override { compute(M); return system_solve(b); }
	using BaseSolver::compute;
	//void compute(const Real lambda_, const UInt nnodes) override { lambda = lambda_; }
};


// Base class for a generic SparseLU solver with a diagonal preconditioner
// Solves the system DMx=Db, where D is a diagonal matrix
//
class BaseDiagPreconditioner : public BaseSolver
{
protected:
	VectorXr prec;  // Stores the diagonal of the preconditioner D
	bool initialized = false;  // Indicates if the preconditioner has been initialized
public:
	BaseDiagPreconditioner() : prec{VectorXr::Ones(1)} {};
	inline VectorXr preconditioner() const { return prec; }   // Get the preconditioner diagonal
	MatrixXr preconditionRHS(const MatrixXr& b) const;		  // Compute D*b
	inline MatrixXr system_solve(const MatrixXr& b) const override { return BaseSolver::system_solve(preconditionRHS(b)); }
	MatrixXr system_solve(const SpMat& M, const MatrixXr& b) override;
	virtual void compute(const SpMat& M)
	{
		if (initialized)
			system_factorize(prec.asDiagonal() * M);
		else
			system_factorize(M);
	}
	using BaseSolver::compute;
};


// The diagonal preconditioner can be set by the user
class DiagonalPreconditioner : public BaseDiagPreconditioner
{
public:
	explicit DiagonalPreconditioner(const VectorXr& p) { prec = p; initialized = true; };
	DiagonalPreconditioner(const SpMat& M, const VectorXr& p) { prec = p; initialized = true; compute(M); }; // Set the preconditioner diagonal
	inline void preconditioner(const VectorXr& p) { prec = p; initialized = true; }
};

// Solves the system DMDx'=Db, x = Dx'
// with D=diag(1,...,1,1/sqrt(lambda),...,1/sqrt(lambda))
//
class LambdaPreconditioner : public BaseDiagPreconditioner
{
	Real lambda;
public:
	LambdaPreconditioner() : lambda(1.0) {};
	LambdaPreconditioner(const Real lambda_, const UInt nnodes) : lambda(lambda_) { compute(lambda_, nnodes); };
	Real getlambda() { return lambda; }
	using BaseSolver::buildSystemMatrix;
	MatrixXr system_solve(const SpMat& M, const MatrixXr& b) override;
	MatrixXr system_solve(const MatrixXr& b) const override;
	void compute(const SpMat& M) override { system_factorize(M); }
	void compute(const Real lambda_);
	void compute(const Real lambda_, const UInt nnodes) override;
 };


// Base class for a generic SparseLU solver with a block diagonal preconditioner with two blocks
// Solves the system BMx=Bb, where B is a block diagonal matrix:  B = | NWblock^-1 |     0      |
//																	  |     0      | SEblock^-1 |
//
class BlockPreconditioner : public BaseSolver
{
protected:
	SpMat SEblock;
	Eigen::SparseLU<SpMat> SEdec;
	bool initialized = false;

public:
	BlockPreconditioner();
	SpMat buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE) override;
	SpMat preconditioner() const;
	MatrixXr preconditionRHS(const MatrixXr& b) const;
	MatrixXr system_solve(const MatrixXr& b) const override;
	MatrixXr system_solve(const SpMat& M, const MatrixXr& b) override;
	using BaseSolver::compute;
};

#endif