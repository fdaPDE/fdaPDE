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
	// Factorization of the system matrix
	Eigen::SparseLU<SpMat> Mdec;
	
	void system_factorize(const SpMat& M)
	{
		Mdec.compute(M);
		decomposed = true; 
	}

	//// Time-dependent elements
	//std::shared_ptr<SpMat> Ptk = NULL;	// additive block corresponding to separate penalization
	//std::shared_ptr<SpMat> LR0k = NULL;	// additive block corresponding to parabolic penalization
	//Real lambdaT = 1.0;					// time smoothing parameter
	//bool parabolic = false;
	//bool timeDependent = false;
	bool iterative = false;

	UInt M_;
	Real lambda = 1.0;      // The space smoothing parameter

	//Only used for BlockPreconditioner
	SpMat SEblock;

public:
	BaseSolver() = default;
	BaseSolver(const SpMat & M) { compute(M); };

	bool decomposed = false;

	//// A method that stores the pointers to the additive blocks corresponding to time penalization
	//virtual void addTimeCorrection(const std::shared_ptr<SpMat> Mat, Real lambda, bool flagParabolic)
	//{
	//	parabolic = flagParabolic;
	//	if (parabolic)
	//		LR0k = Mat;
	//	else
	//		Ptk = Mat;
	//	lambdaT = lambda;
	//	timeDependent = true;
	//}

	// A method assembling the system matrix starting from the space-dependent blocks
	virtual SpMat assembleMatrix(const SpMat & DMat, const SpMat & R0, const SpMat & R1, Real lambdaS);

	// A method assembling the system matrix starting from the four blocks
	SpMat buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE);

	// Factorize the system matrix
	virtual void compute(const SpMat & M) { system_factorize(M); }

	// Solve the system
	virtual MatrixXr system_solve(const MatrixXr & b) const;

	// Iterative method flag setter
	void setIterative(bool flag) { iterative = flag; }

	void setM(UInt M) { M_ = M; }

	// Only used for BlockPreconditioner: A method storing the SE block of the original matrix
	virtual void setSEblock(const SpMat & SE) { SEblock = SE; }

	// A method setting the space smoothing parameter
	void set_lambdaS(Real lambda_) { lambda = lambda_; }
};


// The bottom-right block of the system matrix (i.e. the FE mass matrix) is lumped
// Solving the system Ml x = b, with Ml = | NW  | NE  |
//										  | SW  | SE' |
// SE' = diag(m_ii), m_ii=sum_j (SE_ij)
//
class MassLumping : public BaseSolver
{
protected:
	SpMat lumpMassMatrix(const SpMat& M);    // Compute the diagonal terms of SE'
public:
	MassLumping() = default;
	MassLumping(const SpMat& M) { compute(M); };
	SpMat assembleMatrix(const SpMat& DMat, const SpMat& R0, const SpMat& R1, Real lambdaS) override;
	using BaseSolver::system_solve;
	using BaseSolver::compute;
	//using BaseSolver::addTimeCorrection;
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
	BaseDiagPreconditioner(const SpMat& M) { compute(M); };
	inline VectorXr preconditioner() const { return prec; }   // Get the preconditioner diagonal
	virtual MatrixXr preconditionRHS(const MatrixXr& b) const;		  // Compute D*b
	inline MatrixXr system_solve(const MatrixXr& b) const override { return BaseSolver::system_solve(preconditionRHS(b)); }
	virtual void compute(const SpMat& M)
	{
		if (initialized)
			system_factorize(prec.asDiagonal() * M);
		else
			system_factorize(M);
	}
	using BaseSolver::assembleMatrix;
	//using BaseSolver::addTimeCorrection;
};


// The diagonal preconditioner can be set by the user
class DiagonalPreconditioner : public BaseDiagPreconditioner
{
public:
	explicit DiagonalPreconditioner(const VectorXr& p) { prec = p; initialized = true; };
	DiagonalPreconditioner(const SpMat& M) { compute(M); };
	DiagonalPreconditioner(const SpMat& M, const VectorXr& p) { prec = p; initialized = true; compute(M); }; // Set the preconditioner diagonal
	inline void preconditioner(const VectorXr& p) { prec = p; initialized = true; }
	using BaseSolver::assembleMatrix;
	using BaseDiagPreconditioner::system_solve;
	using BaseDiagPreconditioner::compute;
	//using BaseSolver::addTimeCorrection;
};

// Solves the system DMDx'=Db, x = Dx'
// with D=diag(1,...,1,1/sqrt(lambda),...,1/sqrt(lambda))
//
class LambdaPreconditioner : public BaseDiagPreconditioner
{
protected:
	Real lambda = 1.0;
public:
	LambdaPreconditioner() : lambda(1.0) {};
	LambdaPreconditioner(const SpMat& M) { compute(M); };
	LambdaPreconditioner(Real lambda_) : lambda(lambda_) {};
	Real getlambda() const { return lambda; }
	SpMat assembleMatrix(const SpMat& DMat, const SpMat& R0, const SpMat& R1, Real lambdaS) override;
	MatrixXr preconditionRHS(const MatrixXr& b) const override;
	MatrixXr system_solve(const MatrixXr& b) const override;
	using BaseSolver::compute;
	//using BaseSolver::addTimeCorrection;
 };


// Base class for a generic SparseLU solver with a block diagonal preconditioner with two blocks
// Solves the system BMx=Bb, where B is a block diagonal matrix:  B = | NWblock^-1 |     0      |
//																	  |     0      | SEblock^-1 |
//
class BlockPreconditioner : public BaseSolver
{
protected:
	Real lambda = 1.0;			// The space smoothing parameter

	SpMat SEblock;			// R0_lambda
	SpMat SWblock;			// R0^-1 * R1
	Eigen::SparseLU<SpMat> SEdec;	// LU decomposition of R0
	bool initialized = false;

public:
	BlockPreconditioner();
	BlockPreconditioner(const SpMat& M) { compute(M); };
	//void addTimeCorrection(const std::shared_ptr<SpMat> Mat, Real lambda, bool flagParabolic) override
	//{
	//	// If the time constant does not change, no need for recomputing the factorization and blocks
	//	if (lambda == lambdaT)
	//		recompute = false;
	//	else
	//		recompute = true;
	//	//BaseSolver::addTimeCorrection(Mat, lambda, flagParabolic);
	//}
	SpMat assembleMatrix(const SpMat& DMat, const SpMat& R0, const SpMat& R1, Real lambdaS) override;
	SpMat preconditioner() const; // Get the preconditioner 
	MatrixXr preconditionRHS(const MatrixXr& b) const;
	MatrixXr system_solve(const MatrixXr& b) const override;
	using BaseSolver::compute;
	void setSEblock(const SpMat& SE) override { SEblock = SE; SEdec.compute(SE); initialized = true; }
};

#endif