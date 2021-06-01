#include "../Include/System_solver.h"

// ---------- Base Solver methods ----------

SpMat BaseSolver::buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
{

	UInt nnodes = NW.outerSize();
	// Vector to be filled with the triplets used to build _coeffmatrix (reserved with the right dimension)
	std::vector<coeff> tripletAll;
	tripletAll.reserve(NW.nonZeros() + NE.nonZeros() + SW.nonZeros() + SE.nonZeros());

	// Parsing all matrices, reading the values to be put inside _coeffmatrix, coordinates according to the rules
	for (UInt k = 0; k < NW.outerSize(); ++k)
		for (SpMat::InnerIterator it(NW, k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row(), it.col(), it.value()));
		}
	for (UInt k = 0; k < SE.outerSize(); ++k)
		for (SpMat::InnerIterator it(SE, k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row() + nnodes, it.col() + nnodes, it.value()));
		}
	for (UInt k = 0; k < NE.outerSize(); ++k)
		for (SpMat::InnerIterator it(NE, k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row(), it.col() + nnodes, it.value()));
		}
	for (UInt k = 0; k < SW.outerSize(); ++k)
		for (SpMat::InnerIterator it(SW, k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row() + nnodes, it.col(), it.value()));
		}

	// Define, resize, fill and compress
	SpMat M;
	M.resize(2 * nnodes, 2 * nnodes);
	M.setZero();
	M.setFromTriplets(tripletAll.begin(), tripletAll.end());
	M.makeCompressed();

	Mdec.compute(M);
	decomposed = true;

	return M;
}

MatrixXr BaseSolver::system_solve(const SpMat& M, const MatrixXr& b)
{
	system_factorize(M);
	return Mdec.solve(b);
}

MatrixXr BaseSolver::system_solve(const MatrixXr& b) const
{
	if (!decomposed)
	{
		Rprintf("System matrix not initialized. Returining 0");
		return MatrixXr::Zero(b.rows(), b.cols());
	}
	else
		return Mdec.solve(b);
}


// ---------- Mass lumping methods ----------

MassLumping::MassLumping(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
{
	system_factorize(this->buildSystemMatrix(NW, SE, SW, NE));
}

void MassLumping::lumpMassMatrix(const SpMat& M)
{
	diag.resize(M.outerSize());

	for (UInt k = 0; k < M.outerSize(); ++k)
	{
		Real val = 0.0;
		for (SpMat::InnerIterator it(M, k); it; ++it)
			val += it.value();
		diag(k) = val;
	}

	lumped = true;
}

SpMat MassLumping::buildSystemMatrix(const SpMat& M)
{
	UInt nnodes = M.outerSize() / 2;
	return buildSystemMatrix(M.topLeftCorner(nnodes, nnodes), M.bottomRightCorner(nnodes, nnodes),
		M.bottomLeftCorner(nnodes, nnodes), M.topRightCorner(nnodes, nnodes));
}

SpMat MassLumping::buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
{
	lumpMassMatrix(SE);
	SpMat SEblock(diag.asDiagonal());
	return BaseSolver::buildSystemMatrix(NW,SEblock,SW,NE);
}

MatrixXr MassLumping::system_solve(const SpMat& M, const MatrixXr& b)
{
	compute(M);
	return BaseSolver::system_solve(b);
}

// ---------- Identity Mass methods ----------

SpMat IdentityMass::buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
{
	UInt nnodes = SW.outerSize();
	SpMat Id(nnodes, nnodes);
	Id.setIdentity();

	return BaseSolver::buildSystemMatrix(NW,lambda*Id,SW,NE);
}

// ---------- Lumped Mass Preconditioner methods ----------

SpMat LumpedPreconditioner::buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
{
	lumpMassMatrix(SE); // / lambda);
	UInt nnodes = SE.outerSize();
	SpMat Id(nnodes, nnodes);
	Id.setIdentity();

	return BaseSolver::buildSystemMatrix(NW, Id, diag.cwiseInverse().asDiagonal()*SW, NE);
}

MatrixXr LumpedPreconditioner::system_solve(const MatrixXr& b) const
{
	VectorXr prec = VectorXr::Ones(2 * diag.rows());
	prec.bottomRows(diag.rows()) = diag.cwiseInverse();
	return MassLumping::system_solve(prec.asDiagonal() * b);
}

// ---------- Base diagonal preconditioner methods ----------

MatrixXr BaseDiagPreconditioner::preconditionRHS(const MatrixXr& b) const
{
	if (initialized)
		return prec.asDiagonal() * b;

	Rprintf("Preconditioner not initialized. Using identity");
	return b;
}


MatrixXr BaseDiagPreconditioner::system_solve(const SpMat& M, const MatrixXr& b)
{
	compute(M);
	return system_solve(b);
}


// ---------- Lambda preconditioner methods ----------

void LambdaPreconditioner::compute(const Real lambda_, const UInt nnodes)
{
	lambda = lambda_;
	prec.resize(2*nnodes);
	prec = VectorXr::Ones(2*nnodes);
	prec.bottomRows(nnodes) = prec.bottomRows(nnodes) / sqrt(lambda_);
	initialized = true;
}

void LambdaPreconditioner::compute(const Real lambda_)
{
	if (lambda == lambda_)
		return;

	lambda = lambda_;
	if (initialized)
	{
		UInt nnodes = prec.cols() / 2;
		prec.bottomRows(nnodes) = VectorXr::Ones(nnodes) / sqrt(lambda_);
	}
}

MatrixXr LambdaPreconditioner::system_solve(const SpMat& M, const MatrixXr& b)
{
	compute(M);
	return system_solve(b);
}

MatrixXr LambdaPreconditioner::system_solve(const MatrixXr& b) const
{
	if (initialized)
		return prec.asDiagonal() * BaseDiagPreconditioner::system_solve(b);
	return BaseDiagPreconditioner::system_solve(b);
}

// ---------- Abstract blocks digonal preconditioner methods ----------
BlockPreconditioner::BlockPreconditioner()
{
	SEblock.resize(1,1);
	SEblock.setIdentity();
	SEdec.compute(SEblock);
}

SpMat BlockPreconditioner::buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
{
	UInt nnodes = SE.outerSize();
	SEblock = SE;
	SEdec.compute(SE);
	initialized = true;

	SpMat Id(nnodes, nnodes);
	Id.setIdentity();

	return BaseSolver::buildSystemMatrix(NW, Id, SEdec.solve(SW), NE);
}

SpMat BlockPreconditioner::preconditioner() const
{
	SpMat M;
	if (!initialized)
	{
		Rprintf("Preconditioner not initalized. Returning identity");
		M.resize(1, 1);
		M.setIdentity();
		return M;
	}
	else
	{
		UInt nnodes = SEblock.outerSize();
		MatrixXr M = MatrixXr::Identity(nnodes, nnodes);
		M.bottomRightCorner(nnodes, nnodes) = SEdec.solve(MatrixXr::Identity(nnodes, nnodes));
		return M.sparseView();
	}
}

MatrixXr BlockPreconditioner::preconditionRHS(const MatrixXr& b) const
{
	if (!initialized)
	{
		Rprintf("Preconditioner not initialized. Using identity");
		return b;
	}
	
	UInt nnodes = SEblock.outerSize();
	MatrixXr rhs(b);
	rhs.bottomRows(nnodes) = SEdec.solve(b.bottomRows(nnodes));
	return rhs;
}

MatrixXr BlockPreconditioner::system_solve(const SpMat& M, const MatrixXr& b)
{
	if (initialized)
		return BaseSolver::system_solve(preconditioner() * M, preconditionRHS(b));

	Rprintf("Preconditioner not initialized. Using identity");
	return BaseSolver::system_solve(M, b);
}

MatrixXr BlockPreconditioner::system_solve(const MatrixXr& b) const
{
	if (!initialized)
		Rprintf("Preconditioner not initialized. Using identity");
	return BaseSolver::system_solve(preconditionRHS(b));
}