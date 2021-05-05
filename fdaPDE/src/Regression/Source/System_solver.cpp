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


//// ---------- Jacobi diagonal preconditioner methods ----------
//void JacobiPreconditioner::compute(const SpMat& M)
//{
//	prec.resize(M.outerSize());
//
//	for (UInt k = 0; k < M.outerSize(); ++k)
//	{
//		Real val = 0.0;
//		for (SpMat::InnerIterator it(M, k); it; ++it)
//			val += pow(it.value(), 2);
//		prec(k) = 1 / val;
//	}
//	initialized = true;
//	system_factorize(prec.asDiagonal() * M);
//}
//
//MatrixXr JacobiPreconditioner::system_solve(const SpMat& M, const MatrixXr& b)
//{
//	compute(M);
//	return BaseDiagPreconditioner::system_solve(b);
//}

// ---------- Lambda preconditioner methods ----------
void LambdaPreconditioner::compute(const SpMat& M)
{
	if (initialized)
		system_factorize(prec.asDiagonal() * M * prec.asDiagonal());
	else
	{
		Rprintf("Preconditioner not initialized. Using identity");
		system_factorize(M);
	}
}

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
	if (!initialized)
	{
		Rprintf("Error: system size not specified. Unable to construct preconditioner. Using identity");
		return;
	}

	if (lambda == lambda_)
		return;

	lambda = lambda_;
	UInt nnodes = prec.cols() / 2;
	prec.bottomRows(nnodes) = VectorXr::Ones(nnodes) / sqrt(lambda_);
}

MatrixXr LambdaPreconditioner::system_solve(const SpMat& M, const MatrixXr& b)
{
	compute(M);
	if(initialized)
		return prec.asDiagonal() * BaseDiagPreconditioner::system_solve(b);
	return BaseDiagPreconditioner::system_solve(b);
}


// ---------- Abstract blocks digonal preconditioner methods ----------
BaseBlocksPreconditioner::BaseBlocksPreconditioner()
{
	NWblock.resize(1,1);
	SEblock.resize(1,1);
	NWblock.setIdentity();
	SEblock.setIdentity();
	NWdec.compute(NWblock);
	SEdec.compute(SEblock);
}

SpMat BaseBlocksPreconditioner::preconditioner() const
{
	SpMat M;
	if (!initialized)
	{
		Rprintf("Preconditioner not initalized. Returning identity");
		M.resize(1, 1);
		M.setIdentity();
	}
	else
	{
		UInt nnodes = NWblock.outerSize();
		M.resize(2 * nnodes, 2 * nnodes);
		MatrixXr Mtemp(2 * nnodes, 2 * nnodes);
		Mtemp.topLeftCorner(nnodes, nnodes) = NWdec.solve(MatrixXr::Identity(nnodes, nnodes));
		Mtemp.topRightCorner(nnodes, nnodes) = SEdec.solve(MatrixXr::Identity(nnodes, nnodes));
		M = Mtemp.sparseView();
	}
	return M;
}

MatrixXr BaseBlocksPreconditioner::preconditionRHS(const MatrixXr& b) const
{
	if (!initialized)
	{
		Rprintf("Preconditioner not initialized. Using identity");
		return b;
	}
	
	UInt nnodes = NWblock.outerSize();
	MatrixXr rhs(b.rows(), b.cols());
	rhs.topRows(nnodes) = NWdec.solve(b.topRows(nnodes));
	rhs.bottomRows(nnodes) = SEdec.solve(b.bottomRows(nnodes));
	return rhs;
}

MatrixXr BaseBlocksPreconditioner::system_solve(const SpMat& M, const MatrixXr& b)
{
	if (initialized)
		return BaseSolver::system_solve(preconditioner() * M, preconditionRHS(b));

	Rprintf("Preconditioner not initialized. Using identity");
	return BaseSolver::system_solve(M, b);
}

MatrixXr BaseBlocksPreconditioner::system_solve(const MatrixXr& b) const
{
	if (!initialized)
		Rprintf("Preconditioner not initialized. Using identity");
	return BaseSolver::system_solve(preconditionRHS(b));
}

void BaseBlocksPreconditioner::compute(const SpMat& M)
{
	if (initialized)
	{
		UInt nnodes = M.rows() / 2;
		system_factorize(buildSystemMatrix(NWdec.solve(M.topLeftCorner(nnodes, nnodes)), 
			SEdec.solve(M.bottomRightCorner(nnodes, nnodes)), SEdec.solve(M.bottomLeftCorner(nnodes, nnodes)),
			NWdec.solve(M.topRightCorner(nnodes, nnodes))));
	}
	else
	{
		Rprintf("Preconditioner not initialized. Using identity");
		system_factorize(M);
	}
}


// ---------- Blocks digonal preconditioner methods ----------
void BlocksPreconditioner::compute(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
{
	BaseSolver::system_factorize(this->buildSystemMatrix(NW, SE, SW, NE));
}

SpMat BlocksPreconditioner::buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
{
	NWblock = NW;
	SEblock = SE;

	NWdec.compute(NWblock);
	SEdec.compute(SEblock);
	initialized = true;

	UInt nnodes = NW.outerSize();
	SpMat Id;
	Id.resize(nnodes, nnodes);
	Id.setIdentity();

	return BaseSolver::buildSystemMatrix(Id, Id, SEdec.solve(SW), NWdec.solve(NE));
}

void BlocksPreconditioner::compute(const SpMat& M)
{
	UInt nnodes = M.outerSize() / 2;

	compute(M.topLeftCorner(nnodes, nnodes), M.bottomRightCorner(nnodes, nnodes),
		M.bottomLeftCorner(nnodes, nnodes), M.topRightCorner(nnodes, nnodes));
}

MatrixXr BlocksPreconditioner::system_solve(const SpMat& M, const MatrixXr& b)
{
	compute(M);
	return BaseBlocksPreconditioner::system_solve(b);
}


// --------- Lambda Block Preconditioner methods ---------

//MatrixXr BlockLambdaPreconditioner::system_solve(const SpMat& M, const MatrixXr& b)
//{
//	compute(M);
//	return system_solve(b);
//}
//
//MatrixXr BlockLambdaPreconditioner::system_solve(const MatrixXr& b) const
//{
//	return L.preconditioner().asDiagonal() * BaseBlocksPreconditioner::system_solve(L.preconditionRHS(b));
//}
//
//void BlockLambdaPreconditioner::compute(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
//{
//	system_factorize(buildSystemMatrix(NW, SE, SW, NE));
//}
//
//void BlockLambdaPreconditioner::compute(const SpMat& M)
//{
//	if (L.getlambda() != 1)
//		BlocksPreconditioner::compute(L.preconditioner().asDiagonal() * M * L.preconditioner().asDiagonal());
//	else
//		BlocksPreconditioner::compute(M);
//}
//
//void BlockLambdaPreconditioner::compute(const Real lambda_)
//{
//	if (lambda_ != L.getlambda())
//		L.compute(lambda_);
//}