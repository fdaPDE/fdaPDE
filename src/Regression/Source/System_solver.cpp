#include "../Include/System_solver.h"

SpMat SystemSolver::buildSystemMatrix(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE)
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

	return M;
}

MatrixXr SystemSolver::system_solve(SpMat& M, MatrixXr& b)
{
	Eigen::SparseLU<SpMat> Mdec(M);
	return Mdec.solve(b);
}

MatrixXr SystemSolver::system_solve(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE, MatrixXr& b)
{
	SpMat M = this->buildSystemMatrix(NW, SE, SW, NE);
	return (this->system_solve(M, b));
}


// ------- Mass lumping techniques ------

void MassLumping::builddiag(SpMat& M)
{
	std::vector<coeff> coefficients;
	coefficients.reserve(M.outerSize());

	for (UInt k = 0; k < M.outerSize(); ++k)
	{
		Real val = 0.0;
		for (SpMat::InnerIterator it(M, k); it; ++it)
			val += it.value();
		coefficients.push_back(coeff(k, 0, val));
	}

	diag.resize(M.outerSize(), 1);
	diag.setZero();
	diag.setFromTriplets(coefficients.begin(), coefficients.end());
	diag.makeCompressed();

	isdiagInitialized = true;
}

SpMat MassLumping::buildSystemMatrix(SpMat& M)
{
	if (!isdiagInitialized)
		this->builddiag(M);
	return diag;
}

SpMat MassLumping::buildSystemMatrix(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE)
{
	UInt nnodes = NW.outerSize();
	std::vector<coeff> coefficients;
	coefficients.reserve(4*nnodes);

	for (UInt k = 0; k < nnodes; ++k)
	{
		Real val = 0.0;
		for (SpMat::InnerIterator it(NW, k); it; ++it)
			val += it.value();
		for (SpMat::InnerIterator it(NE, k); it; ++it)
			val += it.value();
		coefficients.push_back(coeff(k, 0, val));
	}
	for (UInt k = 0; k < nnodes; ++k)
	{
		Real val = 0.0;
		for (SpMat::InnerIterator it(SW, k); it; ++it)
			val += it.value();
		for (SpMat::InnerIterator it(SE, k); it; ++it)
			val += it.value();
		coefficients.push_back(coeff(k+nnodes, 0, val));
	}

	diag.resize(nnodes, 1);
	diag.setZero();
	diag.setFromTriplets(coefficients.begin(), coefficients.end());
	diag.makeCompressed();

	this->isdiagInitialized = true;

	return diag;
}

MatrixXr MassLumping::system_solve(SpMat& M, MatrixXr& b)
{
	if (!(this->isdiagInitialized))
		this->builddiag(M);
	return this->system_solve(b);
}

MatrixXr MassLumping::system_solve(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE, MatrixXr& b)
{
	if (!(this->isdiagInitialized))
		this->buildSystemMatrix(NW, SE, SW, NE);
	return this->system_solve(b);
}

// ---------- Diagonal preconditioning techniques ----------
SpMat diagonalPreconditioner::preconditionMatrix(SpMat& M)
{
	if (isdiagInitialized)
		return diag.asDiagonal() * M;
	else
		return M;
}

SpMat diagonalPreconditioner::preconditionMatrix(SpMat& M,VectorXr& prec)
{
	compute(prec);
	return this->preconditionMatrix(M);
}

MatrixXr diagonalPreconditioner::system_solve(SpMat& M, MatrixXr& b, bool preconditioned = false)
{
	if (!preconditioned)
	{
		M = this->preconditionMatrix(M);
		if (isdiagInitialized)
			b = diag.asDiagonal() * b;
	}
	return SystemSolver::system_solve(M, b);
}

MatrixXr diagonalPreconditioner::system_solve(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE, bool preconditioned = false)
{
	SpMat M = SystemSolver::buildSystemMatrix(NW, SE, SW, NE);
	if (!preconditioned)
	{
		M = this->preconditionMatrix(M);
		if (isdiagInitialized)
			b = diag.asDiagonal() * b;
	}
	return SystemSolver::system_solve(M, b);
}

void JacobiPreconditioner::compute(SpMat& M) {
	std::vector<coeff> coefficients;
	coefficients.reserve(M.outerSize());

	for (UInt k = 0; k < M.outerSize(); ++k)
	{
		Real val = 0.0;
		for (SpMat::InnerIterator it(M, k); it; ++it)
			val += pow(it.value(), 2);
		coefficients.push_back(coeff(k, 0, sqrt(val)));
	}

	SpMat prec(M.innerSize(), 1);
	prec.setFromTriplets(coefficients.begin(), coefficients.end());
	prec.makeCompressed();

	diag = VectorXr(prec.cwiseInverse());

	isdiagInitialized = true;
}

SpMat blocksPreconditioner::buildSystemMatrix()
{

	std::vector<coeff> tripletAll;
	tripletAll.reserve(NWblock.nonZeros() + NEblock.nonZeros() + SWblock.nonZeros() + SEblock.nonZeros());

	for (UInt k = 0; k < nnodes; ++k)
		tripletAll.push_back(coeff(k, k, 1));

	for (UInt k = 0; k < SEblock.outerSize(); ++k)
		for (SpMat::InnerIterator it(SEblock, k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row() + nnodes, it.col() + nnodes, it.value()));
		}

	MatrixXr NE = NWdec.solve(NEblock);
	for (UInt k = 0; k < NE.cols(); ++k)
		for (UInt j = 0; j < NE.rows(); ++j)
		{
			Real val = NE(j, k);
			if (val != 0)
				tripletAll.push_back(coeff(j, k + nnodes, val));
		}

	MatrixXr SW = SEdec.solve(SWblock);
	for (UInt k = 0; k < SW.cols(); ++k)
		for (UInt j = 0; j < SW.rows(); ++j)
		{
			Real val = SW(j, k);
			if (val != 0)
				tripletAll.push_back(coeff(j + nnodes, k, val));
		}

	// Define, resize, fill and compress
	SpMat M;
	M.resize(2 * nnodes, 2 * nnodes);
	M.setZero();
	M.setFromTriplets(tripletAll.begin(), tripletAll.end());
	M.makeCompressed();

	return M;

}

SpMat blocksPreconditioner::buildSystemMatrix(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE)
{
	NWblock = NW;
	SEblock = SE;
	SWblock = SW;
	NEblock = NE;
	nnodes = NW.outersize();

	if (!decomposed)
	{
		NWdec.compute(NW);
		SEdec.compute(SE);
		decomposed = true;
	}

	return this->buildSystemMatrix();
}

MatrixXr blocksPreconditioner::preconditionRHS(MatrixXr& b)
{
	MatrixXr rhs(b.rows(), b.cols());
	rhs.topRows(nnodes) = NWdec.solve(b.topRows(nnodes));
	rhs.bottomRows(nnodes) = SEdec.solve(b.bottomRows(nnodes));

	return rhs;
}

MatrixXr blocksPreconditioner::system_solve(SpMat& M, MatrixXr& b, bool preconditioned = false)
{
	if (!preconditioned)
	{
		b = this->preconditionRHS(b);
	}
	SystemSolver::system_solve(M, b);
}

MatrixXr blocksPreconditioner::system_solve(SpMat& NW, SpMat& SE, SpMat& SW, SpMat& NE, MatrixXr& b, bool preconditioned = false)
{
	SpMat M = this->buildSystemMatrix(NW, SE, SW, NE);
	b = this->preconditionRHS(b);

	SystemSolver::system_solve(M, b, preconditioned);
}