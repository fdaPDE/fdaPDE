#include "../Include/System_solver.h"


SpMat SystemSolver::buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
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

MatrixXr SystemSolver::system_solve(const SpMat& M, const MatrixXr& b)
{
	//Eigen::SparseLU<SpMat> Mdec(M);
	Mdec.compute(M);
	return Mdec.solve(b);
}

MatrixXr SystemSolver::system_solve(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE, const MatrixXr& b)
{
	SpMat M = this->buildSystemMatrix(NW, SE, SW, NE);
	return (this->system_solve(M, b));
}

MatrixXr SystemSolver::system_solve(const MatrixXr& b)
{
	if (!decomposed)
	{
		Rprintf("System matrix not initialized");
		return b;
	}

	else
		return Mdec.solve(b);
}


// ------- Mass lumping methods ------

void MassLumping::builddiag(const SpMat& M)
{
	diag.resize(M.outerSize());
	//std::vector<coeff> coefficients;
	//coefficients.reserve(M.outerSize());

	for (UInt k = 0; k < M.outerSize(); ++k)
	{
		Real val = 0.0;
		for (SpMat::InnerIterator it(M, k); it; ++it)
			val += it.value();
		//coefficients.push_back(coeff(k, 0, val));
		diag(k) = val;
	}

	/*diag.resize(M.outerSize(), 1);
	diag.setZero();
	diag.setFromTriplets(coefficients.begin(), coefficients.end());
	diag.makeCompressed();*/

	isdiagInitialized = true;
}

SpMat MassLumping::buildSystemMatrix(const SpMat& M)
{
	this->builddiag(M);
	return SpMat(diag.asDiagonal());
}

SpMat MassLumping::buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
{
	UInt nnodes = NW.outerSize();
	/*std::vector<coeff> coefficients;
	coefficients.reserve(4*nnodes);*/
	diag.resize(2 * nnodes);

	for (UInt k = 0; k < nnodes; ++k)
	{
		Real val = 0.0;
		for (SpMat::InnerIterator it(NW, k); it; ++it)
			val += it.value();
		for (SpMat::InnerIterator it(NE, k); it; ++it)
			val += it.value();
		//coefficients.push_back(coeff(k, 0, val));
		diag(k) = val;
	}
	for (UInt k = 0; k < nnodes; ++k)
	{
		Real val = 0.0;
		for (SpMat::InnerIterator it(SW, k); it; ++it)
			val += it.value();
		for (SpMat::InnerIterator it(SE, k); it; ++it)
			val += it.value();
		//coefficients.push_back(coeff(k+nnodes, 0, val));
		diag(k + nnodes) = val;
	}

	/*diag.resize(nnodes, 1);
	diag.setZero();
	diag.setFromTriplets(coefficients.begin(), coefficients.end());
	diag.makeCompressed();*/

	this->isdiagInitialized = true;

	return SpMat(diag.asDiagonal());
}

MatrixXr MassLumping::system_solve(const MatrixXr& b)
{
	if (!isdiagInitialized)
	{
		Rprintf("System matrix not initialized. Solving for Identity operator");
		return b;
	}
	return diag.cwiseInverse().asDiagonal() * b;
}

MatrixXr MassLumping::system_solve(const SpMat& M, const MatrixXr& b)
{
	this->builddiag(M);
	return this->system_solve(b);
}

MatrixXr MassLumping::system_solve(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE, const MatrixXr& b)
{
	this->buildSystemMatrix(NW, SE, SW, NE);
	return this->system_solve(b);
}

// ---------- Diagonal preconditioner methods ----------
SpMat diagonalPreconditioner::preconditionMatrix(const SpMat& M)
{
	if (isdiagInitialized)
		return diag.asDiagonal() * M;
	else
		return M;
}

SpMat diagonalPreconditioner::preconditionMatrix(const SpMat& M, const VectorXr& prec)
{
	preconditioner(prec);
	return preconditionMatrix(M);
}

MatrixXr diagonalPreconditioner::preconditionRHS(const MatrixXr& b)
{
	if (isdiagInitialized)
		return diag.asDiagonal() * b;

	return b;
}

MatrixXr diagonalPreconditioner::system_solve(const SpMat& M, const MatrixXr& b, bool preconditioned = false)
{
	SpMat Mp(M);
	MatrixXr bp(b);
	if (!preconditioned)
	{
		Mp = preconditionMatrix(M);
		if (isdiagInitialized)
			bp = diag.asDiagonal() * b;
		else
			Rprintf("Using identity preconditioner");
	}
	return SystemSolver::system_solve(Mp, bp);
}

MatrixXr diagonalPreconditioner::system_solve(const SpMat & NW, const SpMat & SE, const SpMat & SW, const SpMat & NE, const MatrixXr & b, bool preconditioned = false)
{
	SpMat Mp = SystemSolver::buildSystemMatrix(NW, SE, SW, NE);
	MatrixXr bp(b);
	if (!preconditioned)
	{
		Mp = preconditionMatrix(Mp);
		bp = preconditionRHS(b);
		/*if (isdiagInitialized)
			b = diag.asDiagonal() * b;*/
	}
	return SystemSolver::system_solve(Mp, bp);
}

MatrixXr diagonalPreconditioner::system_solve(const MatrixXr& b, bool preconditioned = false)
{
	if (!preconditioned)
		return SystemSolver::system_solve(preconditionRHS(b));

	return SystemSolver::system_solve(b);
}

SpMat diagonalPreconditioner::buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
{
	SpMat M = SystemSolver::buildSystemMatrix(NW, SE, SW, NE);
	return preconditionMatrix(M);
}

SpMat diagonalPreconditioner::buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE, const VectorXr& prec)
{
	SpMat M = SystemSolver::buildSystemMatrix(NW, SE, SW, NE);
	return preconditionMatrix(M, prec);
}

// ---------- Jacobi diagonal preconditioner methods ----------
void JacobiPreconditioner::compute(const SpMat& M)
{
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

// ---------- Blocks digonal preconditioner methods ----------
blocksPreconditioner::blocksPreconditioner()
{
	NWblock.resize(1, 1);
	SEblock.resize(1, 1);
	SWblock.resize(1, 1);
	NEblock.resize(1, 1);

	std::vector<coeff> zerovec, onevec;
	zerovec.push_back(coeff(0, 0, 0));
	onevec.push_back(coeff(0, 0, 1));
	NWblock.setFromTriplets(zerovec.begin(), zerovec.end());
	SEblock.setFromTriplets(zerovec.begin(), zerovec.end());
	SWblock.setFromTriplets(onevec.begin(), onevec.end());
	NEblock.setFromTriplets(onevec.begin(), onevec.end());

	nnodes = 1;

	NWdec.compute(NWblock);
	SEdec.compute(SEblock);
}

SpMat blocksPreconditioner::buildSystemMatrix()
{
	if (!blocksDecomposed)
	{
		Rprintf("Blocks not initialized. Return Identity");
		VectorXr onevec(2);
		onevec << 1, 1;
		return SpMat(onevec.asDiagonal());
	}

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

	Mdec.compute(M);
	decomposed = true;

	return M;
}

void blocksPreconditioner::compute(const SpMat& M)
{
	nnodes = M.outerSize() / 2;

	std::vector<coeff> NEtrip, SWtrip, SEtrip, NWtrip;
	NEtrip.reserve(pow(nnodes, 2));
	SWtrip.reserve(pow(nnodes, 2));
	SEtrip.reserve(pow(nnodes, 2));
	NWtrip.reserve(pow(nnodes, 2));

	for (UInt k = 0; k < 2 * nnodes; ++k)
	{
		for (SpMat::InnerIterator it(M, k); it; ++it)
		{
			if (it.col() < nnodes)
			{
				if (it.row() < nnodes)
					NEtrip.push_back(coeff(it.row(), it.col(), it.value()));
				else
					SEtrip.push_back(coeff(it.row() - nnodes, it.col(), it.value()));
			}
			else
			{
				if (it.row() < nnodes)
					NWtrip.push_back(coeff(it.row(), it.col() - nnodes, it.value()));
				else
					SWtrip.push_back(coeff(it.row() - nnodes, it.col() - nnodes, it.value()));
			}
		}
	}

	NEblock.resize(nnodes, nnodes);
	NEblock.setZero();
	NEblock.setFromTriplets(NEtrip.begin(), NEtrip.end());
	NEblock.makeCompressed();

	SWblock.resize(nnodes, nnodes);
	SWblock.setZero();
	SWblock.setFromTriplets(SWtrip.begin(), SWtrip.end());
	SWblock.makeCompressed();

	SEblock.resize(nnodes, nnodes);
	SEblock.setZero();
	SEblock.setFromTriplets(SEtrip.begin(), SEtrip.end());
	SEblock.makeCompressed();

	NWblock.resize(nnodes, nnodes);
	NWblock.setZero();
	NWblock.setFromTriplets(NWtrip.begin(), NWtrip.end());
	NWblock.makeCompressed();

	NWdec.compute(NWblock);
	SEdec.compute(SEblock);

	blocksDecomposed = true;

	buildSystemMatrix();
}

SpMat blocksPreconditioner::buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
{
	NWblock = NW;
	SEblock = SE;
	SWblock = SW;
	NEblock = NE;
	nnodes = NW.outerSize();

	/*if (!decomposed)
	{*/
		NWdec.compute(NW);
		SEdec.compute(SE);
		blocksDecomposed = true;
	/*}*/

	return this->buildSystemMatrix();
}

MatrixXr blocksPreconditioner::preconditionRHS(const MatrixXr& b)
{
	if (!blocksDecomposed)
	{
		Rprintf("Blocks not initialized. Using identity preconditioner");
		return b;
	}

	MatrixXr rhs(b.rows(), b.cols());
	rhs.topRows(nnodes) = NWdec.solve(b.topRows(nnodes));
	rhs.bottomRows(nnodes) = SEdec.solve(b.bottomRows(nnodes));

	return rhs;
}

MatrixXr blocksPreconditioner::system_solve(const MatrixXr& b, bool preconditioned = false)
{
	if (!blocksDecomposed)
	{
		Rprintf("Blocks not initialized. Solving for the identity operator");
		return b;
	}

	MatrixXr bp(b);

	if(!decomposed)
		SpMat Mp = buildSystemMatrix();
	if (!preconditioned)
		bp = preconditionRHS(b);

	return Mdec.solve(bp);
}

MatrixXr blocksPreconditioner::system_solve(const SpMat& M, const MatrixXr& b, bool preconditioned = false)
{
	/*if (!decomposed)
		M = buildSystemMatrix(M);*/

	MatrixXr bp(b);

	if (!preconditioned)
		bp = preconditionRHS(b);

	return SystemSolver::system_solve(M, bp);
}

MatrixXr blocksPreconditioner::system_solve(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE, const MatrixXr& b, bool preconditioned = false)
{
	SpMat M = buildSystemMatrix(NW, SE, SW, NE);

	return system_solve(M, b, preconditioned);
}