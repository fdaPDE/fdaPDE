#include "../Include/System_solver.h"

// ---------- Base Solver methods ----------

SpMat BaseSolver::assembleMatrix(const SpMat& DMat, const SpMat& R0, const SpMat& R1, Real lambdaS)
{
	//Rprintf("Using base solver");
	SpMat R0_lambda = (-lambdaS) * R0; 
	SpMat R1_lambda = (-lambdaS) * R1;

	SpMat R1_lambdaT(R1_lambda.transpose());

	//// Add time-dependent blocks
	//if(timeDependent && parabolic && !iterative)
	//	R1_lambda -= lambdaS * (lambdaT * (*LR0k));
	//if(timeDependent && parabolic && iterative)
	//	R1_lambda = lambdaS * R1 - lambdaT * R0_lambda;

	//if (timeDependent && !parabolic)
	//	return buildSystemMatrix(DMat + lambdaT * (*Ptk), R0_lambda, R1_lambda, R1_lambdaT);
	
	// Build the system matrix from the four blocks
	return buildSystemMatrix(DMat, R0_lambda, R1_lambda, R1_lambdaT);
}

SpMat BaseSolver::buildSystemMatrix(const SpMat& NW, const SpMat& SE, const SpMat& SW, const SpMat& NE)
{

	UInt nnodes = NW.outerSize();
	//	nnodes *= M_;

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

//void BaseSolver::addTimeCorrection(const std::shared_ptr<SpMat> Mat, Real lambda, bool flagParabolic)
//	{
//		parabolic = flagParabolic;
//		if (parabolic)
//			LR0k = Mat;
//		else
//			Ptk = Mat;
//		lambdaT = lambda;
//		timeDependent = true;
//	}


// ---------- Mass lumping methods ----------

SpMat MassLumping::lumpMassMatrix(const SpMat& M)
{
	VectorXr diag(M.outerSize());

	for (UInt k = 0; k < M.outerSize(); ++k)
	{
		Real val = 0.0;
		for (SpMat::InnerIterator it(M, k); it; ++it)
			val += it.value();
		diag(k) = val;
	}

	SpMat lumpedMass(diag.asDiagonal());
	return lumpedMass;
}

SpMat MassLumping::assembleMatrix(const SpMat& DMat, const SpMat& R0, const SpMat& R1, Real lambdaS)
{
	//Rprintf("Using mass lumping");
	SpMat R0_lambda = (-lambdaS) * lumpMassMatrix(R0); 
	SpMat R1_lambda = (-lambdaS) * R1;

	SpMat R1_lambdaT(R1_lambda.transpose());

	// Add time-dependent blocks
	/*if (timeDependent && !parabolic)
		return buildSystemMatrix(DMat + lambdaT * (*Ptk), R0_lambda, R1_lambda, R1_lambdaT);
	else if (parabolic)
		R1_lambda -= lambdaS * (lambdaT * lumpMassMatrix(*LR0k));*/

	// Build the system matrix from the four blocks
	return buildSystemMatrix(DMat, R0_lambda, R1_lambda, R1_lambdaT);
}

// ---------- Base diagonal preconditioner methods ----------

MatrixXr BaseDiagPreconditioner::preconditionRHS(const MatrixXr& b) const
{
	if (initialized)
		return prec.asDiagonal() * b;

	Rprintf("Preconditioner not initialized. Using identity");
	return b;
}

// ---------- Lambda preconditioner methods ----------

SpMat LambdaPreconditioner::assembleMatrix(const SpMat& DMat, const SpMat& R0, const SpMat& R1, Real lambdaS)
{
	//Rprintf("Using lambda preconditioner");
	SpMat R0_lambda = -R0; // build the SouthEast block of the matrix
	SpMat R1_lambda = (-sqrt(lambdaS)) * R1;

	// Store the space smoothing parameter and the diagonal of the preconditioner
	lambda = lambdaS;
	UInt nnodes = R0.outerSize();
	prec.resize(2*nnodes);
	prec = VectorXr::Ones(2*nnodes);
	prec.bottomRows(nnodes) = prec.bottomRows(nnodes) / sqrt(lambdaS);
	initialized = true;

	SpMat R1_lambdaT(R1_lambda.transpose());

	//// Add time-dependent blocks
	//if (timeDependent && !parabolic)
	//	return buildSystemMatrix(DMat + lambdaT * (*Ptk), R0_lambda, R1_lambda, R1_lambdaT);
	//else if (parabolic)
	//	R1_lambda -= sqrt(lambdaS) * (lambdaT * (*LR0k));

	// Build the system matrix from the four blocks
	return buildSystemMatrix(DMat, R0_lambda, R1_lambda, R1_lambdaT);
}

MatrixXr LambdaPreconditioner::system_solve(const MatrixXr& b) const
{
	if (initialized)
	{
		MatrixXr sol = BaseDiagPreconditioner::system_solve(preconditionRHS(b));
		UInt nnodes = b.rows() / 2;
		sol.bottomRows(nnodes) = sol.bottomRows(nnodes) / sqrt(lambda);
		return sol;
	}
	return BaseDiagPreconditioner::system_solve(b);
}

MatrixXr LambdaPreconditioner::preconditionRHS(const MatrixXr& b) const
{
	UInt nnodes = b.rows() / 2;
	MatrixXr rhs(b);
	if (initialized && rhs.bottomRows(nnodes) != MatrixXr::Zero(nnodes, b.cols()))
		rhs.bottomRows(nnodes) = rhs.bottomRows(nnodes) / sqrt(lambda);

	return rhs;
}

// ---------- Block digonal preconditioner methods ----------

SpMat BlockPreconditioner::assembleMatrix(const SpMat& DMat, const SpMat& R0, const SpMat& R1, Real lambdaS)
{
	//Rprintf("Using block preconditioner");
	SpMat R0_lambda = (-lambdaS) * R0;
	SpMat R1_lambda = -R1;

	//// Add time-dependent blocks in the parabolic case
	//if (timeDependent && parabolic)
	//	R1_lambda -= (lambdaT * (*LR0k));

	UInt nnodes = R0_lambda.outerSize();
	
	SEblock = R0_lambda;
	SEdec.compute(R0);
	SWblock = SEdec.solve(R1_lambda);
	initialized = true;

	R1_lambda = lambdaS * R1_lambda;

	SpMat Id(nnodes, nnodes);
	Id.setIdentity();

	SpMat R1_lambdaT(R1_lambda.transpose());

	//// Add time-dependent blocks in the separable case
	//// And build the system matrix from the four blocks
	//if (timeDependent && !parabolic)
	//	return buildSystemMatrix(DMat + lambdaT * (*Ptk), Id, SWblock, R1_lambdaT);

	// Build the system matrix from the four blocks
	return buildSystemMatrix(DMat, Id, SWblock, R1_lambdaT);
}

BlockPreconditioner::BlockPreconditioner()
{
	SEblock.resize(1,1);
	SEblock.setIdentity();
	SEdec.compute(SEblock);
	lambda = 1.0;
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
	if (rhs.bottomRows(nnodes) != MatrixXr::Zero(nnodes, b.cols()))
		rhs.bottomRows(nnodes) = (SEdec.solve(b.bottomRows(nnodes)))/lambda;
	return rhs;
}

MatrixXr BlockPreconditioner::system_solve(const MatrixXr& b) const
{
	if (!initialized)
		Rprintf("Preconditioner not initialized. Using identity");
	MatrixXr res = BaseSolver::system_solve(preconditionRHS(b));
	return res;
}

//void BlockPreconditioner::addTimeCorrection(const std::shared_ptr<SpMat> Mat, Real lambda, bool flagParabolic)
//	{
//		// If the time constant does not change, no need for recomputing the factorization and blocks
//		if (lambda == lambdaT)
//			recompute = false;
//		else
//			recompute = true;
//		//BaseSolver::addTimeCorrection(Mat, lambda, flagParabolic);
//	}
