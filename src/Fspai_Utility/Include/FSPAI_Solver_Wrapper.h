////////////////////////////////////////////////////////////
/// \mainpage FSPAI (Factorized Sparse Approximate Inverses)
/// \section section1 ABOUT
/// FSPAI is the implementation of the known FSPAI algorithm
/// introduced by T. Huckle. \n It generates a preconditioner for large sparse
/// and ill-conditioned symmetric positive \n definite systems
/// of linear equations. FSPAI is inherently parallel and the  computed \n
/// preconditioner approximates the inverse of the Cholesky
/// factor of the system matrix. \n It is the factorized version of
/// the SPAI algorithm. The \e FSPAI was implemented at the \n research unit \n\n
/// \e Informatik \e V -- \e Scientific \e Computing \e in \e Computer \e Science \n
/// \e Technische \e Universität \e München. \n \n \n
/// \b DESIGNED \b BY \n \n
/// Matous Sedlacek <sedlacek@in.tum.de> \n \n \n
/// \b RELEASED \b 2011 \n \n
/// FSPAI is published under the LGPL in year 2011. \n \n \n
/// \b LICENSE \n \n
/// FSPAI: Factorized Sparse Approximate Inverses \n
/// Copyright © 2011 Matous Sedlacek \n
/// Scientific Computing in Computer Science -- Informatics V \n
/// Technische Universität München \n \n
/// FSPAI is free software: you can redistribute it and/or modify it under the \n
/// terms of the GNU Lesser General Public License as published by the Free Software \n
/// Foundation, either version 3 of the License, or (at your option) any later version. \n
/// \n
/// FSPAI is distributed in the hope that it will be useful, but WITHOUT ANY \n
/// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR \n
/// A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details. \n
/// \n
/// You should have received a copy of the GNU Lesser General Public License along with \n
/// FSPAI. If not, see http://www.gnu.org/licenses/. \n
/// \n
/// If you obtain any results with FSPAI we would appreciate that you refer to FSPAI. \n
/////////////////////////////////////////////////////////////

#ifndef FSPAI_SOLVER_WRAPPER_H_
#define FSPAI_SOLVER_WRAPPER_H_


int FSPAI_Solver_Wrapper (int argc, char ** argv);  



#endif
