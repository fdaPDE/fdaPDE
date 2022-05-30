#ifndef FSPAI_SOLVER_WRAPPER_H_
#define FSPAI_SOLVER_WRAPPER_H_


/*!
  This function is the main interface with the imported FSPAI code. The arguments argc and argv have the same meaning of argc and argv arguments in a generic main.cpp function. argc represents the number of inputs passed in argv minus 1. argv is a vector of vectors of characters containing as a first element "FSPAI_Solver_Wrapper". This function is used by FSPAI_Wrapper.h function when non-exact inversion for a sparse symmetric matrix is required.
*/
int FSPAI_Solver_Wrapper (int argc, char ** argv);  

#endif
