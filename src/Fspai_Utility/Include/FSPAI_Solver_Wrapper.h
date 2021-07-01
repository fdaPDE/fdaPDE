#ifndef FSPAI_SOLVER_WRAPPER_H_
#define FSPAI_SOLVER_WRAPPER_H_

//C++/MPI includings
#include <iostream>
#include <stdexcept>

// file includings
#include "ENV_Handler.h"
#include "Macros.h"
#include "Pe_Exception.h"
#include "Command_Parser.h"
#include "Timer.h"
#include "Type_Base_Handler.h"
#include "Double_Handler.h"
#include "Complex_Handler.h"
#include "MMio.h"
#include "Pattern.h"

int FSPAI_Solver_Wrapper (int argc, char ** argv);  

#endif
