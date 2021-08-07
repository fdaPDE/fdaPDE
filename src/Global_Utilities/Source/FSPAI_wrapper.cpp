#include "../Include/FSPAI_Wrapper.h"
#include "../../Fspai_Utility/Include/FSPAI_Solver_Wrapper.h"

bool FSPAI_Wrapper(const SpMat & A, SpMat & A_inv){
  
  std::string Temp_name_write = "Tmp_file_mtx_to_inv_FSPAI"; // Name of temporary file for the matrix to be inverted
  std::string Temp_name_read  = "Tmp_file_mtx_invtd_FSPAI";  // Name of temporary file for inverted matrix
  bool saved_Mat = false;
  bool read_Mat  = false;
  
  saved_Mat=Eigen::saveMarket(A, Temp_name_write, Eigen::Symmetric); // Save the matrix in market format into a temporary file that will be read by FSPAI interface

  if(saved_Mat!=true){
    Rprintf("Internal error: unable to communicate with FSPAI, inference discarded");
    return false;
  }

  // FSPAI LIBRARY BEGINS

  FSPAI_data FSPAI_dat;
  FSPAI_dat.out_File = Temp_name_read;
  
  // Vector of parameters that are needed by FSPAI library for the computation of the inverse
  std::vector<std::string> SPAI_Arguments = {"FSPAI_Solver_Wrapper", Temp_name_write.c_str(), "-diag", "1", "-ep", FSPAI_dat.tol_Inverse.c_str(),
    "-ns", FSPAI_dat.max_Step_Col.c_str(), "-mn", FSPAI_dat.max_New_Nz.c_str(), "-out", FSPAI_dat.out_File.c_str(), "-sol", FSPAI_dat.sol.c_str()}; //ONLY SEQUENTIAL UP TO NOW

  // Since FSPAI_Solver_Wrapper is expecting as inputs: (int argc, char *argv[])
  std::vector<char*> SPAI_Argv;
  for (const auto& arg : SPAI_Arguments){
    SPAI_Argv.push_back((char*)arg.data());
  }
  SPAI_Argv.push_back(nullptr);
  
  int FSPAI_Inverted = FSPAI_Solver_Wrapper(SPAI_Argv.size() - 1, SPAI_Argv.data()); // Run the FSPAI library solver, computes apprixmate inverse of an SPD matrix 

  // FSPAI LIBRARY ENDS

  SpMat Chol_A_inv;

  read_Mat=Eigen::loadMarket(Chol_A_inv, Temp_name_read); // Read the matrix from a temporary file produced by FSPAI in market format

  if(FSPAI_Inverted!=0 || read_Mat!=true){
    Rprintf("Internal error: unable to communicate with FSPAI correctly, inference discarded");
    return false;
  }

  Chol_A_inv.makeCompressed();

  //Recover matrix A_inv via Cholesky factorization
  A_inv = Chol_A_inv * Chol_A_inv.transpose();
  A_inv.makeCompressed();


  remove(Temp_name_write.c_str()); // Remove the temporary file
  remove(Temp_name_read.c_str());  // Remove the temporary file
  
  return true;
}

