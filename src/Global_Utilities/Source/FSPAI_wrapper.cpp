#include "../Include/FSPAI_wrapper.h"

bool FSPAI_wrapper(const SpMat & A, SpMat & A_inv){
  
  std::string Temp_name_write = "Tmp_file_mtx_to_inv_FSPAI"; // Name of temporary file for the matrix to be inverted
  std::string Temp_name_read  = "Tmp_file_mtx_invtd_FSPAI";  // Name of temporary file for inverted matrix
  bool saved_Mat = false;
  bool read_Mat  = false;
  
  saved_Mat=Eigen::saveMarket(A, Temp_name_write); // Save the matrix in market format into a temporary file that will be read by FSPAI interface

  if(saved_Mat!=true){
    Rprintf("Internal error: unable to communicate with FSPAI, inference discarded");
    return false;
  }

  // DO ALL THAT YOU NEED WITH FSAPI

  read_Mat=Eigen::loadMarket(A_inv, Temp_name_read); // Read the matrix from a temporary file produced by FSPAI in Market format (May be Preconditioner (inv) or PCG sol)

  if(read_Mat!=true){
    Rprintf("Internal error: unable to communicate with FSPAI, inference discarded");
    return false;
  }

  remove(Temp_file_write.c_str()); // Remove the temporary file
  remove(Temp_file_read.c_str());  // Remove the temporary file
  
  return true;
}

