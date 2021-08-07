#include"Inference_Carrier.h"

template<typename InputHandler> 
Inference_Carrier<InputHandler>::Inference_Carrier(const InputHandler * Regression_Data_, const MixedFERegressionBase<InputHandler> * model_, const output_Data * out_regression_, const InferenceData * inf_data_, Real lambda_){

  //Setting lambda (optimal)
  setlambda(lambda_);

  //Setting the problem specifications
  setRegData(Regression_Data_);
  setModel(model_);
  setInfData(inf_data_);

  //Setting from Regression_Data_
  setWp(Regression_Data_->getCovariates());
  setN_obs(Regression_Data_->getNumberofObservations());
  setq(Regression_Data_->getCovariates()->cols());
  setZp(Regression_Data_->getObservations());

  //Setting from MixedFERegressionBase
  setN_nodes(model_->getnnodes_());
  setPsip (model_->getpsi_());
  setPsi_tp (model_->getpsi_t_());
  setWtW_decp(model_->getWTW_());
  setR0p(model_->getR0_());
  setR1p(model_->getR1_());
  setHp(model_->getH_());
  setUp(model_->getU_());
  setVp(model_->getV_());
  setEp(model_->getmatrixNoCov_());
  setE_decp(model_->getmatrixNoCovdec_());
  setG_decp(model_->getGdec_());

  //Setting from Output
  setBeta_hatp(&(out_regression_->betas));
  setZ_hatp(&(out_regression_->z_hat));

};
