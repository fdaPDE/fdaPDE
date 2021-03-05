#include"Inference_Carrier.h"

template<InputHandler> Inference_Carrier<InputHandler>::Inference_Carrier(const InputHandler * Regression_Data_, const MixedFERegressionBase<InputHandler> * model_, const output_Data * out_regression, const InferenceData * inf_data_){

//Setting the datasets
setOptData(Regression_Data_);
setModel(model_);
setInfData(inf_data_);

//Setting from Regression_Data_
setWp(&(Regression_Data_->getCovariates()));
setN_obs(Regression_Data_->getNumberofObservations());
setp(Regression_Data_->getCovariates()->cols())
if(std.:dynamic_cast<RegressionDataEllipticSpaceVarying>(Regression_Data_)!=nullptr){ //check space varying
setKp(&(Regression_Data_->getK()));	//Otherwise it is nullptr
}



//Setting from MixedFERegressionBase
setN_nodes(model_->getnnodes_());
setPsip (model_->getpsi_());
setPsi_tp (model_->getpsi_t_());
setR0p (model_->getR0_());
setR1p (model_->getR1_());
setPp(model_->getR_());
setWtW_decp(model_->getWTW_());
setHp(model_->getH_());
setUp(model_->getU_());
setVp(model_->getV_());
setEp(model_->getmatrixNoCov_());
setE_decp(model->getmatrixNoCovdec_());
setG_decp(model_->getGdec_());

//Setting from Output
setLambda(out_regression->lambda_sol);
setBeta_hatp(&(out_regression->betas));
setZp(&(Regression_Data_->observations_));
setZ_hatp(&(out_regression));
setVar_res(out_regression.sigma_hat_sq)

//Last to be executed (needs the other elements to be set in order to work properly)
setF_hatp();
}
