#include"Inference_Carrier.h"

template<InputHandler> Inference_Carrier<InputHandler>::Inference_Carrier(const InputHandler * Regression_Data_, const MixedFERegressionBase<InputHandler> * model_, const output_Data * out_regression, const InferenceData * inf_data_){

//Setting the datasets
setOptData(Regression_Data_);
setModel(model_);
setInfData(inf_data_);

//Setting from Regression_Data_
setWp(&(Regression_Data_->covariates_));
if(std.:dynamic_cast<RegressionDataEllipticSpaceVarying>(Regression_Data_)!=nullptr){ //check space varying
setKp(&(Regression_Data_->K_));	//Otherwise it is nullptr


//Setting from MixedFERegressionBase
setPsip (&(model_->psi_));
setPsi_tp (&(model_->psi_t_));
setR0p (&(model_->R0_lambda));
setR1p (&(model_->R1_lambda));
setPp(&(model_->R_));
setWtW_decp(&(model_->WTW_));
setHp(&(model_->H_));
setUp(&(model_->U_));
setVp(&(model_->V_));
setEp(&(model_->matrixNoCov_));
setE_decp(&(model->matrixNoCovdec_));
setG_decp(&(model_->Gdec_));

//Setting from Output
setLambda(out_regression->lambda_sol);
setBeta_hatp(&(out_regression->betas));
setZp(&(Regression_Data_->observations_));
setZ_hatp(&(out_regression));

//Last to be executed (needs the other elements to be set in order to work properly)
setF_hatp();
}
