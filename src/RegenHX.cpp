#include "RegenHX.h"

RegenHX::RegenHX()
{
	regenModel = new RegeneratorModel();
}

RegenHX::~RegenHX()
{
	delete regenModel;
}

int RegenHX::getDesignSolution()
{
	RegeneratorSolution solution;
	int status = regenModel->getDesignSolution();
	if (status < 0) {
		return status;
	}

	regenModel->getSolution(&solution);

	costHX = solution.costModule;
	T_H_out = solution.T_H_out;
	T_C_out = solution.T_C_out;
	dP_H = solution.dP_H;
	dP_C = solution.dP_C;
	epsilon = solution.epsilon;
	NTU_R_e = solution.NTU_R_e;
	Q_dot_a = solution.Q_dot_a;
	UA = solution.UA;
	L = solution.L;
	D_fr = solution.D_fr;
	AspectRatio = L / D_fr;
	wallThickness = solution.wallThickness;
	m_dot_carryover = solution.m_dot_carryover;
	m_HTR_LP_dP = solution.m_HTR_LP_dP;
	m_HTR_HP_dP = solution.m_HTR_HP_dP;

	costHX *= modulesInParallel;
	if (operationMode == valveDesignOption::VDO1) {
		UA *= modulesInParallel;
		Q_dot_a *= modulesInParallel;
		m_dot_carryover *= modulesInParallel;
	}


	return status;
}

void RegenHX::set_params(int target_1, int target_2, int operation_mode, int valveMode, double target_2_value, double P_0, double D_s, double e_v, double Q_dot_loss)
{
	this->target_1 = static_cast<targetModes::targetModes>(target_1) ;
	this->target_2 = static_cast<targetModes::target2Modes>(target_2);
	this->valveMode = static_cast<valveDesignOption::valveModes>(valveMode);
	this->target_2_value = target_2_value;
	this->operationMode = static_cast<valveDesignOption::valveDesignOption>(operation_mode);
	this->P_0 = P_0;
	this->D_s = D_s;
	this->Q_dot_loss = Q_dot_loss;
	this->e_v = e_v;
}

int RegenHX::getDesignSolution(double* results)
{
	int status = getDesignSolution();

	if (status < 0) {
		results[0] = -1;
		results[1] = -1;
		results[2] = -1;
		results[3] = -1;
		results[4] = -1;
		results[5] = -1;
		results[6] = -1;
		results[7] = -1;
		results[8] = -1;
		results[9] = -1;
		results[10] = -1;
		results[11] = -1;

		return status;
	}

	results[0] = epsilon;
	results[1] = costHX;
	results[2] = UA;
	results[3] = T_H_out;
	results[4] = dP_H;
	results[5] = dP_C;
	results[6] = D_fr;
	results[7] = L;
	results[8] = wallThickness;
	results[9] = m_dot_carryover;
	results[10] = m_HTR_LP_dP;
	results[11] = m_HTR_HP_dP;

	return status;
}

void RegenHX::initialize(double N_sub_hx)
{
}

double RegenHX::od_delta_p_cold(double m_dot_c /*kg/s*/)
{
	return ms_des_solved.m_DP_cold_des*pow(m_dot_c / m_dot_C, 1.75);
}

double RegenHX::od_delta_p_hot(double m_dot_h /*kg/s*/)
{
	return ms_des_solved.m_DP_hot_des*pow(m_dot_h / m_dot_H, 1.75);
}

void RegenHX::design_fix_UA_calc_outlet(double UA_target, double eff_limit, double T_c_in, double P_c_in, double m_dot_c, double P_c_out, double T_h_in, double P_h_in, double m_dot_h, double P_h_out, double & q_dot, double & T_c_out, double & T_h_out)
{
	design_fix_TARGET_calc_outlet(0, UA_target, eff_limit, T_c_in, P_c_in, m_dot_c, P_c_out,
		T_h_in, P_h_in, m_dot_h, P_h_out, q_dot, T_c_out, T_h_out);
}

void RegenHX::design_fix_TARGET_calc_outlet(int targetType /*-*/, double targetValue /*kW/K or $*/, double eff_limit /*-*/, double T_c_in /*K*/, double P_c_in /*kPa*/, double m_dot_c /*kg/s*/, double P_c_out /*kPa*/,
	double T_h_in /*K*/, double P_h_in /*kPa*/, double m_dot_h /*kg/s*/, double P_h_out /*kPa*/,
	double & q_dot /*kWt*/, double & T_c_out /*K*/, double & T_h_out /*K*/)
{
	resetDesignStructure();
	if (targetValue < 1.E-10) {
		q_dot = 0;
		T_c_out = T_c_in;
		T_h_out = T_h_in;

		resetDesignStructure();
		ms_des_solved.m_m_dot_carryover = 0;
		ms_des_solved.m_DP_cold_des = 0;
		ms_des_solved.m_DP_hot_des = 0;
		return;
	}

	targetModes::target2Modes secondTargetMode = target_2;
	double dP_max_Regen = target_2_value;
	

	setParameters(operationMode, valveMode, Q_dot_loss, P_0, D_s, e_v);
	setInletStates(T_h_in, P_h_in, m_dot_h, T_c_in, P_c_in, m_dot_c);


	setDesignTargets(targetModes::UA, secondTargetMode, targetValue, dP_max_Regen);

	int status = getDesignSolution();

	if (status < 0) {
		resetDesignStructure();
		throw("RegenHX::design",
			"Regenerator model failed!");
	}

	q_dot = Q_dot_a;
	T_c_out = T_C_out;
	T_h_out = T_H_out;
}

void RegenHX::resetDesignStructure()
{
	ms_des_solved.m_DP_cold_des = ms_des_solved.m_DP_hot_des = ms_des_solved.m_eff_design = ms_des_solved.m_min_DT_design = ms_des_solved.m_NTU_design =
		ms_des_solved.m_Q_dot_design = ms_des_solved.m_T_c_out = ms_des_solved.m_T_h_out = ms_des_solved.m_UA_design_total = 
		ms_des_solved.m_aUA_design_total = ms_des_solved.m_cost_design_total =
		ms_des_solved.m_m_dot_carryover = std::numeric_limits<double>::quiet_NaN();
	ms_des_solved.m_eff_limited = false;
}

void RegenHX::setInletStates(double T_H_in, double P_H_in, double m_dot_H, double T_C_in, double P_C_in, double m_dot_C)
{
	this->T_C_in = T_C_in;
	this->T_H_in = T_H_in;
	this->m_dot_H = m_dot_H;
	this->m_dot_C = m_dot_C;

	if (this->operationMode == valveDesignOption::VDO1) {
		regenModel->setInletStates(T_H_in, P_H_in, m_dot_H / modulesInParallel, T_C_in, P_C_in, m_dot_C / modulesInParallel);
	}
	else {
		regenModel->setInletStates(T_H_in, P_H_in, m_dot_H, T_C_in, P_C_in,  m_dot_C);
	}
}

void RegenHX::setParameters(valveDesignOption::valveDesignOption operationMode, valveDesignOption::valveModes valveMode, double Q_dot_loss, double P_0, double D_s, double e_v)
{
	this->operationMode = operationMode;
	this->valveMode = valveMode;
	regenModel->setParameters(valveMode, Q_dot_loss, P_0, D_s, e_v);
}

void RegenHX::setDesignTargets(targetModes::targetModes targetMode, targetModes::target2Modes secondTargetMode, double targetParameter, double secondTargetParameter)
{
	if (targetMode != targetModes::EFF) {
		regenModel->setDesignTargets(targetMode, secondTargetMode, targetParameter / modulesInParallel, secondTargetParameter);
	} else {
		regenModel->setDesignTargets(targetMode, secondTargetMode, targetParameter, secondTargetParameter);
	}
}
