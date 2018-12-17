#include "RegenHX.h"

RegenHX::RegenHX()
{
	regenModel = new RegeneratorModel();
	valves = new valve[5];
}

RegenHX::~RegenHX()
{
	delete regenModel;
	delete[] valves;
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
	
	for (int i = 0; i < 4; i++) {
		valves[i] = solution.valves[i];
	}

	costHX *= modulesInParallel;
	if (operationMode == valveDesignOption::VDO1) {
		UA *= modulesInParallel;
		Q_dot_a *= modulesInParallel;
		m_dot_carryover *= modulesInParallel;
	}

	costHX += costValves;

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

	//High Pressure Low Temperature. T_c_in
	/*valves[0].dP = 8.6;
	valves[0].cost = 36000;
	valves[0].cv = 950;

	//High Pressure High Temperature. T_c_out
	valves[1].dP = 13;
	valves[1].cost = 74000;
	valves[1].cv = 1100;

	//Low Pressure High Temperature. T_h_in
	valves[2].dP = 16;
	valves[2].cost = 101000;
	valves[2].cv = 1750;

	//Low Pressure Low Temperature. T_h_out
	valves[3].dP = 29.6;
	
	valves[4].cost = 15000;

	if (operationMode == valveDesignOption::VDO1) {
		valves[3].cost = 36000;
		valves[3].cv = 950;
	}
	else if (operationMode == valveDesignOption::VDO2) {
		valves[3].cost = 59000;
		valves[3].cv = 1600;
	}
	
	costValves = 4 * (valves[0].cost + valves[1].cost + valves[2].cost + valves[3].cost) + 2 * valves[4].cost;

	regenModel->setValves(valves);
	*/
	targetModes::target2Modes secondTargetMode = target_2;
	double dP_max_Regen = target_2_value;
	
	/*double dP_max_Regen;
	if (secondTargetMode == targetModes::dP_max) {
		if (valveMode == valveDesignOption::FIXED_DP) {
			dP_max_Regen = max((P_h_in - P_h_out) / f_dP - valves[2].dP - valves[3].dP, (P_c_in - P_c_out) / f_dP - valves[0].dP - valves[1].dP);
		}
		else if (valveMode == valveDesignOption::FIXED_CV) {
			dP_max_Regen = max((P_h_in - P_h_out) / f_dP, (P_c_in - P_c_out) / f_dP);
		}
	}
	else if (secondTargetMode == targetModes::AR) {
		dP_max_Regen = target_2_value;
	}
	*/

	setParameters(operationMode, valveMode, Q_dot_loss, P_0, D_s, e_v);
	setInletStates(T_h_in, P_h_in, m_dot_h, T_c_in, P_c_in, m_dot_c);

	/*if (valveMode == valveDesignOption::FIXED_DP) {
		setInletStates(T_h_in, P_h_in - valves[2].dP, m_dot_h, T_c_in, P_c_in - valves[0].dP, m_dot_c);
	}
	else if (valveMode == valveDesignOption::FIXED_CV) {
		setInletStates(T_h_in, P_h_in, m_dot_h, T_c_in, P_c_in, m_dot_c);
	}*/

	if (targetType == 0) {
		setDesignTargets(targetModes::UA, secondTargetMode, targetValue, dP_max_Regen);
	}
	else if (targetType == 1) {
		setDesignTargets(targetModes::COST, secondTargetMode, targetValue - costValves, dP_max_Regen);
	}
	else if (targetType == 2) {
		setDesignTargets(targetModes::EFF, secondTargetMode, targetValue, dP_max_Regen);
	}

	int status = getDesignSolution();

	if (status < 0) {
		resetDesignStructure();
		throw("RegenHX::design",
			"Regenerator model failed!");
	}

	/*if (epsilon > eff_limit && targetType != 2) {
		ms_des_solved.m_eff_limited = true;
		double oldUA = UA;
		double oldCostHX = costHX;
		setParameters(operationMode, valveMode, Q_dot_loss, P_0, D_s, e_v);

		if (valveMode == valveDesignOption::FIXED_DP) {
			setInletStates(T_h_in, P_h_in - valves[2].dP, m_dot_h, T_c_in, P_c_in - valves[0].dP, m_dot_c);
		}
		else if (valveMode == valveDesignOption::FIXED_CV) {
			setInletStates(T_h_in, P_h_in, m_dot_h, T_c_in, P_c_in, m_dot_c);
		}

		setDesignTargets(targetModes::EFF, secondTargetMode, eff_limit, dP_max_Regen);

		status = getDesignSolution();

		if (status < 0) {
			resetDesignStructure();
			throw("RegenHX::design",
				"Regenerator model failed!");
		}

		UA = oldUA;
		costHX = oldCostHX;
	}*/

	q_dot = Q_dot_a;
	T_c_out = T_C_out;
	T_h_out = T_H_out;

	/*if (valveMode == valveDesignOption::FIXED_DP) {
		ms_des_solved.m_DP_cold_des = (m_HTR_HP_dP + valves[0].dP + valves[1].dP) * f_dP;
		ms_des_solved.m_DP_hot_des = (m_HTR_LP_dP + valves[2].dP + valves[3].dP) * f_dP;

		ms_des_solved.m_HTR_valve_HTHP_cv = valves[1].cv;
		ms_des_solved.m_HTR_valve_LTHP_cv = valves[0].cv;
		ms_des_solved.m_HTR_valve_HTLP_cv = valves[2].cv;
		ms_des_solved.m_HTR_valve_LTLP_cv = valves[3].cv;
	}
	else if (valveMode == valveDesignOption::FIXED_CV) {
		ms_des_solved.m_DP_cold_des = m_HTR_HP_dP * f_dP;
		ms_des_solved.m_DP_hot_des = m_HTR_LP_dP * f_dP;

		ms_des_solved.m_HTR_valve_HTHP_cv = valves[1].dP;
		ms_des_solved.m_HTR_valve_LTHP_cv = valves[0].dP;
		ms_des_solved.m_HTR_valve_HTLP_cv = valves[2].dP;
		ms_des_solved.m_HTR_valve_LTLP_cv = valves[3].dP;
	}
	ms_des_solved.m_HTR_valve_HTHP_cv = valves[1].dP;
	ms_des_solved.m_HTR_valve_LTHP_cv = valves[0].dP;
	ms_des_solved.m_HTR_valve_HTLP_cv = valves[2].dP;
	ms_des_solved.m_HTR_valve_LTLP_cv = valves[3].dP;
	
	ms_des_solved.m_eff_design = epsilon;
	ms_des_solved.m_min_DT_design = min((T_H_in - T_H_out), (T_C_out - T_C_in));
	ms_des_solved.m_NTU_design = NTU_R_e;
	ms_des_solved.m_Q_dot_design = Q_dot_a;
	ms_des_solved.m_T_c_out = T_C_out;
	ms_des_solved.m_T_h_out = T_H_out;
	ms_des_solved.m_HTR_AR = AspectRatio;
	ms_des_solved.m_HTR_D_fr = D_fr;
	ms_des_solved.m_HTR_L = L;

	if (targetType == 0 || targetType == 2) {
		ms_des_solved.m_UA_design_total = UA;
	}
	else if (targetType == 1) {
		ms_des_solved.m_UA_design_total = costHX;
	}
	

	ms_des_solved.m_aUA_design_total = UA;
	ms_des_solved.m_cost_design_total = costHX;
	ms_des_solved.m_m_dot_carryover = m_dot_carryover;*/
}

void RegenHX::off_design_solution(double T_c_in, double P_c_in, double m_dot_c, double P_c_out, double T_h_in, double P_h_in, double m_dot_h, double P_h_out, double & q_dot, double & T_c_out, double & T_h_out)
{
	setInletStates(T_h_in, P_h_in, m_dot_h, T_c_in, P_c_in, m_dot_c);
	regenModel->setOutletStates(T_c_out, P_h_out, T_c_out, P_c_out);
	

	/*ms_od_solved.m_eff = epsilon;
	ms_od_solved.m_min_DT = min((T_H_in - T_H_out), (T_C_out - T_C_in));
	ms_od_solved.m_NTU = NTU_R_e;
	ms_od_solved.m_P_c_out = P_C_out;
	ms_od_solved.m_P_h_out = P_H_out;
	ms_od_solved.m_q_dot = Q_dot_calc;
	ms_od_solved.m_T_c_out = T_C_out;
	ms_od_solved.m_T_h_out = T_H_out;
	ms_od_solved.m_UA_total = UA;*/
}

void RegenHX::resetDesignStructure()
{
	ms_des_solved.m_DP_cold_des = ms_des_solved.m_DP_hot_des = ms_des_solved.m_eff_design = ms_des_solved.m_min_DT_design = ms_des_solved.m_NTU_design =
		ms_des_solved.m_Q_dot_design = ms_des_solved.m_T_c_out = ms_des_solved.m_T_h_out = ms_des_solved.m_UA_design_total = 
		ms_des_solved.m_aUA_design_total = ms_des_solved.m_cost_design_total =
		ms_des_solved.m_m_dot_carryover = std::numeric_limits<double>::quiet_NaN();
	ms_des_solved.m_eff_limited = false;
}

void RegenHX::setInletStates(double T_H_in, double P_H_in, double m_dot_H, double T_C_in, double P_C, double m_dot_C)
{
	this->T_C_in = T_C_in;
	this->T_H_in = T_H_in;
	this->m_dot_H = m_dot_H;
	this->m_dot_C = m_dot_C;

	if (this->operationMode == valveDesignOption::VDO1) {
		regenModel->setInletStates(T_H_in, P_H_in, m_dot_H / modulesInParallel, T_C_in, P_C, m_dot_C / modulesInParallel);
	}
	else {
		regenModel->setInletStates(T_H_in, P_H_in, m_dot_H, T_C_in, P_C,  m_dot_C);
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
