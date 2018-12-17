#include "RegeneratorModel.h"

string const RegeneratorModel::LOG_FILEPATH = string(BUILD_PATH) + "/RegenHX_LOG.log";
string const RegeneratorModel::PROPERTY_FILES = string(BUILD_PATH) + "/PropertyFiles/";
string const RegeneratorModel::SPHERES_RP_TABLE_PATH = PROPERTY_FILES + "Spheres_RP.csv";
string const RegeneratorModel::BALANCED_REGENERATOR_TABLE_PATH = PROPERTY_FILES + "balanced-regenerator.csv";
string const RegeneratorModel::FATIGUE_TABLE_PATH = PROPERTY_FILES + "fatigue.csv";


RegeneratorModel::RegeneratorModel()
{
	loadTables();

	HeatTransfer = new MonoSolver<RegeneratorModel>();
	HotPressureDrop = new MonoSolver<RegeneratorModel>();
	ColdPressureDrop = new MonoSolver<RegeneratorModel>();
	Length = new MonoSolver<RegeneratorModel>();
	Diameter = new MonoSolver<RegeneratorModel>();
	WallThickness = new MonoSolver<RegeneratorModel>();
	CarryoverMassFlow = new MonoSolver<RegeneratorModel>();

	PressureSplit = new MonoSolver<RegeneratorModel>();
	Valve = new MonoSolver<RegeneratorModel>();

	valves = new valve[4];

	if (spdlog::get("logger") == nullptr) {
		auto logger = spdlog::basic_logger_mt("logger", LOG_FILEPATH);
		spdlog::flush_on(spdlog::level::info);
		spdlog::drop("logger");
		spdlog::register_logger(logger);
		//spdlog::set_pattern("[%H:%M:%S] [%l] %v");
		spdlog::set_pattern("%v");
	}
}

RegeneratorModel::~RegeneratorModel()
{
	spdlog::drop("logger");

	delete bedMaterialTable;
	delete shellMaterialTable;
	delete insulationMaterialTable;
	delete regeneratorTable;
	delete spheresRPTable;
	delete fatigueTable;

	delete HeatTransfer;
	delete HotPressureDrop;
	delete ColdPressureDrop;
	delete Length;
	delete Diameter;
	delete WallThickness;
	delete CarryoverMassFlow;

	delete PressureSplit;
	delete Valve;

	delete[] valves;
}

void RegeneratorModel::packedspheresNdFit(double Re, double * f, double * j_H)
{
	if (Re < 20 || Re >  50000) {
		(*f) = -1;
		(*j_H) = -1;
		throw invalid_argument("Re should be between 20 and 50,000!");
	}

	//Re is coloumn 0, f is coloumn 1
	(*f) = spheresRPTable->getValue("f", "Re", Re);
	(*j_H) = 0.23*pow(Re, -0.3);
}

double RegeneratorModel::hx(double NTU, double C_1, double C_2)
{
	if (C_1 <= 0 || C_2 <= 0) {
		throw invalid_argument("C_1 and C_2 should be greater than 0!");
	}

	double U = C_1 / C_2;	//[-]

	if (U > 1000) {
		throw invalid_argument("No value provided for U > 1000!");
	}

	return regeneratorTable->getValue(NTU, U);
}

void RegeneratorModel::getpropsregenFitCO2(double T, double P, double * rho, double * mu, double * k, double * Pr, double * Cp)
{
	CO2_state CO2States;

	int error = CO2_TP(T, P, &CO2States);

	if (error != 0)
	{
		throw invalid_argument("Error fixing CO2 state with T and P provided!");
	}

	*Cp = CO2States.cp;
	*rho = CO2States.dens;
	*mu = CO2_visc(*rho, T) / 1000000; //Convert from uPa-s to Pa-s
	*k = CO2_cond(*rho, T);

	*Pr = (*Cp)*(*mu) / (*k)* 1000.0;	//Convert MPa to kPa
}

void RegeneratorModel::packedspheresFitCO2(double m_dot, double d, double A_fr, double L, double T, double P, double porosity, 
																			double * f, double * h, double * NTU, double * DP)
{
	double rho, mu, k, Pr, Cp;
	try {
		getpropsregenFitCO2(T, P, &rho, &mu, &k, &Pr, &Cp);
	}
	catch (const invalid_argument& e) {
		throw e;
	}

	double r_h = porosity * d / (6.0*(1 - porosity));	//[m]
	double G = m_dot / (porosity*A_fr);	//[kg/m^2-s]
	double Re = G * (4 * r_h) / mu;	//[-]

	double j_H;	//[-]
	try {
		packedspheresNdFit(Re, f, &j_H);
	}
	catch (const invalid_argument& e) {
		throw e;
	}


	double St = j_H / pow(Pr, (2.0 / 3.0));	//[-]
	*NTU = St * L / r_h;	//[-]
	*h = j_H * G*Cp / pow(Pr, (2.0 / 3.0)) * 1000.0;	//Convert kW to W. [W]


	double V = m_dot / rho / (A_fr);
	Re = d * rho*V / mu;
	double HELPER = (1 - porosity);
	double Re_m = Re / HELPER;
	double f_1L = 136 / pow(HELPER, 0.38);
	double f_1T = 29 / (pow(HELPER, 1.45) * pow(porosity, 2));
	double f_2 = 1.87 * pow(porosity, 0.75) / pow(HELPER, 0.26);
	double q = exp(-pow(porosity, 2) * HELPER / 12.6 * Re_m);
	double f_p = (q * f_1L / Re_m + (1 - q) * (f_2 + f_1T / Re_m)) * HELPER / pow(porosity, 3);
	*DP = f_p * rho * pow(V, 2) * L / d / 1000.0; //Pa to kPa
}

void RegeneratorModel::massflowVariablesInit()
{
	C_dot_H = C_p_H * m_dot_H;
	C_dot_C = C_p_C * m_dot_C;

	C_dot_min = min(C_dot_H, C_dot_C);
	C_dot_max = max(C_dot_H, C_dot_C);

	C_R = C_dot_min / C_dot_max;
}

void RegeneratorModel::calculateModel()
{
	//							[SECTION 1]
	//Pressure at hot exit
	P_H_out = P_H_in - dP_H;	//[kPa]
							//Pressure at cold exit
	P_C_out = P_C - dP_C;	//[kPa]

	if (P_C_out <= N_co2_props::P_lower_limit || P_H_out <= N_co2_props::P_lower_limit || P_C_out >= N_co2_props::P_upper_limit || P_H_out >= N_co2_props::P_upper_limit) {
		throw invalid_argument("Outlet pressures are either below " + std::to_string(N_co2_props::P_lower_limit) + "[kPa] or above " + std::to_string(N_co2_props::P_upper_limit) + "[kPa]!");
	}

	CO2_TP(T_C_in, P_H_out, &CO2State);
	h_H_out_max = CO2State.enth;

	CO2_TP(T_H_in, P_C_out, &CO2State);
	h_C_out_max = CO2State.enth;

	//Maximum heat transfer if hot stream exiting at cold temp
	Q_dot_max_H = m_dot_H * (h_H_in - h_H_out_max);

	//Maximum heat transfer if cold stream exiting at hot temp
	Q_dot_max_C = m_dot_C * (h_C_out_max - h_C_in);

	//Actual maximum possible heat transfer
	Q_dot_max = min(Q_dot_max_H, Q_dot_max_C);

	if (T_H_out <= N_co2_props::T_lower_limit || T_H_out >= N_co2_props::T_upper_limit) {
		throw invalid_argument("T_H_out is either below " + std::to_string(N_co2_props::T_lower_limit) + "[K] or above " + std::to_string(N_co2_props::T_upper_limit) + "[K]!");
	}

	CO2_TP(T_H_out, P_H_out, &CO2State);
	h_H_out = CO2State.enth;

	Q_dot_a = m_dot_H * (h_H_in - h_H_out);	

	h_C_out = h_C_in + Q_dot_a / m_dot_C;
	CO2_PH(P_C_out, h_C_out, &CO2State);

	T_C_out = CO2State.temp;

	Q_dot = Q_dot_a + Q_dot_loss;

	//							[SECTION 2]

	D_fr = cbrt(4 * V_0 / PI / AR);
	L = D_fr * AR;
	//A_fr = V_0 / L;
	//D_fr = sqrt(4 * A_fr / PI);
	A_fr = PI * pow(D_fr, 2) / 4;
	T_H_f = (T_H_in + T_H_out) / 2;
	T_C_f = (T_C_in + T_C_out) / 2;

	try {
		packedspheresFitCO2(m_dot_H, D_s, A_fr, L, T_H_f, P_H_in, e_v, &f_H, &h_H, &NTU_H, &dP_H_calc);
		packedspheresFitCO2(m_dot_C, D_s, A_fr, L, T_C_f, P_C, e_v, &f_C, &h_C, &NTU_C, &dP_C_calc);
	}
	catch (const invalid_argument& e) {
		throw e;
	}

	//V_0 = A_fr * L;

	m_s = V_0 * (1 - e_v) * rho_s;

	A_s = 6.0 * (1.0 - e_v)*V_0 / D_s;

	UA = 2 * A_s * h_H * h_C / (h_H + h_C) / 1000.0;

	NTU_R = UA / C_dot_min;

	C_m = 2.0 * m_s * c_s / P_0 / C_dot_min;

	NTU_R_e = 2.0 * C_R * NTU_R / (1.0 + C_R);

	C_m_e = 2.0 * C_R*C_m / (1.0 + C_R);

	try {
		if (C_m_e < 0.001) {
			throw invalid_argument("No value provided for U > 1000!");
		}

		epsilon_1 = regeneratorTable->getValue(NTU_R_e, 1 / C_m_e);
	}
	catch (exception e) {
		throw e;
	}

	X = (1.0 - pow(C_R, 2)) / (2.0 * C_R)*(epsilon_1 / (1.0 - epsilon_1));

	epsilon = (1.0 - exp(-X)) / (1.0 - C_R * exp(-X));

	/*if (epsilon < 0 || epsilon > 1) {
		throw invalid_argument("Epsilon is not between 0 and 1!");
	}*/

	Q_dot_calc = epsilon * Q_dot_max;
	Q_dot_a_calc = Q_dot_calc - Q_dot_loss;
	dP_max = max(dP_C, dP_H);
}

void RegeneratorModel::setValves(valve * valves)
{
	for (int i = 0; i < 4; i++) {
		this->valves[i] = valves[i];
		
	}
}

void RegeneratorModel::calcValvePressureDrops()
{
	valves[2].m_dot = valves[3].m_dot = m_dot_H;
	valves[1].m_dot = m_dot_C;
	valves[0].m_dot = m_dot_C + m_dot_carryover;
	
	valves[0].P_in = P_C;
	valves[0].T_in = T_C_in;
	valves[1].P_in = P_C_out;
	valves[1].T_in = T_C_out;
	valves[2].P_in = P_H_in;
	valves[2].T_in = T_H_in;
	valves[3].P_in = P_H_out;
	valves[3].T_in = T_H_out;
	
	Valve_SP.solverName = "Valve";
	Valve_SP.target = 0;
	Valve_SP.guessValue1 = 25;		Valve_SP.guessValue2 = 50;
	Valve_SP.lowerBound = 0.1;					Valve_SP.upperBound = P_H_in;
	Valve_SP.tolerance = 0.001;
	Valve_SP.iterationLimit = 50;
	Valve_SP.isErrorRel = false;
	Valve_SP.classInst = this;
	Valve_SP.monoEquation = &RegeneratorModel::Valve_ME;
	
	Valve->setParameters(&Valve_SP);
	
		for (int i = 0; i < 4; i++) {
			valveIndex = i;
			Valve->solve();
		
		}
	}

void RegeneratorModel::calculateValveCvs()
{
	valves[2].m_dot = valves[3].m_dot = m_dot_H;
	valves[1].m_dot = m_dot_C;
	valves[0].m_dot = m_dot_C + m_dot_carryover;

	valves[0].P_in = P_C;
	valves[0].T_in = T_C_in;
	valves[1].P_in = P_C_out;
	valves[1].T_in = T_C_out;
	valves[2].P_in = P_H_in;
	valves[2].T_in = T_H_in;
	valves[3].P_in = P_H_out;
	valves[3].T_in = T_H_out;

	double f_t = 100;
	//"Valve Constants for fully open flow"
	double F_p = 1;
	double x_t0 = 0.24;

	//"Valve parameters"
	double x_t = x_t0 * ((1.0929e-09)*pow(f_t, 5) - (2.9987e-07)*pow(f_t, 4) + (2.4587e-05)*pow(f_t, 3) - (6.4267e-05)*pow(f_t, 2) - (7.0864E-02)*f_t + 3.1976);

	double rho_in, k_in, Y, C_v;
	for (int index = 0; index < 4; index++) {
		CO2_TP(valves[index].T_in, valves[index].P_in, &CO2State);
		rho_in = CO2State.dens;
		k_in = CO2State.cp / CO2State.cv;

		//"Gas expansion factor"
		Y = 1 - (valves[index].dP / valves[index].P_in) / (3 * (k_in / 1.4)*x_t);

		//"Flow rate - pressure drop relationship"
		C_v = sqrt(100 / rho_in / valves[index].dP) * valves[index].m_dot * 3600 / 27.3 / F_p / Y;
		valves[index].cv = C_v / ((4.0064e-10)*pow(f_t, 5) - (1.0249e-07)*pow(f_t, 4) + (7.4964e-06)*pow(f_t, 3) - (8.2991e-05)*pow(f_t, 2) + (5.7679E-03)*f_t);
		//double dP_calc = 1 / rho_in * pow((valves[valveIndex].m_dot * 3600 / 27.3 / F_p / C_v / Y), 2) * 100;
	}
}

void RegeneratorModel::loadTables() {
	bedMaterialTable = new LookupTable_1D(PROPERTY_FILES + bedMaterialName + ".csv");
	shellMaterialTable = new LookupTable_1D(PROPERTY_FILES + shellMaterialName + ".csv");
	insulationMaterialTable = new LookupTable_1D(PROPERTY_FILES + insulationMaterialName + ".csv");
	regeneratorTable = new LookupTable_2D(BALANCED_REGENERATOR_TABLE_PATH);
	spheresRPTable = new LookupTable_1D(SPHERES_RP_TABLE_PATH);
	fatigueTable = new LookupTable_1D(FATIGUE_TABLE_PATH);
}

int RegeneratorModel::calculateCost()
{
	D_shell = D_fr + insulationThickness * 2;
	R_i = D_shell / 2.0;
	
	double tolerance;
	int statusSolver = WallThickness->solve(&tolerance);
	if (statusSolver != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	volumeShell = 2 * (PI * (pow(R_o, 2) - pow(R_i, 2)) * L + 4 / 3.0 * PI * (pow(R_o, 3) - pow(R_i, 3))); // a cylinder plus two half spheres
	volumeInsulation = 2* (PI * (pow(R_i, 2) - pow(R_i - insulationThickness, 2)) * L + 4 / 3.0 * PI * pow(R_i, 3)*insulationParameter);
	volumeBed = 2* V_0 * (1 - e_v);

	double shellMaterialDensity = shellMaterialTable->getValue("rho", "T", T_f);
	double insulationMaterialDensity = insulationMaterialTable->getValue("rho", "T", T_f);
	double bedMaterialDensity = bedMaterialTable->getValue("rho", "T", T_f);

	massShell = shellMaterialDensity * volumeShell;
	massInsulation = insulationMaterialDensity * volumeInsulation;
	massBed = bedMaterialDensity * volumeBed;

	costShellMaterial = specificCostShellMaterial * massShell;
	costInsulationMaterial = specificCostInsulationMaterial * massInsulation;
	costBedMaterial = specificCostBedMaterial * massBed;

	costMaterial = costShellMaterial + costInsulationMaterial + costBedMaterial;
	costModule = costMaterial + 2 * (priceCasting + priceCastingSteel + priceWelding);
	return 0;
}

void RegeneratorModel::setInletStates(double T_H_in, double P_H_in, double m_dot_H, double T_C_in, double P_C, double m_dot_C)
{
	this->T_H_in = T_H_in;
	this->P_H_in = P_H_in;
	this->m_dot_H = m_dot_H;
	this->T_C_in = T_C_in;
	this->P_C = P_C;
	this->m_dot_C = m_dot_C;
	m_dot_carryover = 0;

	if (P_C <= N_co2_props::P_lower_limit || P_H_in <= N_co2_props::P_lower_limit || P_C >= N_co2_props::P_upper_limit || P_H_in >= N_co2_props::P_upper_limit) {
		throw invalid_argument("Inlet pressures are either below " + std::to_string(N_co2_props::P_lower_limit) + "[kPa] or above " + std::to_string(N_co2_props::P_upper_limit) + "[kPa]!");
	}
	if (T_H_in <= N_co2_props::T_lower_limit || T_C_in <= N_co2_props::T_lower_limit || T_H_in >= N_co2_props::T_upper_limit || T_C_in >= N_co2_props::T_upper_limit) {
		throw invalid_argument("Inlet temperatures are either below " + std::to_string(N_co2_props::T_lower_limit) + "[K] or above " + std::to_string(N_co2_props::T_upper_limit) + "[K]!");
	}

	//Enthalpy at hot inlet
	CO2_TP(T_H_in, P_H_in, &CO2State);
	h_H_in = CO2State.enth;

	//Enthalpy at cold inlet
	CO2_TP(T_C_in, P_C, &CO2State);
	h_C_in = CO2State.enth;

	CO2_TP(T_C_in, P_H_in, &CO2State);
	h_H_out_max_p = CO2State.enth;

	CO2_TP(T_H_in, P_C, &CO2State);
	h_C_out_max_p = CO2State.enth;

	C_p_H = (h_H_in - h_H_out_max_p) / (T_H_in - T_C_in);
	C_p_C = (h_C_out_max_p - h_C_in) / (T_H_in - T_C_in);

	T_f = (T_H_in + T_C_in) / 2;
	rho_s = bedMaterialTable->getValue("rho", "T", T_f);
	c_s = bedMaterialTable->getValue("c", "T", T_f);

	massflowVariablesInit();
}

void RegeneratorModel::setOutletStates(double T_H_out, double P_H_out, double T_C_out, double P_C_out)
{
	this->T_H_out = T_H_out;
	this->dP_H = P_H_in - P_H_out;
	this->T_C_out = T_C_out;
	this->dP_C = P_C - P_C_out;
}

void RegeneratorModel::setParameters(valveDesignOption::valveModes valveMode, double Q_dot_loss, double P_0, double D_s, double e_v)
{
	this->valveMode = valveMode;
	this->Q_dot_loss = Q_dot_loss;
	this->D_s = D_s;
	this->e_v = e_v;
	this->P_0 = P_0;

	numberOfCycles = 1.2 * operationYears * 365 * operationHoursPerDay * 60 * 60 / P_0;
	stressAmplitude = fatigueTable->getValue("Sa", "N", numberOfCycles) * 6.89475729; //Convert ksi to MPA
}

void RegeneratorModel::setDesignTargets(targetModes::targetModes targetMode, targetModes::target2Modes secondTargetMode, double targetParameter, double secondTargetParameter)
{
	this->secondTargetMode = secondTargetMode;
	if (secondTargetMode == targetModes::target2Modes::dP_max) {
		this->targetdP_max_Regen = secondTargetParameter;
	}
	else if (secondTargetMode == targetModes::target2Modes::AR) {
		targetAR = secondTargetParameter;
	}
	
	this->targetMode = targetMode;
	this->targetParameter = targetParameter;
}

double RegeneratorModel::densityIntegral(double T_low, double T_high, double P)
{
	double dT = 0.001 * (T_high - T_low) / L;
	double T_int = T_low - dT;
	int N = abs((int)((T_high - T_low) / dT));
	double integral = 0;

	for (int i = 0; i <= N; i++) {
		T_int += dT;
		CO2_TP(T_int, P, &CO2State);
		integral += CO2State.dens;
	}

	return fabs(integral * dT);
}

void RegeneratorModel::calcCarryoverMassFlow()
{
	CO2_TP(T_C_out, P_C, &CO2State);
	double rho_H_H_extra = CO2State.dens;

	CO2_TP(T_C_in, P_C, &CO2State);
	double rho_H_C_extra = CO2State.dens;

	CO2_TP(T_H_in, P_H_in, &CO2State);
	double rho_L_H_extra = CO2State.dens;

	CO2_TP(T_H_out, P_H_in, &CO2State);
	double rho_L_C_extra = CO2State.dens;

	double mass_H = vol_extra / 4 * (rho_H_H_extra - rho_L_H_extra);
	double mass_C = vol_extra / 4 * (rho_H_C_extra - rho_L_C_extra);

	double integral_H, integral_L;
	
	double dT = 32.1165431 - 1.33843643 * (T_C_in - 273.15) + 0.723818378 * (T_H_in - 273.15)
		+ 0.00222683779 * P_C - 0.000485963027 * P_H_in - 17.4313839 * C_m_e;

	integral_H = densityIntegral(T_C_in, T_C_out - dT, P_C) * L / fabs(T_C_out - dT - T_C_in);
	integral_L = densityIntegral(T_H_out + dT, T_H_in, P_H_in) * L / fabs(T_H_in - T_H_out - dT);
	

	double mass = 2 * (e_v * pow(D_fr, 2) / 4 * PI *(integral_H - integral_L) + (mass_C + mass_H)) / CO;
	m_dot_carryover = mass / P_0;
}

void RegeneratorModel::carryoverEnthDrop()
{
	h_H_out = (h_C_in * m_dot_carryover + h_H_out * (m_dot_H - m_dot_carryover)) / m_dot_H;
	CO2_PH(P_H_out, h_H_out, &CO2State);
	T_H_out = CO2State.temp;
}

int RegeneratorModel::HeatTransfer_ME(double T_H_out, double * QdotAsDifference) {
	this->T_H_out = T_H_out;

	if (T_H_out > T_H_in || T_H_out < T_C_in) {
		return -1;
	}

	try {
		calculateModel();
	}
	catch (const invalid_argument& e) {
		return -1;
	}

	*QdotAsDifference = Q_dot_a_calc - Q_dot_a;
	return 0;
}

int RegeneratorModel::HotPressureDrop_ME(double dP_H, double * dP_HsDifference)
{
	this->dP_H = dP_H;

	if (HeatTransfer->solve() != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	*dP_HsDifference = dP_H_calc - this->dP_H;
	return 0;
}

int RegeneratorModel::ColdPressureDrop_ME(double dP_C, double * dP_CsDifference)
{
	this->dP_C = dP_C;

	if (HotPressureDrop->solve() != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	*dP_CsDifference = dP_C_calc - this->dP_C;
	return 0;
}

int RegeneratorModel::Length_ME(double L, double * dP_max)
{
	this->L = L;

	int status = ColdPressureDrop->solve();
	

	if (status != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	*dP_max = this->dP_max;

	return 0;
}

int RegeneratorModel::Diameter_dP_ME(double D_fr, double * targetParameter)
{
	this->D_fr = D_fr;

	int status = Length->solve();

	if (status != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	if (targetMode == targetModes::COST) {
		calculateCost();
		*targetParameter = costModule;
	}
	else if (targetMode == targetModes::EFF) {
		*targetParameter = epsilon;

	}
	else if (targetMode == targetModes::UA) {
		*targetParameter = UA;
	}

	return 0;
}

int RegeneratorModel::Diameter_AR_ME(double D_fr, double * targetParameter)
{
	this->D_fr = D_fr;
	this->L = D_fr * targetAR;

	int status = ColdPressureDrop->solve();

	if (status != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	if (targetMode == targetModes::COST) {
		calculateCost();
		*targetParameter = costModule;
	}
	else if (targetMode == targetModes::EFF) {
		*targetParameter = epsilon;

	}
	else if (targetMode == targetModes::UA) {
		*targetParameter = UA;
	}

	return 0;
}

int RegeneratorModel::CarryoverMassFlow_FIXED_CV_ME(double m_dot_carryover, double *comass_difference)
{
	this->m_dot_C -= m_dot_carryover;
	this->m_dot_H -= m_dot_carryover;

	massflowVariablesInit();

	int status = PressureSplit->solve();

	this->m_dot_C += m_dot_carryover;
	this->m_dot_H += m_dot_carryover;
	if (status != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	calcCarryoverMassFlow();

	*comass_difference = m_dot_carryover - this->m_dot_carryover;
	return 0;
}

int RegeneratorModel::CarryoverMassFlow_FIXED_DP_ME(double m_dot_carryover, double *comass_difference)
{
	this->m_dot_C -= m_dot_carryover;
	this->m_dot_H -= m_dot_carryover;

	massflowVariablesInit();
	
	int status = Diameter->solve();

	this->m_dot_C += m_dot_carryover;
	this->m_dot_H += m_dot_carryover;
	if (status != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	calcCarryoverMassFlow();

	*comass_difference = m_dot_carryover - this->m_dot_carryover;
	return 0;
}

int RegeneratorModel::WallThickness_ME(double th, double * stressAmplitude)
{
	wallThickness = th;

	R_o = R_i + wallThickness;

	double RiSquared = pow(R_i, 2);

	double RoSquared = pow(R_o, 2);

	double RSquareDifference = RoSquared - RiSquared;

	double RSquareProduct = RoSquared * RiSquared;

	double magicPieceRi = RiSquared * RSquareDifference;

	double magicPieceRo = RoSquared * RSquareDifference;

	double magicPieceHigh = RSquareProduct * (Patm - P_C / 1000.0);

	double magicPieceLow = RSquareProduct * (Patm - P_H_out / 1000.0);

	//sigma_a_h = (P_high*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2);		
	double sigma_a_h = (P_C / 1000.0 * RiSquared - Patm * RoSquared) / RSquareDifference; //axial stress calculation

	//sigma_c_1_h = (P_high*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) - r_i ^ 2 * r_o ^ 2 * (P_o - P_high) / (r_i ^ 2 * (r_o ^ 2 - r_i ^ 2));	//hoop stress at inner surface
	double sigma_c_1_h = sigma_a_h - magicPieceHigh / magicPieceRi;	//hoop stress at inner surface

	//sigma_c_2_h = (P_high*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) - r_i ^ 2 * r_o ^ 2 * (P_o - P_high) / (r_o ^ 2 * (r_o ^ 2 - r_i ^ 2));	//hoop stress at outer surface
	double sigma_c_2_h = sigma_a_h - magicPieceHigh / magicPieceRo;	//hoop stress at outer surface

	//sigma_r_1_h = (P_high*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) + r_i ^ 2 * r_o ^ 2 * (P_o - P_high) / (r_i ^ 2 * (r_o ^ 2 - r_i ^ 2));	//radial stress at inner surface
	double sigma_r_1_h = sigma_a_h + magicPieceHigh / magicPieceRi;	//radial stress at inner surface

	//sigma_r_2_h = (P_high*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) + r_i ^ 2 * r_o ^ 2 * (P_o - P_high) / (r_o ^ 2 * (r_o ^ 2 - r_i ^ 2));	//radial stress at outer surface
	double sigma_r_2_h = sigma_a_h + magicPieceHigh / magicPieceRo;	//radial stress at outer surface

	double sigma_prime_1_h = sqrt((pow((sigma_r_1_h - sigma_c_1_h), 2) + pow((sigma_c_1_h - sigma_a_h), 2) + pow((sigma_a_h - sigma_r_1_h), 2)) / 2.0);	//equivalent stress at inner surface
	double sigma_prime_2_h = sqrt((pow((sigma_r_2_h - sigma_c_2_h), 2) + pow((sigma_c_2_h - sigma_a_h), 2) + pow((sigma_a_h - sigma_r_2_h), 2)) / 2.0);	//equivalent stress at outer surface
	double sigma_prime_h = max(sigma_prime_1_h, sigma_prime_2_h);		//largest equivalent stress


	// sigma_a_l = (P_low*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2);	//axial stress calculation
	double sigma_a_l = (P_H_out / 1000.0 * RiSquared - Patm * RoSquared) / RSquareDifference;	//axial stress calculation

	//sigma_c_1_l = (P_low*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) - r_i ^ 2 * r_o ^ 2 * (P_o - P_low) / (r_i ^ 2 * (r_o ^ 2 - r_i ^ 2));	//hoop stress at inner surface
	double sigma_c_1_l = sigma_a_l - magicPieceLow / magicPieceRi;	//hoop stress at inner surface

	//sigma_c_2_l = (P_low*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) - r_i ^ 2 * r_o ^ 2 * (P_o - P_low) / (r_o ^ 2 * (r_o ^ 2 - r_i ^ 2));	//hoop stress at outer surface
	double sigma_c_2_l = sigma_a_l - magicPieceLow / magicPieceRo;	//hoop stress at outer surface

	//sigma_r_1_l = (P_low*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) + r_i ^ 2 * r_o ^ 2 * (P_o - P_low) / (r_i ^ 2 * (r_o ^ 2 - r_i ^ 2));	//radial stress at inner surface
	double sigma_r_1_l = sigma_a_l + magicPieceLow / magicPieceRi;	//radial stress at inner surface

	//sigma_r_2_l = (P_low*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) + r_i ^ 2 * r_o ^ 2 * (P_o - P_low) / (r_o ^ 2 * (r_o ^ 2 - r_i ^ 2));	//radial stress at outer surface
	double sigma_r_2_l = sigma_a_l + magicPieceLow / magicPieceRo;	//radial stress at outer surface

	double sigma_prime_1_l = sqrt((pow((sigma_r_1_l - sigma_c_1_l), 2) + pow((sigma_c_1_l - sigma_a_l), 2) + pow((sigma_a_l - sigma_r_1_l), 2)) / 2.0);	//equivalent stress at inner surface
	double sigma_prime_2_l = sqrt((pow((sigma_r_2_l - sigma_c_2_l), 2) + pow((sigma_c_2_l - sigma_a_l), 2) + pow((sigma_a_l - sigma_r_2_l), 2)) / 2.0);	//equivalent stress at outer surface
	double sigma_prime_l = max(sigma_prime_1_l, sigma_prime_2_l);		//largest equivalent stress

	*stressAmplitude = (sigma_prime_h - sigma_prime_l) / 2.0;
	
	return 0;
}

int RegeneratorModel::Valve_ME(double dP, double *dP_difference)
{
	valves[valveIndex].dP = dP;
	
	double f_t = 100;
	CO2_TP(valves[valveIndex].T_in, valves[valveIndex].P_in, &CO2State);
	double rho_in = CO2State.dens;
	double k_in = CO2State.cp / CO2State.cv;
	
	//"Valve Constants for fully open flow"
	double F_p = 1;
	double x_t0 = 0.24;
	
	//"Valve parameters"
	double C_v = valves[valveIndex].cv * ((4.0064e-10)*pow(f_t, 5) - (1.0249e-07)*pow(f_t, 4) + (7.4964e-06)*pow(f_t, 3) - (8.2991e-05)*pow(f_t, 2) + (5.7679E-03)*f_t);
	double x_t = x_t0 * ((1.0929e-09)*pow(f_t, 5) - (2.9987e-07)*pow(f_t, 4) + (2.4587e-05)*pow(f_t, 3) - (6.4267e-05)*pow(f_t, 2) - (7.0864E-02)*f_t + 3.1976);
	
	//"Gas expansion factor"
	double Y = 1 - (dP / valves[valveIndex].P_in) / (3 * (k_in / 1.4)*x_t);
	
	//"Flow rate - pressure drop relationship"
	double dP_calc = 1 / rho_in * pow((valves[valveIndex].m_dot * 3600 / 27.3 / F_p / C_v / Y), 2) * 100;
	
	*dP_difference = dP - dP_calc;
	
	return 0;
}

int RegeneratorModel::PressureSplit_ME(double regenMaxDrop_guess, double *regenMaxTotalDrop)
{
	targetdP_max_Regen = regenMaxDrop_guess;
	Length->updateTarget(targetdP_max_Regen);
	
	int status = Diameter->solve();
	if (status != C_monotonic_eq_solver::CONVERGED) {
		return -1;
		
	}
	
	calcValvePressureDrops();
	dP_C_total = dP_C + valves[0].dP + valves[1].dP;
	dP_H_total = dP_H + valves[2].dP + valves[3].dP;
	targetdP_total = max(dP_C_total, dP_H_total);
	
	*regenMaxTotalDrop = this->targetdP_total;
	
	return 0;
}

void RegeneratorModel::setSolver(Eigen::VectorXd& params) {
	this->T_H_out = params(0)*T_H_in;
	this->AR = params(1)*max_Size;
	this->dP_C = params(2)*P_C;
	this->V_0 = params(3)*max_Size;
}

void RegeneratorModel::jacobian(int n, Eigen::VectorXd x_n, Eigen::VectorXd dx, Eigen::MatrixXd& Jn) {
	Eigen::VectorXd f_i(n);
	Eigen::VectorXd f_im1(n);
	Eigen::VectorXd f_ip1(n);
	Eigen::VectorXd x_tmp(n);
	evaluate(x_n, f_i);

	x_tmp = x_n;

	// Through variables
	for (int i = 0; i < n; i++) {
		x_tmp(i) -= dx(i);
		evaluate(x_tmp, f_im1);
		x_tmp(i) += 2*dx(i);
		evaluate(x_tmp, f_ip1);
		Jn(i, i) = (f_im1(i) + f_ip1(i) - 2 * f_i(i)) / pow(dx(i), 2);
	}
}

void RegeneratorModel::evaluate(Eigen::VectorXd point, Eigen::VectorXd& f) {
	setSolver(point);
	calculateModel();

	f(0) = (Q_dot_a - Q_dot_a_calc) / (Q_dot_a + Q_dot_a_calc);
	f(1) = (dP_H - dP_H_calc)/P_H_in;
	f(2) = (dP_C - dP_C_calc)/P_C;
	f(3) = (UA - targetParameter) / (targetParameter * 3);
}

int RegeneratorModel::getDesignSolution()
{
	clock_t begin = clock();
	int n = 4;

	Eigen::VectorXd x_0(n);
	Eigen::VectorXd x_1(n);

	dP_H = targetdP_max_Regen;

	x_0(0) = (T_H_in + T_C_in)/2/T_H_in;
	x_0(1) = 1.1/max_Size;
	x_0(2) = 0.3*targetdP_max_Regen/P_C;
	x_0(3) = 0.5/max_Size;

	Eigen::VectorXd dx(n);
	dx(0) = 1/T_H_in;
	dx(1) = 0.1 / max_Size;
	dx(2) = 1/P_C;
	dx(3) = 0.1/max_Size;

	Eigen::MatrixXd Jn(n, n);
	jacobian(n, x_0, dx, Jn);

	
	Eigen::MatrixXd Jn_inv(n, n);
	Eigen::MatrixXd Jn_tmp(n, n);
	Eigen::VectorXd f_0(n);
	Eigen::VectorXd f_1(n);
	Eigen::VectorXd del_f(n);
	Eigen::VectorXd del_x(n);
	evaluate(x_0, f_0);
	
	
	/*spdlog::get("logger")->info("T_H_out/T_H_in = " + to_string(x_0(0)) + ", T_H_out = " + to_string(x_0(0)*T_H_in));
	spdlog::get("logger")->info("AR = " + to_string(x_0(1)*max_Size) + ", L = " + to_string(L) + ", D_fr = " + to_string(D_fr));
	spdlog::get("logger")->info("dP_C/P_C = " + to_string(x_0(2)) + ", dP_C = " + to_string(x_0(2)*P_C));
	spdlog::get("logger")->info("V_0/max_Size = " + to_string(x_0(3)) + ", V_0 = " + to_string(x_0(3)*max_Size));
	
	spdlog::get("logger")->info("(Q_dot_a - Q_dot_a_calc)/(Q_dot_a + Q_dot_a_calc) = " + to_string(f_0(0)) + ", Q_dot_a - Q_dot_a_calc = " + to_string(f_0(0)*(Q_dot_a + Q_dot_a_calc)));
	spdlog::get("logger")->info("(dP_H - dP_H_calc)/P_H_in = " + to_string(f_0(1)) + ", dP_H - dP_H_calc = " + to_string(f_0(1)*P_H_in));
	spdlog::get("logger")->info("(dP_C - dP_C_calc)/P_C = " + to_string(f_0(2)) + ", dP_C - dP_C_calc = " + to_string(f_0(2)*P_C));
	spdlog::get("logger")->info("(UA - targetParameter) /  (targetParameter + UA) = " + to_string(f_0(3)) + ", UA - targetParameter = " + to_string(f_0(3)*((targetParameter + UA))));
	*/

	double eps = 1e-5;
	int iterations = 0;
	bool isInvertable = false;
	//double* distance = new double[1000];
	while (fabs(f_0(0)) > eps || fabs(f_0(1)) > eps || fabs(f_0(2)) > eps || fabs(f_0(3)) > eps) {

		//distance[iterations] = pow(f_0(0),2) + pow(f_0(1), 2) + pow(f_0(2), 2) + pow(f_0(3), 2);

		/*spdlog::get("logger")->info("\n");
		spdlog::get("logger")->info("Jn");
		for (int i = 0; i < n; i++) {
			std::stringstream ss;
			for (int j = 0; j < n; j++) {
				ss << Jn(i, j) << " ";
			}
			spdlog::get("logger")->info(ss.str());
		}
		spdlog::get("logger")->info("\n");*/

		isInvertable = Jn.fullPivLu().isInvertible();

		if (!isInvertable) {
			//spdlog::get("logger")->info("Non-reversable :( \n");
			//spdlog::get("logger")->info("-");
			return 1;
		}
		else {
			Jn_inv = Jn.inverse();
		}
	
		/*spdlog::get("logger")->info("Jn_inv");
		for (int i = 0; i < n; i++) {
			std::stringstream ss;
			for (int j = 0; j < n; j++) {
				ss << Jn_inv(i, j) << " ";
			}
			spdlog::get("logger")->info(ss.str());
		}
		spdlog::get("logger")->info("\n");*/

		x_1 = x_0 - (Jn_inv * f_0);
		del_x = x_1 - x_0;

		
		if (x_1(0) <= T_C_in/T_H_in || x_1(0) >= 1 || x_1(1) < 0 || x_1(1) >= 1 || x_1(2) < 0 || x_1(3) < 0 || x_1(3) >= 1) {
			//spdlog::get("logger")->info("Refreshing jacobian");
			jacobian(n, x_0, dx, Jn_tmp);
			for (int i = 0; i < n; i++) {
				if (x_1(i) < 0 || x_1(i) >= 1) {
					Jn(i, i) = Jn_tmp(i, i);
				}
			}

			if (x_1(0) <= T_C_in / T_H_in) {
				Jn(0, 0) = Jn_tmp(0, 0);
			}

			/*spdlog::get("logger")->info("\n");
			spdlog::get("logger")->info("Jn");
			for (int i = 0; i < n; i++) {
				std::stringstream ss;
				for (int j = 0; j < n; j++) {
					ss << Jn(i, j) << " ";
				}
				spdlog::get("logger")->info(ss.str());
			}
			spdlog::get("logger")->info("\n");*/

			isInvertable = Jn.fullPivLu().isInvertible();

			if (!isInvertable) {
				//spdlog::get("logger")->info("Non-reversable :( \n");
				return 1;
			}
			else {
				Jn_inv = Jn.inverse();
			}

			/*spdlog::get("logger")->info("Jn_inv");
			for (int i = 0; i < n; i++) {
				std::stringstream ss;
				for (int j = 0; j < n; j++) {
					ss << Jn_inv(i, j) << " ";
				}
				spdlog::get("logger")->info(ss.str());
			}
			spdlog::get("logger")->info("\n");*/

			x_1 = x_0 - (Jn_inv * f_0);
			del_x = x_1 - x_0;
		}


		
		
		/*spdlog::get("logger")->info("T_H_out/T_H_in = " + to_string(x_1(0)) + ", T_H_out = " + to_string(x_1(0)*T_H_in));
		spdlog::get("logger")->info("AR = " + to_string(x_1(1)*max_Size) + ", L = " + to_string(L) + ", D_fr = " + to_string(D_fr));
		spdlog::get("logger")->info("dP_C/P_C = " + to_string(x_1(2)) + ", dP_C = " + to_string(x_1(2)*P_C));
		spdlog::get("logger")->info("V_0/max_Size = " + to_string(x_1(3)) + ", V_0 = " + to_string(x_1(3)*max_Size));
		*/

		evaluate(x_1, f_1);
		
		
		/*spdlog::get("logger")->info("(Q_dot_a - Q_dot_a_calc)/(Q_dot_a + Q_dot_a_calc) = " + to_string(f_1(0)) + ", Q_dot_a - Q_dot_a_calc = " + to_string(f_1(0)*(Q_dot_a + Q_dot_a_calc)));
		spdlog::get("logger")->info("(dP_H - dP_H_calc)/P_H_in = " + to_string(f_1(1)) + ", dP_H - dP_H_calc = " + to_string(f_1(1)*P_H_in));
		spdlog::get("logger")->info("(dP_C - dP_C_calc)/P_C = " + to_string(f_1(2)) + ", dP_C - dP_C_calc = " + to_string(f_1(2)*P_C));
		spdlog::get("logger")->info("(UA - targetParameter) / ( (targetParameter + UA)) = " + to_string(f_1(3)) + ", UA - targetParameter = " + to_string(f_1(3)*((targetParameter + UA))));
		*/
		del_f = f_1 - f_0;
		
		Jn = Jn + (del_f - Jn * del_x)*del_x.transpose() / del_x.squaredNorm();

		Jn(0, 1) = 0;
		Jn(0, 2) = 0;
		Jn(0, 3) = 0;

		Jn(1, 0) = 0;
		Jn(1, 2) = 0;
		//Jn(1, 3) = 0;

		Jn(2, 0) = 0;
		Jn(2, 1) = 0;
		Jn(2, 3) = 0;

		Jn(3, 0) = 0;
		Jn(3, 1) = 0;
		Jn(3, 2) = 0;

		x_0 = x_1;
		f_0 = f_1;

		iterations++;
		//spdlog::get("logger")->info(std::to_string(iterations) + "\n");
	}

	clock_t end = clock();
	double elapsed = double(end - begin) / CLOCKS_PER_SEC*1000;
	spdlog::get("logger")->info(to_string(T_H_out) + "," + to_string(L) + "," + to_string(D_fr) + "," + 
		to_string(V_0) + "," + to_string(AR) + "," + to_string(UA) + "," + to_string(dP_max) + "," + to_string(epsilon) + "," +
		to_string(Q_dot_a - Q_dot_a_calc) + "," + to_string(dP_H - dP_H_calc) + "," + to_string(dP_C - dP_C_calc)
		+ "," + to_string(targetParameter) + "," + to_string(targetdP_max_Regen) + "," + to_string(iterations)
		+ "," + to_string(elapsed));

	/*
	spdlog::get("logger")->info("Solved after " + std::to_string(iterations) + " iterations.\n\n");
	spdlog::get("logger")->info("T_H_out/T_H_in = " + to_string(x_0(0)) + ", T_H_out = " + to_string(x_0(0)*T_H_in));
	spdlog::get("logger")->info("AR = " + to_string(x_0(1)*max_Size) + ", L = " + to_string(L) + ", D_fr = " + to_string(D_fr));
	spdlog::get("logger")->info("dP_C/P_C = " + to_string(x_0(2)) + ", dP_C = " + to_string(x_0(2)*P_C));
	spdlog::get("logger")->info("V_0/max_Size = " + to_string(x_0(3)) + ", V_0 = " + to_string(x_0(3)*max_Size));
	
	spdlog::get("logger")->info("(Q_dot_a - Q_dot_a_calc)/(Q_dot_a + Q_dot_a_calc) = " + to_string(f_0(0)) + ", Q_dot_a - Q_dot_a_calc = " + to_string(f_0(0)*(Q_dot_a + Q_dot_a_calc)));
	spdlog::get("logger")->info("(dP_H - dP_H_calc)/P_H_in = " + to_string(f_0(1)) + ", dP_H - dP_H_calc = " + to_string(f_0(1)*P_H_in));
	spdlog::get("logger")->info("(dP_C - dP_C_calc)/P_C = " + to_string(f_0(2)) + ", dP_C - dP_C_calc = " + to_string(f_0(2)*P_C));
	spdlog::get("logger")->info("(UA - targetParameter) /  (targetParameter + UA) = " + to_string(f_0(3)) + ", UA - targetParameter = " + to_string(f_0(3)* (targetParameter + UA)));
	*/
	/*
	for (int i = 0; i < iterations; i++) {
		spdlog::get("logger")->info(to_string(distance[i]));
	}*/
}

int RegeneratorModel::getOffDesignSolution()
{
	calcCarryoverMassFlow();

	try {
		calculateModel();
	}
	catch (const invalid_argument& e) {
		return -1;
	}

	

	return 0;
}

void RegeneratorModel::getSolution(RegeneratorSolution * solution)
{
	solution->dP_C = dP_C;
	solution->dP_H = dP_H;
	solution->epsilon = epsilon;
	solution->T_H_in = T_H_in;
	solution->T_C_in = T_C_in;
	solution->T_C_out = T_C_out;
	solution->T_H_out = T_H_out;
	solution->NTU_R_e = NTU_R_e;
	solution->Q_dot_a = Q_dot_a;
	solution->UA = UA;
	solution->m_dot_H = m_dot_H;
	solution->m_dot_C = m_dot_C;
	solution->m_dot_carryover = m_dot_carryover;
	solution->m_HTR_HP_dP = dP_C;
	solution->m_HTR_LP_dP = dP_H;
	solution->costModule = costModule;
	solution->L = L;
	solution->D_fr = D_fr;
	solution->valves = valves;
	solution->wallThickness = wallThickness;
}
