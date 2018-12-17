#include "RegeneratorModel.h"

std::string SSCDIR1(std::getenv("SSCDIR"));

string const RegeneratorModel::LOG_FILEPATH = "C:/Users/Dmitrii/Documents/ME459/me459-project/RegenHX_LOG.log";
string const RegeneratorModel::PROPERTY_FILES = "C:/Users/Dmitrii/Documents/ME459/me459-project/PropertyFiles/";
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
	this->dP_H = params(1)*P_H_in;
	this->dP_C = params(2)*P_C;
	this->V_0 = params(3)*max_Size;
	//this->AR = params(4)*max_Size;
	//this->L = params(4)*max_Size;
	//this->wallThickness = params(5);
	//this->m_dot_carryover = params(6);
}

void RegeneratorModel::jacobian_n(int n, Eigen::VectorXd& x_nm1, Eigen::VectorXd& x_n, Eigen::VectorXd& f_nm1, Eigen::VectorXd& f_n, Eigen::MatrixXd& Jn) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Jn(i, j) = (f_n(i) - f_nm1(i)) / (x_n(j) - x_nm1(j));
		}
	}

	Jn(0, 1) = 0;
	Jn(0, 2) = 0;
	Jn(0, 3) = 0;
	Jn(0, 4) = 0;

	Jn(1, 0) = 0;
	Jn(1, 2) = 0;
	Jn(1, 3) = 0;
	Jn(1, 4) = 0;

	Jn(2, 0) = 0;
	Jn(2, 1) = 0;
	Jn(2, 3) = 0;
	Jn(2, 4) = 0;

	Jn(3, 0) = 0;
	Jn(3, 1) = 0;
	Jn(3, 2) = 0;

	Jn(4, 0) = 0;
	Jn(4, 1) = 0;
	Jn(4, 2) = 0;
}

void RegeneratorModel::jacobian_first_time(int n, Eigen::VectorXd x_n, Eigen::VectorXd dx, Eigen::MatrixXd& Jn) {
	Eigen::VectorXd f_i(n);
	Eigen::VectorXd f_im1(n);
	Eigen::VectorXd f_ip1(n);
	Eigen::VectorXd x_tmp(n);
	evaluate(x_n, f_i);

	x_tmp = x_n;

	// Through variables
	for (int j = 0; j < n; j++) {
		x_tmp(j) -= dx(j);
		evaluate(x_tmp, f_im1);
		x_tmp(j) += 2*dx(j);
		evaluate(x_tmp, f_ip1);

		// Through functions
		for (int i = 0; i < n; i++) {
			Jn(i, j) = -(f_im1(i) + f_ip1(i) - 2 * f_i(i)) / pow(dx(j), 2);
		}
	}

	Jn(0, 1) = 0;
	Jn(0, 2) = 0;
	Jn(0, 3) = 0;
	//Jn(0, 4) = 0;

	Jn(1, 0) = 0;
	Jn(1, 2) = 0;
	Jn(1, 3) = 0;
	//Jn(1, 4) = 0;

	Jn(2, 0) = 0;
	Jn(2, 1) = 0;
	Jn(2, 3) = 0;
	//Jn(2, 4) = 0;

	Jn(3, 0) = 0;
	Jn(3, 1) = 0;
	Jn(3, 2) = 0;
	//Jn(3, 4) = 0;

	//Jn(4, 0) = 0;
	//Jn(4, 1) = 0;
	//Jn(4, 2) = 0;
	//Jn(4, 3) = 0;
}

void RegeneratorModel::evaluate(Eigen::VectorXd point, Eigen::VectorXd& f) {
	setSolver(point);

	//this->m_dot_C -= m_dot_carryover_guess;
	//this->m_dot_H -= m_dot_carryover_guess;

	//massflowVariablesInit();
	calculateModel();
	//calculateCost();

	//this->m_dot_C += m_dot_carryover;
	//this->m_dot_H += m_dot_carryover;

	//calcCarryoverMassFlow();

	f(0) = (Q_dot_a - Q_dot_a_calc)/(Q_dot_a + Q_dot_a_calc);
	f(1) = (dP_H - dP_H_calc)/P_H_in;
	f(2) = (dP_C - dP_C_calc)/P_C;
	f(3) = (UA - targetParameter)/ (targetParameter + UA);
	//f(4) = (dP_max - targetdP_max_Regen)/targetdP_max_Regen;
	//values(5) = stressAmplitude_calc - targets(5);
	//values(6) = m_dot_carryover_guess - m_dot_carryover;
}

int RegeneratorModel::getDesignSolution()
{
	spdlog::get("logger")->info("Here1");
	int n = 4;

	Eigen::VectorXd x_0(n);
	Eigen::VectorXd x_1(n);
	/*parameters(0) = 490;
	parameters(1) = 215;
	parameters(2) = 60;
	parameters(3) = 1.2;
	parameters(4) = 1.1;*/

	//T_H_out = 450;
	//dP_C = 60;
	//dP_H = 225;
	//D_fr = 1.2;
	//L = 1.1;
	AR = 1.1;

	spdlog::get("logger")->info("Here2");
	x_0(0) = 600/T_H_in;//(T_C_in + T_H_in) / 2;
	x_0(1) = 0.8*targetdP_max_Regen/P_H_in;
	x_0(2) = 0.3*targetdP_max_Regen/P_C;
	x_0(3) = 2.5/max_Size;
	//x_0(4) = 0.7/max_Size;
	spdlog::get("logger")->info("Here3");

	//Eigen::VectorXd targets(n);
	//targets(0) = 0;
	//targets(1) = 0;
	//targets(2) = 0;
	//targets(3) = targetParameter;
	//targets(4) = targetdP_max_Regen;
	//spdlog::get("logger")->info(targetParameter);
	//spdlog::get("logger")->info(to_string(targetdP_max_Regen) + "\n");
	//targets(5) = stressAmplitude;
	//targets(6) = 0;

	Eigen::VectorXd dx(n);
	dx(0) = 1/T_H_in;
	dx(1) = 1/P_H_in;
	dx(2) = 1/P_C;
	dx(3) = 0.1/max_Size;
	//dx(4) = 0.1/max_Size;
	Eigen::MatrixXd Jn(n, n);
	spdlog::get("logger")->info("Here4");
	jacobian_first_time(n, x_0, dx, Jn);
	spdlog::get("logger")->info("Here5");

	//Eigen::VectorXd damping(n);
	//damping(0) = 1;
	//damping(1) = 0.5;
	//damping(2) = 0.5;
	//damping(3) = 1;
	//damping(4) = 1;
	
	Eigen::MatrixXd Jn_inv(n, n);
	Eigen::VectorXd f_0(n);
	Eigen::VectorXd f_1(n);
	Eigen::VectorXd del_f(n);
	Eigen::VectorXd del_x(n);
	evaluate(x_0, f_0);
	

	spdlog::get("logger")->info("T_H_out/T_H_in = " + to_string(x_0(0)) + ", T_H_out = " + to_string(x_0(0)*T_H_in));
	spdlog::get("logger")->info("dP_H/P_H_in = " + to_string(x_0(1)) + ", dP_H = " + to_string(x_0(1)*P_H_in));
	spdlog::get("logger")->info("dP_C/P_C = " + to_string(x_0(2)) + ", dP_C = " + to_string(x_0(2)*P_C));
	spdlog::get("logger")->info("V_0/max_Size = " + to_string(x_0(3)) + ", V_0 = " + to_string(x_0(3)*max_Size));
	//spdlog::get("logger")->info("AR = " + to_string(x_0(4)*max_Size) + ", L = " + to_string(L) + ", D_fr = " + to_string(D_fr));
	spdlog::get("logger")->info("(Q_dot_a - Q_dot_a_calc)/(Q_dot_a + Q_dot_a_calc) = " + to_string(f_0(0)) + ", Q_dot_a - Q_dot_a_calc = " + to_string(f_0(0)*(Q_dot_a + Q_dot_a_calc)));
	spdlog::get("logger")->info("(dP_H - dP_H_calc)/P_H_in = " + to_string(f_0(1)) + ", dP_H - dP_H_calc = " + to_string(f_0(1)*P_H_in));
	spdlog::get("logger")->info("(dP_C - dP_C_calc)/P_C = " + to_string(f_0(2)) + ", dP_C - dP_C_calc = " + to_string(f_0(2)*P_C));
	spdlog::get("logger")->info("(UA - targetParameter) /  (targetParameter + UA) = " + to_string(f_0(3)) + ", UA - targetParameter = " + to_string(f_0(3)*((targetParameter + UA))));
	//spdlog::get("logger")->info("(dP_max - targetdP_max_Regen) / targetdP_max_Regen = " + to_string(f_0(4)) + ", dP_max - targetdP_max_Regen = " + to_string(f_0(4)*targetdP_max_Regen));

	double eps = 1e-5;
	int iterations = 0;
	bool isInvertable = false;
	//double distance = std::numeric_limits<double>::max();
	//double best_distance = std::numeric_limits<double>::max();
	double* distance = new double[1000];
	while (fabs(f_0(0)) > eps || fabs(f_0(1)) > eps || fabs(f_0(2)) > eps || fabs(f_0(3)) > eps) {

		distance[iterations] = fabs(f_0(0)) + fabs(f_0(1)) + fabs(f_0(2)) + fabs(f_0(3));

		spdlog::get("logger")->info("\n");
		spdlog::get("logger")->info("Jn");
		for (int i = 0; i < n; i++) {
			std::stringstream ss;
			for (int j = 0; j < n; j++) {
				ss << Jn(i, j) << " ";
			}
			spdlog::get("logger")->info(ss.str());
		}
		spdlog::get("logger")->info("\n");

		isInvertable = Jn.fullPivLu().isInvertible();

		if (!isInvertable) {
			spdlog::get("logger")->info("Non-reversable :( \n");
			return 0;
		}
		else {
			Jn_inv = Jn.inverse();
			//Jn_inv.normalize();
		}
	
		spdlog::get("logger")->info("Jn_inv");
		for (int i = 0; i < n; i++) {
			std::stringstream ss;
			for (int j = 0; j < n; j++) {
				ss << Jn_inv(i, j) << " ";
			}
			spdlog::get("logger")->info(ss.str());
		}
		spdlog::get("logger")->info("\n");

		x_1 = x_0 - (Jn_inv * f_0);
		del_x = x_1 - x_0;
		/*if (x_1(3) < 0) {
			x_1(3) = 0.1/max_Size;
			del_x = x_1 - x_0;
		}*/
		/*for (int i = 0; i < n; i++) {
			if (fabs(del_x(i) / x_0(i)) > 0.3) {
				del_x(i) = copysign(x_0(i)*0.3, del_x(i));
			}
		}*/

		

		spdlog::get("logger")->info("T_H_out/T_H_in = " + to_string(x_1(0)) + ", T_H_out = " + to_string(x_1(0)*T_H_in));
		spdlog::get("logger")->info("dP_H/P_H_in = " + to_string(x_1(1)) + ", dP_H = " + to_string(x_1(1)*P_H_in));
		spdlog::get("logger")->info("dP_C/P_C = " + to_string(x_1(2)) + ", dP_C = " + to_string(x_1(2)*P_C));
		spdlog::get("logger")->info("V_0/max_Size = " + to_string(x_1(3)) + ", V_0 = " + to_string(x_1(3)*max_Size));
		//spdlog::get("logger")->info("AR = " + to_string(x_1(4)*max_Size) + ", L = " + to_string(L) + ", D_fr = " + to_string(D_fr));

		evaluate(x_1, f_1);
		
		spdlog::get("logger")->info("(Q_dot_a - Q_dot_a_calc)/(Q_dot_a + Q_dot_a_calc) = " + to_string(f_1(0)) + ", Q_dot_a - Q_dot_a_calc = " + to_string(f_1(0)*(Q_dot_a + Q_dot_a_calc)));
		spdlog::get("logger")->info("(dP_H - dP_H_calc)/P_H_in = " + to_string(f_1(1)) + ", dP_H - dP_H_calc = " + to_string(f_1(1)*P_H_in));
		spdlog::get("logger")->info("(dP_C - dP_C_calc)/P_C = " + to_string(f_1(2)) + ", dP_C - dP_C_calc = " + to_string(f_1(2)*P_C));
		spdlog::get("logger")->info("(UA - targetParameter) / ( (targetParameter + UA)) = " + to_string(f_1(3)) + ", UA - targetParameter = " + to_string(f_1(3)*((targetParameter + UA))));
		//spdlog::get("logger")->info("(dP_max - targetdP_max_Regen) / targetdP_max_Regen = " + to_string(f_1(4)) + ", dP_max - targetdP_max_Regen = " + to_string(f_1(4)*targetdP_max_Regen));

		del_f = f_1 - f_0;
		
		Jn = Jn + (del_f - Jn * del_x)*del_x.transpose() / del_x.squaredNorm();

		Jn(0, 1) = 0;
		Jn(0, 2) = 0;
		//Jn(0, 3) = 0;
		//Jn(0, 4) = 0;

		Jn(1, 0) = 0;
		Jn(1, 2) = 0;
		Jn(1, 3) = 0;
		//Jn(1, 4) = 0;

		Jn(2, 0) = 0;
		Jn(2, 1) = 0;
		Jn(2, 3) = 0;
		//Jn(2, 4) = 0;

		Jn(3, 0) = 0;
		Jn(3, 1) = 0;
		Jn(3, 2) = 0;
		//Jn(3, 4) = 0;

		/*Jn(4, 0) = 0;
		Jn(4, 1) = 0;
		Jn(4, 2) = 0;
		Jn(4, 3) = 0;*/

		x_0 = x_1;
		f_0 = f_1;

		
		

		

		/*spdlog::get("logger")->info("Step size: " + std::to_string(steps(0)));
		for (int i = 0; i < n; i++) {
			if (parameters(i) < lowerBounds(i)) {
				parameters(i) = lowerBounds(i);
			}
		}

		std::stringstream ss1, ss2, ss3, ss4, ss5;
		ss1 << parameters;
		spdlog::get("logger")->info("Parameters");
		spdlog::get("logger")->info(ss1.str());

		differences(n, parameters, targets, steps, diff);
		evaluate(parameters, targets, values);
		Solver::jacobian(n, diff, steps, jacobian);

		ss2 << values;
		ss3 << jacobian;
		ss4 << diff;

		spdlog::get("logger")->info("Values");
		spdlog::get("logger")->info(ss2.str());

		spdlog::get("logger")->info("Jacobian");
		spdlog::get("logger")->info(ss3.str());

		distance = sqrt(values(0)*values(0) + values(1)*values(1) + values(2)*values(2) + values(3)*values(3) + values(4)*values(4) + values(5)*values(5) + values(6)*values(6));

		if (distance >= best_distance) {
			for (int i = 0; i < n; i++) {
				steps(i) /= 2;
			}
			continue;
		}
		else {
			best_distance = distance;
		}

		spdlog::get("logger")->info("Values\n");
		spdlog::get("logger")->info(ss2.str());
		spdlog::get("logger")->info("\n\n");

		spdlog::get("logger")->info("Distance: " + std::to_string(distance));

			spdlog::get("logger")->info("Differences\n");
			spdlog::get("logger")->info(ss4.str());
			spdlog::get("logger")->info("\n\n");

			spdlog::get("logger")->info("Jacobian\n");
			spdlog::get("logger")->info(ss3.str());
			spdlog::get("logger")->info("\n\n");

			//jacobian.computeInverseWithCheck(jacobian_inverse, isInvertable);
		isInvertable = jacobian.fullPivLu().isInvertible();
		//	spdlog::get("logger")->info("Jacobian is invertable: " + std::to_string(isInvertable) + "\n");

		if (!isInvertable) {
			return 0;
		}
		else {
			jacobian_inverse = jacobian.inverse();
		}

		parameters = parameters - limiter * jacobian_inverse * values;*/
		iterations++;
		spdlog::get("logger")->info(std::to_string(iterations) + "\n");
	}

	spdlog::get("logger")->info("Solved after " + std::to_string(iterations) + " iterations.\n\n");
	spdlog::get("logger")->info("T_H_out/T_H_in = " + to_string(x_0(0)) + ", T_H_out = " + to_string(x_0(0)*T_H_in));
	spdlog::get("logger")->info("dP_H/P_H_in = " + to_string(x_0(1)) + ", dP_H = " + to_string(x_0(1)*P_H_in));
	spdlog::get("logger")->info("dP_C/P_C = " + to_string(x_0(2)) + ", dP_C = " + to_string(x_0(2)*P_C));
	spdlog::get("logger")->info("V_0/max_Size = " + to_string(x_0(3)) + ", V_0 = " + to_string(x_0(3)*max_Size));
	//spdlog::get("logger")->info("AR = " + to_string(x_0(4)*max_Size) + ", L = " + to_string(L) + ", D_fr = " + to_string(D_fr));
	spdlog::get("logger")->info("(Q_dot_a - Q_dot_a_calc)/(Q_dot_a + Q_dot_a_calc) = " + to_string(f_0(0)) + ", Q_dot_a - Q_dot_a_calc = " + to_string(f_0(0)*(Q_dot_a + Q_dot_a_calc)));
	spdlog::get("logger")->info("(dP_H - dP_H_calc)/P_H_in = " + to_string(f_0(1)) + ", dP_H - dP_H_calc = " + to_string(f_0(1)*P_H_in));
	spdlog::get("logger")->info("(dP_C - dP_C_calc)/P_C = " + to_string(f_0(2)) + ", dP_C - dP_C_calc = " + to_string(f_0(2)*P_C));
	spdlog::get("logger")->info("(UA - targetParameter) /  (targetParameter + UA) = " + to_string(f_0(3)) + ", UA - targetParameter = " + to_string(f_0(3)* (targetParameter + UA)));
	//spdlog::get("logger")->info("(dP_max - targetdP_max_Regen) / targetdP_max_Regen = " + to_string(f_0(4)) + ", dP_max - targetdP_max_Regen = " + to_string(f_0(4)*targetdP_max_Regen));

	for (int i = 0; i < iterations; i++) {
		spdlog::get("logger")->info(to_string(distance[i]));
	}

	/*std::stringstream ss1, ss2;
	ss1 << x_n;
	ss2 << f_n;
	spdlog::get("logger")->info(ss1.str());
	spdlog::get("logger")->info("\n\n");
	spdlog::get("logger")->info(ss2.str());
	spdlog::get("logger")->info("\n\n");*/
}
	/*HeatTransfer_SP.solverName = "Balance Heat Tansfer";
	HeatTransfer_SP.target = 0;
	HeatTransfer_SP.guessValue1 = T_C_in;						HeatTransfer_SP.guessValue2 = (T_C_in + T_H_in) / 2;
	HeatTransfer_SP.lowerBound = N_co2_props::T_lower_limit;	HeatTransfer_SP.upperBound = N_co2_props::T_upper_limit;
	HeatTransfer_SP.tolerance = 0.1;
	HeatTransfer_SP.iterationLimit = 50;
	HeatTransfer_SP.isErrorRel = false;
	HeatTransfer_SP.classInst = this;
	HeatTransfer_SP.monoEquation = &RegeneratorModel::HeatTransfer_ME;
	HeatTransfer->setParameters(&HeatTransfer_SP);

	if (secondTargetMode == targetModes::AR) {
		targetdP_max_Regen = 187.5;
	}
	HotPressureDrop_SP.solverName = "Balance Hot Pressure Drop";
	HotPressureDrop_SP.target = 0;
	HotPressureDrop_SP.guessValue1 = 0.8*targetdP_max_Regen;		HotPressureDrop_SP.guessValue2 = targetdP_max_Regen;
	HotPressureDrop_SP.lowerBound = 0.1;					HotPressureDrop_SP.upperBound = P_H_in;
	HotPressureDrop_SP.tolerance = 0.0001;
	HotPressureDrop_SP.iterationLimit = 50;
	HotPressureDrop_SP.isErrorRel = false;
	HotPressureDrop_SP.classInst = this;
	HotPressureDrop_SP.monoEquation = &RegeneratorModel::HotPressureDrop_ME;
	HotPressureDrop->setParameters(&HotPressureDrop_SP);

	ColdPressureDrop_SP.solverName = "Balance Cold Pressure Drop";
	ColdPressureDrop_SP.target = 0;
	ColdPressureDrop_SP.guessValue1 = 0.3*targetdP_max_Regen - 1;		ColdPressureDrop_SP.guessValue2 = 0.3*targetdP_max_Regen;
	ColdPressureDrop_SP.lowerBound = 0.1;							ColdPressureDrop_SP.upperBound = P_C;
	ColdPressureDrop_SP.tolerance = 0.0001;
	ColdPressureDrop_SP.iterationLimit = 50;
	ColdPressureDrop_SP.isErrorRel = false;
	ColdPressureDrop_SP.classInst = this;
	ColdPressureDrop_SP.monoEquation = &RegeneratorModel::ColdPressureDrop_ME;
	ColdPressureDrop->setParameters(&ColdPressureDrop_SP);
	
	double m_dot_average = (m_dot_C + m_dot_H) / 2;

	Diameter_SP.target = targetParameter;
	Diameter_SP.guessValue1 = 0.8 * m_dot_average / 37;		Diameter_SP.guessValue2 = 1.1 * m_dot_average / 37;
	Diameter_SP.lowerBound = 0.1;		Diameter_SP.upperBound = 10;
	Diameter_SP.tolerance = 0.0001;
	Diameter_SP.iterationLimit = 50;
	Diameter_SP.isErrorRel = true;
	Diameter_SP.classInst = this;

	CarryoverMassFlow_SP.target = 0;
	CarryoverMassFlow_SP.lowerBound = 0;		CarryoverMassFlow_SP.upperBound = min(m_dot_C, m_dot_H);
	CarryoverMassFlow_SP.tolerance = 0.001;
	CarryoverMassFlow_SP.iterationLimit = 50;
	CarryoverMassFlow_SP.isErrorRel = false;
	CarryoverMassFlow_SP.classInst = this;

	if (secondTargetMode == targetModes::dP_max) {
		Length_SP.solverName = "Length Solver";
		Length_SP.target = targetdP_max_Regen;
		Length_SP.guessValue1 = 0.8;		Length_SP.guessValue2 = 1.1;
		Length_SP.lowerBound = 0.1;			Length_SP.upperBound = 10;
		Length_SP.tolerance = 0.0001;
		Length_SP.iterationLimit = 50;
		Length_SP.isErrorRel = true;
		Length_SP.classInst = this;
		Length_SP.monoEquation = &RegeneratorModel::Length_ME;
		Length->setParameters(&Length_SP);

		Diameter_SP.solverName = "Diameter_dP Solver";
		Diameter_SP.monoEquation = &RegeneratorModel::Diameter_dP_ME;
	}
	else if (secondTargetMode == targetModes::AR) {
		Diameter_SP.solverName = "Diameter_AR Solver";
		Diameter_SP.monoEquation = &RegeneratorModel::Diameter_AR_ME;
	}

	Diameter->setParameters(&Diameter_SP);

	if (valveMode == valveDesignOption::valveModes::FIXED_CV) {
		PressureSplit_SP.solverName = "Splits Pressure Drop Between Regenerator and Valves";
		PressureSplit_SP.target = targetdP_max_Regen;
		PressureSplit_SP.guessValue1 = targetdP_max_Regen * 0.8;		PressureSplit_SP.guessValue2 = targetdP_max_Regen * 0.9;
		PressureSplit_SP.lowerBound = 0.1;					PressureSplit_SP.upperBound = P_H_in;
		PressureSplit_SP.tolerance = 0.0001;
		PressureSplit_SP.iterationLimit = 50;
		PressureSplit_SP.isErrorRel = true;
		PressureSplit_SP.classInst = this;
		PressureSplit_SP.monoEquation = &RegeneratorModel::PressureSplit_ME;
		PressureSplit->setParameters(&PressureSplit_SP);

		CarryoverMassFlow_SP.solverName = "Carryover Mass Fixed CV Solver";
		CarryoverMassFlow_SP.monoEquation = &RegeneratorModel::CarryoverMassFlow_FIXED_CV_ME;
	}
	else if (valveMode == valveDesignOption::valveModes::FIXED_DP) {
		CarryoverMassFlow_SP.solverName = "Carryover Mass Fixed DP Solver";
		CarryoverMassFlow_SP.monoEquation = &RegeneratorModel::CarryoverMassFlow_FIXED_DP_ME;
	}

	WallThickness_SP.solverName = "Wall Thickness";
	WallThickness_SP.target = stressAmplitude;
	WallThickness_SP.guessValue1 = 0.04;		WallThickness_SP.guessValue2 = 0.05;
	WallThickness_SP.lowerBound = 0;			WallThickness_SP.upperBound = 1;
	WallThickness_SP.tolerance = 0.001;
	WallThickness_SP.iterationLimit = 50;
	WallThickness_SP.isErrorRel = true;
	WallThickness_SP.classInst = this;
	WallThickness_SP.monoEquation = &RegeneratorModel::WallThickness_ME;
	WallThickness->setParameters(&WallThickness_SP);
	
	int statusSolver;
	statusSolver = Diameter->solve();

	if (statusSolver != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}
	
	calcCarryoverMassFlow();

	CarryoverMassFlow_SP.guessValue1 = m_dot_carryover;
	CarryoverMassFlow_SP.guessValue2 = m_dot_carryover + 1;
	CarryoverMassFlow->setParameters(&CarryoverMassFlow_SP);

	if (D_fr - 0.001 < Diameter_SP.lowerBound) {
		Diameter->updateGuesses(D_fr, D_fr + 0.001);
	}
	else {
		Diameter->updateGuesses(D_fr - 0.001, D_fr);
	}

	if (secondTargetMode == targetModes::dP_max) {
		if (L - 0.001 < Length_SP.lowerBound) {
			Length->updateGuesses(L, L + 0.001);
		}
		else {
			Length->updateGuesses(L - 0.001, L);
		}
	}

	if (T_H_out - 1 < HeatTransfer_SP.lowerBound) {
		HeatTransfer->updateGuesses(T_H_out, T_H_out + 1);
	}
	else {
		HeatTransfer->updateGuesses(T_H_out - 1, T_H_out);
	}

	if (dP_H - 10 < HotPressureDrop_SP.lowerBound) {
		HotPressureDrop->updateGuesses(dP_H, dP_H + 10);
	}
	else {
		HotPressureDrop->updateGuesses(dP_H - 10, dP_H);
	}

	if (dP_C - 10 < ColdPressureDrop_SP.lowerBound) {
		ColdPressureDrop->updateGuesses(dP_C, dP_C+10);
	}
	else {
		ColdPressureDrop->updateGuesses(dP_C - 10, dP_C);
	}

	double tolerance;
	statusSolver = CarryoverMassFlow->solve(&tolerance);
	
	if (statusSolver != C_monotonic_eq_solver::CONVERGED) {
		if (isfinite(tolerance) == false) {
			return -1;
		}

		if (tolerance / min(m_dot_C, m_dot_H) > 0.01)
		{
			return -1;
		}
	}

	if (valveMode == valveDesignOption::valveModes::FIXED_CV) {
		dP_C = dP_C_total;
		dP_H = dP_H_total;
		calculateValveCvs();
	}
	
	carryoverEnthDrop();
	if (calculateCost() != 0) {
		return -1;
	}

	m_dot_H -= m_dot_carryover;
	m_dot_C -= m_dot_carryover;

	return 0;
}*/

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
