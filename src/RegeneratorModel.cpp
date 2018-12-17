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

	if (spdlog::get("logger") == nullptr) {
		auto logger = spdlog::basic_logger_mt("logger", LOG_FILEPATH);
		spdlog::flush_on(spdlog::level::info);
		spdlog::drop("logger");
		spdlog::register_logger(logger);
		spdlog::set_pattern("%v");
	}
}

RegeneratorModel::~RegeneratorModel()
{
	spdlog::drop("logger");

	delete bedMaterialTable;
	delete regeneratorTable;
	delete spheresRPTable;

	delete HeatTransfer;
	delete HotPressureDrop;
	delete ColdPressureDrop;
	delete Length;
	delete Diameter;
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
	P_C_out = P_C_in - dP_C;	//[kPa]

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

	A_fr = PI * pow(D_fr, 2) / 4;
	T_H_f = (T_H_in + T_H_out) / 2;
	T_C_f = (T_C_in + T_C_out) / 2;

	try {
		packedspheresFitCO2(m_dot_H, D_s, A_fr, L, T_H_f, P_H_in, e_v, &f_H, &h_H, &NTU_H, &dP_H_calc);
		packedspheresFitCO2(m_dot_C, D_s, A_fr, L, T_C_f, P_C_in, e_v, &f_C, &h_C, &NTU_C, &dP_C_calc);
	}
	catch (const invalid_argument& e) {
		throw e;
	}

	V_0 = A_fr * L;

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

	Q_dot_calc = epsilon * Q_dot_max;
	Q_dot_a_calc = Q_dot_calc - Q_dot_loss;
	dP_max = max(dP_C, dP_H);
}

void RegeneratorModel::loadTables() {
	bedMaterialTable = new LookupTable_1D(PROPERTY_FILES + bedMaterialName + ".csv");
	regeneratorTable = new LookupTable_2D(BALANCED_REGENERATOR_TABLE_PATH);
	spheresRPTable = new LookupTable_1D(SPHERES_RP_TABLE_PATH);
}

void RegeneratorModel::setInletStates(double T_H_in, double P_H_in, double m_dot_H, double T_C_in, double P_C_in, double m_dot_C)
{
	this->T_H_in = T_H_in;
	this->P_H_in = P_H_in;
	this->m_dot_H = m_dot_H;
	this->T_C_in = T_C_in;
	this->P_C_in = P_C_in;
	this->m_dot_C = m_dot_C;
	m_dot_carryover = 0;

	if (P_C_in <= N_co2_props::P_lower_limit || P_H_in <= N_co2_props::P_lower_limit || P_C_in >= N_co2_props::P_upper_limit || P_H_in >= N_co2_props::P_upper_limit) {
		throw invalid_argument("Inlet pressures are either below " + std::to_string(N_co2_props::P_lower_limit) + "[kPa] or above " + std::to_string(N_co2_props::P_upper_limit) + "[kPa]!");
	}
	if (T_H_in <= N_co2_props::T_lower_limit || T_C_in <= N_co2_props::T_lower_limit || T_H_in >= N_co2_props::T_upper_limit || T_C_in >= N_co2_props::T_upper_limit) {
		throw invalid_argument("Inlet temperatures are either below " + std::to_string(N_co2_props::T_lower_limit) + "[K] or above " + std::to_string(N_co2_props::T_upper_limit) + "[K]!");
	}

	//Enthalpy at hot inlet
	CO2_TP(T_H_in, P_H_in, &CO2State);
	h_H_in = CO2State.enth;

	//Enthalpy at cold inlet
	CO2_TP(T_C_in, P_C_in, &CO2State);
	h_C_in = CO2State.enth;

	CO2_TP(T_C_in, P_H_in, &CO2State);
	h_H_out_max_p = CO2State.enth;

	CO2_TP(T_H_in, P_C_in, &CO2State);
	h_C_out_max_p = CO2State.enth;

	C_p_H = (h_H_in - h_H_out_max_p) / (T_H_in - T_C_in);
	C_p_C = (h_C_out_max_p - h_C_in) / (T_H_in - T_C_in);

	T_f = (T_H_in + T_C_in) / 2;
	rho_s = bedMaterialTable->getValue("rho", "T", T_f);
	c_s = bedMaterialTable->getValue("c", "T", T_f);

	massflowVariablesInit();
}

void RegeneratorModel::setParameters(valveDesignOption::valveModes valveMode, double Q_dot_loss, double P_0, double D_s, double e_v)
{
	this->Q_dot_loss = Q_dot_loss;
	this->D_s = D_s;
	this->e_v = e_v;
	this->P_0 = P_0;

	numberOfCycles = 1.2 * operationYears * 365 * operationHoursPerDay * 60 * 60 / P_0;
}

void RegeneratorModel::setDesignTargets(targetModes::targetModes targetMode, targetModes::target2Modes secondTargetMode, double targetParameter, double secondTargetParameter)
{
	this->secondTargetMode = secondTargetMode;
	this->targetdP_max_Regen = secondTargetParameter;
	
	
	this->targetMode = targetMode;
	this->targetParameter = targetParameter;
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

	*targetParameter = UA;

	return 0;
}

int RegeneratorModel::getDesignSolution()
{
	clock_t begin = clock();
	HeatTransfer_SP.solverName = "Balance Heat Tansfer";
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
	ColdPressureDrop_SP.lowerBound = 0.1;							ColdPressureDrop_SP.upperBound = P_C_in;
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

	Diameter->setParameters(&Diameter_SP);
	
	int statusSolver;
	statusSolver = Diameter->solve();

	if (statusSolver != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	clock_t end = clock();
	double elapsed = double(end - begin) / CLOCKS_PER_SEC * 1000;
	cout << T_H_out << ", " << L << ", " << D_fr << ", " << V_0 << ", " << AR << ", " << UA << ", " <<
		dP_max << ", " << epsilon << ", " << (Q_dot_a - Q_dot_a_calc) << ", " << (dP_H - dP_H_calc) << ", " << (dP_C - dP_C_calc)
		<< ", " << targetParameter << ", " << targetdP_max_Regen << ", "  << ", " << elapsed << endl;

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
	solution->L = L;
	solution->D_fr = D_fr;
}
