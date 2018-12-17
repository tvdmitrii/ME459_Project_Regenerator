#include "RegeneratorModel.h"

string const RegeneratorModel::LOG_FILEPATH = string(BUILD_PATH) + "/RegenHX_LOG.log";
string const RegeneratorModel::PROPERTY_FILES = string(BUILD_PATH) + "/PropertyFiles/";
string const RegeneratorModel::SPHERES_RP_TABLE_PATH = PROPERTY_FILES + "Spheres_RP.csv";
string const RegeneratorModel::BALANCED_REGENERATOR_TABLE_PATH = PROPERTY_FILES + "balanced-regenerator.csv";


RegeneratorModel::RegeneratorModel()
{
	loadTables();

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

	D_fr = cbrt(4 * V_0 / PI / AR);
	L = D_fr * AR;
	//A_fr = V_0 / L;
	//D_fr = sqrt(4 * A_fr / PI);
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
}

void RegeneratorModel::setDesignTargets(targetModes::targetModes targetMode, targetModes::target2Modes secondTargetMode, double targetParameter, double secondTargetParameter)
{
	this->secondTargetMode = secondTargetMode;
	this->targetdP_max_Regen = secondTargetParameter;
	
	this->targetMode = targetMode;
	this->targetParameter = targetParameter;
}


void RegeneratorModel::setPoint(Eigen::VectorXd point) {
	this->T_H_out = point(0)*T_H_in;
	this->AR = point(1)*max_Size;
	this->dP_C = point(2)*P_C_in;
	this->V_0 = point(3)*max_Size;
}

void RegeneratorModel::jacobian(int n, Eigen::VectorXd x_n, Eigen::VectorXd dx, Eigen::MatrixXd& Jn) {
	Eigen::VectorXd f_i(n);
	Eigen::VectorXd f_im1(n);
	Eigen::VectorXd f_ip1(n);
	Eigen::VectorXd x_tmp(n);
	evaluate(x_n, f_i);

	x_tmp = x_n;

	// Second order central sheme
	// Due to stability issues only the diagonal is calculated. More on that in documentation
	for (int i = 0; i < n; i++) {
		x_tmp(i) -= dx(i);
		evaluate(x_tmp, f_im1);
		x_tmp(i) += 2*dx(i);
		evaluate(x_tmp, f_ip1);
		Jn(i, i) = (f_im1(i) + f_ip1(i) - 2 * f_i(i)) / pow(dx(i), 2);
	}
}

void RegeneratorModel::evaluate(Eigen::VectorXd point, Eigen::VectorXd& f) {
	setPoint(point);
	calculateModel();

	// Function values are normalized
	f(0) = (Q_dot_a - Q_dot_a_calc) / (Q_dot_a + Q_dot_a_calc);
	f(1) = (dP_H - dP_H_calc)/P_H_in;
	f(2) = (dP_C - dP_C_calc)/P_C_in;
	f(3) = (UA - targetParameter) / (targetParameter * 3);
}

int RegeneratorModel::getDesignSolution()
{
	clock_t begin = clock();

	// Number of functions/variables
	int n = 4;

	// Previous coordiante
	Eigen::VectorXd x_0(n);

	// Current coordiante
	Eigen::VectorXd x_1(n);

	dP_H = targetdP_max_Regen;

	// Variables are normalized
	x_0(0) = (T_H_in + T_C_in) / 2 / T_H_in;
	x_0(1) = 1.1 / max_Size;
	x_0(2) = 0.3*targetdP_max_Regen / P_C_in;
	x_0(3) = 0.5 / max_Size;

	// Step size for finite-difference
	Eigen::VectorXd dx(n);
	dx(0) = 1 / T_H_in;
	dx(1) = 0.1 / max_Size;
	dx(2) = 1 / P_C_in;
	dx(3) = 0.1 / max_Size;

	Eigen::MatrixXd Jn(n, n);
	Eigen::MatrixXd Jn_inv(n, n);
	Eigen::MatrixXd Jn_tmp(n, n);
	Eigen::VectorXd f_0(n);
	Eigen::VectorXd f_1(n);
	Eigen::VectorXd del_f(n);
	Eigen::VectorXd del_x(n);

	// Calculate Jacobian and function values before entering iteration loop
	jacobian(n, x_0, dx, Jn);
	evaluate(x_0, f_0);

	///Solution tolerance
	double eps = 1e-5;

	int iterations = 0;
	while (fabs(f_0(0)) > eps || fabs(f_0(1)) > eps || fabs(f_0(2)) > eps || fabs(f_0(3)) > eps) {
		
		// Check if Jacobian matrix is invertable
		if (!Jn.fullPivLu().isInvertible()) {
			return 1;
		}
		else {
			Jn_inv = Jn.inverse();
		}

		// Take a single step
		x_1 = x_0 - (Jn_inv * f_0);
		del_x = x_1 - x_0;

		/// If new variable values are not physically possible: recalculate Jacobian using jacobian(), update necessary entries in the matrix, redo the step.
		if (x_1(0) <= T_C_in / T_H_in || x_1(0) >= 1 || x_1(1) < 0 || x_1(1) >= 1 || x_1(2) < 0 || x_1(3) < 0 || x_1(3) >= 1) {
			jacobian(n, x_0, dx, Jn_tmp);
			for (int i = 0; i < n; i++) {
				if (x_1(i) < 0 || x_1(i) >= 1) {
					Jn(i, i) = Jn_tmp(i, i);
				}
			}

			// Hot outlet temperature cannot drop below cold inlet temperature
			if (x_1(0) <= T_C_in / T_H_in) {
				Jn(0, 0) = Jn_tmp(0, 0);
			}

			// Check if Jacobian matrix is invertable
			if (!Jn.fullPivLu().isInvertible()) {
				return 1;
			}
			else {
				Jn_inv = Jn.inverse();
			}

			// Take a single step
			x_1 = x_0 - (Jn_inv * f_0);
			del_x = x_1 - x_0;
		}

		evaluate(x_1, f_1);
		del_f = f_1 - f_0;

		///"Good Broyden" method of updating Jacobian. Calculation of Jacobian is costly and is avoided.
		Jn = Jn + (del_f - Jn * del_x)*del_x.transpose() / del_x.squaredNorm();

		// Which variables should have an effect on which functions
		Jn(0, 1) = 0;
		Jn(0, 2) = 0;
		Jn(0, 3) = 0;

		Jn(1, 0) = 0;
		Jn(1, 2) = 0;

		Jn(2, 0) = 0;
		Jn(2, 1) = 0;
		Jn(2, 3) = 0;

		Jn(3, 0) = 0;
		Jn(3, 1) = 0;
		Jn(3, 2) = 0;

		x_0 = x_1;
		f_0 = f_1;

		iterations++;
	}

	clock_t end = clock();
	double elapsed = double(end - begin) / CLOCKS_PER_SEC * 1000;
	cout << T_H_out << ", " << L << ", " << D_fr << ", " << V_0 << ", " << AR << ", " << UA << ", " <<
		dP_max << ", " << epsilon << ", " << (Q_dot_a - Q_dot_a_calc) << ", " << (dP_H - dP_H_calc) << ", " << (dP_C - dP_C_calc)
		<< ", " << targetParameter << ", " << targetdP_max_Regen << ", " << iterations << ", " << elapsed << endl;

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
	solution->m_HTR_HP_dP = dP_C;
	solution->m_HTR_LP_dP = dP_H;
	solution->L = L;
	solution->D_fr = D_fr;
}
