#include <string>
#include <iostream>
#include "RegenHX.h"

int main(int argc, char **argv)
{
	if (argc != 3) {
		printf("Usage: %s <UA [kW/K]> <dP_max [kPa]>\n", argv[0]);
		return 1;
	}

	double dP_max = atof(argv[2]);
	double UA = atof(argv[1]);
	RegenHX* regenerator = new RegenHX();
	regenerator->initialize(10);
	cout << "T_H_out, L, D_fr, V_0, AR, UA, dP_max, epsilon, Q_dot_diff, dP_H_diff, dP_C_diff, targetParameter, targetdP_max_Regen, Time" << endl;

	double q_dot, T_c_out, T_h_out;
	regenerator->set_params(2, 0, 0, 1, dP_max, 145, 0.003, 0.37, 100);
	try {
		regenerator->design_fix_UA_calc_outlet(UA, 0.98, 444, 25000, 70, 24900, 800, 7600, 70, 7600 - dP_max, q_dot, T_c_out, T_h_out);
	}
	catch (exception e) {
		cout << "Could not converge on solution";
	}

	return 0;
}
