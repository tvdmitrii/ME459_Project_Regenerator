#include <string>
#include <iostream>
#include "RegenHX.h"

int main(int argc, char **argv)
{
	double dP_max = 100;
	double UA = 1000;
	RegenHX* regenerator = new RegenHX();
	regenerator->initialize(10);
	spdlog::get("logger")->info("T_H_out,L,D_fr,V_0,AR,UA,dP_max,epsilon,Q_dot_diff,dP_H_diff,dP_C_diff,targetParameter,targetdP_max_Regen,Time");

	double q_dot, T_c_out, T_h_out;
	while (UA <= 10000) {
		dP_max = 100;
		while (dP_max <= 300) {
			regenerator->set_params(2, 0, 0, 1, dP_max, 145, 0.003, 0.37, 100);
			try {
				regenerator->design_fix_UA_calc_outlet(UA, 0.98, 444, 25000, 70, 24900, 800, 7600, 70, 7600 - dP_max, q_dot, T_c_out, T_h_out);
			}
			catch (exception e) {
				spdlog::get("logger")->info("--" + to_string(UA) + ", " + to_string(dP_max) + "--");
			}
			dP_max += 25;
		}
		UA += 250;
	}
}

