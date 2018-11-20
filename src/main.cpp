#include <string>
#include <iostream>
#include "RegenHX.h"

int main(int argc, char **argv)
{
	RegenHX* regenerator = new RegenHX();
	regenerator->initialize(10);
	regenerator->set_params(2, 0, 0, 1, 215, 145, 0.003, 0.37, 100);
	double q_dot, T_c_out, T_h_out;
	regenerator->design_fix_UA_calc_outlet(7000, 0.98, 444, 25000, 70, 24900, 800, 7600, 70, 7380, q_dot, T_c_out, T_h_out);
	cout << "Q_dot = " << q_dot << ", T_C_out = " << T_c_out << ", T_H_out = " << T_h_out << endl;
	cout << "Efficiency " << regenerator->ms_des_solved.m_eff_design;
}
