#pragma once
#include "numeric_solvers.h"
#include <stdlib.h>
#include <string>

using namespace std;

template <class T>
struct SolverParameters {
	string solverName;
	double target;
	double guessValue1;
	double guessValue2;
	double lowerBound;
	double upperBound;
	double tolerance;
	int iterationLimit;
	bool isErrorRel;
	T* classInst;
	int (T::*monoEquation)(double x, double* y);
};

template <class T>
class MonoSolver
{
private:
	
	string solverName;

	double target;

	double guessValue1;
	double guessValue2;

	double lowerBound;
	double upperBound;

	double tolerance;

	int iterationLimit;

	bool isErrorRel;

	C_monotonic_equation* equation = nullptr;

	C_monotonic_eq_solver* eqSolver = nullptr;

	bool inProgress;
	
public:
	MonoSolver(SolverParameters<T>* params)
	{
		solverName = params->solverName;
		target = params->target;
		guessValue1 = params->guessValue1;
		guessValue2 = params->guessValue2;
		lowerBound = params->lowerBound;
		upperBound = params->upperBound;
		tolerance = params->tolerance;
		iterationLimit = params->iterationLimit;
		isErrorRel = params->isErrorRel;

		equation = new C_member_mono_eq<T>(params->classInst, params->monoEquation);
		eqSolver = new C_monotonic_eq_solver(*equation);
		eqSolver->settings(tolerance, iterationLimit, lowerBound, upperBound, isErrorRel);

		inProgress = false;
	}

	MonoSolver()
	{
		inProgress = false;
	}

	void setParameters(SolverParameters<T>* params)
	{
		if (inProgress) {
			throw exception();
		}

		solverName = params->solverName;
		target = params->target;
		guessValue1 = params->guessValue1;
		guessValue2 = params->guessValue2;
		lowerBound = params->lowerBound;
		upperBound = params->upperBound;
		tolerance = params->tolerance;
		iterationLimit = params->iterationLimit;
		isErrorRel = params->isErrorRel;

		if (equation != nullptr && eqSolver != nullptr) {
			delete equation;
			delete eqSolver;
		}
		equation = new C_member_mono_eq<T>(params->classInst, params->monoEquation);
		eqSolver = new C_monotonic_eq_solver(*equation);
		eqSolver->settings(tolerance, iterationLimit, lowerBound, upperBound, isErrorRel);

		inProgress = false;
	}
	
	void updateGuesses(double guess1, double guess2) {
		if (inProgress) {
			throw exception();
		}

		guessValue1 = guess1;
		guessValue2 = guess2;
	}

	void updateTarget(double target) {
		if (inProgress) {
			throw exception();
		}

		this->target = target;
	}

	void updateTolerance(double tolerance) {
		if (inProgress) {
			throw exception();
		}

		this->tolerance = tolerance;
	}

	int solve(double* tolSolved = nullptr)
	{
		inProgress = true;
		eqSolver->reset();

		double solution, toleranceSolved;
		solution = toleranceSolved = numeric_limits<double>::quiet_NaN();
		int iter_solved = -1;

		int status = eqSolver->solve(guessValue1, guessValue2, target, solution, toleranceSolved, iter_solved);

		if (tolSolved != nullptr && toleranceSolved != numeric_limits<double>::quiet_NaN()) {
			*tolSolved = toleranceSolved;
		}

		inProgress = false;
		return status;
	}

	~MonoSolver() {
		delete equation;
		delete eqSolver;
	};
};
