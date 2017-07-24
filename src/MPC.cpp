#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 20; // anything more than 6 i.e no. of forward points showing in simulator
double dt = .05; // ie. 20Hz

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// TODO: inserted positions
// Both the reference cross track and orientation errors are 0.
// The reference velocity is set to 40 mph ~ 18m/s.
double ref_v = 18.;

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.

// TODO: Set the number of model variables (includes both states and inputs).
// For example: If the state is a 4 element vector, the actuators is a 2
// element vector and there are 10 timesteps. The number of variables is:
// #states * (N_timestamps) + #actuators * (N_timestamps - 1) = 4*10+2*(10-1)=58
// N timesteps == N - 1 actuations
// TODO: Set the number of independent variables, constraints
size_t _X_ = 0;
size_t _Y_ = _X_ + N;
size_t _Psi_ = _Y_ + N;
size_t _V_ = _Psi_ + N;
size_t _cte_ = _V_ + N;
size_t _epsi_ = _cte_ + N;
size_t n_constraints = _epsi_ + N; // N * 6;
size_t _delta_ = n_constraints;
size_t _acc_ = _delta_ + N - 1;
size_t n_vars = _acc_ + N - 1; // N * 6 + (N - 1) * 2;

  typedef CPPAD_TESTVECTOR(double) Dvector;
void setBounds(Dvector& vars_lowerbound, Dvector& vars_upperbound); //called in Solve()
void setConstraints(const Eigen::VectorXd state, Dvector& constraints_lowerbound, Dvector& constraints_upperbound); //called in Solve()

class FG_eval {
 public:
  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  AD<double> getCost(const ADvector& vars); //called in operator()

  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    //fetch different state variables from vars[] at offsets setup at the top
#define ip(A,i) vars[A+i]
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
#define eqn(A,i) fg[1+A+i]
	  
    // Setup model Cost
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
	fg[0] = getCost(vars);

    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.
    // The idea here is to constraint this value to be 0.

    // Initial constraints
	{
      eqn(_X_,0) = ip(_X_,0);
      eqn(_Y_,0) = ip(_Y_,0);
      eqn(_Psi_,0) = ip(_Psi_,0);
      eqn(_V_,0) = ip(_V_,0);
      eqn(_cte_,0) = ip(_cte_,0);
      eqn(_epsi_,0) = ip(_epsi_,0);
	}

    // The rest of the constraints
    // For each state param, set the constraints, as eqn = actual(t) - estimate(t) = 0.
    for (size_t t = 0; t < N-1; t++) {
      // Considers state at time t+1, t, and actuation at time t

      // Frequently used params defined
	  AD<double> x0 = ip(_X_,t);
      AD<double> psi0 = ip(_Psi_,t);
      AD<double> v0 = ip(_V_,t);
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * pow(x0,2) + coeffs[3] * pow(x0,3);
      AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * x0 * coeffs[2] + 3 * pow(x0,2) * coeffs[3] );

      //
      // Recall the equations for the model:
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
      eqn(_X_,t+1) = ip(_X_,t+1) - (ip(_X_,t) + v0 * CppAD::cos(psi0) * dt);
      eqn(_Y_,t+1) = ip(_Y_,t+1) - (ip(_Y_,t) + v0 * CppAD::sin(psi0) * dt);
      eqn(_Psi_,t+1) = ip(_Psi_,t+1) - (psi0 + v0 * ip(_delta_,t) / Lf * dt);
      eqn(_V_,t+1) = ip(_V_,t+1) - (v0 + ip(_acc_,t) * dt);
      eqn(_cte_,t+1) = ip(_cte_,t+1) - ((f0 - ip(_Y_,t)) + (v0 * CppAD::sin(ip(_epsi_,t)) * dt));
      eqn(_epsi_,t+1) = ip(_epsi_,t+1) - ((psi0 - psides0) + v0 * ip(_delta_,t) / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
Eigen::VectorXd wts(10);
MPC::MPC() {
	wts << 40., 12, .1, 3000., 1000., 1., 1., 1., 200., 1.;
	const char *names[] = {"ref_vel", "N fwd", "dt", "cte", "pse", "v", "delta", "acc", "delta_change", "acc_change"};
	for (auto i=0;i<10; i++) {
//		cout << names[i] << ": "; cin >> wts[i];
	} ref_v = wts[0]; N = wts[1]; dt = wts[2];  
}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  //size_t i;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (size_t i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }
  // Set the initial variable values
  ip(_X_,0)   = state[0];
  ip(_Y_,0)   = state[1];
  ip(_Psi_,0) = state[2];
  ip(_V_,0)   = state[3];
  ip(_cte_,0) = state[4];
  ip(_epsi_,0)= state[5];

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.  
  setBounds(vars_lowerbound, vars_upperbound);
  
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  setConstraints(state, constraints_lowerbound, constraints_upperbound);

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << ok << " Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> result;
  
  result.push_back(solution.x[_delta_]);
  result.push_back(solution.x[_acc_]);
  
  for (size_t i = 0; i < N-1; i++)
  {
    result.push_back(solution.x[_X_ + i + 1]);
    result.push_back(solution.x[_Y_ + i + 1]);
  }
  return result;
}

void setBounds(Dvector& vars_lowerbound, Dvector& vars_upperbound) {
  //in SI units of radians, m, s.

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (size_t i = 0; i < n_constraints; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (size_t i = n_constraints; i < _acc_; i++) {
    vars_lowerbound[i] = -25.*M_PI/180.;
    vars_upperbound[i] =  25.*M_PI/180.;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (size_t i = _acc_; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
}

// Our constraints are all equality, so same lower and upper bounds
void setConstraints(const Eigen::VectorXd state, Dvector& constraints_lowerbound, Dvector& constraints_upperbound) {
  for (size_t i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[_X_]   = state[0];
  constraints_lowerbound[_Y_]   = state[1];
  constraints_lowerbound[_Psi_] = state[2];
  constraints_lowerbound[_V_]   = state[3];
  constraints_lowerbound[_cte_] = state[4];
  constraints_lowerbound[_epsi_]= state[5];

  constraints_upperbound[_X_]   = state[0];
  constraints_upperbound[_Y_]   = state[1];
  constraints_upperbound[_Psi_] = state[2];
  constraints_upperbound[_V_]   = state[3];
  constraints_upperbound[_cte_] = state[4];
  constraints_upperbound[_epsi_]= state[5];
}


// TODO: Set cost weights to tune the model for a 'smooth ride'
AD<double> FG_eval::getCost(const ADvector& vars) {
    AD<double> cost = 0;
	
    // The part of the cost based on the reference state.
    for (size_t t = 0; t < N; t++) {
      cost += wts[3]*CppAD::pow(ip( _cte_,t), 2);
      cost += wts[4]*CppAD::pow(ip(_epsi_,t), 2);
      cost += wts[5]*CppAD::pow(ip(_V_,t) - ref_v, 2);
    }

    // Minimize the use of actuators.
    for (size_t t = 0; t < N - 1; t++) {
      cost += wts[6]*CppAD::pow(ip(_delta_,t), 2);
      cost += wts[7]*CppAD::pow(ip( _acc_ ,t), 2);
    }

    // Minimize the value gap between sequential actuations.
    for (size_t t = 1; t < N - 1; t++) {
      cost += wts[8]*CppAD::pow(ip(_delta_,t) - ip(_delta_,t-1), 2);
      cost += wts[9]*CppAD::pow(ip( _acc_ ,t) - ip( _acc_ ,t-1), 2);
    }
	return cost;
}