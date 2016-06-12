model DID "Double Integrator Discrete-time"
  parameter Real p = 1 "gain for input";
  parameter Real y1_start = 1 "start value for first state";
  parameter Real y2_start = 0 "start value for second state";
  input Real u(start = -2);
  output Real y1(start = y1_start, fixed = true);
  output Real y2(start = y2_start, fixed = true);
  Real ud;
equation
  ud = sample(u, Clock(Clock(/*inferred*/), solverMethod = "ImplicitEuler"));
  der(y1) = p * ud;
  der(y2) = previous(y1) + 0.5 * p * ud * interval(ud);
end DID;
