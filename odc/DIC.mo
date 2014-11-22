model DIC "Double Integrator Continuous-time"
  parameter Real p = 1 "gain for input";
  parameter Real y1_start = 1 "start value for first state";
  parameter Real y2_start = 0 "start value for second state";
  input Real u(start = -2);
  output Real y1, y2;
initial equation
  y1 = y1_start;
  y2 = y2_start;
equation
  der(y1) = p * u;
  der(y2) = y1;
end DIC;
