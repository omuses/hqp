model DID "Double Integrator Discrete-time"
  parameter Real p = 1 "gain for input";
  parameter Real y1_start = 1 "start value for first state";
  parameter Real y2_start = 0 "start value for second state";
  input Real u(start = -2);
  Real x1(start = y1_start), x2(start = y2_start);
  output Real y1, y2;
equation
  when Clock() then
    x1 = previous(x1) + p * u * interval(u);
    x2 = previous(x2) + previous(x1) * interval(x1) + 0.5 * p * u * interval(u)^2;
    y1 = previous(x1);
    y2 = previous(x2);
  end when;
end DID;
