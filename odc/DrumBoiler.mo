model DrumBoiler "Process model for dynamic and steady-state optimization examples"
  extends Modelica.Fluid.Examples.DrumBoiler.DrumBoiler(use_inputs = true);
  Modelica.Blocks.Interfaces.RealOutput sigma_D(unit = "N/mm2") "thermal stress of drum" annotation(Placement(transformation(extent = {{100, -66}, {112, -54}}, rotation = 0)));
equation
  sigma_D = (-1e3 * der(evaporator.T_D)) + 1e-5 * evaporator.p;
  annotation(uses(Modelica(version = "3.2.1")), experiment(StopTime = 3600), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics), Documentation(info = "<html>
  <p>This process model serves as simple test case for optimization. It exhibits the following features:</p>
  <p><ul>
  <li>based on the Modelica Standard Library, including Media and Fluid</li>
  <li>simple embedded control using Blocks, including PI controller and discontinuous limiter</li>
  <li>three states: drum boiler pressure, drum boiler water level, integral part of controller</li>
  </ul></p>
  <p>The optimization should care about</p>
  <p><ul>
  <li>usage of regular Modelica simulation model as optimization contraints</li>
  <li>initialization of states from parameters</li>
  <li>treatment of nominal values (e.g. 100 bar is represented as 1e7 Pa in the pressure state, while the valve is fully open at 1)</li>
  <li>physical state/output constraint for thermal stress</li>
  <li>controller state/output constraint for limiter (i.e. avoid negative feedwater flow)</li>
  </ul></p>
 The model extends from Modelica.Fluid.Examples.DrumBoiler.DrumBoiler and adds the new output sigma_D for thermal stress and membrane stress</p>
  <pre>
    sigma_D = -1e3*der(evaporator.T_D) + 1e-5*evaporator.p
  </pre>
  <p>For more information see: </p>
  <p><i>R. Franke, M. Rode, K. Kr&uuml;ger: On-line Optimization of Drum Boiler Startup, Modelica 2003. </i></p>
  <p><a href=\"https://www.modelica.org/events/Conference2003/papers/h29_Franke.pdf\">https://www.modelica.org/events/Conference2003/papers/h29_Franke.pdf</a></p>

  <h4>Simulation model</h4>
  <p>The Modelica model can be translated to an ODE of the form:</p>
  <pre>
    der(x) = f(x,u)
        y  = h(x,u)
  </pre>
  <p>with:</p>
  <pre>
    x = {controller.x, evaporator.p, evaporator.V_l}
    u = {q_F, Y_Valve}
    y = {p_S, qm_S, sigma_D, T_S, V_l}
  </pre>
  <p>The initial states are defined by model parameters:</p>
  <pre>
    x(0) = f0(p)
  </pre>
  <p>with the parameter values:</p>
  <pre>
    evaporator.p_start = 1 bar
    evaporator.V_start = 67 m3
    controller.x_start = 0
  </pre>
  <p>The simulation model serves as basis for multiple optimization fomulations, i.e. experiments.</p>

  <h4>Trajectory optimization (drumboiler.tcl)</font></h4>
  <p>The aim is to obtain an optimal startup control for fuel flow rate and steam flow considering a constraint on termal stress. </p>
  <p>Minimize:</p>
  <pre>
   3600s
     &int; 1e-7*(p_S-110bar)^2 + 1e-8*(qm_S-180kg/s)^2 + 1e-4*(dq_F/dt)^2 dt
     0
  </pre>
  <p>subject to the model:</p>
  <pre>
    der(x) = f(x,u)
      x(0) = f0(p)
  </pre>
  <p>the control bounds:</p>
  <pre>
    0          &LT;= Y_Valve  &LT;= 1
    0          &LT;=   q_F    &LT;= 500 MW
    -24 MW/min &LT;= der(q_F) &LT;= 24 MW/min
    q_F(0) = 0
  </pre>
  <p>and the state/output constraint:</p>
  <pre>
    -150 N/mm2 &LT;= sigma_D
  </pre>
  <p>An appropriate discretization of the control inputs is piecewise linear with a step size of 60s. The model variables can be initialized with the results of an initial-value simulation keeping the valve fully open and constantly ramping up the fuel flow rate, e.g. by 400MW/1h.</p>

  <h4>Set point optimization (drumboiler_sp.tcl)</h4>
  <p>The goal of this optimization is to find a steady state together with values for the fuel flow rate q_F and the opening of the steam valve Y_Valve for which the heat input is minimized, subject to required steam pressure and mass flow rate.</p>
  <p>This is:</p>
  <pre>
    q_F  --&GT;  min
                 x,u
   </pre>
  <p>subject to the steady-state model:</p>
  <pre>
    0 = f(x,u)
  </pre>
  <p>and the constraints:</p>
  <pre>
    0        &LT;= q_F     &LT;= 500 MW
    0        &LT;= Y_Valve &LT;= 1
    100 bar  &LT;= p_S     &LT;= 120 bar
    qm_S = 150 kg/s
  </pre>
  <p>The solution can be found at q_F = 328 MW, Y_Valve = 0.63 with the states evaporator.p = 120 bar, evaporator.V_liquid = 67 m3, and controller.x = 15.</p>
  </html>"));
end DrumBoiler;