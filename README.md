HQP: a solver for sparse nonlinear optimization
===============================================

HQP is a solver for nonlinearly constrained large-scale optimization. It is intended for problems with sufficient regular sparsity structure. Such optimization problems arise e.g. from the numerical treatment of optimal control problems. External interfaces allow the formulation of optimization problems based on widely used model formats.

<h2>Overview</h2>

HQP (Huge Quadratic Programming) consists of mainly two parts: the actual HQP optimizer and the front-end Omuses. Both parts are designed as framework in the programming language C++.

<h3>HQP</h3>

The actual HQP optimizer treats nonlinearly constrained problems with a sequential quadratic programming (SQP) algorithm. An interior-point method is applied to the solution of convex quadratic subproblems.


The implementation is based on sparse matrix codes of the [Meschach library for matrix computations in C](http://www.netlib.org/c/meschach). The matrix library was extended with additional routines for the analysis and direct solution of sparse equation systems.  

The tool command language [Tcl](http://www.tcl.tk) is used for selecting solver modules, configuring parameters and for controlling the execution.

<h3>Omuses</h3>

The front-end Omuses provides additional support for the efficient problem formulation. This is possible thanks to the availability of great software packages that have been integrated with Omuses. 

[ADOL-C](https://projects.coin-or.org/ADOL-C) is exploited for the automatic differentiation and structural analysis of model equations.

Furthermore Omuses provides numerical solvers for differential equations defining constraints in dynamic optimization
problems. Besides own implementations (Dopri5, Euler, GRK4, IMP, OdeTs, SDIRK, RK4), the following additional software packages are currently integrated:

- [DASPK](http://www.engineering.ucsb.edu/~cse/software.html): solution of implicit and stiff differential equations
- [RKsuite](http://www.netlib.org"): solution of explicit differential equations

Please note that the integrated software packages underly copyright restrictions of their respective authors. That is why the DASPK software is not included, though an interface is provided. 

<h2>Interfaces for formulating optimization problems</h2>

<h3>External Model Interfaces</h3>

- Prg_SFunctionOpt/Prg_SFunctionEst: pre-formulated dynamic optimization/estimation problems for a model accessible as [Simulink](http://www.mathworks.com) S-function
- Prg_CUTE: implementation of interface to [CUTE](http://www.netlib.org/utk/misc/sw_survey/urc/html/215.1.html) (Constrained and Unconstrained Testing Environment), providing a large collection of sparse nonlinear optimization problems formulated in the Standard Input Format [SIF](http://www.numerical.rl.ac.uk/lancelot/sif/sifhtml.html) 

<h3>Internal Model Interfaces</h3>

Optimization problems can be formulated natively in C/C++. The following interfaces exist (sorted from high-level to low-level):

- Omu_Program: dynamic optimization problems exploiting the front-end Omuses, including treatment of differential equations and automatic differentiation (see the odc subdirectory for examples). Omuses converts a dynamic optimization problem to a discrete-time optimal control problem
- Hqp_Docp: discrete-time optimal control problem (DOCP) for the Hqp solver, specified as Hqp_DocpSpec (see the hqp_docp subdirectory for an example). A DOCP gets converted to a large-scale non-linear optimization problem
- Hqp_SqpProgram: large-scale nonlinear optimization problem, as actually being treated by the HQP solver

<h2>References</h2>

<h3>Applications</h3>

The most interesting applications, where HQP is running round-the-clock to help optimizing our world, include

- Online optimization of power plants
- Water management in a large canal system

<b>See:</b>

R. Franke and B. Weidmann: Steaming ahead with optimizing power plant boiler startup.
<i>Power Engineering International -- PEi</i>, 15(6), 2007.

R. Franke and L. Vogelbacher: Nonlinear model predictive control for cost optimal startup of steam power plants. <i>at -- Automatisierungstechnik</i>, 54(12):630--637, 2006.

H. Linke, E. Arnold, and R. Franke: Optimal water management of a canal system. In J.C. Refsgaard and E.A. Karalis, editors, <i>Operational Water Management</i>, pages 211--218. Balkema, Rotterdam, 1997.

<h3>Research Applications</h3>

Furthermore, there finds several applications to the research on

- batch process control
- energy systems
- water systems

<b>See:</b>

Z.K. Nagy, B. Mahn, R. Franke, and F. Allg&ouml;wer: Evaluation study of an efficient output feedback nonlinear model predictive control for temperature tracking in an industrial batch reactor. <i>Control Engineering Practice</i>, 15(7):839 -- 850, 2007.

C. Hoffmann and H. Puta: Dynamic optimization of energy supply systems with Modelica models. In <i>Proceedings of the 5th International Modelica Conference.</i> Modelica Association, Vienna, Austria, September 2006.

H. Linke: <i>Wasserbewirtschaftung von Binnenschiffahrtsgew&auml;ssern auf Basis einer modellgest&uuml;tzten Vorhersage des Systemverhaltens.</i> Dissertation, Cuvillier Verlag, G&ouml;ttingen, 2006.

G. Reichl: <i>Optimierte Bewirtschaftung von Kl&auml;ranlagen basierend auf der Modellierung mit Modelica.</i> Dissertation, Cuvillier Verlag, G&ouml;ttingen, 2005.

R. Nystr&ouml;m, R. Franke, I. Harjunkowski, and A. Kroll: Production campaign planning including grade transition sequencing and dynamic optimization. <i>Computers and Chemical Engineering</i>, 29:2163--2179, 2005.

R. Franke: Integrierte dynamische Modellierung und Optimierung von Systemen mit saisonaler W&auml;rmespeicherung, Dissertation, volume 394 of <i>Fortschritt-Berichte VDI</i>, series 6. VDI-Verlag, D&uuml;sseldorf, 1998.


<h2>Directory structure of HQP distribution</h2>

The HQP distribution contains the following sub-directories:

<h3>Main modules</h3>

<b>hqp</b>: The actual solver for sparse nonlinear optimization

<b>omu</b>: Front-end Omuses for treating dynamic optimization problems

<h3>Examples</h3>

<b>odc</b>: Omuses Demo Collection (examples using the front-end Omuses)

<b>hqp_docp</b>: Example for directly using the DOCP interface of HQP

<h3>External model interfaces</h3>

<b>hxi</b>: Hqp eXternal Interfaces (currently S-function)

<b>hqp_cute</b>: Files required to run HQP with CUTE

<h3>Additional differential equation solvers</h3>

<b>daspk</b>: Place holder for DASPK differential-algebraic equation solver

<b>rksuite</b>: RKsuite differential equation solvers (by R.W. Brankin et al)

<h3>Base modules</h3>

<b>meschach</b>: Matrix library (by D.E. Steward and Z. Leyk)

<b>adol-c</b>: Automatic differentiation code (by A. Griewank et al)

<b>iftcl</b>: Interface wrapping Tcl (Tool Command Language; used for configuring solver options and for controlling the execution)

<b>malloc</b>: GNU malloc library (HQP may be configured to use GNU malloc instead of the system malloc)

<h3>Documentation</h3>

<b>doc</b>: Doxygen input file for generating reference documentation. 

Please see also the separately available User Manual (omuses.ps.gz), including explanation of simple examples.

<h2>License</h2>

This software is free according to the conditions of the GNU LIBRARY GENERAL PUBLIC LICENSE, Version 2 (see COPYING.LIB).
