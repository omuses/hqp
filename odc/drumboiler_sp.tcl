#
# DrumBoiler setpoint optimization example
#   providing the model via FMU generated from Modelica.Fluid
#

set mdl_name DrumBoiler

## translate DrumBoiler.mo to DrumBoiler.fmu

puts "Compiling FMU..."
source omc.tcl
if {![file exists $mdl_name.fmu]
    || [file mtime $mdl_name.mo] > [file mtime $mdl_name.fmu]} {
    compileFMU $mdl_name.mo
}

## run optimization

puts "Running optimization..."
package require Omuses

## configure optimization program

# optimization program and model name
prg_name DynamicOpt
mdl_path $mdl_name.fmu

# load and setup model
prg_setup_stages

# steady state condition
# controller.x evaporator.V_l evaporator.p
mdl_x0_active  {1 1 1}
mdl_der_x0_min {0 0 0}
mdl_der_x0_max {0 0 0}

# optimized (active) inputs
#               q_F Y_Valve
mdl_u_active  {   1    1}
mdl_u_max     { 500    1}
mdl_u_min     {   0    0}

# optimization objective
mdl_u_weight1 {   1    0}

# constraints on outputs
#               p_S qm_S sigma_D T_S V_l
mdl_y_max     { 120  150   Inf   Inf  Inf}
mdl_y_min     { 100  150  -Inf  -Inf -Inf}

## setup and solve optimization program

# initial solution
mdl_x0 {0 1e5 65}
mdl_u0 {0 1}

# setup optimization program
prg_setup

# initialize optimizer
sqp_init

# drive solution
catch hqp_solve result

# summarize results
puts "Result   : $result"
puts "Objective: [prg_f]"
puts "Obj-evals: [prg_fbd_evals]"

proc putsVars {names vals} {
    for {set i 0} {$i < [llength $names]} {incr i} {
	puts [format "%15s: %.3g"  [lindex $names $i] [lindex $vals $i]]
    }
}

puts "\nOptimized steady state, inputs, outputs:"
putsVars [set ::fmu::DrumBoiler::stateNames] [mdl_x0]
putsVars [set ::fmu::DrumBoiler::inputNames] [mdl_u0]
putsVars [set ::fmu::DrumBoiler::outputNames] [mdl_y0]
