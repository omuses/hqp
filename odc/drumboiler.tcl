#
# DrumBoiler trajectory optimization example
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

## configure program

# optimization program and model name
prg_name DynamicOpt
mdl_path $mdl_name.fmu

# setup stages and sample periods per stage
prg_KK 60
prg_setup_stages

# model inputs and sample time points
set KK [prg_KK]
set dt [expr 3600.0/$KK]
set us {}
set ts {}
for {set kk 0} {$kk <= $KK} {incr kk} {
    # q_F Y_Valve
    lappend us "[expr 400.0*$kk/$KK] 1.0"
    lappend ts [expr $dt*$kk]
}
mdl_us $us
prg_ts $ts

# define controlled (active) inputs
#                  q_F  Y_Valve
mdl_u_active      {   1    1}
mdl_u_order       {   1    1}
mdl_u_max         { 500    1}
mdl_u_min         {   0    0}
mdl_u0_max        {   0    1}
mdl_u0_min        {   0    0}
mdl_der_u_max     { 0.4  Inf}
mdl_der_u_min     {-0.4 -Inf}
mdl_der_u_weight2 {1e-4    0}

# define optimization objective
#              p_S qm_S sigma_D T_S V_l
mdl_y_nominal { 100  100  100  100  100}
mdl_y_ref     { 110  180  0.0  0.0  0.0}
mdl_y_weight2 {1e-3 1e-4  0.0  0.0  0.0}

# define constraint on thermal stress
mdl_y_min     {-Inf -Inf -150 -Inf -Inf}

## setup solver
prg_integrator IMP
prg_int_stepsize 30
sqp_eps 1e-4

## solve
proc solve {} {
    # setup optimization program
    prg_setup

    # simulate initial solution
    prg_simulate

    # plot results of initial simulation
    plot

    # initialize optimizer
    sqp_init

    # drive solution and plot after each iteration
    for {set iter 1} {$iter <= 20} {incr iter} {
	sqp_max_iters $iter
	catch hqp_solve result
	plot
	if {$result != "iters"} break
    }

    # summarize results
    puts "Result   : $result"
    puts "Objective: [prg_f]"
    puts "Obj-evals: [prg_fbd_evals]"
}

## plot results if Tk and BLT are present
if {[catch {package present Tk}] 
    || ([catch {package require BLT}] ?
	[catch {package require rbc}] : 0)} {
    puts "Skipping plot as no blt::graph available."
    proc plot {} {}
} else {
  proc plot {} {
    # read results and form vector of cumulative times
    set ts [prg_ts]
    set us [mdl_us]
    set ys [mdl_ys]
    set KK [prg_KK]
    set u_names [set ::fmu::DrumBoiler::inputNames]
    set y_names [set ::fmu::DrumBoiler::outputNames]
    foreach u_name $u_names {
	set vars($u_name) {}
    }
    foreach y_name $y_names {
	set vars($y_name) {}
    }
    for {set kk 0} {$kk <= $KK} {incr kk} {
	set ukk [lindex $us $kk]
	set ykk [lindex $ys $kk]
	set i 0
	foreach u_name $u_names {
	    lappend vars($u_name) [lindex $ukk $i]
	    incr i
	}
	set i 0
	foreach y_name $y_names {
	    lappend vars($y_name) [lindex $ykk $i]
	    incr i
	}
    }
    # create blt::graph in a toplevel window and plot data
    destroy .graph
    if {![catch {package present BLT}]} {
	set graph [blt::graph .graph -title "[prg_name] demo"]
    } else {
	set graph [rbc::graph .graph -title "[prg_name] demo"]
    }
    $graph element create q_F -xdata $ts -ydata $vars(q_F) \
	-color black -symbol ""
    $graph element create Y_Valve -xdata $ts -ydata $vars(Y_Valve) \
	-color green -symbol "" -mapy y2
    $graph element create p_S -xdata $ts -ydata $vars(p_S) \
	-color blue -symbol ""
    $graph element create qm_S -xdata $ts -ydata $vars(qm_S) \
	-color yellow -symbol ""
    $graph element create sigma_D -xdata $ts -ydata $vars(sigma_D) \
	-color brown -symbol ""
    $graph axis configure x -title "Time"
    $graph axis configure y2 -hide false -title "Y_Valve"
    pack $graph -fill both -expand true
    wm deiconify .	;# bring window to front
    update idletasks
  }
}

solve
