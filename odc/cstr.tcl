#
# Optimal control of Klatt Engell reactor
#   formulated as discrete-time model with inline integration in Modelica
#
# See: Franke et al: Discrete-time models for control applications with FMI,
#      Proceedings of the 12th International Modelica Conference,
#      May 15-17, 2017, Prague, Czech Republic,
#      DOI 10.3384/ecp17132507

## translate CSTR.mo to CSTR.fmu

source omc.tcl
if {![file exists CSTR.fmu]
    || [file mtime CSTR.mo] > [file mtime CSTR.fmu]} {
    puts "Compiling FMU..."
    compileFMU CSTR.mo
}

## run optimization

puts "Running optimization..."
package require Omuses

## configure problem

# optimization program and model name
prg_name DOCP
mdl_path CSTR.fmu

# setup stages and sample periods per stage
prg_K 150
prg_setup_stages

# model inputs and sample time points
set KK [prg_K]
set dt [expr 20.0]
set us {}
set ts {}
for {set kk 0} {$kk <= $KK} {incr kk} {
    set t [expr $dt*$kk]
    if {$t < 1500} {
        #           QK_flow TF VF_flow
        lappend us {-1113.5 104.9 14.19}
    } elseif {$t < 1700} {
        lappend us "-1113.5 [expr 104.9+5.1*($t-1500)/200] 14.19"
    } elseif {$t < 2500} {
        lappend us {-1113.5 110 14.19}
    } elseif {$t < 2700} {
        lappend us "-1113.5 [expr 110-5.1*($t-2500)/200] 14.19"
    } else {
        lappend us {-1113.5 104.9 14.19}
    }
    lappend ts $t
}
mdl_us $us
prg_ts $ts

# path constraints at all times
#               cA cB TK
mdl_y_ref       {2.14 1.07 112.9}
mdl_y_weight2   {0    1    0}

# nominal values and interpolation order for inputs
#               QK_flow TF VF_flow
mdl_u_nominal 	{1000 100 10}
mdl_u_order 	{0 0 0}

# determine controlled (active) inputs
mdl_u_active 	{    1    0    0}
mdl_u_max 	{    0  Inf  Inf}
mdl_u_min 	{-9000 -Inf -Inf}

#mdl_der_u_weight2 {0.0001 0 0}
mdl_der_u_weight2 {0.001 0 0}

## setup and solve problem
qp_mat_solver LQDOCP

# setup optimization problem
prg_setup

# simulate initial solution
prg_simulate

# initialize optimizer
sqp_init

# drive solution
catch hqp_solve result

# output results
puts "Result   : $result"
puts "Objective: [prg_f]"
puts "Obj-evals: [prg_fbd_evals]"

# plot result
puts "Plotting resulting signals..."
if {[catch {package present Tk}]
    || ([catch {package require BLT}] ?
	[catch {package require rbc}] : 0)} {
    puts "Skipping plot as no blt::graph available."
} else {

    # read results and form vector of cumulative times
    set ts [prg_ts]
    set us [mdl_us]
    set ys [mdl_ys]
    set KK [prg_K]
    set y1 {}
    set y2 {}
    set y3 {}
    set u1 {}
    set u2 {}
    set u3 {}
    for {set kk 0} {$kk <= $KK} {incr kk} {
	lappend y1 [lindex [lindex $ys $kk] 0]
	lappend y2 [lindex [lindex $ys $kk] 1]
	lappend y3 [lindex [lindex $ys $kk] 2]
	lappend y4 [lindex [lindex $ys $kk] 3]
	lappend u1 [lindex [lindex $us $kk] 0]
	lappend u2 [lindex [lindex $us $kk] 1]
	lappend u3 [lindex [lindex $us $kk] 2]
    }
    if {[lindex [mdl_u_order] 0] == 0} {
	set usmooth "step"
    } else {
	set usmooth "linear"
    }

    # create blt::graph in a toplevel window and plot data
    destroy .graph
    if {![catch {package present BLT}]} {
	set graph [blt::graph .graph -title "[prg_name] demo"]
    } else {
	set graph [rbc::graph .graph -title "[prg_name] demo"]
    }
    $graph element create {QK [kJ/h]} -xdata $ts -ydata $u1 -smooth $usmooth \
	-color red -symbol ""
#    $graph element create {TF [°C]} -xdata $ts -ydata $u2 -smooth $usmooth \
#	-color black -symbol ""
#    $graph element create {cA [mol/l]} -xdata $ts -ydata $y1 \
#	-color blue -symbol "" -mapy y2
    $graph element create {cB [mol/l]} -xdata $ts -ydata $y2 \
	-color green -symbol "" -mapy y2
#    $graph element create {TK [°C]} -xdata $ts -ydata $y3 \
#	-color yellow -symbol "" -mapy y2
    $graph axis configure x -title "Time"
    $graph axis configure y -title "Input"
    $graph axis configure y2 -hide false -title "Outputs"
    pack $graph -fill both -expand true
    wm deiconify .	;# bring window to front
}
