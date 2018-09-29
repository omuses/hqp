#
# Parameter and initial state estimation example
#   for Double Integrator (Discrete-time) given as FMU
#

## translate DID.mo to DID.fmu

source omc.tcl
if {![file exists DID.fmu]
    || [file mtime DID.mo] > [file mtime DID.fmu]} {
    puts "Compiling FMU..."
    compileFMU DID.mo
}

## prepare estimation data set
## (calculate reference outputs from known solution and add noise)

# same parameters for all experiments
set p           {2.0}

# individual initial value for the first state in each experiment
set x01s        {{1.5} {0.5}}

# initialize seed of the random number generator used for noise
expr srand(0)

# first experiment
set x011 [lindex $x01s 0]
set us1 {}
set ts1 {}
set ys_ref1 {}
set K1 10
set u1 1.0
set dt1 [expr 2.0/$K1]
for {set k 0} {$k <= $K1} {incr k} {
    lappend us1 $u1
    lappend ts1 [expr $k*$dt1]
    set y11 [expr $x011+$k*$dt1*$p*$u1]
    lappend ys_ref1 [expr ($y11 + 0.5*(rand()-0.5))]
}

# second experiment
set x012 [lindex $x01s 1]
set us2 {}
set ts2 {}
set ys_ref2 {}
set K2 10
set u2 2.0
set dt2 [expr 1.0/$K2]
for {set k 0} {$k <= $K2} {incr k} {
    lappend us2 $u2
    lappend ts2 [expr $k*$dt2]
    set y12 [expr $x012+$k*$dt2*$p*$u2]
    lappend ys_ref2 [expr ($y12 + 0.5*(rand()-0.5))]
}

## run estimation

puts "Running estimation..."
package require Omuses

## configure problem

# optimization program and model name
prg_name DynamicEst
mdl_path DID.fmu

# setup stages and number of experiments
prg_K [expr $K1+$K2+1]
prg_nex 2
prg_setup_stages

# initial guesses for states for all experiments
mdl_x0s { {1 0} {1 0} }

# model inputs and sample time points
mdl_us [concat $us1 $us2]
prg_ts [concat $ts1 $ts2]

# define parameters and initial states to estimate
mdl_p_active    {1 0 0}
mdl_x0_active   {1 0}

# interpolation order for inputs
mdl_u_order 	{0}

# define outputs used for estimation, i.e. use first output
mdl_y_active    {1 0}

## setup and solve problem

# setup estimation problem
prg_setup

# provide reference outputs
prg_ys_ref [concat $ys_ref1 $ys_ref2]

# simulate initial solution
prg_simulate

puts "\nStart"
puts "Parameters    : [mdl_p]"
puts "Initial states: [mdl_x0s]"
puts "Residuum      : [prg_f]"
puts ""

# initialize optimizer
sqp_init

# drive solution
catch hqp_solve result

# output results
puts "Result   : $result"
puts "Objective: [prg_f]"
puts "Obj-evals: [prg_fbd_evals]"

puts "\nEstimated"
puts "Parameters    : [mdl_p] +- [mdl_p_confidence]"
puts "Initial states: [mdl_x0s] +- [mdl_x0_confidence]"
puts "Residuum      : [prg_f]"
puts ""

# plot signals
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
    set ys_ref [prg_ys_ref]
    set K [prg_K]
    set times {}
    set t_offs 0.0
    set y1 {}
    set y2 {}
    set tk [lindex $ts 0]
    for {set k 0} {$k <= $K} {incr k} {
	set tk_last $tk
	set tk [lindex $ts $k]
	if {$tk < $tk_last} {
	    set t_offs [expr $t_offs+$tk_last-$tk]
	}
	lappend times [expr $t_offs+$tk]
	lappend y1 [lindex [lindex $ys $k] 0]
	lappend y2 [lindex [lindex $ys $k] 1]
    }

    # create blt::graph in a toplevel window and plot data
    set graph .g
    destroy $graph
    if {![catch {package present BLT}]} {
	set graph [blt::graph .graph -title "[prg_name] demo"]
    } else {
	set graph [rbc::graph .graph -title "[prg_name] demo"]
    }
    $graph element create u -xdata $times -ydata $us -smooth step \
        -color red -pixels 3
    $graph element create y1 -xdata $times -ydata $y1 -smooth step \
        -color blue -pixels 3
    $graph element create y2 -xdata $times -ydata $y2 -smooth step \
        -color green -pixels 3
    $graph element create y1_ref -xdata $times -ydata $ys_ref \
	-color blue -linewidth 0 -symbol "scross" -pixels 8
    $graph axis configure x -title "Time (cumulative)"
    $graph axis configure y -title "Signals"
    pack $graph -fill both -expand true
}
