#
# Crane optimal control problem
#  (graphical output with BLT 2)
#
# rf, 1/16/97
# <e-arnold> 2014-01-30 try rbc, if BLT not available
#

source omu.tcl

if {[catch {package present Tk}] 
    || ([catch {package require BLT}] ?
	[catch {package require rbc}] : 0)} {
    set plots 0
    puts stderr "Crane plots disabled"
} else {
    set plots 1
    if {![catch {package present BLT}]} {
	namespace import blt::graph
    } else {
	namespace import rbc::graph
	proc Blt_ZoomStack {g} {Rbc_ZoomStack $g}
	proc Blt_Crosshairs {g} {Rbc_Crosshairs $g}
	proc Blt_ActiveLegend {g} {Rbc_ActiveLegend $g}
	proc Blt_ClosestPoint {g} {Rbc_ClosestPoint $g}
    }
    option add *symbol {}
}

proc plot_vars {} {
  set tscale [lindex [prg_x] 0]
  omu_plot .u0 5 $tscale; omu_plot_titles .u0 {u [kN]} {time [s]}
  omu_plot .x1 1 $tscale; omu_plot_titles .x1 {phi [rad]} {time [s]}
  omu_plot .x2 2 $tscale; omu_plot_titles .x2 {omega [rad/s]} {time [s]}
  omu_plot .x3 3 $tscale; omu_plot_titles .x3 {v [m/s]} {time [s]}
  omu_plot .x4 4 $tscale; omu_plot_titles .x4 {s [m]} {time [s]}
  update
}

prg_name Crane

## apply the bounds to the initial solution in each stage
#prg_bound_init 1

prg_tf_guess 10.0
#prg_tf_guess 15.0

prg_K 50
#prg_K 10
#prg_K 250

prg_setup
prg_simulate

if $plots plot_vars

## configure and run the solver
qp_mat_solver LQDOCP

sqp_init
catch hqp_solve result
puts "Result: $result"

if $plots plot_vars

if {$result == "optimal"} {
  set tscale [lindex [prg_x] 0]
  omu_write_plt control.plt $tscale
}

if 0 {
omu_plot_dump .u0 crane_u.eps
omu_plot_dump .x1 crane_phi.eps
omu_plot_dump .x4 crane_s.eps
}
