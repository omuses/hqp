#
# Bio-technological process optimal control problem
# graphical output with BLT
#
# E. Arnold   2003-02-16
#

source omu.tcl

if { [catch {set tk_version}] || [catch {package require BLT}] } {
    set plots 0
    puts stderr "plots disabled"
} else {
    set plots 1
    namespace import blt::*
    option add *symbol {}
}

#------------------------------------------------------------------------------
#
# plot results prg_ts, prg_y into a BLT graph
#
proc omu_plot_y {w sidx {tscale 1.0}} {

    if {![winfo exists $w]} {
	toplevel $w
	wm minsize $w 10 10
	set xtitle time
	set ytitle ""
	regexp {.*\.([^.]*)} $w all ytitle
	graph $w.g \
		-title {} \
		-plotborderwidth 0 \
		-background white \
		-plotbackground white \
		-width 500 \
		-height 200
	$w.g xaxis configure \
		-background white \
		-title $xtitle
	$w.g yaxis configure \
		-background white \
		-loose true \
		-title $ytitle
	$w.g legend configure \
		-background white \
		-hide true

	pack $w.g -fill both -expand true

	Blt_ZoomStack $w.g
	Blt_Crosshairs $w.g
	Blt_ActiveLegend $w.g
	Blt_ClosestPoint $w.g
    }

    set vars [prg_y]
    set k 0
    set idx $sidx
    set sdim [expr [llength $vars]/[llength [prg_ts]]]
    set xdata {}
    set ydata {}
    while { $k < [llength [prg_ts]] } {
	lappend xdata [expr $tscale*[lindex [prg_ts] $k]]
	lappend ydata [lindex $vars $idx]
	incr k
	incr idx $sdim
    }

    set nitems [llength [$w.g element names]]
    $w.g element create g$nitems \
	    -symbol {} \
	    -xdata $xdata \
	    -ydata $ydata

    incr nitems -1
    if {$nitems >= 0} {
	$w.g element configure g$nitems \
		-color gray
    }
}

proc plot_vars {} {
    set tscale 1.0
    omu_plot .u1 2;   omu_plot_titles .u1 {Fs [L/h]} {time [h]}
    omu_plot .x1 0;   omu_plot_titles .x1 {p [g]} {time [h]}
    omu_plot .x2 1;   omu_plot_titles .x2 {INT(Fs) [L]} {time [h]}
    omu_plot_y .y1 0; omu_plot_titles .y1 {cs [g/L]} {time [h]}
    omu_plot_y .y2 1; omu_plot_titles .y2 {cp [g/L]} {time [h]}
    update
}

prg_name Bio

# change model parameters

# number of stages
#prg_K 151
prg_K 51

# initial control trajectory generation: 0 - constant (prg_uinit), 
#                                        1 - nonlinear controller
prg_controller 0
#prg_uinit 0.02

# initial substrate concentration
#prg_cs0 20.0

# optimization horizon
#prg_tf 150
prg_tf 10

# initialize 
prg_setup
prg_simulate

if { $plots } {
    plot_vars
}

# configure and run the solver
qp_mat_solver LQDOCP

sqp_init
catch hqp_solve result
puts "Result: $result"

if { $plots } {
    plot_vars
}

if { $result == "optimal" } {
    omu_write_plt control.plt 
}

if { 0 } {
    omu_plot_dump .u1 bio_fs.eps
    omu_plot_dump .x1 bio_p.eps
    omu_plot_dump .y1 bio_cs.eps
    omu_plot_dump .y2 bio_cp.eps
}
