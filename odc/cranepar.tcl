#
# Crane parameter and initial states estimation
#
# rf, 1/16/97
#

source omu.tcl

proc hqp_exit {reason} {
  puts [prg_x]
  exit 0
}

prg_name CranePar

## uncomment to change to single stage
#prg_multistage 0

## uncomment to switch off automatic differentiation
#prg_ad 0

## read the data file and setup the stages
omu_read_plt record.plt
prg_KK [expr {[llength $data(Base::Time)] - 1}]
prg_setup_stages

## initialize the measurement data and the sample time points for the
## position trajectory [m]; max. deviation to add (+/-) [m]; seed value
prg_ts     $data(Base::Time)
prg_s_ref  $data(Crane.s.value)
prg_maxdev 0.05
prg_seed   1234
prg_disturb

## optimizer stopping limit
sqp_eps 1e-7

## run the optimizer
prg_setup
prg_simulate
sqp_init
catch hqp_solve result

## print the results
puts "Result:         $result"
puts "Obj-evals:      [prg_fbd_evals]"
if {$result == "optimal"} {
  set x [prg_x]
  puts "Estimated mass: [expr 1e3*[lindex $x 0]]"
  puts "Initial states: [lrange $x 1 4]"
} else {
  hqp_exit $result
}

exit
