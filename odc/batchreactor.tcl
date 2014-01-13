#
# BatchReactor optimal control example
#
# Interesting to investigate:
#  - different settings for prg_K
#  - tradeoff between formulation effort and solution time
#    for high-level ADOL-C, compared to low-level bare version
#
# rf, 1/13/14
#

foreach example {BatchReactor BatchReactor_bare} {
  puts "\n$example:"

  ## configure example
  prg_name $example
  prg_K 160
  prg_kinf 0.5

  ## configure program formulation
  prg_integrator IMP 	;# implicit midpoint rule
  prg_int_stepsize Inf 	;# one fixed step per interval K

  ## setup program and generate initial solution
  prg_setup
  prg_simulate

  ## setup and initialize solver
  sqp_qp_solver Mehrotra
  sqp_init

  ## run solver
  set timespec [time {catch hqp_solve result}]
  puts "Result: $result ([expr [lindex $timespec 0]/1000] milliseconds)"
}
