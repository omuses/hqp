#
# demo for Schittkowski's differentiable penalty function
#

prg_name Maratos

## toggle comments
#sqp_solver Powell
sqp_solver Schittkowski

prg_setup
sqp_init
catch hqp_solve result

puts "Result   : $result"
puts "Objective: [prg_f]"
puts "Obj-evals: [prg_fbd_evals]"
puts "Variables: [prg_x]"
exit
