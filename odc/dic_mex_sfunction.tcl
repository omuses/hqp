#
# Double integrator (continuous) example
#   providing the model via MEX S-function
#

## compile sfun_dic.c to a MEX S-function
## (use platform independent file startup.m for communication)
if 1 {
    puts -nonewline "Calling Matlab to compile MEX S-function..."
    flush stdout
    set fp [open startup.m w]
    puts $fp "mex sfun_dic.c; exit"
    close $fp
    exec matlab -nosplash -nodesktop
    file delete startup.m
    puts ""
}

puts "Running optimization..."

## configure problem

# optimization program and model name
prg_name SFunctionOpt
mdl_name sfun_dic

# setup stages and sample periods per stage
prg_KK 60
prg_sps 1
prg_setup_stages

# model inputs and sample time points
set KK [prg_KK]
set dt [expr 1.0/$KK]
set us {}
set ts {}
for {set kk 0} {$kk <= $KK} {incr kk} {
    lappend us -2.0
    lappend ts [expr $dt*$kk]
}
mdl_us $us
prg_ts $ts

# initial states
mdl_x0 		{1.0 0.0}

# path constraints at all times
mdl_y_max 	{Inf 0.1}

# constraints at final time
mdl_yf_min 	{-1.0 0.0}
mdl_yf_max 	{-1.0 0.0}

# determine controlled (active) inputs
mdl_u_active 	1

# define optimization objective
mdl_u_weight2 	1.0

## setup and solve problem

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
