#
# hqp_cute.tcl -- HQP (1.5) configuration and execution for CUTE problems
#
# Copyright(C) 1995-98 Ruediger Franke
#
# uncomment lines or insert statements if appropriate
#

#
# print Meschach version on stdout
#

#m_version

#
# set parameters for a CUTE problem
#

proc hqp_config {} {

  # toggle comments to apply stretching
  prg_name CUTE
#  prg_name CUTE_ST

  # uncomment for analytical Lagrangian Hessian
#  prg_hela 1

  # toggle comments for the SQP algorithm
  sqp_solver Powell
#  sqp_solver Schittkowski

  sqp_max_iters 9999
  qp_max_iters   999

  # toggle comments for the matrix solver
  qp_mat_solver RedSpBKP
#  qp_mat_solver SpBKP

  if [prg_hela] {
    sqp_hela Gerschgorin
  } elseif {[prg_name] == "CUTE_ST"} {
    sqp_hela BFGS
  } else {
    sqp_hela DScale
  }
}

#
# you should not need to edit anything below this line
#

# read env(HQP_OPTIONS)
#
proc hqp_options {{stream stdout}} {
  global env

  set noptions 0

  if {![info exists env(HQP_OPTIONS)]} {
    return 0
  }
  
  puts $stream " HQP_OPTIONS"
  foreach option [split $env(HQP_OPTIONS) ","] {
    if {$noptions > 0} {
       puts $stream ","
    }
    if [regexp {(.*)\=(.*)} $option all name value] {
      if {![catch {$name $value} reason]} {
        puts -nonewline $stream "  $name = $value"

        # adapt prg_hela automatically to sqp_hela
        if {$name == "sqp_hela"} {
	  if {[sqp_hela] == "Gerschgorin"} {
            prg_hela 1
          } else {
            prg_hela 0
          }
          puts -nonewline $stream ", prg_hela = [prg_hela]"
        }
      } else {
        puts -nonewline $stream "  $reason"
      }
    } else {
      puts -nonewline $stream "  $option: not understood"
    }
    incr noptions
  }
  puts $stream "\n"

  return $noptions
}

# exit HQP
#
proc hqp_exit {reason} {

  puts "HQP [hqp_version]: $reason"

  if {[sqp_iter] > 0} {
    catch prg_f f
    puts "f      : $f"
    puts "f_evals: [prg_fbd_evals]"
    prg_write_soln "$reason"
  } else {
    prg_write_soln -nosoln "$reason"
  }      

  exit 0
}

#
# configure HQP
# read environment variable HQP_OPTIONS
#

hqp_config
hqp_options

#
# solve the problem
#

prg_setup
sqp_init

#
# some statistics
#
catch {
  puts [format " Stages/Groups: %6d / %6d" [prg_K] [prg_NEL]]
}
catch {
  puts [format " HQP n/ me/ m : %6d / %6d / %6d" \
        [llength [prg_x]] [llength [sqp_y]] [llength [sqp_z]] ]
}
catch {
  puts [format " Max stagesize: %6d" [sqp_hela_max_bsize]]
}
catch {
  puts [format " Semibandwidth: %6d" [mat_sbw]]
}
puts ""

#set fail [catch {hqp_solve [open /dev/null w]} reason]
set fail [catch hqp_solve reason]

#
# write the solution and exit
#

if {$fail} {
  set source [lindex $errorInfo 3]
  if {[lindex [lindex $source 0] 0] != "error"} {
    hqp_exit "$reason from $source"
  } else {
    hqp_exit $reason
  }
}

hqp_exit "optimal solution"
