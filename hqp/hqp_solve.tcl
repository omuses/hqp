#
# some Tcl procedures for the Hqp package
#
# rf, 9/8/96
#
# rf, 4/6/00
#  - added hqp_solve_hot for hot start of SQP solver
#

proc hqp_exit {reason} {
  exit 0
}

#
# procedure for defining hqp_logging option
#
proc hqp_logging {{value {}}} {
  global hqp_logging
  # default value
  if {![info exists hqp_logging]} {
    set hqp_logging 1
  }
  # process value argument
  if {$value != {}} {
    set hqp_logging $value
  }
  if {$hqp_logging} {
    proc hqp_puts {args} {eval puts $args}
    proc hqp_flush {args} {eval flush $args}
  } else {
    proc hqp_puts {args} {}
    proc hqp_flush {args} {}
  }
  return $hqp_logging
}
# initialize default logging
hqp_logging

#
# solver initialization and iteration counter
#

proc hqp_init {} {
  if {[sqp_solver] != "None"} {
    sqp_init
  } elseif {[mip_solver] != "None"} {
    # mip_init
  } else {
    error "No solver configured. Expecting sqp_solver or mip_solver."
  }
}

proc hqp_iter {{value {}}} {
  if {[sqp_solver] != "None"} {
    return [eval sqp_iter $value]
  } elseif {[mip_solver] != "None"} {
    return 0
  } else {
    error "No solver configured. Expecting sqp_solver or mip_solver."
  }
}

proc hqp_max_iters {{value {}}} {
  if {[sqp_solver] != "None"} {
    return [eval sqp_max_iters $value]
  } elseif {[mip_solver] != "None"} {
    return 1
  } else {
    error "No solver configured. Expecting sqp_solver or mip_solver."
  }
}

#
# hot start routine
#
proc hqp_solve_hot {{stream stdout}} {
  hqp_solve $stream 1
}

#
# general solution routine
#
proc hqp_solve {{stream stdout} {hot 0}} {

  if {[sqp_solver] == "None" && [mip_solver] == "None"} {
      error "No solver configured; require sqp_solver or mip_solver."
  }

  #
  # Call SQP solver if it has been configured
  #

 if {[sqp_solver] != "None"} {

  global qp_iters

  if {![info exists qp_iters] || [sqp_iter] == 0} {
    set qp_iters 0
  }
  set nullsteps 0

  if {[sqp_iter] == 0 || $hot} {
    hqp_puts $stream [format "%3s %12s %10s %10s \[%3s %3.3s\] %10s %10s %8s"\
                  it obj ||inf|| ||grdL||\
                  qp result ||s|| s'Qs stepsize]
  }

  while 1 {

    if {$qp_iters == 0} {
      if {!$hot} {
	# can't re-use higher order information
	sqp_qp_update
      }
      hqp_puts -nonewline $stream [format "%3d %12.6g %10.4g %10.4g " \
         [sqp_iter] [prg_f] [sqp_norm_inf] [sqp_norm_grd_L]]
      hqp_flush $stream
      # check feasibility of model evaluation
      if [catch {expr [prg_f]*[sqp_norm_inf]}] {
        hqp_puts $stream ""
        error evaluation
      }
    } else {
      hqp_puts -nonewline $stream [format "%3d %12.6g %10.4g " \
         [sqp_iter] [prg_f] [sqp_norm_inf]]
      hqp_flush $stream
      # check feasibility of model evaluation
      if [catch {expr [prg_f]*[sqp_norm_inf]}] {
        hqp_puts $stream ""
        error evaluation
      }

      # additional break test as norm_inf changed with taken step
      # (Don't use this brake test after cold start,
      #  as Hessian approximation is not updated anymore!
      #  This update is useful for subsequent hot starts.)
      if {$hot && [sqp_iter] > 0 && [sqp_sQs] >= 0.0 && !$hela_restart} {
        if {[sqp_norm_inf] < [sqp_eps] && [qp_result] == "optimal"} {
          if {[sqp_sQs] < [expr [sqp_eps]*[sqp_eps]]} break
 	  if {[sqp_iter] > 2} {
            if {[sqp_norm_s] < [expr [sqp_eps]*[sqp_norm_x]]
 	        && [sqp_norm_df] < [expr [sqp_eps]*abs([prg_f])]
	        && [sqp_sQs] < [sqp_eps]} break
  	  }
        }
      }

      sqp_qp_update

      hqp_puts -nonewline $stream [format "%10.4g " [sqp_norm_grd_L]]
      hqp_flush $stream
    }

    # test for failed hot start
    if {0 && $hot} {
      if {$qp_iters == 0} {
        set norm_grd_L0 [sqp_norm_grd_L]
      } elseif {[sqp_norm_grd_L] > 1.0 &&
	        [expr [sqp_eps]*[sqp_norm_grd_L]] >= $norm_grd_L0} {
        hqp_puts "\nRestart cold:"
        prg_setup
        prg_simulate
        sqp_init
        return [hqp_solve $stream]
      }
    }
    
    # test for bad Hessian
    if {[sqp_xQx] < 0.0} {
      sqp_hela_restart
      set hela_restart 1
    } else {
      set hela_restart 0
    }
    set xQx_old [sqp_xQx]

    if {[sqp_iter] > 0} {
      if {[sqp_norm_inf] < [sqp_eps]} {
        if {[sqp_norm_grd_L] < [sqp_eps]} break
      }
    }

    sqp_qp_solve
    incr qp_iters [qp_iter]

    hqp_puts -nonewline $stream [format "\[%3d %3.3s\] " \
	[qp_iter] [qp_result]]
    hqp_flush $stream

    if {[qp_iter] == 0} {
	hqp_puts $stream ""
	error [qp_result]
    }

    hqp_puts -nonewline $stream [format "%10.4g %10.4g " \
        [sqp_norm_s] [sqp_sQs]]
    hqp_flush $stream

    if {[sqp_sQs] < 0.0} {
      sqp_hela_restart
    }

    if {[sqp_iter] > 0 && [sqp_sQs] >= 0.0 && !$hela_restart} {
      if {[sqp_norm_inf] < [sqp_eps] && [qp_result] == "optimal"} {
        if {[sqp_sQs] < [expr [sqp_eps]*[sqp_eps]]} break
 	if {[sqp_iter] > 2} {
          if {[sqp_norm_s] < [expr [sqp_eps]*[sqp_norm_x]]
 	      && [sqp_norm_df] < [expr [sqp_eps]*abs([prg_f])]
	      && [sqp_sQs] < [sqp_eps]} break
	}
      }
    }

    sqp_step

    hqp_puts $stream [format "%8.3g" [sqp_alpha]]

    if {[qp_iter] >= [qp_max_iters] && [qp_result] != "feasible"} {
        error subiters
    }
    if {[sqp_iter] >= [sqp_max_iters]} {
        error iters
    }
    if {[sqp_inf_iters] >= [sqp_max_inf_iters]} {
	if {[qp_result] == "suboptimal"} {
	    error infeasible
	} else {
	    error degenerate
	}
    }

    if {[sqp_alpha] < 1e-8
        && [sqp_norm_df] < [expr [sqp_eps]*abs([prg_f])]} {
	incr nullsteps 1
    } else {
	set nullsteps 0
    }
    if {$nullsteps > 5} {
	error stall
    }

    # update quasi-parallel running Tcl stuff
    catch {update}
  }

  hqp_puts $stream [format "\n%43d qp-it" $qp_iters]
  unset qp_iters

 }

  #
  # Call mixed integer solver if it has been configured
  # Note: both SQP and MIP solver may be useful
  #       to first solve for non-linearities and then
  #       for integers using the last linearized problem
  #

  if {[mip_solver] != "None"} {
    mip_solve
    prg_update_fbd 	;# take over solution
    return [mip_result]
  }

  return optimal
}
