#
# some Tcl procedures for Omuses
#
# rf, 1/16/97
#

#
# generate a vector of start/end times of all stages
#
proc omu_k_times {} {

  set ts [prg_ts]
  set k_times {}
  foreach k [prg_ks] {
    lappend k_times [lindex $ts $k]
  }
  return $k_times
}

#
# read an OmSim-like data file into the global array data
#
proc omu_read_plt {fname {tstart all} {tend all} {dtmin 0.0}} {
  global data

  set fp [open $fname r]

  gets $fp line
  set ndata [lindex $line 2]
  set names {}
  for {set i 0} {$i < $ndata} {incr i} {
      gets $fp line
      lappend names $line
      set data($line) {}
      set data(zeros) {}
  }
  set tprev first
  set lastpos -1
  while {1} {
      gets $fp line
      if [eof $fp] break
      set time [lindex $line 0]
      if {$tstart != "all" && $time < $tstart} continue
      if {$tend != "all" && $time > $tend} break
      if {$time == $tprev} {
	  # replace the last point until the time increases
	  for {set i 0} {$i < $ndata} {incr i} {
	      set data([lindex $names $i]) \
		  [lreplace $data([lindex $names $i]) $lastpos $lastpos \
		   [lindex $line $i]]
	  }
      } elseif {$tprev == "first" || $time >= [expr $tprev+$dtmin]} {
	  # append a new point
	  for {set i 0} {$i < $ndata} {incr i} {
	      lappend data([lindex $names $i]) [lindex $line $i]
	  }
	  incr lastpos
	  set tprev $time
      }
  }
  close $fp
}

#
# write current prg_x into an OmSim-like data file
# (this procedure only works for a constant number of variables per stage!)
#
proc omu_write_plt {fname {tscale 1.0}} {

  set fp [open $fname w]

  set vars [prg_x]
  set k_times [omu_k_times]
  set kdim [llength $k_times]
  set xdim [lindex [prg_nxs] 0]
  set udim [lindex [prg_nus] 0]
  set sdim [expr $xdim+$udim]

  # write the file header
  puts $fp "$kdim 0 [expr $sdim+1]"
  puts -nonewline $fp "time"
  for {set i 0} {$i < $xdim} {incr i} {
    puts -nonewline $fp "\nx$i"
  }
  for {set i 0} {$i < $udim} {incr i} {
    puts -nonewline $fp "\nu$i"
  }

  set idx 0
  foreach kt $k_times {
    puts -nonewline $fp "\n[expr $tscale*$kt] "
    puts -nonewline $fp [lrange $vars $idx [expr $idx+$sdim-1]]
    incr idx $sdim
  }

  # rewrite the controls of the last stage at the final time
  incr idx [expr -$sdim-$udim]
  puts $fp " [lrange $vars $idx [expr $idx+$udim-1]]"

  close $fp
}

#
# plot results into a BLT graph
# (this procedure only works for a constant number of variables per stage!)
#
proc omu_plot {w sidx {tscale 1.0}} {

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
     -title $xtitle
    $w.g yaxis configure \
     -loose true \
     -title $ytitle
    $w.g legend configure \
     -mapped false

    pack $w.g -fill both -expand true
  }

  set vars [prg_x]
  set k_times [omu_k_times]
  set sdim [lindex [prg_nxs] 0]
  if {$sidx >= $sdim} {
    set control 1
  } else {
    set control 0
  }
  incr sdim [lindex [prg_nus] 0]
  set k 0
  set K [llength $k_times]; incr K -1
  set xdata {}
  set ydata {}
  set idx $sidx
  set dim [llength $vars]
  set t [expr $tscale*[lindex $k_times 0]]
  while {$idx < $dim} {
      set var [lindex $vars $idx]
      lappend xdata $t
      lappend ydata $var
      if {$idx >= $dim} break
      incr k
      if {$k > $K} break
      set t [expr $tscale*[lindex $k_times $k]]
      if {$control} {
	  # piecewise constant
	  lappend xdata $t
	  lappend ydata $var
      }
      incr idx $sdim
  }

  set nitems [llength [$w.g element names]]
  $w.g element create g$nitems \
    -xdata $xdata \
    -ydata $ydata

  incr nitems -1
  if {$nitems >= 0} {
    $w.g element configure g$nitems \
      -color gray
  }
}

proc omu_plot_titles {w ytitle xtitle} {

  $w.g yaxis configure \
    -title $ytitle
  $w.g xaxis configure \
    -title $xtitle
}

proc omu_plot_dump {w fname} {

  update
  $w.g postscript config -padx 0 -pady 0
  $w.g postscript output $fname -decorations no
}

