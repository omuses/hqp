#
# OpenModelica Compiler utility routines
#

# configuration settings -- may add installation path
set omc {omc}

# load Modelica file and translate contained model to FMU
proc compileFMU {mo_file} {
    global omc
    set model_name [file rootname $mo_file]
    set mo_file_normalized [file normalize $mo_file]
    # create a temporary directory
    set cwd [pwd]
    if [file exists tmp] {
	file delete -force tmp
    }
    file mkdir tmp
    cd tmp
    # create script file for omc
    set fp [open "omc_tcl_commands.mos" w]
    puts $fp "// temporary file with omc commands"
    puts $fp "loadFile(\"$mo_file_normalized\");"
    puts $fp "getErrorString();"
    puts $fp "translateModelFMU($model_name, version=\"2.0\");"
    puts $fp "getErrorString();"
    close $fp
    # call omc, check for errors and clean up
    if {[catch {eval [concat exec $omc "omc_tcl_commands.mos"]} log]} {
	cd $cwd
	error "Could not compile FMU: $log"
    }
    cd $cwd
    set logLines [split $log "\n"]
    set loadFile [lindex $logLines 1]
    if {$loadFile != {""}} {
	error $loadFile
    }
    set translateModelFMULines [lrange $logLines 3 end]
    set translateModelFMU [join $translateModelFMULines "\n"]
    foreach logLine $translateModelFMULines {
	if {[string compare -length 6 "Error:" $logLine] == 0} {
	    error $translateModelFMU
	}
    }
    if {$translateModelFMU != {""}} {
	puts $translateModelFMU
    }
    file copy -force "tmp/$model_name.fmu" .
    file delete -force tmp
    return [lindex $log 2]
}
