#
# OpenModelica Compiler utility routines
#

# configuration settings -- may add installation path
set omc {omc -d=-disableDirectionalDerivatives --simCodeTarget=Cpp --exportClocksInModelDescription}

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
    puts $fp "setCommandLineOptions(\"--std=3.3\");"
    puts $fp "setCommandLineOptions(\"+simCodeTarget=Cpp\");"
    puts $fp "setCommandLineOptions(\"-d=-disableDirectionalDerivatives\");"
    puts $fp "loadFile(\"$mo_file_normalized\");"
    puts $fp "getErrorString();"
    puts $fp "buildModelFMU($model_name, version=\"2.0\");"
    puts $fp "getErrorString();"
    close $fp
    # call omc, check for errors and clean up
    if {[catch {eval [concat exec $omc "omc_tcl_commands.mos"]} log]} {
        cd $cwd
        error "Could not compile FMU: $log"
    }
    cd $cwd
    file copy -force "tmp/$model_name.fmu" .
    file delete -force tmp
    return [file rootname $mo_file_normalized].fmu
}
