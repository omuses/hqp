## @file fmi.tcl
#  Read fmu (zip) files, extract binaries and parse xml model descriptions.
#  rf, 02/22/2014

package provide fmi 2.0

package require vfs::zip	;# used to read FMU files

## hold one sub-namespace with infos for each loaded FMU
namespace eval ::fmu {
}

## fmi functions
namespace eval ::fmi {
    # translate tcl_platform(os) to directory name of FMU binaries
    array set platformSubDir {
	"Darwin" 	darwin
	"Linux" 	linux
	"Windows NT" 	win
    }
}

## Obtain name of FMU from fmuPath (excluding directories and file extension)
#  @return fmuName
proc ::fmi::getName {fmuPath} {
    return [file rootname [file tail $fmuPath]]
}

## Get location of extracted FMU resources
#  @return URI of resources directory
proc ::fmi::getResourcesURI {fmuPath} {
    set dirPath [file rootname $fmuPath]
    return "file:///[string trimleft [file normalize $dirPath] /]/resources"
}

## Get path of FMU binary from attributes exported by readModelDescription
#  @return binary path for respective platform
proc ::fmi::getBinaryPath {fmuPath} {
    set binPath [file rootname $fmuPath]
    set binPath ${binPath}/binaries/$::fmi::platformSubDir($::tcl_platform(os))
    set binPath ${binPath}[expr 8*$::tcl_platform(wordSize)]
    set fmuName [file rootname [file tail $fmuPath]]
    set modelIdentifier [set ::fmu::${fmuName}::attributes(modelIdentifier)]
    return ${binPath}/${modelIdentifier}[info sharedlibextension]
}

## Test if folder with rootname of FMU exists and has same or newer mtime
#  @return true in case of success
proc ::fmi::testExtracted {fmuPath} {
    set dirPath [file rootname $fmuPath]
    if {[file exists $dirPath] &&
        [file mtime $dirPath] >= [file mtime $fmuPath]} {
        return true
    } else {
        return false
    }
}

## Extract all files of FMU, including e.g. binaries and resources
proc ::fmi::extractModel {fmuPath} {
    set fmuTime [file mtime $fmuPath]
    set dirPath [file rootname $fmuPath]
    if [file exists $dirPath] {
	file delete -force $dirPath
    }
    set fp [::vfs::zip::Mount $fmuPath $fmuPath]
    file copy $fmuPath $dirPath
    ::vfs::zip::Unmount $fp $fmuPath
    file mtime $dirPath $fmuTime
}

## Parse modelDescription.xml of extracted model and
#  store infos in namespace array ::fmu::${fmuName}
proc ::fmi::readModelDescription {fmuPath} {

    # read modelDescription.xml
    set dirPath [file rootname $fmuPath]
    set fp [open $dirPath/modelDescription.xml]
    fconfigure $fp -encoding utf-8
    set fmiModelDescription [xml2list [read $fp]]
    close $fp

    # create a namespace to hold infos
    set fmuName [file rootname [file tail $fmuPath]]
    if [namespace exists ::fmu::$fmuName] {
	namespace delete ::fmu::$fmuName
    }
    namespace eval ::fmu::$fmuName {
    }

    # get model attributes, 
    # e.g. fmuAttributes(fmiVersion) or fmuAttributes(guid)
    array set fmuAttributes [lindex $fmiModelDescription 1]
    set fmiVersion $fmuAttributes(fmiVersion)

    # default values for optional elements of model description
    set fmiElements(TypeDefinitions) {}
    # get model description, e.g. fmiElements(ModelVariables)
    foreach element [lindex $fmiModelDescription 2] {
	set fmiAttributes([lindex $element 0]) [lindex $element 1]
	set fmiElements([lindex $element 0]) [lindex $element 2]
    }

    # complement fmuAttributes for ModelExchange
    # e.g. fmuAttributes(modelIdentifier)
    if {$fmiVersion >= 2.0} {
	if {[array names fmiAttributes ModelExchange] == {}} {
	    error "Missing tag ModelExchange in modelDescription" 
	}
	array set fmuAttributes $fmiAttributes(ModelExchange)
    }

    # collect type definitions
    array unset typeDefinitions
    foreach typeDefinition $fmiElements(TypeDefinitions) {
	set typeName [lindex [lindex $typeDefinition 1] 1]
	set typeDefinitions($typeName) [lindex [lindex $typeDefinition 2] 0]
    }

    # identify states
    set stateIndices {}
    if {$fmiVersion < 2.0} {
	if {$fmuAttributes(numberOfContinuousStates) > 0} {
	    error "States are not supported in FMI $fmiVersion"
	}
    } else {
	# find states through ModelStructure and derivative attribute
	foreach element $fmiElements(ModelStructure) {
	    if {[lindex $element 0] == "Derivatives"} {
		foreach derivative [lindex $element 2] {
		    set attributes [lindex $derivative 1]
		    set indexIdx [lsearch $attributes "index"]
		    set dIdx [expr [lindex $attributes [incr indexIdx]]-1]
		    set dVar [lindex $fmiElements(ModelVariables) $dIdx]
		    set typeAttributes [lindex [lindex [lindex $dVar 2] 0] 1]
		    set sIdx [lsearch $typeAttributes "derivative"]
		    lappend stateIndices [lindex $typeAttributes [incr sIdx]]
		}
	    }
	}
    }

    # collect variables per category
    foreach category {"parameter" "state" "input" "output"} {
	set ::fmu::${fmuName}::${category}Names {}
	set ::fmu::${fmuName}::${category}Indices {}
	set ::fmu::${fmuName}::${category}References {}
	set ::fmu::${fmuName}::${category}BaseTypes {}
	set ::fmu::${fmuName}::${category}NominalValues {}
	set ::fmu::${fmuName}::${category}StartValues {}
    }
    set index 1
    foreach modelVariable $fmiElements(ModelVariables) {
	array unset v
	set v(name) ""
	set v(valueReference) ""
	set v(causality) ""
	set v(variability) ""
	set v(start) "NaN"
	set v(fixed) ""
	set v(nominal) "1"
	# fetch attributes from model variable
	array set v [lindex $modelVariable 1]
	# complement with attributes from declared variable type
	set vBaseType [lindex [lindex [lindex $modelVariable 2] 0] 0]
	set vTypeAttributes [lindex [lindex [lindex $modelVariable 2] 0] 1]
	set idx [lsearch $vTypeAttributes declaredType]
	if {$idx >= 0} {
	    set declaredType [lindex $vTypeAttributes [incr idx]]
	    array set v [lindex $typeDefinitions($declaredType) 1]
	}
	# complement/override with attributes from variable type
	array set v $vTypeAttributes

	# detect parameters in FMI 1.0
	# Note: use fixed start value as indicator for unbound parameter
	if {$fmiVersion < 2.0} {
	    if {$v(variability) == "parameter"
		&& $v(start) != "NaN" && $v(fixed) == "true"} {
		set v(causality) "parameter"
	    }
	}

	# obtain category (state, parameter, input/output)
	set category ""
	if {[lsearch $stateIndices $index] >= 0} {
	    set category "state"
	} elseif {[regexp ^(parameter|input|output)$ $v(causality)]} {
	    set category $v(causality)
	}
	# store variable infos
	if {$category != ""} {
	    lappend ::fmu::${fmuName}::${category}Names $v(name)
	    lappend ::fmu::${fmuName}::${category}Indices $index
	    lappend ::fmu::${fmuName}::${category}References $v(valueReference)
	    lappend ::fmu::${fmuName}::${category}BaseTypes $vBaseType
	    lappend ::fmu::${fmuName}::${category}NominalValues $v(nominal)
	    lappend ::fmu::${fmuName}::${category}StartValues $v(start)
	}

	# keep mapping from valueReference to name per base type
	set typeId [string tolower [string index $vBaseType 0]]
	set ::fmu::${fmuName}::${typeId}Names($v(valueReference)) $v(name)

	incr index
    }

    # get model structure
    # use zero based indices counted per category
    set catIdx 0
    foreach category {parameter state input output} {
        foreach index [set ::fmu::${fmuName}::${category}Indices] {
            set mapIndex($index) $catIdx
            incr catIdx
        }
    }
    foreach element $fmiElements(ModelStructure) {
	set structureElements([lindex $element 0]) [lindex $element 2]
    }
    foreach category {output derivative discreteState} {
        set elementName [string replace $category 0 0 \
                             [string toupper [string index $category 0]]]s
        foreach element $structureElements($elementName) {
            array set structureAttributes [lindex $element 1]
            set var $mapIndex($structureAttributes(index))
            set dependencies {}
            foreach dependency $structureAttributes(dependencies) {
                lappend dependencies $mapIndex($dependency)
            }
            set ::fmu::${fmuName}::${category}Dependencies($var) $dependencies
        }
    }

    # export fmuAttributes
    array set ::fmu::${fmuName}::attributes [array get fmuAttributes]
}

## Replace valueReferences of the form #r123# with variable names
#  @return processed message
proc ::fmi::mapNames {fmuName message} {
    set i [string first \# $message]
    while {$i >= 0} {
	set j [string first \# $message [expr $i+1]]
	if {$j < 0} {
	    break
	} elseif {[expr $j-$i] == 1} {
	    set message [string replace $message $i $j \#]
	} else {
	    set ref [string trim [string range $message $i $j] \#]
	    set typeId [string index $ref 0]
	    set valueReference [string range $ref 1 end]
	    set aName "::fmu::${fmuName}::${typeId}Names"
	    if {[lsearch [array names $aName] $valueReference] >= 0} {
		set name [lindex [array get $aName $valueReference] 1]
		set message [string replace $message $i $j $name]
	    }
	}
	set i [string first \# $message [incr i]]
    }
    return $message
}

## Parse XML to a Tcl list.
# See: A little XML parser, http://wiki.tcl.tk/3919
#      by Richard Suchenwirth 2002-08-20
#
# Converts between XML and well-formed Tcl list with:
# - if the first element is #text, the second is the text content;
# - else the first element is the tag, 
#        the second the attributes alternating name value..., and 
#        the third is a list of the child elements.
#
# Equivalent to tDOM $node asList command.
proc ::fmi::xml2list xml {
    regsub -all {<\?.*?\?>} $xml "" xml  ;# remove processing instructions
    regsub -all {<!--.*?-->} $xml "" xml ;# remove comments
    regsub -all {>\s*<} [string trim $xml " \n\t<>"] "\} \{" xml
    set xml [string map {> "\} \{#text \{" < "\}\} \{"}  $xml]

    set res ""   ;# string to collect the result   
    set stack {} ;# track open tags
    set rest {}
    foreach item "{$xml}" {
	switch -regexp -- $item {
	    ^# {append res "{[lrange $item 0 end]} " ;# text item
	    }
	    ^/ {
		regexp {/(.+)} $item -> tagname ;# end tag
		set expected [lindex $stack end]
		if {$tagname!=$expected} {error "$item != $expected"}
		set stack [lrange $stack 0 end-1]
		append res "\}\} "
	    }
	    /$ { # singleton - start and end in one <> group
		regexp {([^ ]+)( (.+))?/$} $item -> tagname - rest
		set rest [lrange [string map {= " "} $rest] 0 end]
		append res "{$tagname [list $rest] {}} "
	    }
	    default {
		set tagname [lindex $item 0] ;# start tag
		set rest [lrange [string map {= " "} $item] 1 end]
		lappend stack $tagname
		append res "\{$tagname [list $rest] \{"
	    }
	}
	if {[llength $rest]%2} {error "att's not paired: $rest"}
    }
    if [llength $stack] {error "unresolved: $stack"}
    string map {"\} \}" "\}\}"} [lindex $res 0]
}

## Convert a Tcl list to XML
proc ::fmi::list2xml list {
    switch -- [llength $list] {
	2 {lindex $list 1}
	3 {
	    foreach {tag attributes children} $list break
	    set res <$tag
	    foreach {name value} $attributes {
		append res " $name=\"$value\""
	    }
	    if [llength $children] {
		append res >
		foreach child $children {
		    append res [list2xml $child]
		}
		append res </$tag>
	    } else {append res />}
	}
	default {error "could not parse $list"}
    }
}
