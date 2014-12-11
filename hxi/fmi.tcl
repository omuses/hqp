## @file fmi.tcl
#  Read fmu (zip) files, extract binaries and parse xml model descriptions.
#  rf, 02/22/2014

package provide fmi 2.0

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

## Obtain directory name of extracted FMU from fmuPath
#  @return dirPath
proc ::fmi::getDirPath {fmuPath} {
    return [file dirname $fmuPath]/.[file rootname [file tail $fmuPath]]
}

## Get location of extracted FMU resources
#  @return URI of resources directory
proc ::fmi::getResourcesURI {fmuPath} {
    set dirPath [::fmi::getDirPath $fmuPath]
    return "file:///[string trimleft [file normalize $dirPath] /]/resources"
}

## Get path of FMU binary from attributes exported by readModelDescription
#  @return binary path for respective platform
proc ::fmi::getBinaryPath {fmuPath} {
    set dirPath [::fmi::getDirPath $fmuPath]
    set binPath ${dirPath}/binaries/$::fmi::platformSubDir($::tcl_platform(os))
    set binPath ${binPath}[expr 8*$::tcl_platform(wordSize)]
    set fmuName [file rootname [file tail $fmuPath]]
    set modelIdentifier [set ::fmu::${fmuName}::attributes(modelIdentifier)]
    return ${binPath}/${modelIdentifier}[info sharedlibextension]
}

## Test if folder with rootname of FMU exists and has same or newer mtime
#  @return true in case of success
proc ::fmi::testExtracted {fmuPath} {
    set dirPath [::fmi::getDirPath $fmuPath]
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
    set dirPath [::fmi::getDirPath $fmuPath]
    if [file exists $dirPath] {
	file delete -force $dirPath
    }
    if {[package vcompare [info tclversion] "8.6"] >= 0} {
	# use Tcl's built-in zlib
	::fmi::unzip $fmuPath
    } else {
	package require vfs::zip
	set fp [::vfs::zip::Mount $fmuPath $fmuPath]
	file copy $fmuPath $dirPath
	::vfs::zip::Unmount $fp $fmuPath
    }
    file mtime $dirPath $fmuTime
}

## Parse modelDescription.xml of extracted model and
#  store infos in namespace array ::fmu::${fmuName}
#  @return fmiVersion string
proc ::fmi::readModelDescription {fmuPath} {

    # extract model if not yet done
    if {![::fmi::testExtracted $fmuPath]} {
        ::fmi::extractModel $fmuPath
    }

    # read modelDescription.xml
    set dirPath [::fmi::getDirPath $fmuPath]
    set fp [open $dirPath/modelDescription.xml]
    fconfigure $fp -encoding utf-8
    set fmiModelDescription [xml2list [read $fp]]
    close $fp

    # create a namespace for parsed model description
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
    set derivativeIndices {}
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
		    lappend derivativeIndices [incr dIdx]
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

	# obtain categories (state, parameter, input/output)
	# note that a variable may be state and output
	set categories {}
	if {[lsearch $stateIndices $index] >= 0} {
	    lappend categories "state"
	}
	if {[regexp ^(parameter|input|output)$ $v(causality)]} {
	    lappend categories $v(causality)
	}
        # unify boolean literals to 1 and 0 for easier treatment as numbers
	if {$categories != {} && $vBaseType == "Boolean"} {
	    set v(start) [string map {true 1 false 0} $v(start)]
	}
	# store variable infos
	foreach category $categories {
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
    if {$fmiVersion >= 2.0} {
        set catIdx 0
        foreach category {parameter state input output} {
            foreach index [set ::fmu::${fmuName}::${category}Indices] {
                set mapIndex($index) $catIdx
                incr catIdx
            }
        }
        foreach category {derivative} {
            foreach index [set ${category}Indices] {
                set mapIndex($index) $catIdx
                incr catIdx
            }
        }
        foreach element $fmiElements(ModelStructure) {
            set structureElements([lindex $element 0]) [lindex $element 2]
        }
	# collect given ModelStructure elements
	set categories {}
	foreach category {output derivative discreteState} {
            set elementName [string replace $category 0 0 \
                                 [string toupper [string index $category 0]]]s
	    if {[lsearch [array names structureElements] $elementName] >= 0} {
		lappend categories $category
	    }
	}
	# process ModelStructure elements
        foreach category $categories {
            set elementName [string replace $category 0 0 \
                                 [string toupper [string index $category 0]]]s
            foreach element $structureElements($elementName) {
                array set structureAttributes [lindex $element 1]
                if {[lsearch \
                     [array names structureAttributes] dependencies] < 0} {
                    # no dependencies given
                    continue
                }
                set var $mapIndex($structureAttributes(index))
                set deps {}
                foreach dependency $structureAttributes(dependencies) {
                    lappend deps $mapIndex($dependency)
                }
                set ::fmu::${fmuName}::${category}Dependencies($var) $deps
            }
        }
    }

    # export fmuAttributes
    array set ::fmu::${fmuName}::attributes [array get fmuAttributes]

    return $fmiVersion
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
		set item [string map {= " "} $item]
		set tagname [lindex $item 0] ;# start tag
		set rest [lrange $item 1 end]
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

## Unzip a zip archive.
# See: http://wiki.tcl.tk/3307
#      http://wiki.tcl.tk/15158
proc ::fmi::unzip fmuPath {
    set fd [open $fmuPath rb]
    # scan directory at end of zip file
    set off -22
    while 1 {
	seek $fd $off end
	binary scan [read $fd 4] i sig
	if {$sig == 0x06054b50} {
	    seek $fd $off end
	    break
	}
	incr off -1
    }
    binary scan [read $fd 22] issssiis sig disk cddisk nrecd nrec \
	dirsize diroff clen
    if {$clen > 0} {
	set comment [read $fd $clen]
    } else {
	set comment ""
    }
    if {$disk != 0} {
	error "multi-file zip not supported"
    }
    seek $fd $diroff
    # collect all directory entries
    for {set i 0} {$i < $nrec} {incr i} {
	binary scan [read $fd 46] issssssiiisssssii \
	    sig ver mver flag method time date crc csz usz n m k d ia ea \
	    off
	if {$sig != 0x02014b50} {
	    error "bad directory entry"
	}
	set name [read $fd $n]
	set extra [read $fd $m]
	if {$k == 0} {
	    set c ""
	} else {
	    set c [read $fd $k]
	}
	set directory($name) [dict create date $date time $time \
			      size $csz disksize $usz offset $off \
			      method $method extra $extra comment $c]
    }
    # extract all files
    set dirPath [::fmi::getDirPath $fmuPath]
    foreach name [array names directory] {
        if {[string match */ $name]} {
            # skip directory names
            continue
        }
        dict with directory($name) {}
	# extract file content
        seek $fd $offset
        binary scan [read $fd 30] isssssiiiss sig - - - - - - - - nlen xlen
        if {$sig != 0x04034b50} {
            error "not a file record"
        }
        seek $fd [expr {$nlen + $xlen}] current
        set data [read $fd $size]
        if {[string length $data] != $size} {
            error "read length mismatch: $size expected"
        }
        if {$method == 8} {
            set data [zlib inflate $data]
        } elseif {$method != 0} {
            error "unsupported method: $method"
        }
	# build directory
	set pathName ${dirPath}/${name}
	if {![file exists [file dirname $pathName]]} {
	    file mkdir [file dirname $pathName]
	}
	# write file
	set fp [open $pathName wb]
	puts -nonewline $fp $data
	close $fp
	# set file modification time from DOS timestamp
	# date: |Y|Y|Y|Y|Y|Y|Y|m| |m|m|m|d|d|d|d|d|
	# time: |H|H|H|H|H|M|M|M| |M|M|M|S|S|S|S|S|.|
	set ymdhms [list [expr (($date & 0xFE00) >> 9)+1980] \
			 [expr (($date & 0x01E0) >> 5)] \
			 [expr  ($date & 0x001F)] \
			 [expr (($time & 0xF800) >> 11)] \
			 [expr (($time & 0x07E0) >> 5)] \
			 [expr (($time & 0x001F) << 1)]]
	file mtime $pathName [clock scan $ymdhms \
				  -format "%Y %m %d %k %M %S" \
				  -timezone :UTC]
    }
    close $fd
}
