# Farhad Ramezanghorbani
# Measure the PEG/W bulk ratio 
# at the vicinity (R = 0.7 nm) 
# of each amino acid

proc vicinRatio {R start end fout} {
    set outfile [open $fout w]
    set sel [atomselect top "all not (resname PEG4 or resname W or resname ION)"]
    set reslist [lsort -unique [$sel get resname]]
    set W [atomselect top "resname W"]
    set P [atomselect top "resname PEG4"]
    set Wtot [$W num]
    set Ptot [$P num]
    set bulk [expr {double($Ptot)/(4*double($Wtot))}]
    
    set nframe [molinfo top get numframes]
    set nframe [expr $end - $start + 1]
    
    foreach res $reslist {
        set WL [atomselect top "resname W and pbwithin $R of resname $res"]
        set PL [atomselect top "resname PEG4 and pbwithin $R of resname $res"]
        set Wcount 0
        set Pcount 0
        set wval 0
        set pval 0    
        for {set frame $start} {$frame  <= $end} {incr frame} {
            $WL frame $frame
            $WL update 
            $PL frame $frame
            $PL update
            set wval [$WL num]
            set pval [$PL num]
            set Wcount [expr {$Wcount + 4*($wval)}]
            set Pcount [expr {$Pcount + $pval}] 
        }
                puts "$res : [expr {double($Pcount)/double($Wcount)/double($bulk)}]"
        puts $outfile "$res : [expr {double($Pcount)/double($Wcount)/double($bulk)}]"
    }
    puts "Bulk ratio : $bulk"
    puts $outfile "Bulk ratio : $bulk"
    close $outfile
}

# return dictionary of data: ratio(beadtype)
proc vicinRatioBeads {R start end List} {
    global rateArray
    set sel [atomselect top "all not (resname PEG or resname W or resname ION)"]
    set namelist [lsort -unique [$sel get name]]     
    set W [atomselect top "resname W"]
    set P [atomselect top "resname PEG"]
    set Wtot [$W num]
    set Ptot [$P num]
    set bulk [expr {double($Ptot)/(4*double($Wtot))}] 

    set nframe [molinfo top get numframes]
    set nframe [expr $end - $start + 1]
    # changed namelist to List and name to index after here
    foreach index $List {
        set WL [atomselect top "resname W and pbwithin $R of index $index"]
        set PL [atomselect top "resname PEG and pbwithin $R of index $index"]
        set Wcount 0
        set Pcount 0
        set wval 0
        set pval 0
        for {set frame $start} {$frame  <= $end} {incr frame} {
            $WL frame $frame
            $WL update
            $PL frame $frame
            $PL update
            set wval [$WL num]
            set pval [$PL num]
            set Wcount [expr {$Wcount + 4*($wval)}]
            set Pcount [expr {$Pcount + $pval}]
        }
        puts "done with $index"
        #puts "$name : [expr {double($Pcount)/double($Wcount)/double($bulk)}]"
        set rateArray($index) [expr {double($Pcount)/double($Wcount)/double($bulk)}]
    }
    #puts "Bulk ratio : $bulk"
    set rateArray(bulk) $bulk
}

