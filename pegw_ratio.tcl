# Farhad Ramezanghorbani
# measure PEG/W bulk ratio

source libs.tcl

set mol [mol new "confout.gro" type gro waitfor all]
mol addfile "traj_comp.xtc" type xtc first 10000 last -1 step 1 waitfor all molid $mol
set nframe [molinfo top get numframes]

# Block Analysis
set diff [expr $nframe / 4]
for {set block 0} {$block < 4} {incr block} {
    set start [expr $block * $diff]
    set finish [expr $start + $diff - 1]
    puts "Block analysis for \[$start : $finish]"
    set fname7 "rate7block$block"
    
    puts "Calc. PEGWAA 0.7 nm cutoff"
    vicinRatio 7.0 $start $finish $fname7
} 

set start 0 
puts "Calc. PEGWAA 0.7 nm cutoff"
vicinRatio 7.0 $start $nframe rate7.out

puts "done"
exit
