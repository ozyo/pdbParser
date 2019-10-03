package require Orient
namespace import Orient::orient
puts [lindex $argv 1]
mol load pdb [lindex $argv 1]

set sel [atomselect top "all"]
set I [draw principalaxes $sel]
puts "align axis is [lindex $I 2]"
set im [measure inertia [atomselect top "all"]]
puts "measured inertia $im"
set A [orient $sel [lindex $I 2] {0 0 1}]
$sel move $A
set sel [atomselect top all]
set chains [lsort -ascii -unique [$sel get chain]]

set coma [measure center [atomselect top "chain [lindex $chains 0]"] weight mass]
set comb [measure center [atomselect top "chain [lindex $chains 1]"] weight mass]
set comc [measure center [atomselect top "chain [lindex $chains 2]"] weight mass]
set comd [measure center [atomselect top "chain [lindex $chains 3]"] weight mass]
set come [measure center [atomselect top "chain [lindex $chains 4]"] weight mass]

set e1 [expr ( [lindex $comb 0] - [lindex $coma 0] ) * ( [lindex $comb 1] + [lindex $coma 1] ) ]
set e2 [expr ( [lindex $comc 0] - [lindex $comb 0] ) * ( [lindex $comc 1] + [lindex $comb 1] ) ]
set e3 [expr ( [lindex $comd 0] - [lindex $comc 0] ) * ( [lindex $comd 1] + [lindex $comc 1] ) ]
set e4 [expr ( [lindex $come 0] - [lindex $comd 0] ) * ( [lindex $come 1] + [lindex $comd 1] ) ]
set e5 [expr ( [lindex $coma 0] - [lindex $come 0] ) * ( [lindex $coma 1] + [lindex $come 1] ) ]

set summed [ expr $e1+$e2+$e3+$e4+$e5]

set com [measure center [atomselect top all] weight mass]
set x1 [measure center [atomselect top "chain A and index 1"] weight mass]

if {[ lindex $x1 2] > [lindex $com 2]} {
    set summed [expr $summed * -1]
    set sign "change"
} else {
    set sign "same"
}

set outf [open "[lindex $argv 0]" a]

puts $outf  "[lindex $argv 1] $summed $sign" 

close $outf

quit