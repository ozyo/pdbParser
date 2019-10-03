touch areas.dat
for pdb in `ls -1 correct*pdb`
do
echo $pdb
$VMD -dispdev text -e align_direction.tcl -args areas.dat $pdb
done