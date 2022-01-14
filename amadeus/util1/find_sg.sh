#!/bin/sh

npop=50
for i in `seq -s " " -f %04g 1 $npop`
do 
head -1 POSCAR_$i  > ztmp
/usr/local/findSpaceGroup/build/findSpg POSCAR_$i |grep spg >> ztmp
cat ztmp |awk 'ORS=NR%2?FS:RS' | awk '{print $1, $2, $3, $5}' > zhead
tail -n +2 POSCAR_$i >> zhead
mv zhead POSCAR_$i 
rm ztmp
done
