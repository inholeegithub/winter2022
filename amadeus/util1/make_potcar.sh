variable=`awk -v line=6 'NR==line' POSCAR`
echo $variable
number=0
for element in $variable ;do
echo $element
number=`expr $number + 1`
if [ $number -eq 1 ];then
cat /TGM/Apps/VASP/POTCAR/1.POTPAW.LDA.54.RECOMMEND/$element/POTCAR  > ./POTCAR
#cat /TGM/Apps/VASP/POTCAR/2.POTPAW.PBE.54.RECOMMEND/$element/POTCAR  > ./POTCAR
else
cat /TGM/Apps/VASP/POTCAR/1.POTPAW.LDA.54.RECOMMEND/$element/POTCAR  >> ./POTCAR
#cat /TGM/Apps/VASP/POTCAR/1.POTPAW.PBE.54.RECOMMEND/$element/POTCAR  >> ./POTCAR
fi
done
