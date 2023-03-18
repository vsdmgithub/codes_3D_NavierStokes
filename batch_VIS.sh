# !/bin/bash
dim_list0=("200" "275")
dim_list1=("080")
for i in ${dim_list1[@]}
do
	old="GHD"
	new=$i
	echo $new
	sed -i "s/$old/$new/g" NSE_code.f90
	sed -i "s/$old/$new/g" makefile
	make
	sed -i "s/$new/$old/g" NSE_code.f90
	sed -i "s/$new/$old/g" makefile
	cp run_VIS.sub run_VIS_$i.sub
	sed -i "s/$old/$new/g" run_VIS_$i.sub
done
# #=================================================
