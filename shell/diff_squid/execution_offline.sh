#!/bin/bash

#podモード数
nm=$3
#POD計算領域数
nd=$4
#並列計算領域数 (=並列数)
np=$5
#基底本数可変の閾値 1.0E-{pa}
pa=$6
#solver type
st=$7

# 実行ディレクトリ
directory="result_diff/${nm}-${np}-${nd}"

#. shell/install.sh

cd solvers/diff

make -f Makefile_ROM_squid clean
make -f Makefile_ROM_squid

fname="execution.sh"
touch $fname
echo "#!/bin/bash" >> $fname
echo "#------- qsub option -----------" >>$fname
echo "#PBS -q SQUID" >> $fname                          # debug == DBG, regular == SQUID
echo "#PBS --group=jh240017" >> $fname
echo "#PBS -b ${N_node}" >>$fname                       # node DBG == 2, SQUID == max 512
echo "#PBS -l elapstim_req=1:00:00" >>$fname            # debug == max 0:10:00, regular == max 24:00:00
echo "#PBS -l cpunum_job=${N_cpu}" >> $fname            # CPU/node == max 72
echo "#PBS -l memsz_job=248GB" >> $fname
echo "#PBS -T intmpi" >> $fname
echo "#------- Program execution -----------" >> $fname
echo "module load BaseCPU/2023" >> $fname
echo "cd \$PBS_O_WORKDIR" >> $fname
echo "" >> $fname

echo "cp -r hlpod_diff_offline ./../../$directory" >> $fname
echo "cp -r hlpod_diff_online ./../../$directory" >> $fname

echo "cd ./../../$directory" >> $fname

echo "mkdir -p {pod_modes_vtk,pod_modes,fem_solver_prm,pod_solver_prm,calctime}" >> $fname
echo "for ((i=0; i<${nd}; i++))" >> "$fname"
echo "do" >> "$fname"
echo '    mkdir -p "pod_modes/subdomain${i}"' >> "$fname"
echo "done" >> "$fname"


echo "mpirun -np $np  ./hlpod_diff_offline ./ -nd $nd -nm $nm -pa $pa -st $st" >> $fname
echo "mpirun -np $np  ./hlpod_diff_online ./ -nd $nd -nm $nm -pa $pa -st $st" >> $fname
echo "cd ../.." >> $fname