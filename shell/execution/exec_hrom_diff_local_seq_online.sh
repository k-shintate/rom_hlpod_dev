#!/bin/bash

#mesh
#一方向分割数
e=20
#解析領域の大きさ
ep=5

#podモード数
num_modes=(20)
#POD計算領域数
num_1stdd=(6)
#並列計算領域数 (=並列数)
num_parallel=(1)
#基底本数可変の閾値 1.0E-{pa}
pa=0
#solver type
st=2

for nm in "${num_modes[@]}"
do
	for nd in "${num_1stdd[@]}"
	do
    	for np in "${num_parallel[@]}"
    	do
	
	    . shell/diff_hrom/execution_online_seq.sh $e $ep $nm $nd $np $pa $st

        done
	done
done

