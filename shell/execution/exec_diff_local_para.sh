#!/bin/bash

#mesh
#一方向分割数
e=30
#解析領域の大きさ
ep=5

#podモード数
num_modes=(10 20)
#POD計算領域数
num_1stdd=(8 16 32)
#並列計算領域数 (=並列数)
num_parallel=(2 4 8)
#基底本数可変の閾値 1.0E-{pa}
pa=0
#solver type
st=3

for nm in "${num_modes[@]}"
do
	for nd in "${num_1stdd[@]}"
	do
    	for np in "${num_parallel[@]}"
    	do
	
        . shell/diff/meshgen.sh $e $ep $nm $nd $np $pa
	    . shell/diff/merge_graph.sh $e $ep $nm $nd $np $pa
	    . shell/diff/execution.sh $e $ep $nm $nd $np $pa $st

        done
	done
done
