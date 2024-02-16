FILES=(`ls instances/undirected_2024/`)
NR_THREADS=1
for j in $( seq 1 $NR_THREADS )
do
    for (( i = $j-1; i < ${#FILES[@]}; i += $NR_THREADS )); do
        if ! ( timeout 600 ./a.out -s pathfbs -t 12 -m 62000 < instances/undirected_2024/${FILES[$i]} 2>>results/latest | xargs echo "${FILES[$i]} $@"  )
        then
            echo "$i timeout"
        fi
    done &
done