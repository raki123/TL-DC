NR_THREADS=1
STRATEGY=treefbs

FILES=(`ls instances/one_pair/`)
for j in $( seq 1 $NR_THREADS )
do
    for (( i = $j-1; i < ${#FILES[@]}; i += $NR_THREADS )); do
        if ! ( timeout 600 ./a.out -s $STRATEGY -t 12 -m 62000 < instances/one_pair/${FILES[$i]} 2>>results/latest | xargs echo "${FILES[$i]} $@"  )
        then
            echo "$i timeout"
        fi
    done 
done

FILES=(`ls instances/undirected_2024/`)
for j in $( seq 1 $NR_THREADS )
do
    for (( i = $j-1; i < ${#FILES[@]}; i += $NR_THREADS )); do
        if ! ( timeout 600 ./a.out -s $STRATEGY -t 12 -m 62000 < instances/undirected_2024/${FILES[$i]} 2>>results/latest | xargs echo "${FILES[$i]} $@"  )
        then
            echo "$i timeout"
        fi
    done 
done

FILES=(`ls instances/directed_2024/`)
for j in $( seq 1 $NR_THREADS )
do
    for (( i = $j-1; i < ${#FILES[@]}; i += $NR_THREADS )); do
        if ! ( timeout 600 ./a.out -d -s $STRATEGY -t 12 -m 62000 < instances/directed_2024/${FILES[$i]} 2>>results/latest | xargs echo "${FILES[$i]} $@"  )
        then
            echo "$i timeout"
        fi
    done 
done
