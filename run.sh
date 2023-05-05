FILES=(`ls one_pair/`)
NR_THREADS=4
for j in $( seq 1 $NR_THREADS )
do
    for (( i = $j-1; i < ${#FILES[@]}; i += $NR_THREADS )); do
        if ! ( timeout 10 ./build/read < one_pair/${FILES[$i]} 2>/dev/null | xargs echo "$i $@"  )
        then
            echo "$i timeout"
        fi
    done &
done