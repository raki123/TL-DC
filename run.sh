for i in `ls one_pair/`
do 
    if ! timeout 10 ./build/read < one_pair/$i 2>/dev/null
    then
        echo "timeout"
    fi
done