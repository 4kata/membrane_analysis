#!/bin/bash

cd programs

program_list=`ls`
echo ${program_list}
filtered_list=$(echo "$program_list" | sed '/bin/d; /mod/d')
echo $filtered_list
IFS=$'\n' read -r -d '' -a filtered_array <<<"$filtered_list"
rm -f ./bin/${filtered_array}

for i_dir in "${filtered_array[@]}"
do 
  cd ./${i_dir}
  rm -f ${i_dir}
  make
  cd ../
done
cd ./bin


