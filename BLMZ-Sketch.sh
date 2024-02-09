#!/bin/bash

w=1000000
n=3000000
cpp_program_path="./ICLR/DP_Sliding_window"

# 4层2倍
lambda_values=(0.256 0.105 0.052 0.0256 0.0105 0.0052 0.00256 0.00105 0.00052 0.000256)
all_budget_values=(0.1 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2)

length=${#lambda_values[@]}
level=3
for (( i=0; i<$length; i++ ))
do
    lambda=${lambda_values[$i]}

    for all_budget in "${all_budget_values[@]}"
    do
        args=("$cpp_program_path" "$w" "$level" "$lambda" "$all_budget" "$n" "./res/ICLR_${level}_${lambda}_${all_budget}.txt" "./res/ICLR_${level}_${lambda}_${all_budget}_time.txt" "./dataset/new_file.txt")
        echo "Running command with lambda=${lambda}, all_budget=${all_budget}: ${args[*]}"
        "${args[@]}"
        echo "Command executed for lambda=${lambda}, all_budget=${all_budget}"
    done
done
