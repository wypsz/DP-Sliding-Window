#!/bin/bash


w=1000000
n=3000000
cpp_program_path="./PCC/DP_Sliding_window"

all_budget_values=(0.1 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2)


lambda_values=(0.025 0.0375 0.05 0.0625 0.075 0.0875 0.1 0.1125 0.125)
true_values=(0.025 0.0375 0.05 0.0625 0.075 0.0875 0.1 0.1125 0.125)
# 获取数组长度
length=${#lambda_values[@]}

for (( i=0; i<$length; i++ ))
do
    lambda=${lambda_values[$i]}
    true_value=${true_values[$i]}

    for all_budget in "${all_budget_values[@]}"
    do
        args=("$cpp_program_path" "$w" "$lambda" "$all_budget" "$n" "$true_value" "./res/PCC_${lambda}_${all_budget}.txt" "./res/PCC_${lambda}_${all_budget}_time.txt" "./dataset/new_file.txt")
        echo "Running command with lambda=${lambda}, true_value=${true_value}, all_budget=${all_budget}: ${args[*]}"
        "${args[@]}"
        echo "Command executed for lambda=${lambda}, true_value=${true_value}, all_budget=${all_budget}"
    done
done