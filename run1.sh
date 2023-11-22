#!/bin/bash

# Define the arrays of values
rho_values=(0.599)
gamma_values=(0.002)
beta_values=(0.00005)
alpha_values=(0.5 0.6 0.8 0.99)

# Constants
w=1000000
sub_num=10
step=2000
all_budget=1
n=10000000
query=20000
cpp_program_path="./build/DP_Sliding_window"

# Counter for unique file paths
run_counter=52

# Nested loops
for rho in "${rho_values[@]}"
do
    for gamma in "${gamma_values[@]}"
    do
        for beta in "${beta_values[@]}"
        do
            for alpha in "${alpha_values[@]}"
            do
                # Generate unique file paths
                file_path_1="./results/output_${run_counter}_1.txt"
                file_path_2="./results/output_${run_counter}_2.txt"
                file_path_3="./results/output_${run_counter}_3.txt"
                run_counter=$((run_counter + 1))

                # Prepare and run the command
                args=("$cpp_program_path" "$w" "$sub_num" "$step" "$all_budget" "$rho" "$gamma" "$beta" "$alpha" "$(bc <<< "1-($rho+0.001)")" "$n" "$query" "$file_path_1" "$file_path_2" "$file_path_3")
                echo "Running command: ${args[*]}"
                "${args[@]}"
                echo "Command executed"
            done
        done
    done
done