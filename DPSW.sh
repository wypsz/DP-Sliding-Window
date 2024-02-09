#!/bin/bash

# 定义变量
rho_values=(0.1 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2)
all_sizes=(33 83 166 333 833 1666 3333 8333 16666 33333)
beta_values=(0.01 0.03 0.05 0.2)
bs=(5 4 3 2)
alpha_values=(0.75)
datasize=10000000

# Constants
w=1000000
sub_num=10
step=10000
n=3000000
q=0.8

cpp_program_path="./DPSWN/DP_Sliding_window"

# Nested loops
for rho in "${rho_values[@]}"
do
    for size in "${all_sizes[@]}"
    do
        # Calculate gamma as 1/(size*4)
        gamma=$(bc -l <<< "scale=7; 1/($size/4)")

        for j in "${!beta_values[@]}"
        do
            beta=${beta_values[j]}
            b=${bs[j]}
            for alpha in "${alpha_values[@]}"
            do
                for i in {1..10}
                do
                    # Generate unique file paths
                    input_path="./dataset/newfile.txt"
                    file_path_1="./res/DPSW_${i}_${sub_num}_${rho}_all${size}xb_1.txt"
                    file_path_2="./res/DPSW_${i}_${sub_num}_${rho}_all${size}xb_2.txt"
                    # Prepare and run the command
                    args=("$cpp_program_path" "$w" "$sub_num" "$step" "$rho" "$q" "$gamma" "$beta" "$alpha" "$n" "$i" "$file_path_1" "$file_path_2" "$input_path" "$datasize")
                    echo "Running command: ${args[*]}"
                    "${args[@]}"
                    echo "Command executed for size=$size, gamma=$gamma"
                done
            done
        done
    done
done
