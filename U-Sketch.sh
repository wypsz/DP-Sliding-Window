#!/bin/bash


w=1000000
sub_num=100
zeta_values=(0.01 0.03 0.05 0.2)
bs=(5 4 3 2)
beta=0.9
n=3000000
step=1000
cpp_program_path="./UP/DP_Sliding_window"
result_dir="./res"
epsilons=(0.1 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2)
all_size=(5 12 25 50 125 250 500 1250 2500 5000)
# 遍历 epsilon
for epsilon in "${epsilons[@]}"
do
    # 使用索引遍历 all_size
    for i in "${!all_size[@]}"
    do
        size="${all_size[i]}"

        for j in "${!zeta_values[@]}"
        do
          zeta="${zeta_values[j]}"
          b="${bs[j]}"
          gamma=$(echo "scale=7; e(1)/($size / $b)" | bc -l)

          # 设置文件名
          result_path="$result_dir/UP_${sub_num_}${epsilon}_${size}x${b}.txt"
          time_path="$result_dir/UP_${sub_num}_${epsilon}_${size}x${b}_time.txt"
          input_path="./dataset/new_file.txt"

          # 构建参数数组
          args=("$cpp_program_path" "$w" "$sub_num" "$step" "$epsilon" "$gamma" "$beta" "$zeta" "$n" "$result_path" "$time_path" "$input_path")

          # 打印要运行的命令
          echo "Running command: ${args[*]}"

          # 执行命令
          "${args[@]}"

          # 表明命令执行完毕
          echo "Command executed for epsilon=$epsilon, size=$size"
          done
    done
done
