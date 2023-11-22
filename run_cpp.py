import subprocess


# 参数列表
rho_values = [0.599,0.699,0.799,0.899]
gamma_values =[0.0001,0.00005,3e-5,2.5e-5,2e-5,1.6e-5,1.4e-5,1.25e-5,1.11e-5,1e-5]
beta_values = [0.1,0.07,0.02,0.01,0.003,0.001,0.0005,0.0001,0.00005,0.00002,0.00001]
alpha_values = [0.2,0.22,0.23,0.24,0.3,0.4,0.45,0.5,0.6,0.8,0.99]


w = 1000000
sub_num = 10
step = 2000
all_budget = 1
n = 10000000
query = 20000

cpp_program_path = "./cmake-build-debug/DP_Sliding_window.exe"
run_counter = 0

for rho in rho_values:
    for gamma in gamma_values:
        for beta in beta_values:
            for alpha in alpha_values:
                # 生成独特的文件路径
                file_path_1 = f'./results/output_{run_counter}_1.txt'
                file_path_2 = f'./results/output_{run_counter}_2.txt'
                file_path_3 = f'./results/output_{run_counter}_3.txt'
                run_counter += 1
                args = [cpp_program_path,str(w),str(sub_num),str(step),str(all_budget),str(rho),str(gamma),str(beta),str(alpha),str((1-(rho+0.001))),str(n),str(query),file_path_1,file_path_2,file_path_3]
                print("Running command:",args)
                subprocess.run(args,check=True)

                print("Command executed")

