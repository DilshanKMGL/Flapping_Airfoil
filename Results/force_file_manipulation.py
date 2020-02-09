# time_step_pl_pl = [0.02, 0.01, 0.005, 0.004, 0.002, 0.001]
# iteration_pl_pl = [250, 500, 1000, 1250, 2500, 5000]

name = '0.0 0.001 5000 force_file'
file1 = open(name + '.txt', 'r')
file2 = open(name + ' updated.txt', 'w')
lines = file1.readlines()
iteration = int(lines[21])
start_line = 28

for index in range(start_line):
    file2.write(lines[index])

step_val = 20
for index in range(int(250)):
    file2.write(lines[start_line + index*step_val])
