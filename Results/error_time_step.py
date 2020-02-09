name1 = '0.0 0.002 2500 force_file'
name2 = '0.0 0.001 5000 force_file'
start_line = 28
file1 = open(name1+'.txt', 'r')
line1 = file1.readlines()
iteration1 = int(line1[21])
time_step1 = float(line1[17])

file2 = open(name2+'.txt', 'r')
line2 = file2.readlines()
iteration2 = int(line2[21])
time_step2 = float(line2[17])

file3 = open(str(time_step1) + ' ' + str(time_step2) + ' updated.txt', 'w')
for index in range(start_line):
    file3.write(line1[index])

for index in range(start_line):
    file3.write(line2[index])

step_val = 5
for index in range(int(500)):
    file3.write(line1[start_line + index*step_val])

step_val = 10
# for index in range(int(iteration2/step_val)):
for index in range(int(500)):
    file3.write(line2[start_line + index*step_val])
