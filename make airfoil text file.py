file_list = ['0006', '0008', '0009', '0010', '0012', '0015', '0018', '0021', '0024', '1408', '1410', '1412', '2408',
             '2410', '2411', '2412', '2414', '2415', '2418', '2421', '2424', '4412', '4415', '4418', '4421', '4424',
             '6409', '6412']

for file in file_list:
    file_name = 'NACA' + file + '.txt'
    file1 = open(file_name, 'r')
    data = file1.readlines()
    file1.close()

    write_data = []
    for line in data:
        row = line[2:-1].split(' ')
        if len(row) == 2:
            write_data.append(str(row[0]) + '+' + str(row[1]) + 'i')
        else:
            write_data.append(str(row[0]) + '+' + str(row[2]) + 'i')

    file_name = 'NACA' + str(file) + '_updated.txt'
    file1 = open(file_name, 'w')
    file1.close()
    file1 = open(file_name, 'a+')

    file1.write('[')
    for i in range(len(write_data)):
        if len(write_data) == (i + 1):
            file1.write(write_data[i])
        else:
            file1.write(write_data[i] + ';\n')
    file1.write(']')
