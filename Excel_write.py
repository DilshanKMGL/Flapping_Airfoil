import xlsxwriter as xl
from _datetime import datetime as dt
import os

# get current date and time
ctime = dt.now().strftime("%Y-%m-%d %H.%M.%S")

try:
    p = os.listdir('Results')
    path_dir = 'Results/' + p[-1] # get the latest folder
    print(path_dir)
    # os.chdir(path_dir) # move to sub directory

    force_file = open(path_dir + '/force_file.txt')
    result_file = open(path_dir + '/result_file.txt')

    force_line = force_file.readlines()
    print(force_line[0][:-1])
    result_line = result_file.readlines()
    print(result_line[0][:-1])

    workbook = xl.Workbook(path_dir + '/Data.xlsx')
    worksheet = workbook.add_worksheet('NACA2412')
    workbook.close()
except FileNotFoundError:
    print('No such directory exists')
