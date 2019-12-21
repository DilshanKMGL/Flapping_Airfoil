import xlsxwriter as xl
from _datetime import datetime as dt
import os

ctime = dt.now().strftime("%Y-%m-%d %H.%M.%S")
if not os.path.exists('Results'):
    os.mkdir('Results')
os.mkdir('Results/Data '+ctime)
airfoil_name = 2412
workbook = xl.Workbook('Data.xlsx')
worksheet = workbook.add_worksheet('NACA2412')
try:
    force_file = open('force_file_NACA'+str(airfoil_name)+'.txt')
    result_file = open('result_file_NACA'+str(airfoil_name)+'.txt')
except FileNotFoundError:
    print('no such file')
workbook.close()

p = os.listdir("Results") # get all the folders name in a list
print(p)
