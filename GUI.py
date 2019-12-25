from tkinter import *
from tkinter.ttk import *
import calculation_DVM as dmv_cal


def button_command():
    # solution.main()
    airfoil = combo.get()
    re_num = float(text4.get())
    density = float(text5.get())
    viscosity = float(text6.get())
    free_aoa = float(text7.get())
    pl_amplitude = float(text9.get())
    pl_frequency = float(text10.get())
    time_step = float(text12.get())
    iteration = int(text13.get())
    new_vor_distance = float(text15.get())
    new_vor_angle = float(text16.get())

    result_file = float(chk_state1.get())
    force_file = float(chk_state2.get())
    excel_file = False  # float(chk_state3.get())

    grph_vortex_move = float(chk_state4.get())
    grph_clt = float(chk_state5.get())
    grph_cdt = float(chk_state6.get())
    grph_pld = float(chk_state7.get())

    dmv_cal.main(airfoil, re_num, density, viscosity, free_aoa, pl_amplitude, pl_frequency, time_step, iteration,
                 new_vor_distance, new_vor_angle, result_file, force_file, excel_file, grph_vortex_move, grph_clt,
                 grph_cdt, grph_pld)


def exit_command():
    exit()


label_width = 25
text_width = 15
window = Tk()
window.title("DVM Application")

# column 1
label0 = Label(window, text=' ')
label0.grid(row=0, column=0)

label1 = Label(window, text='Aerofoil selection', font=('Arial Black', 10), width=label_width, anchor=W)
label1.grid(row=1, column=1)
label2 = Label(window, text='Aerofoil name', width=label_width, anchor=W)
label2.grid(row=2, column=1)
combo = Combobox(window, width=text_width - 3)
combo['values'] = ('NACA0006', 'NACA0008', 'NACA0009', 'NACA0010', 'NACA0012', 'NACA0015', 'NACA0018', 'NACA0021',
                   'NACA0024', 'NACA1408', 'NACA1410', 'NACA1412', 'NACA2408', 'NACA2410', 'NACA2411', 'NACA2412',
                   'NACA2414', 'NACA2415', 'NACA2418', 'NACA2421', 'NACA2424', 'NACA4412', 'NACA4415', 'NACA4418',
                   'NACA4421', 'NACA4424', 'NACA6409', 'NACA6412')
combo.current(15)
combo.grid(row=2, column=2)
# text2 = Entry(window, width=text_width)
# text2.grid(row=2, column=2)

label3 = Label(window, text='Flow parameters', font=('Arial Black', 10), width=label_width, anchor=W)
label3.grid(row=3, column=1)

label4 = Label(window, text='Reynolds number', width=label_width, anchor=W)
label4.grid(row=4, column=1)
text4 = Entry(window, width=text_width)
text4.grid(row=4, column=2)
text4.insert(END, '1e6')

label5 = Label(window, text='Density', width=label_width, anchor=W)
label5.grid(row=5, column=1)
text5 = Entry(window, width=text_width)
text5.grid(row=5, column=2)
text5.insert(END, '1.225')

label6 = Label(window, text='Viscosity', width=label_width, anchor=W)
label6.grid(row=6, column=1)
text6 = Entry(window, width=text_width)
text6.grid(row=6, column=2)
text6.insert(END, '1.789e-5')

label7 = Label(window, text='Angle of attack - degrees', width=label_width, anchor=W)
label7.grid(row=7, column=1)
text7 = Entry(window, width=text_width)
text7.grid(row=7, column=2)
text7.insert(END, '0.0')

label8 = Label(window, text='Plunging parameters', font=('Arial Black', 10), width=label_width, anchor=W)
label8.grid(row=8, column=1)

label9 = Label(window, text='Plunging amplitude - m', width=label_width, anchor=W)
label9.grid(row=9, column=1)
text9 = Entry(window, width=text_width)
text9.grid(row=9, column=2)
text9.insert(END, '1.0')

label10 = Label(window, text='Plunging frequency - Hz', width=label_width, anchor=W)
label10.grid(row=10, column=1)
text10 = Entry(window, width=text_width)
text10.grid(row=10, column=2)
text10.insert(END, '5.0')

label11 = Label(window, text='Time parameters', font=('Arial Black', 10), width=label_width, anchor=W)
label11.grid(row=11, column=1)

label12 = Label(window, text='Time step - s', width=label_width, anchor=W)
label12.grid(row=12, column=1)
text12 = Entry(window, width=text_width)
text12.grid(row=12, column=2)
text12.insert(END, '0.005')

label13 = Label(window, text='Iteration', width=label_width, anchor=W)
label13.grid(row=13, column=1)
text13 = Entry(window, width=text_width)
text13.grid(row=13, column=2)
text13.insert(END, '1000')

label14 = Label(window, text='New vortex placement', font=('Arial Black', 10), width=label_width, anchor=W)
label14.grid(row=14, column=1)

label15 = Label(window, text='Distance - %chord', width=label_width, anchor=W)
label15.grid(row=15, column=1)
text15 = Entry(window, width=text_width)
text15.grid(row=15, column=2)
text15.insert(END, '0.001')

label16 = Label(window, text='Angle - degrees', width=label_width, anchor=W)
label16.grid(row=16, column=1)
text16 = Entry(window, width=text_width)
text16.grid(row=16, column=2)
text16.insert(END, '0.0')

# column 2
label01 = Label(window, text='         ')
label01.grid(row=1, column=3)

# column 3
label17 = Label(window, text='Data files', font=('Arial Black', 10), width=10, anchor=W)
label17.grid(row=1, column=4, sticky=W)

chk_state1 = BooleanVar()
chk_state1.set(True)
checkbtn1 = Checkbutton(window, text='Vortex position file', var=chk_state1)
checkbtn1.grid(row=2, column=4, sticky=W)
chk_state2 = BooleanVar()
chk_state2.set(True)
checkbtn2 = Checkbutton(window, text='Force calculation file', var=chk_state2)
checkbtn2.grid(row=3, column=4, sticky=W)
chk_state3 = BooleanVar()
chk_state3.set(False)
# checkbtn3 = Checkbutton(window, text='Excel file', var=chk_state3)
# checkbtn3.grid(row=4, column=4, sticky=W)

label18 = Label(window, text='Graph plot', font=('Arial Black', 10), width=10, anchor=W)
label18.grid(row=5, column=4, sticky=W)

chk_state4 = BooleanVar()
chk_state4.set(True)
checkbtn4 = Checkbutton(window, text='Vortex movement', var=chk_state4)
checkbtn4.grid(row=6, column=4, sticky=W)
chk_state5 = BooleanVar()
chk_state5.set(True)
checkbtn5 = Checkbutton(window, text='Cl vs time', var=chk_state5)
checkbtn5.grid(row=7, column=4, sticky=W)
chk_state6 = BooleanVar()
chk_state6.set(True)
checkbtn6 = Checkbutton(window, text='Cd vs time', var=chk_state6)
checkbtn6.grid(row=8, column=4, sticky=W)
chk_state7 = BooleanVar()
chk_state7.set(False)
checkbtn7 = Checkbutton(window, text='Plunging distance vs time', var=chk_state7)
checkbtn7.grid(row=9, column=4, sticky=W)

button1 = Button(window, text='Calculate', command=button_command)
button1.grid(row=17, column=4, sticky=E)
button2 = Button(window, text='Stop', command=exit_command)
button2.grid(row=18, column=4, sticky=E)

window.geometry('530x430')  # windwo size customization
window.resizable(0, 0)
window.mainloop()

"""
label1 = Label(window, text='label1', font=('Arial', 10))
label1.grid(column=0, row=0)

text1 = Entry(window, width=10)
text1.grid(column=1, row=0)
text1.focus()

button1 = Button(window, text='Import', command=button_command)
button1.grid(column=1, row=1)

chk_state = BooleanVar()
chk_state.set(True)
checkbtn1 = Checkbutton(window, text='Choose', var=chk_state)
checkbtn1.grid(column=0, row=2)
"""
