import os
from tkinter import *
from tkinter import filedialog
import tkinter.messagebox
from pygame import mixer
import time
from mutagen.mp3 import MP3
import threading

root = Tk()
# root.geometry('300x400')
root.title('Melody')
root.iconbitmap(r'images/melody.ico')


def about_us():
    tkinter.messagebox.showinfo('About Melody', 'This is a music player using Python tkinter')


def browse_file():
    global filename
    filename = filedialog.askopenfilename()


# create menu bar
menubar = Menu(root)
root.config(menu=menubar)
# create submenu
submenu = Menu(menubar, tearoff=0)
menubar.add_cascade(label='File', menu=submenu)
submenu.add_command(label='Open', command=browse_file)
submenu.add_command(label='Exit', command=root.destroy)

submenu = Menu(menubar, tearoff=0)
menubar.add_cascade(label='Help', menu=submenu)
submenu.add_command(label='About Us', command=about_us)

leftFrame = Frame(root)
leftFrame.pack(side=LEFT, padx=30)

lb1 = Listbox(leftFrame)
lb1.insert(0, '1')
lb1.insert(1, '2')
lb1.pack()

bt1 = Button(leftFrame, text = '+ Add')
bt1.pack(side=LEFT)
bt2 = Button(leftFrame, text = '+ Add')
bt2.pack(side=LEFT)

rightFrame = Frame(root)
rightFrame.pack()

topFrame = Frame(rightFrame)
topFrame.pack()

mixer.init()  # initializing the mixer
paused = False

lengthLabel = Label(topFrame, text='Total time - --:--')
lengthLabel.pack(pady=5)
currentTimeLabel = Label(topFrame, text='Current time - --:--', relief=GROOVE)
currentTimeLabel.pack(pady=0)
