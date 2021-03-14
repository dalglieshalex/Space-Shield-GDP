# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 14:29:29 2021

@author: dalgl
"""
from tkinter import *
#from PIL import ImageTK, Image

root = Tk()

Windframe = Frame(root, padx = 10, pady = 10).pack()

reqframe = LabelFrame(Windframe, text = "Select Shield Requirements", padx = 5, pady = 5)
reqframe.pack()

didef = Label(reqframe, text="Max particle diameter").grid(row=0, column=0)
partmat = Label(reqframe, text="Particle Material").grid(row=1, column=0)
Succrt = Label(reqframe, text="Shield success rate").grid(row=2, column=0)



root.mainloop()