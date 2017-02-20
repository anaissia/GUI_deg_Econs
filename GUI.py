# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 08:16:06 2017

@author: anaissia
"""

import Tkinter
import ttk

from econs_gui import econs
"""
adder.py
~~~~~~

Creates a simple GUI for summing two numbers.
"""



class Adder(ttk.Frame):
    """The adders gui and functions."""
    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.root = parent
        self.init_gui()

    def on_quit(self):
        """Exits program."""
        quit()

#    def calculate(self):
#        """Calculates the sum of the two inputted numbers."""
#        num1 = int(self.num1_entry.get())
#        num2 = int(self.num2_entry.get())
#        num3 = num1 + num2
#        self.answer_label['text'] = num3


    def run(self):
        """Calculates the sum of the two inputted numbers."""
        path = str(self.path.get())
        capa = float(self.capa.get())
        unit = float(self.unit.get())
        start_passenger=float(self.start_passenger.get())
        eff_t=float(self.eff_t.get())
        eff_convert=float(self.eff_convert.get())
        reg_brak_eff=float(self.reg_brak_eff.get())
        gear_ratio=float(self.gear_ratio.get())
        radius=float(self.radius.get())
        M=float(self.M.get())
        A=float(self.A.get())
        rho_air=float(self.rho_air.get())
        C_drag=float(self.C_drag.get())
        aux_load=float(self.aux_load.get())
        T_max=float(self.T_max.get())
        P_max=float(self.P_max.get())
        
        econs(path,capa,unit,start_passenger,eff_t,eff_convert,reg_brak_eff,gear_ratio,radius,M,A, rho_air, C_drag, aux_load, T_max, P_max)
        root.quit()
        root.destroy()
    def CloseWindow(self):
        root.quit()  
        root.destroy()
    def init_gui(self):
        """Builds GUI."""
        self.root.title('A.F Ebus code')
        self.root.option_add('*tearOff', 'FALSE')
        
        self.grid(column=0, row=0, sticky='nsew')
        self.root.configure(background='black')
        """ Menu Bar """
        self.menubar = Tkinter.Menu(self.root)

        self.menu_file = Tkinter.Menu(self.menubar)
        self.menu_file.add_command(label='Exit', command=self.on_quit)

        self.menu_edit = Tkinter.Menu(self.menubar)

        self.menubar.add_cascade(menu=self.menu_file, label='File')
        self.menubar.add_cascade(menu=self.menu_edit, label='Edit')

        self.root.config(menu=self.menubar)
        """ All num entries""" 
        self.path = ttk.Entry(self, width=35)
        self.path.grid(column=1, row = 2)
        self.path.insert(0, "Brampton_23_East_Heavy Duty")
        
        self.capa = ttk.Entry(self, width=35)
        self.capa.grid(column=6, row = 2, columnspan=4)
        self.capa.insert(0, "200")
        
        self.unit = ttk.Entry(self, width=35)
        self.unit.grid(column=1, row = 3,columnspan=4)
        self.unit.insert(0, "1")

        self.start_passenger = ttk.Entry(self, width=35)
        self.start_passenger.grid(column=6, row=3,columnspan=4)
        self.start_passenger.insert(0, "61")
        
        self.eff_t = ttk.Entry(self, width=35)
        self.eff_t.grid(column= 1, row=4)
        self.eff_t.insert(0, "0.95")

        self.eff_convert = ttk.Entry(self, width=35)
        self.eff_convert.grid(column=6, row=4,columnspan=4, sticky='ew')
        self.eff_convert.insert(0, "0.97")

        
        self.reg_brak_eff = ttk.Entry(self, width=35)
        self.reg_brak_eff.grid(column= 1, row=5)
        self.reg_brak_eff.insert(0, "0.4")

        self.gear_ratio = ttk.Entry(self, width=35)
        self.gear_ratio.grid(column=6, row=5,columnspan=4, sticky='ew')
        self.gear_ratio.insert(0, "5.67")


        self.radius = ttk.Entry(self, width=35)
        self.radius.grid(column=1, row=6,columnspan=4, sticky='ew')
        self.radius.insert(0, "0.571")
        
        self.M = ttk.Entry(self, width=35)
        self.M.grid(column=6, row=6,columnspan=4, sticky='ew')
        self.M.insert(0, "14864")        

        self.A = ttk.Entry(self, width=35)
        self.A.grid(column=1, row=7,columnspan=4, sticky='ew')
        self.A.insert(0, "7.44")
        
        self.rho_air = ttk.Entry(self, width=35)
        self.rho_air.grid(column=6, row=7,columnspan=4, sticky='ew')
        self.rho_air.insert(0, "1.225")        

        self.C_drag = ttk.Entry(self, width=35)
        self.C_drag.grid(column=1, row=8,columnspan=4, sticky='ew')
        self.C_drag.insert(0, "0.65")
        
        self.aux_load = ttk.Entry(self, width=35)
        self.aux_load.grid(column=6, row=8,columnspan=4, sticky='ew')
        self.aux_load.insert(0, "10000")  

        self.T_max = ttk.Entry(self, width=35)
        self.T_max.grid(column=1, row=9,columnspan=4, sticky='ew')
        self.T_max.insert(0, "2500")
        
        self.P_max = ttk.Entry(self, width=35)
        self.P_max.grid(column=6, row=9,columnspan=4, sticky='ew')
        self.P_max.insert(0, "153000")  
                
        self.run_button = ttk.Button(self, text='Run !',
                command=self.run)
        self.run_button.grid(column=3, row=10, columnspan=4)

        
                
     
        self.quit_text = ttk.Button(self, text='Quit',
                command=self.CloseWindow)
                
        self.quit_text.grid(column=3, row=12, columnspan=4)        

     
#        self.answer_frame = ttk.LabelFrame(self, text='Answer',
#                height=100)
#        self.answer_frame.grid(column=0, row=4, columnspan=4, sticky='nesw')
#
#        self.answer_label = ttk.Label(self.answer_frame, text='')
#        self.answer_label.grid(column=0, row=0)

        # Labels that remain constant throughout execution.
        ttk.Label(self, text='Electric bus energy consumption model').grid(column=0, row=0,
                columnspan=10)
        ttk.Label(self, text='Path').grid(column=0, row=2,
                sticky='w')
        ttk.Label(self, text='Battery Capacity (Ah)').grid(column=5, row=2,
                sticky='w')
        ttk.Label(self, text='Unit (1 for km/hour or m/s defined in import_data or 2 for MPH)').grid(column=0, row=3,
                sticky='w')
        ttk.Label(self, text='Number of passengers').grid(column=5, row=3,
                sticky='w')                
        ttk.Label(self, text='Transmission efficiency').grid(column=0, row=4,
                sticky='w')
        ttk.Label(self, text='Converter efficiency').grid(column=5, row=4,
                sticky='w')              
        ttk.Label(self, text='Regenerative braking efficiency').grid(column=0, row=5,
                sticky='w')
        ttk.Label(self, text='Gear ratio').grid(column=5, row=5,
                sticky='w')      
        ttk.Label(self, text='Wheel radius (m)').grid(column=0, row=6,
                sticky='w')
        ttk.Label(self, text='Vehicle weight (unloaded)').grid(column=5, row=6,
                sticky='w')                 

        ttk.Label(self, text='Frontal Area (m2)').grid(column=0, row=7,
                sticky='w')
        ttk.Label(self, text='Air density (kg/m3)').grid(column=5, row=7,
                sticky='w')   

        ttk.Label(self, text='Drag coefficient').grid(column=0, row=8,
                sticky='w')
        ttk.Label(self, text='Auxilliary load (W) ').grid(column=5, row=8,
                sticky='w')                   
                
        ttk.Label(self, text='Maximum torque (Nm)').grid(column=0, row=9,
                sticky='w')
        ttk.Label(self, text='Maximum power (W) ').grid(column=5, row=9,
                sticky='w')                  
                
        ttk.Separator(self, orient='horizontal').grid(column=0,
                row=1, columnspan=10, sticky='ew')

        for child in self.winfo_children():
            child.grid_configure(padx=10, pady=10)

if __name__ == '__main__':
    root = Tkinter.Tk()
    Adder(root)
    root.mainloop()
