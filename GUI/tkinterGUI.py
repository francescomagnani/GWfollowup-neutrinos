import tkinter as tk
from tkinter import *
from tkinter import filedialog
from PIL import Image, ImageTk

import selection.optimization


class tkGUI:
    def __init__(self,conn,itt):
        self.root = tk.Tk()
        #self.root.attributes('-fullscreen', True)
        self.root.geometry("600x600")
        self.root.title("F.Magnani GW follow-up")
        Grid.rowconfigure(self.root, 0, weight=1)
        Grid.rowconfigure(self.root, 1, weight=1)
        Grid.rowconfigure(self.root, 2, weight=1)
        Grid.rowconfigure(self.root, 3, weight=1)
        Grid.rowconfigure(self.root, 4, weight=1)
        Grid.rowconfigure(self.root, 5, weight=1)
        Grid.rowconfigure(self.root, 6, weight=1)
        Grid.rowconfigure(self.root, 7, weight=1)
        Grid.rowconfigure(self.root, 8, weight=1)
        Grid.rowconfigure(self.root, 9, weight=1)
        Grid.columnconfigure(self.root, 1, weight=1)
        Grid.columnconfigure(self.root, 2, weight=1)
        Grid.columnconfigure(self.root, 3, weight=1)
        Grid.columnconfigure(self.root, 4, weight=1)
        self.label = tk.Label(self.root, text="GW follow-up with your neutrino telescope", font=("Arial",20))
        self.label.grid(row=0,column=0, sticky="nsew")
        self.label.pack()
        # select GW
        self.browser = tk.Frame(self.root)
        self.browser.columnconfigure(0,weight=1)
        self.browser.columnconfigure(1,weight=1)
        self.label = tk.Label(self.browser, text = "Select a GW file", font=('Arial',10))
        self.label.grid(row=1, column=0, sticky="nsew")
        self.brw = tk.Button(self.browser, text="Browse file", command=lambda : self.set_display(conn,itt))
        self.brw.grid(row=1, column=1, sticky="nsew")
        self.browser.pack(side = tk.TOP, anchor = NW)
        # optimization variables
        self.likelihood = tk.StringVar()
        self.tracklength = tk.StringVar()
        self.hits = tk.StringVar()
        self.bdtscore = tk.StringVar()

        print("\n")
        print("Please, press the browse button and upload a '.fits' file that I can analyse (you may find a few files in the 'alerts' folder)")

        # variables
        # display+radio buttons
        self.dispradio = tk.Frame(self.root)
        self.dispradio.columnconfigure(0,weight=1)
        self.dispradio.columnconfigure(1,weight=1)
        self.label_var = tk.Label(self.dispradio, text="Optimization variables:", font=("Arial",10))
        self.label_var.grid(row=2,column=1, sticky="nsew")
        self.ch1 = tk.Checkbutton(self.dispradio, text="likelihood", variable=self.likelihood, onvalue="likelihood", offvalue="")
        self.ch1.grid(row=3,column=1, sticky="nsew")
        self.ch2 = tk.Checkbutton(self.dispradio, text="tracklength", variable=self.tracklength, onvalue="fitinf10", offvalue="")
        self.ch2.grid(row=3,column=2, sticky="nsew")
        self.ch3 = tk.Checkbutton(self.dispradio, text="n. of hits", variable=self.hits, onvalue="fitinf3", offvalue="")
        self.ch3.grid(row=3,column=3, sticky="nsew")
        self.ch4 = tk.Checkbutton(self.dispradio, text="bdt score", variable=self.bdtscore, onvalue="bdt_score", offvalue="")
        self.ch4.grid(row=3,column=4, sticky="nsew")
        self.label_mods = tk.Label(self.dispradio, text="Operation modalities:", font=("Arial",10))
        self.label_mods.grid(row=4,column=1, sticky="nsew")
        self.mod1 = tk.Label(self.dispradio, text="Phase space", font=("Arial",10))
        self.mod1.grid(row=5,column=1, sticky="nsew")
        self.txt1 = tk.Entry(self.dispradio)
        self.txt1.grid(row=5,column=2, sticky="nsew")
        self.mod2 = tk.Label(self.dispradio, text="Orders", font=("Arial",10))
        self.mod2.grid(row=5,column=3, sticky="nsew")
        self.txt2 = tk.Entry(self.dispradio)
        self.txt2.grid(row=5,column=4, sticky="nsew")
        self.dispradio.pack(after = self.browser)

        # start / stop
        self.buttonframe = tk.Frame(self.root)
        self.buttonframe.columnconfigure(0,weight=1)
        self.buttonframe.columnconfigure(1,weight=1)
        self.button1 = tk.Button(self.buttonframe, text="Start", font=('Arial',18), command=lambda : self.set_display2())#, state = DISABLED)
        self.button1.grid(row=6, column=0, sticky="nsew")
        self.button2 = tk.Button(self.buttonframe, text="Stop", font=('Arial',18), command=lambda : self.root.destroy())
        self.button2.grid(row=6, column=1, sticky="nsew")
        self.buttonframe.pack(after = self.dispradio)

        #display
        self.canvas1 = tk.Canvas(self.root, width=400, height=200, bd=0, highlightthickness=0, relief='ridge')
        self.canvas1.pack(after=self.buttonframe)

        # results
        self.lab_on = tk.Label(self.root, text="ON events: ", font=("Arial",12))
        self.lab_on.pack(after=self.buttonframe)
        self.lab_bkg = tk.Label(self.root, text="Expected bkg: ", font=("Arial",12))
        self.lab_bkg.pack(after=self.lab_on)

        self.root.mainloop()

    def stretch_image_1(self,event):
        ratio = event.width / event.height
        ratio_img = self.img1.size[0] / self.img1.size[1]
        if ratio > ratio_img:
            height = int(event.height)
            width = int(height*ratio_img)
        else:
            width = int(event.width)
            height = int(width/ratio_img)
        newimg = self.img1.resize((width,height))
        newimgtk = ImageTk.PhotoImage(newimg)
        self.canvas1.create_image(int(event.width/2),int(event.height/2), image = newimgtk, anchor='nw')


    def set_display(self,configu,iteratio):
        # read file
        filename = filedialog.askopenfilename()
        print("\n")
        print("File selected successfully: ",filename)
        print("Reading metadata... (it might require a few minutes)")
        if filename.split('.')[1] != 'fits':
            raise Exception("Sorry, incorrect file format (only '.fits' files are supported)")
        # create an optimization object
        self.opt_obj = selection.optimization.Opt(configu,iteratio)
        # set GW
        self.opt_obj.read_GW(filename) # create img
        outfigdir = self.opt_obj.outfig_skymap
        self.img1 = Image.open(outfigdir).resize((400,200))
        print("File successfully uploaded.")
        print("\n")
        print("Display legend:")
        print("In the image, you can see the GW shape in red with the 50% and 90% containing areas.")
        print("The greyish area all over the GW defines the ON region: a neutrino coming from here might be a neutrino emitted in corrispondence of the GW!!!")
        print("The blue points are all the events recorded between December 5 and 18, 2022.")
        print("The globe is colored in two (slightly) different colors: the part on top is the down-going sky.")

        self.display1 = ImageTk.PhotoImage(self.img1)
        self.canvas1.create_image(0,0, image = self.display1, anchor=NW)
        self.canvas1.bind('<Configure>', self.stretch_image_1)
        
        print("\n")
        print("Selection optimization:")
        print("Unfortunately, the down-going sky is highly contaminated by muons, so not all of those blue points are neutrinos.")
        print("But we can hope to find some by optimizing our cuts on 4 different variables: likelihood, tracklength, number of hits, and bdt score.")
        print("After having selected your variable(s), you may want to tell me how many combinations do I have to optimize.")
        print("Therefore, you may define the 'Phase space' parameter which tells me the fraction of total cuts to optimize, and the 'order' parameter wich helps me to understand what's the optimization order of the variables.")
        print("\n")
        print("Press 'Start' when you're ready.")

        
    
    @staticmethod
    def is_float(string):
        try:
            float(string)
            return True
        except ValueError:
            return False

    def set_display2(self):
        print("\n")
        text1 = self.txt1.get()
        if self.is_float(text1):
            text1 = float(self.txt1.get())
            if (text1>0.) & (text1<=1):
                text2 = self.txt2.get()
                if self.is_float(text2):
                    text2 = float(self.txt2.get())
                    if (text2>0.) & (text2<=1):
                        # take optimization variables
                        vars = list()
                        if self.hits.get():
                            vars.append("fitinf3")
                        if self.tracklength.get():
                            vars.append("fitinf10")
                        if self.likelihood.get():
                            vars.append("likelihood")
                        if self.bdtscore.get():
                            vars.append("bdt_score")
                        # check if vars
                        if vars:
                            print("Parameters successfully acquired. Starting the optimization...")
                            self.opt_obj.variables = vars #save variables into optimization object
                            self.opt_obj.stpoint = text1
                            self.opt_obj.ord = text2
                            # results
                            non, expbk, path = self.opt_obj.optimization_procedure() # optimize
                            # upload image
                            self.img2 = Image.open(path).resize((400,200))
                            self.display1 = ImageTk.PhotoImage(self.img2)
                            self.canvas1.create_image(0,0, image = self.display1, anchor=NW)
                            self.canvas1.bind('<Configure>', self.stretch_image_1)
                            #self.labeldisp = tk.Label(self.result, text = "gw_event_optimized", image = self.display2)
                            #self.labeldisp.grid(row=7, column=0, rowspan = 2, sticky="nsew")
                            # update labels
                            self.lab_on.config(text="On events: "+str(non))
                            self.lab_bkg.config(text="Expected bkg: "+str(expbk))
                            # self.lab_cond = tk.Label(self.result, text="3sigma significance: "+str(cond), font=("Arial",10))
                            # self.lab_on.grid(row=2, column=1)

                        else:
                            print("No variable selected. You need to select at least one variable before starting the optimization.")
                    else:
                        print("The order values must be within (0,1].")
                else:
                    print("Insert numeric order (0,1], please.")
            else:
                print("The phase space factor has to be within (0,1].")

        else:
            print("Insert numeric phase space (0,1], please")
