# Import everything
import sys, tkFileDialog, os, matplotlib, bisect, math
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
from Tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
from matplotlib.ticker import ScalarFormatter

class MainWindow(Frame):
    '''Create and manage main window'''
    
    def __init__(self, root):
        '''Initialize the frame'''
        Frame.__init__(self, root)
        self.pack()
        
        # Store class variables
        self.root = root
        self.simulFile = None
        
        # Set up buttons
        self.createButtons()
        
        # If provided simulation file with sys.argv, load it now
        if len(sys.argv) > 1:
            self.simulFile = sys.argv[1]
            self.indexFile()
            self.outDir = os.path.split(self.simulFile)[0]
            
        else:
            # Get Snuppat directories
            self.__loadDirectories()
    
    def __loadDirectories(self):
        '''Store Snuppat directories if found'''
        
        # Look for snuppat directory
        parentDir = os.path.split(os.getcwd())[0]
        outDir = os.path.join(parentDir, "Snuppat", "output")
        
        # Store initial directory
        if os.path.exists(outDir):
            self.outDir = outDir
        else:
            self.outDir = os.getcwd()
    
    def createButtons(self):
        '''Set all buttons'''
        
        # Plot
        self.plotButton = Button(self, text = "PLOT", command = self.plotFile,
                                 state = DISABLED)
        self.plotButton.pack(fill = BOTH, expand = True)
        
        # Load
        self.loadButton = Button(self, text = "LOAD", command = self.loadFile)
        self.loadButton.pack(fill = BOTH, expand = True)
        
        # Quit
        self.quitButton = Button(self, text = "QUIT", command = self.quit)
        self.quitButton.pack(fill = BOTH, expand = True)
    
    def loadFile(self):
        '''Load file manager'''
        
        # Get file name
        newFile = tkFileDialog.askopenfilename(initialdir = self.outDir,
                               title = "Select output file")
        
        # If no file loaded, return silently without changing value
        if len(newFile) == 0:
            return
        
        # Remember directory
        self.outDir = os.path.split(newFile)[0]
        
        # Create file index
        self.simulFile = newFile
        self.indexFile()
    
    def indexFile(self):
        '''Create a reference index of the simulation file'''
        
        # Check if in unix
        inUnix = True if os.name == "posix" else False
        
        # Get total of models
        totModels = 0
        if inUnix:
            nL = 5000 # Try with nL lines
            
            # Look for last model
            for line in os.popen("tail -n {} {}".format(nL, self.simulFile)):
                if "Model" in line:
                    totModels = line
            
            # If no model found in nL lines, continue
            if totModels != 0:
                totModels = int(totModels.split()[-1])
        
        # Indicate loading
        print "Loading..."
        print
        
        # Open file
        sortedModels = list(); sortedAges = list()
        fileIndex = {}; iniModel = None; oldPrctg = 0
        with open(self.simulFile, "r") as fread:
            # Look for each line with "model" and get the position
            while True:
                currPos = fread.tell()
                line = fread.readline()
                
                # Exit if EOF
                if len(line) == 0:
                    break
                
                # Store model and position
                if "Model" in line:
                    lnlst = line.split()
                    
                    # Read modNum, mass and age
                    modNum = int(lnlst[-1])
                    mass = float(fread.readline().split()[-1])
                    age = 10**(float(fread.readline().split()[-1]) - 3)
                    
                    fileIndex[modNum] = (currPos, age)
                    sortedModels.append(modNum)
                    sortedAges.append(age)
                    
                    if iniModel is None:
                        iniModel = modNum
                
                # Print progress
                if totModels != iniModel:
                    prctg = float(modNum - iniModel)/(totModels - iniModel)
                    prctg *= 100
                    
                    if totModels > 0 and ((prctg - oldPrctg) > 1):
                        # Move in shell: up one line, back three columns
                        # and erase line
                        print "[1A",; print "[30D",; print "[K",
                        
                        # Write precentage.
                        print "Done {}%".format(int(prctg))
                        oldPrctg = prctg
        
        print "Loaded"
        self.fileIndex = fileIndex
        self.sortedModels = sortedModels
        self.sortedAges = sortedAges
        
        # Enable plot button
        self.plotButton["state"] = NORMAL
    
    def plotFile(self):
        '''The meat of the GUI, here the output is plotted'''
        
        # Create the object
        plotRoot = Tk(className = self.simulFile)
        plotRoot.protocol("WM_DELETE_WINDOW", plotRoot.quit)
        win = PlotWindow(plotRoot, self.fileIndex, self.sortedModels,
                self.sortedAges, self.simulFile)
        
        # Main loop
        plotRoot.mainloop()
        
        # Close
        win.fread.close()
        plotRoot.destroy()

class PlotWindow(Frame):
    '''Create and manage main window'''
    
    def __init__(self, root, fileIndex, sortedModels, sortedAges, simulFile):
        '''Initialize the frame'''
        Frame.__init__(self, root)
        self.pack()
        
        # Store class variables
        self.root = root
        self.fileIndex = fileIndex
        self.atrbFile = "atrb.sav"
        self.sortedModels = sortedModels
        self.sortedAges = sortedAges
        self.simulFile = simulFile
        self.pltAtrb = {}
        self.firstAge = sortedAges[0]
        self.searchFor = "model"
        
        self.currModelii = None
        self.modelStart = None
        self.plotMasses = None
        self.plotBordMass = None
        self.plotTemp = None
        self.plotRad = None
        self.plotListsOfData = None
        self.plotConvRegions = None
        self.ax = None
        self.ax2 = None
        
        # Open output file
        self.fread = open(self.simulFile, "r")
        
        # Load data directory
        self.getFiles()
        
        # Set up plot window and buttons
        self.createPlotWindow()
        self.createButtons()
    
    def createButtons(self):
        '''Set all buttons'''
        
        # Option frame
        optionFrame = Frame(self)
        optionFrame.pack(fill = BOTH, expand = True)
        
        # Axes ranges
        self.xRangeLabel = Label(optionFrame, text = "x range")
        self.xRangeLabel.grid(row = 0, column = 0)
        
        self.xRangeTxt = Entry(optionFrame)
        self.xRangeTxt.grid(row = 0, column = 1)
        
        self.yRangeLabel = Label(optionFrame, text = "y range")
        self.yRangeLabel.grid(row = 1, column = 0)
        
        self.yRangeTxt = Entry(optionFrame)
        self.yRangeTxt.grid(row = 1, column = 1)
        
        # Model and age
        self.modelLabel = Label(optionFrame, text = "Model")
        self.modelLabel.grid(row = 0, column = 2)
        
        self.modelTxt = Entry(optionFrame)
        self.modelTxt.grid(row = 0, column = 3)
        
        self.ageLabel = Label(optionFrame, text = "Age (ky)")
        self.ageLabel.grid(row = 1, column = 2)
        
        self.ageTxt = Entry(optionFrame)
        self.ageTxt.grid(row = 1, column = 3)
        
        # Add or remove element
        self.elementLabel = Label(optionFrame, text = "Element")
        self.elementLabel.grid(row = 0, column = 4)
        
        self.elementTxt = Entry(optionFrame)
        self.elementTxt.grid(row = 0, column = 5)
        
        # Show temperature, density, mesh, convective zones and grid
        self.tempCheck = Checkbutton(optionFrame, text = "Temperature",
                command = self.cTemp)
        self.tempCheck.grid(row = 2, column = 0)
        self.tempCheck.deselect()
        if self.showTemp:
            self.tempCheck.select()
        
        self.rhoCheck = Checkbutton(optionFrame, text = "Neutron density",
                command = self.cRho)
        self.rhoCheck.grid(row = 2, column = 1)
        self.rhoCheck.deselect()
        if self.showRho:
            self.rhoCheck.select()
        
        self.meshCheck = Checkbutton(optionFrame, text = "Mesh",
                command = self.cMesh)
        self.meshCheck.grid(row = 2, column = 2)
        self.meshCheck.deselect()
        if self.showMesh:
            self.meshCheck.select()
        
        self.conveCheck = Checkbutton(optionFrame, text = "Convective zones",
                command = self.cConv)
        self.conveCheck.grid(row = 2, column = 3)
        self.conveCheck.deselect()
        if self.showConve:
            self.conveCheck.select()
        
        self.gridCheck = Checkbutton(optionFrame, text = "Grid",
                command = self.cGrid)
        self.gridCheck.grid(row = 2, column = 4)
        self.gridCheck.deselect()
        if self.showGrid:
            self.gridCheck.select()
        
        # Dashed element lines
        self.dashCheck = Checkbutton(optionFrame, text = "Dashed elements",
                command = self.cDash)
        self.dashCheck.grid(row = 2, column = 5)
        self.dashCheck.deselect()
        if self.dashElements:
            self.dashCheck.select()
        
        # Frame for add and remove buttons
        addRmFrame = Frame(optionFrame)
        addRmFrame.grid(row = 1, column = 5)
        
        self.addButton = Button(addRmFrame, text = "Add", command =
                self.addElem)
        self.addButton.pack(side = LEFT, fill = BOTH, expand = True)
        
        self.rmButton = Button(addRmFrame, text = "Remove", command =
                self.rmElem)
        self.rmButton.pack(side = LEFT, fill = BOTH, expand = True)
        
        # Next and back frame
        prevNextFrame = Frame(self)
        prevNextFrame.pack(fill = BOTH, expand = True)
        
        # Next and back buttons
        self.prevButton = Button(prevNextFrame, text = "Previous model",
                command = self.prevModel)
        self.prevButton.pack(side = LEFT, fill = BOTH, expand = True)
        self.nextButton = Button(prevNextFrame, text = "Next model", command =
                self.nextModel)
        self.nextButton.pack(side = LEFT, fill = BOTH, expand = True)
        
        # Button frame
        butFrame = Frame(self)
        butFrame.pack(fill = BOTH, expand = True)
        
        # Update button
        self.updateButton = Button(butFrame, text = "UPDATE", command =
                self.update)
        self.updateButton.pack(fill = BOTH, expand = True)
        
        # Print button
        self.printButton = Button(butFrame, text = "PRINT", command =
                self.printB)
        self.printButton.pack(fill = BOTH, expand = True)
        
        # Save and erase frame
        saveEraseFrame = Frame(self)
        saveEraseFrame.pack(fill = BOTH, expand = True)
        
        # Save/Erase state buttons
        self.saveButton = Button(saveEraseFrame, text = "SAVE STATE",
                command = self.saveB)
        self.saveButton.pack(side = LEFT, fill = BOTH, expand = True)
        self.eraseButton = Button(saveEraseFrame, text = "ERASE STATE",
                command = self.eraseB)
        self.eraseButton.pack(side = LEFT, fill = BOTH, expand = True)
    
    def createPlotWindow(self):
        '''Create and set the plot window'''
        
        # Figure and initial plot
        self.fig = plt.Figure()
        plt.rcParams.update({"font.size": 14})
        self.setupAndPlot()
        
        # Show canvas
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill = BOTH, expand = True)
        toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        toolbar.update()
        self.canvas._tkcanvas.pack(fill = BOTH, expand = True)
        
        self.cid = self.fig.canvas.mpl_connect("draw_event", self.updateDraw)
    
    def getFiles(self):
        '''Load data directory and relevant files'''
        
        # Look for snuppat directory
        parentDir = os.path.split(os.getcwd())[0]
        dataDir = os.path.join(parentDir, "Snuppat", "data")
        
        # Store initial directory
        if not os.path.exists(dataDir):
            dataDir = tkFileDialog.askdirectory(title = "Select data directory")
        
        self.dataDir = dataDir
        self.speciesFile = os.path.join(self.dataDir, "species.dat")
    
    def __initiatePlotAttributes(self):
        '''Create placeholders for all plot attributes'''
        
        self.showTemp = True
        self.showRho = False
        self.showConve = True
        self.showMesh = True
        self.showGrid = False
        self.dashElements = False
        
        self.pltAtrb["model"] = None
        self.pltAtrb["age"] = None
        self.pltAtrb["xrange"] = None
        self.pltAtrb["yrange"] = None
        self.pltAtrb["elements"] = ["he4", "n14", "c12", "c13", "ne22", "h"]
        
        # If self.atrbFile exists, get values from there
        if os.path.isfile(self.atrbFile):
            fread = open(self.atrbFile, "r")
            lnlst = fread.readline().split()[0]
            self.showTemp = True if lnlst == "True" else False
            lnlst = fread.readline().split()[0]
            self.showRho = True if lnlst == "True" else False
            lnlst = fread.readline().split()[0]
            self.showConve = True if lnlst == "True" else False
            lnlst = fread.readline().split()[0]
            self.showMesh = True if lnlst == "True" else False
            lnlst = fread.readline().split()[0]
            self.showGrid = True if lnlst == "True" else False
            lnlst = fread.readline().split()[0]
            self.dashElements = True if lnlst == "True" else False
            
            self.pltAtrb["model"] = int(fread.readline())
            self.pltAtrb["xrange"] = map(float, fread.readline().split())
            self.pltAtrb["yrange"] = map(float, fread.readline().split())
            self.pltAtrb["elements"] = fread.readline().split()
            fread.close()
    
    def searchModel(self):
        '''Put reading head at the beginning of the requested model'''
        
        # If model and age are not specified, get first model
        if self.pltAtrb["model"] is None and self.pltAtrb["age"] is None:
            ii = 0
            
        elif self.searchFor == "model":
            # Look for next closest model
            ii = bisect.bisect_left(self.sortedModels, self.pltAtrb["model"])
            if ii > len(self.sortedModels) - 1:
                ii = len(self.sortedModels) - 1
            
        elif self.searchFor == "age":
            # Look for next closest age
            ii = bisect.bisect_left(self.sortedAges, self.pltAtrb["age"] +
                    self.firstAge)
            if ii > len(self.sortedAges) - 1:
                ii = len(self.sortedAges) - 1
        
        self.modelStart = self.fileIndex[self.sortedModels[ii]][0]
        self.currModelii = ii
        
        self.fread.seek(self.modelStart)
    
    def addElem(self):
        '''Add element to the plot'''
        
        tempElem = self.elementTxt.get()
        if tempElem in self.pltAtrb["elements"]:
            return
        
        # Check if element in the network
        if len(elementPos(self.speciesFile, tempElem)) == 0:
            return
        
        self.pltAtrb["elements"].append(tempElem)
        
        # Enable rmButton
        self.rmButton["state"] = NORMAL
        
        # Update
        self.ax.clear()
        if self.ax2 is not None: self.ax2.clear()
        self.setupAndPlot()
        self.canvas.draw()
    
    def rmElem(self):
        '''Remove element from the plot'''
        
        tempElem = self.elementTxt.get()
        if tempElem not in self.pltAtrb["elements"]:
            return
        
        self.pltAtrb["elements"].remove(tempElem)
        
        # Disable if only one element
        if len(self.pltAtrb["elements"]) <= 1:
            self.rmButton["state"] = DISABLED
        
        # Update
        self.ax.clear()
        if self.ax2 is not None: self.ax2.clear()
        self.setupAndPlot()
        self.canvas.draw()
    
    def nextModel(self):
        '''Search for next model'''
        
        self.currModelii += 1
        self.prevButton["state"] = NORMAL
        
        if self.currModelii > len(self.sortedModels) - 1:
            self.currModelii -= 1
            self.nextButton["state"] = DISABLED
        
        # Update values
        self.pltAtrb["model"] = self.sortedModels[self.currModelii]
        self.pltAtrb["age"] = self.sortedAges[self.currModelii] - self.firstAge
        
        # Update
        self.ax.clear()
        if self.ax2 is not None: self.ax2.clear()
        self.setupAndPlot()
        self.canvas.draw()
    
    def prevModel(self):
        '''Search for previous model'''
        
        self.currModelii -= 1
        self.nextButton["state"] = NORMAL
        
        if self.currModelii < 0:
            self.currModelii += 1
            self.prevButton["state"] = DISABLED
        
        # Update values
        self.pltAtrb["model"] = self.sortedModels[self.currModelii]
        self.pltAtrb["age"] = self.sortedAges[self.currModelii] - self.firstAge
        
        # Update
        self.ax.clear()
        if self.ax2 is not None: self.ax2.clear()
        self.setupAndPlot()
        self.canvas.draw()
    
    def cTemp(self):
        '''Toggle variable and change plot'''
        
        self.showTemp = not self.showTemp
        
        # Remove ax or deselect density
        if not self.showTemp:
            self.ax2.remove()
            self.ax2 = None
            
        else:
            self.showRho = False
            self.rhoCheck.deselect()
        
        # Update
        self.ax.clear()
        if self.ax2 is not None: self.ax2.clear()
        self.plotElements()
        self.canvas.draw()
    
    def cRho(self):
        '''Toggle variable and change plot'''
        
        self.showRho = not self.showRho
        
        # Remove ax or deselect temperature
        if not self.showRho:
            self.ax2.remove()
            self.ax2 = None
            
        else:
            self.showTemp = False
            self.tempCheck.deselect()
        
        # Update
        self.ax.clear()
        if self.ax2 is not None: self.ax2.clear()
        self.plotElements()
        self.canvas.draw()
    
    def cConv(self):
        '''Toggle variable and change plot'''
        
        self.showConve = not self.showConve
        
        # Update
        self.ax.clear()
        if self.ax2 is not None: self.ax2.clear()
        self.plotElements()
        self.canvas.draw()
    
    def cMesh(self):
        '''Toggle variable and change plot'''
        
        self.showMesh = not self.showMesh
        
        # Update
        self.ax.clear()
        if self.ax2 is not None: self.ax2.clear()
        self.plotElements()
        self.canvas.draw()
    
    def cGrid(self):
        '''Toggle variable and change plot'''
        
        self.showGrid = not self.showGrid
        
        # Update
        self.ax.clear()
        if self.ax2 is not None: self.ax2.clear()
        self.plotElements()
        self.canvas.draw()
    
    def cDash(self):
        '''Toggle variable and change plot'''
        
        self.dashElements = not self.dashElements
        
        # Update
        self.ax.clear()
        if self.ax2 is not None: self.ax2.clear()
        self.plotElements()
        self.canvas.draw()
    
    def updateDraw(self, event):
        '''Capture drawing event'''
        
        # Get current axis values in case they have changed
        self.pltAtrb["xrange"] = self.ax.get_xlim()
        self.pltAtrb["yrange"] = self.ax.get_ylim()
        
        # Write them into the entry
        self.xRangeTxt.delete(0, END)
        self.xRangeTxt.insert(0,
                "{:11.7f} {:11.7f}".format(*self.pltAtrb["xrange"]))
        self.yRangeTxt.delete(0, END)
        self.yRangeTxt.insert(0,
                "{:9.2e} {:9.2e}".format(*self.pltAtrb["yrange"]))
        
        # Now update model number and age
        self.modelTxt.delete(0, END)
        self.modelTxt.insert(0, "{}".format(self.pltAtrb["model"]))
        self.ageTxt.delete(0, END)
        self.ageTxt.insert(0, "{:.5f}".format(self.pltAtrb["age"]))
    
    def update(self):
        '''Read all the options and update the plot'''
        
        # xRange
        xRange = self.xRangeTxt.get().split()
        
        try:
            self.pltAtrb["xrange"] = [float(xRange[0]), float(xRange[1])]
        except ValueError:
            self.pltAtrb["xrange"] = None
        except IndexError:
            self.pltAtrb["xrange"] = None
        except:
            raise
        
        # yRange
        yRange = self.yRangeTxt.get().split()
        
        try:
            self.pltAtrb["yrange"] = [float(yRange[0]), float(yRange[1])]
        except ValueError:
            self.pltAtrb["yrange"] = None
        except IndexError:
            self.pltAtrb["yrange"] = None
        except:
            raise
        
        # Check if age or model has changed
        ageModelChanged = False
        try:
            modNumText = int(self.modelTxt.get())
            ageText = float(self.ageTxt.get())
            
            if modNumText != self.pltAtrb["model"]:
                self.pltAtrb["model"] = modNumText
                ageModelChanged = True
                self.searchFor = "model"
                
            elif abs(ageText - self.pltAtrb["age"]) > 1e-6:
                self.pltAtrb["age"] = ageText
                ageModelChanged = True
                self.searchFor = "age"
            
        except ValueError:
            pass
        except:
            raise
        
        # Update
        self.ax.clear()
        if self.ax2 is not None: self.ax2.clear()
        
        if ageModelChanged:
            self.setupAndPlot()
        else:
            self.plotElements()
        
        self.canvas.draw()
    
    def printB(self):
        '''Print the values currently shown in the window columnwise
        with the mass in the first column, temperature in the second,
        neutron density in the third and then the abundances'''
        
        # xRange and yRange
        xRange = [float(x) for x in self.xRangeTxt.get().split()]
        
        # Values
        mass = self.plotMasses
        temp = self.plotTemp
        rho = self.plotRho
        elements = self.plotListsOfData
        
        # Open output file
        filNam = os.path.split(self.simulFile)[-1] + "PRINT"
        filNam += str(self.pltAtrb["model"]) + str(xRange)
        fout = open(filNam, 'w')
        
        # Compose header
        header = "Mass | Temperature | Neutron Density"
        for elem in self.pltAtrb["elements"]:
            header += " | {}".format(elem)
        header += '\n'
        
        writtenHead = False
        
        # Search for mass coordinate
        for ii in range(len(self.plotMasses)):
            if mass[ii] < xRange[0]:
                continue
            
            if mass[ii] > xRange[1]:
                break
            
            # Write header
            if not writtenHead:
                fout.write(header)
                writtenHead = True
            
            # Compose line
            s = "{} {} {}".format(mass[ii], temp[ii], rho[ii])
            for elem in elements:
                s += " {}".format(elem[ii])
            s += '\n'
            
            # Write
            fout.write(s)
        
        # Close
        fout.close()
    
    def saveB(self):
        '''Save plot state such as elements, model and range'''
        
        # Open file
        atrbFile = open(self.atrbFile, "w")
        
        # Write options
        atrbFile.write("{}\n".format(self.showTemp))
        atrbFile.write("{}\n".format(self.showRho))
        atrbFile.write("{}\n".format(self.showConve))
        atrbFile.write("{}\n".format(self.showMesh))
        atrbFile.write("{}\n".format(self.showGrid))
        atrbFile.write("{}\n".format(self.dashElements))
        
        # Write attributes
        atrbFile.write("{}\n".format(self.pltAtrb["model"]))
        
        for elem in self.pltAtrb["xrange"]:
            atrbFile.write("{} ".format(elem))
        atrbFile.write("\n")
        
        for elem in self.pltAtrb["yrange"]:
            atrbFile.write("{} ".format(elem))
        atrbFile.write("\n")
        
        for elem in self.pltAtrb["elements"]:
            atrbFile.write("{} ".format(elem))
        atrbFile.write("\n")
        
        # Close
        atrbFile.close()
    
    def eraseB(self):
        '''Erase atrb.sav file'''
        
        if os.path.isfile(self.atrbFile):
            os.remove(self.atrbFile)

    def plotElements(self):
        '''Plot elements and convective regions'''
        
        # Plot preprocess
        if self.ax is None:
            self.ax = self.fig.add_subplot(1, 1, 1)
        
        # Manage plot attributes beforehand
        xLimit = self.pltAtrb["xrange"]
        if xLimit is None:
            xLimit = [self.plotMasses[0], self.plotMasses[-1]]
        
        yLimit = self.pltAtrb["yrange"]
        if yLimit is None:
            yLimit = [1e-24, 1]
        
        self.ax.set_xlabel("M/M$_\odot$")
        
        # Plot convective regions
        if self.showConve:
            for region in self.plotConvRegions:
                self.ax.fill_between(region, [yLimit[1], yLimit[1]],
                        color = "none", hatch = "/", edgecolor = "k")
        
        # Set xSpan
        xSpan = xLimit[1] - xLimit[0]
        
        # Plot chemistry
        self.ax.set_ylabel("Mass fraction")
        self.ax.set_yscale("log")
        self.ax.set_ylim(yLimit)
        
        lines = []
        for ii in range(len(self.pltAtrb["elements"])):
            name = self.pltAtrb["elements"][ii]
            capel = name[0].upper() + name[1:]
            if self.dashElements:
                lines.append(self.ax.plot(self.plotMasses,
                    self.plotListsOfData[ii], "--", label = capel, lw = 2))
            else:
                lines.append(self.ax.plot(self.plotMasses,
                    self.plotListsOfData[ii], label = capel, lw = 2))
        
        # Plot mesh
        if self.showMesh:
            constLevel = [math.sqrt(yLimit[0]*yLimit[1]) for x in
                    self.plotBordMass]
            lines.append(self.ax.plot(self.plotBordMass, constLevel, "r.",
                label = "Mesh", lw = 2))
        
        # Plot grid
        if self.showGrid:
            nHz = int(math.log10(yLimit[1]/yLimit[0]))
            for ii in range(nHz):
                val = 10**(int(math.log10(yLimit[0])) + ii)
                self.ax.plot(xLimit, (val, val), "r-")
                val = xLimit[0] + ii*(xLimit[1] - xLimit[0])/(nHz - 1)
                self.ax.plot((val, val), yLimit, "r-")
        
        # Plot temperature
        if self.showTemp:
            if self.ax2 is None:
                self.ax2 = self.ax.twinx()
            
            self.ax2.set_ylabel("Temperature (K)")
            self.ax2.set_yscale("log")
            self.ax2.set_ylim([1e5, 1e9])
            lines.append(self.ax2.plot(self.plotMasses, self.plotTemp, "k--",
                label = "Temperature", lw = 2))
            
        elif self.showRho:
            if self.ax2 is None:
                self.ax2 = self.ax.twinx()
            
            self.ax2.set_ylabel("Neutron density (n/cm$^3$)")
            self.ax2.set_yscale("log")
            lines.append(self.ax2.plot(self.plotMasses, self.plotRho, "k--",
                label = "Neutron density", lw = 2))
            
        
        # Legend
        lins = None
        for li in lines:
            if lins is None:
                lins = li
            else:
                lins += li
        
        labs = [l.get_label() for l in lins]
        self.ax.legend(lins, labs, prop = {"size": 10})
        
        # Xlimit
        self.ax.set_xlim(xLimit)
        self.ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=5))
        self.ax.xaxis.set_major_formatter(ScalarFormatter(useOffset = False))
    
    def setupAndPlot(self):
        '''Plot with self.pltAtrb attributes'''
        
        if len(self.pltAtrb) == 0:
            self.__initiatePlotAttributes()
        
        # Graphic elements and isotopes
        eleList = list()
        for elem in self.pltAtrb["elements"]:
            eleList.append(elementPos(self.speciesFile, elem))
        
        # Get next model
        self.searchModel()
        
        # Read the header
        self.pltAtrb["model"] = int(self.fread.readline().split()[-1])
        mass = float(self.fread.readline().split()[-1])
        self.pltAtrb["age"] = 10**(float(self.fread.readline().split()[-1]) - 3)
        self.pltAtrb["age"] -= self.firstAge
        
        # Initialize variables
        self.plotMasses = list()
        self.plotBordMass = list()
        self.plotTemp = list()
        self.plotRho = list()
        self.plotRad = list()
        
        # Initialize listsOfData
        self.plotListsOfData = list()
        for el in self.pltAtrb["elements"]:
            self.plotListsOfData.append(list())
        
        # Now let's go mass by mass
        lnlst1 = self.fread.readline().split()
        while True:
            line2 = self.fread.readline()
            lnlst2 = line2.split()
            if len(lnlst2) == 0:
                break
            
            # Check if in new model
            if "Model" in line2:
                break
            
            # Calculate and append values
            # Mass
            self.plotMasses.append((float(lnlst1[0]) +
                                    float(lnlst2[0]))*0.5*mass)
            self.plotBordMass.append(float(lnlst1[0])*mass)
            
            # Temperature
            valtemp = (float(lnlst2[1]) + float(lnlst1[1]))*5e8
            self.plotTemp.append(valtemp)
            
            # Neutron density
            valrho = (float(lnlst2[2]) + float(lnlst1[2]))*0.5
            valN = (float(lnlst2[4]) + float(lnlst1[4]))*0.5
            valrho *= 6.022e23*valN
            self.plotRho.append(valrho)
            
            # Nabla radiative - nabla adiabatic
            valRad = float(lnlst2[3]) + float(lnlst1[3])
            self.plotRad.append(valRad)
            
            # For each element or isotope, add all values corresponding to each
            # list. Treat effc13 separately
            for ii in range(len(self.pltAtrb["elements"])):
                val = 0
                if self.pltAtrb["elements"][ii] == "effc13":
                    c13, c13Mass = eleList[ii][0]
                    n14, n14Mass = eleList[ii][1]
                    
                    # Apply the definition XC13Eff = 13*(YC13 - YN14)
                    val = float(lnlst2[c13]) + float(lnlst1[c13])
                    val -= float(lnlst2[n14]) + float(lnlst1[n14])
                    val *= 0.5*c13Mass
                    
                else:
                    for posMass in eleList[ii]:
                        posEl, massEl = posMass
                        val += (float(lnlst2[posEl]) +
                                float(lnlst1[posEl]))*0.5*massEl
                
                self.plotListsOfData[ii].append(val)
            
            lnlst1 = lnlst2
        
        # Get convective regions
        self.plotConvRegions = list()
        inConvec = False; xConvLow = None; xConvHi = None
        for i in range(len(self.plotMasses)):
            # Look for first and last convective
            if self.plotRad[i] > 0 and not inConvec:
                # If more than one convective in a row
                if self.plotRad[i + 1] > 0:
                    xConvLow = self.plotMasses[i]
                    inConvec = True
            elif self.plotRad[i] <= 0 and inConvec:
                # If more than one radiative in a row
                if i < len(self.plotMasses) - 1 and self.plotRad[i + 1] <= 0:
                    xConvHi = self.plotMasses[i - 1]
                    inConvec = False
                    
                    # Add region
                    self.plotConvRegions.append([xConvLow, xConvHi])
        
        # Add last one
        if inConvec:
            self.plotConvRegions.append([xConvLow, self.plotMasses[-1]])
        
        # Plot elements
        self.plotElements()

def elementPos(dataFile, elem):
    '''Get element list of indices and masses'''
    
    # Be careful with effc13
    if elem == "effc13":
        c13 = elementPos(dataFile, "c13")
        n14 = elementPos(dataFile, "n14")
        
        return c13 + n14
    
    # Open file and read
    indicesMass = []
    with open(dataFile, "r") as fread:
        ii = 1
        for line in fread:
            lnlst = line.split()
            stri = "{}{}".format(lnlst[1], lnlst[0])
            elmMass = int(lnlst[0])
            
            # Convert hydrogen and neutrons
            if stri == "p1" or stri == "d2":
                lnlst[1] = "h"
            elif stri == "n1":
                lnlst[1] = "neut"
            
            # Now check if is exactly stri
            if elem == stri:
                indicesMass.append((ii - 1 + 4, elmMass))
                break
            
            # Check if its an element instead of isotope
            if elem == lnlst[1]:
                indicesMass.append((ii - 1 + 4, elmMass))
            
            ii += 1
    
    if len(indicesMass) == 0:
        print "Element not in the list. Please check {}.".format(dataFile)
    
    return indicesMass

def main():
    '''Main program controlling the GUI'''
    
    # Initialize
    root = Tk(className = "Snuppat GUI")
    root.protocol("WM_DELETE_WINDOW", root.quit)
    MainWindow(root)
    
    # Call main loop
    root.mainloop()
    
    # Close
    root.destroy()

if __name__ == "__main__":
    main()
