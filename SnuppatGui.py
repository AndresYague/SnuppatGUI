# Import everything
import sys, os, matplotlib, bisect, math, struct, pickle, random
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
from matplotlib.ticker import ScalarFormatter
from tkinter import filedialog

class readBinaryModels(object):
    '''Class for reading binary models'''

    def __init__(self, fil):
        '''Initialize'''

        super(readBinaryModels, self).__init__()
        self.fread = open(fil, "rb")

    def close(self):
        '''Close file'''

        self.fread.close()

    def readHeader(self, currPos=None):
        '''Return header'''

        if currPos is not None:
            self.fread.seek(currPos)

        head = []
        byte = self.fread.read(4)
        if len(byte) == 0:
            return None

        head.append(*struct.unpack('i', byte))
        head.append(*struct.unpack('d', self.fread.read(8)))
        head.append(*struct.unpack('d', self.fread.read(8)))
        head.append(*struct.unpack('i', self.fread.read(4)))
        head.append(*struct.unpack('i', self.fread.read(4)))

        return head

    def nextModel(self, currPos=None):
        '''Return next model, unpacked'''

        if currPos is not None:
            self.fread.seek(currPos)

        # Read header
        head = self.readHeader()
        if head is None:
            return None

        model = [head]
        for ii in range(head[3]):
            s = []
            for jj in range(head[4]):
                s.append(*struct.unpack('d', self.fread.read(8)))

            model.append(s)

        return model

    def readOnlyHeader(self):
        '''Look only for the header and skip the rest'''

        # Read header
        currPos = self.fread.tell()
        head = self.readHeader()
        if head is None:
            return None, None

        # Skip file
        for ii in range(head[3]):
            for jj in range(head[4]):
                self.fread.read(8)

        return head, currPos

class MainWindow(Frame):
    '''Create and manage main window'''

    def __init__(self, root):
        '''Initialize the frame'''
        Frame.__init__(self, root)
        self.pack(fill=BOTH, expand=True)

        # Store class variables
        self.root = root
        self.fileIndex = None
        self.sortedModels = None
        self.sortedAges = None
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

    def __checkIndex(self, models, fileIndex, sortedModels):
        '''Check some random models to check that everything is ok'''

        # Choose either a 10% of models, 50, or all of them, whatever is
        # the minimum there are
        checkNMods = min([max([len(sortedModels)*0.1, 25]), len(sortedModels)])
        checkNMods = int(checkNMods)

        print(f"Checking {checkNMods} random models")
        for ii in range(checkNMods + 2):
            if ii < checkNMods:
                modNum = random.choice(sortedModels)
            elif ii == checkNMods:
                print("Checking first model")
                modNum = sortedModels[0]
            else:
                print("Checking last model")
                modNum = sortedModels[-1]

            tup = fileIndex[modNum]
            modHead = models.readHeader(tup[0])

            if modNum != modHead[0]:
                return False

        return True

    def createButtons(self):
        '''Set all buttons'''

        # Plot
        self.plotButton = Button(self, text="PLOT", command=self.plotFile,
                                 state=DISABLED)
        self.plotButton.pack(fill=BOTH, expand=True)

        # Load
        self.loadButton = Button(self, text="LOAD", command=self.loadFile)
        self.loadButton.pack(fill=BOTH, expand=True)

        # Reload
        self.reloadButton = Button(self, text="RELOAD", command=self.reload,
                                 state=DISABLED)
        self.reloadButton.pack(fill=BOTH, expand=True)

        # Quit
        self.quitButton = Button(self, text="QUIT", command=self.quit)
        self.quitButton.pack(fill=BOTH, expand=True)

    def loadFile(self):
        '''Load file manager'''

        # Get file name
        newFile = filedialog.askopenfilename(initialdir = self.outDir,
                               title = "Select output file")

        # If no file loaded, return silently without changing value
        if len(newFile) == 0:
            return

        # Remember directory
        self.outDir = os.path.split(newFile)[0]

        # Create file index
        self.simulFile = newFile
        self.indexFile()

    def reload(self):
        '''Reload the last loaded file'''

        self.indexFile()

    def indexFile(self):
        '''Create a reference index of the simulation file'''

        # Indicate loading
        print("Loading...")
        print()

        indDir = "indices"
        if not os.path.exists(indDir):
            os.mkdir(indDir)

        breakFilePath = os.path.split(self.simulFile)
        indicesFile = os.path.join(indDir, breakFilePath[1] + ".sav")

        # If it exists load it directly, otherwise create file
        if os.path.exists(indicesFile):
            with open(indicesFile, "rb") as fread:
                fileIndex = pickle.load(fread)
                sortedModels = pickle.load(fread)
                sortedAges = pickle.load(fread)

            # Check index
            models = readBinaryModels(self.simulFile)
            isFine = self.__checkIndex(models, fileIndex, sortedModels)
            models.close()

            if not isFine:
                # Remove file
                print("Error during check. Removing file and reloading")
                os.remove(indicesFile)

                # Now load the file again
                self.indexFile()
                return

            print("Loaded from file")

        else:

            # Open file
            loaded = 0
            sortedModels = list(); sortedAges = list(); fileIndex = {}
            models = readBinaryModels(self.simulFile)

            # Look for each line with "model" and get the position
            while True:
                head, currPos = models.readOnlyHeader()
                if head is None:
                    break

                # Store model and position
                loaded += 1

                # Assign modNum, mass and age
                modNum = head[0]
                mass = head[1]
                age = 10**(head[2] - 3)

                fileIndex[modNum] = (currPos, age)
                sortedModels.append(modNum)
                sortedAges.append(age)

            with open(indicesFile, "wb") as fwrite:
                pickle.dump(fileIndex, fwrite)
                pickle.dump(sortedModels, fwrite)
                pickle.dump(sortedAges, fwrite)

            models.close()

        print("Loaded")
        self.fileIndex = fileIndex
        self.sortedModels = sortedModels
        self.sortedAges = sortedAges

        # Enable plot and reload button
        self.plotButton["state"] = NORMAL
        self.reloadButton["state"] = NORMAL

    def plotFile(self):
        '''The meat of the GUI, here the output is plotted'''

        # Create the object
        plotRoot = Tk(className=self.simulFile)
        plotRoot.protocol("WM_DELETE_WINDOW", plotRoot.quit)
        win = PlotWindow(plotRoot, self.fileIndex, self.sortedModels,
                self.sortedAges, self.simulFile)

        # Main loop
        plotRoot.mainloop()

        # Close
        if win.loadedModels:
            win.models.close()

        plotRoot.destroy()

class PlotWindow(Frame):
    '''Create and manage main window'''

    def __init__(self, root, fileIndex, sortedModels, sortedAges, simulFile):
        '''Initialize the frame'''
        Frame.__init__(self, root)
        self.pack(fill=BOTH, expand=True)

        # Store class variables
        self.root = root
        self.fileIndex = fileIndex
        self.sortedModels = sortedModels
        self.sortedAges = sortedAges
        self.simulFile = simulFile

        self.atrbFile = "atrb.sav"
        self.pltAtrb = {}
        if sortedAges is not None:
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

        self.models = None
        self.loadedModels = False
        if self.simulFile is not None:
            self.models = readBinaryModels(self.simulFile)
            self.loadedModels = True

        # Load data directory
        self.getFiles()

        # Set up plot window and buttons
        self.createPlotWindow()
        self.createButtons()

    def createButtons(self):
        '''Set all buttons'''

        # Option frame
        optionFrame = Frame(self)
        optionFrame.pack(fill=BOTH, expand=True)

        # Axes ranges
        self.xRangeLabel = Label(optionFrame, text="x range")
        self.xRangeLabel.grid(row=0, column=0)

        self.xRangeTxt = Entry(optionFrame)
        self.xRangeTxt.grid(row=0, column=1)

        self.yRangeLabel = Label(optionFrame, text="y range")
        self.yRangeLabel.grid(row=1, column=0)

        self.yRangeTxt = Entry(optionFrame)
        self.yRangeTxt.grid(row=1, column=1)

        # Model and age
        self.modelLabel = Label(optionFrame, text="Model")
        self.modelLabel.grid(row=0, column=2)

        self.modelTxt = Entry(optionFrame)
        self.modelTxt.grid(row=0, column=3)

        self.ageLabel = Label(optionFrame, text="Age (ky)")
        self.ageLabel.grid(row=1, column=2)

        self.ageTxt = Entry(optionFrame)
        self.ageTxt.grid(row=1, column=3)

        # Add or remove element
        self.elementLabel = Label(optionFrame, text="Element")
        self.elementLabel.grid(row=0, column=4)

        self.elementTxt = Entry(optionFrame)
        self.elementTxt.grid(row=0, column=5)

        # Show temperature, mesh and convective zones
        self.tempCheck = Checkbutton(optionFrame, text = "Temperature",
                command = self.cTemp)
        self.tempCheck.grid(row=2, column=0)
        self.tempCheck.deselect()
        if self.showTemp:
            self.tempCheck.select()

        self.rhoCheck = Checkbutton(optionFrame, text = "Neutron density",
                command = self.cRho)
        self.rhoCheck.grid(row=2, column=1)
        self.rhoCheck.deselect()
        if self.showRho:
            self.rhoCheck.select()

        self.meshCheck = Checkbutton(optionFrame, text = "Mesh",
                command = self.cMesh)
        self.meshCheck.grid(row=2, column=2)
        self.meshCheck.deselect()
        if self.showMesh:
            self.meshCheck.select()

        self.conveCheck = Checkbutton(optionFrame, text = "Convective zones",
                command = self.cConv)
        self.conveCheck.grid(row=2, column=3)
        self.conveCheck.deselect()
        if self.showConve:
            self.conveCheck.select()

        self.gridCheck = Checkbutton(optionFrame, text = "Grid",
                command = self.cGrid)
        self.gridCheck.grid(row=2, column=4)
        self.gridCheck.deselect()
        if self.showGrid:
            self.gridCheck.select()

        # Dashed element lines
        self.dashCheck = Checkbutton(optionFrame, text = "Dashed elements",
                command = self.cDash)
        self.dashCheck.grid(row=2, column=5)
        self.dashCheck.deselect()
        if self.dashElements:
            self.dashCheck.select()

        # Frame for add and remove buttons
        addRmFrame = Frame(optionFrame)
        addRmFrame.grid(row=1, column=5)

        self.addButton = Button(addRmFrame, text = "Add", command =
                self.addElem)
        self.addButton.pack(side=LEFT, fill=BOTH, expand=True)

        self.rmButton = Button(addRmFrame, text = "Remove", command =
                self.rmElem)
        self.rmButton.pack(side=LEFT, fill=BOTH, expand=True)

        # Next and back frame
        prevNextFrame = Frame(self)
        prevNextFrame.pack(fill=BOTH, expand=True)

        # Next and back buttons
        self.prevButton = Button(prevNextFrame, text = "Previous model",
                command = self.prevModel)
        self.prevButton.pack(side=LEFT, fill=BOTH, expand=True)
        self.nextButton = Button(prevNextFrame, text = "Next model", command =
                self.nextModel)
        self.nextButton.pack(side=LEFT, fill=BOTH, expand=True)

        # Button frame
        butFrame = Frame(self)
        butFrame.pack(fill=BOTH, expand=True)

        # Update button
        self.updateButton = Button(butFrame, text = "UPDATE", command =
                self.update)
        self.updateButton.pack(fill=BOTH, expand=True)

        # Print button
        self.printButton = Button(butFrame, text = "PRINT", command =
                self.printB)
        self.printButton.pack(fill=BOTH, expand=True)

        # Save and erase frame
        saveEraseFrame = Frame(self)
        saveEraseFrame.pack(fill=BOTH, expand=True)

        # Save/Erase state buttons
        self.saveButton = Button(saveEraseFrame, text = "SAVE STATE",
                command = self.saveB)
        self.saveButton.pack(side=LEFT, fill=BOTH, expand=True)
        self.eraseButton = Button(saveEraseFrame, text = "ERASE STATE",
                command = self.eraseB)
        self.eraseButton.pack(side=LEFT, fill=BOTH, expand=True)

    def createPlotWindow(self):
        '''Create and set the plot window'''

        # Figure and initial plot
        self.fig = plt.Figure()
        plt.rcParams.update({"font.size": 14})
        self.setupAndPlot()

        # Show canvas
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(self.canvas, self)
        toolbar.update()
        self.canvas._tkcanvas.pack(fill=BOTH, expand=True)

        self.cid = self.fig.canvas.mpl_connect("draw_event", self.updateDraw)

    def getFiles(self):
        '''Load data directory and relevant files'''

        # Look for snuppat directory
        parentDir = os.path.split(os.getcwd())[0]
        dataDir = os.path.join(parentDir, "Snuppat", "data")

        # Store initial directory
        if not os.path.exists(dataDir):
            dataDir = filedialog.askdirectory(title="Select data directory")

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
            self.pltAtrb["xrange"] = list(map(float, fread.readline().split()))
            self.pltAtrb["yrange"] = list(map(float, fread.readline().split()))
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

    def addElem(self):
        '''Add element to the plot'''

        tempElem = self.elementTxt.get()
        if tempElem in self.pltAtrb["elements"]:
            return

        # Check if element in the network
        if len(elementPos(self.speciesFile, tempElem, 4)) == 0:
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

        # Get current axis values in case they have changed
        if self.loadedModels:

            # Now update model number and age
            self.modelTxt.delete(0, END)
            self.modelTxt.insert(0, f"{self.pltAtrb['model']}")
            self.ageTxt.delete(0, END)
            self.ageTxt.insert(0, f"{self.pltAtrb['age']:.5f}")

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

            modAtr = self.pltAtrb["model"]
            ageAtr = self.pltAtrb["age"]

            if modNumText != modAtr:
                modAtr = modNumText
                ageModelChanged = True
                self.searchFor = "model"

            elif abs(ageText - ageAtr) > 1e-6:
                ageAtr = ageText
                ageModelChanged = True
                self.searchFor = "age"

            self.pltAtrb["model"] = modAtr
            self.pltAtrb["age"] = ageAtr

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

        # Open it
        fout = open(filNam, 'w')

        # Compose header
        header = "Mass | Temperature | Neutron Density"
        for elem in self.pltAtrb["elements"]:
            header += f" | {elem}"
        header += '\n'

        writtenHead = False

        # Search for mass coordinate
        for ii in range(len(mass)):
            if mass[ii] < xRange[0]:
                continue

            if mass[ii] > xRange[1]:
                break

            # Write header
            if not writtenHead:
                fout.write(header)
                writtenHead = True

            # Compose line
            s = f"{mass[ii]} {temp[ii]} {rho[ii]}"
            for elem in elements:
                s += f" {elem[ii]}"
            s += '\n'

            # Write
            fout.write(s)

        # Close
        fout.close()

    def saveB(self):
        '''Save plot state such as elements, model and range'''

        # Open file
        atrbFile = open(self.atrbFile, "w")
        temp = self.showTemp
        conve = self.showConve
        mesh = self.showMesh
        dash = self.dashElements
        model = self.pltAtrb["model"]
        xRange = self.pltAtrb["xrange"]
        yRange = self.pltAtrb["yrange"]
        elements = self.pltAtrb["elements"]

        # Write options
        atrbFile.write(f"{self.showTemp} # Show temperature\n")
        atrbFile.write(f"{self.showRho} # Show neutron density\n")
        atrbFile.write(f"{self.showConve} # Show convection\n")
        atrbFile.write(f"{self.showMesh} # Show the mesh\n")
        atrbFile.write(f"{self.showGrid} # Show the grid\n")
        atrbFile.write(f"{self.dashElements} # Dash the elements\n")

        # Write attributes
        atrbFile.write(f"{model}\n")

        for elem in xRange:
            atrbFile.write(f"{elem} ")
        atrbFile.write("\n")

        for elem in yRange:
            atrbFile.write(f"{elem} ")
        atrbFile.write("\n")

        for elem in elements:
            atrbFile.write(f"{elem} ")
        atrbFile.write("\n")

        # Close
        atrbFile.close()

    def eraseB(self):
        '''Erase atrb.sav file'''

        if os.path.isfile(self.atrbFile):
            os.remove(self.atrbFile)

    def plotElements(self):
        '''Plot elements and convective regions'''

        # Initialize parameters
        xLimit = None
        yLimit = None
        lines = []

        # Plot preprocess
        if self.ax is None:
            self.ax = self.fig.add_subplot(1, 1, 1)

        self.ax.set_xlabel("M/M$_\odot$")
        self.ax.set_ylabel("Mass fraction")
        self.ax.set_yscale("log")

        # Manage plot attributes beforehand
        if self.loadedModels:
            xLimit = self.pltAtrb["xrange"]
            if xLimit is None:
                xLimit = [self.plotMasses[0], self.plotMasses[-1]]

            yLimit = self.pltAtrb["yrange"]
            if yLimit is None:
                yLimit = [1e-24, 1]

            # Plot convective regions
            if self.showConve:
                for region in self.plotConvRegions:
                    self.ax.fill_between(region, yLimit[0], yLimit[1],
                            facecolor = "none", hatch = "/", edgecolor = "k")

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

            # Plot neutron density
            elif self.showRho:
                if self.ax2 is None:
                    self.ax2 = self.ax.twinx()

                self.ax2.set_ylabel("Neutron density (n/cm$^3$)")
                self.ax2.set_yscale("log")
                lines.append(self.ax2.plot(self.plotMasses, self.plotRho, "k--",
                    label = "Neutron density", lw = 2))

        # Unpack lines
        lins = []
        for line in lines:
            lins += line
        lines = lins

        labs = [l.get_label() for l in lines]
        self.ax.legend(lines, labs, prop={"size": 10})

        # Limits
        self.ax.set_xlim(xLimit)
        self.ax.set_ylim(yLimit)
        self.ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=5))
        self.ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))

    def setupAndPlot(self):
        '''Plot with self.pltAtrb attributes'''

        if len(self.pltAtrb) == 0:
            self.__initiatePlotAttributes()

        # Graphic elements and isotopes
        eleList = list()
        for elem in self.pltAtrb["elements"]:
            eleList.append(elementPos(self.speciesFile, elem, 4))

        # Get next model and read header
        if self.loadedModels:
            self.searchModel()
            model = self.models.nextModel(self.modelStart)

            header = model[0]
            self.pltAtrb["model"] = header[0]
            mass = header[1]
            self.pltAtrb["age"] = 10**(header[2] - 3)
            self.pltAtrb["age"] -= self.firstAge

        if self.loadedModels:
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
            kk = 1
            while kk < header[3]:
                lin1 = model[kk]
                lin2 = model[kk + 1]

                # Calculate and append values
                # Mass
                self.plotMasses.append((lin1[0] + lin2[0])*0.5*mass)
                self.plotBordMass.append(lin1[0]*mass)

                # Temperature
                valtemp = (lin2[1] + lin1[1])*5e8
                self.plotTemp.append(valtemp)

                # Neutron density
                valrho = (lin2[2] + lin1[2])*0.5
                valN = (lin2[4] + lin1[4])*0.5
                valrho *= 6.022e23*valN
                self.plotRho.append(valrho)

                # Nabla radiative - nabla adiabatic
                valRad = (lin2[3] + lin1[3])*0.5
                self.plotRad.append(valRad)

                # For each element or isotope, add all values corresponding to
                # each list. Treat effc13 separately
                for ii in range(len(self.pltAtrb["elements"])):
                    val = 0
                    if self.pltAtrb["elements"][ii] == "effc13":
                        c13, c13Mass = eleList[ii][0]
                        n14, n14Mass = eleList[ii][1]

                        # Apply the definition XC13Eff = 13*(YC13 - YN14)
                        val = lin2[c13] + lin1[c13]
                        val -= lin2[n14] + lin1[n14]
                        val *= 0.5*c13Mass

                    else:
                        for posMass in eleList[ii]:
                            posEl, massEl = posMass
                            val += (lin2[posEl] + lin1[posEl])*0.5*massEl

                    self.plotListsOfData[ii].append(val)

                kk += 1

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

def elementPos(dataFile, elem, addIndx):
    '''Get element list of indices and masses'''

    # Be careful with effc13
    if elem == "effc13":
        c13 = elementPos(dataFile, "c13", addIndx)
        n14 = elementPos(dataFile, "n14", addIndx)

        return c13 + n14

    # Open file and read
    indicesMass = []
    with open(dataFile, "r") as fread:
        ii = 1
        for line in fread:
            lnlst = line.split()
            stri = f"{lnlst[1]}{lnlst[0]}"
            elmMass = int(lnlst[0])

            # Convert hydrogen and neutrons
            if stri == "p1" or stri == "d2":
                lnlst[1] = "h"
            elif stri == "n1":
                lnlst[1] = "neut"

            # Now check if is exactly stri
            if elem == stri:
                indicesMass.append((ii - 1 + addIndx, elmMass))
                break

            # Check if its an element instead of isotope
            if elem == lnlst[1]:
                indicesMass.append((ii - 1 + addIndx, elmMass))

            ii += 1

    if len(indicesMass) == 0:
        print(f"Element not in the list. Please check {dataFile}.")

    return indicesMass

def main():
    '''Main program controlling the GUI'''

    # Initialize
    root = Tk(className="Snuppat GUI")
    root.protocol("WM_DELETE_WINDOW", root.quit)
    MainWindow(root)

    # Call main loop
    root.mainloop()

    # Close
    root.destroy()

if __name__ == "__main__":
    main()
