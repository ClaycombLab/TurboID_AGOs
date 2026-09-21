from ij import IJ, ImagePlus
from ij import WindowManager as WM
from ij.plugin import ChannelSplitter
from ij.plugin.filter import ParticleAnalyzer as PA
from ij.measure import Measurements
from ij.measure import ResultsTable
from ij.plugin.frame import RoiManager
from ij.plugin import Duplicator
import ij.ImageStack

rm = RoiManager.getRoiManager()
RM = RoiManager()

# Define the path to your image file
strain = ("CSR-1","HRDE-1","WAGO-1","WAGO-4","PRG-1","D2005","PPW-2")
image_path = "C:\\Users\\Samantha\\Claycomb Lab Dropbox\\Claycomb Lab Team Folder\\Sebastian Fuentes\\Grad School\\Reasearch Project\\Images\\GFP_AGO\\TIFFs\\"
thresholds = {"CSR-1":1510,"lCSR-1":546,"hCSR-1":1509, 
		"HRDE-1":594, "lHRDE-1":418, "hHRDE-1":593,
		"WAGO-1":1092,"lWAGO-1":434,"hWAGO-1":1091,
		"WAGO-4":1927,"lWAGO-4":964,"hWAGO-4":1926,
		"PRG-1":1012,"lPRG-1":353,"hPRG-1":1011,
		"D2005":1172,"lD2005":755,"hD2005":1171,
		"PPW-2":1092,"lPPW-2":482,"hPPW-2":1091}


#sets the Measurements we want to quantify
IJ.run("Set Measurements...", "area mean min integrated decimal=3")

for s in range(7):
	for w in range(3):
		#resets the ROI manager so we measure and process in the right places (ie not the ROI from any previous images we may have been working on)
		RM.reset()
		RM.resetMultiMeasureResults()
		IJ.run("Clear Results")
		#opens the image to analyze (same file as in the w loop above but we want to make sure there isn't any carryover edits from that loop so we open it in a fresh window
		imp = IJ.openImage(image_path + strain[s] + "\\" + strain[s] + "_" + str(w+1) + ".tif")
		title = imp.getTitle()
		imp.show()
		dupe=Duplicator().run(imp, 1, 1, 1, 1, 1, 1)
		dupetitle=dupe.getTitle()
		#higher threshold range for granules (or nucleus in the case of HRDE)
		high=Duplicator().run(imp, 1, 1, 1, 1, 1, 1)
		hightitle=high.getTitle()
		high.show()
		# Apply specific threshold to each channel
		IJ.setThreshold(high, thresholds[strain[s]], 65535)
		# Use threshold to make mask of what is a granule
		IJ.run(high, "Convert to Mask", "")
		# Analyze particles: adds mask to each particle, numbers it, and adds ROI
		IJ.run(high, "Analyze Particles...", " show=[Overlay Masks] clear add overlay composite")
		# Save image of mask with numbered particles
		IJ.saveAs(high,"Tiff",image_path + strain[s] + "\\" + strain[s] + "_" + str(w+1) + "_higher.tif")
		high.close()
		# Using ROIs generated from Analyze particles, Measure (uses settings from "Set Measurements" above)
		IJ.run(dupe, "Select All","")
		rm.runCommand(dupe,"Measure")
		# Saves ROI measurements by worm, nucleus and channel
		IJ.saveAs("Results", (image_path + strain[s] + "\\" + strain[s] + "_" + str(w+1) + "_higher.txt"))
		low=Duplicator().run(imp, 1, 1, 1, 1, 1, 1)
		lowtitle=low.getTitle()
		low.show()
		# Apply specific threshold to each channel
		IJ.setThreshold(low, thresholds["l"+strain[s]], thresholds["h"+strain[s]])
		# Use threshold to make mask of what is a granule
		IJ.run(low, "Convert to Mask", "")
		# Analyze particles: adds mask to each particle, numbers it, and adds ROI
		IJ.run(low, "Analyze Particles...", " show=[Overlay Masks] clear add overlay composite")
		# Save image of mask with numbered particles
		IJ.saveAs(low,"Tiff",image_path + strain[s] + "\\" + strain[s] + "_" + str(w+1) + "_lower.tif")
		low.close()
		# Using ROIs generated from Analyze particles, Measure (uses settings from "Set Measurements" above)
		IJ.run(dupe, "Select All","")
		rm.runCommand(dupe,"Measure")
		# Saves ROI measurements by worm, nucleus and channel
		IJ.saveAs("Results", (image_path + strain[s] + "\\" + strain[s] + "_" + str(w+1) + "_lower.txt"))
		imp.close()
		
IJ.selectWindow("Results")
IJ.run("Close")
IJ.selectWindow("ROI Manager")
IJ.run("Close")
