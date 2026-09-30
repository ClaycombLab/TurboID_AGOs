from ij import IJ, ImagePlus
from ij import WindowManager as WM
from ij.plugin import ChannelSplitter
from ij.plugin.filter import ParticleAnalyzer as PA
from ij.measure import Measurements
from ij.measure import ResultsTable
from ij.plugin.frame import RoiManager
from ij.plugin import Duplicator
from ij import ImageStack

rm = RoiManager.getRoiManager()
if rm is None:
    rm = RoiManager()
    
strain = ("wtGFP","deltac","deltan")
image_path = ("C:\\Users\\Samantha\\Claycomb Lab Dropbox\\Claycomb Lab Team Folder\\Sebastian Fuentes\\Grad School\\Reasearch Project\\Images\\Mutant D2005 Images\\Confocal\\TIFFs\\July21\\")
thresholds = {"high":997, "low":482, "mid":996}
wormcount = (3, 6, 8)

#sets the Measurements we want to quantify
IJ.run("Set Measurements...", "area mean min integrated display label integrated decimal=3")

for s in range(3):
	for w in range(wormcount[s]):
		rm.reset()
		#opens the image to analyze (same file as in the w loop above but we want to make sure there isn't any carryover edits from that loop so we open it in a fresh window
		imp = IJ.openImage(image_path + strain[s] + "\\" + strain[s] + "_" + str(w+1) + ".tif")
		nSlices = imp.getNSlices()
		highStack = ImageStack(imp.getWidth(),imp.getHeight())
		lowStack = ImageStack(imp.getWidth(),imp.getHeight())
		#resets the ROI manager so we measure and process in the right places (ie not the ROI from any previous images we may have been working on)
		rm.reset()
		rm.resetMultiMeasureResults()
		IJ.run("Clear Results")
		title = imp.getTitle()
		#imp.show()
		for z in range(1, nSlices + 1):
			dupe=Duplicator().run(imp, 1, 1, z, z, 1, 1)
			dupe.setTitle("%s_z%02d" % (strain[s], z))
			#higher threshold range for granules (or nucleus in the case of HRDE)
			high=Duplicator().run(imp, 1, 1, z, z, 1, 1)
			hightitle=high.getTitle()
			#high.show()
			# Apply specific threshold to each channel
			IJ.setThreshold(high, thresholds["high"], 65535)
			# Use threshold to make mask of what is a granule
			IJ.run(high, "Convert to Mask", "")
			highStack.addSlice(high.getProcessor())
			# Analyze particles: adds mask to each particle, numbers it, and adds ROI
			IJ.run(high, "Analyze Particles...", " show=[Overlay Masks] add overlay composite")
			high.close()
			# Using ROIs generated from Analyze particles, Measure (uses settings from "Set Measurements" above)
			IJ.run(dupe, "Select All", "")
			rm.runCommand(dupe, "Measure")
		IJ.saveAs("Results", (image_path + strain[s] + "\\" + strain[s] + "_" + str(w+1) + "_higher.txt"))
		# Create ImagePlus from mask stack
		highMaskImp = ImagePlus(strain[s] + "_" + str(w+1) + "_higher_mask_stack",highStack)
		# Save mask stack
		IJ.saveAs(highMaskImp,"Tiff",image_path + strain[s] + "\\" + strain[s] + "_" + str(w+1) +"_higher_mask_stack.tif")
		highMaskImp.close()
		for z in range(1, nSlices + 1):
			dupe=Duplicator().run(imp, 1, 1, z, z, 1, 1)
			dupe.setTitle("%s_z%02d" % (strain[s], z))
			#higher threshold range for granules (or nucleus in the case of HRDE)
			low=Duplicator().run(imp, 1, 1, z, z, 1, 1)
			lowtitle=low.getTitle()
			#low.show()
			# Apply specific threshold to each channel
			IJ.setThreshold(low, thresholds["low"], thresholds["mid"])
			# Use threshold to make mask of what is a granule
			IJ.run(low, "Convert to Mask", "")
			lowStack.addSlice(low.getProcessor())
			# Analyze particles: adds mask to each particle, numbers it, and adds ROI
			IJ.run(low, "Analyze Particles...", " show=[Overlay Masks] add overlay composite")
			low.close()
			# Using ROIs generated from Analyze particles, Measure (uses settings from "Set Measurements" above)
			IJ.run(dupe, "Select All", "")
			rm.runCommand(dupe, "Measure")
		IJ.saveAs("Results", (image_path + strain[s] + "\\" + strain[s] + "_" + str(w+1) + "_lower.txt"))
		# Create ImagePlus from mask stack
		lowMaskImp = ImagePlus(strain[s] + "_" + str(w+1) + "_lower_mask_stack",lowStack)
		# Save mask stack
		IJ.saveAs(lowMaskImp,"Tiff",image_path + strain[s] + "\\" + strain[s] + "_" + str(w+1) +"_lower_mask_stack.tif")
		lowMaskImp.close()
IJ.selectWindow("Results")
IJ.run("Close")
IJ.selectWindow("ROI Manager")
IJ.run("Close")