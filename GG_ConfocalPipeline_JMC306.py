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
strain = "JMC306"
RNAi = ("L4440","C27","D","F13","F46","Y39","Y110","Y111")
date = "Mar21"
image_path = "C:\\Users\\Samantha\\Claycomb Lab Dropbox\\Claycomb Lab Team Folder\\Julie Claycomb\\TurboID GENETICS MS\\Primary Data\\Confocal\\tetraspeck\\Mar21\\processing\\aligned\\"
#this should match the order that the images were taken (ie. which channels correspond to which proteins)
cnames = ("pgl-1","wago-4","znfx-1")
thresholds = [1550, 1700, 500] #Make sure you change the values in thresholds if you decide the ones we have set are too stringent or not stringent enough.
#They should match the thresholds used in the colocalization analysis below!

#sets the Measurements we want to quantify
IJ.run("Set Measurements...", "area mean min decimal=3")
#s loops through RNAi - currently have 8 strains to analyze corresponding to the 8 RNAis listed in the RNAi variable above. 
#Change the range value if the RNAi list gets shorter or longer
for s in range(8):
	#w loops through each worm imaged - should have 3 germlines to process for each strain
    for w in range(3):
    	#resets the ROI manager so we measure and process in the right places (ie not the ROI from any previous images we may have been working on)
    	RM.reset()
    	#opens the image to analyze. Make sure the path matches the naming convention you used for your .tif files from step 8 of the Benchling protocol
    	imp = IJ.openImage(image_path + RNAi[s] + "\\" + strain + "_" + RNAi[s] + "_" + date + "_w" + str(w+1) +".tif")
    	title = imp.getTitle()
    	imp.show()
    	#this takes the ROIs in the overlay and adds them to the ROI Manager for processing each nuclei individually
    	IJ.run("To ROI Manager", "")
    	#Runs the JACoP (Just Another Colocalization Plugin) tool from BIOP using manually set thresholds.
    	#Make sure you change the values in "manual_threshold_a" or "manual_threshold_b" if you decide the ones we have set are too stringent or not stringent enough. 
    	#They should match the thresholds variable above! 
    	#Compare channels 1 and 3
    	IJ.run(imp, "BIOP JACoP", "channel_a=1 channel_b=3 threshold_for_channel_a=[Use Manual Threshold Below] threshold_for_channel_b=[Use Manual Threshold Below] manual_threshold_a=1550 manual_threshold_b=500 crop_rois get_pearsons get_manders costes_block_size=5 costes_number_of_shuffling=100")
    	#saves the thresholded colocalization comparison images per nucleus
    	IJ.saveAs("Tiff", image_path+ RNAi[s] + "\\" + RNAi[s] +  "_" + strain + "_w" + str(w+1) + "_n1_C1_C3.tif")
    	IJ.run("Close")
    	IJ.saveAs("Tiff", image_path+ RNAi[s] + "\\" + RNAi[s] +  "_" + strain + "_w" + str(w+1) + "_n2_C1_C3.tif")
    	IJ.run("Close")
    	IJ.saveAs("Tiff", image_path+ RNAi[s] + "\\" + RNAi[s] +  "_" + strain + "_w" + str(w+1) + "_n3_C1_C3.tif")
    	IJ.run("Close")
    	#Runs JACoP for channels 1 and 2 and saves the comparison per nucleus
    	IJ.run(imp, "BIOP JACoP", "channel_a=1 channel_b=2 threshold_for_channel_a=[Use Manual Threshold Below] threshold_for_channel_b=[Use Manual Threshold Below] manual_threshold_a=1550 manual_threshold_b=1700 crop_rois get_pearsons get_manders costes_block_size=5 costes_number_of_shuffling=100")
    	IJ.saveAs("Tiff", image_path+ RNAi[s] + "\\" + RNAi[s] +  "_" + strain + "_w" + str(w+1) + "_n1_C1_C2.tif")
    	IJ.run("Close")
    	IJ.saveAs("Tiff", image_path+ RNAi[s] + "\\" + RNAi[s] +  "_" + strain + "_w" + str(w+1) + "_n2_C1_C2.tif")
    	IJ.run("Close")
    	IJ.saveAs("Tiff", image_path+ RNAi[s] + "\\" + RNAi[s] +  "_" + strain + "_w" + str(w+1) + "_n3_C1_C2.tif")
    	IJ.run("Close")
    	#Runs JACoP for channels 2 and 3 and saves the comparison per nucleus
        IJ.run(imp, "BIOP JACoP", "channel_a=2 channel_b=3 threshold_for_channel_a=[Use Manual Threshold Below] threshold_for_channel_b=[Use Manual Threshold Below] manual_threshold_a=1700 manual_threshold_b=500 crop_rois get_pearsons get_manders costes_block_size=5 costes_number_of_shuffling=100")
        IJ.saveAs("Tiff", image_path+ RNAi[s] + "\\" + RNAi[s] +  "_" + strain + "_w" + str(w+1) + "_n1_C2_C3.tif")
    	IJ.run("Close")
        IJ.saveAs("Tiff", image_path+ RNAi[s] + "\\" + RNAi[s] +  "_" + strain + "_w" + str(w+1) + "_n2_C2_C3.tif")
    	IJ.run("Close")
    	IJ.saveAs("Tiff", image_path+ RNAi[s] + "\\" + RNAi[s] +  "_" + strain + "_w" + str(w+1) + "_n3_C2_C3.tif")
    	IJ.run("Close")
    	#Saves the results of all the colocalization comparisons for the RNAi strain analyzed
        IJ.saveAs("Results", image_path+ RNAi[s] + "\\" + RNAi[s] +  "_" + strain + "_w" + str(w+1) + "_Colocs.txt")
        IJ.selectWindow("Results")
        IJ.run("Close")
        IJ.selectWindow(title)
        IJ.run("Close")
        #n loops through nuclei - 3 nuclei per germline
        for n in range(3):
        	#remove any ROIs in the ROI Manager
            RM.reset()
            #opens the image to analyze (same file as in the w loop above but we want to make sure there isn't any carryover edits from that loop so we open it in a fresh window
            imp = IJ.openImage(image_path + RNAi[s] + "\\" + strain + "_" + RNAi[s] + "_" + date + "_w" + str(w+1) +".tif")
            imp.show()
            #this takes the ROIs in the overlay and adds them to the ROI Manager for processing each nuclei individually
            IJ.run("To ROI Manager", "")
            #this essentially crops each ROI into an individual image
            rm.select(n)
            dupe=Duplicator().run(imp)
            #save the cropped ROIs
            IJ.saveAs(dupe, "Tiff", image_path + RNAi[s] + "\\" + RNAi[s] +  "_" + strain + "_w" + str(w+1) + "_n" + str(n+1) + ".tif")
            dupetitle=dupe.getTitle()
            dupe.show()
            # Split the channels
            channels = ChannelSplitter.split(dupe)
            raws = ChannelSplitter.split(dupe)
            #i loops through each channel (3 since there are 3 fluorescent proteins in JMC306)
            for i in range(3):
            	#clears the measurements window (removes any previous analysis)
                RM.resetMultiMeasureResults()
                RM.reset()
                og = raws[i]
                # Apply specific threshold to each channel
                IJ.setThreshold(channels[i], thresholds[i], 65535)
                # Use threshold to make mask of what is a granule
                IJ.run(channels[i], "Convert to Mask", "")
                # Analyze particles: adds mask to each particle, numbers it, and adds ROI
                p = IJ.run(channels[i], "Analyze Particles...", " show=[Overlay Masks] clear add overlay composite")
                # Save image of mask with numbered particles
                IJ.saveAs(channels[i],"Tiff",image_path + "\\" + RNAi[s] + "\\" + RNAi[s] + "_Overlay_Worm" + str(w+1) + "_n" + str(n+1) + "_" + cnames[i])
                # Using ROIs generated from Analyze particles, Measure (uses settings from "Set Measurements" above)
                IJ.selectWindow(dupetitle)
                IJ.run(og, "Select All","")
                rm.runCommand(og,"Measure")
                # Saves ROI measurements by worm, nucleus and channel
                IJ.saveAs("Results", (image_path + "\\" + RNAi[s] + "\\" + RNAi[s] + "_Measurements_Worm" + str(w+1) + "_n" + str(n+1) + "_" + cnames[i] + ".txt"))
                #just closing all the windows we opened
                IJ.selectWindow("Results")
                IJ.run("Close")
            IJ.selectWindow(dupetitle)
            IJ.run("Close")
            IJ.selectWindow(title)
            IJ.run("Close")
IJ.selectWindow("ROI Manager")
IJ.run("Close")



