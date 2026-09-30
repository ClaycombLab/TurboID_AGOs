import os
from ij import IJ, ImagePlus
from ij import WindowManager as WM
from ij.plugin import ChannelSplitter
from ij.measure import Measurements
from ij.measure import ResultsTable
from ij.plugin.frame import RoiManager
from ij.plugin import Duplicator
from ij.gui import ShapeRoi
from ij.macro import Interpreter
import ij.ImageStack

# rm holds the Inner/Outer concentric-circle ROIs loaded from each image's
# overlay, and is also used to invoke generic commands (Measure, XOR).
rm = RoiManager.getRoiManager()

# Batch mode suppresses image windows from popping up/flashing on screen
# while the script runs; everything still processes normally underneath.
Interpreter.batchMode = True

image_path = "C:\\Users\\Samantha\\Claycomb Lab Dropbox\\Claycomb Lab Team Folder\\Julie Claycomb\\TURBO ID MANUSCRIPT 04_2026\\Primary Data\\WAGO-4_Imaging\\withROIs\\"

# Per-strain thresholds (the two strains were imaged with different
# settings, so their intensity ranges aren't comparable on one scale)
THRESHOLDS = {
	"wago-4": {"high": 14392, "low_min": 9766, "low_max": 14391},
	"wago-4_3xfg": {"high": 9966, "low_min": 4000, "low_max": 9965},
}

# Sets the Measurements we want to quantify
IJ.run("Set Measurements...", "area mean min integrated decimal=3")

# --- Discover all raw source images in the folder ---
# Raw images are named Strain_Treatment_Worm.tif (Strain may itself contain
# an underscore, e.g. "wago-4_3xfg"). Generated/output files from a previous
# run of this script contain "_higher" or "_lower" in the name, so they are
# excluded here to avoid re-processing them as if they were source images.
all_files = os.listdir(image_path)
image_files = [
	f for f in all_files
	if f.lower().endswith(".tif")
	and "_higher" not in f.lower()
	and "_lower" not in f.lower()
]

print("Found " + str(len(image_files)) + " image(s) to process:")
for f in image_files:
	print("  " + f)

for image_file in image_files:
	print("Processing " + image_file + " ...")

	# Reset the ROI manager so we measure/process in the right places
	# (i.e. not ROIs left over from a previous image)
	rm.reset()
	rm.resetMultiMeasureResults()
	# ResultsTable.reset() clears the table directly, rather than going
	# through "Clear Results" (a window-based command that can silently
	# fail to find/clear the Results window while running in batch mode)
	ResultsTable.getResultsTable().reset()

	# --- Parse Strain, Treatment, Worm from the filename, reading from the
	# end (since Strain can itself contain an underscore) ---
	file_stem = os.path.splitext(image_file)[0]
	parts = file_stem.split("_")
	n = len(parts)
	worm = parts[n - 1]
	treatment = parts[n - 2]
	strain = "_".join(parts[0:(n - 2)])

	# --- Look up this image's strain-specific thresholds ---
	if strain not in THRESHOLDS:
		print("  WARNING: no thresholds defined for strain '" + strain + "' - skipping " + image_file)
		continue
	THRESH_HIGH = THRESHOLDS[strain]["high"]
	THRESH_LOW_MIN = THRESHOLDS[strain]["low_min"]
	THRESH_LOW_MAX = THRESHOLDS[strain]["low_max"]

	# Opens the image to analyze
	imp = IJ.openImage(image_path + image_file)
	title = imp.getTitle()
	imp.show()

	# Loads the Inner/Outer concentric-circle ROIs from the image's overlay
	# into the ROI Manager (rm). ASSUMES ROIs alternate Inner, Outer, Inner,
	# Outer... in the order they were drawn (index 0 = Inner pair 1, index 1
	# = Outer pair 1, index 2 = Inner pair 2, etc.)
	IJ.run(imp, "To ROI Manager", "")

	roi_count = rm.getCount()
	print("  found " + str(roi_count) + " ROI(s) in overlay")
	num_pairs = roi_count // 2

	# Working duplicate used for the plain Inner measurement and for
	# computing the XOR selection
	dupe = Duplicator().run(imp, 1, 1, 1, 1, 1, 1)
	dupetitle = dupe.getTitle()

	for p in range(num_pairs):
		inner_idx = p * 2
		outer_idx = p * 2 + 1
		pair_num = p + 1
		print("  pair " + str(pair_num) + " of " + str(num_pairs))

		inner_roi = rm.getRoi(inner_idx)

		# --- Measure the Inner ROI directly (no thresholding) ---
		# Deselect anything selected in the ROI Manager first: if a
		# previous pair's XOR step left ROIs selected there, "Measure"
		# via the Manager would re-measure THOSE instead of dupe's
		# actual current selection. IJ.run(dupe, "Measure", "") measures
		# the image's active selection directly, sidestepping that.
		rm.deselect()
		dupe.setRoi(inner_roi)
		IJ.run(dupe, "Measure", "")
		rt = ResultsTable.getResultsTable()
		last_row = rt.size() - 1
		rt.setValue("Pair", last_row, pair_num)
		rt.setLabel("Inner", last_row)
		rt.setValue("Threshold", last_row, float('nan'))
		rt.updateResults()

		# --- The second ROI in each pair IS already the XOR/annulus shape
		# (drawn that way directly), so no XOR computation is needed -
		# just use it as-is ---
		xor_roi = rm.getRoi(outer_idx)
		dupe.setRoi(xor_roi)
		annulus_stats = dupe.getStatistics()
		full_annulus_area = annulus_stats.area
		print("    full annulus area: " + str(full_annulus_area))
		print("    annulus intensity - min: " + str(annulus_stats.min) + " mean: " + str(annulus_stats.mean) + " max: " + str(annulus_stats.max))

		# --- HIGHER threshold: select all thresholded pixels (whole image),
		# then intersect with the XOR annulus, then measure the REAL
		# (un-thresholded) intensity of that intersected region once ---
		high = Duplicator().run(imp, 1, 1, 1, 1, 1, 1)
		high.show()
		IJ.setThreshold(high, THRESH_HIGH, 65535)
		IJ.run(high, "Create Selection", "")
		thresh_roi = high.getRoi()
		high.close()

		rt = ResultsTable.getResultsTable()
		combined_roi = None
		if thresh_roi is not None:
			combined_roi = ShapeRoi(thresh_roi).and(ShapeRoi(xor_roi))
		if combined_roi is not None and combined_roi.getFloatWidth() > 0 and combined_roi.getFloatHeight() > 0:
			rm.deselect()
			dupe.setRoi(combined_roi)
			print("    HIGHER filtered area: " + str(dupe.getStatistics().area) + " (of " + str(full_annulus_area) + ")")
			IJ.run(dupe, "Measure", "")
			last_row = ResultsTable.getResultsTable().size() - 1
			rt.setValue("Pair", last_row, pair_num)
			rt.setLabel("XOR", last_row)
			rt.setValue("Threshold", last_row, 1)  # 1 = higher
		else:
			print("    HIGHER filtered area: 0 (of " + str(full_annulus_area) + ")")
			# No pixels above this threshold within the XOR annulus -
			# record a zero row so every pair still has exactly 3 rows
			# (Inner, XOR-higher, XOR-lower) for consistent downstream parsing
			rt.incrementCounter()
			last_row = rt.size() - 1
			rt.setValue("Area", last_row, 0)
			rt.setValue("Mean", last_row, 0)
			rt.setValue("Min", last_row, 0)
			rt.setValue("Max", last_row, 0)
			rt.setValue("IntDen", last_row, 0)
			rt.setValue("RawIntDen", last_row, 0)
			rt.setValue("Pair", last_row, pair_num)
			rt.setLabel("XOR", last_row)
			rt.setValue("Threshold", last_row, 1)
		rt.updateResults()


		# --- LOWER threshold: same pattern ---
		low = Duplicator().run(imp, 1, 1, 1, 1, 1, 1)
		low.show()
		IJ.setThreshold(low, THRESH_LOW_MIN, THRESH_LOW_MAX)
		IJ.run(low, "Create Selection", "")
		thresh_roi = low.getRoi()
		low.close()

		rt = ResultsTable.getResultsTable()
		combined_roi = None
		if thresh_roi is not None:
			combined_roi = ShapeRoi(thresh_roi).and(ShapeRoi(xor_roi))
		if combined_roi is not None and combined_roi.getFloatWidth() > 0 and combined_roi.getFloatHeight() > 0:
			rm.deselect()
			dupe.setRoi(combined_roi)
			print("    LOWER filtered area: " + str(dupe.getStatistics().area) + " (of " + str(full_annulus_area) + ")")
			IJ.run(dupe, "Measure", "")
			last_row = ResultsTable.getResultsTable().size() - 1
			rt.setValue("Pair", last_row, pair_num)
			rt.setLabel("XOR", last_row)
			rt.setValue("Threshold", last_row, 0)  # 0 = lower
		else:
			print("    LOWER filtered area: 0 (of " + str(full_annulus_area) + ")")
			rt.incrementCounter()
			last_row = rt.size() - 1
			rt.setValue("Area", last_row, 0)
			rt.setValue("Mean", last_row, 0)
			rt.setValue("Min", last_row, 0)
			rt.setValue("Max", last_row, 0)
			rt.setValue("IntDen", last_row, 0)
			rt.setValue("RawIntDen", last_row, 0)
			rt.setValue("Pair", last_row, pair_num)
			rt.setLabel("XOR", last_row)
			rt.setValue("Threshold", last_row, 0)
		rt.updateResults()

	# --- Save ONE combined Results file for this image: 3 rows per pair
	# (Inner, XOR-higher, XOR-lower), matching the Strain_Treatment_Worm
	# _Measurements naming your R script already parses ---
	IJ.saveAs("Results", (image_path + strain + "_" + treatment + "_" + worm + "_Measurements.txt"))
	print("  saved " + strain + "_" + treatment + "_" + worm + "_Measurements.txt")

	dupe.close()
	imp.close()

print("All " + str(len(image_files)) + " image(s) processed.")

# Turn batch mode back off and clean up any leftover windows
Interpreter.batchMode = False

# Results/ROI Manager windows may not exist if nothing was ever shown in
# batch mode, so these are wrapped defensively
try:
	IJ.selectWindow("Results")
	IJ.run("Close")
except Exception:
	pass
try:
	IJ.selectWindow("ROI Manager")
	IJ.run("Close")
except Exception:
	pass