Images and code for sulcal plotting
Download the sulcal plotting folder, and then source it in the top of your R or R markdown file such as with:
source("~/Downloads/sulcal_plotting/view_sulci.R")

Then you can use the function:
sulcal_map <- combine_and_color_sulci_parallel_snyder_2024(sulci_values, color_palette, scale_width = "9000x", border_width = 0, bounds = c(lower_bound,upper_bound)))

Where sulci_values is a 40-length vector of the values you want to plot, but must be in the order from the sulci_for_analysis.txt file set as default in the SPN container (copied into this folder as well for reference, should allow rerranging when reading data with either English or French sulcal labels).
Color palette is a choice from "magma", "bwr" (blue white red), "viridis", "plasma", "mako", "cbp" (cyan blue purple), and "ppy" (purple pink yellow) but you may also add in more color palettes to suit your scientific figure aesthetic needs.
scale_width set at 9000x approximately retains full size rendering, which I recommend, also with zero border width. Increasing border width can give a more cartoon-y aesthetic by putting black borders around the sulci plotted, and can also help colors pop more when the colors are light.
Bounds sets the color scale bounds for the data so that you can have control over an interpretable color scale given the range of values in your dataset / sulci_values

You may need to install a few R packages to get this script to work, namely imagemagick / magick which does most of the heavy lifting in rendering / combining images.
Also note that this function is parallelized and by default will take all available cores - 1, and runs locally for me on a Mac within ~30s to a minute. 

combine_and_color_sulci_alternate_parallel does the same thing but uses a slightly different sulcal parcellation (see README.txt in sulci_images_alternate)
