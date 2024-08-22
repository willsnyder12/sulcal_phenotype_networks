
#install.packages("ggplot2")
library("ggplot2")

library("grid")
library("magrittr")
library("viridis")
#install.packages("raster")
library("raster")
library("readxl")
library(ggplot2)
library(viridisLite)
library(grid)
library("cowplot")
library("cowplot")
#install.packages("scales")
library("scales")

if (!requireNamespace("magick", quietly = TRUE)) {
  install.packages("magick")
}
library("magick")

library(foreach)
library(doParallel)



recolor_sulcus <- function(sulcus_image, color, fuzz = 0){
  return(image_transparent(image_background(image_transparent(sulcus_image, 'black',fuzz = fuzz), color), "white", fuzz = fuzz))
}
find_closest_index <- function(value, vector2) {
  abs_diff <- abs(value - vector2)
  return(which.min(abs_diff))
}
combine_and_color_sulci_alternate <- function(sulci_values, color_palette, scale_width, bounds = NA, border_width = 0){
  
  base_image_path <- "~/Downloads/sulcal_plotting/sulci_images_alternate/"
  new_names <- c("Lateral Fissure", "Anterior Cingulate Sulcus","Posterior Cingulate Sulcus", "Calcarine Fissure", "Collateral Fissure",
                 "Intra-Parietal Fissure", "Parieto-Occipital Fissure", "Occipital Sulcus", "Central Sulcus", "Inferior Frontal Sulcus", "Paracingulate Sulcus", "Intermediate Frontal Sulcus", "Superior Frontal Sulcus", "Occipital-Temporal Lateral Sulcus", "Orbitofrontal Sulcus", "Medial Parietal Sulcus", "Pre-Central Sulcus", "Post-Central Sulcus", "Inferior Temporal Sulcus", "Superior Temporal Sulcus"  )
  new_names_2 <- paste(rep(c("Left", "Right"), 20),rep(new_names, each = 2))
  new_names_2 <- new_names_2[-c(14)]
  new_names_2[40] <- "Right Parieto-Occipital Fissure"
  nomenclature <- new_names_2
  all_bg <- image_read("~/Downloads/sulcal_plotting/sulci_images_alternate/all_brains.png")
  sulci_files <- c(list.files("~/Downloads/sulcal_plotting/sulci_images_alternate", pattern = "^Left", full.names = TRUE),
                   list.files("~/Downloads/sulcal_plotting/sulci_images_alternate", pattern = "^Right", full.names = TRUE))
  composite.img <- all_bg
  
  remapped_sulci_values <- round(scales::rescale(sulci_values, to = c(1,1000)))
  if(!is.na(sum(bounds))){
    bounded_seq <- seq(from = bounds[1], to = bounds[2], length.out = 1000)
    remapped_sulci_values <- sapply(sulci_values, find_closest_index, bounded_seq)
  }
  
  if(color_palette == "magma"){
    hex_codes <- magma(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "bwr"){
    hex_codes <- colorRampPalette(c("blue","white", "red"))(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "viridis"){
    hex_codes <- viridis(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "plasma"){
    hex_codes <- plasma(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "mako"){
    hex_codes <- mako(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "cbp"){
    hex_codes <- colorRampPalette(c("cyan","blue", "purple"))(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "ppy"){
    hex_codes <- colorRampPalette(c("purple","deeppink", "yellow"))(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  
  for(i in 1:length(nomenclature)){
    tmp.img <- image_read(paste(base_image_path, nomenclature[i], ".png", sep = ""))
    tmp.img <- recolor_sulcus(tmp.img, sulci_hex_codes[i])
    if(border_width > 0){
      border.img <- image_morphology(image = tmp.img, operation = "dilate",
                                     kernel = "Disk", radius = border_width)
      tmp.img <- image_composite(border.img, tmp.img)
    }
    composite.img <- image_composite(composite.img,tmp.img)
  }
  return(image_scale(image_transparent(composite.img, "white"), scale_width))
}



combine_and_color_sulci_alternate_parallel <- function(sulci_values, color_palette, scale_width, bounds = NA, border_width = 0){
  
  base_image_path <- "~/Downloads/sulcal_plotting/sulci_images_alternate/"
  all_bg <- image_read("~/Downloads/sulcal_plotting/sulci_images_alternate/all_brains.png")
  new_names <- c("Lateral Fissure", "Anterior Cingulate Sulcus","Posterior Cingulate Sulcus", "Calcarine Fissure", "Collateral Fissure",
                 "Intra-Parietal Fissure", "Parieto-Occipital Fissure", "Occipital Sulcus", "Central Sulcus", "Inferior Frontal Sulcus", "Paracingulate Sulcus", "Intermediate Frontal Sulcus", "Superior Frontal Sulcus", "Occipital-Temporal Lateral Sulcus", "Orbitofrontal Sulcus", "Medial Parietal Sulcus", "Pre-Central Sulcus", "Post-Central Sulcus", "Inferior Temporal Sulcus", "Superior Temporal Sulcus"  )
  new_names_2 <- paste(rep(c("Left", "Right"), 20),rep(new_names, each = 2))
  new_names_2 <- new_names_2[-c(14)]
  new_names_2[40] <- "Right Parieto-Occipital Fissure"
  nomenclature <- new_names_2 # Assume this function generates your nomenclature
  remapped_sulci_values <- round(scales::rescale(sulci_values, to = c(1,1000)))
  if(!is.na(sum(bounds))){
    bounded_seq <- seq(from = bounds[1], to = bounds[2], length.out = 1000)
    remapped_sulci_values <- sapply(sulci_values, find_closest_index, bounded_seq)
  }
  
  if(color_palette == "magma"){
    hex_codes <- magma(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "bwr"){
    hex_codes <- colorRampPalette(c("blue","white", "red"))(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "viridis"){
    hex_codes <- viridis(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "plasma"){
    hex_codes <- plasma(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "mako"){
    hex_codes <- mako(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "cbp"){
    hex_codes <- colorRampPalette(c("cyan","blue", "purple"))(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "ppy"){
    hex_codes <- colorRampPalette(c("purple","deeppink", "yellow"))(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  no_cores <- detectCores() - 1
  registerDoParallel(cores=no_cores)
  # Use foreach for parallel processing
  tmp_imgs <- foreach(i = 1:length(nomenclature), .packages = "magick") %dopar% {
    file_path <- paste(base_image_path, nomenclature[i], ".png", sep = "")
    tmp_img <- image_read(file_path)
    tmp_img <- recolor_sulcus(tmp_img, sulci_hex_codes[i])
    if(border_width > 0){
      border_img <- image_morphology(image = tmp_img, operation = "dilate", kernel = "Disk", radius = border_width)
      tmp_img <- image_composite(border_img, tmp_img)
    }
    image_write(tmp_img, format = "png")  # Write to a temporary file and return the file path
  }
  
  # Composite all images onto the base image
  composite_img <- all_bg
  for (img_path in tmp_imgs) {
    tmp_img <- image_read(img_path)
    composite_img <- image_composite(composite_img, tmp_img)
  }
  
  # Stop the parallel backend
  stopImplicitCluster()
  
  return(image_scale(image_transparent(composite_img, "white"), scale_width))
}




combine_and_color_sulci_parallel_snyder_2024 <- function(sulci_values, color_palette, scale_width, bounds = NA, border_width = 0){
  
  base_image_path <- "~/Downloads/sulcal_plotting/sulci_images_snyder_2024/"
  all_bg <- image_read("~/Downloads/sulcal_plotting/sulci_images_snyder_2024/all_brains.png")
  new_names <- c("Lateral Fissure", "Anterior Cingulate Sulcus","Posterior Cingulate Sulcus", "Calcarine Fissure", "Collateral Fissure",
                 "Intra-Parietal Fissure", "Parieto-Occipital Fissure", "Occipital Sulcus", "Central Sulcus", "Inferior Frontal Sulcus", "Paracingulate Sulcus", "Intermediate Frontal Sulcus", "Superior Frontal Sulcus", "Occipital-Temporal Lateral Sulcus", "Orbitofrontal Sulcus", "Medial Parietal Sulcus", "Pre-Central Sulcus", "Post-Central Sulcus", "Inferior Temporal Sulcus", "Superior Temporal Sulcus"  )
  new_names_2 <- paste(rep(c("Left", "Right"), 20),rep(new_names, each = 2))
  new_names_2 <- new_names_2[-c(14)]
  new_names_2[40] <- "Right Parieto-Occipital Fissure"
  nomenclature <- new_names_2 # Assume this function generates your nomenclature
  remapped_sulci_values <- round(scales::rescale(sulci_values, to = c(1,1000)))
  if(!is.na(sum(bounds))){
    bounded_seq <- seq(from = bounds[1], to = bounds[2], length.out = 1000)
    remapped_sulci_values <- sapply(sulci_values, find_closest_index, bounded_seq)
  }
  
  if(color_palette == "magma"){
    hex_codes <- magma(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "bwr"){
    hex_codes <- colorRampPalette(c("blue","white", "red"))(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "viridis"){
    hex_codes <- viridis(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "plasma"){
    hex_codes <- plasma(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "mako"){
    hex_codes <- mako(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "cbp"){
    hex_codes <- colorRampPalette(c("cyan","blue", "purple"))(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  if(color_palette == "ppy"){
    hex_codes <- colorRampPalette(c("purple","deeppink", "yellow"))(1000)
    sulci_hex_codes <- hex_codes[remapped_sulci_values]
  }
  no_cores <- detectCores() - 1
  registerDoParallel(cores=no_cores)
  # Use foreach for parallel processing
  tmp_imgs <- foreach(i = 1:length(nomenclature), .packages = "magick") %dopar% {
    file_path <- paste(base_image_path, nomenclature[i], ".png", sep = "")
    tmp_img <- image_read(file_path)
    tmp_img <- recolor_sulcus(tmp_img, sulci_hex_codes[i])
    if(border_width > 0){
      border_img <- image_morphology(image = tmp_img, operation = "dilate", kernel = "Disk", radius = border_width)
      tmp_img <- image_composite(border_img, tmp_img)
    }
    image_write(tmp_img, format = "png")  # Write to a temporary file and return the file path
  }
  
  # Composite all images onto the base image
  composite_img <- all_bg
  for (img_path in tmp_imgs) {
    tmp_img <- image_read(img_path)
    composite_img <- image_composite(composite_img, tmp_img)
  }
  
  # Stop the parallel backend
  stopImplicitCluster()
  
  return(image_scale(image_transparent(composite_img, "white"), scale_width))
}
#sulcal_map <- combine_and_color_sulci_parallel(1:40, "ppy", scale_width = "9000x", border_width = 0, bounds = c(1,40))
#ggplot()  + ggtitle("Sulcal plot in R!") + annotation_custom(rasterGrob(sulcal_map)) + theme_minimal() 
