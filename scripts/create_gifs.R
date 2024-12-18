library(magick)

## list file names and read in
dir_out <- grep("CULT", dir("plots/landuse"), value = TRUE)
dir_out <- paste0("plots/landuse/CULT")
imgs <- list.files(dir_out, full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 10)

## view animated image
img_animated

## save to disk
image_write(image = img_animated,
            path = "CULT.gif")

## list file names and read in
dir_out <- paste0("plots/landuse/PASTPLANT")
imgs <- list.files(dir_out, full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 10)

## view animated image
img_animated

## save to disk
image_write(image = img_animated,
            path = "PASTPLANT.gif")
