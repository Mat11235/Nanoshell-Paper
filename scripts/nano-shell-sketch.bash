#!/bin/bash
if [ "$1" == "" ]; then
    echo "Usage: bash $0 <rho>"
    exit
    fi
radius_ratio=$1
g++ ../src/eV2ex.cxx -o ../bin/eV2ex
g++ -Wall -I/usr/include/ -L/usr/local/lib ../src/rho2ome_sp.cxx -o ../bin/rho2ome_sp -lgsl -lgslcblas -lm -larmadillo
resu=($(./../bin/rho2ome_sp $1))


inner_circle_color=`./../bin/eV2ex ${resu[0]}`
echo ${resu[0]} $inner_circle_color
# Set the colors for the circles (in hex format)
outer_circle_color="#C0C0C0CC"    # Silver
# inner_circle_color="#00FF00"    # Green
contour_color="#000000"         # Black
light_color="#FFFFFF"           # White

# Set the image dimensions
image_width=600
image_height=500

# Calculate the necessary canvas size to accommodate the circles
canvas_width=$((image_width > image_height ? image_width : image_height))
canvas_height=$canvas_width

# Calculate the radii based on the given ratio
outer_radius=$(awk "BEGIN{ print $canvas_width / 2 }")
inner_radius=$(awk "BEGIN{ print $outer_radius * $radius_ratio }")

# Calculate the offset for centering the circles on the canvas
offset_x=$(($canvas_width / 2))
offset_y=$(($canvas_height / 2))

# Define the light radius and sigma for the first spot
light_radius=5
light_sigma=50

# Define the parameters for the second spot
second_spot_radius=$((outer_radius / 6))
second_spot_blur_radius=$((second_spot_radius / 2))

second_spot_x=$(($canvas_width  - 4*$second_spot_radius))
second_spot_y=$(($canvas_height - 3*$second_spot_radius))

# Generate the outer circle with lighting and contour
convert -size ${canvas_width}x${canvas_height} xc:none \
  -fill "$outer_circle_color" \
  -stroke "$contour_color" -strokewidth 3 \
  -draw "translate $offset_x,$offset_y circle 0,0 0,$outer_radius" \
  \( +clone -background "$light_color" -blur 0x80 -shadow ${light_radius}x${light_sigma}+0+0 +repage \) \
  -composite \
  outer_circle.png

# Generate the inner circle
convert -size ${canvas_width}x${canvas_height} xc:none \
  -fill "$inner_circle_color" \
  -stroke "$contour_color" -strokewidth 3 \
  -draw "translate $offset_x,$offset_y circle 0,0 0,$inner_radius" \
  inner_circle.png


# Generate the second spot as a solid white circle with blurred borders
convert -size ${canvas_width}x${canvas_height} xc:none \
  -fill "$light_color" \
  -draw "translate $second_spot_x,$second_spot_y circle 0,0 $second_spot_radius,0" \
  -blur 0x$second_spot_blur_radius \
  second_spot.png
  
# Composite the two circles and the second spot into a single image
# composite -gravity center inner_circle.png outer_circle.png output1.png

convert outer_circle.png '(' outer_circle.png -fuzz 1% ')' -channel Alpha -evaluate Multiply 4 -composite output1.png
convert second_spot.png '(' second_spot.png -blur 0x$second_spot_blur_radius -fuzz 1% ')' -channel Alpha -evaluate Multiply 3 -composite output2.png
composite -gravity center output2.png output1.png output0.png
convert output0.png '(' inner_circle.png  -fuzz 20% ')' -channel Alpha -evaluate Multiply .5 -composite output.png
composite second_spot.png output.png ultima.png
composite -gravity center ultima.png flare.png ../img/nanoshell.png

# Cleanup temporary files
 rm outer_circle.png inner_circle.png second_spot.png output2.png output1.png output0.png ultima.png output.png

echo "Image generated successfully!"
