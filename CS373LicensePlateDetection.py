from cProfile import label
import math
import sys
from pathlib import Path

from matplotlib import pyplot
from matplotlib.patches import Rectangle

# import our basic, light-weight png reader library
import imageIO.png

# this function reads an RGB color png file and returns width, height, as well as pixel arrays for r,g,b
def readRGBImageToSeparatePixelArrays(input_filename):

    image_reader = imageIO.png.Reader(filename=input_filename)
    # png reader gives us width and height, as well as RGB data in image_rows (a list of rows of RGB triplets)
    (image_width, image_height, rgb_image_rows, rgb_image_info) = image_reader.read()

    print("read image width={}, height={}".format(image_width, image_height))

    # our pixel arrays are lists of lists, where each inner list stores one row of greyscale pixels
    pixel_array_r = []
    pixel_array_g = []
    pixel_array_b = []

    for row in rgb_image_rows:
        pixel_row_r = []
        pixel_row_g = []
        pixel_row_b = []
        r = 0
        g = 0
        b = 0
        for elem in range(len(row)):
            # RGB triplets are stored consecutively in image_rows
            if elem % 3 == 0:
                r = row[elem]
            elif elem % 3 == 1:
                g = row[elem]
            else:
                b = row[elem]
                pixel_row_r.append(r)
                pixel_row_g.append(g)
                pixel_row_b.append(b)

        pixel_array_r.append(pixel_row_r)
        pixel_array_g.append(pixel_row_g)
        pixel_array_b.append(pixel_row_b)

    return (image_width, image_height, pixel_array_r, pixel_array_g, pixel_array_b)


# a useful shortcut method to create a list of lists based array representation for an image, initialized with a value
def createInitializedGreyscalePixelArray(image_width, image_height, initValue = 0):

    new_array = [[initValue for x in range(image_width)] for y in range(image_height)]
    return new_array

# --------------------------------------------------------------------------------------------------------------Helper Methods 

class Queue:
    def __init__(self):
        self.items = []

    def isEmpty(self):
        return self.items == []

    def enqueue(self, item):
        self.items.insert(0,item)

    def dequeue(self):
        return self.items.pop()

    def size(self):
        return len(self.items)

def computeRGBToGreyscale(pixel_array_r, pixel_array_g, pixel_array_b, image_width, image_height):
    greyscale_pixel_array = createInitializedGreyscalePixelArray(image_width, image_height)
    
    for row in range(image_height):
        for col in range(image_width):
            red = 0.299*pixel_array_r[row][col]
            green = 0.587*pixel_array_g[row][col]
            blue = 0.114*pixel_array_b[row][col]
            greyscale_pixel_array[row][col] = round(red+green+blue)
            
    return greyscale_pixel_array

def scaleTo0And255AndQuantize(pixel_array, image_width, image_height):
    returnArray = createInitializedGreyscalePixelArray(image_width, image_height)
    minVal = 255
    maxVal = 0

    for width in range(0,image_width):
        for height in range(image_height):
            if pixel_array[height][width] < minVal:
                minVal = pixel_array[height][width]
            elif pixel_array[height][width] > maxVal:
                maxVal = pixel_array[height][width]

    for width in range(0,image_width):
        for height in range(image_height):
            if maxVal-minVal != 0:

                sout = (pixel_array[height][width]-minVal) * (((255-0)/(maxVal-minVal)))

                if sout < 0:
                    returnArray[height][width] = 0
                elif sout > 255:
                    returnArray[height][width] = 255
                else:
                    returnArray[height][width] = round(sout)
    
    return returnArray

def computeStandardDeviationImage5x5(pixel_array, image_width, image_height):
    returnArray = createInitializedGreyscalePixelArray(image_width, image_height)

    for row in range(2,image_height-2):
        for col in range(2, image_width-2):
            valueList = []
            for horoCycle in range(-2,3):
                for vertCycle in range(-2,3):
                    pixelValue = pixel_array[row+vertCycle][horoCycle+col]
                    valueList.append(pixelValue)
            valueList = sorted(valueList)
            median = 0
            for i in range(len(valueList)):
                median += valueList[i]
            median = median/len(valueList)

            stdeviation = 0
            for i in range(len(valueList)):
                stdeviation += (valueList[i]-median)*(valueList[i]-median)
            
            stdeviation = math.sqrt((stdeviation/25))
            returnArray[row][col]=stdeviation
    return returnArray

def computeThresholdGE(pixel_array, threshold_value, image_width, image_height):
    
    for width in range(image_width):
        for height in range(image_height):
            pixel = pixel_array[height][width]
            if (pixel < threshold_value):
                pixel_array[height][width] = 0
            else:
                pixel_array[height][width] = 255
        
    return pixel_array

def computeDilation8Nbh3x3FlatSE(pixel_array, image_width, image_height):
    outputMatrix = [([0]*(image_width)) for i in range(image_height)]
    for row in range(1, image_height-1):
        for col in range(1, image_width-1):
            hit = 0
            for horozontalProbe in range(-1,2):
                for verticalProbe in range(-1,2):
                    if (pixel_array[row+horozontalProbe][col+verticalProbe] > 0):
                        hit = 1
            outputMatrix[row][col] = hit
    return outputMatrix

def computeErosion8Nbh3x3FlatSE(pixel_array, image_width, image_height):
    outputMatrix = [([0]*(image_width)) for i in range(image_height)]
    for row in range(1, image_height-1):
        for col in range(1, image_width-1):
            fit = 1
            for horozontalProbe in range(-1,2):
                for verticalProbe in range(-1,2):
                    if (pixel_array[row+horozontalProbe][col+verticalProbe] == 0):
                        fit = 0
            outputMatrix[row][col] = fit
    return outputMatrix

def breadthFirstSearch(pixel_array, outputMatrix, outputDict, visited,  i, j, image_height, image_width, label):
    queue = Queue()
    queue.enqueue([i, j])
    visited[i][j] = True
    outputDict[label] += 1

    x = [-1, 1, 0, 0]
    y = [0, 0, -1, 1]

    while (not queue.isEmpty()):
        location = queue.dequeue()
        row = location[0]
        col = location[1]
        outputMatrix[row][col] = label

        for k in range(4):
            vertProbe = row + x[k]
            horoProbe = col + y[k]

            if ((vertProbe >= 0) and (horoProbe >= 0) and (vertProbe < image_height) and (horoProbe < image_width) and (visited[vertProbe][horoProbe] == False) and (pixel_array[vertProbe][horoProbe] > 0)):
                queue.enqueue([vertProbe, horoProbe])
                visited[vertProbe][horoProbe] = True
                outputDict[label] += 1

def computeConnectedComponentLabeling(pixel_array, image_width, image_height):
    outputMatrix = [([0]*(image_width)) for i in range(image_height)]
    outputDict = {}
    visited = [[False for x in range(image_width)] for y in range(image_height)] 
    label = 0

    for i in range(image_height):
        for j in range(image_width):
            if ((visited[i][j] == False) & (pixel_array[i][j] > 0)):
                label += 1
                outputDict[label] = 0
                breadthFirstSearch(pixel_array, outputMatrix,
                                   outputDict, visited,  i, j, image_height, image_width, label)

    return (outputMatrix, outputDict)

# This is our code skeleton that performs the license plate detection.
# Feel free to try it on your own images of cars, but keep in mind that with our algorithm developed in this lecture,
# we won't detect arbitrary or difficult to detect license plates!
def main():

    command_line_arguments = sys.argv[1:]

    SHOW_DEBUG_FIGURES = True

    # this is the default input image filename
    input_filename = "numberplate1.png"

    if command_line_arguments != []:
        input_filename = command_line_arguments[0]
        SHOW_DEBUG_FIGURES = False

    output_path = Path("output_images")
    if not output_path.exists():
        # create output directory
        output_path.mkdir(parents=True, exist_ok=True)

    output_filename = output_path / Path(input_filename.replace(".png", "_output.png"))
    if len(command_line_arguments) == 2:
        output_filename = Path(command_line_arguments[1])


    # we read in the png file, and receive three pixel arrays for red, green and blue components, respectively
    # each pixel array contains 8 bit integer values between 0 and 255 encoding the color values
    (image_width, image_height, px_array_r, px_array_g, px_array_b) = readRGBImageToSeparatePixelArrays(input_filename)

    # setup the plots for intermediate results in a figure
    fig1, axs1 = pyplot.subplots(2, 2)
    axs1[0, 0].set_title('Input red channel of image')
    axs1[0, 0].imshow(px_array_r, cmap='gray')
    axs1[0, 1].set_title('Input green channel of image')
    axs1[0, 1].imshow(px_array_g, cmap='gray')
    axs1[1, 0].set_title('Input blue channel of image')
    axs1[1, 0].imshow(px_array_b, cmap='gray')


    # STUDENT IMPLEMENTATION here

    px_array = computeRGBToGreyscale(px_array_r, px_array_g, px_array_b, image_width, image_height)
    compute_array = scaleTo0And255AndQuantize(px_array, image_width, image_height)
    compute_array = computeStandardDeviationImage5x5(compute_array, image_width, image_height)
    compute_array = computeThresholdGE(compute_array, 60, image_width, image_height)
    for i in range(5):
        compute_array = computeDilation8Nbh3x3FlatSE(compute_array, image_width, image_height)
    for i in range(5):
        compute_array = computeErosion8Nbh3x3FlatSE(compute_array, image_width, image_height)
    (px_array_lables, px_size) = computeConnectedComponentLabeling(compute_array, image_width, image_height)
    sizeDict = sorted(px_size.items(), key=lambda item: item[1], reverse=True)

    top = 0
    bottom = image_height
    left = image_width
    right = 0
    detectPlate = True
    labelToDetect = 0

    while detectPlate:
        for y in range(image_height):
            for x in range(image_width):
                if px_array_lables[y][x] == sizeDict[labelToDetect][0]:
                    top = max(top, y)
                    bottom = min(bottom, y)
                    left = min(left, x)
                    right = max(right, x)
        if ((right-left)/(top-bottom)) > 1.5 and ((right-left)/(top-bottom)) < 5:
            detectPlate = False
        else:
            labelToDetect += 1


    # print(top,bottom,left,right,sizeDict[0][0])
    # for i in range(len(sizeDict)):



    # compute a dummy bounding box centered in the middle of the input image, and with as size of half of width and height
    bbox_min_x = left
    bbox_max_x = right
    bbox_min_y = bottom
    bbox_max_y = top





    # Draw a bounding box as a rectangle into the input image
    axs1[1, 1].set_title('Final image of detection')
    axs1[1, 1].imshow(px_array, cmap='gray')
    rect = Rectangle((bbox_min_x, bbox_min_y), bbox_max_x - bbox_min_x, bbox_max_y - bbox_min_y, linewidth=1,
                     edgecolor='g', facecolor='none')
    axs1[1, 1].add_patch(rect)



    # write the output image into output_filename, using the matplotlib savefig method
    extent = axs1[1, 1].get_window_extent().transformed(fig1.dpi_scale_trans.inverted())
    pyplot.savefig(output_filename, bbox_inches=extent, dpi=600)

    if SHOW_DEBUG_FIGURES:
        # plot the current figure
        pyplot.show()


if __name__ == "__main__":
    main()