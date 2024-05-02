import os
from array import array

def writeBOV(outDir, variableName, BOVname, formatstring, ndarray):
    dataSizeY = len(ndarray);
    dataSizeX = len(ndarray[0])
    dataSizeZ = len(ndarray[0][0])

    if formatstring=='i':
        dataFormat = "INT"
    elif formatstring=='f':
        dataFormat = "FLOAT"

    divideBrick = 0

    # write header file
    fid = open(os.path.join(outDir, BOVname + '.bov'), 'w')
    fid.write('TIME: 00.00\n')
    fid.write('DATA_FILE: ' + BOVname + '.dat \n')
    fid.write('# The data size corresponds to NX,NY,NZ \n')
    fid.write('DATA_SIZE: ' + str(dataSizeY) + ' ' + str(dataSizeX) + ' ' + str(dataSizeZ) + '\n')
    fid.write('# Allowable values for DATA_FORMAT are: BYTE,SHORT,INT,FLOAT,DOUBLE \n')
    fid.write('DATA_FORMAT: ' + dataFormat + '\n')
    fid.write('VARIABLE: ' + variableName + '\n')
    fid.write('# Endian representation of the computer that created the data.\n')
    fid.write('# Intel is LITTLE, many other processors are BIG.\n')
    fid.write('DATA_ENDIAN: LITTLE\n')
    fid.write('# Centering refers to how the data is distributed in a cell. If you \n')
    fid.write('# give zonal then its 1 data value per zone. Otherwise the data \n')
    fid.write('# will be centered at the nodes. \n')
    fid.write('CENTERING: zonal \n')
    fid.write('# BRICK_ORIGIN lets you specify a new coordinate system origin for \n')
    fid.write('# the mesh that will be created to suit your data. \n')
    fid.write('BRICK_ORIGIN: 0. 0. 0.\n')
    fid.write('# BRICK_SIZE lets you specify the size of the brick. \n')
    fid.write('BRICK_SIZE: ' + str(dataSizeY) + ' ' + str(dataSizeX) + ' ' + str(dataSizeZ) + '\n')
    fid.write('# BYTE_OFFSET: is optional and lets you specify some number of \n')
    fid.write('# bytes to skip at the front of the file. This can be useful for \n')
    fid.write('# skipping the 4-byte header that Fortran tends to write to files. \n')
    fid.write('# If your file does not have a header then DO NOT USE BYTE_OFFSET. \n')
    fid.write('# BYTE_OFFSET: 4 \n')
    fid.write('# DIVIDE_BRICK: is optional and can be set to true or false. \n')
    fid.write('# When DIVIDE_BRICK is true, the BOV reader uses the values stored \n')
    fid.write('# in DATA_BRICKLETS to divide the data into chunks that can be \n')
    fid.write('# processed in parallel. \n')
    fid.write('DIVIDE_BRICK: false \n')
    fid.write('# DATA_BRICKLETS: is optional and requires you to specify 3 integers \n')
    fid.write('# that indicate the size of the bricklets to create when you have \n')
    fid.write('# also specified the DIVIDE_BRICK option. The values chosen for \n')
    fid.write('# DATA_BRICKLETS must be factors of the numbers used for DATA_SIZE. \n')

    fid.write('#DATA_BRICKLETS: ' + str(dataSizeY / 10) + ' ' + str(dataSizeX / 10) + ' ' + str(dataSizeZ / 10) + '\n')
    fid.write('# DATA_COMPONENTS: is optional and tells the BOV reader how many \r\n')
    fid.write('# components your data has. 1=scalar, 2=complex number, 3=vector, \r\n')
    fid.write('# 4 and beyond indicate an array variable. You can use COMPLEX \r\n')
    fid.write('# instead of 2 for complex numbers. When your data consists of \r\n')
    fid.write('# multiple components, all components for a cell or node are written \r\n')
    fid.write('# sequentially to the file before going to the next cell or node. \r\n')
    fid.write('# DATA_COMPONENTS: 1 \r\n')

    fid.close()
    arr = array(formatstring, ndarray.flatten('F'))
    f = open(os.path.join(outDir, BOVname + '.dat'), 'wb')
    arr.tofile(f)
    f.close()