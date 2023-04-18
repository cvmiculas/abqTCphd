# initialization
from abaqus import *
from abaqusConstants import *
import time
import math


def create_abqTCphd(nameModel, nameCol, nameSF, bC, hC, L, tC, radE, hS, tS, dS, meshSizeGlobal, meshNoCorners, meshNoThickness, meshSizeBiasMin, meshSizeBiasMax, nameMat, dataElasticity, dataPlasticity, loadFaceA_YN, loadFaceA, loadFaceB_YN, loadFaceB, loadFaceC_YN, loadFaceC, loadFaceD_YN, loadFaceD, nlGeom, maxDispVal, DOF, columnFace, maxLPFVal, maxNumIncVal, initialArcIncVal, minArcIncVal, maxArcIncVal, totalArcLengthVal, numCpusVal, numGPUsVal):
    start_timeALL = time.time()
    
    # INITIAL CALCULATION
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ## internal radius
    if radE > 0:
        radI = radE-tC;
    else:
        radI = 0;

    # GEOMETRY VERIFICATION
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    if radI <= 0:
        print('= = = E R R O R = = => "radE - tC < = 0".\n')
        return
    
    if dS+2*tS > bC-2*radE:
        print('= = = E R R O R = = => "dS+2*tS > bC-2*radE".\n')
        return
    
    ## names for each socket faces
    nameSF1   = nameSF + '-1';
    nameSF2   = nameSF + '-2';
    nameSF3   = nameSF + '-3';
    nameSF4   = nameSF + '-4';
    
    ## section name
    nameSectionMat = 'Section-1';
    
    ## steps name
    nameStepInitial = 'Initial'
    nameStepRiks    = 'Riks'
    
    ## some lengths
    L1 = L/2-hS/2;
    L2 = L/2;
    L3 = L/2+hS/2;
    
    lengths = (0.0, L1, L2, L3, L)
    lengths_L = (L1/2, L3+L1/2)
    
    ## coordinates for Datum Points (DPs) used for partition, mesh, etc.
        ## X-dir: bottom edge
    dpX_1 = (-(bC/2-tC-radI), -hC/2, L)
    dpX_2 = (-(dS/2+tS), -hC/2, L)
    dpX_4 = (0.0, -hC/2, L)
    dpX_6 = ( (dS/2+tS), -hC/2, L)
    dpX_7 = ( (bC/2-tC-radI), -hC/2, L)
        
        ## Y-dir: right edge
    dpY_1 = ( bC/2, -(hC/2-tC-radI), L)
    dpY_2 = ( bC/2, -(dS/2+tS), L)
    dpY_4 = ( bC/2, 0.0, L)
    dpY_6 = ( bC/2,  (dS/2+tS), L)
    dpY_7 = ( bC/2,  (hC/2-tC-radI), L)
        
        ## Z-dir: top face
    dpZ_1 = (0.0, hC/2, L/2-hS/2)
    dpZ_2 = (0.0, hC/2, L/2)
    dpZ_3 = (0.0, hC/2, L/2+hS/2)
    
        ## create listS with DPs of interest
    if dS > 0:
        dpX_3  = (-(dS/2)   , -hC/2, L)
        dpX_5  = ( (dS/2)   , -hC/2, L)
        
        dpY_3  = ( bC/2, -(dS/2)   , L)
        dpY_5  = ( bC/2,  (dS/2)   , L)
        dpX = (dpX_1, dpX_2, dpX_3, dpX_4, dpX_5, dpX_6, dpX_7)
        dpY = (dpY_1, dpY_2, dpY_3, dpY_4, dpY_5, dpY_6, dpY_7)
        
    else:
        dpX = (dpX_1, dpX_2, dpX_4, dpX_6, dpX_7)
        dpY = (dpY_1, dpY_2, dpY_4, dpY_6, dpY_7)
    
    dpZ = (dpZ_1, dpZ_2, dpZ_3)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    
    
    # MODEL GENERATION
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    # # C R E A T E   C O L U M N   P A R T   * * * * * * * * * * * * * * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    ## create mdb model
    mdb.Model(name=nameModel, modelType=STANDARD_EXPLICIT)
    
    ## set current model
    modelCurrent = mdb.models[nameModel]
    
    # set current sketch
    sketchC = modelCurrent.ConstrainedSketch(name='model sketch', sheetSize=200.0)
    
    ## create sketch
        # rectangle (exterior edges)
    sketchC.rectangle(point1=(-bC/2, -hC/2), point2=(bC/2, hC/2))
        
        # rectangle (interior edges)
    sketchC.rectangle(point1=(-(bC-2*tC)/2, -(hC-2*tC)/2), point2=((bC-2*tC)/2, (hC-2*tC)/2))
    
    ## create part
    partC = modelCurrent.Part(name=nameCol, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    
    # make extrude of part
    partC.BaseSolidExtrude(sketch=sketchC, depth=L)
    
    ## create rounded corners
    if radE > 0:
        ## find external edges and round them
        edgeE1 = partC.edges.findAt(( bC/2, hC/2,L/2),)
        edgeE2 = partC.edges.findAt((-bC/2,-hC/2,L/2),)
        edgeE3 = partC.edges.findAt((-bC/2, hC/2,L/2),)
        edgeE4 = partC.edges.findAt(( bC/2,-hC/2,L/2),)
        partC.Round(radius=radE, edgeList=(edgeE1, edgeE2, edgeE3, edgeE4,))
        
        ## find internal edges and round them
        edgeI1 = partC.edges.findAt(( (bC-2*tC)/2, (hC-2*tC)/2,L/2),)
        edgeI2 = partC.edges.findAt((-(bC-2*tC)/2,-(hC-2*tC)/2,L/2),)
        edgeI3 = partC.edges.findAt((-(bC-2*tC)/2, (hC-2*tC)/2,L/2),)
        edgeI4 = partC.edges.findAt(( (bC-2*tC)/2,-(hC-2*tC)/2,L/2),)
        partC.Round(radius=radI, edgeList=(edgeI1, edgeI2, edgeI3, edgeI4,))
  
    ## create column model part
    partC = mdb.models[nameModel].parts[nameCol]
    
    ## set current viewport on column model part
    session.viewports['Viewport: 1'].setValues(displayedObject=partC)
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    
    
    # # C R E A T E   D A T U M   P O I N T S   (DPs)   (column)  * * * * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    for iX in dpX:
        partC.DatumPointByCoordinate(coords=iX)
        
    for iY in dpY:
        partC.DatumPointByCoordinate(coords=iY)
        
    for iZ in dpZ:
        partC.DatumPointByCoordinate(coords=iZ)
    
    
    
    # # C R E A T E   P A R T I T I O N S
        ## using DPs  (XY-plane)
    coord_Z_line = ((bC/2-tC-radI), hC/2, (L-hS)/3.5)
    for iZ in dpZ:
        edge_Z_int = partC.edges.findAt(coord_Z_line,)
        cells2Partition = partC.cells[:]
        partC.PartitionCellByPlanePointNormal(point=iZ, normal=edge_Z_int, cells=cells2Partition)
        
        ## using DPs  (ZY-plane)
    coord_X_line = ( (dS+tS)/2, -hC/2+tC, L)
    for iX in dpX:
        edge_X_int = partC.edges.findAt(coord_X_line,)
        cells2Partition = partC.cells[:]
        partC.PartitionCellByPlanePointNormal(point=iX, normal=edge_X_int, cells=cells2Partition)
     
        ## using DPs  (ZX-plane)
    coord_Y_line = (bC/2-tC, (dS+tS)/2, L)
    for iY in dpY:
        edge_Y_int = partC.edges.findAt(coord_Y_line,)
        cells2Partition = partC.cells[:]
        partC.PartitionCellByPlanePointNormal(point=iY, normal=edge_Y_int, cells=cells2Partition)
    
        ## at 45 deg in the rounded corners
    x_sign = [1,  1, -1, -1]
    y_sign = [1, -1, -1,  1]
    
    pt1 = ( (bC/2-tC-radI),  (hC/2-tC-radI), 0,  (bC/2)        ,  (hC/2)        , L)
    pt2 = ( (bC/2-tC-radI), -(hC/2)        , 0,  (bC/2)        , -(hC/2-tC-radI), L)
    pt3 = (-(bC/2)        , -(hC/2)        , 0, -(bC/2-tC-radI), -(hC/2-tC-radI), L)
    pt4 = (-(bC/2)        ,  (hC/2-tC-radI), 0, -(bC/2-tC-radI),  (hC/2)        , L)
    pts = (pt1,pt2,pt3,pt4)
    
    for ii in range(0,4):
        signX = x_sign[ii]
        signY = y_sign[ii]
       
        pt1_C = (signX*(bC/2-tC-radI + radE*sin(math.pi/4)), signY*(hC/2-tC-radI + radE*cos(math.pi/4)), 0)
        pt2_C = (signX*(bC/2-tC-radI + radE*sin(math.pi/4)), signY*(hC/2-tC-radI + radE*cos(math.pi/4)), L)
        pt3_C = (signX*(bC/2-tC-radI + radI*sin(math.pi/4)), signY*(hC/2-tC-radI + radI*cos(math.pi/4)), L)
        
        cells_all = partC.cells.getByBoundingBox(xMin= pts[ii][0], yMin=pts[ii][1], zMin=pts[ii][2], xMax=pts[ii][3], yMax=pts[ii][4], zMax=pts[ii][5])
        
        partC.PartitionCellByPlaneThreePoints(cells=cells_all, point1=pt1_C, point2=pt2_C, point3=pt3_C)
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    
    
    # # C R E A T E   S E T   * * * * * * * * * * * * * * * * * * * * * * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        ## all column part
    cellsAll = partC.cells[:];
    setName_ALL = 'set_ALL';
    partC.Set(cells=cellsAll, name=setName_ALL)
    
        ## Thickness for seed edges by number
    list_T = []
    for ii in range(0,len(lengths)):
        myLength = lengths[ii]
        
        list_X_pos = []
        list_X_neg = []

        for iX in range(0,len(dpX)):
            line_X_pos = (dpX[iX][0],  1*(dpX[iX][1] + tC/2), myLength)
            list_X_pos.append(partC.edges.findAt((line_X_pos,) ))
            
            line_X_neg = (dpX[iX][0], -1*(dpX[iX][1] + tC/2), myLength)
            list_X_neg.append(partC.edges.findAt((line_X_neg,) ))
        
        list_Y_pos = []
        list_Y_neg = []

        for iY in range(0,len(dpY)):
            line_Y_pos = ( 1*(dpY[iY][0] - tC/2), dpY[iY][1], myLength)
            list_Y_pos.append(partC.edges.findAt((line_Y_pos,) ))

            line_Y_neg = (-1*(dpY[iY][0] - tC/2), dpY[iY][1], myLength)
            list_Y_neg.append(partC.edges.findAt((line_Y_neg,) ))

        list_iter = list_X_pos + list_X_neg + list_Y_pos + list_Y_neg
        list_T = list_T + list_iter
        
    setName_T = 'set_T'
    
    partC.Set(edges=list_T, name=setName_T)
    
        ## Length for seed edges by single-bias
        ## L0 (first part)
    edges_L0_all = partC.edges.getByBoundingBox(xMin=-bC/2, yMin=-hC/2, zMin=0, xMax=bC/2, yMax=hC/2, zMax=(L-hS)/2)
    
    edges_L0_int = []
    
    for e_L0 in edges_L0_all:
        e_size = e_L0.getSize(0)
        if e_size == (L-hS)/2:
            edges_L0_int.append(partC.edges.findAt((e_L0.pointOn[0],) ))
            
    setName_L0 = 'set_L0'
    partC.Set(edges=edges_L0_int, name=setName_L0)
    
        ## L2 (second part)
    edges_L1_all = partC.edges.getByBoundingBox(xMin=-bC/2, yMin=-hC/2, zMin=(L+hS)/2, xMax=bC/2, yMax=hC/2, zMax=L)
    
    edges_L1_int = []
    
    for e_L1 in edges_L1_all:
        e_size = e_L1.getSize(0)
        if e_size == (L-hS)/2:
            edges_L1_int.append(partC.edges.findAt((e_L1.pointOn[0],) ))
            
    setName_L1 = 'set_L1'
    partC.Set(edges=edges_L1_int, name=setName_L1)
    
    
        ## Rounded Corners for seed edges by number
    edges_C_int = []
    for kk in range(0,len(lengths)):
        myLength3 = lengths[kk]
        
        edges_C = partC.edges.getByBoundingBox(xMin=-bC/2, yMin=-hC/2, zMin=0.99*myLength3, xMax=bC/2, yMax=hC/2, zMax=1.01*myLength3)
   
        for e_C in edges_C:
            e_size = e_C.getSize(0)
            if (abs(e_size - radE*math.pi/4) < radE/1E5) or (abs(e_size - radI*math.pi/4) < radE/1E5):
                edges_C_int.append(partC.edges.findAt((e_C.pointOn[0],) ))

    setName_C = 'set_RC'
    partC.Set(edges=edges_C_int, name=setName_C)
    
    
        ## RC for boundary conditions (BCs)
    lC_1 = ( (bC/2-tC-radI + radE*sin(math.pi/4)), (hC/2-tC-radI + radE*cos(math.pi/4)))
    lC_2 = (-(bC/2-tC-radI + radE*sin(math.pi/4)), (hC/2-tC-radI + radE*cos(math.pi/4)))
    lC_3 = (-(bC/2-tC-radI + radE*sin(math.pi/4)),-(hC/2-tC-radI + radE*cos(math.pi/4)))
    lC_4 = ( (bC/2-tC-radI + radE*sin(math.pi/4)),-(hC/2-tC-radI + radE*cos(math.pi/4)))
    
    edgeBC_1 = partC.edges.getByBoundingBox(xMin=0.99*lC_1[0], yMin=0.99*lC_1[1], zMin=0, xMax=1.01*lC_1[0], yMax=1.01*lC_1[1], zMax=L)
    
    edgeBC_2 = partC.edges.getByBoundingBox(xMin=1.01*lC_2[0], yMin=0.99*lC_2[1], zMin=0, xMax=0.99*lC_2[0], yMax=1.01*lC_2[1], zMax=L)
    
    edgeBC_3 = partC.edges.getByBoundingBox(xMin=1.01*lC_3[0], yMin=1.01*lC_3[1], zMin=0, xMax=0.99*lC_3[0], yMax=0.99*lC_3[1], zMax=L)
    
    edgeBC_4 = partC.edges.getByBoundingBox(xMin=0.99*lC_4[0], yMin=1.01*lC_4[1], zMin=0, xMax=1.01*lC_4[0], yMax=0.99*lC_4[1], zMax=L)
    
    partC.Set(edges=[edgeBC_1,edgeBC_2,edgeBC_3,edgeBC_4], name='set_BC')
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    
    
    # # C R E A T E   S U R F A C E   * * * * * * * * * * * * * * * * * * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        ## surface Socket Face # 1
    sf_1_1 = (-(dS/2+tS/2), -hC/2, L/2-hS/4);
    sf_1_2 = (-(dS/2+tS/2), -hC/2, L/2+hS/4);
    sf_1_3 = (+(dS/2+tS/2), -hC/2, L/2-hS/4);
    sf_1_4 = (+(dS/2+tS/2), -hC/2, L/2+hS/4);

    sf_1 = partC.faces.findAt((sf_1_1, ), (sf_1_2,), (sf_1_3,), (sf_1_4,))
    partC.Surface(side1Faces=sf_1, name='Surf-1')

        ## surface Socket Face # 3
    sf_3_1 = (-(dS/2+tS/2), hC/2, L/2-hS/4);
    sf_3_2 = (-(dS/2+tS/2), hC/2, L/2+hS/4);
    sf_3_3 = (+(dS/2+tS/2), hC/2, L/2-hS/4);
    sf_3_4 = (+(dS/2+tS/2), hC/2, L/2+hS/4);

    sf_3 = partC.faces.findAt((sf_3_1, ), (sf_3_2,), (sf_3_3,), (sf_3_4,))
    partC.Surface(side1Faces=sf_3, name='Surf-3')
    
        # ## surface Socket Face # 2
    sf_2_1 = ( bC/2,-(dS/2+tS/2), L/2-hS/4);
    sf_2_2 = ( bC/2,-(dS/2+tS/2), L/2+hS/4);
    sf_2_3 = ( bC/2,+(dS/2+tS/2), L/2-hS/4);
    sf_2_4 = ( bC/2,+(dS/2+tS/2), L/2+hS/4);

    sf_2 = partC.faces.findAt((sf_2_1, ), (sf_2_2,), (sf_2_3,), (sf_2_4,))
    partC.Surface(side1Faces=sf_2, name='Surf-2')
    
        ## surface Socket Face # 4
    sf_4_1 = (-bC/2,-(dS/2+tS/2), L/2-hS/4);
    sf_4_2 = (-bC/2,-(dS/2+tS/2), L/2+hS/4);
    sf_4_3 = (-bC/2,+(dS/2+tS/2), L/2-hS/4);
    sf_4_4 = (-bC/2,+(dS/2+tS/2), L/2+hS/4);

    sf_4 = partC.faces.findAt((sf_4_1, ), (sf_4_2,), (sf_4_3,), (sf_4_4,))
    partC.Surface(side1Faces=sf_4, name='Surf-4')
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *




    # # C R E A T E   M A T E R I A L   * * * * * * * * * * * * * * * * * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        ## define material name
    mdb.models[nameModel].Material(name=nameMat)
    
        ## define material type
    if   (len(dataPlasticity) == 0):
        mdb.models[nameModel].materials[nameMat].Elastic(table=dataElasticity)
    elif (len(dataPlasticity) != 0):
        mdb.models[nameModel].materials[nameMat].Elastic(table=dataElasticity)
        mdb.models[nameModel].materials[nameMat].Plastic(table=dataPlasticity)
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    
    # # C R E A T E   +   A S S I G N  S E C T I O N   T O  C E L L S   * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    mdb.models[nameModel].HomogeneousSolidSection(name=nameSectionMat, material=nameMat, thickness=None)
    regMat = partC.sets[setName_ALL]
    partC.SectionAssignment(region=regMat, sectionName=nameSectionMat, offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    
    
    # # C R E A T E   M E S H   * * * * * * * * * * * * * * * * * * * * * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        ## STEP 1: assign global seeds
    partC.seedPart(size=meshSizeGlobal, minSizeFactor=0.1)

    
        ## STEP 2: assign local seeds
            ## Rounded Corners
    pickedEdgesC = partC.sets[setName_C].edges
    partC.seedEdgeByNumber(edges=pickedEdgesC, number=meshNoCorners, constraint=FINER)
    
            ## Thickness
    pickedEdgesT = partC.sets[setName_T].edges
    partC.seedEdgeByNumber(edges=pickedEdgesT, number=meshNoThickness, constraint=FINER)

            ## Length using single-biasing direction - L0
    pickedEdges1 = partC.sets[setName_L0].edges
    edges_dir_1 = []
    edges_dir_2 = []

    for ii in pickedEdges1:
        verticesVal = ii.getVertices()
        # print(verticesVal)
            
        coord_Z_vert1 = partC.vertices[verticesVal[0]].pointOn
        # print(coord_Z_vert1)
        
        coord_Z_vert2 = partC.vertices[verticesVal[1]].pointOn
        # print(coord_Z_vert2)            
                    
        if (coord_Z_vert1[0][2] < coord_Z_vert2[0][2]):            
            edges_dir_1.append(ii)            
        elif (coord_Z_vert1[0][2] > coord_Z_vert2[0][2]):
            edges_dir_2.append(ii)
            
    partC.seedEdgeByBias(biasMethod=SINGLE, end1Edges=edges_dir_2 , end2Edges=edges_dir_1, minSize=5.0, maxSize=50.0, constraint=FINER)

            ## Length using single-biasing direction - L1
    pickedEdges2 = partC.sets[setName_L1].edges
    edges_dir_1 = []
    edges_dir_2 = []

    for ii in pickedEdges2:
        verticesVal = ii.getVertices()
        # print(verticesVal)
            
        coord_Z_vert1 = partC.vertices[verticesVal[0]].pointOn
        # print(coord_Z_vert1)
        
        coord_Z_vert2 = partC.vertices[verticesVal[1]].pointOn
        # print(coord_Z_vert2)            
                    
        if (coord_Z_vert1[0][2] < coord_Z_vert2[0][2]):            
            edges_dir_1.append(ii)            
        elif (coord_Z_vert1[0][2] > coord_Z_vert2[0][2]):
            edges_dir_2.append(ii)
            
    partC.seedEdgeByBias(biasMethod=SINGLE, end1Edges=edges_dir_1 , end2Edges=edges_dir_2, minSize=meshSizeBiasMin, maxSize=meshSizeBiasMax, constraint=FINER)

    
        ## STEP 4: generate mesh
    partC.generateMesh()
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    

    # # C R E A T E    A S S E M B  L Y   * * * * * * * * * * * * * * * * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    assembly = mdb.models[nameModel].rootAssembly
    assembly.Instance(name=nameCol, part=partC, dependent=ON)
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *



    # # A D D    B C s   T O   C O L U M N    * * * * * * * * * * * * * * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        ## add BCs to Rounded Corners
    regCol = assembly.instances[nameCol].sets['set_BC']
    mdb.models[nameModel].DisplacementBC(name='BC_rc', createStepName=nameStepInitial, region=regCol, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *



    # # C R E A T E   S T A T I C   R I K S   S T E  P  * * * * * * * * * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    mdb.models[nameModel].StaticRiksStep(name=nameStepRiks, previous=nameStepInitial)
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *



    # # C R E A T E   F I E L D   O U T P U T   * * * * * * * * * * * * * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    mdb.models[nameModel].FieldOutputRequest(name='F-Output-1', createStepName=nameStepRiks, variables=PRESELECT)
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    
    
    # # C R E A T E   H I S T O R Y   O U T P U T   * * * * * * * * * * * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    mdb.models[nameModel].HistoryOutputRequest(name='H-Output-1', createStepName=nameStepRiks, variables=PRESELECT)
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


    # # C R E A T E   S O C K E T   F A C E S   * * * * * * * * * * * * * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    existanceSF = [loadFaceA_YN, loadFaceB_YN, loadFaceC_YN, loadFaceD_YN]
    
    if not any(existanceSF):
        print('= = = E R R O R = = => No load defined. Load required on at least one face.\n')
        return

    else:
    
        if loadFaceA_YN:
            
            if loadFaceA == 0:
                print('= = = E R R O R = = => Load value on face A is zero. Define a non-zero value.\n')
                return
                
            
            else:
                ## CREATE SOCKET FACE 1/A
                    ## create sketch and draw
                sketchSF1 = modelCurrent.ConstrainedSketch(name='model sf 1', sheetSize=200.0)
                sketchSF1.Line(point1=( -(dS/2+tS), -hC/2), point2=(-(dS/2),    -hC/2) )
                sketchSF1.Line(point1=(  (dS/2),    -hC/2), point2=( (dS/2+tS), -hC/2) )
                partSF1 = modelCurrent.Part(name=nameSF1, dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
                
                partSF1 = mdb.models[nameModel].parts[nameSF1]
                    
                    ## create Datum Point
                partSF1.BaseShellExtrude(sketch=sketchSF1, depth=hS)
                partSF1.DatumPointByCoordinate(coords=(0.0, -hC/2, hS/2))
                    
                    ## assign a Referrence Point to the Datum Point
                sf1_DP_key = partSF1.datums.keys()
                partSF1.ReferencePoint(point=partSF1.datums[sf1_DP_key[0]])
                
                    ## create set-1 for RP (to be used in L/BC app)
                r_key = partSF1.referencePoints.keys()
                partSF1.Set(referencePoints=(partSF1.referencePoints[r_key[0]], ), name='set-1')
                
                    ## create surface for contact
                list_1 = []
                list_2 = []
                for ii in range(0,len(partSF1.faces)):
                    if   (partSF1.faces[ii].getNormal() == (0.0,1.0,0.0)):
                        centroid = partSF1.faces[ii].getCentroid()
                        list_1.append(partSF1.faces.findAt(centroid,))
                    elif (partSF1.faces[ii].getNormal() == (0.0,-1.0,0.0)):
                        centroid = partSF1.faces[ii].getCentroid()
                        list_2.append(partSF1.faces.findAt(centroid,))
                
                partSF1.Surface(side1Faces=list_1, side2Faces=list_2, name='Surf-1')
                    
                    ## mesh part
                partSF1.seedPart(size=5.0, deviationFactor=0.1, minSizeFactor=0.1)
                partSF1.generateMesh()
                
                    ## assembly
                    ## add, translate, assing regions, create constraint
                assembly.Instance(name=nameSF1, part=partSF1, dependent=ON)
                assembly.translate(instanceList=(nameSF1, ), vector=(0.0, 0.0, 1440.0))
                region1=assembly.instances[nameCol].surfaces['Surf-1']
                region2=assembly.instances[nameSF1].surfaces['Surf-1']
                mdb.models[nameModel].Tie(name='Constraint-1', master=region2, slave=region1, 
                positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
                
                    ## create set for BC
                setRP1 = assembly.allInstances[nameSF1].sets['set-1']
                
                    ## add BC
                mdb.models[nameModel].DisplacementBC(name='BC_sf1', createStepName=nameStepInitial, region=setRP1, u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
                
                    ## add load
                mdb.models[nameModel].ConcentratedForce(name='Load-1', createStepName=nameStepRiks, region=setRP1, cf2=loadFaceA, distributionType=UNIFORM, field='', localCsys=None)
                
                    ## create history output from RP for socket-face-1
                mdb.models[nameModel].HistoryOutputRequest(name='HO_SF1', createStepName=nameStepRiks, variables=('CF2', 'RF2', 'U2'), region=setRP1, sectionPoints=DEFAULT, rebar=EXCLUDE)



        if loadFaceB_YN:
        
            if loadFaceB == 0:
                print('= = = E R R O R = = => Load value on face B is zero. Define a non-zero value.\n')
                return
        
            else:
                ## CREATE SOCKET FACE 2/B
                    ## create sketch and draw
                sketchSF2 = modelCurrent.ConstrainedSketch(name='model sf 2', sheetSize=200.0)
                sketchSF2.Line(point1=( bC/2, -(dS/2+tS),), point2=(bC/2, -(dS/2),) )
                sketchSF2.Line(point1=( bC/2,  (dS/2),),    point2=(bC/2,  (dS/2+tS),) )
                partSF2 = modelCurrent.Part(name=nameSF2, dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
                
                partSF2 = mdb.models[nameModel].parts[nameSF2]
                    
                    ## create Datum Point
                partSF2.BaseShellExtrude(sketch=sketchSF2, depth=hS)
                partSF2.DatumPointByCoordinate(coords=(bC/2, 0.0, hS/2))
                    
                    ## assign a Referrence Point to the Datum Point
                sf2_DP_key = partSF2.datums.keys()
                partSF2.ReferencePoint(point=partSF2.datums[sf2_DP_key[0]])
                
                    ## create set-1 for RP (to be used in L/BC app)
                r_key = partSF2.referencePoints.keys()
                partSF2.Set(referencePoints=(partSF2.referencePoints[r_key[0]], ), name='set-1')
                
                    ## create surface for contact
                list_1 = []
                list_2 = []
                for ii in range(0,len(partSF2.faces)):
                    if   (partSF2.faces[ii].getNormal() == (1.0,0.0,0.0)):
                        centroid = partSF2.faces[ii].getCentroid()
                        list_1.append(partSF2.faces.findAt(centroid,))
                    elif (partSF2.faces[ii].getNormal() == (-1.0,0.0,0.0)):
                        centroid = partSF2.faces[ii].getCentroid()
                        list_2.append(partSF2.faces.findAt(centroid,))
                
                partSF2.Surface(side1Faces=list_2, side2Faces=list_1, name='Surf-2')
                    
                    ## mesh part
                partSF2.seedPart(size=5.0, deviationFactor=0.1, minSizeFactor=0.1)
                partSF2.generateMesh()
                
                    ## assembly
                    ## add, translate, assing regions, create constraint
                assembly.Instance(name=nameSF2, part=partSF2, dependent=ON)
                assembly.translate(instanceList=(nameSF2, ), vector=(0.0, 0.0, 1440.0))
                region1=assembly.instances[nameCol].surfaces['Surf-2']
                region2=assembly.instances[nameSF2].surfaces['Surf-2']
                mdb.models[nameModel].Tie(name='Constraint-2', master=region2, slave=region1, 
                positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
                
                    ## create set for BC
                setRP2 = assembly.allInstances[nameSF2].sets['set-1']
                
                    ## add BC
                mdb.models[nameModel].DisplacementBC(name='BC_sf2', createStepName=nameStepInitial, region=setRP2, u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
                
                    ## add load
                mdb.models[nameModel].ConcentratedForce(name='Load-2', createStepName=nameStepRiks, region=setRP2, cf1=loadFaceB, distributionType=UNIFORM, field='', localCsys=None)
                
                    ## create history output from RP for socket-face-2
                mdb.models[nameModel].HistoryOutputRequest(name='HO_SF2', createStepName=nameStepRiks, variables=('CF1', 'RF1', 'U1'), region=setRP2, sectionPoints=DEFAULT, rebar=EXCLUDE)


            
        if loadFaceC_YN:
        
            if loadFaceC == 0:
                print('= = = E R R O R = = => Load value on face C is zero. Define a non-zero value.\n')
                return
                
            else:
                ## CREATE SOCKET FACE 3/C
                    ## create sketch and draw
                sketchSF3 = modelCurrent.ConstrainedSketch(name='model sf 3', sheetSize=200.0)
                sketchSF3.Line(point1=( -(dS/2+tS), hC/2), point2=(-(dS/2),    hC/2) )
                sketchSF3.Line(point1=(  (dS/2),    hC/2), point2=( (dS/2+tS), hC/2) )
                partSF3 = modelCurrent.Part(name=nameSF3, dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
                
                partSF3 = mdb.models[nameModel].parts[nameSF3]
                    
                    ## create Datum Point
                partSF3.BaseShellExtrude(sketch=sketchSF3, depth=hS)
                partSF3.DatumPointByCoordinate(coords=(0.0, hC/2, hS/2))
                    
                    ## assign a Referrence Point to the Datum Point
                sf1_DP_key = partSF3.datums.keys()
                partSF3.ReferencePoint(point=partSF3.datums[sf1_DP_key[0]])
                
                    ## create set-1 for RP (to be used in L/BC app)
                r_key = partSF3.referencePoints.keys()
                partSF3.Set(referencePoints=(partSF3.referencePoints[r_key[0]], ), name='set-1')
                
                    ## create surface for contact
                list_1 = []
                list_2 = []
                for ii in range(0,len(partSF3.faces)):
                    if   (partSF3.faces[ii].getNormal() == (0.0,1.0,0.0)):
                        centroid = partSF3.faces[ii].getCentroid()
                        list_1.append(partSF3.faces.findAt(centroid,))
                    elif (partSF3.faces[ii].getNormal() == (0.0,-1.0,0.0)):
                        centroid = partSF3.faces[ii].getCentroid()
                        list_2.append(partSF3.faces.findAt(centroid,))
                
                partSF3.Surface(side1Faces=list_2, side2Faces=list_1, name='Surf-3')
                    
                    ## mesh part
                partSF3.seedPart(size=5.0, deviationFactor=0.1, minSizeFactor=0.1)
                partSF3.generateMesh()
                
                    ## assembly
                    ## add, translate, assing regions, create constraint
                assembly.Instance(name=nameSF3, part=partSF3, dependent=ON)
                assembly.translate(instanceList=(nameSF3, ), vector=(0.0, 0.0, 1440.0))
                region1=assembly.instances[nameCol].surfaces['Surf-3']
                region2=assembly.instances[nameSF3].surfaces['Surf-3']
                mdb.models[nameModel].Tie(name='Constraint-3', master=region2, slave=region1, 
                positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON) 
                
                    ## create set for BC
                setRP3 = assembly.allInstances[nameSF3].sets['set-1']
                
                    ## add BC
                mdb.models[nameModel].DisplacementBC(name='BC_sf3', createStepName=nameStepInitial, region=setRP3, u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
                
                    ## add load
                mdb.models[nameModel].ConcentratedForce(name='Load-3', createStepName=nameStepRiks, region=setRP3, cf2=loadFaceC, distributionType=UNIFORM, field='', localCsys=None)
                    
                    ## create history output from RP for socket-face-3
                mdb.models[nameModel].HistoryOutputRequest(name='HO_SF3', createStepName=nameStepRiks, variables=('CF2', 'RF2', 'U2'), region=setRP3, sectionPoints=DEFAULT, rebar=EXCLUDE)
            

            
        if loadFaceD_YN:
        
            if loadFaceD == 0:
                print('= = = E R R O R = = => Load value on face D is zero. Define a non-zero value.\n')
                return
                
            else:
                
                ## CREATE SOCKET FACE 4/D
                    ## create sketch and draw
                sketchSF4 = modelCurrent.ConstrainedSketch(name='model sf 4', sheetSize=200.0)
                sketchSF4.Line(point1=(-bC/2, -(dS/2+tS),), point2=(-bC/2, -(dS/2),) )
                sketchSF4.Line(point1=(-bC/2,  (dS/2),),    point2=(-bC/2,  (dS/2+tS),) )
                partSF4 = modelCurrent.Part(name=nameSF4, dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
                
                partSF4 = mdb.models[nameModel].parts[nameSF4]
                    
                    ## create Datum Point
                partSF4.BaseShellExtrude(sketch=sketchSF4, depth=hS)
                partSF4.DatumPointByCoordinate(coords=(-bC/2, 0.0, hS/2))
                    
                    ## assign a Referrence Point to the Datum Point
                sf4_DP_key = partSF4.datums.keys()
                partSF4.ReferencePoint(point=partSF4.datums[sf4_DP_key[0]])
                
                    ## create set-1 for RP (to be used in L/BC app)
                r_key = partSF4.referencePoints.keys()
                partSF4.Set(referencePoints=(partSF4.referencePoints[r_key[0]], ), name='set-1')
                
                    ## create surface for contact
                list_1 = []
                list_2 = []
                for ii in range(0,len(partSF4.faces)):
                    if   (partSF4.faces[ii].getNormal() == (1.0,0.0,0.0)):
                        centroid = partSF4.faces[ii].getCentroid()
                        list_1.append(partSF4.faces.findAt(centroid,))
                    elif (partSF4.faces[ii].getNormal() == (-1.0,0.0,0.0)):
                        centroid = partSF4.faces[ii].getCentroid()
                        list_2.append(partSF4.faces.findAt(centroid,))
                
                partSF4.Surface(side1Faces=list_1, side2Faces=list_2, name='Surf-4')
                    
                    ## mesh part
                partSF4.seedPart(size=5.0, deviationFactor=0.1, minSizeFactor=0.1)
                partSF4.generateMesh()
                
                    ## assembly
                    ## add, translate, assing regions, create constraint
                assembly.Instance(name=nameSF4, part=partSF4, dependent=ON)
                assembly.translate(instanceList=(nameSF4, ), vector=(0.0, 0.0, 1440.0))
                region1=assembly.instances[nameCol].surfaces['Surf-4']
                region2=assembly.instances[nameSF4].surfaces['Surf-4']
                mdb.models[nameModel].Tie(name='Constraint-4', master=region2, slave=region1, 
                positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
                
                    ## create set for BC
                setRP4 = assembly.allInstances[nameSF4].sets['set-1']
                
                    ## on socket face 4
                mdb.models[nameModel].DisplacementBC(name='BC_sf4', createStepName=nameStepInitial, region=setRP4, u1=UNSET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
                
                    ## add load
                mdb.models[nameModel].ConcentratedForce(name='Load-4', createStepName=nameStepRiks, region=setRP4, cf1=loadFaceD, distributionType=UNIFORM, field='', localCsys=None)
                
                    ## create history output from RP for socket-face-4
                mdb.models[nameModel].HistoryOutputRequest(name='HO_SF4', createStepName=nameStepRiks, variables=('CF1', 'RF1', 'U1'), region=setRP4, sectionPoints=DEFAULT, rebar=EXCLUDE)
                
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    
    
    # # U P D A T E   S T E  P   R I K S  * * * * * * * * * * * * * * * * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    ## geoetric nonlinearity
    if nlGeom:
        mdb.models[nameModel].steps[nameStepRiks].setValues(nlgeom=ON)
    else:
        mdb.models[nameModel].steps[nameStepRiks].setValues(nlgeom=OFF)
        
    ## basic 
    mdb.models[nameModel].steps[nameStepRiks].setValues(maxLPF=maxLPFVal)
    
    ## stopping criteria
    if (maxDispVal != 0) and (DOF != '0') and (columnFace != '0'):
        
        if   (columnFace == 'A') and (loadFaceA_YN):
            regionMaxDisp = mdb.models[nameModel].rootAssembly.allInstances[nameSF1].sets['set-1']
            
            mdb.models[nameModel].steps[nameStepRiks].setValues(nodeOn=ON, maximumDisplacement=maxDispVal, region=regionMaxDisp, dof=int(DOF))
            
        elif (columnFace == 'B') and (loadFaceB_YN):
            regionMaxDisp = mdb.models[nameModel].rootAssembly.allInstances[nameSF2].sets['set-1']
            
            mdb.models[nameModel].steps[nameStepRiks].setValues(nodeOn=ON, maximumDisplacement=maxDispVal, region=regionMaxDisp, dof=int(DOF))
            
        elif (columnFace == 'C') and (loadFaceC_YN):
            regionMaxDisp = mdb.models[nameModel].rootAssembly.allInstances[nameSF3].sets['set-1']
            
            mdb.models[nameModel].steps[nameStepRiks].setValues(nodeOn=ON, maximumDisplacement=maxDispVal, region=regionMaxDisp, dof=int(DOF))
            
        elif (columnFace == 'D') and (loadFaceD_YN):
            regionMaxDisp = mdb.models[nameModel].rootAssembly.allInstances[nameSF4].sets['set-1']
            
            mdb.models[nameModel].steps[nameStepRiks].setValues(nodeOn=ON, maximumDisplacement=maxDispVal, region=regionMaxDisp, dof=int(DOF))
        
    ## incrementation
    mdb.models[nameModel].steps[nameStepRiks].setValues(maxNumInc=maxNumIncVal)
    mdb.models[nameModel].steps[nameStepRiks].setValues(initialArcInc=initialArcIncVal)
    mdb.models[nameModel].steps[nameStepRiks].setValues(maxArcInc=maxArcIncVal)
    mdb.models[nameModel].steps[nameStepRiks].setValues(minArcInc=minArcIncVal)
    mdb.models[nameModel].steps[nameStepRiks].setValues(totalArcLength=totalArcLengthVal)
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    
    # # C R E A T E   J O B   * * * * * * * * * * * * * * * * * * * * * * * * * * *
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    ## CREATE JOB
    jobName = 'j_' + nameModel
    
    if numCpusVal == 1:
        mdb.Job(name=jobName, model=nameModel, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=numCpusVal, numGPUs=numGPUsVal)
    elif numCpusVal > 1 :
        mdb.Job(name=jobName, model=nameModel, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',     scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=numCpusVal, numDomains=numCpusVal, numGPUs=0)
        
    # # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    
    print("- - - - - %s seconds - - - - -\n" % (time.time() - start_timeALL))
    print("= = = = = Model created successfully! = = = = = \n")