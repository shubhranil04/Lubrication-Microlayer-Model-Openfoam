#### import the simple module from the paraview
from paraview.simple import *
import os

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get current working directory
cwd = os.getcwd()

# create a new 'OpenFOAMReader'
foamReader = OpenFOAMReader(
    FileName=os.path.join(cwd, "system", "controlDict")
)
foamReader.MeshRegions = ['/fluid/internalMesh']
foamReader.CellArrays = ['interfaceCentre.liquid']

# get animation scene
animationScene = GetAnimationScene()

# get the time-keeper
timeKeeper = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene.UpdateAnimationUsingDataTimeSteps()

timeSteps = timeKeeper.TimestepValues

# loop through time steps
for t in timeSteps:
    animationScene.AnimationTime = t
    UpdatePipeline(time=t)
    
    # create output filename
    filename = os.path.join(
        cwd, "postProcessing", "interface", f"interface_t_{t:0.4f}.csv"
    )
    
    SaveData(
        filename,
        proxy=foamReader,
        ChooseArraysToWrite=1,
        CellDataArrays=['interfaceCentre.liquid'],
        FieldAssociation='Cell Data',
        Precision=6,
        AddTime=1,
        UseScientificNotation=1
    )
