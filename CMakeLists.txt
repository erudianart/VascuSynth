PROJECT( VascuSynth )

FIND_PACKAGE( ITK )
IF( ITK_FOUND )
	INCLUDE( ${ITK_USE_FILE} )
ENDIF (ITK_FOUND )

ADD_EXECUTABLE(VascuSynth VascuSynth.cpp SupplyMap.cpp OxygenationMap.cpp NodeTable.cpp VascularTree.cpp TreeDrawer.cpp)

TARGET_LINK_LIBRARIES( VascuSynth ITKCommon ITKIO)