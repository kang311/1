# Microsoft Developer Studio Project File - Name="RMC_POT" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=RMC_POT - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "RMC_POT.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "RMC_POT.mak" CFG="RMC_POT - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "RMC_POT - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "RMC_POT - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "RMC_POT - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /Od /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /D "_MULTI" /FR /YX /FD /Zm200 /c
# ADD BASE RSC /l 0x40e /d "NDEBUG"
# ADD RSC /l 0x40e /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib pthreadVC2.lib pthreadVSE2.lib /nologo /subsystem:console /machine:I386 /out:"Release/RMC_POT_multi_noper.exe"
# SUBTRACT LINK32 /pdb:none

!ELSEIF  "$(CFG)" == "RMC_POT - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "_MULTI" /FR /YX /FD /GZ /Zm200 /c
# ADD BASE RSC /l 0x40e /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib pthreadVC2.lib pthreadVSE2.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "RMC_POT - Win32 Release"
# Name "RMC_POT - Win32 Debug"
# Begin Group "Source"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=.\altern.cpp
# End Source File
# Begin Source File

SOURCE=.\AvCoordConst.cpp
# End Source File
# Begin Source File

SOURCE=.\CalcData.cpp
# End Source File
# Begin Source File

SOURCE=.\CalcPart.cpp
# End Source File
# Begin Source File

SOURCE=.\ChiSquared.cpp
# End Source File
# Begin Source File

SOURCE=.\CoordNumbConst.cpp
# End Source File
# Begin Source File

SOURCE=.\CosDistrConst.cpp
# End Source File
# Begin Source File

SOURCE=.\DataMat.cpp
# End Source File
# Begin Source File

SOURCE=.\ExptsData.cpp
# End Source File
# Begin Source File

SOURCE=.\FNC_POT.cpp
# End Source File
# Begin Source File

SOURCE=.\History.cpp
# End Source File
# Begin Source File

SOURCE=.\HistoSet.cpp
# End Source File
# Begin Source File

SOURCE=.\makemove.cpp
# End Source File
# Begin Source File

SOURCE=.\makemovecus.cpp
# End Source File
# Begin Source File

SOURCE=.\Move.cpp
# End Source File
# Begin Source File

SOURCE=.\NeighbourList.cpp
# End Source File
# Begin Source File

SOURCE=.\PPCFSet.cpp
# End Source File
# Begin Source File

SOURCE=.\RMC_POT.cpp
# End Source File
# Begin Source File

SOURCE=.\RunParams.cpp
# End Source File
# Begin Source File

SOURCE=.\SimpleCfg.cpp
# End Source File
# Begin Source File

SOURCE=.\Threads.cpp
# End Source File
# Begin Source File

SOURCE=.\Topology.cpp
# End Source File
# Begin Source File

SOURCE=.\utilities.cpp
# End Source File
# End Group
# Begin Group "Header"

# PROP Default_Filter "*.h"
# Begin Source File

SOURCE=.\altern.h
# End Source File
# Begin Source File

SOURCE=.\classes1.h
# End Source File
# Begin Source File

SOURCE=.\classes2.h
# End Source File
# Begin Source File

SOURCE=.\global.h
# End Source File
# Begin Source File

SOURCE=.\Move.h
# End Source File
# Begin Source File

SOURCE=.\NeighbourList.h
# End Source File
# Begin Source File

SOURCE=.\Threads.h
# End Source File
# Begin Source File

SOURCE=.\Topology.h
# End Source File
# Begin Source File

SOURCE=.\units.h
# End Source File
# Begin Source File

SOURCE=.\utilities.h
# End Source File
# End Group
# End Target
# End Project
