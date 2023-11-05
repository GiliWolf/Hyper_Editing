#!/usr/bin/env bash
BUILD_DIR=/home/alu/hillelr/scripts/Releases/pyInstaller/RNAEditingIndex
pyinstaller --name RNAEditingIndex \
            --onedir \
            --clean \
            --add-data /home/alu/hillelr/scripts/GGPS/Tools/EditingIndex/Configs:./Configs \
            --add-data /home/alu/hillelr/scripts/GGPS/Tools/EditingIndex/path_locator.py:. \
            --paths=/home/alu/hillelr/scripts/GGPS:/private/common/Software/PythonModules \
            --workpath=${BUILD_DIR}/build \
            --distpath=${BUILD_DIR}/dist\
            /home/alu/hillelr/scripts/GGPS/Tools/EditingIndex/A2IEditingIndex.py


#pyinstaller --name RNAEditingIndexSummary \
#            --onedir \
#            --clean \
#            --add-data /home/alu/hillelr/scripts/GGPS/Tools/EditingIndex/Configs:./Configs \
#            --add-data /home/alu/hillelr/scripts/GGPS/Tools/EditingIndex/path_locator.py:. \
#            --paths=/home/alu/hillelr/scripts/GGPS:/private/common/Software/PythonModules \
#            --workpath=${BUILD_DIR}/build \
#            --distpath=${BUILD_DIR}/dist\
#            /home/alu/hillelr/scripts/GGPS/Tools/EditingIndex/A2IEditingIndexSummary.py
