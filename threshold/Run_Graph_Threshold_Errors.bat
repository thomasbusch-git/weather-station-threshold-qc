@echo off
echo ================================================
echo  Hawaii Mesonet Threshold Graph Generator v1.0
echo ================================================
echo.
echo Generating threshold violation plots...
echo.
python "%~dp0Graph_Threshold_Errors.py" -i "%~dp0..\data\stationdata" -o "%~dp0output\graphs"
echo.
if %ERRORLEVEL% EQU 0 (
    echo Done! Opening output folder...
    explorer "%~dp0output\graphs"
) else (
    echo ERROR: Script failed. Check the output above.
    pause
)
