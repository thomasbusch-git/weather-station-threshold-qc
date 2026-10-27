@echo off
echo ================================================
echo  Hawaii Mesonet Threshold Error Detector v1.0
echo ================================================
echo.
echo Running analysis...
echo.
python "%~dp0Threshold_Errors.py" -i "%~dp0..\data\stationdata" -o "%~dp0output"
echo.
if %ERRORLEVEL% EQU 0 (
    echo Analysis complete! Opening report...
    start "" "%~dp0output\threshold_output.html"
) else (
    echo ERROR: Analysis failed. Check the output above.
    pause
)
