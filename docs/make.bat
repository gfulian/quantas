@echo off
setlocal

rem Build the Quantas Sphinx documentation on Windows.
rem Usage from the repository root:
rem     docs\make.bat html
rem If no target is provided, HTML is built by default.

pushd "%~dp0"

if "%PYTHON%"=="" set "PYTHON=python"
if "%~1"=="" (
    set "TARGET=html"
) else (
    set "TARGET=%~1"
)

%PYTHON% -c "import sphinx" >nul 2>nul
if errorlevel 1 (
    echo.
    echo Sphinx is not available in the selected Python environment.
    echo Install the documentation dependencies first:
    echo     %PYTHON% -m pip install -e "..[docs]"
    echo.
    set "EXITCODE=1"
    goto :finish
)

if /I "%TARGET%"=="assets" (
    %PYTHON% tools\generate_elasticity_assets.py
    set "EXITCODE=%ERRORLEVEL%"
    goto :finish
)

%PYTHON% -m sphinx -M "%TARGET%" source _build -W --keep-going %SPHINXOPTS%
set "EXITCODE=%ERRORLEVEL%"

if "%EXITCODE%"=="0" if /I "%TARGET%"=="html" (
    echo.
    echo Documentation written to:
    echo     %CD%\_build\html\index.html
)

:finish
popd
endlocal & exit /b %EXITCODE%
