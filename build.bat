@echo off
setlocal

:: Compiler and flags
set CXX=g++
:: set CXXFLAGS=-O2 -Wall -std=c++17

:: Target executable
set TARGET=susy.exe
set OBJDIR=obj

:: Source and object files
set SOURCES=susy.c functions.cpp variables.cpp readwrite.cpp electroweak.cpp loop.cpp complex.cpp radiative.cpp initcond.cpp numericx.cpp higgs.cpp
set OBJECTS=%OBJDIR%\susy.o %OBJDIR%\functions.o %OBJDIR%\variables.o %OBJDIR%\readwrite.o %OBJDIR%\electroweak.o %OBJDIR%\loop.o %OBJDIR%\complex.o %OBJDIR%\radiative.o %OBJDIR%\initcond.o %OBJDIR%\numericx.o %OBJDIR%\higgs.o

if not exist "%OBJDIR%" mkdir "%OBJDIR%"

echo.
echo Compiling source files...
for %%F in (%SOURCES%) do (
    echo Compiling %%F...
    %CXX% %CXXFLAGS% -c %%F -o "%OBJDIR%\%%~nF.o"
    if errorlevel 1 goto :error
)

echo.
echo Linking %TARGET%...
%CXX% -o %TARGET% %OBJECTS%
if errorlevel 1 goto :error

echo.
echo Build successful!
goto :eof

:error
echo.
echo Build failed!
exit /b 1
