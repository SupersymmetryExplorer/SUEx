@echo off
setlocal

:: Compiler and flags
set CXX=g++
:: set CXXFLAGS=-O2 -Wall -std=c++17

:: Target executable
set TARGET=susy.exe

:: Source and object files
set SOURCES=susy.c functions.cpp variables.cpp readwrite.cpp electroweak.cpp loop.cpp complex.cpp radiative.cpp initcond.cpp numericx.cpp higgs.cpp
set OBJECTS=susy.o functions.o variables.o readwrite.o electroweak.o loop.o complex.o radiative.o initcond.o numericx.o higgs.o

echo.
echo Compiling source files...
for %%F in (%SOURCES%) do (
    echo Compiling %%F...
    %CXX% %CXXFLAGS% -c %%F
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
