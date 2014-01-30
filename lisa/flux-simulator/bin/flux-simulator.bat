@echo off
setlocal ENABLEDELAYEDEXPANSION
::CONFIGURATION
IF "%FLUX_MEM%" == "" (
    if "%PROCESSOR_ARCHITECTURE%" == "amd64" (
        ::64 bit machine
        set FLUX_MEM=3G
    ) else (
        if "%PROCESSOR_ARCHITEW6432%" == "amd64" (
            ::32 bit app on 64 bit machine
            set FLUX_MEM=3G
        ) else (
            ::32 bit machine
            set FLUX_MEM=1400M
        )
    )
)

::add java_home to path if its set
IF defined JAVA_HOME set PATH=%JAVA_HOME%\bin;%Path%

:: test for a java installation in the path
for %%X in (java.exe) do (set FOUND=%%~$PATH:X)
if not defined FOUND (
   echo No Java installation found, please install Java ^>= 1.6 and make sure it is in your PATH.
   exit /B 1
)

set parent=%~dp0%..\lib
FOR /R "%parent%" %%G IN (*.jar) DO set CLASSPATH=!CLASSPATH!;"%%G"
java -Xmx%FLUX_MEM% %JAVA_OPTS% -cp %CLASSPATH% -Dflux.tool=simulator -Dflux.app="barna.simulator" barna.commons.launcher.Flux %*