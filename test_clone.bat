@echo off
setlocal

rem Clone a committed SciFort repository into a temporary directory and test it.
rem Usage: test_clone.bat [repository-url-or-path]

set "source_repo=%~1"
if not defined source_repo set "source_repo=%~dp0."

set "clone_dir=%TEMP%\scifort-clone-test-%RANDOM%-%RANDOM%"

where git >nul 2>nul
if errorlevel 1 goto missing_git

where fpm >nul 2>nul
if errorlevel 1 goto missing_fpm

echo Cloning "%source_repo%" into "%clone_dir%"...
git clone --quiet "%source_repo%" "%clone_dir%"
if errorlevel 1 goto clone_failed

pushd "%clone_dir%"
if errorlevel 1 goto directory_failed

echo Building cloned repository...
fpm build
if errorlevel 1 goto build_failed

echo Testing cloned repository...
fpm test
if errorlevel 1 goto test_failed

popd
rmdir /s /q "%clone_dir%"
echo Clone build and tests passed.
exit /b 0

:build_failed
echo ERROR: fpm build failed. 1>&2
popd
goto test_cleanup

:test_failed
echo ERROR: fpm test failed. 1>&2
popd
goto test_cleanup

:directory_failed
echo ERROR: Could not enter cloned repository "%clone_dir%". 1>&2
goto test_cleanup

:clone_failed
echo ERROR: git clone failed. 1>&2
goto test_cleanup

:missing_git
echo ERROR: git was not found on PATH. 1>&2
exit /b 1

:missing_fpm
echo ERROR: fpm was not found on PATH. 1>&2
exit /b 1

:test_cleanup
if exist "%clone_dir%" rmdir /s /q "%clone_dir%"
exit /b 1
