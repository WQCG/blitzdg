# Assumes MinGW64 installation exists.
$ErrorActionPreference = 'Stop'

Invoke-RestMethod -Uri https://s3.amazonaws.com/dsteinmo-libs/include-deps.zip -OutFile include-deps.zip
Invoke-RestMethod -Uri https://s3.amazonaws.com/dsteinmo-libs/deps-win-mingw64.zip -OutFile deps-win-mingw64.zip
Expand-Archive -Path include-deps.zip -DestinationPath .\include\
Expand-Archive -Path deps-win-mingw64.zip -DestinationPath .\lib\