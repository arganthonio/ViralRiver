@echo off
echo Starting ViralRiver Graphical Interface... Please wait a few seconds.
:: Run Docker with 2GB of shared memory (--shm-size=2g) to accelerate Kraken2
start /b docker run --rm -p 8501:8501 --shm-size=2g -v "%cd%":/workspace ajriverar/viralriver:latest
timeout /t 5 /nobreak > NUL
start http://localhost:8501