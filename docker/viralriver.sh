#!/bin/bash
echo "Starting ViralRiver Graphical Interface... Please wait a few seconds."

docker run -d --rm \
  -p 8501:8501 \
  --shm-size=2g \
  -v "$(pwd)":/workspace \
  --pull missing \
  ajriverar/viralriver:latest

sleep 5

if command -v xdg-open > /dev/null; then
    xdg-open http://localhost:8501
elif command -v open > /dev/null; then
    open http://localhost:8501
else
    echo "Open http://localhost:8501 in your browser."
fi
