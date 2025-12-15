#!/bin/bash
set -e

echo "Installing Python dependencies..."
pip install -r requirements.txt

echo "Downloading dataset from GitHub Release..."
# Download the GSE161340.zip from the release
wget https://github.com/uniyeh/vis-final-scRNA/releases/download/dataset/GSE161340.zip -O GSE161340.zip

echo "Extracting dataset..."
unzip -q GSE161340.zip

echo "Cleaning up zip file..."
rm GSE161340.zip

echo "Build complete!"
