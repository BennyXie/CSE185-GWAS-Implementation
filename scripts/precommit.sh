#!/bin/bash
# Get current version from __version__.py
CURRENT_VERSION=$(awk -F'=' '/__version__/{print $2}' gwas_tools/__version__.py | tr -d '"' | tr -d ' ')

# Increment the version number
IFS='.' read -ra VERSION_PARTS <<< "$CURRENT_VERSION"
MINOR_VERSION=${VERSION_PARTS[1]}
MINOR_VERSION=$((MINOR_VERSION + 1))
NEW_VERSION="${VERSION_PARTS[0]}.$MINOR_VERSION"

# Update __version__.py with the new version
sed -i "s/__version__=\"$CURRENT_VERSION\"/__version__=\"$NEW_VERSION\"/" gwas_tools/__version__.py
