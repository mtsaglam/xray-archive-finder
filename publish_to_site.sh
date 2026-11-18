#!/bin/bash
# Copies the UI into the website repo, pointing it at the hosted API.
# Usage: ./publish_to_site.sh https://xray-archive-finder.onrender.com [site_repo_dir]
set -euo pipefail
API="${1:?API base URL required, e.g. https://xray-archive-finder.onrender.com}"
SITE="${2:-$HOME/Documents/mtsaglam.github.io}"
API="${API%/}/"
mkdir -p "$SITE/xray-finder"
sed "s#<!--API_BASE-->#<script>window.XRAY_API_BASE = \"$API\";</script>#" \
  "$(dirname "$0")/static/index.html" > "$SITE/xray-finder/index.html"
echo "Wrote $SITE/xray-finder/index.html (API: $API)"
