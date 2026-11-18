#!/bin/bash
# Çift tıklayınca: gerekirse ortamı kurar, sunucuyu başlatır ve tarayıcıyı açar.
cd "$(dirname "$0")" || exit 1
if [ ! -x venv/bin/python ]; then
  echo "İlk kurulum yapılıyor..."
  python3 -m venv venv && venv/bin/python -m pip install -q --upgrade pip && venv/bin/python -m pip install -q -r requirements.txt || exit 1
fi
PORT=${PORT:-5050}
(sleep 2 && open "http://127.0.0.1:$PORT") &
PORT=$PORT venv/bin/python app.py
