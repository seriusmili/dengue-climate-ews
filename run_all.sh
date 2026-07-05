#!/bin/bash
# ============================================================
#  Jalankan SEMUA analisis dengan satu perintah.
#  Cara pakai (di Terminal):   bash run_all.sh
# ============================================================
set -e
cd "$(dirname "$0")"   # pindah ke folder tempat file ini berada

echo ""
echo ">>> Langkah 1/3: memasang paket Python yang dibutuhkan (sekali saja)..."
python3 -m pip install -r requirements.txt || python3 -m pip install --user -r requirements.txt

echo ""
echo ">>> Langkah 2/3: menyiapkan data (memakai iklim demo dulu)..."
cd scripts
python3 00_make_synthetic_climate.py
python3 02_build_dataset.py

echo ""
echo ">>> Langkah 3/3: menjalankan model & membuat grafik..."
python3 03_dlnm_model.py
python3 04_multiprovince.py
python3 05_forecast_eval.py

echo ""
echo "============================================================"
echo "  SELESAI. Buka folder 'figures' untuk melihat grafiknya."
echo "  Ringkasan angka ada di folder 'writeup'."
echo "============================================================"
