"""
provinces.py
------------
Shared configuration: the provinces analysed and a representative climate
coordinate for each (used by both the real Open-Meteo fetch and the synthetic
demo generator so the two stay consistent).
"""

# province name (matching OpenDengue after Yogyakarta merge) -> (lat, lon)
PROVINCES = {
    "JAWA TIMUR":    (-7.50, 112.50),   # Surabaya / Malang corridor
    "JAWA BARAT":    (-6.90, 107.60),   # Bandung
    "JAWA TENGAH":   (-7.00, 110.40),   # Semarang
    "DKI JAKARTA":   (-6.20, 106.85),   # Jakarta
    "BANTEN":        (-6.12, 106.15),   # Serang
    "BALI":          (-8.65, 115.20),   # Denpasar
    "DI YOGYAKARTA": (-7.80, 110.37),   # Yogyakarta
}

PRIMARY = "JAWA TIMUR"   # province used for the single-province deep-dive

START = "2004-01-01"
END = "2024-07-31"
