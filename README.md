[glycopeptide-analyzer-v12-en_README.md](https://github.com/user-attachments/files/33844445/glycopeptide-analyzer-v12-en_README.md)
# 🔬 GlycoPeptide Designer & Analyzer v12.0

A comprehensive, research-grade bioinformatics platform for protein glycosylation site prediction, in-silico enzymatic digestion, theoretical glycopeptide mass/isotopic envelope calculation, and targeted mass spectrometry export (optimized for **Skyline** and **Thermo Orbitrap PRM**).

---

## 🌟 Key Features

1. **Automated UniProt Integration**
   - Fetch protein sequences automatically via UniProt Accession ID (e.g., `P01857`).
   - Automatic signal peptide cleavage to analyze mature chains directly.

2. **Advanced In-Silico Digestion**
   - Supports standard and combined proteases: Trypsin, Chymotrypsin, Lys-C, Glu-C, Asp-N, Pepsin, Pronase (non-specific cleavage), and enzyme combinations (`Trypsin/Lys-C`, `Trypsin/Glu-C`).
   - Configurable missed cleavages and peptide length filters (Min/Max AA length).

3. **Post-Translational Modifications (PTMs)**
   - **Static Modifications:** Carbamidomethylation of Cysteine (`+57.02 Da`).
   - **Variable Modifications:** Oxidation of Methionine (`+15.99 Da`).
   - Combinatorial variant generation for accurate mass matching.

4. **Dynamic N-Glycan Library Generator**
   - Builds N-glycan compositions dynamically based on Core rules (Minimum 3 Hex, 2 HexNAc) with fully customizable user-defined ranges (`Hex`, `HexNAc`, `Fuc`, `Neu5Ac`).

5. **Isotopic Envelope Simulation**
   - Calculates exact atomic carbon counts (including PTM and glycan carbons) and simulates natural $^{13}	ext{C}$ isotopic distributions (M to M+5) using binomial approximation for charge states $z = +2, +3, +4$.

6. **Targeted MS Export (Skyline & Orbitrap PRM)**
   - **Skyline CSV Export:** Formatted as Small Molecules with exact Chemical Formulas to bypass sequence-to-mass mismatch errors.
   - **Orbitrap PRM Inclusion List (CSV):** Ready-to-use inclusion list for Thermo Xcalibur method editors.
   - **Detailed Excel Report:** Full analytical overview including neutral masses, motif positions, and m/z values.

---

## 🛠️ Technology Stack

- **Backend:** Python 3.11, FastAPI, Uvicorn
- **Bioinformatics & Math:** Pyteomics, Pandas, NumPy, SQLAlchemy
- **Frontend:** HTML5, Bootstrap 5, JavaScript (Vanilla), Chart.js, SheetJS (XLSX)
- **Containerization:** Docker & Docker Compose

---

## 🚀 Quick Start (Docker - Recommended)

Make sure you have [Docker Desktop](https://www.docker.com/) installed and running on your system.

1. Open a terminal in the project root directory (`glycopeptide-analyzer`).
2. Run the following command:
   ```bash
   docker compose up -d --build
   ```
3. Open your web browser and navigate to:
   **`http://localhost:8080`**

---

## 💻 Manual Setup (Anaconda / Python)

If you prefer running it locally without Docker:

1. Create and activate a clean Python environment:
   ```bash
   conda create -n glyco_env python=3.11 -y
   conda activate glyco_env
   ```
2. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Start the FastAPI backend server:
   ```bash
   uvicorn app.main:app --reload
   ```
4. Open `index.html` directly in your web browser.

---

## 📊 Skyline Import Instructions

When importing the exported **Skyline (CSV)** transition list into Skyline:
1. Go to **Edit > Insert > Transition List**.
2. **Crucial:** Select the **"Molecules"** radio button at the bottom-left of the import window. 
3. Skyline will use the embedded **Chemical Formula** column to calculate exact precursor masses effortlessly.

---

## 📄 License
This project is open-source and developed for academic and clinical research purposes.
