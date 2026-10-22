import re
from typing import Dict, List, Any
from app.services.mass_calculator import MassCalculator

class FragmentationService:
    AA_MASSES = {
        'A': 71.037114, 'R': 156.101111, 'N': 114.042927, 'D': 115.026943,
        'C': 103.009185, 'E': 129.042593, 'Q': 128.058578, 'G': 57.021464,
        'H': 137.058912, 'I': 113.084064, 'L': 113.084064, 'K': 128.094963,
        'M': 131.040485, 'F': 147.068414, 'P': 97.052764, 'S': 87.032028,
        'T': 101.047679, 'W': 186.079313, 'Y': 163.063329, 'V': 99.068414
    }
    PROTON_MASS = 1.007276
    H2O_MASS = 18.010565

    # Oxonium Ions with realistic relative HCD intensities
    OXONIUM_IONS = [
        {"name": "Hex", "mz": 163.060, "intensity": 25},
        {"name": "HexNAc-2H2O", "mz": 168.065, "intensity": 30},
        {"name": "HexNAc-H2O", "mz": 186.076, "intensity": 60},
        {"name": "HexNAc", "mz": 204.086, "intensity": 100}, # Usually Base Peak
        {"name": "Neu5Ac-H2O", "mz": 274.092, "requires": "Neu5Ac", "intensity": 80},
        {"name": "Neu5Ac", "mz": 292.102, "requires": "Neu5Ac", "intensity": 90},
        {"name": "HexHexNAc", "mz": 366.139, "intensity": 45},
        {"name": "Fuc", "mz": 147.065, "requires": "Fuc", "intensity": 20},
        {"name": "HexNAcFuc", "mz": 350.145, "requires": "Fuc", "intensity": 35},
        {"name": "HexNAcHexNeu5Ac", "mz": 657.234, "requires": "Neu5Ac", "intensity": 15}
    ]

    @classmethod
    def _parse_ptms(cls, sequence: str, peptide_mods: str) -> Dict[int, float]:
        ptm_map = {}
        if "Carbamidomethyl" in peptide_mods:
            for i, aa in enumerate(sequence):
                if aa == 'C': ptm_map[i] = 57.021464
        if "Oxidation" in peptide_mods:
            for i, aa in enumerate(sequence):
                if aa == 'M': 
                    ptm_map[i] = ptm_map.get(i, 0.0) + 15.994915
                    break 
        return ptm_map

    @classmethod
    def calculate_ms2_spectrum(cls, sequence: str, glycan: str, peptide_mods: str, precursor_charge: int) -> Dict[str, Any]:
        ptm_map = cls._parse_ptms(sequence, peptide_mods)
        glyco_comp = MassCalculator.parse_glycan_composition(glycan)
        
        b_ions = []
        y_ions = []
        seq_length = len(sequence)
        
        b_mass = 0.0
        y_mass = cls.H2O_MASS

        # Realistic backbone ion intensities are generally low for intact glycopeptides (10-30%)
        for i in range(seq_length - 1):
            n_term_aa = sequence[i]
            b_mass += cls.AA_MASSES.get(n_term_aa, 0) + ptm_map.get(i, 0.0)
            
            c_term_idx = seq_length - 1 - i
            c_term_aa = sequence[c_term_idx]
            y_mass += cls.AA_MASSES.get(c_term_aa, 0) + ptm_map.get(c_term_idx, 0.0)

            for z in range(1, min(precursor_charge, 3) + 1):
                b_ions.append({"label": f"b{i+1}^{z}+", "mz": round((b_mass + z*cls.PROTON_MASS)/z, 4), "type": "b", "intensity": 15.0})
                y_ions.append({"label": f"y{i+1}^{z}+", "mz": round((y_mass + z*cls.PROTON_MASS)/z, 4), "type": "y", "intensity": 25.0})

        # Peptide Y-ions (glycan losses)
        naked_pep_mass = cls.H2O_MASS + sum(cls.AA_MASSES.get(aa, 0) for aa in sequence) + sum(ptm_map.values())
        pep_y_ions = []
        
        for z in range(1, precursor_charge + 1):
            pep_y_ions.append({"label": f"Y0^{z}+ (Peptide)", "mz": round((naked_pep_mass + z*cls.PROTON_MASS)/z, 4), "type": "pep_y", "intensity": 65.0})
            if glyco_comp.get("HexNAc", 0) >= 1:
                pep_y_ions.append({"label": f"Y1^{z}+ (Pep+HexNAc)", "mz": round((naked_pep_mass + 203.079373 + z*cls.PROTON_MASS)/z, 4), "type": "pep_y", "intensity": 85.0})
            if glyco_comp.get("HexNAc", 0) >= 2:
                pep_y_ions.append({"label": f"Y2^{z}+ (Pep+2HexNAc)", "mz": round((naked_pep_mass + 2*203.079373 + z*cls.PROTON_MASS)/z, 4), "type": "pep_y", "intensity": 50.0})
            if glyco_comp.get("HexNAc", 0) >= 1 and glyco_comp.get("Fuc", 0) >= 1:
                pep_y_ions.append({"label": f"Y1F^{z}+ (Pep+HexNAc+Fuc)", "mz": round((naked_pep_mass + 203.079373 + 146.057909 + z*cls.PROTON_MASS)/z, 4), "type": "pep_y", "intensity": 70.0})

        # Oxonium Ions
        oxoniums = []
        for oxo in cls.OXONIUM_IONS:
            if "requires" in oxo and glyco_comp.get(oxo["requires"], 0) == 0:
                continue
            oxoniums.append({"label": oxo["name"], "mz": oxo["mz"], "type": "oxo", "intensity": oxo["intensity"]})

        return {
            "oxonium_ions": oxoniums,
            "b_ions": b_ions,
            "y_ions": y_ions,
            "pep_y_ions": pep_y_ions
        }
