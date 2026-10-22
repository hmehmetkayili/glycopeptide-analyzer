import re
from typing import Dict, List, Any
from pyteomics import mass

class MassCalculator:
    MONOSACCHARIDE_MASSES = {
        "Hex": 162.052824, "HexNAc": 203.079373, "Fuc": 146.057909, "Neu5Ac": 291.095417
    }
    PROTON_MASS = 1.007276
    C13_PROBABILITY = 0.0107

    @classmethod
    def generate_n_glycan_library(cls, ranges: Dict[str, Dict[str, int]]) -> List[str]:
        glycans = []
        for h in range(ranges['Hex']['min'], ranges['Hex']['max'] + 1):
            for n in range(ranges['HexNAc']['min'], ranges['HexNAc']['max'] + 1):
                for f in range(ranges['Fuc']['min'], ranges['Fuc']['max'] + 1):
                    for s in range(ranges['Neu5Ac']['min'], ranges['Neu5Ac']['max'] + 1):
                        glycans.append(f"H{h}N{n}F{f}S{s}")
        return glycans

    @classmethod
    def parse_glycan_composition(cls, comp_str: str) -> Dict[str, int]:
        comp = {"Hex": 0, "HexNAc": 0, "Fuc": 0, "Neu5Ac": 0}
        for match in re.finditer(r'([HNFS])(\d+)', comp_str.upper()):
            block, count = match.groups()
            if block == 'H': comp["Hex"] = int(count)
            elif block == 'N': comp["HexNAc"] = int(count)
            elif block == 'F': comp["Fuc"] = int(count)
            elif block == 'S': comp["Neu5Ac"] = int(count)
        return comp

    @classmethod
    def calculate_mz(cls, neutral_mass: float, charge: int) -> float:
        return round((neutral_mass + (charge * cls.PROTON_MASS)) / charge, 5)

    @classmethod
    def get_chemical_formula(cls, sequence: str, glycan_comp: Dict[str, int], c_mod_count: int, m_mod_count: int) -> str:
        try:
            comp = mass.Composition(sequence=sequence)
        except Exception:
            return "Unknown"

        if c_mod_count > 0:
            comp += mass.Composition({'C': 2, 'H': 3, 'N': 1, 'O': 1}) * c_mod_count

        if m_mod_count > 0:
            comp += mass.Composition({'O': 1}) * m_mod_count

        if glycan_comp.get('Hex', 0) > 0:
            comp += mass.Composition({'C': 6, 'H': 10, 'O': 5}) * glycan_comp['Hex']
        if glycan_comp.get('HexNAc', 0) > 0:
            comp += mass.Composition({'C': 8, 'H': 13, 'N': 1, 'O': 5}) * glycan_comp['HexNAc']
        if glycan_comp.get('Fuc', 0) > 0:
            comp += mass.Composition({'C': 6, 'H': 10, 'O': 4}) * glycan_comp['Fuc']
        if glycan_comp.get('Neu5Ac', 0) > 0:
            comp += mass.Composition({'C': 11, 'H': 17, 'N': 1, 'O': 8}) * glycan_comp['Neu5Ac']

        formula_str = ""
        for element in ['C', 'H', 'N', 'O', 'P', 'S']:
            count = comp.get(element, 0)
            if count > 0:
                formula_str += f"{element}{int(count)}"
        return formula_str

    @classmethod
    def count_carbons(cls, sequence: str, glycan_comp: Dict[str, int], c_mod_count: int = 0) -> int:
        try:
            pep_comp = mass.Composition(sequence=sequence)
            pep_carbons = pep_comp.get('C', 0)
        except:
            pep_carbons = len(sequence) * 5
            
        pep_carbons += (c_mod_count * 2)
        glycan_carbons = (glycan_comp.get("Hex", 0) * 6 + 
                          glycan_comp.get("HexNAc", 0) * 8 + 
                          glycan_comp.get("Fuc", 0) * 6 + 
                          glycan_comp.get("Neu5Ac", 0) * 11)
        return pep_carbons + glycan_carbons

    @classmethod
    def get_isotopic_envelope(cls, mono_mz: float, charge: int, carbon_count: int) -> List[Dict[str, float]]:
        isotopes = []
        p = cls.C13_PROBABILITY
        q = 1 - p
        n = carbon_count
        
        p0 = q ** n
        isotopes.append({"mz": mono_mz, "intensity": 100.0, "label": "M"})
        if p0 == 0: p0 = 1e-10 
        
        for k in range(1, 6):
            prev_p = (q ** n) if k == 1 else p_k
            coeff = (n - k + 1) / k
            p_k = prev_p * coeff * (p / q) if k > 1 else n * (q**(n-1)) * p
            
            rel_int = (p_k / p0) * 100.0
            isotopes.append({"mz": mono_mz + (k * 1.00335 / charge), "intensity": round(rel_int, 2), "label": f"M+{k}"})

        return isotopes

    @classmethod
    def generate_theoretical_glycopeptides(cls, mapped_peptides, glycan_compositions):
        spectra = []
        for pep in mapped_peptides:
            for glyco_str in glycan_compositions:
                comp_dict = cls.parse_glycan_composition(glyco_str)
                glycan_mass = sum([cls.MONOSACCHARIDE_MASSES[k] * v for k, v in comp_dict.items() if k in cls.MONOSACCHARIDE_MASSES])
                
                glycopep_mass = pep["mass_da"] + glycan_mass
                carbons = cls.count_carbons(pep["sequence"], comp_dict, pep.get("c_mod_count", 0))
                target_site = pep["glycosylation_sites"][0]
                
                formula = cls.get_chemical_formula(pep["sequence"], comp_dict, pep.get("c_mod_count", 0), pep.get("m_mod_count", 0))
                
                result_item = {
                    "sequence": pep["sequence"],
                    "length": pep["length"],
                    "start_pos": pep["start_pos"],
                    "end_pos": pep["end_pos"],
                    "motif": target_site["motif"],
                    "peptide_mods": pep["peptide_mods"],
                    "glycan": glyco_str,
                    "neutral_mass": round(glycopep_mass, 5),
                    "formula": formula,
                    "isotopes": {}
                }
                
                for z in [2, 3, 4]: 
                    mz_val = cls.calculate_mz(glycopep_mass, z)
                    result_item[f"mz_{z}"] = mz_val
                    result_item["isotopes"][f"z{z}"] = cls.get_isotopic_envelope(mz_val, z, carbons)
                
                spectra.append(result_item)
        return spectra
