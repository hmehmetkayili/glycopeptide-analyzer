import re
from pyteomics import parser, mass
from typing import List, Dict, Any

class DigestionService:
    def __init__(self):
        self.rules = parser.expasy_rules.copy()
        self.rules['pepsin'] = r'(?<=[FLWAEQ])(?!P)'
        self.rules['glu-c'] = r'(?<=[E])(?!P)'
        self.rules['lys-c'] = r'(?<=[K])(?!P)'
        self.rules['arg-c'] = r'(?<=[R])(?!P)'
        self.rules['asp-n'] = r'(?=[D])'
        self.rules['trypsin_lys-c'] = r'(?<=[KR])(?!P)'
        self.rules['trypsin_glu-c'] = r'(?<=[KRE])(?!P)'

        self.MASS_CARBAMIDOMETHYL = 57.021464
        self.MASS_OXIDATION = 15.994915

    def _apply_ptms(self, pep: str, base_mass: float, static_mods: List[str], variable_mods: List[str]) -> List[Dict[str, Any]]:
        variants = []
        c_count = pep.count('C')
        m_count = pep.count('M')
        
        static_mass_add = 0.0
        mod_texts = []
        c_mod_count = 0
        
        if 'carbamidomethyl_c' in static_mods and c_count > 0:
            static_mass_add = c_count * self.MASS_CARBAMIDOMETHYL
            mod_texts.append(f"{c_count}xCarbamidomethyl(C)")
            c_mod_count = c_count
            
        fixed_mass = base_mass + static_mass_add
        
        variants.append({
            "mass": fixed_mass,
            "mod_string": " + ".join(mod_texts) if mod_texts else "None",
            "c_mod_count": c_mod_count,
            "m_mod_count": 0
        })
        
        if 'oxidation_m' in variable_mods and m_count > 0:
            for i in range(1, m_count + 1):
                v_texts = mod_texts.copy()
                v_texts.append(f"{i}xOxidation(M)")
                variants.append({
                    "mass": fixed_mass + (i * self.MASS_OXIDATION),
                    "mod_string": " + ".join(v_texts),
                    "c_mod_count": c_mod_count,
                    "m_mod_count": i
                })
                
        return variants

    def in_silico_digest(self, sequence: str, enzyme: str = 'trypsin', 
                         max_missed_cleavages: int = 2, min_len: int = 4, max_len: int = 30,
                         static_mods: List[str] = None, variable_mods: List[str] = None) -> List[Dict[str, Any]]:
                         
        if static_mods is None: static_mods = []
        if variable_mods is None: variable_mods = []
                         
        enzyme_lower = enzyme.lower().replace("+", "_").replace(" ", "")
        results = []
        
        if enzyme_lower == 'pronase':
            for i in range(len(sequence)):
                for j in range(i + min_len, min(i + max_len + 1, len(sequence) + 1)):
                    pep = sequence[i:j]
                    pep_mass = mass.calculate_mass(sequence=pep, ion_type='M', charge=0)
                    
                    variants = self._apply_ptms(pep, pep_mass, static_mods, variable_mods)
                    for var in variants:
                        results.append({
                            "sequence": pep, "length": len(pep), 
                            "start_pos": i + 1, "end_pos": j, 
                            "mass_da": round(var["mass"], 5), "missed_cleavages": 0,
                            "peptide_mods": var["mod_string"], "c_mod_count": var["c_mod_count"],
                            "m_mod_count": var["m_mod_count"]
                        })
            return results
            
        if enzyme_lower not in self.rules:
            raise ValueError(f"Unsupported enzyme: {enzyme}")
            
        rule = self.rules[enzyme_lower]
        unique_peptides = parser.cleave(sequence, rule, missed_cleavages=max_missed_cleavages)
        
        for pep in unique_peptides:
            if not (min_len <= len(pep) <= max_len):
                continue
            actual_missed = parser.num_sites(pep, rule)
            pep_mass = mass.calculate_mass(sequence=pep, ion_type='M', charge=0)
            
            variants = self._apply_ptms(pep, pep_mass, static_mods, variable_mods)
            starts = [m.start() + 1 for m in re.finditer(f'(?={re.escape(pep)})', sequence)]
            
            for start in starts:
                end = start + len(pep) - 1
                for var in variants:
                    results.append({
                        "sequence": pep, "length": len(pep), 
                        "start_pos": start, "end_pos": end, 
                        "mass_da": round(var["mass"], 5), "missed_cleavages": actual_missed,
                        "peptide_mods": var["mod_string"], "c_mod_count": var["c_mod_count"],
                        "m_mod_count": var["m_mod_count"]
                    })
        
        results.sort(key=lambda x: (x["start_pos"], x["length"]))
        return results
