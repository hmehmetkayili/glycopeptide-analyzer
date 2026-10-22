import re
from typing import List, Dict, Any

class GlycosylationService:
    N_GLYCO_MOTIF = re.compile(r'(?=N[^P][ST])')
    
    @classmethod
    def find_n_linked_sites(cls, sequence: str, context_size: int = 5) -> List[Dict[str, Any]]:
        sites = []
        for match in cls.N_GLYCO_MOTIF.finditer(sequence):
            pos = match.start()
            motif = sequence[pos:pos+3]
            sites.append({"position": pos + 1, "glyco_type": "N-LINKED", "motif": motif})
        return sites

    @classmethod
    def map_peptides_to_sites(cls, peptides: List[Dict[str, Any]], sites: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        mapped_results = []
        for pep in peptides:
            pep_sites = []
            for site in sites:
                if pep["start_pos"] <= site["position"] <= pep["end_pos"]:
                    pep_sites.append(site)
            if pep_sites:
                pep_copy = pep.copy()
                pep_copy["glycosylation_sites"] = pep_sites
                mapped_results.append(pep_copy)
        return mapped_results
