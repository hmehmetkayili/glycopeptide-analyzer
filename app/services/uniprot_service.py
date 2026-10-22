import httpx
import logging
from typing import Optional, Dict, Any, List

logger = logging.getLogger(__name__)

class UniProtService:
    def __init__(self):
        self.base_url = "https://rest.uniprot.org/uniprotkb/"
        self.experimental_eco_codes = ["ECO:0000269", "ECO:0000303"]

    async def fetch_protein_data(self, accession: str) -> Optional[Dict[str, Any]]:
        url = f"{self.base_url}{accession}.json"
        async with httpx.AsyncClient() as client:
            try:
                response = await client.get(url, timeout=15.0)
                if response.status_code != 200:
                    return None
                raw_data = response.json()
                return self._parse_uniprot_json(raw_data)
            except Exception as e:
                logger.error(f"UniProt Error: {e}")
        return None

    def _parse_uniprot_json(self, data: Dict[str, Any]) -> Dict[str, Any]:
        accession = data.get("primaryAccession", "")
        sequence_data = data.get("sequence", {})
        full_sequence = sequence_data.get("value", "")
        
        features = data.get("features", [])
        signal_end = 0
        mature_sequence = full_sequence
        
        for feature in features:
            if feature.get("type") == "Signal":
                try:
                    signal_end = int(feature["location"]["end"]["value"])
                except:
                    pass
                    
        if signal_end > 0 and signal_end < len(full_sequence):
            mature_sequence = full_sequence[signal_end:]

        return {
            "accession": accession,
            "full_sequence": full_sequence,
            "mature_sequence": mature_sequence,
            "signal_cleaved": signal_end > 0
        }
