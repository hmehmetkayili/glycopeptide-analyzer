from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field
from typing import List, Dict

from app.services.uniprot_service import UniProtService
from app.services.digestion_service import DigestionService
from app.services.glycosylation_service import GlycosylationService
from app.services.mass_calculator import MassCalculator

app = FastAPI(title="GlycoPeptide Analyzer v12")

@app.get("/", include_in_schema=False)
async def home():
    return FileResponse("/app/index.html")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_credentials=True, allow_methods=["*"], allow_headers=["*"])

uniprot_service = UniProtService()
digestion_service = DigestionService()

class RangeLimit(BaseModel):
    min: int
    max: int

class AnalysisRequest(BaseModel):
    sequence: str
    enzyme: str = "trypsin"
    missed_cleavages: int = 2
    pep_min_len: int = 4
    pep_max_len: int = 30
    static_mods: List[str] = ["carbamidomethyl_c"]
    variable_mods: List[str] = ["oxidation_m"]
    glycan_ranges: Dict[str, RangeLimit] = {
        "Hex": RangeLimit(min=3, max=5),
        "HexNAc": RangeLimit(min=2, max=4),
        "Fuc": RangeLimit(min=0, max=1),
        "Neu5Ac": RangeLimit(min=0, max=1)
    }

@app.get("/api/uniprot/{accession}")
async def fetch_protein_data(accession: str):
    data = await uniprot_service.fetch_protein_data(accession)
    if not data:
        raise HTTPException(status_code=404, detail="UniProt Accession not found.")
    return data

@app.post("/api/analyze")
async def run_full_analysis(req: AnalysisRequest):
    try:
        peptides = digestion_service.in_silico_digest(
            sequence=req.sequence, 
            enzyme=req.enzyme, 
            max_missed_cleavages=req.missed_cleavages,
            min_len=req.pep_min_len,
            max_len=req.pep_max_len,
            static_mods=req.static_mods,
            variable_mods=req.variable_mods
        )
        n_sites = GlycosylationService.find_n_linked_sites(req.sequence)
        mapped_peptides = GlycosylationService.map_peptides_to_sites(peptides, n_sites)
        
        ranges = {k: {"min": v.min, "max": v.max} for k, v in req.glycan_ranges.items()}
        glycan_library = MassCalculator.generate_n_glycan_library(ranges)
        
        theoretical_results = MassCalculator.generate_theoretical_glycopeptides(mapped_peptides, glycan_library)
        
        return {
            "status": "success",
            "metadata": {
                "enzyme": req.enzyme,
                "peptide_length_filter": f"{req.pep_min_len}-{req.pep_max_len}",
                "glycan_count_generated": len(glycan_library),
                "glycopeptides_found": len(theoretical_results)
            },
            "results": theoretical_results
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
