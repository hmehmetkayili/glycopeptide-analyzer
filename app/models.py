from sqlalchemy import Column, Integer, String, Text, Boolean, ForeignKey, Enum, Float, Table
from sqlalchemy.orm import declarative_base, relationship
import enum

Base = declarative_base()

class ProteinType(enum.Enum):
    SERUM = "Serum Glycoprotein"
    MAB = "Monoclonal Antibody"
    OTHER = "Other"

class GylcoType(enum.Enum):
    N_LINKED = "N-Linked"
    O_LINKED = "O-Linked"

site_literature_link = Table('site_literature_link', Base.metadata,
    Column('site_id', Integer, ForeignKey('glycosylation_sites.id'), primary_key=True),
    Column('literature_id', Integer, ForeignKey('literature.id'), primary_key=True)
)

class Protein(Base):
    __tablename__ = "proteins"
    id = Column(Integer, primary_key=True, index=True)
    accession = Column(String(50), unique=True, index=True, nullable=True)
    name = Column(String(200), nullable=False)
    gene_name = Column(String(100))
    organism = Column(String(100), default="Homo sapiens (9606)")
    protein_type = Column(Enum(ProteinType), default=ProteinType.SERUM)
    chains = relationship("ProteinChain", back_populates="protein", cascade="all, delete-orphan")
    sites = relationship("GlycosylationSite", back_populates="protein", cascade="all, delete-orphan")

class ProteinChain(Base):
    __tablename__ = "protein_chains"
    id = Column(Integer, primary_key=True, index=True)
    protein_id = Column(Integer, ForeignKey("proteins.id"), nullable=False)
    chain_name = Column(String(100), nullable=False)
    sequence = Column(Text, nullable=False)
    molecular_weight_da = Column(Float, nullable=True)
    protein = relationship("Protein", back_populates="chains")
    sites = relationship("GlycosylationSite", back_populates="chain", cascade="all, delete-orphan")

class GlycosylationSite(Base):
    __tablename__ = "glycosylation_sites"
    id = Column(Integer, primary_key=True, index=True)
    protein_id = Column(Integer, ForeignKey("proteins.id"), nullable=False)
    chain_id = Column(Integer, ForeignKey("protein_chains.id"), nullable=True)
    position = Column(Integer, nullable=False)
    glyco_type = Column(Enum(GylcoType), nullable=False)
    motif = Column(String(10))
    local_sequence = Column(String(20))
    is_verified = Column(Boolean, default=False)
    evidence_type = Column(String(100))
    source = Column(String(100))
    protein = relationship("Protein", back_populates="sites")
    chain = relationship("ProteinChain", back_populates="sites")
    literature = relationship("Literature", secondary=site_literature_link, back_populates="sites")

class Literature(Base):
    __tablename__ = "literature"
    id = Column(Integer, primary_key=True, index=True)
    pubmed_id = Column(String(50), unique=True, nullable=True)
    doi = Column(String(100), unique=True, nullable=True)
    title = Column(Text, nullable=False)
    authors = Column(Text)
    year = Column(Integer)
    sites = relationship("GlycosylationSite", secondary=site_literature_link, back_populates="literature")
