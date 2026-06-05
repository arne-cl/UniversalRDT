"""Preprocess the iML1515 metabolic model for RDT atom mapping.

Reads the iML1515 model (gzipped JSON, BiGG format), resolves SMILES and
InChIKeys for each metabolite via external APIs, and writes per-reaction
input folders that AraCore/run_rdt.py expects.
"""

import argparse
import gzip
import json
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

MODELSEED_URL = "https://modelseed.org/solr/compounds/select"
PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey"
CHEBI_URL = "https://www.ebi.ac.uk/chebi/backend/api/public/compound"


def bigg_to_bracket(bigg_id: str) -> tuple[str, str]:
    """Convert a BiGG metabolite ID to (base_name, bracket_notation).

    BiGG IDs use underscore-delimited compartments (dhap_c). RDT requires
    bracket-enclosed compartments (dhap[c]). Double underscores in the base
    name (e.g. cysi__L_c) are preserved as-is since the split is on the
    *last* underscore.

    Returns:
        Tuple of (base_name, bracket_id), e.g. ("dhap", "dhap[c]").
    """
    base, sep, compartment = bigg_id.rpartition("_")
    if not sep:
        return (bigg_id, bigg_id)
    return (base, f"{base}[{compartment}]")


def is_exchange_reaction(rxn: dict) -> bool:
    """True if reaction is an exchange or single-metabolite boundary reaction.

    Exchange reactions (EX_ prefix) and single-metabolite boundary reactions
    have no meaningful atom mapping and should be excluded.
    """
    rid = rxn["id"]
    if rid.startswith("EX_"):
        return True
    if len(rxn.get("metabolites", {})) == 1:
        return True
    return False


def is_biomass_reaction(rxn: dict) -> bool:
    """True if the reaction is a biomass pseudoreaction."""
    return "BIOMASS" in rxn["id"]


def has_all_smiles(rxn: dict, smiles_map: dict) -> bool:
    """True if every metabolite in the reaction has a SMILES entry."""
    return all(mid in smiles_map for mid in rxn.get("metabolites", {}))


def filter_reaction(rxn: dict, smiles_map: dict) -> bool:
    """True if the reaction should be included in RDT processing.

    Excludes exchange reactions, biomass pseudoreactions, and reactions
    where any metabolite lacks a SMILES string.
    """
    if is_exchange_reaction(rxn):
        return False
    if is_biomass_reaction(rxn):
        return False
    if not has_all_smiles(rxn, smiles_map):
        return False
    return True


def split_by_sign(metabolites: dict[str, float]) -> tuple[dict[str, float], dict[str, float]]:
    """Split reaction metabolites into reactants (negative coef) and products (positive coef).

    In the BiGG/COBRA convention, negative coefficients denote reactants
    and positive coefficients denote products. Coefficients are returned as
    absolute values.
    """
    reactants = {}
    products = {}
    for mid, coeff in metabolites.items():
        if coeff < 0:
            reactants[mid] = abs(coeff)
        elif coeff > 0:
            products[mid] = coeff
    return reactants, products


def flatten_coefficients(metabolites: dict[str, float]) -> list[str]:
    """Repeat each metabolite ID according to its stoichiometric coefficient.

    Non-unit stoichiometry (e.g. `pi_c: 2.0`) is flattened by repeating the
    species ID, since RDT needs each molecule instance as a separate entry.
    Coefficients are rounded to the nearest integer.
    """
    result = []
    for mid, coeff in metabolites.items():
        count = int(round(abs(coeff)))
        result.extend([mid] * count)
    return result


def build_smiles_string(reactant_smiles: list[str], product_smiles: list[str]) -> str:
    """Join reactant/product SMILES into the rxn.smiles format.

    Reactant SMILES are dot-joined on the left, product SMILES on the right,
    separated by >> (regardless of reaction reversibility).
    """
    left = ".".join(reactant_smiles)
    right = ".".join(product_smiles)
    return f"{left}>>{right}"


def build_smiles_from_reaction(rxn: dict, smiles_map: dict) -> str:
    """Build the full SMILES string for a reaction with flattened stoichiometry.

    Splits by coefficient sign, flattens any non-unit stoichiometry, and
    joins the resulting SMILES into the standard rxn.smiles format.
    """
    mets = rxn["metabolites"]
    reactants, products = split_by_sign(mets)
    flat_reactants = flatten_coefficients(reactants)
    flat_products = flatten_coefficients(products)
    reactant_smiles = [smiles_map[mid][0] for mid in flat_reactants]
    product_smiles = [smiles_map[mid][0] for mid in flat_products]
    return build_smiles_string(reactant_smiles, product_smiles)


def get_from_to_species(rxn: dict) -> tuple[list[str], list[str]]:
    """Get reactant (from) and product (to) species IDs in bracket notation.

    Returns flattened lists that may contain duplicates from non-unit
    stoichiometry. Callers typically deduplicate with set().
    """
    mets = rxn["metabolites"]
    reactants, products = split_by_sign(mets)
    flat_reactants = flatten_coefficients(reactants)
    flat_products = flatten_coefficients(products)
    from_species = [bigg_to_bracket(mid)[1] for mid in flat_reactants]
    to_species = [bigg_to_bracket(mid)[1] for mid in flat_products]
    return from_species, to_species


def build_inchikey_map(rxn: dict, smiles_map: dict) -> dict[str, str]:
    """Build a mapping from base name to InChIKey for all metabolites in a reaction.

    Base names exclude the compartment suffix. Each unique base name appears
    at most once even if the metabolite appears in multiple compartments.
    """
    result = {}
    for mid in rxn["metabolites"]:
        if mid in smiles_map:
            base = bigg_to_bracket(mid)[0]
            if base not in result:
                result[base] = smiles_map[mid][1]
    return result


def write_reaction_output(rxn: dict, smiles_map: dict, output_dir: Path) -> None:
    """Write the four input files that run_rdt.py expects for a single reaction.

    Creates a subdirectory named after the reaction ID containing:
      - rxn.smiles: SMILES string with flattened stoichiometry
      - from_species_with_cmp: sorted unique reactant species IDs (bracket)
      - to_species_with_cmp: sorted unique product species IDs (bracket)
      - species_id_inchikey.txt: tab-separated base_name -> InChIKey

    Reactions that fail filtering are silently skipped.
    """
    if not filter_reaction(rxn, smiles_map):
        return
    rxn_id = rxn["id"]
    rxn_dir = Path(output_dir) / rxn_id
    rxn_dir.mkdir(parents=True, exist_ok=True)

    rxn_smiles = build_smiles_from_reaction(rxn, smiles_map)
    (rxn_dir / "rxn.smiles").write_text(rxn_smiles + "\n")

    from_species, to_species = get_from_to_species(rxn)
    from_unique = sorted(set(from_species))
    to_unique = sorted(set(to_species))
    (rxn_dir / "from_species_with_cmp").write_text("\n".join(from_unique) + "\n")
    (rxn_dir / "to_species_with_cmp").write_text("\n".join(to_unique) + "\n")

    inchikey_map = build_inchikey_map(rxn, smiles_map)
    ik_lines = [f"{base}\t{inchikey}" for base, inchikey in sorted(inchikey_map.items())]
    (rxn_dir / "species_id_inchikey.txt").write_text("\n".join(ik_lines) + "\n")


def process_model(model: dict, smiles_map: dict, output_dir: Path) -> dict:
    """Run the full preprocessing pipeline over all model reactions.

    Filters each reaction, writes output for included ones, and returns a
    stats dict with counts for total, included, exchange-excluded,
    biomass-excluded, and no-SMILES-excluded reactions.
    """
    stats = {
        "total": len(model["reactions"]),
        "excluded_exchange": 0,
        "excluded_biomass": 0,
        "excluded_no_smiles": 0,
        "included": 0,
    }
    for rxn in model["reactions"]:
        if is_exchange_reaction(rxn):
            stats["excluded_exchange"] += 1
            continue
        if is_biomass_reaction(rxn):
            stats["excluded_biomass"] += 1
            continue
        if not has_all_smiles(rxn, smiles_map):
            stats["excluded_no_smiles"] += 1
            continue
        write_reaction_output(rxn, smiles_map, output_dir)
        stats["included"] += 1
    return stats


def load_model(path: str | Path) -> dict:
    """Load a gzipped BiGG JSON model from disk."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Model file not found: {path}")
    with gzip.open(path, "rt") as f:
        return json.load(f)


def load_cache(path: str | Path) -> dict:
    """Load the SMILES/InChIKey cache from a JSON file, or return empty dict."""
    path = Path(path)
    if not path.exists():
        return {}
    with open(path) as f:
        return json.load(f)


def save_cache(data: dict, path: str | Path) -> None:
    """Persist the SMILES/InChIKey cache to a JSON file (creates parent dirs)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def _json_get(url: str, params: dict | None = None) -> dict | None:
    """GET a JSON response from a URL, returning None on any failure.

    Appends query parameters, sets a User-Agent header, and applies a 10s
    timeout. Silently returns None for network, HTTP, or parse errors.
    """
    full_url = f"{url}?{urlencode(params)}" if params else url
    req = Request(full_url, headers={"User-Agent": "UniversalRDT-preprocessor/1.0"})
    try:
        with urlopen(req, timeout=10) as resp:
            return json.loads(resp.read())
    except (URLError, HTTPError, OSError, json.JSONDecodeError):
        return None


def fetch_smiles_modelseed(seed_id: str) -> dict | None:
    """Look up SMILES and InChIKey for a ModelSEED compound ID via Solr API.

    Primary source. Modelseed covers most of the iML1515 metabolites.
    """
    params = {"q": f"id:{seed_id}", "fl": "id,smiles,inchikey", "wt": "json"}
    data = _json_get(MODELSEED_URL, params)
    if data is None:
        return None
    docs = data.get("response", {}).get("docs", [])
    if docs and docs[0].get("smiles"):
        doc = docs[0]
        return {
            "smiles": doc["smiles"],
            "inchikey": doc.get("inchikey", ""),
            "source": "ModelSEED",
        }
    return None


def fetch_smiles_pubchem(inchi_key: str) -> dict | None:
    """Look up SMILES and InChIKey from PubChem REST API using an InChIKey.

    Fallback source for metabolites not found via ModelSEED.
    """
    url = f"{PUBCHEM_URL}/{inchi_key}/property/IsomericSMILES,InChIKey/JSON"
    data = _json_get(url)
    if data is None:
        return None
    props = data.get("PropertyTable", {}).get("Properties", [])
    if props and props[0].get("IsomericSMILES"):
        return {
            "smiles": props[0]["IsomericSMILES"],
            "inchikey": props[0].get("InChIKey", inchi_key),
            "source": "PubChem",
        }
    return None


def fetch_smiles_chebi(chebi_id: str) -> dict | None:
    """Look up SMILES and InChIKey from ChEBI REST API using a ChEBI ID.

    Last-resort fallback. Recovers 1 unique metabolite (moco_c) beyond
    PubChem.
    """
    url = f"{CHEBI_URL}/{chebi_id}/"
    data = _json_get(url)
    if data is None:
        return None
    struct = data.get("default_structure", {})
    if struct.get("smiles"):
        return {
            "smiles": struct["smiles"],
            "inchikey": struct.get("inchikey", ""),
            "source": "ChEBI",
        }
    return None


def resolve_metabolite(met: dict, cache: dict) -> dict | None:
    """Resolve SMILES and InChIKey for a single metabolite, trying sources in order.

    Priority chain (stop at first success):
      1. ModelSEED Solr API (via seed.compound annotation)
      2. PubChem REST API (via inchi_key annotation)
      3. ChEBI REST API (via chebi annotation)

    Checks the cache first. Successful results are written into cache so
    callers can persist it later.
    """
    mid = met["id"]
    ann = met.get("annotation", {})

    if mid in cache:
        cached = cache[mid]
        if cached.get("smiles") and cached.get("inchikey"):
            return cached

    seed_ids = ann.get("seed.compound", [])
    for seed_id in seed_ids:
        time.sleep(2)
        result = fetch_smiles_modelseed(seed_id)
        if result:
            cache[mid] = result
            return result

    inchi_keys = ann.get("inchi_key", [])
    for inchi_key in inchi_keys:
        time.sleep(2)
        result = fetch_smiles_pubchem(inchi_key)
        if result:
            cache[mid] = result
            return result

    chebi_ids = ann.get("chebi", [])
    for chebi_id in chebi_ids:
        time.sleep(2)
        result = fetch_smiles_chebi(chebi_id)
        if result:
            cache[mid] = result
            return result

    return None


def resolve_all_metabolites(model: dict, cache: dict) -> dict:
    """Build a smiles_map (met_id -> (smiles, inchikey)) for all model metabolites.

    Uses cached entries where available, otherwise resolves via the API
    priority chain. The cache dict is mutated in-place with any new results.
    """
    smiles_map = {}
    for met in model["metabolites"]:
        mid = met["id"]
        if mid in cache and cache[mid].get("smiles"):
            smiles_map[mid] = (cache[mid]["smiles"], cache[mid].get("inchikey", ""))
            continue
        result = resolve_metabolite(met, cache)
        if result:
            smiles_map[mid] = (result["smiles"], result.get("inchikey", ""))
    return smiles_map


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments, falling back to defaults for model path, output, and cache."""
    parser = argparse.ArgumentParser(
        description="Preprocess iML1515 model for RDT atom mapping"
    )
    default_model = Path.home() / "repos/atnlib/data/BiGG/iML1515.json.gz"
    parser.add_argument(
        "--model-path",
        type=Path,
        default=default_model,
        help="Path to gzipped iML1515 JSON model",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("iML1515/reaction_intermediates"),
        help="Output directory for per-reaction folders",
    )
    parser.add_argument(
        "--cache-file",
        type=Path,
        default=Path("iML1515/smiles_cache.json"),
        help="JSON cache file for SMILES/InChIKey results",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity (repeatable)",
    )
    return parser.parse_args(argv)


def main() -> None:
    """CLI entry point: load model, resolve SMILES, write per-reaction output."""
    args = parse_args()

    if args.verbose:
        print(f"Loading model from {args.model_path}...", file=sys.stderr)
    model = load_model(args.model_path)

    if args.verbose:
        print(f"Loading cache from {args.cache_file}...", file=sys.stderr)
    cache = load_cache(args.cache_file)

    if args.verbose:
        print("Resolving metabolite SMILES...", file=sys.stderr)
    smiles_map = resolve_all_metabolites(model, cache)

    if args.verbose:
        print(f"Saving cache ({len(cache)} entries)...", file=sys.stderr)
    save_cache(cache, args.cache_file)

    covered = len(smiles_map)
    total_mets = len(model["metabolites"])
    uncovered = total_mets - covered
    print(f"Metabolites: {covered}/{total_mets} covered, {uncovered} uncovered", file=sys.stderr)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stats = process_model(model, smiles_map, output_dir)
    print(
        f"Reactions: {stats['total']} total, "
        f"{stats['included']} included, "
        f"{stats['excluded_exchange']} exchange, "
        f"{stats['excluded_biomass']} biomass, "
        f"{stats['excluded_no_smiles']} no SMILES",
        file=sys.stderr,
    )
    print(f"Output written to {output_dir.resolve()}", file=sys.stderr)


if __name__ == "__main__":
    main()
