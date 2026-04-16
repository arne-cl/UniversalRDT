.PHONY: install docker \
        prepare-aracore prepare-metacyc \
        run-aracore run-metacyc \
        run-aracore-docker run-metacyc-docker \
        postprocess-metacyc compare-metacyc \
        clean clean-aracore clean-metacyc

RDT_JAR     := rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar
DOCKER_IMG  := universal-rdt

# ── Install (native, Ubuntu 24.04) ────────────────────────────────

install:
	sudo apt-get update
	sudo apt-get install -y --no-install-recommends \
		openjdk-21-jre-headless openbabel unzip

# ── Docker ─────────────────────────────────────────────────────────

docker: Dockerfile $(RDT_JAR)
	docker build -t $(DOCKER_IMG) .

# ── Data extraction ───────────────────────────────────────────────

prepare-aracore: AraCore/reaction_intermediates.zip
	unzip -o -O UTF-8 $< -d AraCore/

MetaCyc/reaction_intermediates.zip: MetaCyc/reaction_intermediates.zip.01 MetaCyc/reaction_intermediates.zip.02
	cat $^ > $@

prepare-metacyc: MetaCyc/reaction_intermediates.zip
	unzip -o -u -O UTF-8 $< -d MetaCyc/
	unzip -o -u -O UTF-8 MetaCyc/atom_mappings.zip -d MetaCyc/

# ── Run pipeline (native) ─────────────────────────────────────────

run-aracore: prepare-aracore
	cd AraCore && RDT_JAR=$(CURDIR)/$(RDT_JAR) bash run_rdt.sh
	cd AraCore && bash unite_mappings.sh

run-metacyc: prepare-metacyc
	cd MetaCyc && RDT_JAR=$(CURDIR)/$(RDT_JAR) bash run_rdt_metacyc.sh

postprocess-metacyc:
	cd MetaCyc && bash create_inchi_equivalent_atoms.sh
	cd MetaCyc && bash convert_existing_mappings.sh
	cd MetaCyc && bash sort_mappings.sh
	cd MetaCyc && bash unite_converted_mappings.sh
	cd MetaCyc && bash find_mapping_issues.sh
	cd MetaCyc && bash diff_converted_mappings.sh

compare-metacyc: run-metacyc postprocess-metacyc

# ── Run pipeline (Docker) ─────────────────────────────────────────

DOCKER_RUN = docker run --rm -v "$(CURDIR):/data" -w /data \
	-e RDT_JAR=/opt/rdt/$(RDT_JAR) \
	-e DISPLAY=:99 \
	$(DOCKER_IMG)

run-aracore-docker: docker prepare-aracore
	$(DOCKER_RUN) xvfb-run -a bash -c 'cd AraCore && bash run_rdt.sh && bash unite_mappings.sh'

run-metacyc-docker: docker prepare-metacyc
	$(DOCKER_RUN) xvfb-run -a bash -c 'cd MetaCyc && bash run_rdt_metacyc.sh'

# ── Clean ──────────────────────────────────────────────────────────

clean: clean-aracore clean-metacyc

clean-aracore:
	rm -rf AraCore/reaction_intermediates/
	rm -f AraCore/all_mapping.txt AraCore/all_mapping.sorted.txt \
	      AraCore/all_mapping.N.sorted.txt \
	      AraCore/all_atoms.N.sorted.txt AraCore/all_atoms.N.count.txt \
	      AraCore/all_atoms.N.count.histo \
	      AraCore/all_rxn_N_count.txt AraCore/all_rxn_N_count.histo
	rm -f AraCore/inchikey_diagnostic_results.csv

clean-metacyc:
	rm -rf MetaCyc/reaction_intermediates/ MetaCyc/atom_mappings/
	rm -f MetaCyc/reaction_intermediates.zip
