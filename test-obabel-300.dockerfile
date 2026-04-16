# Test Open Babel 3.0.0 vs 3.1.1 InChIKey generation
FROM ubuntu:20.04

# Install Open Babel 3.0.0 (Ubuntu 20.04 default)
RUN apt-get update && apt-get install -y openbabel && rm -rf /var/lib/apt/lists/*

# Copy test script
COPY test_inchikey_mismatch.sh /test.sh
RUN chmod +x /test.sh

# Run the test
CMD ["/bin/bash", "/test.sh"]
