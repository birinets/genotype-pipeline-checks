# The image that workflows/genotype_checks.wdl runs the check in: the kit installed from this
# repository, so the workflow needs no network access when it runs. .github/workflows/image.yml
# builds and pushes it for release tags only.
#
# python:3.12-slim, pinned by digest (3.12.16-slim-trixie on 8 November 2026).
FROM python:3.12-slim@sha256:5e2bdb1b0d542d5039cffb8181c65930093949194a65c365dfb1b0fe13a4fa5b

LABEL org.opencontainers.image.source="https://github.com/birinets/genotype-pipeline-checks" \
      org.opencontainers.image.description="genotype-pipeline-checks: allele-identity and strand checks for genotype pipelines" \
      org.opencontainers.image.licenses="Apache-2.0"

COPY pyproject.toml README.md LICENSE NOTICE /tmp/kit/
COPY src/ /tmp/kit/src/
RUN python -m pip install --no-cache-dir --disable-pip-version-check /tmp/kit \
 && rm -rf /tmp/kit \
 && genotype-checks selftest > /dev/null

CMD ["genotype-checks", "--help"]
