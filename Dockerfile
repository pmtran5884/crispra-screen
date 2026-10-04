# Turnkey image: `docker run` needs nothing on the host but Docker itself.
#
#   docker build -t crispra-screen .
#   docker run --rm -v "$PWD":/data crispra-screen qc --fastq /data/lib_R1.fastq.gz --outdir /data/qc_out
#
# Mount your working directory at /data; outputs land back on the host there.
FROM python:3.11-slim

LABEL org.opencontainers.image.title="crispra-screen" \
      org.opencontainers.image.source="https://github.com/pmtran5884/crispra-screen"

WORKDIR /opt/crispra-screen
COPY pyproject.toml README.md ./
COPY src ./src
RUN pip install --no-cache-dir .

# run as a non-root user and default the working directory to the bind mount
RUN useradd -m runner
USER runner
WORKDIR /data

ENTRYPOINT ["crispra-screen"]
CMD ["--help"]
