# syntax=docker/dockerfile:1.7

FROM debian:bookworm-slim AS downloader

ARG MODULE_VERSION
LABEL module.version="${MODULE_VERSION}"
RUN echo "Building version ${MODULE_VERSION}"

COPY images/pandepth/VERSION /tmp/pandepth.VERSION
COPY images/pandepth/SHA256SUMS /tmp/pandepth.SHA256SUMS

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        curl && \
    rm -rf /var/lib/apt/lists/* && \
    PANDEPTH_VERSION="$(tr -d '[:space:]' < /tmp/pandepth.VERSION)" && \
    case "${PANDEPTH_VERSION}" in ""|*[!0-9.]*) echo "Invalid pandepth version: ${PANDEPTH_VERSION}" >&2; exit 1 ;; esac && \
    PANDEPTH_TARBALL="PanDepth-${PANDEPTH_VERSION}-Linux-x86_64.tar.gz" && \
    curl -fsSL "https://github.com/HuiyangYu/PanDepth/releases/download/v${PANDEPTH_VERSION}/${PANDEPTH_TARBALL}" \
        -o "/tmp/${PANDEPTH_TARBALL}" && \
    echo "$(cat /tmp/pandepth.SHA256SUMS)  /tmp/${PANDEPTH_TARBALL}" | sha256sum -c && \
    tar -xzf "/tmp/${PANDEPTH_TARBALL}" -C /opt && \
    rm "/tmp/${PANDEPTH_TARBALL}"

FROM debian:bookworm-slim

COPY --from=downloader /opt/PanDepth-2.26-Linux-x86_64/pandepth /opt/pandepth/bin/pandepth
COPY images/pandepth/VERSION /opt/pandepth/VERSION

ENV PATH="/opt/pandepth/bin:${PATH}"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libz1g \
        libstdc++6 \
        procps && \
    rm -rf /var/lib/apt/lists/* && \
    pandepth --help > /tmp/pandepth.help 2>&1 && \
    grep -qi "pandepth" /tmp/pandepth.help && \
    rm /tmp/pandepth.help

WORKDIR /data

CMD ["bash"]
