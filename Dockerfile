# Environnement reproductible pour construire, tester et essayer SymAlgo++.
#
#   docker build -t symalgopp .            # compile et lance tous les tests
#   docker run --rm symalgopp              # démo
#   docker run --rm symalgopp exemple_ajustement
#   docker build --target construction -t symalgopp-dev .   # image de développement (g++, cmake)

FROM debian:12-slim AS construction
RUN apt-get update \
 && apt-get install -y --no-install-recommends g++ cmake ninja-build pkg-config libgmp-dev libeigen3-dev ca-certificates \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY . .
# la construction échoue si un test ou un cas d'usage échoue
RUN cmake -S . -B /build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/symalgopp \
 && cmake --build /build \
 && ctest --test-dir /build --output-on-failure \
 && cmake --install /build \
 && mkdir -p /opt/symalgopp/bin \
 && cp /build/demo /build/exemple_* /opt/symalgopp/bin/

FROM debian:12-slim AS execution
RUN apt-get update \
 && apt-get install -y --no-install-recommends libgmpxx4ldbl \
 && rm -rf /var/lib/apt/lists/*
COPY --from=construction /opt/symalgopp/bin /opt/symalgopp/bin
ENV PATH="/opt/symalgopp/bin:${PATH}"
CMD ["demo"]
