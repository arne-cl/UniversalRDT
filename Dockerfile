FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-21-jre-headless \
    openbabel \
    unzip \
    xvfb \
    && rm -rf /var/lib/apt/lists/*

ENV RDT_JAR=/opt/rdt/rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar

COPY rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar /opt/rdt/rdt-2.5.0-SNAPSHOT-jar-with-dependencies.jar

WORKDIR /data
