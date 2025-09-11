# Dockerfile for unit testing calibratedcpp
FROM debian:stable

# Install dependencies
RUN apt-get update && \
    apt-get install -y openjdk-21-jdk openjfx ant junit4 && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /root

# Copy all source code
ADD . ./

# Default entrypoint: run tests
ENTRYPOINT JAVA_FX_HOME=/usr/share/java/ ant test