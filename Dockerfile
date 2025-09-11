# Dockerfile for unit testing calibratedcpp
FROM debian:stable

# Install dependencies
RUN apt-get update && \
    apt-get install -y openjdk-21-jdk openjfx ant junit4 unzip wget && \
    rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /root

# Copy source, test, and build.xml
COPY src ./src
COPY test ./test
COPY build.xml ./build.xml

# Default entrypoint: run tests
ENTRYPOINT JAVA_FX_HOME=/usr/share/java/ant ant test
