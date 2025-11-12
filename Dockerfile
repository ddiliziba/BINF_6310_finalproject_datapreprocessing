
FROM mambaorg/micromamba:1.5.8

USER root
RUN apt-get update && apt-get install -y pigz && rm -rf /var/lib/apt/lists/*

# bring the env file into the image
COPY environment.yml /tmp/environment.yml

SHELL ["/bin/bash", "-lc"]
# build the env (create + install from the YAML)
RUN micromamba create -y -n hap_phase -f /tmp/environment.yml && micromamba clean --all -y

# auto-activate the env in interactive shells
ENV MAMBA_DOCKERFILE_ACTIVATE=1
RUN echo "micromamba activate hap_phase" >> ~/.bashrc

WORKDIR /data
ENTRYPOINT ["/bin/bash", "-lc"]
CMD ["bash"]
