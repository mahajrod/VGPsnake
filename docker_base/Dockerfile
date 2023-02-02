FROM mambaorg/micromamba:0.15.3
COPY --chown=micromamba:micromamba conda_spec.yaml /tmp/conda_spec.yaml
COPY pip_spec.txt /tmp/pip_spec.txt
RUN micromamba install -y -n base -f /tmp/conda_spec.yaml
RUN pip install -r /tmp/pip_spec.txt
