FROM interpace_array_base

USER root

RUN apt-get update
RUN apt-get install --assume-yes libicu63 openssl

WORKDIR /workdir
RUN chown micromamba:micromamba /workdir

USER micromamba

COPY --chown=micromamba:micromamba ./ ./


ENV PATH="/workdir/tools/:${PATH}"
ENV PATH="/workdir/resources/soft/:${PATH}"

#Set permissions
USER root
RUN chmod -R g+w /opt/conda/etc
RUN chmod -R g+w /home/micromamba
RUN chmod g+w /workdir
#bcl-convert binary MUST be present in resources/soft directory
#RUN chmod a+x /workdir/resources/soft/bcl-convert

USER micromamba

#ENTRYPOINT ["/workdir/entrypoint.sh"]
