FROM maiolino_doc_v3

RUN mkdir scripts

COPY ./scripts /scripts

CMD ["bin/bash"]