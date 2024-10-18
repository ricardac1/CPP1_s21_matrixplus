# корневой image образ для сборки - чистый ubuntu
FROM ubuntu

# обновляем apt пакеты и устанавливаем по очереди нужные зависимости
RUN apt-get update && apt-get install -y \
    g++ \
    make \
    libgtest-dev \
    libsubunit-dev \
    valgrind \
    vim \
    cmake \
    git \
    clang \
    llvm \
    lcov \
    clang-format \
    libpcre3 \
    libpcre3-dev \
    build-essential \
    && apt-get clean

RUN git clone https://github.com/google/googletest.git /googletest && \
    cd /googletest && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install

COPY . /code

# Устанавливаем 
WORKDIR /code

CMD /bin/bash