cmake_minimum_required(VERSION 3.0.2)
project(SIDR2.0 C)

set(CMAKE_C_STANDARD 99)

message(STATUS "Build: Release")

set(SRC
        htslib/libhts.a
        src/kmer_hash.h src/kmer_hash.c
	src/param.h src/param.c
	src/pipeline.h src/pipeline.c
        src/pylink.h src/pylink.c
        src/sequence.h src/sequence.c
        src/sequencing_data.h src/sequencing_data.c
        src/taxonomy.h src/taxonomy.c
	src/main.c
   )

find_library(HTSLIB hts htslib)
find_package (Python3 COMPONENTS Development REQUIRED)
include_directories(src htslib ${Python3_INCLUDE_DIR})
add_executable(SIDR2.0 ${SRC})
target_link_libraries(SIDR2.0 PUBLIC ${HTSLIB} ${Python3_LIBRARY} pthread z m bz2 dl util)
target_compile_options(SIDR2.0 PUBLIC -O2 -g)

install(TARGETS SIDR2.0 DESTINATION .)
