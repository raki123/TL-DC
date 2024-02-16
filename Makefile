all: utils/htd utils/bw utils/flow_cutter src/nauty2_8_6/nauty.a src/jemalloc/lib/libjemalloc.a src/tamaki/tw/heuristic/MainDecomposer$$1.class
	cmake -S . -B build
	cmake --build build -j 12
	cp build/tldc ./a.out

src/nauty2_8_6/nauty.a:
	cd src/nauty2_8_6/; CFLAGS="-Wno-incompatible-pointer-types -O4" ./configure --enable-tls --disable-interrupt; make -j12 nauty.a

src/jemalloc/lib/libjemalloc.a: src/jemalloc/autogen.sh
	cd src/jemalloc/; ./autogen.sh; make -j12

utils/htd: src/htd_opt/CMakeLists.txt
	cmake -S src/htd_opt/ -B src/htd_opt/build -DBUILD_SHARED_LIBS=OFF
	cmake --build src/htd_opt/build -j 12
	cp src/htd_opt/build/bin/htd_main-1.2.0 utils/htd

utils/bw:
	cd src/hicks/; make -j 12
	cp src/hicks/bw utils/bw

utils/flow_cutter: src/flow_cutter/build.sh
	cd src/flow_cutter; bash build.sh
	cp src/flow_cutter/flow_cutter_pace17 utils/flow_cutter

src/tamaki/tw/heuristic/MainDecomposer$$1.class: src/tamaki/tw/heuristic/MainDecomposer.class
	cd src/tamaki; make

src/htd_opt/CMakeLists.txt:
	git submodule update --init --recursive

src/flow_cutter/build.sh:
	git submodule update --init --recursive

src/jemalloc/autogen.sh:
	git submodule update --init --recursive

src/tamaki/tw/heuristic/MainDecomposer.class:
	git submodule update --init --recursive;
