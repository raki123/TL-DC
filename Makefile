all: utils/htd src/nauty2_8_6/nauty.a src/clingo/CMakeLists.txt
	cmake -S . -B build
	cmake --build build -j 12
	cp build/tldc ./a.out

src/nauty2_8_6/nauty.a:
	cd src/nauty2_8_6/; ./configure --enable-tls --disable-interrupt; make -j12 nauty.a

utils/htd: src/htd_opt/CMakeLists.txt
	cmake -S src/htd_opt/ -B src/htd_opt/build -DBUILD_SHARED_LIBS=OFF
	cmake --build src/htd_opt/build -j 12
	cp src/htd_opt/build/bin/htd_main-1.2.0 utils/htd

src/htd_opt/CMakeLists.txt:
	git submodule update --init --recursive

src/clingo/CMakeLists.txt:
	git submodule update --init --recursive

