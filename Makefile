all: utils/htd src/nauty2_8_6/nauty.a
	cmake -S . -B build -DCMAKE_CXX_COMPILER=g++
	cmake --build build -j 12
	cp build/read ./a.out

src/nauty2_8_6/nauty.a:
	cd src/nauty2_8_6/; ./configure --enable-tls --disable-interrupt; make -j12 nauty.a

utils/htd:
	cmake -S src/htd_opt/ -B src/htd_opt/build -DBUILD_SHARED_LIBS=OFF
	cmake --build src/htd_opt/build -j 12
	cp src/htd_opt/build/bin/htd_main-1.2.0 utils/htd


#utils/htd:
#	cmake -S src/htd/ -B src/htd/build
#	cmake --build src/htd/build -j 12
#	cp src/htd/build/bin/htd_main-1.2.0 utils/htd