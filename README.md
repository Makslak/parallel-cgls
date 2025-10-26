# parallel-cgls
C++ 20
Release
Visual studio

	cd C:\path\to\source_code

build:
	Способ 1:
	
	cl /O2 /std:c++20 /EHsc main.cpp ParallelBlockNormalCG.cpp ^
  	/I"%MSMPI_INC%" /Fe:solver.exe ^
	/link /LIBPATH:"%MSMPI_LIB64%" msmpi.lib

Если не работает 1:
	Узнаем нужные пути:	
		Cпособ 1:
		
	echo %MSMPI_INC%
	echo %MSMPI_LIB64%
			
Вывод у меня:
	
	echo %MSMPI_INC%
	C:\Program Files (x86)\Microsoft SDKs\MPI\Include\

	echo %MSMPI_LIB64%
	C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64\
		
Способ 2:
	
	where /r "C:\Program Files (x86)" mpi.h
	where /r "C:\Program Files (x86)" msmpi.lib
		
Вывод у меня:

	where /r "C:\Program Files (x86)" mpi.h
	C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h

	where /r "C:\Program Files (x86)" msmpi.lib
	C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64\msmpi.lib
		
Узнав эти пути:
	%MSMPI_INC% = C:\Program Files (x86)\Microsoft SDKs\MPI\Include\
	%MSMPI_LIB64% = C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64\
Подставляем в команду из "Способ 1" и получаем:

	cl /O2 /std:c++20 /EHsc main.cpp ParallelBlockNormalCG.cpp ^
	/I"C:\Program Files (x86)\Microsoft SDKs\MPI\Include" /Fe:solver.exe ^
 	/link /LIBPATH:"C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64" msmpi.lib

Получаем solver.exe, и запускаем:

	mpiexec -n 4 solver.exe -M 1000 -N 1000 -alpha 1e-3 -iters 1000


Примеры:
C:\path\to\source_code>mpiexec -n 1 solver.exe -M 1000 -N 1000 -alpha 1e-3 -iters 1000
 M=1000 N=1000 alpha=0.001 iters=1000 procs=1 solve_ms=904.761

C:\path\to\source_code>mpiexec -n 4 solver.exe -M 1000 -N 1000 -alpha 1e-3 -iters 1000
 M=1000 N=1000 alpha=0.001 iters=1000 procs=4 solve_ms=251.225

C:\path\to\source_code>mpiexec -n 9 solver.exe -M 1000 -N 1000 -alpha 1e-3 -iters 1000
 M=1000 N=1000 alpha=0.001 iters=1000 procs=9 solve_ms=151.774

C:\path\to\source_code>mpiexec -n 16 solver.exe -M 1000 -N 1000 -alpha 1e-3 -iters 1000
 M=1000 N=1000 alpha=0.001 iters=1000 procs=16 solve_ms=3952.7

(Ryzen 5 3600 6 ядер / 12 потоков)
