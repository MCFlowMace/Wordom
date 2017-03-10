%module wordom

%include fileio.i
%include tools.i
%include qcprot.i
%include xdrfile.i
%include xdrfile_xtc.i

%inline %{
void About()
{
  printf("PyWordom, version 0.23\n");
  return;
}
%}

/*
swig -features autodoc=1 -python wordom.i
gcc -fpic -D LAPACK -c wordom_wrap.c fileio.c xdrfile.c xdrfile_xtc.c tools.c -I`python-config --includes`
gcc -shared -llapack -lblas -lm wordom_wrap.o fileio.o xdrfile.o xdrfile_xtc.o tools.o -o _wordom.so
*/
